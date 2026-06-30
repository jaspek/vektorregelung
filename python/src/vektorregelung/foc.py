"""Indirekte feldorientierte Regelung (IFOC) der Asynchronmaschine.

Indirect field-oriented control of an induction motor. The controller performs,
once per sample:

1. Clarke + Park transform of the measured stator phase currents into the
   rotor-flux frame, using the *estimated* field angle ``theta_e``.
2. An outer speed PI loop that produces the torque-current reference ``iq*``.
3. A current-model rotor-flux estimator and the slip relation, which together
   give the field angle by integration (this is what makes the scheme
   *indirect* -- no flux sensor is required).
4. Two inner current PI loops (d and q) with cross-coupling and back-EMF
   feed-forward decoupling, producing the stator voltage command.
5. Voltage-vector limiting, inverse Park and space-vector PWM to obtain the
   three half-bridge duty cycles.

The ``d`` axis is aligned with the rotor flux, so ``id`` sets the magnetising
flux and ``iq`` sets the torque -- exactly as in a separately excited DC machine.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

from .induction_motor import InductionMotorParameters
from .pi_controller import PIController
from .svpwm import max_voltage, space_vector_pwm
from .transforms import clarke, inverse_park, park, wrap_angle


@dataclass
class FocOutput:
    """Everything a single controller step produces (handy for logging)."""

    duty_a: float
    duty_b: float
    duty_c: float
    v_alpha: float
    v_beta: float
    v_d: float
    v_q: float
    i_d: float
    i_q: float
    id_ref: float
    iq_ref: float
    theta_e: float
    omega_slip: float
    flux_estimate: float


@dataclass
class FocLimits:
    """Operating limits enforced by the controller."""

    v_dc: float          # DC-link voltage [V]
    i_max: float         # maximum stator current magnitude [A]
    flux_min: float = 1e-3  # floor for the flux estimate to avoid divide-by-zero


class IndirectFOC:
    """Indirect field-oriented controller for an induction machine."""

    def __init__(
        self,
        motor: InductionMotorParameters,
        dt: float,
        limits: FocLimits,
        rotor_flux_ref: float,
        id_pi: PIController,
        iq_pi: PIController,
        speed_pi: PIController,
    ):
        self.motor = motor
        self.dt = dt
        self.limits = limits
        self.rotor_flux_ref = rotor_flux_ref

        self.id_pi = id_pi
        self.iq_pi = iq_pi
        self.speed_pi = speed_pi

        # Estimator / integrator state.
        self.theta_e = 0.0          # estimated field angle [rad]
        self.flux_estimate = 0.0    # current-model rotor flux estimate [Wb]

        # Magnetising-current reference that establishes the rotor flux.
        # id* = psi_r* / Lm  (steady-state relation lambda_r = Lm * id).
        self.id_ref = rotor_flux_ref / motor.lm

    # --------------------------------------------------------------- factory
    @classmethod
    def from_motor(
        cls,
        motor: InductionMotorParameters,
        dt: float,
        v_dc: float,
        rotor_flux_ref: float,
        i_max: float,
        current_bandwidth_hz: float = 200.0,
        speed_bandwidth_hz: float = 15.0,
    ) -> "IndirectFOC":
        """Build a controller with PI gains tuned from the motor parameters.

        Current loops are tuned by internal-model control (pole placement on the
        stator transient circuit); the speed loop places its zero a decade below
        the crossover for a robust phase margin.
        """
        # --- Current loop: plant 1/(sigmaLs*s + Rsigma) -> bandwidth w_ci -----
        w_ci = 2.0 * math.pi * current_bandwidth_hz
        kp_i = motor.sigma_ls * w_ci
        ki_i = motor.r_sigma * w_ci

        v_lim = max_voltage(v_dc)
        id_pi = PIController(kp_i, ki_i, out_min=-v_lim, out_max=v_lim)
        iq_pi = PIController(kp_i, ki_i, out_min=-v_lim, out_max=v_lim)

        # --- Speed loop: Te = Kt*iq, plant Kt/(J*s) -------------------------
        w_cs = 2.0 * math.pi * speed_bandwidth_hz
        kt = motor.torque_constant(rotor_flux_ref)
        kp_s = motor.inertia * w_cs / kt
        ki_s = kp_s * w_cs / 10.0
        speed_pi = PIController(kp_s, ki_s, out_min=-i_max, out_max=i_max)

        limits = FocLimits(v_dc=v_dc, i_max=i_max)
        return cls(
            motor=motor,
            dt=dt,
            limits=limits,
            rotor_flux_ref=rotor_flux_ref,
            id_pi=id_pi,
            iq_pi=iq_pi,
            speed_pi=speed_pi,
        )

    # ------------------------------------------------------------------ reset
    def reset(self) -> None:
        self.theta_e = 0.0
        self.flux_estimate = 0.0
        self.id_pi.reset()
        self.iq_pi.reset()
        self.speed_pi.reset()

    # ------------------------------------------------------------------- step
    def update(
        self,
        speed_ref: float,
        ia: float,
        ib: float,
        ic: float,
        omega_m: float,
    ) -> FocOutput:
        """Run one control step.

        Parameters
        ----------
        speed_ref:
            Mechanical speed setpoint [rad/s].
        ia, ib, ic:
            Measured stator phase currents [A].
        omega_m:
            Measured mechanical shaft speed [rad/s].
        """
        motor = self.motor
        dt = self.dt

        # 1) Measure: abc -> alpha/beta -> d/q (rotor-flux frame).
        i_alpha, i_beta = clarke(ia, ib, ic)
        i_d, i_q = park(i_alpha, i_beta, self.theta_e)

        # 2) Rotor-flux estimator (current model):
        #    d(psi)/dt = (Lm*id - psi) / Tr.
        self.flux_estimate += dt * (motor.lm * i_d - self.flux_estimate) / motor.tr
        flux = max(self.flux_estimate, self.limits.flux_min)

        # 3) Slip relation and field-angle integration.
        #    omega_slip = Lm * iq / (Tr * psi_r).
        omega_slip = motor.lm * i_q / (motor.tr * flux)
        omega_e = motor.pole_pairs * omega_m + omega_slip
        self.theta_e = wrap_angle(self.theta_e + dt * omega_e)

        # 4) Outer speed loop -> torque-current reference.
        iq_ref = self.speed_pi.update(speed_ref - omega_m, dt)
        id_ref = self.id_ref
        iq_ref = self._limit_iq(id_ref, iq_ref)

        # 5) Inner current loops.
        v_d = self.id_pi.update(id_ref - i_d, dt)
        v_q = self.iq_pi.update(iq_ref - i_q, dt)

        # 6) Cross-coupling / back-EMF feed-forward decoupling.
        #    vd_ff = -omega_e * sigmaLs * iq
        #    vq_ff =  omega_e * (sigmaLs * id + (Lm/Lr) * psi_r)
        sigma_ls = motor.sigma_ls
        lm_over_lr = motor.lm / motor.lr
        v_d += -omega_e * sigma_ls * i_q
        v_q += omega_e * (sigma_ls * i_d + lm_over_lr * flux)

        # 7) Limit the voltage vector to the SVPWM linear region.
        v_d, v_q = self._limit_voltage(v_d, v_q)

        # 8) Inverse Park + space-vector PWM.
        v_alpha, v_beta = inverse_park(v_d, v_q, self.theta_e)
        duty_a, duty_b, duty_c = space_vector_pwm(v_alpha, v_beta, self.limits.v_dc)

        return FocOutput(
            duty_a=duty_a,
            duty_b=duty_b,
            duty_c=duty_c,
            v_alpha=v_alpha,
            v_beta=v_beta,
            v_d=v_d,
            v_q=v_q,
            i_d=i_d,
            i_q=i_q,
            id_ref=id_ref,
            iq_ref=iq_ref,
            theta_e=self.theta_e,
            omega_slip=omega_slip,
            flux_estimate=self.flux_estimate,
        )

    # --------------------------------------------------------------- helpers
    def _limit_iq(self, id_ref: float, iq_ref: float) -> float:
        """Clamp iq* so that the total current vector stays within ``i_max``.

        The magnetising current id* has priority (flux must be maintained), so
        the torque current gets whatever headroom is left.
        """
        i_max = self.limits.i_max
        head = i_max * i_max - id_ref * id_ref
        if head <= 0.0:
            return 0.0
        iq_limit = math.sqrt(head)
        if iq_ref > iq_limit:
            return iq_limit
        if iq_ref < -iq_limit:
            return -iq_limit
        return iq_ref

    def _limit_voltage(self, v_d: float, v_q: float) -> tuple[float, float]:
        """Scale the (vd, vq) vector down to the maximum realisable magnitude."""
        v_max = max_voltage(self.limits.v_dc)
        mag = math.hypot(v_d, v_q)
        if mag > v_max and mag > 0.0:
            scale = v_max / mag
            return v_d * scale, v_q * scale
        return v_d, v_q

"""Dynamisches Modell einer Drehstrom-Asynchronmaschine (Käfigläufer).

Dynamic model of a three-phase squirrel-cage induction motor in the *stationary*
(alpha/beta) reference frame. State vector:

    x = [i_s_alpha, i_s_beta, psi_r_alpha, psi_r_beta, omega_m, theta_m]

where the stator currents and rotor fluxes are stationary-frame space-vector
components, ``omega_m`` is the mechanical angular speed [rad/s] and ``theta_m``
the mechanical rotor angle [rad]. Integration uses a fixed-step RK4 scheme.

The model is intended for closed-loop simulation of the field-oriented
controller; it is the "plant" the controller does not have direct access to
(the controller only sees stator phase currents and the measured shaft speed).
"""

from __future__ import annotations

from dataclasses import dataclass

from .transforms import inverse_clarke


@dataclass(frozen=True)
class InductionMotorParameters:
    """Per-phase equivalent-circuit parameters of an induction machine.

    All rotor quantities are referred to the stator side.
    """

    rs: float          # stator resistance [Ohm]
    rr: float          # rotor resistance [Ohm]
    ls: float          # stator inductance  Ls = Lls + Lm [H]
    lr: float          # rotor inductance   Lr = Llr + Lm [H]
    lm: float          # mutual (magnetising) inductance [H]
    pole_pairs: int    # number of pole pairs
    inertia: float     # rotor inertia J [kg*m^2]
    friction: float = 0.0  # viscous friction coefficient B [N*m*s]

    # --- Derived quantities -------------------------------------------------
    @property
    def sigma(self) -> float:
        """Total leakage factor sigma = 1 - Lm^2 / (Ls*Lr)."""
        return 1.0 - (self.lm * self.lm) / (self.ls * self.lr)

    @property
    def sigma_ls(self) -> float:
        """Transient (leakage) inductance sigma*Ls seen by the current loop."""
        return self.sigma * self.ls

    @property
    def tr(self) -> float:
        """Rotor time constant Tr = Lr / Rr [s]."""
        return self.lr / self.rr

    @property
    def r_sigma(self) -> float:
        """Equivalent series resistance for current-loop tuning.

        Rsigma = Rs + Rr*(Lm/Lr)^2 is the effective damping resistance of the
        stator transient circuit.
        """
        ratio = self.lm / self.lr
        return self.rs + self.rr * ratio * ratio

    def torque_constant(self, rotor_flux: float) -> float:
        """Torque per q-axis current Kt = 1.5*p*(Lm/Lr)*psi_r [N*m/A]."""
        return 1.5 * self.pole_pairs * (self.lm / self.lr) * rotor_flux


class InductionMotor:
    """Stationary-frame induction-motor simulation model."""

    def __init__(self, params: InductionMotorParameters):
        self.p = params
        # State: stator current (alpha/beta), rotor flux (alpha/beta),
        # mechanical speed, mechanical angle.
        self.i_alpha = 0.0
        self.i_beta = 0.0
        self.psi_r_alpha = 0.0
        self.psi_r_beta = 0.0
        self.omega_m = 0.0
        self.theta_m = 0.0

    # ------------------------------------------------------------------ state
    @property
    def state(self) -> tuple[float, float, float, float, float, float]:
        return (
            self.i_alpha,
            self.i_beta,
            self.psi_r_alpha,
            self.psi_r_beta,
            self.omega_m,
            self.theta_m,
        )

    def phase_currents(self) -> tuple[float, float, float]:
        """Stator currents (ia, ib, ic) reconstructed from the alpha/beta state."""
        return inverse_clarke(self.i_alpha, self.i_beta)

    @property
    def rotor_flux_magnitude(self) -> float:
        return (self.psi_r_alpha ** 2 + self.psi_r_beta ** 2) ** 0.5

    @property
    def electromagnetic_torque(self) -> float:
        p = self.p
        return (
            1.5
            * p.pole_pairs
            * (p.lm / p.lr)
            * (self.psi_r_alpha * self.i_beta - self.psi_r_beta * self.i_alpha)
        )

    # ------------------------------------------------------------- dynamics
    def _derivatives(
        self,
        state: tuple[float, float, float, float, float, float],
        v_alpha: float,
        v_beta: float,
        load_torque: float,
    ) -> tuple[float, float, float, float, float, float]:
        i_a, i_b, psi_a, psi_b, omega_m, _theta_m = state
        p = self.p

        sigma_ls = p.sigma_ls
        tr_inv = p.rr / p.lr            # 1 / Tr
        lm_over_lr = p.lm / p.lr
        omega_r = p.pole_pairs * omega_m  # electrical rotor speed

        # Effective stator resistance term: Rs + Rr*(Lm/Lr)^2.
        r_eff = p.rs + p.rr * lm_over_lr * lm_over_lr

        # Stator current dynamics (stationary frame).
        # di/dt = (1/sigmaLs) * [ -r_eff*i + (Lm/Lr)*(1/Tr)*psi
        #                          + (Lm/Lr)*omega_r*psi_perp + v ]
        di_a = (
            -r_eff * i_a
            + lm_over_lr * tr_inv * psi_a
            + lm_over_lr * omega_r * psi_b
            + v_alpha
        ) / sigma_ls
        di_b = (
            -r_eff * i_b
            + lm_over_lr * tr_inv * psi_b
            - lm_over_lr * omega_r * psi_a
            + v_beta
        ) / sigma_ls

        # Rotor flux dynamics.
        dpsi_a = p.lm * tr_inv * i_a - tr_inv * psi_a - omega_r * psi_b
        dpsi_b = p.lm * tr_inv * i_b - tr_inv * psi_b + omega_r * psi_a

        # Electromagnetic torque and mechanical dynamics.
        torque = 1.5 * p.pole_pairs * lm_over_lr * (psi_a * i_b - psi_b * i_a)
        domega = (torque - load_torque - p.friction * omega_m) / p.inertia
        dtheta = omega_m

        return di_a, di_b, dpsi_a, dpsi_b, domega, dtheta

    def step(self, v_alpha: float, v_beta: float, load_torque: float, dt: float) -> None:
        """Advance the model by ``dt`` seconds using classical RK4."""
        s0 = self.state

        k1 = self._derivatives(s0, v_alpha, v_beta, load_torque)
        s1 = _axpy(s0, k1, dt * 0.5)
        k2 = self._derivatives(s1, v_alpha, v_beta, load_torque)
        s2 = _axpy(s0, k2, dt * 0.5)
        k3 = self._derivatives(s2, v_alpha, v_beta, load_torque)
        s3 = _axpy(s0, k3, dt)
        k4 = self._derivatives(s3, v_alpha, v_beta, load_torque)

        new_state = tuple(
            s0[i] + (dt / 6.0) * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i])
            for i in range(6)
        )
        (
            self.i_alpha,
            self.i_beta,
            self.psi_r_alpha,
            self.psi_r_beta,
            self.omega_m,
            self.theta_m,
        ) = new_state


def _axpy(base, deriv, h):
    """Return base + h*deriv element-wise for the 6-element state tuple."""
    return tuple(base[i] + h * deriv[i] for i in range(6))

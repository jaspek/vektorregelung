"""Closed-loop simulation harness: IndirectFOC driving an InductionMotor.

Couples the controller and the plant at a fixed control rate and records the
trajectory. Reference speed and load torque are supplied as callables of time so
arbitrary test profiles (steps, ramps, reversals) can be expressed.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable, List

from .foc import IndirectFOC
from .induction_motor import InductionMotor


@dataclass
class SimulationLog:
    """Time series recorded during a simulation run."""

    t: List[float] = field(default_factory=list)
    speed_ref: List[float] = field(default_factory=list)
    speed: List[float] = field(default_factory=list)
    torque: List[float] = field(default_factory=list)
    load_torque: List[float] = field(default_factory=list)
    i_d: List[float] = field(default_factory=list)
    i_q: List[float] = field(default_factory=list)
    iq_ref: List[float] = field(default_factory=list)
    flux_actual: List[float] = field(default_factory=list)
    flux_estimate: List[float] = field(default_factory=list)
    theta_e: List[float] = field(default_factory=list)
    duty_a: List[float] = field(default_factory=list)
    duty_b: List[float] = field(default_factory=list)
    duty_c: List[float] = field(default_factory=list)


def simulate(
    controller: IndirectFOC,
    motor: InductionMotor,
    speed_ref: Callable[[float], float],
    load_torque: Callable[[float], float],
    duration: float,
    dt: float,
) -> SimulationLog:
    """Run a closed-loop simulation.

    Parameters
    ----------
    controller, motor:
        The controller and plant. ``controller.dt`` should equal ``dt``.
    speed_ref, load_torque:
        Functions of time ``t`` returning the speed setpoint [rad/s] and load
        torque [N*m].
    duration, dt:
        Total simulated time and fixed control/integration step [s].
    """
    log = SimulationLog()
    n_steps = int(round(duration / dt))

    for k in range(n_steps):
        t = k * dt
        w_ref = speed_ref(t)
        t_load = load_torque(t)

        ia, ib, ic = motor.phase_currents()
        out = controller.update(w_ref, ia, ib, ic, motor.omega_m)

        # The averaged inverter reproduces the commanded (alpha, beta) voltage,
        # so we apply it directly to the plant.
        motor.step(out.v_alpha, out.v_beta, t_load, dt)

        log.t.append(t)
        log.speed_ref.append(w_ref)
        log.speed.append(motor.omega_m)
        log.torque.append(motor.electromagnetic_torque)
        log.load_torque.append(t_load)
        log.i_d.append(out.i_d)
        log.i_q.append(out.i_q)
        log.iq_ref.append(out.iq_ref)
        log.flux_actual.append(motor.rotor_flux_magnitude)
        log.flux_estimate.append(out.flux_estimate)
        log.theta_e.append(out.theta_e)
        log.duty_a.append(out.duty_a)
        log.duty_b.append(out.duty_b)
        log.duty_c.append(out.duty_c)

    return log

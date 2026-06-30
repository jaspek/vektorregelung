import math

import pytest

from vektorregelung import (
    IndirectFOC,
    InductionMotor,
    InductionMotorParameters,
    simulate,
)

MOTOR = InductionMotorParameters(
    rs=1.5,
    rr=1.2,
    ls=0.14,
    lr=0.14,
    lm=0.135,
    pole_pairs=2,
    inertia=0.02,
    friction=0.001,
)

DT = 1.0e-4


def _build():
    controller = IndirectFOC.from_motor(
        MOTOR, dt=DT, v_dc=600.0, rotor_flux_ref=1.0, i_max=20.0
    )
    motor = InductionMotor(MOTOR)
    return controller, motor


def test_motor_parameter_derivations():
    # sigma = 1 - Lm^2/(Ls*Lr)
    expected_sigma = 1.0 - 0.135**2 / (0.14 * 0.14)
    assert MOTOR.sigma == pytest.approx(expected_sigma)
    assert MOTOR.tr == pytest.approx(0.14 / 1.2)
    assert MOTOR.sigma_ls == pytest.approx(expected_sigma * 0.14)


def test_flux_establishes_during_magnetisation():
    controller, motor = _build()
    log = simulate(
        controller,
        motor,
        speed_ref=lambda t: 0.0,
        load_torque=lambda t: 0.0,
        duration=0.6,
        dt=DT,
    )
    # After several rotor time constants the flux should reach the reference.
    assert log.flux_actual[-1] == pytest.approx(1.0, rel=0.05)
    # The current-model estimate must track the real flux closely.
    assert log.flux_estimate[-1] == pytest.approx(log.flux_actual[-1], rel=0.05)
    # Speed must stay at zero with no torque command.
    assert abs(log.speed[-1]) < 1.0


def test_speed_tracking_step():
    controller, motor = _build()
    log = simulate(
        controller,
        motor,
        speed_ref=lambda t: 0.0 if t < 0.3 else 120.0,
        load_torque=lambda t: 0.0,
        duration=1.0,
        dt=DT,
    )
    assert log.speed[-1] == pytest.approx(120.0, abs=1.0)


def test_load_rejection():
    controller, motor = _build()
    load = 8.0

    def load_profile(t):
        return load if t >= 0.8 else 0.0

    log = simulate(
        controller,
        motor,
        speed_ref=lambda t: 0.0 if t < 0.3 else 120.0,
        load_torque=load_profile,
        duration=1.6,
        dt=DT,
    )
    # Speed recovers to the setpoint despite the load.
    assert log.speed[-1] == pytest.approx(120.0, abs=1.5)
    # Steady-state torque balances the applied load (plus a little friction).
    assert log.torque[-1] == pytest.approx(load, abs=1.0)


def test_current_limit_respected():
    controller, motor = _build()
    log = simulate(
        controller,
        motor,
        # Aggressive step forces the q-current into the limit during accel.
        speed_ref=lambda t: 0.0 if t < 0.3 else 200.0,
        load_torque=lambda t: 0.0,
        duration=0.8,
        dt=DT,
    )
    i_max = 20.0
    peak = max(math.hypot(d, q) for d, q in zip(log.i_d, log.i_q))
    # Allow a small numerical/transient margin above the limit.
    assert peak <= i_max * 1.05


def test_speed_reversal():
    controller, motor = _build()
    log = simulate(
        controller,
        motor,
        speed_ref=lambda t: 0.0 if t < 0.3 else (100.0 if t < 1.0 else -100.0),
        load_torque=lambda t: 0.0,
        duration=1.8,
        dt=DT,
    )
    assert log.speed[-1] == pytest.approx(-100.0, abs=2.0)

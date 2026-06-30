import pytest

from vektorregelung import PIController


def test_proportional_only():
    pi = PIController(kp=2.0, ki=0.0)
    assert pi.update(3.0, 0.001) == pytest.approx(6.0)
    # No integral action -> output depends only on the present error.
    assert pi.update(0.0, 0.001) == pytest.approx(0.0)


def test_integral_accumulates_to_remove_steady_error():
    pi = PIController(kp=0.0, ki=10.0)
    dt = 0.01
    out = 0.0
    for _ in range(100):
        out = pi.update(1.0, dt)
    # integral of error over 1 s with ki=10 -> ~10.
    assert out == pytest.approx(10.0, rel=1e-6)


def test_output_saturation():
    pi = PIController(kp=5.0, ki=0.0, out_min=-1.0, out_max=1.0)
    assert pi.update(100.0, 0.001) == pytest.approx(1.0)
    assert pi.update(-100.0, 0.001) == pytest.approx(-1.0)


def test_anti_windup_prevents_integrator_buildup():
    # With anti-windup the integrator must not keep growing while saturated, so
    # the output recovers immediately once the error reverses sign.
    pi = PIController(kp=1.0, ki=1000.0, out_min=-1.0, out_max=1.0)
    dt = 0.001
    for _ in range(1000):  # drive hard into the upper limit
        out = pi.update(10.0, dt)
    assert out == pytest.approx(1.0)

    # Now command a negative error; output should leave saturation within a few
    # steps (it would take far longer if the integrator had wound up).
    recovered = False
    for _ in range(5):
        out = pi.update(-10.0, dt)
        if out < 1.0:
            recovered = True
            break
    assert recovered


def test_reset():
    pi = PIController(kp=0.0, ki=1.0)
    pi.update(1.0, 1.0)
    assert pi.integrator != 0.0
    pi.reset()
    assert pi.integrator == 0.0

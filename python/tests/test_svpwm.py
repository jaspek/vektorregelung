import math

import pytest

from vektorregelung import duty_to_phase_voltages, max_voltage, space_vector_pwm
from vektorregelung.transforms import clarke


def test_zero_command_gives_50_percent_duty():
    da, db, dc = space_vector_pwm(0.0, 0.0, 540.0)
    assert (da, db, dc) == pytest.approx((0.5, 0.5, 0.5))


def test_duties_within_bounds():
    v_dc = 540.0
    v_max = max_voltage(v_dc)
    for k in range(360):
        theta = math.radians(k)
        v_alpha = v_max * math.cos(theta)
        v_beta = v_max * math.sin(theta)
        da, db, dc = space_vector_pwm(v_alpha, v_beta, v_dc)
        for d in (da, db, dc):
            assert 0.0 <= d <= 1.0


def test_command_reconstructed_from_duties():
    # Within the linear region, the duties must reproduce the commanded
    # (alpha, beta) voltage when mapped back via the inverter + Clarke.
    v_dc = 540.0
    v_max = max_voltage(v_dc) * 0.95  # stay safely inside the linear region
    for k in range(0, 360, 7):
        theta = math.radians(k)
        v_alpha = v_max * math.cos(theta)
        v_beta = v_max * math.sin(theta)
        da, db, dc = space_vector_pwm(v_alpha, v_beta, v_dc)
        va, vb, vc = duty_to_phase_voltages(da, db, dc, v_dc)
        a_rec, b_rec = clarke(va, vb, vc)
        assert a_rec == pytest.approx(v_alpha, abs=1e-6)
        assert b_rec == pytest.approx(v_beta, abs=1e-6)


def test_max_voltage_is_vdc_over_sqrt3():
    assert max_voltage(600.0) == pytest.approx(600.0 / math.sqrt(3.0))

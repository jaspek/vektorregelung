import math

import pytest

from vektorregelung import transforms as tf


def test_clarke_amplitude_invariant():
    # Balanced three-phase set, peak amplitude 10 at angle 0.
    a = 10.0
    b = 10.0 * math.cos(-2.0 * math.pi / 3.0)
    c = 10.0 * math.cos(2.0 * math.pi / 3.0)
    alpha, beta = tf.clarke(a, b, c)
    assert alpha == pytest.approx(10.0)
    assert beta == pytest.approx(0.0, abs=1e-12)


def test_clarke_balanced_matches_full_clarke():
    a, b = 3.0, -1.5
    c = -a - b
    full = tf.clarke(a, b, c)
    reduced = tf.clarke_balanced(a, b)
    assert reduced == pytest.approx(full)


def test_clarke_inverse_roundtrip_balanced():
    alpha, beta = 4.0, -2.0
    a, b, c = tf.inverse_clarke(alpha, beta)
    assert a + b + c == pytest.approx(0.0, abs=1e-12)
    alpha2, beta2 = tf.clarke(a, b, c)
    assert (alpha2, beta2) == pytest.approx((alpha, beta))


@pytest.mark.parametrize("theta", [0.0, 0.5, 1.234, math.pi, 5.0])
def test_park_inverse_roundtrip(theta):
    d, q = 7.0, -3.0
    alpha, beta = tf.inverse_park(d, q, theta)
    d2, q2 = tf.park(alpha, beta, theta)
    assert (d2, q2) == pytest.approx((d, q))


def test_park_preserves_magnitude():
    alpha, beta = 6.0, 8.0  # magnitude 10
    for theta in (0.0, 1.0, 2.5, 4.0):
        d, q = tf.park(alpha, beta, theta)
        assert math.hypot(d, q) == pytest.approx(10.0)


def test_park_aligns_vector_on_d_axis():
    # A vector at angle theta maps onto the +d axis (q == 0).
    theta = 0.9
    mag = 5.0
    alpha = mag * math.cos(theta)
    beta = mag * math.sin(theta)
    d, q = tf.park(alpha, beta, theta)
    assert d == pytest.approx(mag)
    assert q == pytest.approx(0.0, abs=1e-12)


def test_wrap_angle():
    two_pi = 2.0 * math.pi
    assert tf.wrap_angle(0.0) == pytest.approx(0.0)
    assert tf.wrap_angle(two_pi + 1.0) == pytest.approx(1.0)
    assert tf.wrap_angle(-1.0) == pytest.approx(two_pi - 1.0)
    assert 0.0 <= tf.wrap_angle(100.0) < two_pi

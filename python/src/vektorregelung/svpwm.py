"""Raumzeiger-Modulation (Space-Vector PWM).

Computes the three half-bridge duty cycles that synthesise a commanded stator
voltage vector (alpha/beta) from a DC-link voltage ``v_dc``.

This implementation uses the equivalent and very efficient *min/max common-mode
injection* method: the three sinusoidal phase references are shifted by a common
zero-sequence component so that the modulation index can reach the full
``v_dc / sqrt(3)`` of the linear SVPWM region. The resulting duty cycles are
identical to the classical sector/timing-based SVPWM but require no sector logic.
"""

from __future__ import annotations

from .transforms import inverse_clarke

# Maximum length of the voltage vector that SVPWM can synthesise without
# overmodulation, expressed as a fraction of v_dc:  |v| <= v_dc / sqrt(3).
import math

MAX_MODULATION = 1.0 / math.sqrt(3.0)


def max_voltage(v_dc: float) -> float:
    """Largest realisable stator voltage vector magnitude for a given v_dc."""
    return MAX_MODULATION * v_dc


def space_vector_pwm(
    v_alpha: float, v_beta: float, v_dc: float
) -> tuple[float, float, float]:
    """Return phase duty cycles ``(da, db, dc)`` in [0, 1].

    Parameters
    ----------
    v_alpha, v_beta:
        Commanded stator voltage vector in the stationary frame [V].
    v_dc:
        DC-link voltage [V].
    """
    if v_dc <= 0.0:
        return 0.5, 0.5, 0.5

    va, vb, vc = inverse_clarke(v_alpha, v_beta)

    v_max = max(va, vb, vc)
    v_min = min(va, vb, vc)
    # Zero-sequence (common-mode) injection centres the references in the
    # available DC-link window.
    v_offset = 0.5 * (v_max + v_min)

    inv_vdc = 1.0 / v_dc
    da = 0.5 + (va - v_offset) * inv_vdc
    db = 0.5 + (vb - v_offset) * inv_vdc
    dc = 0.5 + (vc - v_offset) * inv_vdc

    return _clamp01(da), _clamp01(db), _clamp01(dc)


def duty_to_phase_voltages(
    da: float, db: float, dc: float, v_dc: float
) -> tuple[float, float, float]:
    """Reconstruct the line-to-neutral voltages produced by given duty cycles.

    Useful for verification: feeding the SVPWM output back through this function
    must reproduce the commanded (alpha, beta) vector via :func:`clarke`.
    """
    # Leg (pole) voltages relative to the negative DC rail.
    leg_a = da * v_dc
    leg_b = db * v_dc
    leg_c = dc * v_dc
    common = (leg_a + leg_b + leg_c) / 3.0
    return leg_a - common, leg_b - common, leg_c - common


def _clamp01(x: float) -> float:
    if x < 0.0:
        return 0.0
    if x > 1.0:
        return 1.0
    return x

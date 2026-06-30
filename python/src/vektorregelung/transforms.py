"""Koordinatentransformationen für die Vektorregelung.

Clarke- und Park-Transformationen (und ihre Inversen) bilden das Fundament
jeder feldorientierten Regelung. Sie überführen die drei Strangsgrößen einer
Drehfeldmaschine in ein zweiachsiges, mit dem Fluss mitrotierendes
Koordinatensystem, in dem sich die Maschine wie eine Gleichstrommaschine regeln
lässt.

Coordinate transforms for field-oriented control. Pure standard library so the
module is portable and dependency free (the same math is mirrored in the C
library under ``c/``).

Konventionen / conventions
--------------------------
* Clarke transform is *amplitude invariant* (factor 2/3): a balanced three-phase
  set of peak amplitude ``A`` maps to an (alpha, beta) vector of magnitude ``A``.
* Angle ``theta`` is the electrical field angle in radians.
"""

from __future__ import annotations

import math

# sqrt(3)/2 and 1/sqrt(3) precomputed for clarity and speed.
_SQRT3_2 = math.sqrt(3.0) / 2.0
_ONE_OVER_SQRT3 = 1.0 / math.sqrt(3.0)


def clarke(a: float, b: float, c: float) -> tuple[float, float]:
    """Amplitude-invariant Clarke transform (abc -> alpha/beta).

    For a balanced system (``a + b + c == 0``) this reduces to ``alpha = a``.
    The general form is kept so unbalanced/zero-sequence inputs are handled
    correctly.
    """
    alpha = (2.0 / 3.0) * (a - 0.5 * b - 0.5 * c)
    beta = _ONE_OVER_SQRT3 * (b - c)
    return alpha, beta


def clarke_balanced(a: float, b: float) -> tuple[float, float]:
    """Clarke transform assuming a balanced system (``c = -a - b``).

    Only two current sensors are needed in many drives; this overload avoids
    measuring the third phase.
    """
    alpha = a
    beta = _ONE_OVER_SQRT3 * (a + 2.0 * b)
    return alpha, beta


def inverse_clarke(alpha: float, beta: float) -> tuple[float, float, float]:
    """Inverse Clarke transform (alpha/beta -> abc)."""
    a = alpha
    b = -0.5 * alpha + _SQRT3_2 * beta
    c = -0.5 * alpha - _SQRT3_2 * beta
    return a, b, c


def park(alpha: float, beta: float, theta: float) -> tuple[float, float]:
    """Park transform (alpha/beta -> d/q) for field angle ``theta``."""
    cos_t = math.cos(theta)
    sin_t = math.sin(theta)
    d = alpha * cos_t + beta * sin_t
    q = -alpha * sin_t + beta * cos_t
    return d, q


def inverse_park(d: float, q: float, theta: float) -> tuple[float, float]:
    """Inverse Park transform (d/q -> alpha/beta) for field angle ``theta``."""
    cos_t = math.cos(theta)
    sin_t = math.sin(theta)
    alpha = d * cos_t - q * sin_t
    beta = d * sin_t + q * cos_t
    return alpha, beta


def wrap_angle(theta: float) -> float:
    """Wrap an angle into the interval [0, 2*pi).

    Keeping the integrated field angle bounded avoids floating point precision
    loss over long runs (relevant for the integrator in indirect FOC).
    """
    two_pi = 2.0 * math.pi
    theta = math.fmod(theta, two_pi)
    if theta < 0.0:
        theta += two_pi
    return theta

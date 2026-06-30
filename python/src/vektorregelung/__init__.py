"""Vektorregelung -- indirect field-oriented control of an induction motor.

Public API:

* :class:`~vektorregelung.induction_motor.InductionMotorParameters`
* :class:`~vektorregelung.induction_motor.InductionMotor`
* :class:`~vektorregelung.pi_controller.PIController`
* :class:`~vektorregelung.foc.IndirectFOC`,
  :class:`~vektorregelung.foc.FocLimits`,
  :class:`~vektorregelung.foc.FocOutput`
* :func:`~vektorregelung.simulation.simulate`,
  :class:`~vektorregelung.simulation.SimulationLog`
* transform and SVPWM helpers
"""

from .foc import FocLimits, FocOutput, IndirectFOC
from .induction_motor import InductionMotor, InductionMotorParameters
from .pi_controller import PIController
from .simulation import SimulationLog, simulate
from .svpwm import (
    duty_to_phase_voltages,
    max_voltage,
    space_vector_pwm,
)
from .transforms import (
    clarke,
    clarke_balanced,
    inverse_clarke,
    inverse_park,
    park,
    wrap_angle,
)

__version__ = "0.1.0"

__all__ = [
    "InductionMotorParameters",
    "InductionMotor",
    "PIController",
    "IndirectFOC",
    "FocLimits",
    "FocOutput",
    "simulate",
    "SimulationLog",
    "clarke",
    "clarke_balanced",
    "inverse_clarke",
    "park",
    "inverse_park",
    "wrap_angle",
    "space_vector_pwm",
    "duty_to_phase_voltages",
    "max_voltage",
]

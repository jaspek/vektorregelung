"""Zeitdiskreter PI-Regler mit Anti-Windup.

Discrete-time PI controller with output saturation and conditional-integration
anti-windup. Used for the current loops (d and q) and the outer speed loop of
the field-oriented controller.

The integrator only accumulates while the output is unsaturated, or while
further integration would move the output *back* into the valid range. This is
the parameter-free "clamping"/conditional-integration anti-windup scheme, which
is robust and maps cleanly onto the C implementation.
"""

from __future__ import annotations

from dataclasses import dataclass, field

_INF = float("inf")


@dataclass
class PIController:
    """Parallel-form PI controller: ``u = kp*e + ki*integral(e)``.

    Parameters
    ----------
    kp, ki:
        Proportional and integral gains.
    out_min, out_max:
        Saturation limits applied to the controller output.
    """

    kp: float
    ki: float
    out_min: float = -_INF
    out_max: float = _INF
    # Accumulated integral term (already scaled, i.e. it stores ki * integral(e)).
    integrator: float = field(default=0.0)

    def reset(self, integrator: float = 0.0) -> None:
        """Reset the integrator state (e.g. when the drive is disabled)."""
        self.integrator = integrator

    def update(self, error: float, dt: float) -> float:
        """Advance the controller by one step.

        Parameters
        ----------
        error:
            Setpoint minus measurement.
        dt:
            Sample period in seconds.

        Returns
        -------
        float
            Saturated controller output.
        """
        proportional = self.kp * error
        integ_step = self.ki * error * dt
        tentative_integ = self.integrator + integ_step
        unsaturated = proportional + tentative_integ

        # Saturate the output.
        if unsaturated > self.out_max:
            output = self.out_max
        elif unsaturated < self.out_min:
            output = self.out_min
        else:
            output = unsaturated

        # Conditional integration (anti-windup): commit the integral step only
        # if the output is not saturated, or if this step pushes the output back
        # toward the linear region (i.e. error opposes the saturation excess).
        saturated = output != unsaturated
        excess = unsaturated - output  # >0 when clamped at the upper limit
        if not saturated or error * excess < 0.0:
            self.integrator = tentative_integ

        return output

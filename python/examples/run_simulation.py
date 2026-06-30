#!/usr/bin/env python3
"""Demo: indirect field-oriented control of a ~3 kW induction motor.

Runs a closed-loop simulation with three phases:

* 0.00-0.30 s: magnetisation -- the d-axis current establishes the rotor flux
  while the speed reference is held at zero.
* 0.30 s: speed reference steps to 150 rad/s (~1430 rpm).
* 1.00 s: a 10 N*m load torque is applied.

The script prints a summary of the steady-state tracking performance and, if
matplotlib is available, writes a plot to ``foc_simulation.png``.

Run from the ``python`` directory:

    PYTHONPATH=src python3 examples/run_simulation.py
"""

from __future__ import annotations

import os
import sys

# Allow running straight from the repo without installing the package.
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

from vektorregelung import (  # noqa: E402
    IndirectFOC,
    InductionMotor,
    InductionMotorParameters,
    simulate,
)

# --- Plant: typical 3 kW / 4-pole / 400 V induction machine -----------------
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

DT = 1.0e-4            # 10 kHz control / integration step
DURATION = 1.6        # seconds
V_DC = 600.0          # DC-link voltage [V]
ROTOR_FLUX_REF = 1.0  # rated rotor flux [Wb]
I_MAX = 20.0          # current limit [A]

SPEED_SETPOINT = 150.0   # rad/s
SPEED_STEP_TIME = 0.30   # s
LOAD_STEP_TIME = 1.00    # s
LOAD_TORQUE = 10.0       # N*m


def speed_reference(t: float) -> float:
    return SPEED_SETPOINT if t >= SPEED_STEP_TIME else 0.0


def load_profile(t: float) -> float:
    return LOAD_TORQUE if t >= LOAD_STEP_TIME else 0.0


def main() -> int:
    controller = IndirectFOC.from_motor(
        MOTOR,
        dt=DT,
        v_dc=V_DC,
        rotor_flux_ref=ROTOR_FLUX_REF,
        i_max=I_MAX,
    )
    motor = InductionMotor(MOTOR)

    log = simulate(
        controller=controller,
        motor=motor,
        speed_ref=speed_reference,
        load_torque=load_profile,
        duration=DURATION,
        dt=DT,
    )

    _print_summary(log)
    _maybe_plot(log)
    return 0


def _print_summary(log) -> None:
    final_speed = log.speed[-1]
    final_ref = log.speed_ref[-1]
    final_flux = log.flux_actual[-1]
    final_torque = log.torque[-1]

    # Settling: index where speed first reaches 95 % of the setpoint after the
    # speed step.
    target = SPEED_SETPOINT
    settle_t = None
    for t, w in zip(log.t, log.speed):
        if t >= SPEED_STEP_TIME and w >= 0.95 * target:
            settle_t = t - SPEED_STEP_TIME
            break

    print("=== Indirect FOC simulation summary ===")
    print(f"  steps simulated     : {len(log.t)}")
    print(f"  speed setpoint      : {final_ref:8.2f} rad/s")
    print(f"  final speed         : {final_speed:8.2f} rad/s")
    print(f"  speed error         : {final_ref - final_speed:8.4f} rad/s")
    if settle_t is not None:
        print(f"  95% settling time   : {settle_t * 1e3:8.1f} ms")
    print(f"  rotor flux (ref {ROTOR_FLUX_REF:.2f}) : {final_flux:8.4f} Wb")
    print(f"  final torque        : {final_torque:8.2f} N*m  "
          f"(load {LOAD_TORQUE:.1f} N*m)")


def _maybe_plot(log) -> None:
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("\n(matplotlib not installed -- skipping plot)")
        return

    fig, axes = plt.subplots(3, 1, figsize=(9, 9), sharex=True)

    axes[0].plot(log.t, log.speed_ref, "k--", label="reference")
    axes[0].plot(log.t, log.speed, "b", label="speed")
    axes[0].set_ylabel("speed [rad/s]")
    axes[0].legend(loc="lower right")
    axes[0].grid(True)

    axes[1].plot(log.t, log.torque, "g", label="Te")
    axes[1].plot(log.t, log.load_torque, "r--", label="load")
    axes[1].set_ylabel("torque [N*m]")
    axes[1].legend(loc="upper right")
    axes[1].grid(True)

    axes[2].plot(log.t, log.i_d, "m", label="id")
    axes[2].plot(log.t, log.i_q, "c", label="iq")
    axes[2].plot(log.t, log.flux_actual, "orange", label="rotor flux [Wb]")
    axes[2].set_ylabel("current [A] / flux [Wb]")
    axes[2].set_xlabel("time [s]")
    axes[2].legend(loc="upper right")
    axes[2].grid(True)

    fig.suptitle("Indirect Field-Oriented Control of an Induction Motor")
    fig.tight_layout()
    out_path = os.path.join(os.path.dirname(__file__), "foc_simulation.png")
    fig.savefig(out_path, dpi=110)
    print(f"\nPlot written to {out_path}")


if __name__ == "__main__":
    raise SystemExit(main())

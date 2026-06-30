/*
 * Closed-loop demo of the C Vektorregelung library: indirect FOC driving the
 * induction-motor model. Mirrors python/examples/run_simulation.py.
 *
 * Prints a steady-state summary and, when given a path argument, writes the
 * full trajectory as CSV (useful to diff against the Python reference).
 *
 * Build with the host (double-precision) configuration so results match the
 * Python implementation:  make demo   (see ../Makefile).
 */
#include <stdio.h>

#include "vektorregelung/foc.h"
#include "vektorregelung/induction_motor.h"

#define DT             ((vr_float)1.0e-4)
#define DURATION       ((vr_float)1.6)
#define V_DC           ((vr_float)600.0)
#define ROTOR_FLUX_REF ((vr_float)1.0)
#define I_MAX          ((vr_float)20.0)

#define SPEED_SETPOINT   ((vr_float)150.0)
#define SPEED_STEP_TIME  ((vr_float)0.30)
#define LOAD_STEP_TIME   ((vr_float)1.00)
#define LOAD_TORQUE      ((vr_float)10.0)

static vr_float speed_reference(vr_float t)
{
    return (t >= SPEED_STEP_TIME) ? SPEED_SETPOINT : (vr_float)0.0;
}

static vr_float load_profile(vr_float t)
{
    return (t >= LOAD_STEP_TIME) ? LOAD_TORQUE : (vr_float)0.0;
}

int main(int argc, char **argv)
{
    vr_motor_params_t params;
    vr_foc_t foc;
    vr_induction_motor_t motor;
    FILE *csv = NULL;
    long n_steps, k;
    vr_float final_speed = 0, final_flux = 0, final_torque = 0;
    vr_float settle_t = -1.0;

    /* 3 kW / 4-pole / 400 V induction machine. */
    params.rs = (vr_float)1.5;
    params.rr = (vr_float)1.2;
    params.ls = (vr_float)0.14;
    params.lr = (vr_float)0.14;
    params.lm = (vr_float)0.135;
    params.pole_pairs = 2;
    params.inertia = (vr_float)0.02;
    params.friction = (vr_float)0.001;

    vr_foc_init_tuned(&foc, &params, DT, V_DC, ROTOR_FLUX_REF, I_MAX,
                      (vr_float)0.0, (vr_float)0.0);
    vr_motor_init(&motor, &params);

    if (argc > 1) {
        csv = fopen(argv[1], "w");
        if (csv) {
            fprintf(csv, "t,speed_ref,speed,torque,load,id,iq,iq_ref,"
                         "flux_actual,flux_estimate\n");
        }
    }

    n_steps = (long)(DURATION / DT + 0.5);
    for (k = 0; k < n_steps; ++k) {
        vr_float t = (vr_float)k * DT;
        vr_float w_ref = speed_reference(t);
        vr_float t_load = load_profile(t);
        vr_abc_t i = vr_motor_phase_currents(&motor);
        vr_foc_output_t out = vr_foc_update(&foc, w_ref, i.a, i.b, i.c,
                                            motor.omega_m);

        vr_motor_step(&motor, out.v_ab.alpha, out.v_ab.beta, t_load, DT);

        if (settle_t < 0.0 && t >= SPEED_STEP_TIME
            && motor.omega_m >= (vr_float)0.95 * SPEED_SETPOINT) {
            settle_t = (double)(t - SPEED_STEP_TIME);
        }

        if (csv) {
            fprintf(csv,
                    "%.5f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f\n",
                    (double)t, (double)w_ref, (double)motor.omega_m,
                    (double)vr_motor_torque(&motor), (double)t_load,
                    (double)out.i_dq.d, (double)out.i_dq.q, (double)out.iq_ref,
                    (double)vr_motor_flux_magnitude(&motor),
                    (double)out.flux_estimate);
        }

        if (k == n_steps - 1) {
            final_speed = motor.omega_m;
            final_flux = vr_motor_flux_magnitude(&motor);
            final_torque = vr_motor_torque(&motor);
        }
    }

    if (csv) {
        fclose(csv);
    }

    printf("=== Indirect FOC simulation summary (C) ===\n");
    printf("  steps simulated     : %ld\n", n_steps);
    printf("  speed setpoint      : %8.2f rad/s\n", (double)SPEED_SETPOINT);
    printf("  final speed         : %8.2f rad/s\n", (double)final_speed);
    printf("  speed error         : %8.4f rad/s\n",
           (double)(SPEED_SETPOINT - final_speed));
    if (settle_t >= 0.0) {
        printf("  95%% settling time   : %8.1f ms\n", settle_t * 1e3);
    }
    printf("  rotor flux (ref %.2f) : %8.4f Wb\n",
           (double)ROTOR_FLUX_REF, (double)final_flux);
    printf("  final torque        : %8.2f N*m  (load %.1f N*m)\n",
           (double)final_torque, (double)LOAD_TORQUE);

    return 0;
}

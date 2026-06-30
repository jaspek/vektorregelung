/*
 * Induction-motor parameters, derived quantities, and a stationary-frame
 * dynamic model for host-side simulation.
 * Mirrors python/src/vektorregelung/induction_motor.py.
 *
 * The model (vr_induction_motor_t) is intended for simulation/testing on a
 * host; the parameter struct and derived helpers are also used by the
 * controller for tuning and feed-forward decoupling on the target.
 */
#ifndef VEKTORREGELUNG_INDUCTION_MOTOR_H
#define VEKTORREGELUNG_INDUCTION_MOTOR_H

#include "vektorregelung/types.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    vr_float rs;        /* stator resistance [Ohm] */
    vr_float rr;        /* rotor resistance [Ohm] */
    vr_float ls;        /* stator inductance Ls = Lls + Lm [H] */
    vr_float lr;        /* rotor inductance Lr = Llr + Lm [H] */
    vr_float lm;        /* mutual inductance Lm [H] */
    int pole_pairs;     /* number of pole pairs */
    vr_float inertia;   /* rotor inertia J [kg*m^2] */
    vr_float friction;  /* viscous friction B [N*m*s] */
} vr_motor_params_t;

/* Derived parameters. */
vr_float vr_motor_sigma(const vr_motor_params_t *p);     /* leakage factor */
vr_float vr_motor_sigma_ls(const vr_motor_params_t *p);  /* transient inductance */
vr_float vr_motor_tr(const vr_motor_params_t *p);        /* rotor time constant */
vr_float vr_motor_r_sigma(const vr_motor_params_t *p);   /* effective resistance */
vr_float vr_motor_torque_constant(const vr_motor_params_t *p, vr_float flux);

/* Stationary-frame simulation model. State: stator current (alpha/beta),
 * rotor flux (alpha/beta), mechanical speed and angle. */
typedef struct {
    vr_motor_params_t params;
    vr_float i_alpha;
    vr_float i_beta;
    vr_float psi_r_alpha;
    vr_float psi_r_beta;
    vr_float omega_m;
    vr_float theta_m;
} vr_induction_motor_t;

void vr_motor_init(vr_induction_motor_t *m, const vr_motor_params_t *params);
vr_abc_t vr_motor_phase_currents(const vr_induction_motor_t *m);
vr_float vr_motor_torque(const vr_induction_motor_t *m);
vr_float vr_motor_flux_magnitude(const vr_induction_motor_t *m);

/* Advance the model by dt seconds (classical RK4). */
void vr_motor_step(vr_induction_motor_t *m, vr_float v_alpha, vr_float v_beta,
                   vr_float load_torque, vr_float dt);

#ifdef __cplusplus
}
#endif

#endif /* VEKTORREGELUNG_INDUCTION_MOTOR_H */

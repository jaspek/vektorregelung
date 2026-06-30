/*
 * Indirect field-oriented controller (IFOC) for an induction motor.
 * Mirrors python/src/vektorregelung/foc.py.
 *
 * One call to vr_foc_update() performs the full chain: Clarke/Park of the
 * measured currents, current-model flux estimation and slip computation,
 * outer speed PI loop, inner d/q current PI loops with decoupling
 * feed-forward, voltage limiting, inverse Park and space-vector PWM.
 */
#ifndef VEKTORREGELUNG_FOC_H
#define VEKTORREGELUNG_FOC_H

#include "vektorregelung/types.h"
#include "vektorregelung/pi_controller.h"
#include "vektorregelung/induction_motor.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    vr_float v_dc;     /* DC-link voltage [V] */
    vr_float i_max;    /* stator current magnitude limit [A] */
    vr_float flux_min; /* floor for the flux estimate [Wb] */
} vr_foc_limits_t;

typedef struct {
    vr_abc_t duty;       /* half-bridge duty cycles in [0, 1] */
    vr_ab_t v_ab;        /* commanded stator voltage (alpha/beta) [V] */
    vr_dq_t v_dq;        /* commanded stator voltage (d/q) [V] */
    vr_dq_t i_dq;        /* measured stator current (d/q) [A] */
    vr_float id_ref;     /* magnetising-current reference [A] */
    vr_float iq_ref;     /* torque-current reference [A] */
    vr_float theta_e;    /* estimated field angle [rad] */
    vr_float omega_slip; /* slip frequency [rad/s] */
    vr_float flux_estimate; /* estimated rotor flux [Wb] */
} vr_foc_output_t;

typedef struct {
    vr_motor_params_t motor;
    vr_float dt;
    vr_foc_limits_t limits;
    vr_float rotor_flux_ref;
    vr_float id_ref;            /* = rotor_flux_ref / Lm */

    vr_pi_t id_pi;
    vr_pi_t iq_pi;
    vr_pi_t speed_pi;

    vr_float theta_e;          /* integrated field angle */
    vr_float flux_estimate;    /* current-model rotor flux estimate */
} vr_foc_t;

/* Initialise with explicit, pre-built PI controllers. */
void vr_foc_init(vr_foc_t *foc, const vr_motor_params_t *motor, vr_float dt,
                 const vr_foc_limits_t *limits, vr_float rotor_flux_ref,
                 vr_pi_t id_pi, vr_pi_t iq_pi, vr_pi_t speed_pi);

/* Initialise with PI gains automatically tuned from the motor parameters.
 * Pass current_bw_hz <= 0 / speed_bw_hz <= 0 to use the defaults (200 / 15). */
void vr_foc_init_tuned(vr_foc_t *foc, const vr_motor_params_t *motor, vr_float dt,
                       vr_float v_dc, vr_float rotor_flux_ref, vr_float i_max,
                       vr_float current_bw_hz, vr_float speed_bw_hz);

/* Reset estimators and PI integrators (e.g. on drive enable). */
void vr_foc_reset(vr_foc_t *foc);

/* Run one control step. */
vr_foc_output_t vr_foc_update(vr_foc_t *foc, vr_float speed_ref,
                              vr_float ia, vr_float ib, vr_float ic,
                              vr_float omega_m);

#ifdef __cplusplus
}
#endif

#endif /* VEKTORREGELUNG_FOC_H */

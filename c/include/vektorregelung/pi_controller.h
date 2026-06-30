/*
 * Discrete-time PI controller with conditional-integration anti-windup.
 * Mirrors python/src/vektorregelung/pi_controller.py.
 */
#ifndef VEKTORREGELUNG_PI_CONTROLLER_H
#define VEKTORREGELUNG_PI_CONTROLLER_H

#include "vektorregelung/types.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    vr_float kp;
    vr_float ki;
    vr_float out_min;
    vr_float out_max;
    vr_float integrator; /* stores ki * integral(error) */
} vr_pi_t;

/* Initialise gains and limits; clears the integrator. */
void vr_pi_init(vr_pi_t *pi, vr_float kp, vr_float ki,
                vr_float out_min, vr_float out_max);

/* Reset the integrator to a chosen value (commonly 0). */
void vr_pi_reset(vr_pi_t *pi, vr_float integrator);

/* Advance one step; returns the saturated output. */
vr_float vr_pi_update(vr_pi_t *pi, vr_float error, vr_float dt);

#ifdef __cplusplus
}
#endif

#endif /* VEKTORREGELUNG_PI_CONTROLLER_H */

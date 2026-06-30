#include "vektorregelung/pi_controller.h"

void vr_pi_init(vr_pi_t *pi, vr_float kp, vr_float ki,
                vr_float out_min, vr_float out_max)
{
    pi->kp = kp;
    pi->ki = ki;
    pi->out_min = out_min;
    pi->out_max = out_max;
    pi->integrator = (vr_float)0.0;
}

void vr_pi_reset(vr_pi_t *pi, vr_float integrator)
{
    pi->integrator = integrator;
}

vr_float vr_pi_update(vr_pi_t *pi, vr_float error, vr_float dt)
{
    vr_float proportional = pi->kp * error;
    vr_float integ_step = pi->ki * error * dt;
    vr_float tentative_integ = pi->integrator + integ_step;
    vr_float unsaturated = proportional + tentative_integ;
    vr_float output;
    vr_float excess;
    int saturated;

    if (unsaturated > pi->out_max) {
        output = pi->out_max;
    } else if (unsaturated < pi->out_min) {
        output = pi->out_min;
    } else {
        output = unsaturated;
    }

    /* Conditional integration anti-windup: only commit the integral step when
     * unsaturated, or when it pushes the output back toward the linear range. */
    saturated = (output != unsaturated);
    excess = unsaturated - output; /* >0 when clamped at the upper limit */
    if (!saturated || (error * excess < (vr_float)0.0)) {
        pi->integrator = tentative_integ;
    }

    return output;
}

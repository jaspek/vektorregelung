#include "vektorregelung/foc.h"
#include "vektorregelung/transforms.h"
#include "vektorregelung/svpwm.h"

static vr_float limit_iq(const vr_foc_t *foc, vr_float id_ref, vr_float iq_ref)
{
    /* Magnetising current has priority; torque current uses the remaining
     * headroom so the total current vector stays within i_max. */
    vr_float i_max = foc->limits.i_max;
    vr_float head = i_max * i_max - id_ref * id_ref;
    vr_float iq_limit;

    if (head <= (vr_float)0.0) {
        return (vr_float)0.0;
    }
    iq_limit = VR_SQRT(head);
    if (iq_ref > iq_limit) {
        return iq_limit;
    }
    if (iq_ref < -iq_limit) {
        return -iq_limit;
    }
    return iq_ref;
}

static void limit_voltage(const vr_foc_t *foc, vr_float *v_d, vr_float *v_q)
{
    vr_float v_max = vr_max_voltage(foc->limits.v_dc);
    vr_float mag = vr_magnitude(*v_d, *v_q);
    if (mag > v_max && mag > (vr_float)0.0) {
        vr_float scale = v_max / mag;
        *v_d *= scale;
        *v_q *= scale;
    }
}

void vr_foc_init(vr_foc_t *foc, const vr_motor_params_t *motor, vr_float dt,
                 const vr_foc_limits_t *limits, vr_float rotor_flux_ref,
                 vr_pi_t id_pi, vr_pi_t iq_pi, vr_pi_t speed_pi)
{
    foc->motor = *motor;
    foc->dt = dt;
    foc->limits = *limits;
    foc->rotor_flux_ref = rotor_flux_ref;
    foc->id_ref = rotor_flux_ref / motor->lm;
    foc->id_pi = id_pi;
    foc->iq_pi = iq_pi;
    foc->speed_pi = speed_pi;
    foc->theta_e = (vr_float)0.0;
    foc->flux_estimate = (vr_float)0.0;
}

void vr_foc_init_tuned(vr_foc_t *foc, const vr_motor_params_t *motor, vr_float dt,
                       vr_float v_dc, vr_float rotor_flux_ref, vr_float i_max,
                       vr_float current_bw_hz, vr_float speed_bw_hz)
{
    vr_foc_limits_t limits;
    vr_pi_t id_pi, iq_pi, speed_pi;
    vr_float w_ci, kp_i, ki_i, v_lim;
    vr_float w_cs, kt, kp_s, ki_s;

    if (current_bw_hz <= (vr_float)0.0) current_bw_hz = (vr_float)200.0;
    if (speed_bw_hz <= (vr_float)0.0) speed_bw_hz = (vr_float)15.0;

    /* Current loops: plant 1/(sigmaLs*s + Rsigma), internal-model tuning. */
    w_ci = (vr_float)2.0 * VR_PI * current_bw_hz;
    kp_i = vr_motor_sigma_ls(motor) * w_ci;
    ki_i = vr_motor_r_sigma(motor) * w_ci;
    v_lim = vr_max_voltage(v_dc);
    vr_pi_init(&id_pi, kp_i, ki_i, -v_lim, v_lim);
    vr_pi_init(&iq_pi, kp_i, ki_i, -v_lim, v_lim);

    /* Speed loop: Te = Kt*iq, plant Kt/(J*s); zero a decade below crossover. */
    w_cs = (vr_float)2.0 * VR_PI * speed_bw_hz;
    kt = vr_motor_torque_constant(motor, rotor_flux_ref);
    kp_s = motor->inertia * w_cs / kt;
    ki_s = kp_s * w_cs / (vr_float)10.0;
    vr_pi_init(&speed_pi, kp_s, ki_s, -i_max, i_max);

    limits.v_dc = v_dc;
    limits.i_max = i_max;
    limits.flux_min = (vr_float)1e-3;

    vr_foc_init(foc, motor, dt, &limits, rotor_flux_ref, id_pi, iq_pi, speed_pi);
}

void vr_foc_reset(vr_foc_t *foc)
{
    foc->theta_e = (vr_float)0.0;
    foc->flux_estimate = (vr_float)0.0;
    vr_pi_reset(&foc->id_pi, (vr_float)0.0);
    vr_pi_reset(&foc->iq_pi, (vr_float)0.0);
    vr_pi_reset(&foc->speed_pi, (vr_float)0.0);
}

vr_foc_output_t vr_foc_update(vr_foc_t *foc, vr_float speed_ref,
                              vr_float ia, vr_float ib, vr_float ic,
                              vr_float omega_m)
{
    const vr_motor_params_t *motor = &foc->motor;
    vr_float dt = foc->dt;
    vr_ab_t i_ab;
    vr_dq_t i_dq;
    vr_float flux, omega_slip, omega_e;
    vr_float iq_ref, id_ref, v_d, v_q;
    vr_float sigma_ls, lm_over_lr, tr;
    vr_ab_t v_ab;
    vr_abc_t duty;
    vr_foc_output_t out;

    tr = vr_motor_tr(motor);

    /* 1) Measure: abc -> alpha/beta -> d/q (rotor-flux frame). */
    i_ab = vr_clarke(ia, ib, ic);
    i_dq = vr_park(i_ab.alpha, i_ab.beta, foc->theta_e);

    /* 2) Rotor-flux estimator (current model): dpsi/dt = (Lm*id - psi)/Tr. */
    foc->flux_estimate += dt * (motor->lm * i_dq.d - foc->flux_estimate) / tr;
    flux = foc->flux_estimate;
    if (flux < foc->limits.flux_min) {
        flux = foc->limits.flux_min;
    }

    /* 3) Slip relation and field-angle integration. */
    omega_slip = motor->lm * i_dq.q / (tr * flux);
    omega_e = (vr_float)motor->pole_pairs * omega_m + omega_slip;
    foc->theta_e = vr_wrap_angle(foc->theta_e + dt * omega_e);

    /* 4) Outer speed loop -> torque-current reference. */
    iq_ref = vr_pi_update(&foc->speed_pi, speed_ref - omega_m, dt);
    id_ref = foc->id_ref;
    iq_ref = limit_iq(foc, id_ref, iq_ref);

    /* 5) Inner current loops. */
    v_d = vr_pi_update(&foc->id_pi, id_ref - i_dq.d, dt);
    v_q = vr_pi_update(&foc->iq_pi, iq_ref - i_dq.q, dt);

    /* 6) Cross-coupling / back-EMF feed-forward decoupling. */
    sigma_ls = vr_motor_sigma_ls(motor);
    lm_over_lr = motor->lm / motor->lr;
    v_d += -omega_e * sigma_ls * i_dq.q;
    v_q += omega_e * (sigma_ls * i_dq.d + lm_over_lr * flux);

    /* 7) Limit the voltage vector to the SVPWM linear region. */
    limit_voltage(foc, &v_d, &v_q);

    /* 8) Inverse Park + space-vector PWM. */
    v_ab = vr_inverse_park(v_d, v_q, foc->theta_e);
    duty = vr_svpwm(v_ab.alpha, v_ab.beta, foc->limits.v_dc);

    out.duty = duty;
    out.v_ab = v_ab;
    out.v_dq.d = v_d;
    out.v_dq.q = v_q;
    out.i_dq = i_dq;
    out.id_ref = id_ref;
    out.iq_ref = iq_ref;
    out.theta_e = foc->theta_e;
    out.omega_slip = omega_slip;
    out.flux_estimate = foc->flux_estimate;
    return out;
}

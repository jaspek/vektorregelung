#include "vektorregelung/induction_motor.h"
#include "vektorregelung/transforms.h"

/* State vector layout for the RK4 integrator. */
enum {
    S_IA = 0,   /* stator current alpha */
    S_IB,       /* stator current beta */
    S_PA,       /* rotor flux alpha */
    S_PB,       /* rotor flux beta */
    S_W,        /* mechanical speed */
    S_TH,       /* mechanical angle */
    S_N
};

vr_float vr_motor_sigma(const vr_motor_params_t *p)
{
    return (vr_float)1.0 - (p->lm * p->lm) / (p->ls * p->lr);
}

vr_float vr_motor_sigma_ls(const vr_motor_params_t *p)
{
    return vr_motor_sigma(p) * p->ls;
}

vr_float vr_motor_tr(const vr_motor_params_t *p)
{
    return p->lr / p->rr;
}

vr_float vr_motor_r_sigma(const vr_motor_params_t *p)
{
    vr_float ratio = p->lm / p->lr;
    return p->rs + p->rr * ratio * ratio;
}

vr_float vr_motor_torque_constant(const vr_motor_params_t *p, vr_float flux)
{
    return (vr_float)1.5 * (vr_float)p->pole_pairs * (p->lm / p->lr) * flux;
}

void vr_motor_init(vr_induction_motor_t *m, const vr_motor_params_t *params)
{
    m->params = *params;
    m->i_alpha = (vr_float)0.0;
    m->i_beta = (vr_float)0.0;
    m->psi_r_alpha = (vr_float)0.0;
    m->psi_r_beta = (vr_float)0.0;
    m->omega_m = (vr_float)0.0;
    m->theta_m = (vr_float)0.0;
}

vr_abc_t vr_motor_phase_currents(const vr_induction_motor_t *m)
{
    return vr_inverse_clarke(m->i_alpha, m->i_beta);
}

vr_float vr_motor_torque(const vr_induction_motor_t *m)
{
    const vr_motor_params_t *p = &m->params;
    return (vr_float)1.5 * (vr_float)p->pole_pairs * (p->lm / p->lr)
           * (m->psi_r_alpha * m->i_beta - m->psi_r_beta * m->i_alpha);
}

vr_float vr_motor_flux_magnitude(const vr_induction_motor_t *m)
{
    return vr_magnitude(m->psi_r_alpha, m->psi_r_beta);
}

/* Compute the state derivative for a given (sub-)state. */
static void derivatives(const vr_motor_params_t *p, const vr_float s[S_N],
                        vr_float v_alpha, vr_float v_beta, vr_float load_torque,
                        vr_float ds[S_N])
{
    vr_float sigma_ls = vr_motor_sigma_ls(p);
    vr_float tr_inv = p->rr / p->lr;          /* 1 / Tr */
    vr_float lm_over_lr = p->lm / p->lr;
    vr_float omega_r = (vr_float)p->pole_pairs * s[S_W]; /* electrical speed */
    vr_float r_eff = p->rs + p->rr * lm_over_lr * lm_over_lr;
    vr_float torque;

    ds[S_IA] = (-r_eff * s[S_IA]
                + lm_over_lr * tr_inv * s[S_PA]
                + lm_over_lr * omega_r * s[S_PB]
                + v_alpha) / sigma_ls;
    ds[S_IB] = (-r_eff * s[S_IB]
                + lm_over_lr * tr_inv * s[S_PB]
                - lm_over_lr * omega_r * s[S_PA]
                + v_beta) / sigma_ls;

    ds[S_PA] = p->lm * tr_inv * s[S_IA] - tr_inv * s[S_PA] - omega_r * s[S_PB];
    ds[S_PB] = p->lm * tr_inv * s[S_IB] - tr_inv * s[S_PB] + omega_r * s[S_PA];

    torque = (vr_float)1.5 * (vr_float)p->pole_pairs * lm_over_lr
             * (s[S_PA] * s[S_IB] - s[S_PB] * s[S_IA]);
    ds[S_W] = (torque - load_torque - p->friction * s[S_W]) / p->inertia;
    ds[S_TH] = s[S_W];
}

void vr_motor_step(vr_induction_motor_t *m, vr_float v_alpha, vr_float v_beta,
                   vr_float load_torque, vr_float dt)
{
    const vr_motor_params_t *p = &m->params;
    vr_float s0[S_N], stmp[S_N];
    vr_float k1[S_N], k2[S_N], k3[S_N], k4[S_N];
    int i;

    s0[S_IA] = m->i_alpha;
    s0[S_IB] = m->i_beta;
    s0[S_PA] = m->psi_r_alpha;
    s0[S_PB] = m->psi_r_beta;
    s0[S_W] = m->omega_m;
    s0[S_TH] = m->theta_m;

    derivatives(p, s0, v_alpha, v_beta, load_torque, k1);
    for (i = 0; i < S_N; ++i) stmp[i] = s0[i] + (vr_float)0.5 * dt * k1[i];
    derivatives(p, stmp, v_alpha, v_beta, load_torque, k2);
    for (i = 0; i < S_N; ++i) stmp[i] = s0[i] + (vr_float)0.5 * dt * k2[i];
    derivatives(p, stmp, v_alpha, v_beta, load_torque, k3);
    for (i = 0; i < S_N; ++i) stmp[i] = s0[i] + dt * k3[i];
    derivatives(p, stmp, v_alpha, v_beta, load_torque, k4);

    for (i = 0; i < S_N; ++i) {
        s0[i] += (dt / (vr_float)6.0) * (k1[i] + (vr_float)2.0 * k2[i]
                                         + (vr_float)2.0 * k3[i] + k4[i]);
    }

    m->i_alpha = s0[S_IA];
    m->i_beta = s0[S_IB];
    m->psi_r_alpha = s0[S_PA];
    m->psi_r_beta = s0[S_PB];
    m->omega_m = s0[S_W];
    m->theta_m = s0[S_TH];
}

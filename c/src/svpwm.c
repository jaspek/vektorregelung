#include "vektorregelung/svpwm.h"
#include "vektorregelung/transforms.h"

#define VR_ONE_OVER_SQRT3 ((vr_float)0.57735026918962576451)

static vr_float vr_clamp01(vr_float x)
{
    if (x < (vr_float)0.0) {
        return (vr_float)0.0;
    }
    if (x > (vr_float)1.0) {
        return (vr_float)1.0;
    }
    return x;
}

vr_float vr_max_voltage(vr_float v_dc)
{
    return VR_ONE_OVER_SQRT3 * v_dc;
}

vr_abc_t vr_svpwm(vr_float v_alpha, vr_float v_beta, vr_float v_dc)
{
    vr_abc_t duty;
    vr_abc_t v;
    vr_float v_max, v_min, v_offset, inv_vdc;

    if (v_dc <= (vr_float)0.0) {
        duty.a = duty.b = duty.c = (vr_float)0.5;
        return duty;
    }

    v = vr_inverse_clarke(v_alpha, v_beta);

    v_max = v.a;
    if (v.b > v_max) v_max = v.b;
    if (v.c > v_max) v_max = v.c;

    v_min = v.a;
    if (v.b < v_min) v_min = v.b;
    if (v.c < v_min) v_min = v.c;

    /* Zero-sequence injection centres the references in the DC-link window. */
    v_offset = (vr_float)0.5 * (v_max + v_min);
    inv_vdc = (vr_float)1.0 / v_dc;

    duty.a = vr_clamp01((vr_float)0.5 + (v.a - v_offset) * inv_vdc);
    duty.b = vr_clamp01((vr_float)0.5 + (v.b - v_offset) * inv_vdc);
    duty.c = vr_clamp01((vr_float)0.5 + (v.c - v_offset) * inv_vdc);

    return duty;
}

vr_abc_t vr_duty_to_phase_voltages(vr_float da, vr_float db, vr_float dc,
                                   vr_float v_dc)
{
    vr_abc_t out;
    vr_float leg_a = da * v_dc;
    vr_float leg_b = db * v_dc;
    vr_float leg_c = dc * v_dc;
    vr_float common = (leg_a + leg_b + leg_c) / (vr_float)3.0;
    out.a = leg_a - common;
    out.b = leg_b - common;
    out.c = leg_c - common;
    return out;
}

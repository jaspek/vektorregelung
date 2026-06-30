#include "vektorregelung/transforms.h"

/* sqrt(3)/2 and 1/sqrt(3). */
#define VR_SQRT3_2       ((vr_float)0.86602540378443864676)
#define VR_ONE_OVER_SQRT3 ((vr_float)0.57735026918962576451)
#define VR_TWO_THIRDS    ((vr_float)0.66666666666666666667)

vr_ab_t vr_clarke(vr_float a, vr_float b, vr_float c)
{
    vr_ab_t out;
    out.alpha = VR_TWO_THIRDS * (a - (vr_float)0.5 * b - (vr_float)0.5 * c);
    out.beta = VR_ONE_OVER_SQRT3 * (b - c);
    return out;
}

vr_ab_t vr_clarke_balanced(vr_float a, vr_float b)
{
    vr_ab_t out;
    out.alpha = a;
    out.beta = VR_ONE_OVER_SQRT3 * (a + (vr_float)2.0 * b);
    return out;
}

vr_abc_t vr_inverse_clarke(vr_float alpha, vr_float beta)
{
    vr_abc_t out;
    out.a = alpha;
    out.b = -(vr_float)0.5 * alpha + VR_SQRT3_2 * beta;
    out.c = -(vr_float)0.5 * alpha - VR_SQRT3_2 * beta;
    return out;
}

vr_dq_t vr_park(vr_float alpha, vr_float beta, vr_float theta)
{
    vr_float cos_t = VR_COS(theta);
    vr_float sin_t = VR_SIN(theta);
    vr_dq_t out;
    out.d = alpha * cos_t + beta * sin_t;
    out.q = -alpha * sin_t + beta * cos_t;
    return out;
}

vr_ab_t vr_inverse_park(vr_float d, vr_float q, vr_float theta)
{
    vr_float cos_t = VR_COS(theta);
    vr_float sin_t = VR_SIN(theta);
    vr_ab_t out;
    out.alpha = d * cos_t - q * sin_t;
    out.beta = d * sin_t + q * cos_t;
    return out;
}

vr_float vr_wrap_angle(vr_float theta)
{
    vr_float two_pi = (vr_float)2.0 * VR_PI;
    theta = VR_FMOD(theta, two_pi);
    if (theta < (vr_float)0.0) {
        theta += two_pi;
    }
    return theta;
}

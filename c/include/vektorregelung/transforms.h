/*
 * Clarke/Park coordinate transforms (Koordinatentransformationen).
 *
 * Mirrors python/src/vektorregelung/transforms.py. The Clarke transform is
 * amplitude invariant (factor 2/3).
 */
#ifndef VEKTORREGELUNG_TRANSFORMS_H
#define VEKTORREGELUNG_TRANSFORMS_H

#include "vektorregelung/types.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Amplitude-invariant Clarke transform (abc -> alpha/beta). */
vr_ab_t vr_clarke(vr_float a, vr_float b, vr_float c);

/* Clarke transform assuming a balanced system (c = -a - b). */
vr_ab_t vr_clarke_balanced(vr_float a, vr_float b);

/* Inverse Clarke transform (alpha/beta -> abc). */
vr_abc_t vr_inverse_clarke(vr_float alpha, vr_float beta);

/* Park transform (alpha/beta -> d/q) for field angle theta. */
vr_dq_t vr_park(vr_float alpha, vr_float beta, vr_float theta);

/* Inverse Park transform (d/q -> alpha/beta) for field angle theta. */
vr_ab_t vr_inverse_park(vr_float d, vr_float q, vr_float theta);

/* Wrap an angle into [0, 2*pi). */
vr_float vr_wrap_angle(vr_float theta);

#ifdef __cplusplus
}
#endif

#endif /* VEKTORREGELUNG_TRANSFORMS_H */

/*
 * Common types and math configuration for the Vektorregelung C library.
 *
 * The library is written around a single configurable real type `vr_float`.
 * By default it is `float`, which matches the single-precision FPU of most
 * motor-control microcontrollers. Define VR_REAL_DOUBLE to use `double`
 * (e.g. for host-side simulation that must match the Python reference).
 *
 * The library performs no dynamic memory allocation and depends only on the
 * standard C math library, so it is suitable for bare-metal embedded targets.
 */
#ifndef VEKTORREGELUNG_TYPES_H
#define VEKTORREGELUNG_TYPES_H

#include <math.h>

#ifdef VR_REAL_DOUBLE
typedef double vr_float;
#define VR_COS  cos
#define VR_SIN  sin
#define VR_SQRT sqrt
#define VR_FMOD fmod
#define VR_FABS fabs
#define VR_PI   3.14159265358979323846
#else
typedef float vr_float;
#define VR_COS  cosf
#define VR_SIN  sinf
#define VR_SQRT sqrtf
#define VR_FMOD fmodf
#define VR_FABS fabsf
#define VR_PI   3.14159265358979323846f
#endif

/* Stationary-frame (alpha/beta) space vector. */
typedef struct {
    vr_float alpha;
    vr_float beta;
} vr_ab_t;

/* Rotating-frame (d/q) space vector. */
typedef struct {
    vr_float d;
    vr_float q;
} vr_dq_t;

/* Three-phase (a/b/c) quantity. */
typedef struct {
    vr_float a;
    vr_float b;
    vr_float c;
} vr_abc_t;

/* Magnitude of a 2-D vector without pulling in hypot() (embedded friendly). */
static inline vr_float vr_magnitude(vr_float x, vr_float y)
{
    return VR_SQRT(x * x + y * y);
}

#endif /* VEKTORREGELUNG_TYPES_H */

/*
 * Space-vector PWM (Raumzeiger-Modulation) via min/max common-mode injection.
 * Mirrors python/src/vektorregelung/svpwm.py.
 */
#ifndef VEKTORREGELUNG_SVPWM_H
#define VEKTORREGELUNG_SVPWM_H

#include "vektorregelung/types.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Largest realisable voltage-vector magnitude: v_dc / sqrt(3). */
vr_float vr_max_voltage(vr_float v_dc);

/* Phase duty cycles in [0, 1] for a commanded (alpha, beta) voltage. */
vr_abc_t vr_svpwm(vr_float v_alpha, vr_float v_beta, vr_float v_dc);

/* Reconstruct line-to-neutral voltages produced by given duty cycles. */
vr_abc_t vr_duty_to_phase_voltages(vr_float da, vr_float db, vr_float dc,
                                   vr_float v_dc);

#ifdef __cplusplus
}
#endif

#endif /* VEKTORREGELUNG_SVPWM_H */

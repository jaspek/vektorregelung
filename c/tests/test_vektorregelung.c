/*
 * Lightweight self-contained test suite for the C library (no framework).
 * Exits non-zero if any check fails.
 */
#include <math.h>
#include <stdio.h>

#include "vektorregelung/foc.h"
#include "vektorregelung/induction_motor.h"
#include "vektorregelung/pi_controller.h"
#include "vektorregelung/svpwm.h"
#include "vektorregelung/transforms.h"

/* M_PI is not part of strict ISO C99; define a local constant instead. */
#define TEST_PI 3.14159265358979323846

static int g_failures = 0;
static int g_checks = 0;

#define CHECK(cond)                                                          \
    do {                                                                     \
        g_checks++;                                                          \
        if (!(cond)) {                                                       \
            g_failures++;                                                    \
            printf("  FAIL %s:%d  %s\n", __FILE__, __LINE__, #cond);         \
        }                                                                    \
    } while (0)

#define CLOSE(a, b, tol) (fabs((double)(a) - (double)(b)) <= (double)(tol))

static void test_transforms(void)
{
    /* Amplitude-invariant Clarke on a balanced set. */
    vr_float a = 10.0f;
    vr_float b = 10.0f * (vr_float)cos(-2.0 * TEST_PI / 3.0);
    vr_float c = 10.0f * (vr_float)cos(2.0 * TEST_PI / 3.0);
    vr_ab_t ab = vr_clarke(a, b, c);
    CHECK(CLOSE(ab.alpha, 10.0, 1e-4));
    CHECK(CLOSE(ab.beta, 0.0, 1e-4));

    /* Park/inverse-Park round trip preserves the vector. */
    {
        vr_float theta = 1.234f;
        vr_ab_t r = vr_inverse_park(7.0f, -3.0f, theta);
        vr_dq_t d = vr_park(r.alpha, r.beta, theta);
        CHECK(CLOSE(d.d, 7.0, 1e-4));
        CHECK(CLOSE(d.q, -3.0, 1e-4));
    }

    /* Park preserves magnitude. */
    {
        vr_dq_t d = vr_park(6.0f, 8.0f, 1.0f);
        CHECK(CLOSE(vr_magnitude(d.d, d.q), 10.0, 1e-4));
    }

    /* wrap_angle. */
    CHECK(CLOSE(vr_wrap_angle((vr_float)(2.0 * TEST_PI + 1.0)), 1.0, 1e-4));
    CHECK(vr_wrap_angle(-1.0f) > 0.0f);
}

static void test_pi(void)
{
    /* Proportional only. */
    vr_pi_t pi;
    vr_pi_init(&pi, 2.0f, 0.0f, -1e30f, 1e30f);
    CHECK(CLOSE(vr_pi_update(&pi, 3.0f, 0.001f), 6.0, 1e-5));

    /* Integral accumulates to remove steady error. */
    vr_pi_init(&pi, 0.0f, 10.0f, -1e30f, 1e30f);
    {
        int i;
        vr_float out = 0.0f;
        for (i = 0; i < 100; ++i) out = vr_pi_update(&pi, 1.0f, 0.01f);
        CHECK(CLOSE(out, 10.0, 1e-3));
    }

    /* Saturation. */
    vr_pi_init(&pi, 5.0f, 0.0f, -1.0f, 1.0f);
    CHECK(CLOSE(vr_pi_update(&pi, 100.0f, 0.001f), 1.0, 1e-6));

    /* Anti-windup: recovers quickly after reversing the error. */
    vr_pi_init(&pi, 1.0f, 1000.0f, -1.0f, 1.0f);
    {
        int i, recovered = 0;
        for (i = 0; i < 1000; ++i) vr_pi_update(&pi, 10.0f, 0.001f);
        for (i = 0; i < 5; ++i) {
            if (vr_pi_update(&pi, -10.0f, 0.001f) < 1.0f) { recovered = 1; break; }
        }
        CHECK(recovered);
    }
}

static void test_svpwm(void)
{
    vr_float v_dc = 540.0f;
    vr_float v_max = vr_max_voltage(v_dc) * 0.95f;
    int k;

    CHECK(CLOSE(vr_max_voltage(600.0f), 600.0 / sqrt(3.0), 1e-3));

    {
        vr_abc_t duty = vr_svpwm(0.0f, 0.0f, v_dc);
        CHECK(CLOSE(duty.a, 0.5, 1e-6));
        CHECK(CLOSE(duty.b, 0.5, 1e-6));
        CHECK(CLOSE(duty.c, 0.5, 1e-6));
    }

    for (k = 0; k < 360; k += 7) {
        vr_float th = (vr_float)(k * TEST_PI / 180.0);
        vr_float va = v_max * (vr_float)cos(th);
        vr_float vb = v_max * (vr_float)sin(th);
        vr_abc_t duty = vr_svpwm(va, vb, v_dc);
        vr_abc_t vph;
        vr_ab_t rec;

        CHECK(duty.a >= 0.0f && duty.a <= 1.0f);
        CHECK(duty.b >= 0.0f && duty.b <= 1.0f);
        CHECK(duty.c >= 0.0f && duty.c <= 1.0f);

        vph = vr_duty_to_phase_voltages(duty.a, duty.b, duty.c, v_dc);
        rec = vr_clarke(vph.a, vph.b, vph.c);
        CHECK(CLOSE(rec.alpha, va, 1e-2));
        CHECK(CLOSE(rec.beta, vb, 1e-2));
    }
}

/* Run a closed-loop simulation and return the final speed. */
static vr_float run_closed_loop(vr_float speed_setpoint, vr_float step_time,
                                vr_float load, vr_float load_time,
                                vr_float duration, vr_float *flux_out,
                                vr_float *torque_out)
{
    vr_motor_params_t p;
    vr_foc_t foc;
    vr_induction_motor_t motor;
    vr_float dt = 1.0e-4f;
    long n, k;

    p.rs = 1.5f; p.rr = 1.2f; p.ls = 0.14f; p.lr = 0.14f; p.lm = 0.135f;
    p.pole_pairs = 2; p.inertia = 0.02f; p.friction = 0.001f;

    vr_foc_init_tuned(&foc, &p, dt, 600.0f, 1.0f, 20.0f, 0.0f, 0.0f);
    vr_motor_init(&motor, &p);

    n = (long)(duration / dt + 0.5);
    for (k = 0; k < n; ++k) {
        vr_float t = (vr_float)k * dt;
        vr_float w_ref = (t >= step_time) ? speed_setpoint : 0.0f;
        vr_float t_load = (t >= load_time) ? load : 0.0f;
        vr_abc_t i = vr_motor_phase_currents(&motor);
        vr_foc_output_t out = vr_foc_update(&foc, w_ref, i.a, i.b, i.c,
                                            motor.omega_m);
        vr_motor_step(&motor, out.v_ab.alpha, out.v_ab.beta, t_load, dt);
    }
    if (flux_out) *flux_out = vr_motor_flux_magnitude(&motor);
    if (torque_out) *torque_out = vr_motor_torque(&motor);
    return motor.omega_m;
}

static void test_closed_loop(void)
{
    vr_float flux = 0.0f, torque = 0.0f;

    /* Speed tracking. */
    vr_float w = run_closed_loop(120.0f, 0.3f, 0.0f, 10.0f, 1.0f, &flux, &torque);
    CHECK(CLOSE(w, 120.0, 1.5));
    CHECK(CLOSE(flux, 1.0, 0.05)); /* rotor flux established */

    /* Load rejection: speed holds, torque balances the load. */
    w = run_closed_loop(120.0f, 0.3f, 8.0f, 0.8f, 1.6f, &flux, &torque);
    CHECK(CLOSE(w, 120.0, 1.5));
    CHECK(CLOSE(torque, 8.0, 1.0));

    /* Speed reversal. */
    w = run_closed_loop(-100.0f, 0.3f, 0.0f, 10.0f, 1.4f, &flux, &torque);
    CHECK(CLOSE(w, -100.0, 2.5));
}

int main(void)
{
    printf("transforms...\n");  test_transforms();
    printf("pi...\n");          test_pi();
    printf("svpwm...\n");       test_svpwm();
    printf("closed loop...\n"); test_closed_loop();

    printf("\n%d checks, %d failures\n", g_checks, g_failures);
    return g_failures == 0 ? 0 : 1;
}

#include "llr_16qam.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

/* ── colour output ───────────────────────────────────────────── */
#define GRN  "\033[32m"
#define RED  "\033[31m"
#define YLW  "\033[33m"
#define RST  "\033[0m"

static int g_pass = 0, g_fail = 0;

static void chk_near(const char *label, float got, float ref, float tol)
{
    float err = fabsf(got - ref);
    int ok = (err <= tol);
    printf("  %-44s  %+8.4f  (ref %+8.4f  err %.1e)  [%s%s%s]\n",
           label, got, ref, err,
           ok ? GRN : RED, ok ? "PASS" : "FAIL", RST);
    ok ? g_pass++ : g_fail++;
}

static void chk_sign(const char *label, float llr, int tx_bit)
{
    int dec = (llr >= 0.0f) ? 0 : 1;
    int ok  = (dec == tx_bit);
    printf("  %-44s  LLR=%+7.3f  dec=%d  tx=%d  [%s%s%s]\n",
           label, llr, dec, tx_bit,
           ok ? GRN : RED, ok ? "PASS" : "FAIL", RST);
    ok ? g_pass++ : g_fail++;
}

static void chk_sign_quiet(float llr, int tx_bit)
{
    int ok = (((llr >= 0.0f) ? 0 : 1) == tx_bit);
    ok ? g_pass++ : g_fail++;
}

static void banner(const char *s)
{
    printf("\n%s=== %s ===%s\n", YLW, s, RST);
}

/* ──────────────── level / bit constants ─────────────────────── */

#define A    0.31622776601683793f   /*  1/sqrt(10) */
#define A2   0.63245553203367586f   /*  2/sqrt(10) */
#define A3   0.94868329805051381f   /*  3/sqrt(10) */

#define A64  0.15430334996209194f   /*  1/sqrt(42) */
#define A256 0.07669649888473704f   /* 1/sqrt(170) */

/* PAM4 symbols (16-QAM) */
static const struct { float si, sq; int b[4]; } C[16] = {
    { A3,  A3, {0,0,0,0}},  { A3,  A,  {0,0,0,1}},
    { A3, -A,  {0,0,1,1}},  { A3, -A3, {0,0,1,0}},
    { A,   A3, {0,1,0,0}},  { A,   A,  {0,1,0,1}},
    { A,  -A,  {0,1,1,1}},  { A,  -A3, {0,1,1,0}},
    {-A,   A3, {1,1,0,0}},  {-A,   A,  {1,1,0,1}},
    {-A,  -A,  {1,1,1,1}},  {-A,  -A3, {1,1,1,0}},
    {-A3,  A3, {1,0,0,0}},  {-A3,  A,  {1,0,0,1}},
    {-A3, -A,  {1,0,1,1}},  {-A3, -A3, {1,0,1,0}},
};

/* PAM8 levels and Gray bits (64-QAM) */
static const float P8L[8] = {
    7.0f*A64, 5.0f*A64, 3.0f*A64,      A64,
       -A64, -3.0f*A64, -5.0f*A64, -7.0f*A64
};
static const int P8B[8][3] = {
    {0,0,0},{0,0,1},{0,1,1},{0,1,0},
    {1,1,0},{1,1,1},{1,0,1},{1,0,0}
};

/* PAM16 levels and Gray bits (256-QAM) */
static const float P16L[16] = {
    15.0f*A256, 13.0f*A256, 11.0f*A256,  9.0f*A256,
     7.0f*A256,  5.0f*A256,  3.0f*A256,       A256,
        -A256, -3.0f*A256, -5.0f*A256, -7.0f*A256,
    -9.0f*A256,-11.0f*A256,-13.0f*A256,-15.0f*A256
};
static const int P16B[16][4] = {
    {0,0,0,0},{0,0,0,1},{0,0,1,1},{0,0,1,0},
    {0,1,1,0},{0,1,1,1},{0,1,0,1},{0,1,0,0},
    {1,1,0,0},{1,1,0,1},{1,1,1,1},{1,1,1,0},
    {1,0,1,0},{1,0,1,1},{1,0,0,1},{1,0,0,0}
};

/* ====================================================================
 * TEST 1 — 16-QAM Worked Example
 * ==================================================================== */
static void test_worked_example(void)
{
    banner("Test 1: 16-QAM Worked Example from Tutorial");

    cf32_t y   = { -0.28622776f,  0.87868330f };
    cf32_t h   = {  0.8f,         0.6f        };
    float  s2n = 0.1f;

    const int   tx[4]  = { 0, 1, 0, 0 };
    const float ref[4] = { 3.7723f, -4.2277f, 14.1279f, 3.0640f };
    const float tol    = 5e-3f;

    printf("  y      = %.5f + %.5fj\n", y.re, y.im);
    printf("  h      = %.1f + %.1fj   |h|^2 = %.4f\n",
           h.re, h.im, h.re*h.re + h.im*h.im);
    printf("  sigma2n = %.2f   =>  sigma2_eff = %.2f\n\n", s2n, s2n);

    float llr_ml[4], llr_lin[4];
    llr_16qam_maxlog(y, h, s2n, llr_ml);
    llr_16qam_linear(y, h, s2n, llr_lin);

    printf("  Max-Log method:\n");
    for (int k = 0; k < 4; k++) {
        char lbl[64];
        snprintf(lbl, sizeof lbl, "Lambda(b%d)  [tx=%d]", k, tx[k]);
        chk_near(lbl, llr_ml[k], ref[k], tol);
    }

    printf("\n  Piecewise-linear method:\n");
    for (int k = 0; k < 4; k++) {
        char lbl[64];
        snprintf(lbl, sizeof lbl, "Lambda(b%d)  [tx=%d]", k, tx[k]);
        chk_near(lbl, llr_lin[k], ref[k], tol);
    }
}

/* ====================================================================
 * TEST 2 — All 16 Noiseless 16-QAM Symbols
 * ==================================================================== */
static void test_all_symbols(void)
{
    banner("Test 2: All 16 Noiseless 16-QAM Symbols (h=1, no noise)");

    cf32_t h = { 1.0f, 0.0f };
    float  s2n = 0.001f;

    for (int i = 0; i < 16; i++) {
        cf32_t y = { C[i].si, C[i].sq };
        float llr[4];
        llr_16qam_maxlog(y, h, s2n, llr);

        for (int k = 0; k < 4; k++) {
            char lbl[64];
            snprintf(lbl, sizeof lbl, "sym[%02d]  b%d  (b0b1b2b3=%d%d%d%d)",
                     i, k, C[i].b[0], C[i].b[1], C[i].b[2], C[i].b[3]);
            chk_sign(lbl, llr[k], C[i].b[k]);
        }
    }
}

/* ====================================================================
 * TEST 3 — 16-QAM Decision Boundaries (LLR = 0 exactly)
 * ==================================================================== */
static void test_boundaries(void)
{
    banner("Test 3: 16-QAM Decision Boundaries (LLR must equal 0)");

    cf32_t h = { 1.0f, 0.0f };
    float  s2n = 0.1f;
    float  tol = 1e-6f;
    float  llr[4];

    { cf32_t y = { 0.0f,  A  }; llr_16qam_maxlog(y, h, s2n, llr);
      chk_near("b0  r_I =  0.0000  (sign boundary)",    llr[0], 0.0f, tol); }
    { cf32_t y = {  A2,   A  }; llr_16qam_maxlog(y, h, s2n, llr);
      chk_near("b1  r_I = +2A     (inner/outer +)",     llr[1], 0.0f, tol); }
    { cf32_t y = { -A2,   A  }; llr_16qam_maxlog(y, h, s2n, llr);
      chk_near("b1  r_I = -2A     (inner/outer -)",     llr[1], 0.0f, tol); }
    { cf32_t y = {  A,   0.0f}; llr_16qam_maxlog(y, h, s2n, llr);
      chk_near("b2  r_Q =  0.0000  (sign boundary)",    llr[2], 0.0f, tol); }
    { cf32_t y = {  A,    A2 }; llr_16qam_maxlog(y, h, s2n, llr);
      chk_near("b3  r_Q = +2A     (inner/outer +)",     llr[3], 0.0f, tol); }
    { cf32_t y = {  A,   -A2 }; llr_16qam_maxlog(y, h, s2n, llr);
      chk_near("b3  r_Q = -2A     (inner/outer -)",     llr[3], 0.0f, tol); }
}

/* ====================================================================
 * TEST 4 — 16-QAM Max-Log vs Piecewise-Linear Equivalence
 * ==================================================================== */
static void test_equivalence(void)
{
    banner("Test 4: 16-QAM Max-Log vs Piecewise-Linear Equivalence");

    cf32_t h  = { 0.7f, -0.5f };   /* |h|^2 = 0.74 */
    float  s2n = 0.15f;
    float  tol = 1e-4f;

    static const struct { float ri, rq; const char *desc; } pts[] = {
        { 0.15f,  0.80f, "I:inner  Q:outer" },
        {-0.50f,  0.30f, "I:inner  Q:inner" },
        { 0.75f, -0.20f, "I:outer  Q:inner" },
        {-0.85f, -0.70f, "I:outer  Q:outer" },
        { A2,    A,      "I:+2A boundary"   },
        {-A2,   -A,      "I:-2A boundary"   },
        { A,     0.0f,   "Q:0 boundary"     },
    };
    int np = (int)(sizeof pts / sizeof pts[0]);

    for (int i = 0; i < np; i++) {
        float ri = pts[i].ri, rq = pts[i].rq;
        cf32_t y = {
            h.re * ri - h.im * rq,
            h.im * ri + h.re * rq
        };

        float llr_ml[4], llr_lin[4];
        llr_16qam_maxlog(y, h, s2n, llr_ml);
        llr_16qam_linear(y, h, s2n, llr_lin);

        printf("  [%d] %s\n", i, pts[i].desc);
        for (int k = 0; k < 4; k++) {
            char lbl[64];
            snprintf(lbl, sizeof lbl, "    b%d  maxlog=%+8.4f", k, llr_ml[k]);
            chk_near(lbl, llr_lin[k], llr_ml[k], tol);
        }
    }
}

/* ====================================================================
 * TEST 5 — 16-QAM Bit Error Case
 * ==================================================================== */
static void test_error_case(void)
{
    banner("Test 5: 16-QAM Bit Error — LLR tracks reliability");

    cf32_t h  = { 1.0f, 0.0f };
    float  s2n = 0.1f;
    float  si  = A;

    printf("  Tx: b0=0  s_I = +A = %.5f   boundary at r_I = 0\n\n", si);
    printf("  %-10s  %-10s  %-12s  %-10s  %s\n",
           "noise_I", "r_I", "Lambda(b0)", "decision", "outcome");
    printf("  %s\n", "------------------------------------------------------------");

    float offsets[] = { 0.0f, -0.15f, -0.30f, -0.35f, -0.50f, -0.70f };
    int n = (int)(sizeof offsets / sizeof offsets[0]);

    int ok_all = 1;
    for (int i = 0; i < n; i++) {
        float noise_i = offsets[i];
        float ri = si + noise_i;
        cf32_t y  = { ri, A };
        float llr[4];
        llr_16qam_maxlog(y, h, s2n, llr);

        int dec = (llr[0] >= 0.0f) ? 0 : 1;
        int err = (dec != 0);

        printf("  %+8.4f    %+8.4f    %+10.4f    bit=%d       %s\n",
               noise_i, ri, llr[0], dec, err ? "BIT ERROR" : "correct");

        if (fabsf(ri) < 0.05f && fabsf(llr[0]) > 1.0f) ok_all = 0;
    }

    printf("\n  %-44s  [%s%s%s]\n",
           "LLR magnitude small near boundary",
           ok_all ? GRN : RED, ok_all ? "PASS" : "FAIL", RST);
    ok_all ? g_pass++ : g_fail++;
}

/* ====================================================================
 * TEST 6 — 16-QAM Channel Weighting
 * ==================================================================== */
static void test_channel_weighting(void)
{
    banner("Test 6: 16-QAM Channel Weighting — |LLR| scales with |h|");

    float s2n = 0.1f;
    float ri  = 0.20f;

    printf("  sigma2n = %.2f,  r_I = %.2f\n\n", s2n, ri);
    printf("  %-6s  %-12s  %-12s  %s\n",
           "|h|", "sigma2_eff", "Lambda(b0)", "direction");
    printf("  %s\n", "----------------------------------------------");

    float mags[] = { 0.3f, 0.5f, 0.8f, 1.0f, 1.5f, 2.0f };
    int   nm     = (int)(sizeof mags / sizeof mags[0]);
    float prev   = 0.0f;
    int   ok     = 1;

    for (int i = 0; i < nm; i++) {
        float hmag = mags[i];
        cf32_t h   = { hmag, 0.0f };
        cf32_t y   = { hmag * ri, hmag * A };
        float llr[4];
        llr_16qam_maxlog(y, h, s2n, llr);

        float sigma2_eff = s2n / (hmag * hmag);
        printf("  %-6.1f  %-12.5f  %+10.4f", hmag, sigma2_eff, llr[0]);

        if (i > 0) {
            if (llr[0] > prev) { printf("   up\n"); }
            else               { printf("   NOT monotone\n"); ok = 0; }
        } else {
            printf("\n");
        }
        prev = llr[0];
    }

    printf("\n  %-44s  [%s%s%s]\n",
           "|LLR| strictly increases with |h|",
           ok ? GRN : RED, ok ? "PASS" : "FAIL", RST);
    ok ? g_pass++ : g_fail++;
}

/* ====================================================================
 * TEST 7 — 64-QAM Worked Example
 *
 *  Tx symbol: I = +5a, Q = +1a  (a = 1/sqrt(42))
 *  Tx bits:   I → Gray 001 (b0=0,b1=0,b2=1)
 *             Q → Gray 010 (b0=0,b1=1,b2=0)
 *  Channel  : h = 1+0j,  σ²_n = 0.05,  σ²_eff = 0.05
 *
 *  Reference LLRs  (scale = 4a²/σ² = 4/(42·0.05) = 40/21):
 *    b0_I = 9 · 40/21 =  120/7   ≈ +17.143
 *    b1_I = 1 · 40/21 =   40/21  ≈  +1.905
 *    b2_I = –1 · 40/21 = –40/21  ≈  –1.905
 *    b0_Q = 1 · 40/21 =   40/21  ≈  +1.905
 *    b1_Q = –4 · 40/21 = –160/21 ≈  –7.619
 *    b2_Q = 1 · 40/21 =   40/21  ≈  +1.905
 * ==================================================================== */
static void test_64qam_worked(void)
{
    banner("Test 7: 64-QAM Worked Example");

    cf32_t y   = { 5.0f*A64,  A64 };   /* noiseless, h=1 */
    cf32_t h   = { 1.0f, 0.0f };
    float  s2n = 0.05f;

    const int   tx[6]  = { 0,0,1, 0,1,0 };
    const float ref[6] = {
        120.0f/7.0f, 40.0f/21.0f, -40.0f/21.0f,
         40.0f/21.0f, -160.0f/21.0f, 40.0f/21.0f
    };
    const float tol = 5e-3f;

    printf("  Tx I=+5a, Q=+1a  (a=1/sqrt(42))   sigma2n=%.2f\n\n", s2n);

    float llr_ml[6], llr_lin[6];
    llr_64qam_maxlog(y, h, s2n, llr_ml);
    llr_64qam_linear(y, h, s2n, llr_lin);

    static const char *bname[6] = {"b0_I","b1_I","b2_I","b0_Q","b1_Q","b2_Q"};
    printf("  Max-Log method:\n");
    for (int k = 0; k < 6; k++) {
        char lbl[64];
        snprintf(lbl, sizeof lbl, "Lambda(%s)  [tx=%d]", bname[k], tx[k]);
        chk_near(lbl, llr_ml[k], ref[k], tol);
    }
    printf("\n  Piecewise-linear method:\n");
    for (int k = 0; k < 6; k++) {
        char lbl[64];
        snprintf(lbl, sizeof lbl, "Lambda(%s)  [tx=%d]", bname[k], tx[k]);
        chk_near(lbl, llr_lin[k], ref[k], tol);
    }
}

/* ====================================================================
 * TEST 8 — All 64 Noiseless 64-QAM Symbols
 * ==================================================================== */
static void test_64qam_all_symbols(void)
{
    banner("Test 8: All 64 Noiseless 64-QAM Symbols");

    cf32_t h = { 1.0f, 0.0f };
    float  s2n = 0.001f;

    for (int i = 0; i < 8; i++) {
        for (int j = 0; j < 8; j++) {
            cf32_t y = { P8L[i], P8L[j] };
            float llr[6];
            llr_64qam_maxlog(y, h, s2n, llr);

            int local_ok = 1;
            for (int k = 0; k < 3; k++) {
                int ok_i = (((llr[k]   >= 0.0f) ? 0 : 1) == P8B[i][k]);
                int ok_q = (((llr[3+k] >= 0.0f) ? 0 : 1) == P8B[j][k]);
                g_pass += (ok_i + ok_q);
                g_fail += (!ok_i + !ok_q);
                local_ok = local_ok && ok_i && ok_q;
            }
            printf("  [%d%d] I=%+6.4f Q=%+6.4f  "
                   "I-bits=%d%d%d  Q-bits=%d%d%d  [%s%s%s]\n",
                   i, j, P8L[i], P8L[j],
                   P8B[i][0], P8B[i][1], P8B[i][2],
                   P8B[j][0], P8B[j][1], P8B[j][2],
                   local_ok ? GRN : RED, local_ok ? "PASS" : "FAIL", RST);
        }
    }
}

/* ====================================================================
 * TEST 9 — 64-QAM Decision Boundaries (LLR = 0 exactly)
 *
 *  b0 boundary:  r = 0
 *  b1 boundary:  |r| = 4a
 *  b2 boundaries:|r| = 2a, 6a
 * ==================================================================== */
static void test_64qam_boundaries(void)
{
    banner("Test 9: 64-QAM Decision Boundaries (LLR must equal 0)");

    cf32_t h = { 1.0f, 0.0f };
    float  s2n = 0.1f;
    float  tol = 1e-5f;
    float  llr[6];
    float  a = A64;

    /* I-axis */
    { cf32_t y = {      0.0f, a }; llr_64qam_maxlog(y, h, s2n, llr);
      chk_near("b0_I  r_I=0       (sign boundary)",     llr[0], 0.0f, tol); }
    { cf32_t y = {  4.0f*a,   a }; llr_64qam_maxlog(y, h, s2n, llr);
      chk_near("b1_I  r_I=+4a     (b1 boundary)",       llr[1], 0.0f, tol); }
    { cf32_t y = { -4.0f*a,   a }; llr_64qam_maxlog(y, h, s2n, llr);
      chk_near("b1_I  r_I=-4a     (b1 boundary)",       llr[1], 0.0f, tol); }
    { cf32_t y = {  2.0f*a,   a }; llr_64qam_maxlog(y, h, s2n, llr);
      chk_near("b2_I  r_I=+2a     (b2 inner bndry)",    llr[2], 0.0f, tol); }
    { cf32_t y = { -2.0f*a,   a }; llr_64qam_maxlog(y, h, s2n, llr);
      chk_near("b2_I  r_I=-2a     (b2 inner bndry)",    llr[2], 0.0f, tol); }
    { cf32_t y = {  6.0f*a,   a }; llr_64qam_maxlog(y, h, s2n, llr);
      chk_near("b2_I  r_I=+6a     (b2 outer bndry)",    llr[2], 0.0f, tol); }
    { cf32_t y = { -6.0f*a,   a }; llr_64qam_maxlog(y, h, s2n, llr);
      chk_near("b2_I  r_I=-6a     (b2 outer bndry)",    llr[2], 0.0f, tol); }

    /* Q-axis */
    { cf32_t y = { a,      0.0f }; llr_64qam_maxlog(y, h, s2n, llr);
      chk_near("b3_Q  r_Q=0       (sign boundary)",     llr[3], 0.0f, tol); }
    { cf32_t y = { a,  4.0f*a   }; llr_64qam_maxlog(y, h, s2n, llr);
      chk_near("b4_Q  r_Q=+4a     (b1 boundary)",       llr[4], 0.0f, tol); }
    { cf32_t y = { a,  2.0f*a   }; llr_64qam_maxlog(y, h, s2n, llr);
      chk_near("b5_Q  r_Q=+2a     (b2 inner bndry)",    llr[5], 0.0f, tol); }
    { cf32_t y = { a,  6.0f*a   }; llr_64qam_maxlog(y, h, s2n, llr);
      chk_near("b5_Q  r_Q=+6a     (b2 outer bndry)",    llr[5], 0.0f, tol); }
}

/* ====================================================================
 * TEST 10 — 64-QAM Max-Log vs Piecewise-Linear Equivalence
 * ==================================================================== */
static void test_64qam_equivalence(void)
{
    banner("Test 10: 64-QAM Max-Log vs Piecewise-Linear Equivalence");

    cf32_t h  = { 0.6f, -0.7f };   /* |h|^2 = 0.85 */
    float  s2n = 0.12f;
    float  tol = 1e-4f;
    float  a   = A64;

    const struct { float ri, rq; const char *desc; } pts[] = {
        { 0.05f,   0.35f,  "I:very_inner  Q:mid"      },
        {-0.20f,   0.55f,  "I:inner       Q:outer"    },
        { 0.42f,  -0.12f,  "I:mid         Q:inner"    },
        {-0.70f,  -0.55f,  "I:outer       Q:mid"      },
        { 4.0f*a,     a,   "I:b1_boundary Q:+1a"      },
        { 2.0f*a, 6.0f*a,  "I:b2_low      Q:b2_high"  },
        {-6.0f*a, 4.0f*a,  "I:b2_hi_neg   Q:b1_bndry" },
    };
    int np = (int)(sizeof pts / sizeof pts[0]);

    for (int i = 0; i < np; i++) {
        float ri = pts[i].ri, rq = pts[i].rq;
        cf32_t y = {
            h.re * ri - h.im * rq,
            h.im * ri + h.re * rq
        };

        float llr_ml[6], llr_lin[6];
        llr_64qam_maxlog(y, h, s2n, llr_ml);
        llr_64qam_linear(y, h, s2n, llr_lin);

        printf("  [%d] %s\n", i, pts[i].desc);
        static const char *bname[6] = {"b0_I","b1_I","b2_I","b0_Q","b1_Q","b2_Q"};
        for (int k = 0; k < 6; k++) {
            char lbl[64];
            snprintf(lbl, sizeof lbl, "    %s  maxlog=%+8.4f", bname[k], llr_ml[k]);
            chk_near(lbl, llr_lin[k], llr_ml[k], tol);
        }
    }
}

/* ====================================================================
 * TEST 11 — 256-QAM Worked Example
 *
 *  Tx symbol: I = +5a, Q = −9a  (a = 1/sqrt(170))
 *  Tx bits:   I → 0111 (b0=0,b1=1,b2=1,b3=1)
 *             Q → 1010 (b0=1,b1=0,b2=1,b3=0)
 *  Channel  : h = 1+0j,  σ²_n = 0.03
 *
 *  Reference LLRs  (scale = 4a²/σ² = 4/(170·0.03) = 40/51 ≈ 0.784):
 *    b0_I = 9·(40/51) =  120/17 ≈ +7.059
 *    b1_I = −4·(40/51) = −160/51 ≈ −3.137
 *    b2_I = b3_I = −1·(40/51) = −40/51 ≈ −0.784
 *    b0_Q = −25·(40/51) = −1000/51 ≈ −19.608
 *    b1_Q = b3_Q = 1·(40/51) =  40/51 ≈ +0.784
 *    b2_Q = −4·(40/51) = −160/51 ≈ −3.137
 * ==================================================================== */
static void test_256qam_worked(void)
{
    banner("Test 11: 256-QAM Worked Example");

    cf32_t y   = { 5.0f*A256, -9.0f*A256 };
    cf32_t h   = { 1.0f, 0.0f };
    float  s2n = 0.03f;

    const int   tx[8]  = { 0,1,1,1, 1,0,1,0 };
    const float ref[8] = {
        120.0f/17.0f,  -160.0f/51.0f, -40.0f/51.0f, -40.0f/51.0f,
       -1000.0f/51.0f,   40.0f/51.0f,-160.0f/51.0f,  40.0f/51.0f
    };
    const float tol = 5e-3f;

    printf("  Tx I=+5a, Q=-9a  (a=1/sqrt(170))   sigma2n=%.2f\n\n", s2n);

    float llr_ml[8], llr_lin[8];
    llr_256qam_maxlog(y, h, s2n, llr_ml);
    llr_256qam_linear(y, h, s2n, llr_lin);

    static const char *bname[8] = {
        "b0_I","b1_I","b2_I","b3_I","b0_Q","b1_Q","b2_Q","b3_Q"
    };
    printf("  Max-Log method:\n");
    for (int k = 0; k < 8; k++) {
        char lbl[64];
        snprintf(lbl, sizeof lbl, "Lambda(%s)  [tx=%d]", bname[k], tx[k]);
        chk_near(lbl, llr_ml[k], ref[k], tol);
    }
    printf("\n  Piecewise-linear method:\n");
    for (int k = 0; k < 8; k++) {
        char lbl[64];
        snprintf(lbl, sizeof lbl, "Lambda(%s)  [tx=%d]", bname[k], tx[k]);
        chk_near(lbl, llr_lin[k], ref[k], tol);
    }
}

/* ====================================================================
 * TEST 12 — All 256 Noiseless 256-QAM Symbols (quiet)
 * ==================================================================== */
static void test_256qam_all_symbols(void)
{
    banner("Test 12: All 256 Noiseless 256-QAM Symbols (quiet, 2048 checks)");

    cf32_t h = { 1.0f, 0.0f };
    float  s2n = 0.001f;

    int prev_pass = g_pass, prev_fail = g_fail;
    for (int i = 0; i < 16; i++) {
        for (int j = 0; j < 16; j++) {
            cf32_t y = { P16L[i], P16L[j] };
            float llr[8];
            llr_256qam_maxlog(y, h, s2n, llr);
            for (int k = 0; k < 4; k++) {
                chk_sign_quiet(llr[k],   P16B[i][k]);
                chk_sign_quiet(llr[4+k], P16B[j][k]);
            }
        }
    }
    int np = g_pass - prev_pass;
    int nf = g_fail - prev_fail;
    printf("  256×256×8 = 2048 sign checks:  "
           "%d passed, %d failed  [%s%s%s]\n",
           np, nf,
           (nf == 0) ? GRN : RED, (nf == 0) ? "PASS" : "FAIL", RST);
}

/* ====================================================================
 * TEST 13 — 256-QAM Decision Boundaries (LLR = 0 exactly)
 *
 *  b0:  |r| = 0
 *  b1:  |r| = 8a
 *  b2:  |r| = 4a,  12a
 *  b3:  |r| = 2a,  6a,  10a,  14a
 * ==================================================================== */
static void test_256qam_boundaries(void)
{
    banner("Test 13: 256-QAM Decision Boundaries (LLR must equal 0)");

    cf32_t h = { 1.0f, 0.0f };
    float  s2n = 0.1f;
    float  tol = 1e-4f;
    float  llr[8];
    float  a = A256;

    { cf32_t y = {       0.0f, a }; llr_256qam_maxlog(y, h, s2n, llr);
      chk_near("b0_I  r_I=0        (sign boundary)",     llr[0], 0.0f, tol); }
    { cf32_t y = {   8.0f*a,   a }; llr_256qam_maxlog(y, h, s2n, llr);
      chk_near("b1_I  r_I=+8a      (b1 boundary)",       llr[1], 0.0f, tol); }
    { cf32_t y = {  -8.0f*a,   a }; llr_256qam_maxlog(y, h, s2n, llr);
      chk_near("b1_I  r_I=-8a      (b1 boundary)",       llr[1], 0.0f, tol); }
    { cf32_t y = {   4.0f*a,   a }; llr_256qam_maxlog(y, h, s2n, llr);
      chk_near("b2_I  r_I=+4a      (b2 inner bndry)",    llr[2], 0.0f, tol); }
    { cf32_t y = {  12.0f*a,   a }; llr_256qam_maxlog(y, h, s2n, llr);
      chk_near("b2_I  r_I=+12a     (b2 outer bndry)",    llr[2], 0.0f, tol); }
    { cf32_t y = {   2.0f*a,   a }; llr_256qam_maxlog(y, h, s2n, llr);
      chk_near("b3_I  r_I=+2a      (b3 bndry 1)",        llr[3], 0.0f, tol); }
    { cf32_t y = {   6.0f*a,   a }; llr_256qam_maxlog(y, h, s2n, llr);
      chk_near("b3_I  r_I=+6a      (b3 bndry 2)",        llr[3], 0.0f, tol); }
    { cf32_t y = {  10.0f*a,   a }; llr_256qam_maxlog(y, h, s2n, llr);
      chk_near("b3_I  r_I=+10a     (b3 bndry 3)",        llr[3], 0.0f, tol); }
    { cf32_t y = {  14.0f*a,   a }; llr_256qam_maxlog(y, h, s2n, llr);
      chk_near("b3_I  r_I=+14a     (b3 bndry 4)",        llr[3], 0.0f, tol); }
    { cf32_t y = { a,      0.0f  }; llr_256qam_maxlog(y, h, s2n, llr);
      chk_near("b4_Q  r_Q=0        (sign boundary)",     llr[4], 0.0f, tol); }
    { cf32_t y = { a,   8.0f*a   }; llr_256qam_maxlog(y, h, s2n, llr);
      chk_near("b5_Q  r_Q=+8a      (b1 boundary)",       llr[5], 0.0f, tol); }
    { cf32_t y = { a,  12.0f*a   }; llr_256qam_maxlog(y, h, s2n, llr);
      chk_near("b6_Q  r_Q=+12a     (b2 outer bndry)",    llr[6], 0.0f, tol); }
    { cf32_t y = { a,   6.0f*a   }; llr_256qam_maxlog(y, h, s2n, llr);
      chk_near("b7_Q  r_Q=+6a      (b3 bndry 2)",        llr[7], 0.0f, tol); }
}

/* ====================================================================
 * TEST 14 — 256-QAM Max-Log vs Piecewise-Linear Equivalence
 * ==================================================================== */
static void test_256qam_equivalence(void)
{
    banner("Test 14: 256-QAM Max-Log vs Piecewise-Linear Equivalence");

    cf32_t h  = { 0.5f, -0.8f };   /* |h|^2 = 0.89 */
    float  s2n = 0.08f;
    float  tol = 1e-4f;
    float  a   = A256;

    const struct { float ri, rq; const char *desc; } pts[] = {
        { 0.03f,    0.25f,  "I:very_inner Q:mid"         },
        {-0.12f,    0.60f,  "I:inner      Q:outer"       },
        { 0.55f,   -0.08f,  "I:mid-outer  Q:inner"       },
        {-0.80f,   -0.70f,  "I:outer      Q:outer"       },
        { 8.0f*a,   3.0f*a, "I:b1_bndry   Q:inner"       },
        { 4.0f*a,  12.0f*a, "I:b2_inner   Q:b2_outer"    },
        { 2.0f*a,  10.0f*a, "I:b3_bndry1  Q:b3_bndry3"  },
    };
    int np = (int)(sizeof pts / sizeof pts[0]);

    for (int i = 0; i < np; i++) {
        float ri = pts[i].ri, rq = pts[i].rq;
        cf32_t y = {
            h.re * ri - h.im * rq,
            h.im * ri + h.re * rq
        };

        float llr_ml[8], llr_lin[8];
        llr_256qam_maxlog(y, h, s2n, llr_ml);
        llr_256qam_linear(y, h, s2n, llr_lin);

        printf("  [%d] %s\n", i, pts[i].desc);
        static const char *bname[8] = {
            "b0_I","b1_I","b2_I","b3_I","b0_Q","b1_Q","b2_Q","b3_Q"
        };
        for (int k = 0; k < 8; k++) {
            char lbl[64];
            snprintf(lbl, sizeof lbl, "    %s  maxlog=%+9.4f", bname[k], llr_ml[k]);
            chk_near(lbl, llr_lin[k], llr_ml[k], tol);
        }
    }
}

/* ================================================================== */
int main(void)
{
    printf("16/64/256-QAM LLR Demapper — Test Suite\n");
    printf("=========================================\n");

    /* ── 16-QAM ── */
    test_worked_example();
    test_all_symbols();
    test_boundaries();
    test_equivalence();
    test_error_case();
    test_channel_weighting();

    /* ── 64-QAM ── */
    test_64qam_worked();
    test_64qam_all_symbols();
    test_64qam_boundaries();
    test_64qam_equivalence();

    /* ── 256-QAM ── */
    test_256qam_worked();
    test_256qam_all_symbols();
    test_256qam_boundaries();
    test_256qam_equivalence();

    printf("\n=========================================\n");
    if (g_fail == 0)
        printf(GRN "All %d tests passed.\n" RST, g_pass);
    else
        printf(RED "%d/%d tests FAILED.\n" RST, g_fail, g_pass + g_fail);
    printf("=========================================\n");

    return (g_fail > 0) ? EXIT_FAILURE : EXIT_SUCCESS;
}

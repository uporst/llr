#include "llr_16qam.h"
#include <math.h>
#include <float.h>

/* ── amplitude constants ──────────────────────────────────────── */
#define A    0.31622776601683793f   /*  1/sqrt(10)  — 16-QAM */
#define A2   0.63245553203367586f   /*  2/sqrt(10) */
#define A64  0.15430334996209194f   /*  1/sqrt(42)  — 64-QAM */
#define A256 0.07669649888473704f   /*  1/sqrt(170) — 256-QAM */

/* ── generic PAM descriptor ───────────────────────────────────── */
typedef struct {
    int          nlev;
    int          nbits;
    const float *level;
    const int   *bit;   /* row-major [nlev * nbits]: bit[i*nbits + b] */
} pam_t;

/* ── PAM4 tables (16-QAM) ─────────────────────────────────────── */
/* Gray codes: +3a=00, +1a=01, −1a=11, −3a=10 */
static const float pam4_level[4] = {
    3.0f*A,  A,  -A, -3.0f*A
};
static const int pam4_bit[4*2] = {
    0,0,  0,1,  1,1,  1,0
};

/* ── PAM8 tables (64-QAM) ─────────────────────────────────────── */
/* Gray codes: +7a=000, +5a=001, +3a=011, +1a=010,
 *             -1a=110, -3a=111, -5a=101, -7a=100 */
static const float pam8_level[8] = {
    7.0f*A64,  5.0f*A64,  3.0f*A64,       A64,
       -A64, -3.0f*A64, -5.0f*A64, -7.0f*A64
};
static const int pam8_bit[8*3] = {
    0,0,0,  0,0,1,  0,1,1,  0,1,0,
    1,1,0,  1,1,1,  1,0,1,  1,0,0
};

/* ── PAM16 tables (256-QAM) ───────────────────────────────────── */
/* Gray codes: +15a=0000, +13a=0001, +11a=0011, +9a=0010,
 *              +7a=0110,  +5a=0111,  +3a=0101, +1a=0100,
 *              -1a=1100,  -3a=1101,  -5a=1111, -7a=1110,
 *              -9a=1010, -11a=1011, -13a=1001,-15a=1000 */
static const float pam16_level[16] = {
    15.0f*A256, 13.0f*A256, 11.0f*A256,  9.0f*A256,
     7.0f*A256,  5.0f*A256,  3.0f*A256,       A256,
        -A256, -3.0f*A256, -5.0f*A256, -7.0f*A256,
    -9.0f*A256,-11.0f*A256,-13.0f*A256,-15.0f*A256
};
static const int pam16_bit[16*4] = {
    0,0,0,0,  0,0,0,1,  0,0,1,1,  0,0,1,0,
    0,1,1,0,  0,1,1,1,  0,1,0,1,  0,1,0,0,
    1,1,0,0,  1,1,0,1,  1,1,1,1,  1,1,1,0,
    1,0,1,0,  1,0,1,1,  1,0,0,1,  1,0,0,0
};

static const pam_t pam4d  = {  4, 2, pam4_level,  pam4_bit  };
static const pam_t pam8d  = {  8, 3, pam8_level,  pam8_bit  };
static const pam_t pam16d = { 16, 4, pam16_level, pam16_bit };

/* ── ZF equalization ──────────────────────────────────────────── */
static void zf_equalize(cf32_t y, cf32_t h,
                         float *ri, float *rq, float *sigma2_eff,
                         float sigma2n)
{
    float h2    = h.re*h.re + h.im*h.im;
    *ri         = (y.re*h.re + y.im*h.im) / h2;
    *rq         = (y.im*h.re - y.re*h.im) / h2;
    *sigma2_eff = sigma2n / h2;
}

/* ── Generic max-log LLR for one PAM axis ─────────────────────── */
static void axis_llr_maxlog_n(const pam_t *p, float r,
                               float inv_sig2, float *llr)
{
    float min_d2[4][2] = {
        {FLT_MAX, FLT_MAX}, {FLT_MAX, FLT_MAX},
        {FLT_MAX, FLT_MAX}, {FLT_MAX, FLT_MAX}
    };
    for (int i = 0; i < p->nlev; i++) {
        float d2 = (r - p->level[i]) * (r - p->level[i]);
        for (int b = 0; b < p->nbits; b++) {
            int v = p->bit[i * p->nbits + b];
            if (d2 < min_d2[b][v])
                min_d2[b][v] = d2;
        }
    }
    for (int b = 0; b < p->nbits; b++)
        llr[b] = (min_d2[b][1] - min_d2[b][0]) * inv_sig2;
}

/* ── 16-QAM piecewise-linear closed form ─────────────────────── *
 *
 *  MSB (b0, sign bit):
 *    r  < −2a :  8a(r + a)          kink at ±2a
 *   |r| ≤  2a :  4a·r               linear through 0
 *    r  >  2a :  8a(r − a)
 *
 *  LSB (b1, magnitude bit):
 *    4a(|r| − 2a)    — zero at the inner/outer boundary ±2a
 * ─────────────────────────────────────────────────────────────── */
static void axis_llr_linear(float r, float inv_sig2, float llr[2])
{
    float absr = fabsf(r);
    float l0;
    if      (r < -A2) l0 = 8.0f * A * (r + A);
    else if (r <=  A2) l0 = 4.0f * A *  r;
    else               l0 = 8.0f * A * (r - A);
    float l1 = 4.0f * A * (absr - A2);
    llr[0] = l0 * inv_sig2;
    llr[1] = l1 * inv_sig2;
}

/* ── 64-QAM piecewise-linear closed form ─────────────────────── *
 *
 *  b0 (sign):   generic formula — see axis_llr_256qam_linear
 *
 *  b1 (second):
 *    |r| < 2a :  8a(|r| − 3a)   negative in inner  → bit=1
 *   2a ≤ |r| ≤ 6a :  4a(|r| − 4a)  zero at ±4a boundary
 *    |r| > 6a :  8a(|r| − 5a)   positive in outer  → bit=0
 *
 *  b2 (LSB):
 *   |r| ≤ 4a :  4a(2a − |r|)   zero at ±2a and ±6a
 *   |r|  > 4a :  4a(|r| − 6a)
 * ─────────────────────────────────────────────────────────────── */
static void axis_llr_64qam_linear(float r, float inv_sig2, float llr[3])
{
    float absr = fabsf(r);
    float a    = A64;
    float a2   = 2.0f * a;
    float a4   = 4.0f * a;
    float a6   = 6.0f * a;

    /* b0: generic sign-bit formula for 4 level-pairs */
    {
        int   k   = (int)(absr / a2);  if (k > 3) k = 3;
        float sgn = (r >= 0.0f) ? 1.0f : -1.0f;
        llr[0] = sgn * 4.0f * a * (float)(k + 1)
                     * (absr - (float)k * a) * inv_sig2;
    }
    /* b1 */
    {
        float l1;
        if      (absr <  a2) l1 = 8.0f * a * (absr - 3.0f * a);
        else if (absr <= a6) l1 = 4.0f * a * (absr - 4.0f * a);
        else                 l1 = 8.0f * a * (absr - 5.0f * a);
        llr[1] = l1 * inv_sig2;
    }
    /* b2 */
    {
        float l2;
        if (absr <= a4) l2 = 4.0f * a * (2.0f * a - absr);
        else            l2 = 4.0f * a * (absr      - 6.0f * a);
        llr[2] = l2 * inv_sig2;
    }
}

/* ── 256-QAM piecewise-linear closed form ────────────────────── *
 *
 *  b0 (sign):
 *    k = min(⌊|r|/(2a)⌋, 7)
 *    Λ = sign(r)·4a·(k+1)·(|r| − k·a) / σ²_eff
 *
 *  b1 (second):
 *    Inner half (|r| < 8a):  k₁ = min(⌊|r|/(2a)⌋, 3)
 *      Λ = 4(4−k₁)a·(|r| − (5+k₁)a)       — neg inside, zero at 8a
 *    Outer half (|r| ≥ 8a):  k₂ = min(⌊(|r|−8a)/(2a)⌋, 3)
 *      Λ = 4(k₂+1)a·(|r| − (8+k₂)a)       — pos outside, zero at 8a
 *
 *  b2:
 *    |r| <  2a :  8a(3a − |r|)
 *   2a ≤ |r| ≤  6a :  4a(4a − |r|)        — zero at ±4a and ±12a
 *    6a < |r| <  8a :  8a(5a − |r|)
 *    8a ≤ |r| < 10a :  8a(|r| − 11a)
 *   10a ≤ |r| ≤ 14a :  4a(|r| − 12a)
 *        |r| > 14a :  8a(|r| − 13a)
 *
 *  b3 (LSB):
 *    |r| <  4a :  4a(2a − |r|)             — zero at ±2a, ±6a, ±10a, ±14a
 *   4a ≤ |r| <  8a :  4a(|r| − 6a)
 *   8a ≤ |r| < 12a :  4a(10a − |r|)
 *       |r| ≥ 12a :  4a(|r| − 14a)
 * ─────────────────────────────────────────────────────────────── */
static void axis_llr_256qam_linear(float r, float inv_sig2, float llr[4])
{
    float absr = fabsf(r);
    float a    = A256;
    float a2   = 2.0f * a;
    float a4   = 4.0f * a;
    float a6   = 6.0f * a;
    float a8   = 8.0f * a;

    /* b0: sign bit, 8 level-pairs */
    {
        int   k   = (int)(absr / a2);  if (k > 7) k = 7;
        float sgn = (r >= 0.0f) ? 1.0f : -1.0f;
        llr[0] = sgn * 4.0f * a * (float)(k + 1)
                     * (absr - (float)k * a) * inv_sig2;
    }
    /* b1 */
    {
        float l1;
        if (absr < a8) {
            int k1 = (int)(absr / a2);  if (k1 > 3) k1 = 3;
            l1 = 4.0f * (float)(4 - k1) * a
                      * (absr - (float)(5 + k1) * a);
        } else {
            int k2 = (int)((absr - a8) / a2);  if (k2 > 3) k2 = 3;
            l1 = 4.0f * (float)(k2 + 1) * a
                      * (absr - (float)(8 + k2) * a);
        }
        llr[1] = l1 * inv_sig2;
    }
    /* b2 */
    {
        float l2;
        if      (absr <   a2)          l2 =  8.0f * a * (3.0f * a - absr);
        else if (absr <=   a6)         l2 =  4.0f * a * (4.0f * a - absr);
        else if (absr <    a8)         l2 =  8.0f * a * (5.0f * a - absr);
        else if (absr < 10.0f * a)     l2 =  8.0f * a * (absr - 11.0f * a);
        else if (absr <= 14.0f * a)    l2 =  4.0f * a * (absr - 12.0f * a);
        else                            l2 =  8.0f * a * (absr - 13.0f * a);
        llr[2] = l2 * inv_sig2;
    }
    /* b3 */
    {
        float l3;
        if      (absr <   a4)          l3 = 4.0f * a * (2.0f * a - absr);
        else if (absr <   a8)          l3 = 4.0f * a * (absr      - a6);
        else if (absr < 12.0f * a)     l3 = 4.0f * a * (10.0f * a - absr);
        else                            l3 = 4.0f * a * (absr - 14.0f * a);
        llr[3] = l3 * inv_sig2;
    }
}

/* ── Public 16-QAM ───────────────────────────────────────────── */

void llr_16qam_maxlog(cf32_t y, cf32_t h, float sigma2n, float llr[4])
{
    float ri, rq, sigma2_eff;
    zf_equalize(y, h, &ri, &rq, &sigma2_eff, sigma2n);
    float inv = 1.0f / sigma2_eff;
    float li[2], lq[2];
    axis_llr_maxlog_n(&pam4d, ri, inv, li);
    axis_llr_maxlog_n(&pam4d, rq, inv, lq);
    llr[0] = li[0];  llr[1] = li[1];
    llr[2] = lq[0];  llr[3] = lq[1];
}

void llr_16qam_linear(cf32_t y, cf32_t h, float sigma2n, float llr[4])
{
    float ri, rq, sigma2_eff;
    zf_equalize(y, h, &ri, &rq, &sigma2_eff, sigma2n);
    float inv = 1.0f / sigma2_eff;
    float li[2], lq[2];
    axis_llr_linear(ri, inv, li);
    axis_llr_linear(rq, inv, lq);
    llr[0] = li[0];  llr[1] = li[1];
    llr[2] = lq[0];  llr[3] = lq[1];
}

/* ── Public 64-QAM ───────────────────────────────────────────── */

void llr_64qam_maxlog(cf32_t y, cf32_t h, float sigma2n, float llr[6])
{
    float ri, rq, sigma2_eff;
    zf_equalize(y, h, &ri, &rq, &sigma2_eff, sigma2n);
    float inv = 1.0f / sigma2_eff;
    float li[3], lq[3];
    axis_llr_maxlog_n(&pam8d, ri, inv, li);
    axis_llr_maxlog_n(&pam8d, rq, inv, lq);
    llr[0] = li[0];  llr[1] = li[1];  llr[2] = li[2];
    llr[3] = lq[0];  llr[4] = lq[1];  llr[5] = lq[2];
}

void llr_64qam_linear(cf32_t y, cf32_t h, float sigma2n, float llr[6])
{
    float ri, rq, sigma2_eff;
    zf_equalize(y, h, &ri, &rq, &sigma2_eff, sigma2n);
    float inv = 1.0f / sigma2_eff;
    float li[3], lq[3];
    axis_llr_64qam_linear(ri, inv, li);
    axis_llr_64qam_linear(rq, inv, lq);
    llr[0] = li[0];  llr[1] = li[1];  llr[2] = li[2];
    llr[3] = lq[0];  llr[4] = lq[1];  llr[5] = lq[2];
}

/* ── Public 256-QAM ──────────────────────────────────────────── */

void llr_256qam_maxlog(cf32_t y, cf32_t h, float sigma2n, float llr[8])
{
    float ri, rq, sigma2_eff;
    zf_equalize(y, h, &ri, &rq, &sigma2_eff, sigma2n);
    float inv = 1.0f / sigma2_eff;
    float li[4], lq[4];
    axis_llr_maxlog_n(&pam16d, ri, inv, li);
    axis_llr_maxlog_n(&pam16d, rq, inv, lq);
    llr[0]=li[0]; llr[1]=li[1]; llr[2]=li[2]; llr[3]=li[3];
    llr[4]=lq[0]; llr[5]=lq[1]; llr[6]=lq[2]; llr[7]=lq[3];
}

void llr_256qam_linear(cf32_t y, cf32_t h, float sigma2n, float llr[8])
{
    float ri, rq, sigma2_eff;
    zf_equalize(y, h, &ri, &rq, &sigma2_eff, sigma2n);
    float inv = 1.0f / sigma2_eff;
    float li[4], lq[4];
    axis_llr_256qam_linear(ri, inv, li);
    axis_llr_256qam_linear(rq, inv, lq);
    llr[0]=li[0]; llr[1]=li[1]; llr[2]=li[2]; llr[3]=li[3];
    llr[4]=lq[0]; llr[5]=lq[1]; llr[6]=lq[2]; llr[7]=lq[3];
}

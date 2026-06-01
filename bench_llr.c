/* bench_llr.c — RDTSC cycle-count benchmark for all LLR demapper variants
 *
 * Build:   make bench
 * Run:     ./bench_llr
 *
 * For stable results, pin the process to one core and disable frequency
 * scaling:  sudo cpupower frequency-set -g performance
 * Then run: taskset -c 0 ./bench_llr
 *
 * Requires x86 / x86-64 (RDTSC / RDTSCP).
 */
#include "llr_16qam.h"
#include <stdio.h>
#include <stdint.h>

#define WARMUP   100000u   /* iterations before timing starts */
#define ITERS   1000000u   /* timed iterations per function   */
#define NVEC        16u    /* input pool size (must be 2^k)   */

/* ── RDTSC helpers ────────────────────────────────────────────── */
static inline uint64_t tsc_start(void)
{
    uint32_t lo, hi;
    /* lfence: wait for all prior loads/stores to complete */
    __asm__ volatile ("lfence\n\trdtsc"
                      : "=a"(lo), "=d"(hi) :: "memory");
    return ((uint64_t)hi << 32) | lo;
}

static inline uint64_t tsc_end(void)
{
    uint32_t lo, hi, aux;
    /* rdtscp: serialising read; lfence: keep later loads in order */
    __asm__ volatile ("rdtscp\n\tlfence"
                      : "=a"(lo), "=d"(hi), "=c"(aux) :: "memory");
    return ((uint64_t)hi << 32) | lo;
}

/* ── 16 received-signal vectors (range ≈ ±1, h=0.8+0.6j implied) ─ */
static const cf32_t SIG[NVEC] = {
    { 0.25f,  0.75f}, { 0.75f,  0.25f}, {-0.50f,  0.50f}, { 0.50f, -0.50f},
    { 1.00f,  0.10f}, {-0.10f, -1.00f}, { 0.10f,  1.00f}, {-1.00f,  0.10f},
    { 0.33f,  0.67f}, {-0.67f,  0.33f}, { 0.80f, -0.40f}, {-0.40f, -0.80f},
    {-0.20f,  0.90f}, { 0.60f, -0.70f}, {-0.85f, -0.15f}, { 0.45f,  0.55f},
};
static const cf32_t H   = { 0.8f,  0.6f };
static const float  S2N = 0.05f;

/* ── global sink: forces all llr[] writes to stay live ────────── */
static float g_sink;

typedef void (*llr_fn_t)(cf32_t, cf32_t, float, float *);

/* Returns average TSC ticks per call. */
static double bench_one(llr_fn_t fn)
{
    float llr[8];
    float sink = 0.0f;

    /* warm up instruction cache, branch predictor, and TLB */
    for (unsigned i = 0; i < WARMUP; i++) {
        fn(SIG[i & (NVEC-1)], H, S2N, llr);
        sink += llr[0];
    }

    uint64_t t0 = tsc_start();
    for (unsigned i = 0; i < ITERS; i++) {
        fn(SIG[i & (NVEC-1)], H, S2N, llr);
        sink += llr[0];   /* read llr[0] to keep calls live */
    }
    uint64_t t1 = tsc_end();

    g_sink += sink;       /* prevent hoisting of the loop */
    return (double)(t1 - t0) / ITERS;
}

int main(void)
{
    static const struct { const char *name; llr_fn_t fn; } F[] = {
        { "16-QAM  max-log",  (llr_fn_t)llr_16qam_maxlog  },
        { "16-QAM  linear",   (llr_fn_t)llr_16qam_linear  },
        { "64-QAM  max-log",  (llr_fn_t)llr_64qam_maxlog  },
        { "64-QAM  linear",   (llr_fn_t)llr_64qam_linear  },
        { "256-QAM max-log",  (llr_fn_t)llr_256qam_maxlog },
        { "256-QAM linear",   (llr_fn_t)llr_256qam_linear },
    };
    int nf = (int)(sizeof F / sizeof F[0]);

    printf("LLR demapper — cycle benchmark\n");
    printf("  ITERS=%u   warmup=%u   pool=%u vectors\n",
           ITERS, WARMUP, NVEC);
    printf("  H = %.1f+%.1fj   sigma2n = %.2f\n\n",
           H.re, H.im, S2N);
    printf("  %-26s  %10s  %8s\n", "function", "cy/call", "ratio");
    printf("  %-26s  %10s  %8s\n",
           "--------------------------", "----------", "--------");

    double cy_base = 0.0;
    for (int i = 0; i < nf; i++) {
        double cy = bench_one(F[i].fn);
        if (i == 0) cy_base = cy;
        printf("  %-26s  %10.1f  %8.2f×\n",
               F[i].name, cy, cy / cy_base);
    }

    /* consume sink so the compiler can't discard the whole benchmark */
    printf("\n  (checksum %.4g)\n", g_sink);
    return 0;
}

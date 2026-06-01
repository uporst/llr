#ifndef LLR_16QAM_H
#define LLR_16QAM_H

/*
 * LLR demapper for 16-QAM, 64-QAM, and 256-QAM
 * Gray-coded, 5G NR bit-to-level convention.
 *
 * ZF equalization applied first:  yhat = y * conj(h) / |h|^2
 *   sigma2_eff = sigma2_n / |h|^2
 * Max-Log approximation:
 *   Lambda(bk) = [min_{S1} d^2 - min_{S0} d^2] / sigma2_eff
 * LLR sign convention:  positive => bit=0,  negative => bit=1
 *
 * Bit-to-level mapping (same structure on I and Q axes):
 *
 *  16-QAM  PAM4,  a = 1/sqrt(10):
 *    b0 b1 | level      b0 b1 | level
 *     0  0 | +3a         1  1 | -1a
 *     0  1 | +1a         1  0 | -3a
 *
 *  64-QAM  PAM8,  a = 1/sqrt(42):
 *    b0 b1 b2 | level     b0 b1 b2 | level
 *     0  0  0 | +7a        1  1  0 | -1a
 *     0  0  1 | +5a        1  1  1 | -3a
 *     0  1  1 | +3a        1  0  1 | -5a
 *     0  1  0 | +1a        1  0  0 | -7a
 *
 *  256-QAM  PAM16,  a = 1/sqrt(170):
 *    b0b1b2b3 | level     b0b1b2b3 | level
 *    0 0 0 0  | +15a      1 1 0 0  |  -1a
 *    0 0 0 1  | +13a      1 1 0 1  |  -3a
 *    0 0 1 1  | +11a      1 1 1 1  |  -5a
 *    0 0 1 0  |  +9a      1 1 1 0  |  -7a
 *    0 1 1 0  |  +7a      1 0 1 0  |  -9a
 *    0 1 1 1  |  +5a      1 0 1 1  | -11a
 *    0 1 0 1  |  +3a      1 0 0 1  | -13a
 *    0 1 0 0  |  +1a      1 0 0 0  | -15a
 *
 * Output ordering:  I-axis bits first, then Q-axis bits.
 *   16-QAM  : llr[0..3]  (llr[0-1] from I, llr[2-3] from Q)
 *   64-QAM  : llr[0..5]  (llr[0-2] from I, llr[3-5] from Q)
 *  256-QAM  : llr[0..7]  (llr[0-3] from I, llr[4-7] from Q)
 */

typedef struct { float re, im; } cf32_t;

/* -- 16-QAM  (4 output LLRs) ---------------------------------------- */
void llr_16qam_maxlog(cf32_t y, cf32_t h, float sigma2n, float llr[4]);
void llr_16qam_linear(cf32_t y, cf32_t h, float sigma2n, float llr[4]);

/* -- 64-QAM  (6 output LLRs) ---------------------------------------- */
void llr_64qam_maxlog(cf32_t y, cf32_t h, float sigma2n, float llr[6]);
void llr_64qam_linear(cf32_t y, cf32_t h, float sigma2n, float llr[6]);

/* -- 256-QAM  (8 output LLRs) --------------------------------------- */
void llr_256qam_maxlog(cf32_t y, cf32_t h, float sigma2n, float llr[8]);
void llr_256qam_linear(cf32_t y, cf32_t h, float sigma2n, float llr[8]);

#endif /* LLR_16QAM_H */

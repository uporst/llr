# llr

LLR demapper for 16-QAM, 64-QAM, and 256-QAM in C.

Gray-coded, 5G NR bit-to-level convention. Each modulation has two variants:

- **max-log** — generic engine driven by a PAM descriptor table, exact under
  the max-log approximation
- **linear** — closed-form piecewise-linear expressions, algebraically
  equivalent to max-log but ~3× faster

## API

```c
typedef struct { float re, im; } cf32_t;

void llr_16qam_maxlog (cf32_t y, cf32_t h, float sigma2n, float llr[4]);
void llr_16qam_linear (cf32_t y, cf32_t h, float sigma2n, float llr[4]);

void llr_64qam_maxlog (cf32_t y, cf32_t h, float sigma2n, float llr[6]);
void llr_64qam_linear (cf32_t y, cf32_t h, float sigma2n, float llr[6]);

void llr_256qam_maxlog(cf32_t y, cf32_t h, float sigma2n, float llr[8]);
void llr_256qam_linear(cf32_t y, cf32_t h, float sigma2n, float llr[8]);
```

`y` is the received complex sample, `h` the channel coefficient, `sigma2n`
the noise variance. ZF equalization (`yhat = y*conj(h)/|h|^2`) is applied
internally before demapping.

Output ordering: I-axis bits first, then Q-axis bits.

LLR sign convention: **positive => bit 0, negative => bit 1**.

## Build

```
make          # builds test_llr
make bench    # builds bench_llr
```

Requires gcc and `-lm`. No other dependencies.

## Test

```
./test_llr
```

14 test cases, 2690 checks: worked examples, all noiseless symbols, boundary
conditions, and max-log vs linear equivalence for all three modulations.

## Benchmark

```
taskset -c 0 ./bench_llr
```

Sample results (x86-64, O2):

| Function        | cy/call (avg 5 runs) | vs 16-QAM max-log |
|-----------------|---------------------:|------------------:|
| 16-QAM  max-log |                   71 |             1.00× |
| 16-QAM  linear  |                   14 |             0.19× |
| 64-QAM  max-log |                  133 |             1.88× |
| 64-QAM  linear  |                   33 |             0.46× |
| 256-QAM max-log |                  306 |             4.32× |
| 256-QAM linear  |                   45 |             0.63× |

For stable results, pin to one core and set the CPU governor to performance:

```
sudo cpupower frequency-set -g performance
taskset -c 0 ./bench_llr
```

## Bit mapping

Same Gray-coded structure on I and Q axes. Example for 64-QAM (PAM8,
`a = 1/sqrt(42)`):

```
b0 b1 b2 | level      b0 b1 b2 | level
 0  0  0  | +7a         1  1  0  | -1a
 0  0  1  | +5a         1  1  1  | -3a
 0  1  1  | +3a         1  0  1  | -5a
 0  1  0  | +1a         1  0  0  | -7a
```

16-QAM uses PAM4 (`a = 1/sqrt(10)`), 256-QAM uses PAM16 (`a = 1/sqrt(170)`).
See `llr_16qam.h` for the full tables.

CC      = gcc
CFLAGS  = -Wall -Wextra -O2 -std=c11
LDFLAGS = -lm

TARGET  = test_llr
BENCH   = bench_llr
HDRS    = llr_16qam.h

all: $(TARGET)

$(TARGET): llr_16qam.c test_llr_16qam.c $(HDRS)
	$(CC) $(CFLAGS) -o $@ llr_16qam.c test_llr_16qam.c $(LDFLAGS)

$(BENCH): llr_16qam.c bench_llr.c $(HDRS)
	$(CC) $(CFLAGS) -o $@ llr_16qam.c bench_llr.c $(LDFLAGS)

bench: $(BENCH)

clean:
	rm -f $(TARGET) $(BENCH)

.PHONY: all bench clean

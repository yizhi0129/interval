CC = gcc
CFLAGS = -Wall -O2 -std=c99 -frounding-math -ffp-contract=off
LDFLAGS = -lm
LDFLAGS2 = -lm -lopenblas
LDFLAGS3 = -lm -lmpfi -lmpfr -lgmp
LDFLAGS4 = -lm -lopenblas -lmpfi -lmpfr -lgmp

TARGETS = test_convert test_functions test_newton test_compression test_compression2 test_GS test_GS_CSR

test_convert_SRCS = test_convert.c convert.c
test_functions_SRCS = test_functions.c functions.c convert.c
test_newton_SRCS = test_newton.c functions.c convert.c
test_compression_SRCS = test_compression.c convert.c mpfr_interval.c
test_compression2_SRCS = test_compression2.c convert.c
test_GS_SRCS = test_GS.c convert.c functions.c mpfr_interval.c
test_GS_CSR_SRCS = test_GS_CSR.c convert.c functions.c mpfr_interval.c

test_convert_OBJS = $(test_convert_SRCS:.c=.o)
test_functions_OBJS = $(test_functions_SRCS:.c=.o)
test_newton_OBJS = $(test_newton_SRCS:.c=.o)
test_compression_OBJS = $(test_compression_SRCS:.c=.o)
test_compression2_OBJS = $(test_compression2_SRCS:.c=.o)
test_GS_OBJS = $(test_GS_SRCS:.c=.o)
test_GS_CSR_OBJS = $(test_GS_CSR_SRCS:.c=.o)

all: $(TARGETS)

test_convert: $(test_convert_OBJS)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

test_functions: $(test_functions_OBJS)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS2)

test_newton: $(test_newton_OBJS)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS2)

test_compression: $(test_compression_OBJS)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS3)

test_compression2: $(test_compression2_OBJS)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS3)

test_GS: $(test_GS_OBJS)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS4)

test_GS_CSR: $(test_GS_CSR_OBJS)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS4)

clean:
	rm -f *.o $(TARGETS)

.PHONY: all clean

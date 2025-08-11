CC = gcc
CFLAGS = -Wall -O2 -std=c99 -frounding-math -ffp-contract=off -Ilib
LDFLAGS = -lm
LDFLAGS2 = -lm -lopenblas
LDFLAGS3 = -lm -lmpfi -lmpfr -lgmp
LDFLAGS4 = -lm -lopenblas -lmpfi -lmpfr -lgmp
LDFLAGS5 = -lm -Llib -lsoftposit

TARGETS = test_convert test_functions test_newton test_compression test_compression2 test_GS test_GS_CSR test_posit_convert

test_convert_SRCS = test/test_convert.c src/convert.c
test_functions_SRCS = test/test_functions.c src/functions.c src/convert.c
test_newton_SRCS = test/test_newton.c src/functions.c src/convert.c
test_compression_SRCS = test/test_compression.c src/convert.c src/mpfr_interval.c
test_compression2_SRCS = test/test_compression2.c src/convert.c
test_GS_SRCS = test/test_GS.c src/convert.c src/functions.c src/mpfr_interval.c
test_GS_CSR_SRCS = test/test_GS_CSR.c src/convert.c src/functions.c src/mpfr_interval.c
test_posit_convert_SRCS = test/test_posit_convert.c src/posit_interval.c src/convert.c

test_convert_OBJS = $(test_convert_SRCS:.c=.o)
test_functions_OBJS = $(test_functions_SRCS:.c=.o)
test_newton_OBJS = $(test_newton_SRCS:.c=.o)
test_compression_OBJS = $(test_compression_SRCS:.c=.o)
test_compression2_OBJS = $(test_compression2_SRCS:.c=.o)
test_GS_OBJS = $(test_GS_SRCS:.c=.o)
test_GS_CSR_OBJS = $(test_GS_CSR_SRCS:.c=.o)
test_posit_convert_OBJS = $(test_posit_convert_SRCS:.c=.o)

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

test_posit_convert: $(test_posit_convert_OBJS)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS5)


clean:
	rm -f src/*.o test/*.o $(TARGETS)

.PHONY: all clean

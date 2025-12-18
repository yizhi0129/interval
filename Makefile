CC = gcc
CFLAGS = -Wall -O2 -std=c99 -fno-fast-math -ffloat-store -frounding-math -ffp-contract=off -Ilib
LDFLAGS = -lm
LDFLAGS2 = -lm -lopenblas
LDFLAGS3 = -lm -lmpfi -lmpfr -lgmp
LDFLAGS4 = -lm -lopenblas -lmpfi -lmpfr -lgmp
LDFLAGS5 = -lm -Llib -lsoftposit

TARGETS = test_convert test_matprod test_newton test_mpfr_convert test_storage test_GS test_GS_CSR test_posit_convert test_posit_app test_round

test_convert_SRCS = test/test_convert.c src/convert.c
test_matprod_SRCS = test/test_matprod.c src/functions.c src/convert.c
test_newton_SRCS = test/test_newton.c src/functions.c src/convert.c
test_mpfr_convert_SRCS = test/test_mpfr_convert.c src/convert.c src/mpfr_interval.c
test_storage_SRCS = test/test_storage.c src/convert.c
test_GS_SRCS = test/test_GS.c src/convert.c src/functions.c src/mpfr_interval.c
test_GS_CSR_SRCS = test/test_GS_CSR.c src/convert.c src/functions.c src/mpfr_interval.c
test_posit_convert_SRCS = test/test_posit_convert.c src/posit_interval.c src/convert.c
test_posit_app_SRCS = test/test_posit_app.c src/posit_interval.c src/convert.c
test_round_SRCS = test/test_round.c src/posit_interval.c src/convert.c

test_convert_OBJS = $(test_convert_SRCS:.c=.o)
test_matprod_OBJS = $(test_matprod_SRCS:.c=.o)
test_newton_OBJS = $(test_newton_SRCS:.c=.o)
test_mpfr_convert_OBJS = $(test_mpfr_convert_SRCS:.c=.o)
test_storage_OBJS = $(test_storage_SRCS:.c=.o)
test_GS_OBJS = $(test_GS_SRCS:.c=.o)
test_GS_CSR_OBJS = $(test_GS_CSR_SRCS:.c=.o)
test_posit_convert_OBJS = $(test_posit_convert_SRCS:.c=.o)
test_posit_app_OBJS = $(test_posit_app_SRCS:.c=.o)
test_round_OBJS = $(test_round_SRCS:.c=.o)

all: $(TARGETS)

test_convert: $(test_convert_OBJS)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

test_matprod: $(test_matprod_OBJS)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS2)

test_newton: $(test_newton_OBJS)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS2)

test_mpfr_convert: $(test_mpfr_convert_OBJS)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS3)

test_storage: $(test_storage_OBJS)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS3)

test_GS: $(test_GS_OBJS)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS4)

test_GS_CSR: $(test_GS_CSR_OBJS)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS4)

test_posit_convert: $(test_posit_convert_OBJS)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS5)

test_posit_app: $(test_posit_app_OBJS)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS5)

test_round: $(test_round_OBJS)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS5)

clean:
	rm -f src/*.o test/*.o $(TARGETS)

.PHONY: all clean

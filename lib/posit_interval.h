#ifndef _POSIT_INTERVAL_H_
#define _POSIT_INTERVAL_H_

#include <stdint.h>
#include <stdio.h>
#include <stdint.h>

#include "softposit.h"

#define CONST_IPB_8 set_pow2(-32)
#define CONST_IPB_16 set_pow2(-64)

typedef struct {
    posit16_t posit16;
    uint8_t prec;
} p16_FP_INT_P;

typedef struct {
    posit32_t posit32;
    uint8_t prec;
} p32_FP_INT_P;

void p8_print_binary(posit8_t x);
void p16_print_binary(posit16_t x);
void p32_print_binary(posit32_t x);

posit8_t p8_read_3r(posit8_t c);
posit16_t p16_read_3r(posit16_t c);
posit32_t p32_read_3r(posit32_t c);

p16_FP_INT_P p8_compression(posit8_t c, posit8_t r);
p32_FP_INT_P p16_compression(posit16_t c, posit16_t r);

void p8_read_power(posit8_t x, bool * sign, int8_t * k, uint8_t * reg, uint8_t * n_frac);
void p16_read_power(posit16_t x, bool * sign, int8_t * k, uint8_t * reg, int8_t * exp, uint8_t * n_frac);
void p32_read_power(posit32_t x, bool * sign, int8_t * k, uint8_t * reg, int8_t * exp, uint8_t * n_frac);

posit8_t p8_truncate(posit8_t x, int index);
posit16_t p16_truncate(posit16_t x, int index);
posit32_t p32_truncate(posit32_t x, int index);

#endif
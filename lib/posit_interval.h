#ifndef _POSIT_INTERVAL_H_
#define _POSIT_INTERVAL_H_

#include <stdint.h>
#include <stdio.h>

#include "softposit.h"

#define CONST_IPB_8 1.0 / (2 * 2 * 8)
#define CONST_IPB_16 1.0 / (2 * 2 * 16)

#define MAXPOS8 0x7F
#define MAXPOS16 0x7FFF
#define MAXPOS32 0x7FFFFFFF

#define MAXPOS64 0x7FFFFFFFFFFFFFFF

#define MINNEG8 0x81
#define MINNEG16 0x8001
#define MINNEG32 0x80000001

#define NAR8 0x80
#define NAR16 0x8000
#define NAR32 0x80000000

#define NAR64 0x8000000000000000

#define Fm1_16 convertDoubleToP16(-1.0)
#define F0_16 convertDoubleToP16(0.0)
#define F2_16 convertDoubleToP16(2.0)

#define Fm1_32 convertDoubleToP32(-1.0)
#define F0_32 convertDoubleToP32(0.0)
#define F2_32 convertDoubleToP32(2.0)


typedef struct {
    posit8_t center;
    posit8_t radius;
} C_R_p8;

typedef struct {
    posit16_t center;
    posit16_t radius;
} C_R_p16;

typedef struct {
    posit32_t center;
    posit32_t radius;
} C_R_p32;

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


// find next
posit8_t p8_nextup(posit8_t x);
posit16_t p16_nextup(posit16_t x);
posit32_t p32_nextup(posit32_t x);

posit8_t p8_nextdown(posit8_t x);
posit16_t p16_nextdown(posit16_t x);
posit32_t p32_nextdown(posit32_t x); 


// round up
posit8_t softposit_addMagsP8_up(uint_fast8_t uiA, uint_fast8_t uiB);
posit16_t softposit_addMagsP16_up(uint_fast16_t uiA, uint_fast16_t uiB);
posit32_t softposit_addMagsP32_up(uint_fast32_t uiA, uint_fast32_t uiB);

posit8_t softposit_subMagsP8_up(uint_fast8_t uiA, uint_fast8_t uiB);
posit16_t softposit_subMagsP16_up(uint_fast16_t uiA, uint_fast16_t uiB);
posit32_t softposit_subMagsP32_up(uint_fast32_t uiA, uint_fast32_t uiB);

posit8_t p8_add_up(posit8_t a, posit8_t b);
posit16_t p16_add_up(posit16_t a, posit16_t b);
posit32_t p32_add_up(posit32_t a, posit32_t b);

posit8_t p8_sub_up(posit8_t a, posit8_t b);
posit16_t p16_sub_up(posit16_t a, posit16_t b);
posit32_t p32_sub_up(posit32_t a, posit32_t b);

posit8_t p8_mul_up(posit8_t a, posit8_t b);
posit16_t p16_mul_up(posit16_t a, posit16_t b);
posit32_t p32_mul_up(posit32_t a, posit32_t b);

posit8_t p8_div_up(posit8_t a, posit8_t b);
posit16_t p16_div_up(posit16_t a, posit16_t b);
posit32_t p32_div_up(posit32_t a, posit32_t b);


// round down
posit8_t softposit_addMagsP8_down(uint_fast8_t uiA, uint_fast8_t uiB);
posit16_t softposit_addMagsP16_down(uint_fast16_t uiA, uint_fast16_t uiB);
posit32_t softposit_addMagsP32_down(uint_fast32_t uiA, uint_fast32_t uiB);

posit8_t softposit_subMagsP8_down(uint_fast8_t uiA, uint_fast8_t uiB);
posit16_t softposit_subMagsP16_down(uint_fast16_t uiA, uint_fast16_t uiB);
posit32_t softposit_subMagsP32_down(uint_fast32_t uiA, uint_fast32_t uiB);

posit8_t p8_add_down(posit8_t a, posit8_t b);
posit16_t p16_add_down(posit16_t a, posit16_t b);
posit32_t p32_add_down(posit32_t a, posit32_t b);

posit8_t p8_sub_down(posit8_t a, posit8_t b);
posit16_t p16_sub_down(posit16_t a, posit16_t b);
posit32_t p32_sub_down(posit32_t a, posit32_t b);

posit8_t p8_mul_down(posit8_t a, posit8_t b);
posit16_t p16_mul_down(posit16_t a, posit16_t b);
posit32_t p32_mul_down(posit32_t a, posit32_t b);

posit8_t p8_div_down(posit8_t a, posit8_t b);
posit16_t p16_div_down(posit16_t a, posit16_t b);
posit32_t p32_div_down(posit32_t a, posit32_t b);


// interval tools
C_R_p8 reciprocal_p8(C_R_p8 x);
C_R_p16 reciprocal_p16(C_R_p16 x);
C_R_p32 reciprocal_p32(C_R_p32 x);

C_R_p8 intersect_p8(C_R_p8 x, C_R_p8 y);
C_R_p16 intersect_p16(C_R_p16 x,C_R_p16 y);
C_R_p32 intersect_p32(C_R_p32 x, C_R_p32 y);

C_R_p8 interval_add_p8(C_R_p8 x, C_R_p8 y);
C_R_p16 interval_add_p16(C_R_p16 x, C_R_p16 y);
C_R_p32 interval_add_p32(C_R_p32 x, C_R_p32 y);

C_R_p8 interval_sub_p8(C_R_p8 x, C_R_p8 y);
C_R_p16 interval_sub_p16(C_R_p16 x, C_R_p16 y);
C_R_p32 interval_sub_p32(C_R_p32 x, C_R_p32 y);

C_R_p8 interval_mul_p8(C_R_p8 x, C_R_p8 y);
C_R_p16 interval_mul_p16(C_R_p16 x, C_R_p16 y);
C_R_p32 interval_mul_p32(C_R_p32 x, C_R_p32 y);

// tri-diagonal Gauss-Seidel 
void GS_tridiag_p8(int n, C_R_p8 *A, C_R_p8 *b, C_R_p8 *x);
void GS_tridiag_p16(int n, C_R_p16 *A, C_R_p16 *b, C_R_p16 *x);
void GS_tridiag_p32(int n, C_R_p32 *A, C_R_p32 *b, C_R_p32 *x);

// CSR sparse Gauss-Seidel
void GS_csr_p8(int n, C_R_p8 *A, int *idx, int *col_id, C_R_p8 *b, C_R_p8 *x);
void GS_csr_p16(int n, C_R_p16 *A, int *idx, int *col_id, C_R_p16 *b, C_R_p16 *x);
void GS_csr_p32(int n, C_R_p32 *A, int *idx, int *col_id, C_R_p32 *b, C_R_p32 *x);

#endif
#ifndef _CONVERT_H_
#define _CONVERT_H_

#define DOUBLE_ULS 1023
#define DOUBLE_E 52
#define DOUBLE_M 11

#define EPS set_pow2(-53)
#define REALMIN set_pow2(-1022)

#define ITERMAX 15
#define TOLERANCE set_pow2(-52)

#include <ieee754.h>
#include <time.h>
#include <sys/time.h>
#include <stdint.h>
#include <stdio.h>

typedef double FP_INT;

typedef struct {
    double inf;
    double sup;
} INF_SUP;

typedef struct {
    double center;
    double radius;
} C_R;

typedef struct {
    FP_INT center;
    int precision;
} FP_INT_PREC;

static inline int int_min(int a, int b) 
{
    return (a < b) ? a : b;
}

static inline int int_max(int a, int b) 
{
    return (a > b) ? a : b;
}

int read_m_bit(double x, int index);

double set_m_bit1(double x, int index);

double set_m_bit1_ieee(double x, int index);

double truncate_m(double x, int index);

double truncate_m_ieee(double x, int index);

double set_pow2(int exp); 

int get_sign_bit(double x); 

double set_sign_bit(double x, int sign);

FP_INT CR_FP1(C_R int_cr);
FP_INT CR_FP2(C_R int_cr);
FP_INT CR_FP3(C_R int_cr);
FP_INT CR_FP4(C_R int_cr);
FP_INT CR_FP5(C_R int_cr);
FP_INT CR_FP6(C_R int_cr);
FP_INT CR_FP1b(C_R int_cr);
FP_INT CR_FP1_adj(C_R int_cr);
FP_INT CR_FP1_adjbis(C_R int_cr);
FP_INT CR_FP1_adjter(C_R int_cr);

FP_INT_PREC CR_FP_p(C_R int_cr);
FP_INT_PREC CR_FP_mpfr(C_R int_cr);

FP_INT IS_FP1(INF_SUP int_is);
FP_INT IS_FP2(INF_SUP int_is);
FP_INT IS_FP3(INF_SUP int_is);
FP_INT IS_FP4(INF_SUP int_is);
FP_INT IS_FP5(INF_SUP int_is);
FP_INT IS_FP6(INF_SUP int_is);

C_R FP_CR(FP_INT c_tilde);

C_R FP_CR3(FP_INT c_tilde);

void print_binary(double x);

void fprint_binary(FILE *fp, double x);

void print_binary_u32(uint32_t x);

void fprint_binary_u32(FILE *fp, uint32_t x);

void print_binary_u8(uint8_t x);

void fprint_binary_u8(FILE *fp, uint8_t x);

typedef union {
    double d;
    uint64_t u64;
} cast;

typedef struct D_F {
    double *d;
    float *f;
} D_F;

typedef struct MP {
    double *d;
    uint32_t *u32;
    uint16_t *u16;
    uint8_t *u8;
} MP;

// mask: 0111 precsion (36, 44], 0110 precision (28, 36], 0101 precision (20, 28], 
// 0100 precision (12, 20], 0011 precision (4, 12], 0010 precision [0, 4], default 0000
static inline uint8_t precision_to_mask3(int prec) 
{
    if (prec <= 4)       return 0b0010;
    else if (prec <= 12) return 0b0011;
    else if (prec <= 20) return 0b0100;
    else if (prec <= 28) return 0b0101;
    else if (prec <= 36) return 0b0110;
    else if (prec <= 44) return 0b0111;
    else                 return 0b1000;
}

double get_time_ms(); 

uint32_t * FP_u32(FP_INT_PREC * c_prec, int nnz, uint8_t * mask, int n_mask, int * count32);
D_F * FP_fd(FP_INT_PREC * c_prec, int nnz, uint8_t * mask, int n_mask, int * count_f);
MP * FP_mixed(FP_INT_PREC * c_prec, int nnz, uint8_t * mask, int n_mask, int * count_d, int * count32, int * count16, int * count8);

FP_INT * u32_FP(uint32_t * u32, int nnz, uint8_t * mask, int n_mask);
FP_INT * fd_FP(D_F * d_f, int nnz, uint8_t * mask, int n_mask);
FP_INT * mixed_FP(MP * mixed, int nnz, uint8_t * mask, int n_mask);

#endif
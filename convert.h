#ifndef _CONVERT_H_
#define _CONVERT_H_

#define DOUBLE_ULS 1023
#define DOUBLE_E 52
#define DOUBLE_M 11

#include <ieee754.h>

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

void print_binary(double x);


#endif
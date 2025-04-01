#ifndef _CONVERT_H_
#define _CONVERT_H_

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

static inline int int_min(int a, int b) 
{
    return (a < b) ? a : b;
}

int read_m_bit(double x, int index);

double set_m_bit(double x, int index, int value);

double truncate_m(double x, int index);

double set_pow2(int exp); 

FP_INT CR_FP1(C_R int_cr);

FP_INT CR_FP2(C_R int_cr);

FP_INT IS_FP(INF_SUP int_is);

C_R FP_CR(FP_INT c_tilde);

#endif
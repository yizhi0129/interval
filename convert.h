#ifndef _CONVERT_H_
#define _CONVERT_H_


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

int get_sign_bit(double x); 

FP_INT CR_FP1(C_R int_cr);
FP_INT CR_FP2(C_R int_cr);
FP_INT CR_FP3(C_R int_cr);
FP_INT CR_FP4(C_R int_cr);
FP_INT CR_FP5(C_R int_cr);

FP_INT IS_FP1(INF_SUP int_is);
FP_INT IS_FP2(INF_SUP int_is);
FP_INT IS_FP3(INF_SUP int_is);
FP_INT IS_FP4(INF_SUP int_is);
FP_INT IS_FP5(INF_SUP int_is);

C_R FP_CR(FP_INT c_tilde);

#endif
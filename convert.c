#pragma STDC FENV_ACCESS ON

#include "convert.h"
#include <stdbool.h>
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <fenv.h>

#define DOUBLE_ULS 1023
#define DOUBLE_E 52
#define DOUBLE_M 11


int read_m_bit(double x, int index)
{
    uint64_t bits;
    memcpy(&bits, &x, sizeof(double)); 
    return (bits >> index) & 1;
}

double set_m_bit(double x, int index, int value)
{
    union {
        double d;
        uint64_t u;
    } u = { .d = x };
    if (value)
        u.u |= (1ULL << index);       
    else
        u.u &= ~(1ULL << index);     
    return u.d;
}


double truncate_m(double x, int index) 
{
    union {
        double d;
        uint64_t u;
    } u = { .d = x };
    uint64_t mantissa_mask = ~((1ULL << index) - 1);
    mantissa_mask &= 0x000FFFFFFFFFFFFFULL;
    u.u &= 0xFFF0000000000000ULL | (mantissa_mask & 0x000FFFFFFFFFFFFFULL); 
    return u.d;
}

double set_pow2(int exp) 
{
    union {
        double d;
        uint64_t u;
    } val;

    int biased_exp = exp + DOUBLE_ULS;  
    val.u = (uint64_t)biased_exp << DOUBLE_E;  

    return val.d;
}

// pure bitwise operations
FP_INT CR_FP1(C_R int_cr)
{
    union ieee754_double u1, u2;
    u1.d = int_cr.center;
    u2.d = int_cr.radius;
    int e_c = u1.ieee.exponent - DOUBLE_ULS;
    int e_r = u2.ieee.exponent - DOUBLE_ULS;
    int p = int_min(e_c - e_r - 1, DOUBLE_E);
    FP_INT c_tilde = int_cr.center;
    double r_tilde = 0.0;
    while (true)
    {
        if (p < 1)
        {
            printf("1 Error: No enough precision to convert!\n");
            printf("c = %lf, r = %lf\n", int_cr.center, int_cr.radius);
            return 1;
        }
        int index = DOUBLE_E - p;
        c_tilde = truncate_m(c_tilde, index);
        c_tilde = set_m_bit(c_tilde, index, 1);
        r_tilde = set_pow2(e_c - p);

        //compare
        int b = read_m_bit(int_cr.center, index);
        if ((2*b-1)*int_cr.center + int_cr.radius <= (2*b-1)*c_tilde + r_tilde) 
        {
            break;
        }
        else
        {
            p --;
        }
    }     
    return c_tilde;
}


// floating point operations
FP_INT CR_FP2(C_R int_cr)
{
    double c = int_cr.center;
    double r = int_cr.radius;
    int e_c = (int)floor(log2(fabs(c)));
    int e_r = (int)floor(log2(r));
    int p = int_min(e_c - e_r - 1, DOUBLE_E);
    FP_INT c_tilde = int_cr.center;
    double r_tilde = 0.0;
    while (true)
    {
        if (p < 1)
        {
            printf("2 Error: No enough precision to convert!\n");
            printf("c = %lf, r = %lf\n", c, r);
            return 1;
        }
        int index = DOUBLE_E - p;
        c_tilde = truncate_m(c, index);
        c_tilde = set_m_bit(c_tilde, index, 1);
        r_tilde = set_pow2(e_c - p);

        // compare
        fesetround(FE_DOWNWARD);
        double inf_new = c_tilde - r_tilde;
        fesetround(FE_UPWARD);
        double sup_new = c_tilde + r_tilde;
        double a = inf_new + r;
        double b = c + r;
        fesetround(FE_TONEAREST);
        if (a <= c && b <= sup_new) 
        {
            break;
        }
        else
        {
            p --;
        }
    } 
    return c_tilde;
}


FP_INT CR_FP3(C_R int_cr)
{
    union ieee754_double u1, u2;
    u1.d = int_cr.center;
    u2.d = int_cr.radius;
    int e_c = u1.ieee.exponent - DOUBLE_ULS;
    int e_r = u2.ieee.exponent - DOUBLE_ULS;
    int p = int_min(e_c - e_r - 1, DOUBLE_E);
    double c = int_cr.center;
    double r = int_cr.radius;
    FP_INT c_tilde = c;
    double r_tilde = 0.0;
    while (true)
    {
        if (p < 1)
        {
            printf("3 Error: No enough precision to convert!\n");
            printf("c = %lf, r = %lf\n", c, r);
            return 1;
        }
        int index = DOUBLE_E - p;
        c_tilde = truncate_m(c_tilde, index);
        c_tilde = set_m_bit(c_tilde, index, 1);
        r_tilde = set_pow2(e_c - p);

        fesetround(FE_DOWNWARD);
        double inf_new = c_tilde - r_tilde;
        fesetround(FE_UPWARD);
        double sup_new = c_tilde + r_tilde;
        double a = inf_new + r;
        double b = c + r;
        fesetround(FE_TONEAREST);
        if (a <= c && b <= sup_new) 
        {
            break;
        }
        else
        {
            p --;
        }
    }     
    return c_tilde;
}


FP_INT CR_FP4(C_R int_cr)
{
    double c = int_cr.center;
    double r = int_cr.radius;
    int e_c = (int)floor(log2(fabs(c)));
    int e_r = (int)floor(log2(r));
    int p = int_min(e_c - e_r - 1, DOUBLE_E);
    FP_INT c_tilde = int_cr.center;
    double r_tilde = 0.0;
    bool cond = true;
    while (cond)
    {
        if (p < 1)
        {
            printf("4 Error: No enough precision to convert!\n");
            printf("c = %lf, r = %lf\n", c, r);
            return 1;
        }
        int index = DOUBLE_E - p;
        c_tilde = truncate_m(c, index);
        c_tilde = set_m_bit(c_tilde, index, 1);
        r_tilde = set_pow2(e_c - p);

        int b = read_m_bit(c, index);
        if ((2*b-1) * c + r <= (2*b-1) * c_tilde + r_tilde) 
        {
            cond = false;
            break;
        }
        else
        {
            p --;
        }
    } 
    return c_tilde;
}


FP_INT IS_FP1(INF_SUP int_is)
{
    double x = int_is.inf;
    double y = int_is.sup;
    C_R int_cr;
    fesetround(FE_UPWARD);
    int_cr.center = x + 0.5 * (y - x);
    int_cr.radius = int_cr.center - x;
    fesetround(FE_TONEAREST);
    return CR_FP1(int_cr);
}

FP_INT IS_FP2(INF_SUP int_is)
{
    double x = int_is.inf;
    double y = int_is.sup;
    C_R int_cr;
    fesetround(FE_UPWARD);
    int_cr.center = x + 0.5 * (y - x);
    int_cr.radius = int_cr.center - x;
    fesetround(FE_TONEAREST);
    return CR_FP2(int_cr);
}

FP_INT IS_FP3(INF_SUP int_is)
{
    double x = int_is.inf;
    double y = int_is.sup;
    C_R int_cr;
    fesetround(FE_UPWARD);
    int_cr.center = x + 0.5 * (y - x);
    int_cr.radius = int_cr.center - x;
    fesetround(FE_TONEAREST);
    return CR_FP3(int_cr);
}

FP_INT IS_FP4(INF_SUP int_is)
{
    double x = int_is.inf;
    double y = int_is.sup;
    C_R int_cr;
    fesetround(FE_UPWARD);
    int_cr.center = x + 0.5 * (y - x);
    int_cr.radius = int_cr.center - x;
    fesetround(FE_TONEAREST);
    return CR_FP4(int_cr);
}

C_R FP_CR(FP_INT c_tilde)
{
    union ieee754_double u;
    u.d = c_tilde;
    int e_c = u.ieee.exponent - DOUBLE_ULS;
    int index = 0;
    for (int i = 0; i < DOUBLE_E; i ++)
    {
        int b = read_m_bit(c_tilde, i);
        if (b == 1)
        {
            index = i;
            break;
        }
    }
    int p = DOUBLE_E - index;
    double r_tilde = set_pow2(e_c - p);
    C_R result = {c_tilde, r_tilde};
    return result;
}
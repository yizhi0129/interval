#include "convert.h"
#include <stdbool.h>
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <fenv.h>


int get_sign_bit(double x) 
{
    union {
        double d;
        uint64_t u;
    } u = { .d = x };

    return (u.u >> 63) & 1;  
}


double set_sign_bit(double x, int sign) 
{
    uint64_t* bits = (uint64_t*)&x;
    *bits &= 0x7FFFFFFFFFFFFFFF;      
    *bits |= ((uint64_t)sign << 63);           
    return x;
}


// check valid zone
int read_m_bit(double x, int index)
{
    uint64_t bits;
    memcpy(&bits, &x, sizeof(double)); 
    return (bits >> index) & 1;
}

double set_m_bit1(double x, int index)
{
    union {
        double d;
        uint64_t u;
    } u = { .d = x };
    u.u |= (1ULL << index);         
    return u.d;
}

double set_m_bit1_ieee(double x, int index)
{
    union ieee754_double u;
    u.d = x;

    if (index < 0 || index >= 52) return x;

    if (index < 32) 
    {    
        u.ieee.mantissa1 |= (1U << index); 
    } 
    else 
    {
        int shifted = index - 32;
        u.ieee.mantissa0 |= (1U << shifted); 
    }
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

double truncate_m_ieee(double x, int index)
{
    union ieee754_double u;
    u.d = x;
    if (index == 0) 
    {
        u.ieee.mantissa0 = 0;
        u.ieee.mantissa1 = 0;
    } 
    else if (index < 32) 
    {
        uint32_t mask = ~((1U << index) - 1); 
        u.ieee.mantissa1 &= mask;
    } 
    else 
    {
        int hi_bits = index - 32;
        uint32_t mask = ~((1U << hi_bits) - 1); 
        u.ieee.mantissa0 &= mask;
        u.ieee.mantissa1 = 0;
    }
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
    int p = int_min(e_c - e_r, DOUBLE_E);
    FP_INT c_tilde = int_cr.center;
    double r_tilde = 0.0;
    int s = get_sign_bit(int_cr.center);
    int b = 1 ^ s;
    double u = set_pow2(e_c - 52);
    while (true)
    {
        if (p > 0)
        {
            int index = DOUBLE_E - p;
            c_tilde = truncate_m(c_tilde, index);
            c_tilde = set_m_bit1(c_tilde, index);
            r_tilde = set_pow2(e_c - p);
            b = read_m_bit(int_cr.center, index) ^ s;          
        }
        else
        {
            c_tilde = set_pow2(e_c - p);
            c_tilde = set_sign_bit(c_tilde, s);
            r_tilde = set_pow2(e_c - p);
        }
        
        //compare      
        if ((2*b-1)*int_cr.center + int_cr.radius + u <= (2*b-1)*c_tilde + r_tilde) 
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
    int p = int_min(e_c - e_r, DOUBLE_E);
    FP_INT c_tilde = int_cr.center;
    double r_tilde = 0.0;
    int s = get_sign_bit(c);
    while (true)
    {
        if (p > 0)
        {
            int index = DOUBLE_E - p;
            c_tilde = truncate_m(c, index);
            c_tilde = set_m_bit1(c_tilde, index);
            r_tilde = set_pow2(e_c - p);
        }
        else
        {
            c_tilde = set_pow2(e_c - p);
            c_tilde = set_sign_bit(c_tilde, s);
            r_tilde = set_pow2(e_c - p);
        }

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
    int p = int_min(e_c - e_r, DOUBLE_E);
    double c = int_cr.center;
    double r = int_cr.radius;
    FP_INT c_tilde = c;
    double r_tilde = 0.0;
    int s = get_sign_bit(c);
    while (true)
    {
        if (p > 0)
        {
            int index = DOUBLE_E - p;
            c_tilde = truncate_m(c_tilde, index);
            c_tilde = set_m_bit1(c_tilde, index);
            r_tilde = set_pow2(e_c - p);
        }
        else
        {
            c_tilde = set_pow2(e_c - p);
            c_tilde = set_sign_bit(c_tilde, s);
            r_tilde = set_pow2(e_c - p);
        }

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
    int p = int_min(e_c - e_r, DOUBLE_E);
    FP_INT c_tilde = int_cr.center;
    double r_tilde = 0.0;
    int s = get_sign_bit(c);
    int b = 1 ^ s;
    double u = set_pow2(e_c - 52);
    while (true)
    {
        if (p > 0)
        {
            int index = DOUBLE_E - p;
            c_tilde = truncate_m(c_tilde, index);
            c_tilde = set_m_bit1(c_tilde, index);
            r_tilde = set_pow2(e_c - p);
            b = read_m_bit(c, index) ^ s;
        }
        else
        {
            c_tilde = set_pow2(e_c - p);
            c_tilde = set_sign_bit(c_tilde, s);
            r_tilde = set_pow2(e_c - p);
        }

        if ((2*b-1) * c + r + u <= (2*b-1) * c_tilde + r_tilde) 
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

FP_INT CR_FP5(C_R int_cr)
{
    double c = int_cr.center;
    double r = int_cr.radius;
    int e_c, e_r;
    frexp(c, &e_c);
    frexp(r, &e_r); 
    e_c -= 1;
    e_r -= 1;
    int p = int_min(e_c - e_r, DOUBLE_E);
    FP_INT c_tilde = int_cr.center;
    double r_tilde = 0.0;
    int s = get_sign_bit(c);
    int b = 1 ^ s;
    double u = set_pow2(e_c - 52);
    while (true)
    {
        if (p > 0)
        {
            int index = DOUBLE_E - p;
            c_tilde = truncate_m(c_tilde, index);
            c_tilde = set_m_bit1(c_tilde, index);
            r_tilde = set_pow2(e_c - p);
            b = read_m_bit(c, index) ^ s;
        }
        else
        {
            c_tilde = set_pow2(e_c - p);
            c_tilde = set_sign_bit(c_tilde, s);
            r_tilde = set_pow2(e_c - p);
        }

        if ((2*b-1) * c + r + u <= (2*b-1) * c_tilde + r_tilde) 
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

FP_INT CR_FP6(C_R int_cr)
{
    union ieee754_double u1, u2;
    u1.d = int_cr.center;
    u2.d = int_cr.radius;
    int e_c = u1.ieee.exponent - DOUBLE_ULS;
    int e_r = u2.ieee.exponent - DOUBLE_ULS;
    int p = int_min(e_c - e_r, DOUBLE_E);
    FP_INT c_tilde = int_cr.center;
    double r_tilde = 0.0;
    int s = get_sign_bit(int_cr.center);
    int b = 1 ^ s;
    double u = set_pow2(e_c - 52);
    while (true)
    {
        if (p > 0)
        {
            int index = DOUBLE_E - p;
            c_tilde = truncate_m_ieee(c_tilde, index);
            c_tilde = set_m_bit1_ieee(c_tilde, index);
            r_tilde = set_pow2(e_c - p);
            b = read_m_bit(int_cr.center, index) ^ s;
        }
        else
        {
            c_tilde = set_pow2(e_c - p);
            c_tilde = set_sign_bit(c_tilde, s);
            r_tilde = set_pow2(e_c - p);
        }

        if ((2*b-1) * int_cr.center + int_cr.radius + u <= (2*b-1) * c_tilde + r_tilde) 
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

FP_INT CR_FP1b(C_R int_cr)
{
    union ieee754_double u1, u2;
    u1.d = int_cr.center;
    u2.d = int_cr.radius;
    int e_c = u1.ieee.exponent - DOUBLE_ULS;
    int e_r = u2.ieee.exponent - DOUBLE_ULS;
    int p = int_min(e_c - e_r, DOUBLE_E);
    FP_INT c_tilde = int_cr.center;
    double r_tilde = 0.0;
    int a = read_m_bit(int_cr.radius, DOUBLE_E + 1 - e_c + e_r) ^ 1;
    int s = get_sign_bit(int_cr.center);
    int b = 1 ^ s;
    double u = set_pow2(e_c - 52);
    while (true)
    {
        if (p > 0)
        {
            int index = DOUBLE_E - p;
            c_tilde = truncate_m(c_tilde, index);
            c_tilde = set_m_bit1(c_tilde, index);
            r_tilde = set_pow2(e_c - p);
            b = read_m_bit(int_cr.center, index) ^ s;          
        }
        else
        {
            c_tilde = set_pow2(e_c - p);
            c_tilde = set_sign_bit(c_tilde, s);
            r_tilde = set_pow2(e_c - p);
        }
        
        //compare      
        if ((2*b-1)*int_cr.center + nextafter(int_cr.radius, INFINITY) + a * u <= (2*b-1)*c_tilde + r_tilde) 
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


FP_INT CR_FP1_adj(C_R int_cr)
{
    union ieee754_double u1, u2;
    u1.d = int_cr.center;
    u2.d = int_cr.radius;
    int e_c = u1.ieee.exponent - DOUBLE_ULS;
    int e_r = u2.ieee.exponent - DOUBLE_ULS;
    int p = int_min(e_c - e_r, DOUBLE_E);
    FP_INT c_tilde = int_cr.center;
    double r_tilde = 0.0;
    int s = get_sign_bit(int_cr.center);
    int b = 1 ^ s;
    double u = set_pow2(e_c - 52);
    while (true)
    {
        if (p > 0)
        {
            int index = DOUBLE_E - p;
            c_tilde = truncate_m(c_tilde, index);
            c_tilde = set_m_bit1(c_tilde, index);
            r_tilde = set_pow2(e_c - p);
            b = read_m_bit(int_cr.center, index) ^ s;          
        }
        else
        {
            c_tilde = set_pow2(e_c - p);
            c_tilde = set_sign_bit(c_tilde, s);
            r_tilde = set_pow2(e_c - p);
        }
        
        //compare      
        if ((2*b-1)*int_cr.center + int_cr.radius + u <= (2*b-1)*c_tilde + 2 * r_tilde) 
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

FP_INT IS_FP5(INF_SUP int_is)
{
    double x = int_is.inf;
    double y = int_is.sup;
    C_R int_cr;
    fesetround(FE_UPWARD);
    int_cr.center = x + 0.5 * (y - x);
    int_cr.radius = int_cr.center - x;
    fesetround(FE_TONEAREST);
    return CR_FP5(int_cr);
}

FP_INT IS_FP6(INF_SUP int_is)
{
    double x = int_is.inf;
    double y = int_is.sup;
    C_R int_cr;
    fesetround(FE_UPWARD);
    int_cr.center = x + 0.5 * (y - x);
    int_cr.radius = int_cr.center - x;
    fesetround(FE_TONEAREST);
    return CR_FP6(int_cr);
}

C_R FP_CR(FP_INT c_tilde)
{
    union ieee754_double u;
    u.d = c_tilde;
    int e_c = u.ieee.exponent - DOUBLE_ULS;
    int p = 0;

    for (int i = 0; i <= DOUBLE_E; i++)
    {
        if (read_m_bit(c_tilde, i) == 1)
        {
            p = DOUBLE_E - i;
            break;
        }
    }

    double r_tilde = set_pow2(e_c - p);
    C_R result = {c_tilde, r_tilde};
    return result;
}

void print_binary(double x)
{
    union ieee754_double u;
    u.d = x;
    int exp = u.ieee.exponent - DOUBLE_ULS;
    int sign = get_sign_bit(x);
    (sign == 1) ? printf("- ") : printf("+ ");
    printf("1.");
    for (int i = DOUBLE_E - 1; i >= 0; i --)
    {
        int b = read_m_bit(x, i);
        printf("%d", b);
    }
    printf(" 2^%d", exp);
}

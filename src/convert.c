#include "convert.h"
#include <stdbool.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
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


FP_INT CR_FP1_adjbis(C_R int_cr)
{
    union ieee754_double u1, u2;
    u1.d = int_cr.center;
    u2.d = int_cr.radius;
    int e_c = u1.ieee.exponent - DOUBLE_ULS;
    int e_r = u2.ieee.exponent - DOUBLE_ULS;
    int p = int_min(e_c - e_r, DOUBLE_E);
    FP_INT c_tilde = int_cr.center;
    int s = get_sign_bit(c_tilde);
    if (p > 0)
    {
        int index = DOUBLE_E - p;
        c_tilde = truncate_m(c_tilde, index);
        c_tilde = set_m_bit1(c_tilde, index);        
    }
    else
    {
        c_tilde = set_pow2(e_c - p);
        c_tilde = set_sign_bit(c_tilde, s);        
    }     
    return c_tilde;
}


FP_INT CR_FP1_adjter(C_R int_cr)
{
    union ieee754_double u1, u2;
    u1.d = int_cr.center;
    u2.d = int_cr.radius;
    int e_c = u1.ieee.exponent - DOUBLE_ULS;
    int e_r = u2.ieee.exponent - DOUBLE_ULS;
    int p = int_min(e_c - e_r, DOUBLE_E);
    FP_INT c_tilde = int_cr.center;
    int s = get_sign_bit(c_tilde);
    if (p > 0)
    {
        int index = DOUBLE_E - p;
        c_tilde = truncate_m(c_tilde, index);
        c_tilde = set_m_bit1(c_tilde, index);        
    }
    else
    {
        c_tilde = set_pow2(e_c - p);
        c_tilde = set_sign_bit(c_tilde, s);        
    }     
    return c_tilde;
}


FP_INT_PREC CR_FP_p(C_R int_cr)
{
    union ieee754_double u1, u2;
    u1.d = int_cr.center;
    u2.d = int_cr.radius;
    int e_c = u1.ieee.exponent - DOUBLE_ULS;
    int e_r = u2.ieee.exponent - DOUBLE_ULS;
    int p = int_min(e_c - e_r, DOUBLE_E);
    FP_INT c_tilde = int_cr.center;
    int s = get_sign_bit(c_tilde);
    if (p > 0)
    {
        int index = DOUBLE_E - p;
        c_tilde = truncate_m(c_tilde, index);
        c_tilde = set_m_bit1(c_tilde, index);        
    }
    else
    {
        c_tilde = set_pow2(e_c - p);
        c_tilde = set_sign_bit(c_tilde, s);        
    }     
    FP_INT_PREC result = {c_tilde, p};
    return result;
}

FP_INT_PREC CR_FP_mpfr(C_R int_cr)
{
    union ieee754_double u1, u2;
    u1.d = int_cr.center;
    u2.d = int_cr.radius;
    int e_c = u1.ieee.exponent - DOUBLE_ULS;
    int e_r = u2.ieee.exponent - DOUBLE_ULS;
    int p = int_min(e_c - e_r, DOUBLE_E);
    FP_INT c_tilde = int_cr.center;
    int s = get_sign_bit(c_tilde);
    if (p > 0)
    {
        int index = DOUBLE_E - p;
        // c_tilde = truncate_m(c_tilde, index);
        c_tilde = set_m_bit1(c_tilde, index);        
    }
    else
    {
        c_tilde = set_pow2(e_c - p);
        c_tilde = set_sign_bit(c_tilde, s);        
    } 
    int prec = int_max(p + 1, 2);    
    FP_INT_PREC result = {c_tilde, prec};
    return result;
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

// read r_tiled = uls(c_tilde)
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

// read r_tiled = 3 * uls(c_tilde)
C_R FP_CR3(FP_INT c_tilde)
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

    double r_tilde = 3.0 * set_pow2(e_c - p);
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
    printf(" 2^%d\n", exp);
}

void fprint_binary(FILE *fp, double x)
{
    union ieee754_double u;
    u.d = x;
    int exp = u.ieee.exponent - DOUBLE_ULS;
    int sign = get_sign_bit(x);
    (sign == 1) ? fprintf(fp, "- ") : fprintf(fp, "+ ");
    fprintf(fp, "1.");
    for (int i = DOUBLE_E - 1; i >= 0; i --)
    {
        int b = read_m_bit(x, i);
        fprintf(fp, "%d", b);
    }
    fprintf(fp, " 2^%d\n", exp);
}

void print_binary_u32(uint32_t x)
{
    for (int i = 31; i >= 0; i --)
    {
        int b = (x >> i) & 1;
        printf("%d", b);
    }
}

void fprint_binary_u32(FILE *fp, uint32_t x)
{
    for (int i = 31; i >= 0; i --)
    {
        int b = (x >> i) & 1;
        fprintf(fp, "%d", b);
    }
}

void print_binary_u8(uint8_t x)
{
    for (int i = 7; i >= 0; i --)
    {
        int b = (x >> i) & 1;
        printf("%d", b);
    }
}

void fprint_binary_u8(FILE *fp, uint8_t x)
{
    for (int i = 7; i >= 0; i --)
    {
        int b = (x >> i) & 1;
        fprintf(fp, "%d", b);
    }
}

double get_time_ms() 
{
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return (double)tv.tv_sec * 1000.0 + (double)tv.tv_usec / 1000.0;
}


// convert FP_INT to uint32_t
uint32_t * FP_u32(FP_INT_PREC * c_prec, int nnz, uint8_t * mask, int n_mask, int * count32)
{
    int low_count = 0;
    for (int i = 0; i < n_mask; i ++)
    {
        // set mask bits
        mask[i] = 0;
        for (int j = 0; j < 8; j ++)
        {
            int idx = i * 8 + j;
            if (idx >= nnz) break;

            if (c_prec[idx].precision <= 20)
            {
                mask[i] |= (1 << (7 - j));
                low_count ++;
            }
        }
    }


    (*count32) = 2 * nnz - low_count;
    uint32_t *result = malloc((*count32) * sizeof(uint32_t));
    if (!result) return NULL;

    int k = 0;  // res32 index

    for (int i = 0; i < n_mask; i ++) 
    {
        for (int j = 0; j < 8; j ++) 
        {
            int idx = i * 8 + j;
            if (idx >= nnz) break;

            cast dc = { .u64 = 0};
            dc.d = c_prec[idx].center;

            if ((mask[i] >> (7 - j)) & 1) 
            {
                result[k ++] = (uint32_t)(dc.u64 >> 32);
            } 
            else 
            {
                result[k ++] = (uint32_t)(dc.u64 >> 32);       // high part
                result[k ++] = (uint32_t)(dc.u64 & 0xFFFFFFFF); // low part
            }
        } 
    }
    return result;
}

// convert FP_INT to float/double
D_F * FP_fd(FP_INT_PREC * c_prec, int nnz, uint8_t * mask, int n_mask, int * count_f)
{
    for (int i = 0; i < n_mask; i ++)
    {
        // set mask bits
        mask[i] = 0;
        for (int j = 0; j < 8; j ++)
        {
            int idx = i * 8 + j;
            if (idx >= nnz) break;

            if (c_prec[idx].precision <= 23)
            {
                mask[i] |= (1 << j);
            }
        }
    }

    // count number of float
    *count_f = 0;
    for (int i = 0; i < n_mask; i ++)
    {
        for (int j = 0; j < 8; j ++) 
        {
            int idx = i * 8 + j;
            if (idx >= nnz) break;
            *count_f += (mask[i] >> j) & 1;
        }
    }

    D_F *result = malloc(sizeof(D_F));
    result->d = malloc((nnz - *count_f) * sizeof(double));
    result->f = malloc(*count_f * sizeof(float));
    
    if (!result->d || !result->f) return NULL;

    int id_f = 0; // index for float array
    int id_d = 0; // index for double array

    for (int i = 0; i < n_mask; i ++) 
    {
        for (int j = 7; j >= 0; j --) 
        {
            int idx = i * 8 + j;
            if (idx >= nnz) continue;

            if ((mask[i] >> j) & 1) 
            {
                result->f[id_f ++] = (float)c_prec[idx].center;
            } 
            else 
            {
                result->d[id_d ++] = c_prec[idx].center;
            }
        } 
    }   
    return result;
}

// convert FP_INT to mixed precision
MP * FP_mixed(FP_INT_PREC * c_prec, int nnz, uint8_t * mask, int n_mask, int * count_d, int * count32, int * count16, int * count8)
{
    // set mask bits
    for (int i = 0; i < n_mask; i ++)
    {
        int idx = i * 2;
        int idx2 = i * 2 + 1;
        uint8_t m1 = precision_to_mask3(c_prec[idx].precision);
        uint8_t m2 = 0;
        if (idx2 < nnz) 
        {
            m2 = precision_to_mask3(c_prec[idx2].precision);
        }
        mask[i] = (m1 << 4) | m2;
    }

    // count number of double, uint32_t, uint16_t, uint8_t
    *count_d = 0;
    *count32 = 0;
    *count16 = 0;
    *count8 = 0;
    for (int i = 0; i < n_mask; i ++)
    {
        uint8_t m = mask[i];
        *count8 += m & 1;
        *count16 += (m >> 1) & 1;
        *count32 += (m >> 2) & 1;
        *count_d += (m >> 3) & 1;
        *count8 += (m >> 4) & 1;
        *count16 += (m >> 5) & 1;
        *count32 += (m >> 6) & 1;
        *count_d += (m >> 7) & 1;
    }
    MP *result = malloc(sizeof(MP));
    if (!result) return NULL;

    result->d = malloc(*count_d * sizeof(double));
    result->u32 = malloc(*count32 * sizeof(uint32_t));
    result->u16 = malloc(*count16 * sizeof(uint16_t));
    result->u8 = malloc(*count8 * sizeof(uint8_t));
    if (!result->d || !result->u32 || !result->u16 || !result->u8) 
    {
        free(result->d);
        free(result->u32);
        free(result->u16);
        free(result->u8);
        free(result);
        return NULL;
    }

    int id_d = 0; // index for double array
    int id_u32 = 0; // index for uint32_t array
    int id_u16 = 0; // index for uint16_t array
    int id_u8 = 0; // index for uint8_t array

    for (int i = 0; i < n_mask; i ++) 
    {
        uint8_t m2 = mask[i] & 0x0F;
        uint8_t m1 = mask[i] >> 4;
        int ind1 = i * 2;
        int ind2 = i * 2 + 1;
        cast dc;

        dc.d = c_prec[ind1].center;

        if (m1 == 0x08) 
        {
            result->d[id_d ++] = dc.d;
        }  
        else 
        {
            if (m1 & 0x04) 
            {  // u32
                result->u32[id_u32 ++] = (uint32_t)(dc.u64 >> 32);
                if (m1 & 0x02) 
                {  // u16
                    result->u16[id_u16 ++] = (uint16_t)(dc.u64 >> 16);
                    if (m1 & 0x01) 
                    {  // u8
                        result->u8[id_u8 ++] = (uint8_t)(dc.u64 >> 8);
                    }
                } 
                else if (m1 & 0x01) 
                {  // u32 + u8
                    result->u8[id_u8 ++] = (uint8_t)(dc.u64 >> 24);
                }
            } 
            else if (m1 & 0x02) 
            {  // u16 only
                result->u16[id_u16 ++] = (uint16_t)(dc.u64 >> 48);
                if (m1 & 0x01) 
                {  // u16 + u8
                    result->u8[id_u8 ++] = (uint8_t)(dc.u64 >> 40);
                }
            } 
            else if (m1 & 0x01) 
            {  // u8 only
                result->u8[id_u8 ++] = (uint8_t)(dc.u64 >> 56);
            }
        }

        if (ind2 < nnz)
        {
            dc.d = c_prec[ind2].center;

            if (m2 == 0x08) 
            {
                result->d[id_d ++] = dc.d;
            } 
            else 
            {
                if (m2 & 0x04) 
                {
                    result->u32[id_u32 ++] = (uint32_t)(dc.u64 >> 32);
                    if (m2 & 0x02) 
                    {
                        result->u16[id_u16 ++] = (uint16_t)(dc.u64 >> 16);
                        if (m2 & 0x01) 
                        {
                            result->u8[id_u8 ++] = (uint8_t)(dc.u64 >> 8);
                        }
                    } 
                    else if (m2 & 0x01) 
                    {
                        result->u8[id_u8 ++] = (uint8_t)(dc.u64 >> 24);
                    }
                }    
                else if (m2 & 0x02) 
                {
                    result->u16[id_u16 ++] = (uint16_t)(dc.u64 >> 48);
                    if (m2 & 0x01) 
                    {
                        result->u8[id_u8 ++] = (uint8_t)(dc.u64 >> 40);
                    }
                } 
                else if (m2 & 0x01) 
                {
                    result->u8[id_u8++] = (uint8_t)(dc.u64 >> 56);
                }
            }
        }
    }
    return result;
}

// convert uint32_t to FP_INT
FP_INT * u32_FP(uint32_t * u32, int nnz, uint8_t * mask, int n_mask)
{
    FP_INT * result = malloc(nnz * sizeof(FP_INT));
    if (!result) return NULL;

    int k = 0;

    for (int i = 0; i < n_mask; i ++) 
    {
        for (int j = 0; j < 8; j ++) 
        {
            int idx = i * 8 + j;
            if (idx >= nnz) break;

            uint64_t high = (uint64_t)u32[k ++] << 32;
            uint64_t low = 0;

            if (!((mask[i] >> (7 - j)) & 1)) 
            {
                low = (uint64_t)u32[k ++];
            }

            uint64_t combined = high | low;

            cast dc = { .u64 = combined };

            result[idx] = dc.d;
        }
    }
    return result;
}

// convert D_F to FP_INT
FP_INT * fd_FP(D_F * d_f, int nnz, uint8_t * mask, int n_mask)
{
    FP_INT * result = malloc(nnz * sizeof(FP_INT));
    if (!result) return NULL;

    int k_d = 0; // index for double array
    int k_f = 0; // index for float array

    for (int i = 0; i < n_mask; i ++) 
    {
        for (int j = 7; j >= 0; j --) 
        {
            int idx = i * 8 + j;
            if (idx >= nnz) continue;

            if ((mask[i] >> j) & 1) 
            {
                result[idx] = (FP_INT)d_f->f[k_f ++];
            } 
            else 
            {
                result[idx] = d_f->d[k_d ++];
            }
        }
    }
    return result;
}

// convert mixed precision MP to FP_INT
FP_INT * mixed_FP(MP * mp, int nnz, uint8_t * mask, int n_mask)
{
    FP_INT * result = malloc(nnz * sizeof(FP_INT));
    if (!result) return NULL;

    int k_d = 0; // index for double array
    int k_u32 = 0; // index for uint32_t array
    int k_u16 = 0; // index for uint16_t array
    int k_u8 = 0; // index for uint8_t array

    for (int i = 0; i < n_mask; i ++) 
    {
        uint8_t m2 = mask[i] & 0x0F;
        uint8_t m1 = mask[i] >> 4;
        int ind1 = i * 2;
        int ind2 = i * 2 + 1;

        cast dc;

        if (m1 == 0x08) 
        {
            result[ind1] = mp->d[k_d ++];
        } 
        else 
        {
            dc.u64 = 0;
            if (m1 & 0x04) 
            {
                dc.u64 |= (uint64_t)mp->u32[k_u32 ++] << 32;
                if (m1 & 0x02) 
                {
                    dc.u64 |= (uint64_t)mp->u16[k_u16 ++] << 16;
                    if (m1 & 0x01)
                        dc.u64 |= (uint64_t)mp->u8[k_u8 ++] << 8;
                } 
                else if (m1 & 0x01) 
                {
                    dc.u64 |= (uint64_t)mp->u8[k_u8 ++] << 24;
                }
            } 
            else if (m1 & 0x02) 
            {
                dc.u64 |= (uint64_t)mp->u16[k_u16 ++] << 48;
                if (m1 & 0x01)
                    dc.u64 |= (uint64_t)mp->u8[k_u8 ++] << 40;
            } 
            else if (m1 & 0x01) 
            {
                dc.u64 |= (uint64_t)mp->u8[k_u8 ++] << 56;
            }
            result[ind1] = dc.d;
        }

        if (ind2 < nnz) 
        {
            if (m2 == 0x08) 
            {
                result[ind2] = mp->d[k_d ++];
            } 
            else 
            {
                dc.u64 = 0;
                if (m2 & 0x04) 
                {
                    dc.u64 |= (uint64_t)mp->u32[k_u32 ++] << 32;
                    if (m2 & 0x02) 
                    {
                        dc.u64 |= (uint64_t)mp->u16[k_u16 ++] << 16;
                        if (m2 & 0x01)
                            dc.u64 |= (uint64_t)mp->u8[k_u8 ++] << 8;
                    } 
                    else if (m2 & 0x01) 
                    {
                        dc.u64 |= (uint64_t)mp->u8[k_u8 ++] << 24;
                    }
                } 
                else if (m2 & 0x02) 
                {
                    dc.u64 |= (uint64_t)mp->u16[k_u16 ++] << 48;
                    if (m2 & 0x01)
                        dc.u64 |= (uint64_t)mp->u8[k_u8 ++] << 40;
                }    
                else if (m2 & 0x01) 
                {
                    dc.u64 |= (uint64_t)mp->u8[k_u8 ++] << 56;
                }
                result[ind2] = dc.d;
            }
        }
    }
    return result;
}

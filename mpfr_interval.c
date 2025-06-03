#include "convert.h"
#include "mpfr_interval.h"

void set_mpfr_mantissa_bit(mpfr_t center, int index) 
{
    mpfr_prec_t prec = mpfr_get_prec(center);
    index = int_min(index, prec - 1);

    mpfr_exp_t e_c = mpfr_get_exp(center);

    mpfr_t uls;
    mpfr_init2(uls, prec);
    mpfr_set_ui_2exp(uls, 1, e_c - index - 1, MPFR_RNDN);  // uls = 2^(e_c - index)

    mpfr_t shifted;
    mpfr_init2(shifted, prec + index + 4);

    mpfr_abs(shifted, center, MPFR_RNDN);
    mpfr_mul_2si(shifted, shifted, -(e_c - index - 1), MPFR_RNDN);

    mpfr_floor(shifted, shifted);
    mpz_t z;
    mpz_init(z);
    mpfr_get_z(z, shifted, MPFR_RNDZ);
    int bit = mpz_odd_p(z);
    mpz_clear(z);

    int sign = mpfr_sgn(center);

    // 0 at index
    if (bit == 0) 
    {
        if (sign > 0)
        {
            mpfr_add(center, center, uls, MPFR_RNDA);
        }
        else if (sign < 0)
        {
            mpfr_sub(center, center, uls, MPFR_RNDA);
        }        
    }
    mpfr_clear(uls);
    mpfr_clear(shifted);
}

int mpfr_compress(MPFR_C_R int_cr)
{
    int e_c = mpfr_get_exp(int_cr.center);
    int e_r = mpfr_get_exp(int_cr.radius);
    int prec = mpfr_get_prec(int_cr.center);
    int p = int_min(e_c - e_r, prec - 1);
    set_mpfr_mantissa_bit(int_cr.center, p);
    p += 1;
    return p;
}
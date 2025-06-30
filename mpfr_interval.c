#include "convert.h"
#include "mpfr_interval.h"

#include <stdlib.h>


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

MPFR_C_R read_mpfr_uls(mpfr_t c_tilde)
{
    MPFR_C_R int_cr;
    mpfr_prec_t prec = mpfr_get_prec(c_tilde);
    mpfr_exp_t e_c = mpfr_get_exp(c_tilde);
    mpfr_init2(int_cr.center, prec);
    mpfr_init2(int_cr.radius, 1);

    mpfr_set(int_cr.center, c_tilde, MPFR_RNDN);
    mpfr_set_ui_2exp(int_cr.radius, 1, e_c - prec, MPFR_RNDN);

    return int_cr;
}

void cr_lr(MPFR_C_R int_cr, mpfi_t int_lr)
{
    mpfr_t lo, hi;
    mpfr_inits2(mpfr_get_prec(int_cr.center), lo, hi, NULL);

    mpfr_sub(lo, int_cr.center, int_cr.radius, MPFR_RNDD);
    mpfr_add(hi, int_cr.center, int_cr.radius, MPFR_RNDU);

    mpfi_interv_fr(int_lr, lo, hi); 

    mpfr_clears(lo, hi, NULL);
}

void lr_cr(mpfi_t int_lr, MPFR_C_R int_cr) 
{
    mpfr_init2(int_cr.center, mpfi_get_prec(int_lr));
    mpfr_init2(int_cr.radius, mpfi_get_prec(int_lr));

    mpfr_t a, b, tmp;
    mpfr_inits2(mpfi_get_prec(int_lr), a, b, tmp, NULL);

    mpfi_get_left(a, int_lr);
    mpfi_get_right(b, int_lr);

    mpfr_sub(tmp, b, a, MPFR_RNDU);
    mpfr_div_ui(int_cr.radius, tmp, 2, MPFR_RNDU);

    mpfr_add(int_cr.center, a, int_cr.radius, MPFR_RNDU);

    mpfr_clears(a, b, tmp, NULL);
}



// A in band storage format diagonal n, subdiagonal n-1, superdiagonal n-1
// b is a vector of size n
void interval_GS_tridiag_mpfi(mpfi_t *A, mpfi_t *b, mpfi_t *x, int n) 
{
    mpfi_t sum, tmp, prod;
    mpfi_inits(sum, tmp, prod, NULL);

    mpfi_t *x_prev = malloc(n * sizeof(mpfi_t));
    for (int i = 0; i < n; i ++) 
    {
        mpfi_init2(x_prev[i], mpfi_get_prec(x[i]));
    }

    for (;;) 
    {
        for (int i = 0; i < n; i ++) 
        {
            mpfi_set(x_prev[i], x[i]);
        }

        for (int i = 0; i < n; i ++) 
        {
            mpfi_set(sum, b[i]);

            if (i > 0) 
            {
                mpfi_mul(prod, A[n + i - 1], x[i - 1]);
                mpfi_sub(sum, sum, prod);
            }

            if (i < n - 1) 
            {
                mpfi_mul(prod, A[2 * n - 1 + i], x_prev[i + 1]);
                mpfi_sub(sum, sum, prod);
            }

            mpfi_div(x[i], sum, A[i]);
        }

        int converged = check_convergence(x, x_prev, n, 1e-6);   
        if (converged) 
        {
            break;
        }
    }

    for (int i = 0; i < n; i ++) 
    {
        mpfi_clear(x_prev[i]);
    }

    free(x_prev);

    mpfi_clears(sum, tmp, prod, NULL);
}

int check_convergence(mpfi_t *x, mpfi_t *x_prev, int n, double tol) 
{
    mp_prec_t prec = mpfi_get_prec(x[0]); 
    mpfr_t w_curr, w_prev, diff;
    mpfr_inits2(prec, w_curr, w_prev, diff, NULL);  

    int converged = 1;

    for (int i = 0; i < n; i ++) 
    {
        mpfi_diam_abs(w_curr, x[i]);
        mpfi_diam_abs(w_prev, x_prev[i]);

        mpfr_sub(diff, w_curr, w_prev, MPFR_RNDN);
        mpfr_abs(diff, diff, MPFR_RNDN);

        if (mpfr_cmp_d(diff, tol) > 0) 
        {
            converged = 0;
            break;
        }
    }

    mpfr_clears(w_curr, w_prev, diff, NULL);
    return converged;
}

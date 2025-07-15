#include <mpfr.h>
#include <mpfi.h>

typedef struct {
    mpfr_t center;
    mpfr_t radius;
} MPFR_C_R;

void set_mpfr_mantissa_bit(mpfr_t center, int index);

int mpfr_compress(MPFR_C_R int_cr);

MPFR_C_R read_mpfr_uls3(mpfr_t c_tilde);

void interval_GS_tridiag_mpfi(mpfi_t *A, mpfi_t *b, mpfi_t *x, int n);

void interval_GS_CSR_mpfi(mpfi_t *A, int *idx, int *col_id, mpfi_t *b, mpfi_t *x, int n);

int check_convergence(mpfi_t *x, mpfi_t *x_prev, int n, double tol);

void cr_lr(MPFR_C_R int_cr, mpfi_t int_lr); 
void lr_cr(mpfi_t int_lr, MPFR_C_R * int_cr);  

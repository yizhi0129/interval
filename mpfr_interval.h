#include <mpfr.h>

typedef struct {
    mpfr_t center;
    mpfr_t radius;
} MPFR_C_R;

void set_mpfr_mantissa_bit(mpfr_t center, int index);

int mpfr_compress(MPFR_C_R int_cr);
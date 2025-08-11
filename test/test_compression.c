#include "convert.h"

#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include <math.h>
#include "mpfr_interval.h"


#define N 20000



void generate_test_intervals(MPFR_C_R *test_int, int precision) 
{
    const double c_exp_min = -6.0;
    const double c_exp_max = 9.0;
    const double c_exp_range = c_exp_max - c_exp_min;
    const double r_exp_min = (double) - precision;
    const double r_exp_max = -1.0;
    const double r_exp_range = r_exp_max - r_exp_min;

    mpfr_prec_t prec = precision;

    mpfr_t rand1, rand2, rand3, c_power, c_range, c_val, r_power, r_range, r_val, abs_c;
    mpfr_inits2(prec, rand1, rand2, rand3, c_power, c_range, c_val, r_power, r_range, r_val, abs_c, NULL);

    for (int i = 0; i < N; i ++) 
    {
        gmp_randstate_t rstate;
        gmp_randinit_default(rstate);
        gmp_randseed_ui(rstate, (unsigned long) get_time_ms());

        mpfr_inits2(prec, test_int[i].center, test_int[i].radius, NULL);
        mpfr_inits2(prec, test_int[i + N].center, test_int[i + N].radius, NULL);

        mpfr_urandomb(rand1, rstate); // [0,1)
        mpfr_urandomb(rand2, rstate); // [0,1)
        mpfr_urandomb(rand3, rstate); // [0,1)

        mpfr_set_d(c_power, c_exp_min, MPFR_RNDN);
        mpfr_set_d(c_range, c_exp_range, MPFR_RNDN);
        mpfr_fma(c_power, rand1, c_range, c_power, MPFR_RNDN); 
        mpfr_ui_pow(c_val, 10, c_power, MPFR_RNDN); // c_val = 10^c_power
        mpfr_mul(c_val, c_val, rand2, MPFR_RNDN); // c_val = rand * 10^c_power

        mpfr_set_d(r_power, r_exp_min, MPFR_RNDN);
        mpfr_set_d(r_range, r_exp_range, MPFR_RNDN);
        mpfr_fma(r_power, rand3, r_range, r_power, MPFR_RNDN);
        mpfr_ui_pow(r_val, 2, r_power, MPFR_RNDN); // two_pow = 2^r_power
        mpfr_abs(abs_c, c_val, MPFR_RNDN);
        mpfr_mul(r_val, r_val, abs_c, MPFR_RNDN);

        mpfr_set(test_int[i].center, c_val, MPFR_RNDN);
        mpfr_set(test_int[i].radius, r_val, MPFR_RNDN);

        mpfr_neg(c_val, c_val, MPFR_RNDN);

        mpfr_set(test_int[i + N].center, c_val, MPFR_RNDN);
        mpfr_set(test_int[i + N].radius, r_val, MPFR_RNDN);

        gmp_randclear(rstate);
    }
    mpfr_clears(rand1, rand2, rand3, c_power, c_range, c_val, r_power, r_range, r_val, abs_c, NULL);
}


int main(int argc, char **argv)
{
    int precision = atoi(argv[1]);

    char time[64];
    sprintf(time, "mpfr_time_%d.txt", precision);

    MPFR_C_R *test_int = malloc(2 * N * sizeof(MPFR_C_R));

    mpfr_t *array = malloc(2 * N * sizeof(mpfr_t));
    int *p = malloc(2 * N * sizeof(int));

    generate_test_intervals(test_int, precision);

    double start = get_time_ms();
    for (int i = 0; i < 2 * N; i ++)
    {
        p[i] = mpfr_compress(test_int[i]);
        mpfr_init2(array[i], p[i]);
        mpfr_set(array[i], test_int[i].center, MPFR_RNDZ);
    }
    double end = get_time_ms();

    FILE *fp = fopen(time, "a");
    fprintf(fp, "%.10f ms\n", end - start);

    for (int i = 0; i < 2 * N; i ++)
    {
        printf("Precision: %d\n", p[i]);
    }

    for (int i = 0; i < 2 * N; i ++)
    {
        mpfr_clear(array[i]);
        mpfr_clear(test_int[i].center);
        mpfr_clear(test_int[i].radius);
    }

    free(p);
    free(array);
    free(test_int);
    mpfr_free_cache();

    fclose(fp);
    return 0;
}
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
    int total = 2 * N;
    int precision = atoi(argv[1]);

    double factor_d = (double) 0.25 / (1 + 64 + precision);

    char time[64], file_ipb[64], file_ipb_rate[64];
    sprintf(time, "mpfr_time_%d.txt", precision);
    sprintf(file_ipb, "mpfr_ipb_%d.txt", precision);
    sprintf(file_ipb_rate, "mpfr_ipb_rate_%d.txt", precision);

    MPFR_C_R *test_int = malloc(total * sizeof(MPFR_C_R));

    mpfr_t *array = malloc(total * sizeof(mpfr_t));
    int *p = malloc(total * sizeof(int));

    MPFR_C_R *res_int = malloc(total * sizeof(MPFR_C_R));

    generate_test_intervals(test_int, precision);

    double start1 = get_time_ms();
    for (int i = 0; i < total; i ++)
    {
        p[i] = mpfr_compress(test_int[i]);
        mpfr_init2(array[i], p[i]);
        mpfr_set(array[i], test_int[i].center, MPFR_RNDZ);
    }
    double end1 = get_time_ms();

    double start2 = get_time_ms();
    for (int i = 0; i < total; i ++)
    {
        res_int[i] = read_mpfr_uls3(array[i]);
    }
    double end2 = get_time_ms();

    FILE *fp = fopen(time, "a");
    fprintf(fp, "%.17e\t%.17e\n", end1 - start1, end2 - start2);
    fclose(fp);


    for (int i = 0; i < total; i ++)
    {
        printf("Precision: %d\n", p[i]);
    }


    mpfr_t avg_ipb, max_ipb, min_ipb, avg_ipb2, max_ipb2, min_ipb2, avg_rate, max_rate, min_rate, temp1, temp2, temp3;
    mpfr_inits2(precision, avg_ipb, max_ipb, min_ipb, avg_ipb2, max_ipb2, min_ipb2, avg_rate, max_rate, min_rate, temp1, temp2, temp3, NULL);
    mpfr_set_d(avg_ipb, 0.0, MPFR_RNDN);
    mpfr_set_d(max_ipb, 0.0, MPFR_RNDN);
    mpfr_set_d(min_ipb, INFINITY, MPFR_RNDN);
    mpfr_set_d(avg_ipb2, 0.0, MPFR_RNDN);
    mpfr_set_d(max_ipb2, 0.0, MPFR_RNDN);
    mpfr_set_d(min_ipb2, INFINITY, MPFR_RNDN);
    mpfr_set_d(avg_rate, 0.0, MPFR_RNDN);
    mpfr_set_d(max_rate, 0.0, MPFR_RNDN);
    mpfr_set_d(min_rate, INFINITY, MPFR_RNDN);

    for (int i = 0; i < total; i ++)
    {
        mpfr_d_div(temp1, factor_d, test_int[i].radius, MPFR_RNDN);
        int n_bits = 1 + 64 + ((p[i] / 64) + !(p[i] % 64) * 1) * 64;
        mpfr_d_div(temp2, 0.5, res_int[i].radius, MPFR_RNDN);
        mpfr_div_d(temp2, temp2, (double) n_bits, MPFR_RNDN);
        mpfr_div(temp3, temp2, temp1, MPFR_RNDN);

        mpfr_max(max_ipb, max_ipb, temp1, MPFR_RNDN);
        mpfr_max(max_ipb2, max_ipb2, temp2, MPFR_RNDN);
        mpfr_max(max_rate, max_rate, temp3, MPFR_RNDN);

        mpfr_min(min_ipb, min_ipb, temp1, MPFR_RNDN);
        mpfr_min(min_ipb2, min_ipb2, temp2, MPFR_RNDN);
        mpfr_min(min_rate, min_rate, temp3, MPFR_RNDN);

        mpfr_add(avg_ipb, avg_ipb, temp1, MPFR_RNDN);
        mpfr_max(avg_ipb2, avg_ipb2, temp2, MPFR_RNDN);
        mpfr_max(avg_rate, avg_rate, temp3, MPFR_RNDN);

        mpfr_clear(array[i]);
        mpfr_clear(test_int[i].center);
        mpfr_clear(test_int[i].radius);
        mpfr_clear(res_int[i].center);
        mpfr_clear(res_int[i].radius);
    }
    mpfr_div_ui(avg_ipb, avg_ipb, total, MPFR_RNDN);
    mpfr_div_ui(avg_ipb2, avg_ipb2, total, MPFR_RNDN);
    mpfr_div_ui(avg_rate, avg_rate, total, MPFR_RNDN);

    FILE *fp_ipb = fopen(file_ipb, "a");
    mpfr_out_str(fp_ipb, 10, 17, avg_ipb, MPFR_RNDN); 
    fputc('\t', fp_ipb);
    mpfr_out_str(fp_ipb, 10, 17, min_ipb, MPFR_RNDN);
    fputc('\t', fp_ipb);
    mpfr_out_str(fp_ipb, 10, 17, max_ipb, MPFR_RNDN);
    fputc('\t', fp_ipb);
    mpfr_out_str(fp_ipb, 10, 17, avg_ipb2, MPFR_RNDN);
    fputc('\t', fp_ipb);
    mpfr_out_str(fp_ipb, 10, 17, min_ipb2, MPFR_RNDN);
    fputc('\t', fp_ipb);
    mpfr_out_str(fp_ipb, 10, 17, max_ipb2, MPFR_RNDN);
    fputc('\n', fp_ipb);
    fclose(fp_ipb);

    FILE *fp_rate = fopen(file_ipb_rate, "a");
    mpfr_out_str(fp_rate, 10, 17, avg_rate, MPFR_RNDN); 
    fputc('\t', fp_rate);
    mpfr_out_str(fp_rate, 10, 17, min_rate, MPFR_RNDN);
    fputc('\t', fp_rate);
    mpfr_out_str(fp_rate, 10, 17, max_rate, MPFR_RNDN);
    fputc('\n', fp_rate);
    fclose(fp_rate);

    free(p);
    free(array);
    free(test_int);
    free(res_int);
    mpfr_free_cache();
    
    return 0;
}
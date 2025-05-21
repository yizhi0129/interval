#include "convert.h"
#include <time.h>
#include <sys/time.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include <math.h>
#include <mpfr.h>


#define N 20000

typedef union {
    double d;
    uint64_t u64;
} cast;

double get_time_ms() 
{
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return (double)tv.tv_sec * 1000.0 + (double)tv.tv_usec / 1000.0;
}

int main(int argc, char ** argv)
{
    C_R *test_int2 = malloc(2 * N * sizeof(C_R));

    int n_mask = (2 * N + 7) >> 3; // ceil(2 * N / 8)
    uint8_t *mask = malloc(n_mask * sizeof(uint8_t));
    int count_1 = 0;

    FP_INT_PREC *cp_32bits = malloc(2 * N * sizeof(FP_INT_PREC));
    double *recovered_vals = malloc(2 * N * sizeof(double)); 

    FP_INT_PREC *cp_mpfr = malloc(2 * N * sizeof(FP_INT_PREC));
    mpfr_t *mpfr_c = malloc(2 * N * sizeof(mpfr_t));

    const double c_exp_min = -6.0;
    const double c_exp_max = 9.0;
    const double r_exp_min = -53.0;
    const double r_exp_max = -1.0;

    srand(get_time_ms());  

    for (int i = 0; i < N; i ++)
    {
        double c_power1 = c_exp_min + ((double)rand() / RAND_MAX) * (c_exp_max - c_exp_min);
        double c1 = ((double)rand() / RAND_MAX) * pow(10, c_power1);   
        double r_power1 = r_exp_min + ((double)rand() / RAND_MAX) * (r_exp_max - r_exp_min);
        double r1 = pow(2, r_power1) * fabs(c1);    
        test_int2[i].center = c1;
        test_int2[i].radius = r1;

        double c_power2 = c_exp_min + ((double)rand() / RAND_MAX) * (c_exp_max - c_exp_min);
        double c2 =  - ((double)rand() / RAND_MAX) * pow(10, c_power2);   
        double r_power2 = r_exp_min + ((double)rand() / RAND_MAX) * (r_exp_max - r_exp_min);
        double r2 = pow(2, r_power2) * fabs(c2);     
        test_int2[i + N].center = c2;
        test_int2[i + N].radius = r2;
    }

    // double to mpfr_t converssion
    double start1 = get_time_ms();
    for (int i = 0; i < 2 * N; i ++)
    {
        cp_mpfr[i] = CR_FP_mpfr(test_int2[i]);
        mpfr_init2(mpfr_c[i], cp_mpfr[i].precision);
        mpfr_set_d(mpfr_c[i], cp_mpfr[i].center, MPFR_RNDZ);
    }
    double end1 = get_time_ms();

    // double converssion
    double start2 = get_time_ms();
    for (int i = 0; i < 2 * N; i ++)
    {
        cp_32bits[i] = CR_FP_p(test_int2[i]);
    }
    double end2 = get_time_ms();

    // set mask bits: 1 if precision <= 20, 0 otherwise
    for (int i = 0; i < n_mask; i ++)
    {
        mask[i] = 0;
        for (int j = 0; j < 8; j ++)
        {
            int idx = i * 8 + j;
            if (idx >= 2 * N) break;

            if (cp_32bits[i * 8 + j].precision <= 20)
            {
                mask[i] |= (1 << j);
            }
        }
    }

    // count number of 1s in mask
    for (int i = 0; i < n_mask; i ++)
    {
        for (int j = 0; j < 8; j ++) 
        {
            int idx = i * 8 + j;
            if (idx >= 2 * N) break;
            count_1 += (mask[i] >> j) & 1;
        }
    }

    // allocate necessary memory
    int n_32 = N * 4 - count_1;
    uint32_t *res32 = malloc(n_32 * sizeof(uint32_t));

    // copy reduced results: double to uint32_t
    int k = 0;  // res32 index

    double start3 = get_time_ms();
    for (int i = 0; i < n_mask; i ++) 
    {
        for (int j = 7; j >= 0; j --) 
        {
            int idx = i * 8 + j;
            if (idx >= 2 * N) continue;

            cast dc;
            dc.d = cp_32bits[idx].center;

            if ((mask[i] >> j) & 1) 
            {
                res32[k ++] = (uint32_t)(dc.u64 >> 32);
            } 
            else 
            {
                res32[k ++] = (uint32_t)(dc.u64 >> 32);       // high part
                res32[k ++] = (uint32_t)(dc.u64 & 0xFFFFFFFF); // low part
            }
        }
    }
    double end3 = get_time_ms();

    // read uint32_t to double
    k = 0;

    double start4 = get_time_ms();
    for (int i = 0; i < n_mask; i ++) 
    {
        for (int j = 7; j >= 0; j --) 
        {
            int idx = i * 8 + j;
            if (idx >= 2 * N) continue;

            uint64_t high = (uint64_t)res32[k++] << 32;
            uint64_t low = 0;

            if (!((mask[i] >> j) & 1)) 
            {
                low = (uint64_t)res32[k++];
            }

            uint64_t combined = high | low;

            cast cast;
            cast.u64 = combined;

            recovered_vals[idx] = cast.d;
        }
    }
    double end4 = get_time_ms();

    
    FILE *fp = fopen("compression_time.txt", "a");
    fprintf(fp, "%.10f\t%.10f\t%.10f\t%.10f\tms\n", end1 - start1, end2 - start2, end3 - start3, end4 - start4);
    fclose(fp);

    for (int i = 0; i < 2 * N; i ++)
    {
        if (fabs(recovered_vals[i] - cp_32bits[i].center) > 1e-10)
        {
            printf("Error: %d\n", i);
            printf("Original: %.10e\n", cp_32bits[i].center);
            printf("Recovered: %.10e\n", recovered_vals[i]);
        }
    }

    printf("%d\t%d\n", n_mask, count_1);

    for (int i = 0; i < 2 * N; i ++)
    {
        mpfr_clear(mpfr_c[i]);
    }

    free(cp_mpfr);
    free(test_int2);
    free(mpfr_c);
    free(mask);
    free(res32);
    free(cp_32bits);

    mpfr_free_cache();

    return 0;
}
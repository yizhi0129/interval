#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <fenv.h>

#include "convert.h"
#include "functions.h"
#include "mpfr_interval.h"


void generate_ieee(C_R *A, C_R *b, int n)
{
    const double r_exp_min = -53.0;
    const double r_exp_max = -5.0; 

    srand(get_time_ms());

    for (int i = 0; i < n; i ++)
    {
        A[i].center = 4.0;
        double r_power1 = r_exp_min + ((double)rand() / RAND_MAX) * (r_exp_max - r_exp_min);
        A[i].radius = pow(2, r_power1) * fabs(A[i].center);

        b[i].center = 1.0;
        double r_power2 = r_exp_min + ((double)rand() / RAND_MAX) * (r_exp_max - r_exp_min);
        b[i].radius = pow(2, r_power2) * fabs(b[i].center);
    }

    for (int i = 0; i < n - 1; i ++)
    {
        A[i + n].center = -1.0;
        double r_power1 = r_exp_min + ((double)rand() / RAND_MAX) * (r_exp_max - r_exp_min);
        A[i + n].radius = pow(2, r_power1) * fabs(A[i + n].center);

        A[i + 2 * n - 1].center = -1.0;
        double r_power2 = r_exp_min + ((double)rand() / RAND_MAX) * (r_exp_max - r_exp_min);
        A[i + 2 * n - 1].radius = pow(2, r_power2) * fabs(A[i + 2 * n - 1].center);
    }
}

void generate_mpfr(MPFR_C_R *A, MPFR_C_R *b, int n, int precision)
{
    if (precision < 1) precision = 53;

    const double r_exp_min = (double) - precision;
    const double r_exp_max = -5.0; // 0.03125
    const double r_exp_range = r_exp_max - r_exp_min;

    mpfr_prec_t prec = (mpfr_prec_t) precision;

    gmp_randstate_t rstate;
    gmp_randinit_default(rstate);
    gmp_randseed_ui(rstate, (unsigned long) get_time_ms());

    mpfr_t rand1, rand2, r_power, r_range, r_val, abs_c, c_val;
    mpfr_inits2(prec, rand1, rand2, r_power, r_range, r_val, abs_c, c_val, NULL);
    mpfr_set_d(r_range, r_exp_range, MPFR_RNDN);

    // diagonal and right-hand side
    for (int i = 0; i < n; i ++) 
    {
        mpfr_inits2(prec, A[i].center, A[i].radius, NULL);
        mpfr_inits2(prec, b[i].center, b[i].radius, NULL);

        mpfr_urandomb(rand1, rstate);
        mpfr_set_d(c_val, 4.0, MPFR_RNDN);
        mpfr_abs(abs_c, c_val, MPFR_RNDN);

        mpfr_set_d(r_power, r_exp_min, MPFR_RNDN);
        mpfr_fma(r_power, rand1, r_range, r_power, MPFR_RNDN); // r_exp_min + rand1 * r_exp_range
        mpfr_ui_pow(r_val, 2, r_power, MPFR_RNDN); // 2 ^ r_power
        mpfr_mul(r_val, r_val, abs_c, MPFR_RNDN);  // r_val *= abs_c

        mpfr_set(A[i].center, c_val, MPFR_RNDN);
        mpfr_set(A[i].radius, r_val, MPFR_RNDN);

        mpfr_urandomb(rand2, rstate);
        mpfr_set_d(c_val, 1.0, MPFR_RNDN);
        mpfr_abs(abs_c, c_val, MPFR_RNDN);

        mpfr_set_d(r_power, r_exp_min, MPFR_RNDN);
        mpfr_fma(r_power, rand2, r_range, r_power, MPFR_RNDN);
        mpfr_ui_pow(r_val, 2, r_power, MPFR_RNDN);
        mpfr_mul(r_val, r_val, abs_c, MPFR_RNDN);

        mpfr_set(b[i].center, c_val, MPFR_RNDN);
        mpfr_set(b[i].radius, r_val, MPFR_RNDN);
    }

    // subdiagonal and superdiagonal
    for (int i = 0; i < n - 1; i ++) 
    {
        int idx1 = i + n;
        int idx2 = i + 2 * n - 1;

        mpfr_inits2(prec, A[idx1].center, A[idx1].radius, NULL);
        mpfr_inits2(prec, A[idx2].center, A[idx2].radius, NULL);

        // subdiagonal
        mpfr_urandomb(rand1, rstate);
        mpfr_set_d(c_val, -1.0, MPFR_RNDN);
        mpfr_abs(abs_c, c_val, MPFR_RNDN);

        mpfr_set_d(r_power, r_exp_min, MPFR_RNDN);
        mpfr_fma(r_power, rand1, r_range, r_power, MPFR_RNDN);
        mpfr_ui_pow(r_val, 2, r_power, MPFR_RNDN);
        mpfr_mul(r_val, r_val, abs_c, MPFR_RNDN);

        mpfr_set(A[idx1].center, c_val, MPFR_RNDN);
        mpfr_set(A[idx1].radius, r_val, MPFR_RNDN);

        // superdiagonal
        mpfr_urandomb(rand2, rstate);
        mpfr_set_d(c_val, -1.0, MPFR_RNDN);
        mpfr_abs(abs_c, c_val, MPFR_RNDN);

        mpfr_set_d(r_power, r_exp_min, MPFR_RNDN);
        mpfr_fma(r_power, rand2, r_range, r_power, MPFR_RNDN);
        mpfr_ui_pow(r_val, 2, r_power, MPFR_RNDN);
        mpfr_mul(r_val, r_val, abs_c, MPFR_RNDN);

        mpfr_set(A[idx2].center, c_val, MPFR_RNDN);
        mpfr_set(A[idx2].radius, r_val, MPFR_RNDN);
    }

    mpfr_clears(rand1, rand2, r_power, r_range, r_val, abs_c, c_val, NULL);
    gmp_randclear(rstate);
}


int main(int argc, char** argv)
{
    int n = atoi(argv[1]);

    int precision = atoi(argv[2]);

    char time[64];
    char matrix[64];
    char res[64];
    char rate_A[64];
    char rate_x[64];

    switch (precision) 
    {
        case 0: // mixed precision storage
        {
            sprintf(time, "GS_mixed_time_%d.txt", n);
            sprintf(matrix, "GS_mixed_matrix_%d.txt", n);
            sprintf(res, "GS_mixed_res_%d.txt", n);

            int nnz = 3 * n - 2; // number of non-zero elements in the banded matrix

            C_R *A = malloc(nnz * sizeof(C_R));
            C_R *b = malloc(n * sizeof(C_R));
            C_R *x1 = malloc(n * sizeof(C_R));
            C_R *x2 = malloc(n * sizeof(C_R));

            FP_INT_PREC *A_prec = malloc(nnz * sizeof(FP_INT_PREC));

            C_R * A_new = malloc(nnz * sizeof(C_R));

            generate_ieee(A, b, n);
            //printf("generated A b\n");

            // initialize x
            for (int i = 0; i < n; i ++)
            {
                x1[i].center = 0.5;
                x1[i].radius = 0.5;
                x2[i].center = 0.5;
                x2[i].radius = 0.5;
            }
            //printf("initialized x\n");

            // GS ref 
            interval_GS_tridiag(A, b, x1, n);
            //printf("GS ref done\n");

            // Convert C_R to FP_INT_PREC
            double start1 = get_time_ms();
            for (int i = 0; i < nnz; i ++)
            {
                A_prec[i] = CR_FP_p(A[i]);
            }
            double end1 = get_time_ms();

            //printf("converted C_R to FP_INT_PREC\n");

            int n_mask = (nnz + 1) >> 1; // ceil(N / 2)
            uint8_t *mask = malloc(n_mask * sizeof(uint8_t));

            int count_d = 0, count_u32 = 0, count_u16 = 0, count_u8 = 0;

            // convert FP_INT_PREC to mixed precision storage
            double start2 = get_time_ms();
            MP *A_mixed = FP_mixed(A_prec, nnz, mask, n_mask, &count_d, &count_u32, &count_u16, &count_u8);
            double end2 = get_time_ms();
            //printf("converted FP_INT_PREC to mixed precision storage\n");

            if (!A_mixed) 
            {
                fprintf(stderr, "Error: memory allocation failed for mixed precision storage.\n");
                free(A);
                free(b);
                free(x1);
                free(x2);
                free(A_prec);
                free(mask);
                return 1;
            }

            // Write matrix A to file
            double start3 = get_time_ms();
            FILE *fp_matrix = fopen(matrix, "w");
            
            fprintf(fp_matrix, "%d\t%d\n", n, n_mask);
            for (int i = 0; i < n_mask; i ++)
            {
                fprintf(fp_matrix, "%d\n", mask[i]);
            }
            for (int i = 0; i < count_d; i ++)
            {
                fprintf(fp_matrix, "%.17e\n", A_mixed->d[i]);
            }
            for (int i = 0; i < count_u32; i ++)
            {
                fprintf(fp_matrix, "%u\n", A_mixed->u32[i]);   
            }
            for (int i = 0; i < count_u16; i ++)
            {
                fprintf(fp_matrix, "%hu\n", A_mixed->u16[i]);
            }
            for (int i = 0; i < count_u8; i ++)
            {
                fprintf(fp_matrix, "%hhu\n", A_mixed->u8[i]);
            }
            
            fclose(fp_matrix);
            double end3 = get_time_ms();
            //printf("wrote matrix A to file\n");
            

            // read matrix A from file
            double start4 = get_time_ms();
            FILE *fp_matrix_read = fopen(matrix, "r");
            
            int n_read, n_mask_read;
            if (fscanf(fp_matrix_read, "%d\t%d\n", &n_read, &n_mask_read) != 2) 
            {
                fprintf(stderr, "Error: failed to read n_read and n_mask_read\n");
                fclose(fp_matrix_read);
                exit(EXIT_FAILURE);
            }

            if (n_read != n || n_mask_read != n_mask)
            {
                fprintf(stderr, "Error: matrix size mismatch.\n");
                fclose(fp_matrix_read);
                free(A);
                free(b);
                free(x1);
                free(x2);
                free(A_prec);
                free(mask);
                free(A_mixed);
                return 1;
            }   

            uint8_t *mask_read = malloc(n_mask_read * sizeof(uint8_t));

            // read mask: next n_mask_read lines
            for (int i = 0; i < n_mask_read; i ++)
            {
                if (fscanf(fp_matrix_read, "%hhu\n", &mask_read[i]) != 1) 
                {
                    fprintf(stderr, "Error: failed to read mask[%d]\n", i);
                    fclose(fp_matrix_read);
                    exit(EXIT_FAILURE);
                }
            }

            int count_d_read = 0;
            int count_32_read = 0;
            int count_16_read = 0;
            int count_8_read = 0;

            for (int i = 0; i < n_mask_read; i ++)
            {
                uint8_t m = mask_read[i];
                count_8_read += m & 1; 
                count_16_read += (m >> 1) & 1;
                count_32_read += (m >> 2) & 1;
                count_d_read += (m >> 3) & 1;
                count_8_read += (m >> 4) & 1;
                count_16_read += (m >> 5) & 1;
                count_32_read += (m >> 6) & 1;
                count_d_read += (m >> 7) & 1;
            }

            MP *A_read = malloc(sizeof(MP));
            A_read->d = malloc(count_d_read * sizeof(double));
            A_read->u32 = malloc(count_32_read * sizeof(uint32_t));
            A_read->u16 = malloc(count_16_read * sizeof(uint16_t));
            A_read->u8 = malloc(count_8_read * sizeof(uint8_t));

            // read mixed precision storage
            for (int i = 0; i < count_d_read; i ++)
            {
                if (fscanf(fp_matrix_read, "%lf\n", &A_read->d[i]) != 1) 
                {
                    fprintf(stderr, "Error: failed to read A_read->d[%d]\n", i);
                    fclose(fp_matrix_read);
                    exit(EXIT_FAILURE);
                }
            }
            for (int i = 0; i < count_32_read; i ++)
            {       
                if (fscanf(fp_matrix_read, "%u\n", &A_read->u32[i]) != 1) 
                {
                    fprintf(stderr, "Error: failed to read A_read->u32[%d]\n", i);
                    fclose(fp_matrix_read);
                    exit(EXIT_FAILURE);
                }
            }
            for (int i = 0; i < count_16_read; i ++)
            {
                if (fscanf(fp_matrix_read, "%hu\n", &A_read->u16[i]) != 1) 
                {
                    fprintf(stderr, "Error: failed to read A_read->u16[%d]\n", i);
                    fclose(fp_matrix_read);
                    exit(EXIT_FAILURE);
                }
            }
            for (int i = 0; i < count_8_read; i ++)
            {
                if (fscanf(fp_matrix_read, "%hhu\n", &A_read->u8[i]) != 1) 
                {
                    fprintf(stderr, "Error: failed to read A_read->u8[%d]\n", i);
                    fclose(fp_matrix_read);
                    exit(EXIT_FAILURE);
                }
            }

            fclose(fp_matrix_read);
            double end4 = get_time_ms();
            //printf("read matrix A from file\n");

            // Convert mixed precision storage to C_R
            FP_INT *A_tilde = NULL;

            double start5 = get_time_ms();
            A_tilde = mixed_FP(A_read, nnz, mask_read, n_mask_read);

            for (int i = 0; i < nnz; i ++)
            {
                A_new[i] = FP_CR3(A_tilde[i]);
            }
            double end5 = get_time_ms();

            //printf("converted mixed precision storage to C_R\n");

            // GS
            double start_GS = get_time_ms();
            interval_GS_tridiag(A_new, b, x2, n_read);
            double end_GS = get_time_ms();
            //printf("GS done\n");

            // check inclusion inf(x2) <= inf(x1), sup(x1) <= sup(x2)
            fesetround(FE_UPWARD);
            for (int i = 0; i < n; i ++)
            {
                if (x1[i].center + x1[i].radius > x2[i].center + x2[i].radius)
                {
                    fprintf(stderr, "Sup[%d] not included!\n", i);
                }
            }
            fesetround(FE_DOWNWARD);
            for (int i = 0; i < n; i ++)
            {
                if (x1[i].center - x1[i].radius < x2[i].center - x2[i].radius)
                {
                    fprintf(stderr, "Inf[%d] not included!\n", i);
                }
            }
            fesetround(FE_TONEAREST);


            // Write accuracy to file
            FILE *fp_res = fopen(res, "a");
            double bias = 0.0, dilat = 0.0, dilat2 = 0.0;
            double max_bias = 0.0, max_dilat = 0.0, max_dilat2 = 0.0;
            for (int i = 0; i < n_read; i ++)
            {
                double temp1 = fabs((x2[i].center - x1[i].center) / x1[i].center);
                double temp2 = fabs(x2[i].radius / x1[i].radius);
                double temp3 = fabs((x2[i].radius / x2[i].center) / (x1[i].radius / x1[i].center));
                bias += temp1;
                dilat += temp2;
                dilat2 += temp3;
                max_bias = fmax(max_bias, temp1);
                max_dilat = fmax(max_dilat, temp2);
                max_dilat2 = fmax(max_dilat2, temp3);
            }
            bias /= n_read;
            dilat /= n_read;
            dilat2 /= n_read;
            fprintf(fp_res, "%.17e\t%.17e\t%.17e\t%.17e\t%.17e\t%.17e\n", bias, max_bias, dilat, max_dilat, dilat2, max_dilat2);
            fclose(fp_res);

            // Write timing information to file
            FILE *fp_time = fopen(time, "a");
            double t_convert = end1 - start1,
                t_store = end2 - start2,
                t_write = end3 - start3,
                t_read = end4 - start4,
                t_recover = end5 - start5,
                t_compute = end_GS - start_GS;
            double perc_convert = t_convert / t_compute * 100,
                perc_store = t_store / t_compute * 100,
                perc_write = t_write / t_compute * 100,
                perc_read = t_read / t_compute * 100,
                perc_recover = t_recover / t_compute * 100;

            fprintf(fp_time, "%.17e\t%.17e\t%.17e\t%.17e\t%.17e\t%.17e\t%.17e\t%.17e\t%.17e\t%.17e\t%.17e\n", 
                t_convert, t_store, t_write, t_read, t_recover, t_compute,
                perc_convert, perc_store, perc_write, perc_read, perc_recover);
            fclose(fp_time);

            free(A);
            free(A_new);
            free(b);
            free(x1);
            free(x2);
            free(A_prec);
            free(mask);
            free(A_mixed->d);
            free(A_mixed->u32);
            free(A_mixed->u16);
            free(A_mixed->u8);
            free(A_mixed);
            free(mask_read);
            free(A_read->d);
            free(A_read->u32);
            free(A_read->u16);
            free(A_read->u8);
            free(A_read);
            free(A_tilde);
            break;
        }


        case 1: // float/double storage
        {
            sprintf(time, "GS_fd_time_%d.txt", n);
            sprintf(matrix, "GS_fd_matrix_%d.txt", n);
            sprintf(res, "GS_fd_res_%d.txt", n);
    
            int nnz = 3 * n - 2;

            C_R *A = malloc(nnz * sizeof(C_R));
            C_R *b = malloc(n * sizeof(C_R));
            C_R *x1 = malloc(n * sizeof(C_R));
            C_R *x2 = malloc(n * sizeof(C_R));

            int n_mask = (nnz + 7) >> 3; // ceil(nnz / 8)
            uint8_t *mask = malloc(n_mask * sizeof(uint8_t));
            int count = 0;

            generate_ieee(A, b, n);
            //printf("generated A b\n");

            // initialize x
            for (int i = 0; i < n; i ++)
            {
                x1[i].center = 0.5;
                x1[i].radius = 0.5;
                x2[i].center = 0.5;
                x2[i].radius = 0.5;
            }
            //printf("initialized x\n");

            // GS ref
            interval_GS_tridiag(A, b, x1, n);
            //printf("GS ref done\n");

            FP_INT_PREC *A_prec = malloc(nnz * sizeof(FP_INT_PREC));

            // Convert C_R to FP_INT_PREC
            double start1 = get_time_ms();
            for (int i = 0; i < nnz; i ++)
            {
                A_prec[i] = CR_FP_p(A[i]);
            }
            double end1 = get_time_ms();

            //printf("converted C_R to FP_INT_PREC\n");

            // Convert FP_INT_PREC to float/double storage
            D_F * A_df = NULL;
            int count_f = 0;
            double start2 = get_time_ms();
            A_df = FP_fd(A_prec, nnz, mask, n_mask, &count_f);
            double end2 = get_time_ms();
            //printf("converted FP_INT_PREC to float/double storage\n");

            // Write matrix A to file
            int count_d = nnz - count_f;
            double start3 = get_time_ms();
            FILE *fp_matrix = fopen(matrix, "w");
            
            fprintf(fp_matrix, "%d\t%d\n", n, n_mask);
            for (int i = 0; i < n_mask; i ++)
            {
                fprintf(fp_matrix, "%d\n", mask[i]);
            }
            for (int i = 0; i < count_f; i ++)
            {
                fprintf(fp_matrix, "%.8e\n", A_df->f[i]);
            }
            for (int i = 0; i < count_d; i ++)
            {
                fprintf(fp_matrix, "%.17e\n", A_df->d[i]);
            }
            
            fclose(fp_matrix);
            double end3 = get_time_ms();
            //printf("wrote matrix A to file\n");

            // read matrix A from file
            int n_read, n_mask_read;
            double start4 = get_time_ms();
            FILE *fp_matrix_read = fopen(matrix, "r");
            
            if (fscanf(fp_matrix_read, "%d\t%d\n", &n_read, &n_mask_read) != 2)
            {
                fprintf(stderr, "Error: failed to read n_read and n_mask_read\n");
                fclose(fp_matrix_read);
                exit(EXIT_FAILURE);
            }
            if (n_read != n || n_mask_read != n_mask)
            {
                fprintf(stderr, "Error: matrix size mismatch.\n");
                fclose(fp_matrix_read);
                free(A);
                free(b);
                free(x1);
                free(x2);
                free(A_prec);
                free(mask);
                free(A_df);
                return 1;
            }

            uint8_t *mask_read = malloc(n_mask_read * sizeof(uint8_t));
            D_F *A_read = malloc(sizeof(D_F));

            // read mask: next n_mask_read lines
            for (int i = 0; i < n_mask_read; i ++)
            {
                if (fscanf(fp_matrix_read, "%hhu\n", &mask_read[i]) != 1)
                {
                    fprintf(stderr, "Error: failed to read mask[%d]\n", i);
                    fclose(fp_matrix_read);
                    exit(EXIT_FAILURE);
                }
            }

            for (int i = 0; i < n_mask_read; i ++)
            {
                for (int j = 0; j < 8; j ++)
                {
                    int idx = i * 8 + j;
                    if (idx >= nnz) break;
                    count += (mask_read[i] >> j) & 1;
                }
            }

            int nf_read = count; // number of float needed
            int nd_read = 3 * n_read - 2 - count; // number of double needed

            A_read->f = malloc(nf_read * sizeof(float));
            A_read->d = malloc(nd_read * sizeof(double));

            // read A_f and A_d: next nf_read + nd_read lines
            for (int i = 0; i < nf_read; i ++)
            {
                if (fscanf(fp_matrix_read, "%f\n", &A_read->f[i]) != 1)
                {
                    fprintf(stderr, "Error: failed to read A_read->f[%d]\n", i);
                    fclose(fp_matrix_read);
                    exit(EXIT_FAILURE);
                }
            }
            for (int i = 0; i < nd_read; i ++)
            {
                if (fscanf(fp_matrix_read, "%lf\n", &A_read->d[i]) != 1)
                {
                    fprintf(stderr, "Error: failed to read A_read->d[%d]\n", i);
                    fclose(fp_matrix_read);
                    exit(EXIT_FAILURE);
                }
            }
            
            fclose(fp_matrix_read);
            double end4 = get_time_ms();
            //printf("read matrix A from file\n");
            
            // Convert mixed precision storage to C_R
            FP_INT *A_tilde = NULL;
            double start5 = get_time_ms();
            A_tilde = fd_FP(A_df, nnz, mask_read, n_mask_read);

            for (int i = 0; i < nnz; i ++)
            {
                A[i] = FP_CR3(A_tilde[i]);
            }
            double end5 = get_time_ms();

            //printf("converted float/double storage to C_R\n");

            // GS
            double start_GS = get_time_ms();
            interval_GS_tridiag(A, b, x2, n_read);
            double end_GS = get_time_ms();
            //printf("GS done\n");

            // check inclusion inf(x2) <= inf(x1), sup(x1) <= sup(x2)
            fesetround(FE_UPWARD);
            for (int i = 0; i < n; i ++)
            {
                if (x1[i].center + x1[i].radius > x2[i].center + x2[i].radius)
                {
                    fprintf(stderr, "Sup[%d] not included!\n", i);
                }
            }
            fesetround(FE_DOWNWARD);
            for (int i = 0; i < n; i ++)
            {
                if (x1[i].center - x1[i].radius < x2[i].center - x2[i].radius)
                {
                    fprintf(stderr, "Inf[%d] not included!\n", i);
                }
            }
            fesetround(FE_TONEAREST);

            // Write results to file
            FILE *fp_res = fopen(res, "a");
            double bias = 0.0, dilat = 0.0, dilat2 = 0.0; 
            double max_bias = 0.0, max_dilat = 0.0, max_dilat2 = 0.0;
            for (int i = 0; i < n_read; i ++)
            {
                double temp1 = fabs((x2[i].center - x1[i].center) / x1[i].center);
                double temp2 = fabs(x2[i].radius / x1[i].radius);
                double temp3 = fabs((x2[i].radius / x2[i].center) / (x1[i].radius / x1[i].center));
                bias += temp1;
                dilat += temp2;
                dilat2 += temp3;
                max_bias = fmax(max_bias, temp1);
                max_dilat = fmax(max_dilat, temp2);
                max_dilat2 = fmax(max_dilat2, temp3);
            }
            bias /= n_read;
            dilat /= n_read;
            dilat2 /= n_read;
            fprintf(fp_res, "%.17e\t%.17e\t%.17e\t%.17e\t%.17e\t%.17e\n", bias, max_bias, dilat, max_dilat, dilat2, max_dilat2);
            fclose(fp_res);

            // Write timing information to file
            FILE *fp_time = fopen(time, "a");
            double t_convert = end1 - start1,
                t_store = end2 - start2,
                t_write = end3 - start3,
                t_read = end4 - start4,
                t_recover = end5 - start5,
                t_compute = end_GS - start_GS;
            double perc_convert = t_convert / t_compute * 100,
                perc_store = t_store / t_compute * 100,
                perc_write = t_write / t_compute * 100,
                perc_read = t_read / t_compute * 100,
                perc_recover = t_recover / t_compute * 100;

            fprintf(fp_time, "%.17e\t%.17e\t%.17e\t%.17e\t%.17e\t%.17e\t%.17e\t%.17e\t%.17e\t%.17e\t%.17e\n", 
                t_convert, t_store, t_write, t_read, t_recover, t_compute,
                perc_convert, perc_store, perc_write, perc_read, perc_recover);
            fclose(fp_time);

            free(A);
            free(b);
            free(x1);
            free(x2);
            free(A_prec);
            free(mask);
            free(mask_read);
            free(A_read->d);
            free(A_read->f);
            free(A_read);
            free(A_df->d);
            free(A_df->f);
            free(A_df);
            free(A_tilde);

            break;
        }


        case 2: // uint32_t storage
        {
            sprintf(time, "GS_u32_time_%d.txt", n);
            sprintf(matrix, "GS_u32_matrix_%d.txt", n);
            sprintf(res, "GS_u32_res_%d.txt", n);

            int nnz = 3 * n - 2;

            C_R *A = malloc(nnz * sizeof(C_R));
            C_R *b = malloc(n * sizeof(C_R));
            C_R *x1 = malloc(n * sizeof(C_R));
            C_R *x2 = malloc(n * sizeof(C_R));

            int n_mask = (nnz + 7) >> 3;
            uint8_t *mask = malloc(n_mask * sizeof(uint8_t));
            int count = 0;

            generate_ieee(A, b, n);
            //printf("generated A b\n");

            // initialize x
            for (int i = 0; i < n; i ++)
            {
                x1[i].center = 0.5;
                x1[i].radius = 0.5;
                x2[i].center = 0.5;
                x2[i].radius = 0.5;
            }

            // GS ref
            interval_GS_tridiag(A, b, x1, n);
            //printf("GS ref done\n");

            FP_INT_PREC *A_prec = malloc(nnz * sizeof(FP_INT_PREC));

            // Convert C_R to FP_INT_PREC
            double start1 = get_time_ms();
            for (int i = 0; i < nnz; i ++)
            {
                A_prec[i] = CR_FP_p(A[i]);
            }
            double end1 = get_time_ms();

            //printf("converted C_R to FP_INT_PREC\n");

            uint32_t *A32 = NULL;
            int count32 = 0;

            // Convert FP_INT_PREC to uint32_t storage
            double start2 = get_time_ms();
            A32 = FP_u32(A_prec, nnz, mask, n_mask, &count32);
            double end2 = get_time_ms();

            //printf("converted FP_INT_PREC to uint32_t storage\n");
            //printf("count32 = %d\n", count32);

            // Write matrix A to file      
            FILE *fp_matrix = fopen(matrix, "w");

            double start3 = get_time_ms();
            fprintf(fp_matrix, "%d\t%d\n", n, n_mask);
            for (int i = 0; i < n_mask; i ++)
            {
                fprintf(fp_matrix, "%d\n", mask[i]);
            }
            for (int i = 0; i < count32; i ++)
            {
                fprintf(fp_matrix, "%u\n", A32[i]);
            }
            double end3 = get_time_ms();

            fclose(fp_matrix);
            //printf("wrote matrix A to file\n");
    
            // read matrix A from file
            int n_read, n_mask_read;
            double start4 = get_time_ms();
            FILE *fp_matrix_read = fopen(matrix, "r");

            if (fscanf(fp_matrix_read, "%d\t%d\n", &n_read, &n_mask_read) != 2)
            {
                fprintf(stderr, "Error: failed to read n_read and n_mask_read\n");
                fclose(fp_matrix_read);
                exit(EXIT_FAILURE);
            }
            if (n_read != n || n_mask_read != n_mask)
            {
                fprintf(stderr, "Error: matrix size mismatch.\n");
                fclose(fp_matrix_read);
                free(A);
                free(b);
                free(x1);
                free(x2);
                free(A_prec);
                free(mask);
                free(A32);
                return 1;
            }
            
            uint8_t *mask_read = malloc(n_mask_read * sizeof(uint8_t));

            // read mask: next n_mask_read lines
            for (int i = 0; i < n_mask_read; i ++)
            {
                if (fscanf(fp_matrix_read, "%hhu\n", &mask_read[i]) != 1)
                {
                    fprintf(stderr, "Error: failed to read mask[%d]\n", i);
                    fclose(fp_matrix_read);
                    exit(EXIT_FAILURE);
                }
            }
            
            count = 0;
            for (int i = 0; i < n_mask_read; i ++)
            {
                for (int j = 0; j < 8; j ++)
                {

                    count += (mask_read[i] >> j & 1);
                }
            }

            int n32_read = 2 * (3 * n_read - 2) - count;

            uint32_t *A32_read = malloc(n32_read * sizeof(uint32_t));

            // read matrix A: next n32_read lines
            for (int i = 0; i < n32_read; i ++)
            {
                if (fscanf(fp_matrix_read, "%u\n", &A32_read[i]) != 1)
                {
                    fprintf(stderr, "Error: failed to read A32_read[%d]\n", i);
                    fclose(fp_matrix_read);
                    exit(EXIT_FAILURE);
                }
            }
            
            fclose(fp_matrix_read);
            double end4 = get_time_ms();
            //printf("read matrix A from file\n");

            // Convert uint32_t storage to C_R
            FP_INT *A_tilde = NULL;

            double start5 = get_time_ms();
            A_tilde = u32_FP(A32_read, nnz, mask_read, n_mask_read);

            for (int i = 0; i < nnz; i ++)
            {
                A[i] = FP_CR3(A_tilde[i]);
            }
            double end5 = get_time_ms();

            //printf("converted uint32_t storage to C_R\n");


            // GS
            double start_GS = get_time_ms();
            interval_GS_tridiag(A, b, x2, n_read);
            double end_GS = get_time_ms();
            //printf("GS done\n");

            // check inclusion inf(x2) <= inf(x1), sup(x1) <= sup(x2)
            fesetround(FE_UPWARD);
            for (int i = 0; i < n; i ++)
            {
                if (x1[i].center + x1[i].radius > x2[i].center + x2[i].radius)
                {
                    fprintf(stderr, "Sup[%d] not included!\n", i);
                }
            }
            fesetround(FE_DOWNWARD);
            for (int i = 0; i < n; i ++)
            {
                if (x1[i].center - x1[i].radius < x2[i].center - x2[i].radius)
                {
                    fprintf(stderr, "Inf[%d] not included!\n", i);
                }
            }
            fesetround(FE_TONEAREST);

            // Write accuracy to file
            FILE *fp_res = fopen(res, "a");
            double bias = 0.0, dilat = 0.0, dilat2 = 0.0;
            double max_bias = 0.0, max_dilat = 0.0, max_dilat2 = 0.0;
            for (int i = 0; i < n_read; i ++)
            {
                double temp1 = fabs((x2[i].center - x1[i].center) / x1[i].center);
                double temp2 = fabs(x2[i].radius / x1[i].radius);
                double temp3 = fabs((x2[i].radius / x2[i].center) / (x1[i].radius / x1[i].center));
                bias += temp1;
                dilat += temp2;
                dilat2 += temp3;
                max_bias = fmax(max_bias, temp1);
                max_dilat = fmax(max_dilat, temp2);
                max_dilat2 = fmax(max_dilat2, temp3);
            }
            bias /= n_read;
            dilat /= n_read;
            dilat2 /= n_read;
            fprintf(fp_res, "%.17e\t%.17e\t%.17e\t%.17e\t%.17e\t%.17e\n", bias, max_bias, dilat, max_dilat, dilat2, max_dilat2);
            fclose(fp_res);

            // Write timing information to file
            FILE *fp_time = fopen(time, "a");
            double t_convert = end1 - start1,
                t_store = end2 - start2,
                t_write = end3 - start3,
                t_read = end4 - start4,
                t_recover = end5 - start5,
                t_compute = end_GS - start_GS;
            double perc_convert = t_convert / t_compute * 100,
                perc_store = t_store / t_compute * 100,
                perc_write = t_write / t_compute * 100,
                perc_read = t_read / t_compute * 100,
                perc_recover = t_recover / t_compute * 100;

            fprintf(fp_time, "%.17e\t%.17e\t%.17e\t%.17e\t%.17e\t%.17e\t%.17e\t%.17e\t%.17e\t%.17e\t%.17e\n", 
                t_convert, t_store, t_write, t_read, t_recover, t_compute,
                perc_convert, perc_store, perc_write, perc_read, perc_recover);
            fclose(fp_time);

            free(A32_read);
            free(mask_read);
            free(mask);
            free(A32);
            free(A);
            free(b);
            free(x1);
            free(x2);
            free(A_prec);
            free(A_tilde);
            break;
        }


        default: // mpfi multiprecision
        {
            sprintf(time, "GS_mpfi_time_%d_%d.txt", n, precision);
            sprintf(matrix, "GS_mpfi_matrix_%d_%d.txt", n, precision);
            sprintf(res, "GS_mpfi_res_%d_%d.txt", n, precision);
            sprintf(rate_A, "GS_mpfi_A_%d_%d.txt", n, precision);
            sprintf(rate_x, "GS_mpfi_x_%d_%d.txt", n, precision);

            int nnz = 3 * n - 2;

            MPFR_C_R *A = calloc(nnz, sizeof(MPFR_C_R));
            mpfr_prec_t *pA = malloc(nnz * sizeof(mpfr_prec_t));
            mpfr_t * A_tilde = calloc(nnz, sizeof(mpfr_t));
            MPFR_C_R *b = calloc(n, sizeof(MPFR_C_R));

            generate_mpfr(A, b, n, precision);
            
            //printf("generated A b\n");

            mpfi_t *A_mpfi = calloc(nnz, sizeof(mpfi_t));
            mpfi_t *A_mpfi2 = calloc(nnz, sizeof(mpfi_t));
            mpfi_t *b_mpfi = calloc(n, sizeof(mpfi_t));
            mpfi_t *x_mpfi = calloc(n, sizeof(mpfi_t));
            mpfi_t *x_mpfi2 = calloc(n, sizeof(mpfi_t));
            MPFR_C_R *x1 = calloc(n, sizeof(MPFR_C_R));
            MPFR_C_R *x2 = calloc(n, sizeof(MPFR_C_R));

            for (int i = 0; i < n; i ++) 
            {
                mpfr_inits2(precision, x1[i].center, x1[i].radius, x2[i].center, x2[i].radius, NULL);
            }

            for (int i = 0; i < nnz; i ++)
            {
                mpfi_init2(A_mpfi[i], precision);
                cr_lr(&A[i], A_mpfi[i]);
            } 

            for (int i = 0; i < n; i ++)
            {
                mpfi_init2(b_mpfi[i], precision);
                cr_lr(&b[i], b_mpfi[i]);
                mpfi_init2(x_mpfi[i], precision);
                mpfi_interv_d(x_mpfi[i], 0.0, 1.0); // Initialize x to [0, 1]
                mpfi_init2(x_mpfi2[i], precision);
                mpfi_interv_d(x_mpfi2[i], 0.0, 1.0); // Initialize x2 to [0, 1]
            }

            // GS ref
            interval_GS_tridiag_mpfi(A_mpfi, b_mpfi, x_mpfi, n);
            //printf("GS ref done\n");

            // conversion
            double start1 = get_time_ms();
            for (int i = 0; i < nnz; i ++)
            {
                pA[i] = (mpfr_prec_t) mpfr_compress(A[i]);
                mpfr_init2(A_tilde[i], pA[i]);
                mpfr_set(A_tilde[i], A[i].center, MPFR_RNDZ);
            }
            double end1 = get_time_ms();
            //printf("converted C_R to mpfr_t\n");


            // write matrix A to file
            double start2 = get_time_ms();
            FILE *fp_matrix = fopen(matrix, "w");
            fprintf(fp_matrix, "%d\t%d\n", n, nnz);
            for (int i = 0; i < nnz; i ++)
            {
                fprintf(fp_matrix, "%d\n", (int)pA[i]);
            }
            for (int i = 0; i < nnz; i ++)
            {
                mpfr_out_str(fp_matrix, 2, 0, A_tilde[i], MPFR_RNDN);
                fprintf(fp_matrix, "\n");
            }
            fclose(fp_matrix);
            double end2 = get_time_ms();
            //printf("wrote matrix A to file\n");


            // read matrix A from file
            int n_read, nnz_read;
            double start3 = get_time_ms();
            FILE *fp_matrix_read = fopen(matrix, "r");
            
            if (fscanf(fp_matrix_read, "%d\t%d\n", &n_read, &nnz_read) != 2)
            {
                fprintf(stderr, "Error: failed to read n_read and nnz_read\n");
                fclose(fp_matrix_read);
                exit(EXIT_FAILURE);
            }
            if (n_read != n || nnz_read != nnz)
            {
                fprintf(stderr, "Error: matrix size mismatch.\n");
                fclose(fp_matrix_read);
                free(A);
                free(b);
                free(pA);
                free(A_tilde);
                free(A_mpfi);
                free(b_mpfi);
                free(x_mpfi);
                free(x_mpfi2);
                return 1;
            }
            int *pA_read = malloc(nnz_read * sizeof(int));
            mpfr_t *A_tilde_read = calloc(nnz_read, sizeof(mpfr_t));
            for (int i = 0; i < nnz_read; i ++)
            {
                if (fscanf(fp_matrix_read, "%d\n", &pA_read[i]) != 1)
                {
                    fprintf(stderr, "Error: failed to read pA_read[%d]\n", i);
                    fclose(fp_matrix_read);
                    exit(EXIT_FAILURE);
                }
            }
            for (int i = 0; i < nnz_read; i ++)
            {
                mpfr_init2(A_tilde_read[i], (mpfr_prec_t)pA_read[i]);
                mpfr_inp_str(A_tilde_read[i], fp_matrix_read, 2, MPFR_RNDN);
            }
            fclose(fp_matrix_read);
            double end3 = get_time_ms();
            //printf("read matrix A from file\n");

            // read radius from A_tilde, convert to mpfi_t
            double start4 = get_time_ms();
            for (int i = 0; i < nnz; i ++)
            {
                MPFR_C_R tmp;
                mpfi_init2(A_mpfi2[i], (mpfr_prec_t)pA_read[i]);
                tmp = read_mpfr_uls3(A_tilde_read[i]);
                cr_lr(&tmp, A_mpfi2[i]);
                mpfr_clear(tmp.center);
                mpfr_clear(tmp.radius);
            }
            double end4 = get_time_ms();    
            //printf("read radius from A_tilde\n");    


            // Gauss-Seidel
            double start_GS = get_time_ms();
            interval_GS_tridiag_mpfi(A_mpfi2, b_mpfi, x_mpfi2, n);
            double end_GS = get_time_ms();
            //printf("GS done\n");

            // check inclusion
            for (int i = 0; i < n; i ++)
            {
                if (!mpfi_is_inside(x_mpfi[i], x_mpfi2[i]))
                {
                    fprintf(stderr, "[%d] Not included!\n", i);
                }
            }

            // accuracy check
            for (int i = 0; i < n; i ++)
            {
                lr_cr(x_mpfi[i], &x1[i]);
                lr_cr(x_mpfi2[i], &x2[i]);
            }

            FILE *fp_res = fopen(res, "a");
            FILE *fp_A = fopen(rate_A, "w");
            FILE *fp_x = fopen(rate_x, "w");
            mpfr_t bias, dilat, dilat2, max_bias, max_dilat, max_dilat2, temp1, temp2, temp3, temp3b;
            mpfr_inits2(precision, bias, dilat, dilat2, max_bias, max_dilat, max_dilat2, temp1, temp2, temp3, temp3b, NULL);
            mpfr_set_d(bias, 0.0, MPFR_RNDN);
            mpfr_set_d(dilat, 0.0, MPFR_RNDN);
            mpfr_set_d(dilat2, 0.0, MPFR_RNDN);
            mpfr_set_d(max_bias, 0.0, MPFR_RNDN);
            mpfr_set_d(max_dilat, 0.0, MPFR_RNDN);
            mpfr_set_d(max_dilat2, 0.0, MPFR_RNDN);
            mpfr_set_d(temp1, 0.0, MPFR_RNDN);
            mpfr_set_d(temp2, 0.0, MPFR_RNDN);
            mpfr_set_d(temp3, 0.0, MPFR_RNDN);
            mpfr_set_d(temp3b, 0.0, MPFR_RNDN);

            for (int i = 0; i < n_read; i ++)
            {
                mpfr_sub(temp1, x2[i].center, x1[i].center, MPFR_RNDN);
                mpfr_div(temp1, temp1, x1[i].center, MPFR_RNDN);
                mpfr_abs(temp1, temp1, MPFR_RNDN);

                mpfr_div(temp2, x2[i].radius, x1[i].radius, MPFR_RNDN);
                mpfr_abs(temp2, temp2, MPFR_RNDN);

                fprintf(fp_x, "%d\t", i);
                mpfr_out_str(fp_x, 10, 17, temp2, MPFR_RNDN); 
                fputc('\n', fp_x);

                mpfr_div(temp3, x2[i].radius, x2[i].center, MPFR_RNDN);
                mpfr_div(temp3b, x1[i].radius, x1[i].center, MPFR_RNDN);
                mpfr_div(temp3, temp3, temp3b, MPFR_RNDN);
                mpfr_abs(temp3, temp3, MPFR_RNDN);

                mpfr_max(max_bias, max_bias, temp1, MPFR_RNDN);
                mpfr_max(max_dilat, max_dilat, temp2, MPFR_RNDN);
                mpfr_max(max_dilat2, max_dilat2, temp3, MPFR_RNDN);

                mpfr_add(bias, bias, temp1, MPFR_RNDN);
                mpfr_add(dilat, dilat, temp2, MPFR_RNDN);
                mpfr_add(dilat2, dilat2, temp3, MPFR_RNDN);
            }
            for (int i = 0; i < n_read; i ++)
            {
                mpfr_sub(temp1, &A_mpfi2[i]->right, &A_mpfi2[i]->left, MPFR_RNDU);
                mpfr_div_d(temp1, temp1, 2.0, MPFR_RNDU);
                mpfr_div(temp2, temp1, A[i].radius, MPFR_RNDN);
                mpfr_abs(temp2, temp2, MPFR_RNDN);

                fprintf(fp_A, "%d\t%d\t", i + 1, i + 1);
                mpfr_out_str(fp_A, 10, 17, temp2, MPFR_RNDN); 
                fputc('\n', fp_A);
            }
            for (int i = 0; i < n_read - 1; i ++)
            {
                int i1 = i + n_read;
                int i2 = i + 2 * n_read - 1;
                mpfr_sub(temp1, &A_mpfi2[i1]->right, &A_mpfi2[i1]->left, MPFR_RNDU);
                mpfr_div_d(temp1, temp1, 2.0, MPFR_RNDU);
                mpfr_div(temp2, temp1, A[i1].radius, MPFR_RNDN);
                mpfr_abs(temp2, temp2, MPFR_RNDN);

                fprintf(fp_A, "%d\t%d\t", i + 2, i + 1);
                mpfr_out_str(fp_A, 10, 17, temp2, MPFR_RNDN); 
                fputc('\n', fp_A);

                mpfr_sub(temp1, &A_mpfi2[i2]->right, &A_mpfi2[i2]->left, MPFR_RNDU);
                mpfr_div_d(temp1, temp1, 2.0, MPFR_RNDU);
                mpfr_div(temp2, temp1, A[i2].radius, MPFR_RNDN);
                mpfr_abs(temp2, temp2, MPFR_RNDN);

                fprintf(fp_A, "%d\t%d\t", i + 1, i + 2);
                mpfr_out_str(fp_A, 10, 17, temp2, MPFR_RNDN); 
                fputc('\n', fp_A);
            }

            mpfr_clears(temp1, temp2, temp3, temp3b, NULL);
            fclose(fp_A);
            fclose(fp_x);

            mpfr_div_ui(bias, bias, n_read, MPFR_RNDN);
            mpfr_div_ui(dilat, dilat, n_read, MPFR_RNDN);
            mpfr_div_ui(dilat2, dilat2, n_read, MPFR_RNDN);

            mpfr_out_str(fp_res, 10, 17, bias, MPFR_RNDN); 
            fputc('\t', fp_res);
            mpfr_out_str(fp_res, 10, 17, max_bias, MPFR_RNDN);
            fputc('\t', fp_res);
            mpfr_out_str(fp_res, 10, 17, dilat, MPFR_RNDN);
            fputc('\t', fp_res);
            mpfr_out_str(fp_res, 10, 17, max_dilat, MPFR_RNDN);
            fputc('\t', fp_res);
            mpfr_out_str(fp_res, 10, 17, dilat2, MPFR_RNDN);
            fputc('\t', fp_res);
            mpfr_out_str(fp_res, 10, 17, max_dilat2, MPFR_RNDN);
            fputc('\n', fp_res);
            fclose(fp_res);
           
            // Write timing information to file
            FILE *fp_time = fopen(time, "a");
            double t_convert = end1 - start1,
                t_write = end2 - start2,
                t_read = end3 - start3,
                t_recover = end4 - start4,
                t_compute = end_GS - start_GS;
            double perc_convert = t_convert / t_compute * 100,
                perc_write = t_write / t_compute * 100,
                perc_read = t_read / t_compute * 100,
                perc_recover = t_recover / t_compute * 100;

            fprintf(fp_time, "%.17e\t%.17e\t%.17e\t%.17e\t%.17e\t%.17e\t%.17e\t%.17e\t%.17e\n", 
                t_convert, t_write, t_read, t_recover, t_compute,
                perc_convert, perc_write, perc_read, perc_recover);
            fclose(fp_time);

            for (int i = 0; i < n; i ++)
            {
                mpfi_clear(x_mpfi[i]);
                mpfi_clear(b_mpfi[i]);
                mpfr_clear(b[i].center);
                mpfr_clear(b[i].radius);
                mpfi_clear(x_mpfi2[i]);
                mpfr_clears(x2[i].center, x2[i].radius, x1[i].center, x1[i].radius, NULL);
            }

            for (int i = 0; i < nnz; i ++)
            {
                mpfi_clear(A_mpfi[i]);
                mpfi_clear(A_mpfi2[i]);
                mpfr_clear(A[i].center);
                mpfr_clear(A[i].radius);
                mpfr_clear(A_tilde_read[i]);
                mpfr_clear(A_tilde[i]);
            }

            mpfr_clears(bias, dilat, dilat2, max_bias, max_dilat, max_dilat2, NULL);

            mpfr_free_cache();


            free(x_mpfi);
            free(A_mpfi);
            free(b_mpfi);
            free(A);
            free(b);
            free(pA);
            free(A_tilde);
            free(A_tilde_read);
            free(pA_read);
            free(x_mpfi2);
            free(x1);
            free(x2);
            free(A_mpfi2);

            break;
        }
    }

    return 0;
}
#include "softposit.h"
#include "posit_interval.h"

#include "convert.h"

#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define HALF_16 convertDoubleToP16(0.5)
#define HALF_32 convertDoubleToP32(0.5)

#define S_16 convertDoubleToP16(0.7)
#define S_32 convertDoubleToP32(0.7)

void generate_rhs_p8(int n, C_R_p8 *rhs_p8)
{
    const double r_exp_min = -7.0;
    const double r_exp_max = -5.0;

    srand(get_time_ms());

    for (int i = 0; i < n; i ++)
    {
        double c = 1.0;
        double r_power = r_exp_min + ((double)rand() / RAND_MAX) * (r_exp_max - r_exp_min);
        double r = pow(2, r_power) * fabs(c);
        rhs_p8[i].center = convertDoubleToP8(c);
        rhs_p8[i].radius = convertDoubleToP8(r);
    }
}

void generate_rhs_p16(int n, C_R_p16 *rhs_p16)
{
    const double r_exp_min = -30.0;
    const double r_exp_max = -5.0;

    srand(get_time_ms());

    for (int i = 0; i < n; i ++)
    {
        double c = 1.0;
        double r_power = r_exp_min + ((double)rand() / RAND_MAX) * (r_exp_max - r_exp_min);
        double r = pow(2, r_power) * fabs(c);
        rhs_p16[i].center = convertDoubleToP16(c);
        rhs_p16[i].radius = convertDoubleToP16(r);
    }
}

void generate_matrix_tridiag_p8(int n, C_R_p8 *A)
{
    const double r_exp_min = -7.0;
    const double r_exp_max = -5.0; 

    srand(get_time_ms());

    for (int i = 0; i < n; i ++)
    {
        double c = 4.0;
        A[i].center = convertDoubleToP8(c);
        double r_power1 = r_exp_min + ((double)rand() / RAND_MAX) * (r_exp_max - r_exp_min);
        double r = pow(2, r_power1) * fabs(c);
        A[i].radius = convertDoubleToP8(r);
    }

    for (int i = 0; i < n - 1; i ++)
    {
        double c = -1.0;
        A[i + n].center = convertDoubleToP8(c);
        double r_power1 = r_exp_min + ((double)rand() / RAND_MAX) * (r_exp_max - r_exp_min);
        double r = pow(2, r_power1) * fabs(c);
        A[i + n].radius = convertDoubleToP8(r);

        A[i + 2 * n - 1].center = convertDoubleToP8(c);
        double r_power2 = r_exp_min + ((double)rand() / RAND_MAX) * (r_exp_max - r_exp_min);
        double r2 = pow(2, r_power2) * fabs(c);
        A[i + 2 * n - 1].radius = convertDoubleToP8(r2);
    }
}

void generate_matrix_tridiag_p16(int n, C_R_p16 *A)
{
    const double r_exp_min = -30.0;
    const double r_exp_max = -5.0; 

    srand(get_time_ms());

    for (int i = 0; i < n; i ++)
    {
        double c = 4.0;
        A[i].center = convertDoubleToP16(c);
        double r_power1 = r_exp_min + ((double)rand() / RAND_MAX) * (r_exp_max - r_exp_min);
        double r = pow(2, r_power1) * fabs(c);
        A[i].radius = convertDoubleToP16(r);
    }

    for (int i = 0; i < n - 1; i ++)
    {
        double c = -1.0;
        A[i + n].center = convertDoubleToP16(c);
        double r_power1 = r_exp_min + ((double)rand() / RAND_MAX) * (r_exp_max - r_exp_min);
        double r = pow(2, r_power1) * fabs(c);
        A[i + n].radius = convertDoubleToP16(r);

        A[i + 2 * n - 1].center = convertDoubleToP16(c);
        double r_power2 = r_exp_min + ((double)rand() / RAND_MAX) * (r_exp_max - r_exp_min);
        double r2 = pow(2, r_power2) * fabs(c);
        A[i + 2 * n - 1].radius = convertDoubleToP16(r2);
    }
}

void generate_matrix_csr_p8(int n, C_R_p8 ** A, int **idx, int **col_id)
{
    const double r_exp_min = -7.0;
    const double r_exp_max = -5.0; 

    *idx = malloc((n + 1) * sizeof(int));
    (*idx)[0] = 0;
    int pos = 0;

    srand(get_time_ms());

    int max_nnz = 5 * n;

    C_R_p8 *A_tmp = (C_R_p8 *)malloc(max_nnz * sizeof(C_R_p8));
    int *col_id_tmp = (int *)malloc(max_nnz * sizeof(int));

    for (int i = 0; i < n; i ++)
    {
        for (int j = 0; j < n; j ++)
        {
            if (j == i)
            {
                double c = 4.0;
                A_tmp[pos].center = convertDoubleToP8(c);
                double r_power = r_exp_min + ((double)rand() / RAND_MAX) * (r_exp_max - r_exp_min);
                double r = pow(2, r_power) * fabs(c);
                A_tmp[pos].radius = convertDoubleToP8(r);
                col_id_tmp[pos] = j;
                pos ++;
            }
            else if ((double)rand() / RAND_MAX < (1.0 / n))
            {
                double c = -1.0;
                A_tmp[pos].center = convertDoubleToP8(c);
                double r_power = r_exp_min + ((double)rand() / RAND_MAX) * (r_exp_max - r_exp_min);
                double r = pow(2, r_power) * fabs(c);
                A_tmp[pos].radius = convertDoubleToP8(r);
                col_id_tmp[pos] = j;
                pos ++;
            }
        }
        (*idx)[i + 1] = pos;
    }
    
    *A = (C_R_p8 *)realloc(A_tmp, pos * sizeof(C_R_p8));
    *col_id = (int *)realloc(col_id_tmp, pos * sizeof(int));
}

void generate_matrix_csr_p16(int n, C_R_p16 ** A, int **idx, int **col_id)
{
    const double r_exp_min = -30.0;
    const double r_exp_max = -5.0; 

    *idx = malloc((n + 1) * sizeof(int));
    (*idx)[0] = 0;
    int pos = 0;

    srand(get_time_ms());

    int max_nnz = 5 * n;

    C_R_p16 *A_tmp = (C_R_p16 *)malloc(max_nnz * sizeof(C_R_p16));
    int *col_id_tmp = (int *)malloc(max_nnz * sizeof(int));

    for (int i = 0; i < n; i ++)
    {
        for (int j = 0; j < n; j ++)
        {
            if (j == i)
            {
                double c = 4.0;
                A_tmp[pos].center = convertDoubleToP16(c);
                double r_power = r_exp_min + ((double)rand() / RAND_MAX) * (r_exp_max - r_exp_min);
                double r = pow(2, r_power) * fabs(c);
                A_tmp[pos].radius = convertDoubleToP16(r);
                col_id_tmp[pos] = j;
                pos ++;
            }
            else if ((double)rand() / RAND_MAX < (1.0 / n))
            {
                double c = -1.0;
                A_tmp[pos].center = convertDoubleToP16(c);
                double r_power = r_exp_min + ((double)rand() / RAND_MAX) * (r_exp_max - r_exp_min);
                double r = pow(2, r_power) * fabs(c);
                A_tmp[pos].radius = convertDoubleToP16(r);
                col_id_tmp[pos] = j;
                pos ++;
            }
        }
        (*idx)[i + 1] = pos;
    }

    *A = (C_R_p16 *)realloc(A_tmp, pos * sizeof(C_R_p16));
    *col_id = (int *)realloc(col_id_tmp, pos * sizeof(int));
}


int main(int argc, char **argv)
{
    printf("HALF_16 = ");
    p16_print_binary(HALF_16);
    printf("HALF_32 = ");
    p32_print_binary(HALF_32);

    printf("S_16 = ");
    p16_print_binary(S_16);
    printf("%lf\n", convertP16ToDouble(S_16));
    printf("S_32 = ");
    p32_print_binary(S_32);
    printf("%lf\n", convertP32ToDouble(S_32));

    int n = atoi(argv[1]);

    char time1_16[64], time1_32[64], time2_16[64], time2_32[64];

    char matrix1_16[64], matrix1_32[64], matrix2_16[64], matrix2_32[64];

    char res1_16[64], res1_32[64], res2_16[64], res2_32[64];

    sprintf(time1_16, "GS_tridiag_p16_time_%d.txt", n);
    sprintf(time1_32, "GS_tridiag_p32_time_%d.txt", n);
    sprintf(time2_16, "GS_csr_p16_time_%d.txt", n);
    sprintf(time2_32, "GS_csr_p32_time_%d.txt", n);
    
    sprintf(matrix1_16, "GS_tridiag_p16_matrix_%d.txt", n);
    sprintf(matrix1_32, "GS_tridiag_p32_matrix_%d.txt", n);
    sprintf(matrix2_16, "GS_csr_p16_matrix_%d.txt", n);
    sprintf(matrix2_32, "GS_csr_p32_matrix_%d.txt", n);

    sprintf(res1_16, "GS_tridiag_p16_res_%d.txt", n);
    sprintf(res1_32, "GS_tridiag_p32_res_%d.txt", n);
    sprintf(res2_16, "GS_csr_p16_res_%d.txt", n);
    sprintf(res2_32, "GS_csr_p32_res_%d.txt", n);

    int nnz1 = 3 * n - 2;

    C_R_p8 *A1_8 = malloc(nnz1 * sizeof(C_R_p8));
    C_R_p16 *A1_8_extend = malloc(nnz1 * sizeof(C_R_p16));
    p16_FP_INT_P *A1_16_compress = malloc(nnz1 * sizeof(p16_FP_INT_P)); 
    C_R_p16 *A1_8_ext_new = malloc(nnz1 * sizeof(C_R_p16));
    C_R_p8 *b1_8 = malloc(n * sizeof(C_R_p8));
    C_R_p16 *b1_8_extend = malloc(n * sizeof(C_R_p16));
    C_R_p16 *x_ref1_16 = malloc(n * sizeof(C_R_p16));
    C_R_p16 *x_res1_16 = malloc(n * sizeof(C_R_p16));

    C_R_p16 *A1_16 = malloc(nnz1 * sizeof(C_R_p16));
    C_R_p32 *A1_16_extend = malloc(nnz1 * sizeof(C_R_p32));
    p32_FP_INT_P *A1_32_compress = malloc(nnz1 * sizeof(p32_FP_INT_P)); 
    C_R_p32 *A1_16_ext_new = malloc(nnz1 * sizeof(C_R_p32));
    C_R_p16 *b1_16 = malloc(n * sizeof(C_R_p16));
    C_R_p32 *b1_16_extend = malloc(n * sizeof(C_R_p32));
    C_R_p32 *x_ref1_32 = malloc(n * sizeof(C_R_p32));
    C_R_p32 *x_res1_32 = malloc(n * sizeof(C_R_p32));

    generate_matrix_tridiag_p8(n, A1_8);
    printf("Matrix A1_8 generated\n");

    generate_rhs_p8(n, b1_8);
    printf("RHS b1_8 generated\n");


    generate_matrix_tridiag_p16(n, A1_16);
    printf("Matrix A1_16 generated\n");

    generate_rhs_p16(n, b1_16);
    printf("RHS b1_16 generated\n");

    for (int i = 0; i < n; i ++)
    {
        b1_8_extend[i].center = p8_to_p16(b1_8[i].center);
        b1_8_extend[i].radius = p8_to_p16(b1_8[i].radius);
        b1_16_extend[1].center = p16_to_p32(b1_16[i].center);
        b1_16_extend[1].radius = p16_to_p32(b1_16[i].radius);
    }

    for (int i = 0; i < nnz1; i ++)
    {
        A1_8_extend[i].center = p8_to_p16(A1_8[i].center);
        A1_8_extend[i].radius = p8_to_p16(A1_8[i].radius);
        A1_16_extend[i].center = p16_to_p32(A1_16[i].center);
        A1_16_extend[i].radius = p16_to_p32(A1_16[i].radius);
    }

    for (int i = 0; i < n; i ++)
    {
        x_ref1_16[i].center = HALF_16;
        x_ref1_16[i].radius = HALF_16;
        x_res1_16[i].center = HALF_16;
        x_res1_16[i].radius = HALF_16;

        x_ref1_32[i].center = HALF_32;
        x_ref1_32[i].radius = HALF_32;
        x_res1_32[i].center = HALF_32;
        x_res1_32[i].radius = HALF_32;
    }

    GS_tridiag_p16(n, A1_8_extend, b1_8_extend, x_ref1_16);
    printf("ref GS_tridiag_p16 completed\n");

    GS_tridiag_p32(n, A1_16_extend, b1_16_extend, x_ref1_32);
    printf("ref GS_tridiag_p32 completed\n");

    double start1_16 = get_time_ms();
    for (int i = 0; i < nnz1; i ++)
    {
        A1_16_compress[i] = p8_compression(A1_8[i].center, A1_8[i].radius);
    }
    double end1_16 = get_time_ms();
    printf("A1_16_compress completed\n");

    double start1_32 = get_time_ms();
    for (int i = 0; i < nnz1; i ++)
    {
        A1_32_compress[i] = p16_compression(A1_16[i].center, A1_16[i].radius);
    }
    double end1_32 = get_time_ms();
    printf("A1_32_compress completed\n");

    double start2_16 = get_time_ms();
    for (int i = 0; i < nnz1; i ++)
    {
        A1_8_ext_new[i].center = A1_16_compress[i].posit16;
        A1_8_ext_new[i].radius = p16_read_3r(A1_16_compress[i].posit16);
    }
    double end2_16 = get_time_ms();
    printf("A1_8_ext_new completed\n");

    double start2_32 = get_time_ms();
    for (int i = 0; i < nnz1; i ++)
    {
        A1_16_ext_new[i].center = A1_32_compress[i].posit32;
        A1_16_ext_new[i].radius = p32_read_3r(A1_32_compress[i].posit32);
    }
    double end2_32 = get_time_ms();
    printf("A1_16_ext_new completed\n");

    double start3_16 = get_time_ms();
    GS_tridiag_p16(n, A1_8_ext_new, b1_8_extend, x_res1_16);
    double end3_16 = get_time_ms();
    printf("GS_tridiag_p16 completed\n");

    double start3_32 = get_time_ms();
    GS_tridiag_p32(n, A1_16_ext_new, b1_16_extend, x_res1_32);
    double end3_32 = get_time_ms();
    printf("GS_tridiag_p32 completed\n");

    FILE *fp_res1 = fopen(res1_16, "a");
    double bias = 0.0, dilat = 0.0;
    double max_bias = 0.0, max_dilat = 0.0;
    for (int i = 0; i < n; i ++)
    {
        double c1_16 = convertP16ToDouble(x_ref1_16[i].center);
        double r1_16 = convertP16ToDouble(x_ref1_16[i].radius);
        double c2_16 = convertP16ToDouble(x_res1_16[i].center);
        double r2_16 = convertP16ToDouble(x_res1_16[i].radius);
        double temp1 = fabs((c2_16 - c1_16) / c1_16);
        double temp2 = fabs(r2_16 / r1_16);
        bias += temp1;
        dilat += temp2;
        max_bias = fmax(max_bias, temp1);
        max_dilat = fmax(max_dilat, temp2);              
    }
    bias /= n;
    dilat /= n;
    fprintf(fp_res1, "%.17e\t%.17e\t%.17e\t%.17e\t\n", bias, max_bias, dilat, max_dilat);
    fclose(fp_res1);

    FILE *fp_res2 = fopen(res1_32, "a");
    bias = 0.0, dilat = 0.0;
    max_bias = 0.0, max_dilat = 0.0;
    for (int i = 0; i < n; i ++)
    {
        double c1_32 = convertP32ToDouble(x_ref1_32[i].center);
        double r1_32 = convertP32ToDouble(x_ref1_32[i].radius);
        double c2_32 = convertP32ToDouble(x_res1_32[i].center);
        double r2_32 = convertP32ToDouble(x_res1_32[i].radius);
        double temp1 = fabs((c2_32 - c1_32) / c1_32);
        double temp2 = fabs(r2_32 / r1_32);
        bias += temp1;
        dilat += temp2;
        max_bias = fmax(max_bias, temp1);
        max_dilat = fmax(max_dilat, temp2);              
    }
    bias /= n;
    dilat /= n;
    fprintf(fp_res2, "%.17e\t%.17e\t%.17e\t%.17e\t\n", bias, max_bias, dilat, max_dilat);
    fclose(fp_res2);

    FILE *fp_time1 = fopen(time1_16, "a");
    fprintf(fp_time1, "%.17e\t%.17e\t%.17e\t(ms)\n", 
        end1_16 - start1_16, end2_16 - start2_16, end3_16 - start3_16);
    fclose(fp_time1);

    FILE *fp_time2 = fopen(time1_32, "a");
    fprintf(fp_time2, "%.17e\t%.17e\t%.17e\t(ms)\n", 
        end1_32 - start1_32, end2_32 - start2_32, end3_32 - start3_32);
    fclose(fp_time2);


    C_R_p8 *A2_8 = NULL;
    int *idx_16 = malloc((n + 1) * sizeof(int));
    int *col_id_16 = NULL;
    C_R_p8 *b2_8 = malloc(n * sizeof(C_R_p8));
    C_R_p16 *b2_8_extend = malloc(n * sizeof(C_R_p16));
    C_R_p16 *x_ref2_16 = malloc(n * sizeof(C_R_p16));
    C_R_p16 *x_res2_16 = malloc(n * sizeof(C_R_p16));

    C_R_p16 *A2_16 = NULL;
    int *idx_32 = malloc((n + 1) * sizeof(int));
    int *col_id_32 = NULL;
    C_R_p16 *b2_16 = malloc(n * sizeof(C_R_p16));
    C_R_p32 *b2_16_extend = malloc(n * sizeof(C_R_p32));
    C_R_p32 *x_ref2_32 = malloc(n * sizeof(C_R_p32));
    C_R_p32 *x_res2_32 = malloc(n * sizeof(C_R_p32));
    
    generate_matrix_csr_p8(n, &A2_8, &idx_16, &col_id_16);
    printf("Matrix A2_8 generated\n");

    generate_rhs_p8(n, b2_8);
    printf("RHS b2_8 generated\n");

    int nnz2_16 = idx_16[n] - idx_16[0];
    C_R_p16 *A2_8_extend = malloc(nnz2_16 * sizeof(C_R_p16));
    p16_FP_INT_P *A2_16_compress = malloc(nnz2_16 * sizeof(p16_FP_INT_P)); 
    C_R_p16 *A2_8_ext_new = malloc(nnz2_16 * sizeof(C_R_p16));

    generate_matrix_csr_p16(n, &A2_16, &idx_32, &col_id_32);
    printf("Matrix A2_16 generated\n");

    generate_rhs_p16(n, b2_16);
    printf("RHS b2_16 generated\n");

    int nnz2_32 = idx_32[n] - idx_32[0];
    C_R_p32 *A2_16_extend = malloc(nnz2_32 * sizeof(C_R_p32));
    p32_FP_INT_P *A2_32_compress = malloc(nnz2_32 * sizeof(p32_FP_INT_P)); 
    C_R_p32 *A2_16_ext_new = malloc(nnz2_32 * sizeof(C_R_p32));

    for (int i = 0; i < n; i ++)
    {
        b2_8_extend[i].center = p8_to_p16(b2_8[i].center); 
        b2_8_extend[i].radius = p8_to_p16(b2_8[i].radius); 

        b2_16_extend[i].center = p16_to_p32(b2_16[i].center); 
        b2_16_extend[i].radius = p16_to_p32(b2_16[i].radius); 
    }
    printf("RHS b2_8_extend and b2_16_extend generated\n");

    for (int i = 0; i < nnz2_16; i ++)
    {
        A2_8_extend[i].center = p8_to_p16(A2_8[i].center);
        A2_8_extend[i].radius = p8_to_p16(A2_8[i].radius);
    }
    printf("Matrix A2_8_extend generated\n");

    for (int i = 0; i < nnz2_32; i ++)
    {
        A2_16_extend[i].center = p16_to_p32(A2_16[i].center);
        A2_16_extend[i].radius = p16_to_p32(A2_16[i].radius);
    }
    printf("Matrix A2_16_extend generated\n");

    for (int i = 0; i < n; i ++)
    {
        x_ref2_16[i].center = S_16;
        x_ref2_16[i].radius = S_16;
        x_res2_16[i].center = S_16;
        x_res2_16[i].radius = S_16;

        x_ref2_32[i].center = S_32;
        x_ref2_32[i].radius = S_32;
        x_res2_32[i].center = S_32;
        x_res2_32[i].radius = S_32;        
    }

    GS_csr_p16(n, A2_8_extend, idx_16, col_id_16, b2_8_extend, x_ref2_16);
    printf("ref GS_csr_p16 completed\n");

    GS_csr_p32(n, A2_16_extend, idx_32, col_id_32, b2_16_extend, x_ref2_32);
    printf("ref GS_csr_p32 completed\n");

    double start11_16 = get_time_ms();
    for (int i = 0; i < nnz2_16; i ++)
    {
        A2_16_compress[i] = p8_compression(A2_8[i].center, A2_8[i].radius);
    }
    double end11_16 = get_time_ms();
    printf("A2_16_compress completed\n");

    double start11_32 = get_time_ms();
    for (int i = 0; i < nnz2_32; i ++)
    {
        A2_32_compress[i] = p16_compression(A2_16[i].center, A2_16[i].radius);
    }
    double end11_32 = get_time_ms();
    printf("A2_32_compress completed\n");

    double start22_16 = get_time_ms();
    for (int i = 0; i < nnz2_16; i ++)
    {
        A2_8_ext_new[i].center = A2_16_compress[i].posit16;
        A2_8_ext_new[i].radius = p16_read_3r(A2_16_compress[i].posit16);
    }
    double end22_16 = get_time_ms();
    printf("A2_8_ext_new completed\n");

    double start22_32 = get_time_ms();
    for (int i = 0; i < nnz2_32; i ++)
    {
        A2_16_ext_new[i].center = A2_32_compress[i].posit32;
        A2_16_ext_new[i].radius = p32_read_3r(A2_32_compress[i].posit32);
    }
    double end22_32 = get_time_ms();
    printf("A2_16_ext_new completed\n");

    double start33_16 = get_time_ms();
    GS_csr_p16(n, A2_8_ext_new, idx_16, col_id_16, b2_8_extend, x_res2_16);
    double end33_16 = get_time_ms();
    printf("GS_csr_p16 completed\n");

    double start33_32 = get_time_ms();
    GS_csr_p32(n, A2_16_ext_new, idx_32, col_id_32, b2_16_extend, x_res2_32);
    double end33_32 = get_time_ms();
    printf("GS_csr_p32 completed\n");

    FILE *fp_res3 = fopen(res2_16, "a");
    bias = 0.0, dilat = 0.0;
    max_bias = 0.0, max_dilat = 0.0;
    for (int i = 0; i < n; i ++)
    {
        double c1_16 = convertP16ToDouble(x_ref2_16[i].center);
        double r1_16 = convertP16ToDouble(x_ref2_16[i].radius);
        double c2_16 = convertP16ToDouble(x_res2_16[i].center);
        double r2_16 = convertP16ToDouble(x_res2_16[i].radius);
        double temp1 = fabs((c2_16 - c1_16) / c1_16);
        double temp2 = fabs(r2_16 / r1_16);
        bias += temp1;
        dilat += temp2;
        max_bias = fmax(max_bias, temp1);
        max_dilat = fmax(max_dilat, temp2);              
    }
    bias /= n;
    dilat /= n;
    fprintf(fp_res3, "%.17e\t%.17e\t%.17e\t%.17e\t\n", bias, max_bias, dilat, max_dilat);
    fclose(fp_res3);

    FILE *fp_res4 = fopen(res2_32, "a");
    bias = 0.0, dilat = 0.0;
    max_bias = 0.0, max_dilat = 0.0;
    for (int i = 0; i < n; i ++)
    {
        double c1_32 = convertP32ToDouble(x_ref2_32[i].center);
        double r1_32 = convertP32ToDouble(x_ref2_32[i].radius);
        double c2_32 = convertP32ToDouble(x_res2_32[i].center);
        double r2_32 = convertP32ToDouble(x_res2_32[i].radius);
        double temp1 = fabs((c2_32 - c1_32) / c1_32);
        double temp2 = fabs(r2_32 / r1_32);
        bias += temp1;
        dilat += temp2;
        max_bias = fmax(max_bias, temp1);
        max_dilat = fmax(max_dilat, temp2);              
    }
    bias /= n;
    dilat /= n;
    fprintf(fp_res4, "%.17e\t%.17e\t%.17e\t%.17e\t\n", bias, max_bias, dilat, max_dilat);
    fclose(fp_res4);

    FILE *fp_time3 = fopen(time2_16, "a");
    fprintf(fp_time3, "%.17e\t%.17e\t%.17e\t(ms)\n", 
        end11_16 - start11_16, end22_16 - start22_16, end33_16 - start33_16);
    fclose(fp_time3);

    FILE *fp_time4 = fopen(time2_32, "a");
    fprintf(fp_time4, "%.17e\t%.17e\t%.17e\t(ms)\n", 
        end11_32 - start11_32, end22_32 - start22_32, end33_32 - start33_32);
    fclose(fp_time4);

    free(A1_8);
    free(A1_16_compress);
    free(A1_8_extend);
    free(A1_8_ext_new);

    free(A2_8);
    free(A2_16_compress);
    free(A2_8_extend);
    free(A2_8_ext_new);

    free(A1_16);
    free(A1_32_compress);
    free(A1_16_extend);
    free(A1_16_ext_new);

    free(A2_16);
    free(A2_32_compress);
    free(A2_16_extend);
    free(A2_16_ext_new);

    free(idx_16);
    free(idx_32);
    free(col_id_16);
    free(col_id_32);

    free(b1_8);
    free(b1_8_extend);
    free(b1_16);
    free(b1_16_extend);

    free(b2_8);
    free(b2_8_extend);
    free(b2_16);
    free(b2_16_extend);

    free(x_ref1_16);
    free(x_ref1_32);
    free(x_res1_16);
    free(x_res1_32);

    free(x_ref2_16);
    free(x_ref2_32);
    free(x_res2_16);
    free(x_res2_32);

    return 0;
}

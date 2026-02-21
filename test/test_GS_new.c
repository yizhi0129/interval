#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <stdbool.h>
#include <fenv.h>

#include "convert.h"
#include "functions.h"
#include "mpfr_interval.h"



void generate_x_ieee(C_R *x, int n)
{
    const double r_exp_min = -52.0;
    const double r_exp_max = -5.0; 

    srand(get_time_ms());

    for (int i = 0; i < n; i ++)
    {
        x[i].center = 1.0;
        double r_power2 = r_exp_min + ((double)rand() / RAND_MAX) * (r_exp_max - r_exp_min);
        x[i].radius = pow(2, r_power2) * fabs(x[i].center);
    }
}

void generate_matrix_CSR_ieee(C_R **A, int *idx, int **col_id, int n)
{
    const double r_exp_min = -52.0;
    const double r_exp_max = -5.0; 
    
    idx[0] = 0;
    int pos = 0;

    srand(get_time_ms());

    // suppose no more than 5 non-zero entries per row
    int max_nnz = 5 * n;

    C_R *A_tmp = (C_R *)malloc(max_nnz * sizeof(C_R));
    int *col_id_tmp = (int *)malloc(max_nnz * sizeof(int));

    for (int i = 0; i < n; i ++)
    {
        for (int j = 0; j < n; j ++)
        {
            if (j == i)
            {
                A_tmp[pos].center = 4.0;
                double r_power = r_exp_min + ((double)rand() / RAND_MAX) * (r_exp_max - r_exp_min);
                A_tmp[pos].radius = pow(2, r_power) * fabs(A_tmp[pos].center);
                col_id_tmp[pos] = j;
                pos ++;
            }
            else if ((double)rand() / RAND_MAX < (1.0 / n))
            {
                A_tmp[pos].center = -1.0;
                double r_power = r_exp_min + ((double)rand() / RAND_MAX) * (r_exp_max - r_exp_min);
                A_tmp[pos].radius = pow(2, r_power) * fabs(A_tmp[pos].center);
                col_id_tmp[pos] = j;
                pos ++;
            }
        }
        idx[i + 1] = pos;
    }
    *A = (C_R *)realloc(A_tmp, pos * sizeof(C_R));
    *col_id = (int *)realloc(col_id_tmp, pos * sizeof(int));
}

// A * x = b 
void generate_rhs_ieee(int n, const C_R *x, C_R *b, const C_R *A, const int *idx, const int *col_id)
{
    for (int i = 0; i < n; i ++) 
    {
        C_R sum = (C_R){0.0, 0.0};
        for (int j = idx[i]; j < idx[i + 1]; j++) 
        {
            int k = col_id[j];
            sum = interval_add(sum, interval_mult(A[j], x[k]));
        }
        b[i] = sum;
    }
}

int main(int argc, char** argv)
{
    int n = atoi(argv[1]);

    char res[64];
    sprintf(res, "GS_csr_new_res_%d.txt", n);


    C_R *A = NULL;
    int *idx = malloc((n + 1) * sizeof(int));
    int *col_id = NULL;
    C_R *b = malloc(n * sizeof(C_R));
    C_R *x = malloc(n * sizeof(C_R));
    C_R *x1 = malloc(n * sizeof(C_R));
    C_R *x2 = malloc(n * sizeof(C_R));

    generate_x_ieee(x, n);
    generate_matrix_CSR_ieee(&A, idx, &col_id, n);
    generate_rhs_ieee(n, x, b, A, idx, col_id);

    for (int i = 0; i < n; i ++)
    {
        x1[i] = (C_R) {1.0, 1.0};
        x2[i] = (C_R) {1.0, 1.0};
    }

    interval_GS_CSR(A, idx, col_id, b, x1, n);

    int nnz = idx[n] - idx[0];
    FP_INT *A_comp = malloc(nnz * sizeof(FP_INT));
    
    for (int i = 0; i < nnz; i ++)
    {
        A_comp[i] = CR_FP1_adjbis(A[i]);
    }


    for (int i = 0; i < nnz; i ++)
    {
        A[i] = FP_CR3(A_comp[i]);
    }
    
    interval_GS_CSR(A, idx, col_id, b, x2, n);

    fesetround(FE_UPWARD);
    for (int i = 0; i < n; i ++)
    {
        double sup1 = x1[i].center + x1[i].radius;
        double sup2 = x2[i].center + x2[i].radius;
        if (sup1 > sup2)
        {
            fprintf(stderr, "Sup[%d] not included! sup1 = %.17e, sup2 = %.17e\n", i, sup1, sup2);
        }
    }
    fesetround(FE_DOWNWARD);
    for (int i = 0; i < n; i ++)
    {
        double inf1 = x1[i].center - x1[i].radius;
        double inf2 = x2[i].center - x2[i].radius;
        if (inf1 < inf2)
        {
            fprintf(stderr, "Inf[%d] not included! inf1 = %.17e, inf2 = %.17e\n", i, inf1, inf2);
        }
    }
    fesetround(FE_TONEAREST);

    FILE *fp_res = fopen(res, "a");
    double dilat_r1 = 0.0, max_dilat_r1 = 0.0;
    double dilat_e1 = 0.0, max_dilat_e1 = 0.0;
    double dilat_r2 = 0.0, max_dilat_r2 = 0.0;
    double dilat_e2 = 0.0, max_dilat_e2 = 0.0;

    for (int i = 0; i < n; i ++)
    {
        double temp01 = fabs(x1[i].radius / x[i].radius);
        double temp02 = fabs((x1[i].radius / x1[i].center) / (x[i].radius / x[i].center));

        dilat_r1 += temp01;
        dilat_e1 += temp02;

        max_dilat_r1 = fmax(temp01, max_dilat_r1);
        max_dilat_e1 = fmax(temp02, max_dilat_e1);

        double temp1 = fabs(x2[i].radius / x1[i].radius);
        double temp2 = fabs((x2[i].radius / x2[i].center) / (x1[i].radius / x1[i].center));

        dilat_r2 += temp1;
        dilat_e2 += temp2;

        max_dilat_r2 = fmax(temp1, max_dilat_r2);
        max_dilat_e2 = fmax(temp2, max_dilat_e2);
    }
    dilat_r1 /= n;
    dilat_e1 /= n;
    dilat_r2 /= n;
    dilat_e2 /= n;
    fprintf(fp_res, "%.17e\t%.17e\t%.17e\t%.17e\n", dilat_r1, max_dilat_r1, dilat_e1, max_dilat_e1);
    fprintf(fp_res, "%.17e\t%.17e\t%.17e\t%.17e\n", dilat_r2, max_dilat_r2, dilat_e2, max_dilat_e2);
    fclose(fp_res);


    free(A);
    free(A_comp);
    free(idx);
    free(col_id); 
    free(b);
    free(x1);
    free(x2);

    return 0;
}
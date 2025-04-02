#include "convert.h"
#include <time.h>
#include <sys/time.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>


#define N 1000

double get_time_ms() 
{
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return (double)tv.tv_sec * 1000.0 + (double)tv.tv_usec / 1000.0;
}

int main(int argc, char **argv)
{
    srand((unsigned int)time(NULL));

    C_R *test_int = malloc(2 * N * sizeof(C_R));
    INF_SUP *test_int2 = malloc(2 * N * sizeof(INF_SUP));

    FP_INT *res1 = malloc(2 * N * sizeof(FP_INT));
    FP_INT *res2 = malloc(2 * N * sizeof(FP_INT));
    FP_INT *res3 = malloc(2 * N * sizeof(FP_INT));
    FP_INT *res4 = malloc(2 * N * sizeof(FP_INT));

    C_R *res11 = malloc(2 * N * sizeof(C_R));
    C_R *res22 = malloc(2 * N * sizeof(C_R));
    C_R *res33 = malloc(2 * N * sizeof(C_R));
    C_R *res44 = malloc(2 * N * sizeof(C_R));

    FP_INT *res111 = malloc(2 * N * sizeof(FP_INT));
    FP_INT *res222 = malloc(2 * N * sizeof(FP_INT));
    FP_INT *res333 = malloc(2 * N * sizeof(FP_INT));
    FP_INT *res444 = malloc(2 * N * sizeof(FP_INT));

    C_R *res1111 = malloc(2 * N * sizeof(C_R));
    C_R *res2222 = malloc(2 * N * sizeof(C_R));
    C_R *res3333 = malloc(2 * N * sizeof(C_R));
    C_R *res4444 = malloc(2 * N * sizeof(C_R));

    if (!test_int || !res1 || !res2 || !res3 || !res4 || !res11 || !res22 || !res33 || !res44) 
    {
        perror("malloc failed!\n");
        return 1;
    }

    const double r_ratio_min = 1e-20;
    const double r_ratio_max = 1e-3;

    for (int i = 0; i < N; i ++) 
    {
        double c1 = ((double)rand() / RAND_MAX) * 1e5  * (2 * i + 1);   
        double r_ratio1 = (rand() % 10 != 0) ? r_ratio_min + ((double)rand() / RAND_MAX) * (r_ratio_max - r_ratio_min) : 1e-10;
        double r1 = r_ratio1 * fabs(c1);    
        test_int[i].center = c1;
        test_int[i].radius = r1;

        double c2 =  - ((double)rand() / RAND_MAX) * 1e5  * (2 * i + 1);   
        double r_ratio2 = (rand() % 10 != 0) ? r_ratio_min + ((double)rand() / RAND_MAX) * (r_ratio_max - r_ratio_min) : 1e-10;
        double r2 = r_ratio2 * fabs(c2);     
        test_int[i + N].center = c2;
        test_int[i + N].radius = r2;

        double inf1 = ((double)rand() / RAND_MAX) * 1e5  * (2 * i + 1);
        double sup1 = inf1 + 2 * r1;
        test_int2[i].inf = inf1;
        test_int2[i].sup = sup1;

        double sup2 = - ((double)rand() / RAND_MAX) * 1e5  * (2 * i + 1);
        double inf2 = sup2 - 2 * r2;
        test_int2[i + N].inf = inf2;    
        test_int2[i + N].sup = sup2;
    }


    double start1 = get_time_ms();
    for (int i = 0; i < 2 * N; i ++)
    {
        res1[i] = CR_FP1(test_int[i]);
    }
    double end1 = get_time_ms();

    double start2 = get_time_ms();
    for (int i = 0; i < 2 * N; i ++)
    {
        res2[i] = CR_FP2(test_int[i]);
    }
    double end2 = get_time_ms();

    double start3 = get_time_ms();
    for (int i = 0; i < 2 * N; i ++)
    {
        res3[i] = CR_FP3(test_int[i]);
    }
    double end3 = get_time_ms();

    double start4 = get_time_ms();
    for (int i = 0; i < 2 * N; i ++)
    {
        res4[i] = CR_FP4(test_int[i]);
    }
    double end4 = get_time_ms();

    double start5 = get_time_ms();
    for (int i = 0; i < 2 * N; i ++)
    {
        res11[i] = FP_CR(res1[i]);
    }
    double end5 = get_time_ms();


    for (int i = 0; i < 2 * N; i ++)
    {
        res22[i] = FP_CR(res2[i]);
        res33[i] = FP_CR(res3[i]);
        res44[i] = FP_CR(res4[i]);
    }
  
    double start111 = get_time_ms();
    for (int i = 0; i < 2 * N; i ++)
    {
        res111[i] = IS_FP1(test_int2[i]);
    }
    double end111 = get_time_ms();

    double start222 = get_time_ms();
    for (int i = 0; i < 2 * N; i ++)
    {
        res222[i] = IS_FP2(test_int2[i]);
    }
    double end222 = get_time_ms();

    double start333 = get_time_ms();
    for (int i = 0; i < 2 * N; i ++)
    {
        res333[i] = IS_FP3(test_int2[i]);
    }
    double end333 = get_time_ms();

    double start444 = get_time_ms();
    for (int i = 0; i < 2 * N; i ++)
    {
        res444[i] = IS_FP4(test_int2[i]);
    }
    double end444 = get_time_ms();

    for (int i = 0; i < 2 * N; i ++)
    {
        res1111[i] = FP_CR(res111[i]);
        res2222[i] = FP_CR(res222[i]);
        res3333[i] = FP_CR(res333[i]);
        res4444[i] = FP_CR(res444[i]);
    }

    FILE *fp = fopen("convert_time.txt", "a"); 
    if (!fp) {
        perror("Failed to open convert_time.txt");
        return 1;
    }
    fprintf(fp, "%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t(ms)\n", 
        end1 - start1, end2 - start2, end3 - start3, end4 - start4, end5 - start5, 
        end111 - start111, end222 - start222, end333 - start333, end444 - start444); 
    fclose(fp);  

    FILE *fp2 = fopen("CR_FPINT_res.txt", "a");
    if (!fp2) {
        perror("Failed to open convert_res.txt");
        return 1;
    }

    for (int i = 0; i < 2 * N; i ++)
    {
        fprintf(fp2, "test_c[%d] = %lf\t", i, test_int[i].center);
        fprintf(fp2, "test_r[%d] = %lf\n", i, test_int[i].radius);
        fprintf(fp2, "res1[%d] = %lf\t", i, res1[i]);
        fprintf(fp2, "res1_r[%d] = %lf\n", i, res11[i].radius);
        fprintf(fp2, "res2[%d] = %lf\t", i, res2[i]);
        fprintf(fp2, "res2_r[%d] = %lf\n", i, res22[i].radius);
        fprintf(fp2, "res3[%d] = %lf\t", i, res3[i]);
        fprintf(fp2, "res3_r[%d] = %lf\n", i, res33[i].radius);
        fprintf(fp2, "res4[%d] = %lf\t", i, res4[i]);
        fprintf(fp2, "res4_r[%d] = %lf\n", i, res44[i].radius);
        for (int j = 51; j >= 0; j --)
        {
            int b0 = read_m_bit(test_int[i].center, j);
            fprintf(fp2, "%d", b0);
        }
        fprintf(fp2, "\n");
        for (int j = 51; j >= 0; j --)
        {
            int b1 = read_m_bit(res1[i], j);
            fprintf(fp2, "%d", b1);
        }
        fprintf(fp2, "\n");
        for (int j = 51; j >= 0; j --)
        {
            int b2 = read_m_bit(res2[i], j);
            fprintf(fp2, "%d", b2);
        }
        fprintf(fp2, "\n");
        for (int j = 51; j >= 0; j --)
        {
            int b3 = read_m_bit(res3[i], j);
            fprintf(fp2, "%d", b3);
        }
        fprintf(fp2, "\n");
        for (int j = 51; j >= 0; j --)
        {
            int b4 = read_m_bit(res4[i], j);
            fprintf(fp2, "%d", b4);
        }
        fprintf(fp2, "\n");       
    }
    fclose(fp2);

    FILE *fp3 = fopen("IS_FPINT_res.txt", "a");
    if (!fp3) {
        perror("Failed to open IS_FPINT_res.txt");
        return 1;
    }

    for (int i = 0; i < 2 * N; i ++)
    {
        fprintf(fp3, "test_inf[%d] = %lf\t", i, test_int2[i].inf);
        fprintf(fp3, "test_sup[%d] = %lf\n", i, test_int2[i].sup);
        fprintf(fp3, "res1_inf[%d] = %lf\t", i, res111[i] - res1111[i].radius);
        fprintf(fp3, "res1_sup[%d] = %lf\n", i, res111[i] + res1111[i].radius);
        fprintf(fp3, "res2_inf[%d] = %lf\t", i, res222[i] - res2222[i].radius);
        fprintf(fp3, "res2_sup[%d] = %lf\n", i, res222[i] + res2222[i].radius);
        fprintf(fp3, "res3_inf[%d] = %lf\t", i, res333[i] - res3333[i].radius);
        fprintf(fp3, "res3_sup[%d] = %lf\n", i, res333[i] + res3333[i].radius);
        fprintf(fp3, "res4_inf[%d] = %lf\t", i, res444[i] - res4444[i].radius);
        fprintf(fp3, "res4_sup[%d] = %lf\n", i, res444[i] + res4444[i].radius);
        fprintf(fp3, "\n");
        for (int j = 51; j >= 0; j --)
        {
            int b1 = read_m_bit(res111[i], j);
            fprintf(fp3, "%d", b1);
        }
        fprintf(fp3, "\n");
        for (int j = 51; j >= 0; j --)
        {
            int b2 = read_m_bit(res222[i], j);
            fprintf(fp3, "%d", b2);
        }
        fprintf(fp3, "\n");
        for (int j = 51; j >= 0; j --)
        {
            int b3 = read_m_bit(res333[i], j);
            fprintf(fp3, "%d", b3);
        }
        fprintf(fp3, "\n");
        for (int j = 51; j >= 0; j --)
        {
            int b4 = read_m_bit(res444[i], j);
            fprintf(fp3, "%d", b4);
        }
        fprintf(fp3, "\n");       
    }
    fclose(fp3);

    free(test_int);
    free(res1);
    free(res2);
    free(res3);
    free(res4);
    free(res11);
    free(res22);
    free(res33);
    free(res44);
    free(res111);
    free(res222);
    free(res333);
    free(res444);
    free(res1111);
    free(res2222);
    free(res3333);
    free(res4444);

    return 0;
}

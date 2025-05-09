#include "functions.h"
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h> 
#include <math.h>
#include <sys/time.h>
#include <fenv.h>


double get_time_ms() 
{
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return (double)tv.tv_sec * 1000.0 + (double)tv.tv_usec / 1000.0;
}

void fprint_binary(FILE *fp, double x)
{
    union ieee754_double u;
    u.d = x;
    int exp = u.ieee.exponent - DOUBLE_ULS;
    int sign = get_sign_bit(x);
    (sign == 1) ? fprintf(fp, "- ") : fprintf(fp, "+ ");
    fprintf(fp, "1.");
    for (int i = DOUBLE_E - 1; i >= 0; i --)
    {
        int b = read_m_bit(x, i);
        fprintf(fp, "%d", b);
    }
    fprintf(fp, " 2^%d\n", exp);
}

int main(int argc, char **argv)
{
    int N = atoi(argv[1]);
    if (N <= 0)
    {
        printf("Please input a positive integer.\n");
        return 0;
    }

    C_R *test = (C_R*)malloc(sizeof(C_R) * N);
    FP_INT *res = (FP_INT *)malloc(sizeof(FP_INT) * N);
    C_R *ref = (C_R *)malloc(sizeof(C_R) * N);

    char file1[64];
    char file2[64];

    sprintf(file1, "newton_res_%d.txt", N);
    sprintf(file2, "newton_time_%d.txt", N);

    srand(get_time_ms());
    for (int i = 0; i < N; i ++)
    {
        double inf = 0.2 + 1.2 * ((double)rand() / RAND_MAX); // [0.2, 1.4]
        double sup = 1.5 + 0.8 * ((double)rand() / RAND_MAX); // [1.5, 2.3]
        fesetround(FE_UPWARD);
        test[i].center = (inf + sup) / 2;
        test[i].radius = test[i].center - inf;
        fesetround(FE_TONEAREST);
    }

    for (int i = 0; i < N; i ++)
    {
        ref[i] = newton_pr(test[i]);
    }

    double start1 = get_time_ms();
    for (int i = 0; i < N; i ++)
    {
        ref[i] = newton_res(test[i]);
    }
    double end1 = get_time_ms();

    double start2 = get_time_ms();
    for (int i = 0; i < N; i ++)
    {
        res[i] = CR_FP1_adjbis(ref[i]);
    }
    double end2 = get_time_ms();
    
    FILE *fp1 = fopen(file1, "a");
    fprintf(fp1, "sqrt(2) =\t");
    fprint_binary(fp1, sqrt(2));
    for (int i = 0; i < N; i ++)
    {
        fprintf(fp1, "res[%d] =\t", i);
        fprint_binary(fp1, res[i]);
    }
    fclose(fp1);

    FILE *fp2 = fopen(file2, "a");
    fprintf(fp2, "%.10e\t%.10e\tms\n", end1 - start1, end2 - start2);
    fclose(fp2);

    free(test);
    free(res);
    free(ref);

    return 0;
}

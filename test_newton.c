#include "functions.h"
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h> 
#include <math.h>
#include <sys/time.h>
#include <fenv.h>


#define N 100

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
    C_R *test = (C_R*)malloc(sizeof(C_R) * N);
    FP_INT *res = (FP_INT *)malloc(sizeof(FP_INT) * N);
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

    double start = get_time_ms();
    for (int i = 0; i < N; i ++)
    {
        res[i] = newton(test[i]);
    }
    double end = get_time_ms();
    
    FILE *fp = fopen("newton_res.txt", "a");
    fprintf(fp, "sqrt(2) =\t");
    fprint_binary(fp, sqrt(2));
    for (int i = 0; i < N; i ++)
    {
        fprintf(fp, "res[%d] =\t", i);
        fprint_binary(fp, res[i]);
    }
    fclose(fp);

    FILE *fp2 = fopen("newton_time.txt", "a");
    fprintf(fp2, "%.10e ms\n", end - start);
    fclose(fp2);

    free(test);
    free(res);

    return 0;
}

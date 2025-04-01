#include "convert.h"
#include <time.h>
#include <sys/time.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>


#define N 50

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
    FP_INT *res1 = malloc(2 * N * sizeof(FP_INT));
    FP_INT *res2 = malloc(2 * N * sizeof(FP_INT));
    if (!test_int || !res1 || !res2)
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
    }


    double start1 = get_time_ms();
    for (int i = 0; i < 2 * N; i ++)
    {
        res1[i] = CR_FP1(test_int[i]);
        printf("res1[%d] = %lf\n", i, res1[i]);
    }
    double end1 = get_time_ms();

    double start2 = get_time_ms();
    for (int i = 0; i < 2 * N; i ++)
    {
        res2[i] = CR_FP2(test_int[i]);
        printf("res2[%d] = %lf\n", i, res2[i]);
    }
    double end2 = get_time_ms();

    FILE *fp = fopen("convert_time.txt", "a"); 
    if (!fp) {
        perror("Failed to open convert_time.txt");
        return 1;
    }
    fprintf(fp, "%.6e\t%.6e\t(ms)\n", end1 - start1, end2 - start2); 
    fclose(fp);  

    free(test_int);
    free(res1);
    free(res2);

    return 0;
}

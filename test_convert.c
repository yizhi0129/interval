#include "convert.h"
#include <time.h>
#include <sys/time.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>


#define N 20000

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
    C_R *test_int = malloc(2 * N * sizeof(C_R));
    INF_SUP *test_int2 = malloc(2 * N * sizeof(INF_SUP));

    FP_INT *CR_c1 = malloc(2 * N * sizeof(FP_INT));
    FP_INT *CR_c2 = malloc(2 * N * sizeof(FP_INT));
    FP_INT *CR_c3 = malloc(2 * N * sizeof(FP_INT));
    FP_INT *CR_c4 = malloc(2 * N * sizeof(FP_INT));
    FP_INT *CR_c5 = malloc(2 * N * sizeof(FP_INT));
    FP_INT *CR_c6 = malloc(2 * N * sizeof(FP_INT));

    C_R *CR_cr1 = malloc(2 * N * sizeof(C_R));
    C_R *CR_cr2 = malloc(2 * N * sizeof(C_R));
    C_R *CR_cr3 = malloc(2 * N * sizeof(C_R));
    C_R *CR_cr4 = malloc(2 * N * sizeof(C_R));
    C_R *CR_cr5 = malloc(2 * N * sizeof(C_R));
    C_R *CR_cr6 = malloc(2 * N * sizeof(C_R));

    FP_INT *IS_c1 = malloc(2 * N * sizeof(FP_INT));
    FP_INT *IS_c2 = malloc(2 * N * sizeof(FP_INT));
    FP_INT *IS_c3 = malloc(2 * N * sizeof(FP_INT));
    FP_INT *IS_c4 = malloc(2 * N * sizeof(FP_INT));
    FP_INT *IS_c5 = malloc(2 * N * sizeof(FP_INT));
    FP_INT *IS_c6 = malloc(2 * N * sizeof(FP_INT));

    C_R *IS_cr1 = malloc(2 * N * sizeof(C_R));
    C_R *IS_cr2 = malloc(2 * N * sizeof(C_R));
    C_R *IS_cr3 = malloc(2 * N * sizeof(C_R));
    C_R *IS_cr4 = malloc(2 * N * sizeof(C_R));
    C_R *IS_cr5 = malloc(2 * N * sizeof(C_R));
    C_R *IS_cr6 = malloc(2 * N * sizeof(C_R));

    const double c_exp_min = -6.0;
    const double c_exp_max = 9.0;
    const double r_exp_min = -10.0;
    const double r_exp_max = -3.0;

    srand(get_time_ms());    

    for (int i = 0; i < N; i ++) 
    {
        double c_power1 = c_exp_min + ((double)rand() / RAND_MAX) * (c_exp_max - c_exp_min);
        double c1 = ((double)rand() / RAND_MAX) * pow(10, c_power1);   
        double r_power1 = r_exp_min + ((double)rand() / RAND_MAX) * (r_exp_max - r_exp_min);
        double r1 = pow(10, r_power1) * fabs(c1);    
        test_int[i].center = c1;
        test_int[i].radius = r1;

        double c_power2 = c_exp_min + ((double)rand() / RAND_MAX) * (c_exp_max - c_exp_min);
        double c2 =  - ((double)rand() / RAND_MAX) * pow(10, c_power2);   
        double r_power2 = r_exp_min + ((double)rand() / RAND_MAX) * (r_exp_max - r_exp_min);
        double r2 = pow(10, r_power2) * fabs(c2);     
        test_int[i + N].center = c2;
        test_int[i + N].radius = r2;

        double inf1 = (0.5 * (double)rand() / RAND_MAX + 0.8) * c1;
        double sup1 = inf1 + (1 + 2 * (double)rand() / RAND_MAX) * r1;
        test_int2[i].inf = inf1;
        test_int2[i].sup = sup1;

        double sup2 = (0.5 * (double)rand() / RAND_MAX + 0.8) * c2;
        double inf2 = sup2 - (1 + 2 * (double)rand() / RAND_MAX) * r2;
        test_int2[i + N].inf = inf2;    
        test_int2[i + N].sup = sup2;
    }


    double start1 = get_time_ms();
    for (int i = 0; i < 2 * N; i ++)
    {
        CR_c1[i] = CR_FP1(test_int[i]);
    }
    double end1 = get_time_ms();

    double start2 = get_time_ms();
    for (int i = 0; i < 2 * N; i ++)
    {
        CR_c2[i] = CR_FP2(test_int[i]);
    }
    double end2 = get_time_ms();

    double start3 = get_time_ms();
    for (int i = 0; i < 2 * N; i ++)
    {
        CR_c3[i] = CR_FP3(test_int[i]);
    }
    double end3 = get_time_ms();

    double start4 = get_time_ms();
    for (int i = 0; i < 2 * N; i ++)
    {
        CR_c4[i] = CR_FP4(test_int[i]);
    }
    double end4 = get_time_ms();

    double start5 = get_time_ms();
    for (int i = 0; i < 2 * N; i ++)
    {
        CR_c5[i] = CR_FP5(test_int[i]);
    }
    double end5 = get_time_ms();

    double start6 = get_time_ms();
    for (int i = 0; i < 2 * N; i ++)
    {
        CR_c6[i] = CR_FP6(test_int[i]);
    }
    double end6 = get_time_ms();
  
    double start111 = get_time_ms();
    for (int i = 0; i < 2 * N; i ++)
    {
        IS_c1[i] = IS_FP1(test_int2[i]);
    }
    double end111 = get_time_ms();

    double start222 = get_time_ms();
    for (int i = 0; i < 2 * N; i ++)
    {
        IS_c2[i] = IS_FP2(test_int2[i]);
    }
    double end222 = get_time_ms();

    double start333 = get_time_ms();
    for (int i = 0; i < 2 * N; i ++)
    {
        IS_c3[i] = IS_FP3(test_int2[i]);
    }
    double end333 = get_time_ms();

    double start444 = get_time_ms();
    for (int i = 0; i < 2 * N; i ++)
    {
        IS_c4[i] = IS_FP4(test_int2[i]);
    }
    double end444 = get_time_ms();

    double start555 = get_time_ms();
    for (int i = 0; i < 2 * N; i ++)
    {
        IS_c5[i] = IS_FP5(test_int2[i]);
    }
    double end555 = get_time_ms();

    double start666 = get_time_ms();
    for (int i = 0; i < 2 * N; i ++)
    {
        IS_c6[i] = IS_FP6(test_int2[i]);
    }
    double end666 = get_time_ms();

    double start0 = get_time_ms();
    for (int i = 0; i < 2 * N; i ++)
    {
        CR_cr1[i] = FP_CR(CR_c1[i]);
        CR_cr2[i] = FP_CR(CR_c2[i]);
        CR_cr3[i] = FP_CR(CR_c3[i]);
        CR_cr4[i] = FP_CR(CR_c4[i]);
        CR_cr5[i] = FP_CR(CR_c5[i]);
        CR_cr6[i] = FP_CR(CR_c6[i]);
        IS_cr1[i] = FP_CR(IS_c1[i]);
        IS_cr2[i] = FP_CR(IS_c2[i]);
        IS_cr3[i] = FP_CR(IS_c3[i]);
        IS_cr4[i] = FP_CR(IS_c4[i]);
        IS_cr5[i] = FP_CR(IS_c5[i]);
        IS_cr6[i] = FP_CR(IS_c6[i]);
    }
    double end0 = get_time_ms();

    FILE *fp = fopen("convert_time.txt", "a"); 
    if (!fp) 
    {
        perror("Failed to open convert_time.txt");
        return 1;
    }
    fprintf(fp, "%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t(ms)\n", 
        end1 - start1, end2 - start2, end3 - start3, end4 - start4, (end0 - start0) / 12.0, 
        end111 - start111, end222 - start222, end333 - start333, end444 - start444, end5 - start5, end555 - start555, end6 - start6, end666 - start666); 
    fclose(fp);  

    FILE *fp2 = fopen("CR_FPINT_res.txt", "a");
    if (!fp2) 
    {
        perror("Failed to open convert_res.txt");
        return 1;
    }

    for(int i = 0; i < 2 * N; i += 1000)
    {
        fprintf(fp2, "test_int[%d] = [%.10e,\t%.10e]\n", i, test_int[i].center - test_int[i].radius, test_int[i].center + test_int[i].radius);
        fprint_binary(fp2, test_int[i].center);
        fprint_binary(fp2, test_int[i].radius);
        fprintf(fp2, "CR1:\n");
        fprint_binary(fp2, CR_c1[i]);
        fprint_binary(fp2, CR_cr1[i].center);
        fprint_binary(fp2, CR_cr1[i].radius);
        fprintf(fp2, "CR2:\n");
        fprint_binary(fp2, CR_c2[i]);
        fprint_binary(fp2, CR_cr2[i].center);
        fprint_binary(fp2, CR_cr2[i].radius);
        fprintf(fp2, "CR3:\n");
        fprint_binary(fp2, CR_c3[i]);
        fprint_binary(fp2, CR_cr3[i].center);
        fprint_binary(fp2, CR_cr3[i].radius);
        fprintf(fp2, "CR4:\n");
        fprint_binary(fp2, CR_c4[i]);
        fprint_binary(fp2, CR_cr4[i].center);
        fprint_binary(fp2, CR_cr4[i].radius);
        fprintf(fp2, "CR5:\n");
        fprint_binary(fp2, CR_c5[i]);
        fprint_binary(fp2, CR_cr5[i].center);
        fprint_binary(fp2, CR_cr5[i].radius);
        fprintf(fp2, "CR6:\n");
        fprint_binary(fp2, CR_c6[i]);
        fprint_binary(fp2, CR_cr6[i].center);
        fprint_binary(fp2, CR_cr6[i].radius);
        fprintf(fp2, "\n");
    }

    for (int i = 0; i < 2 * N; i ++)
    {   
        if (CR_c1[i] - CR_cr1[i].radius > test_int[i].center - test_int[i].radius || 
            CR_c1[i] + CR_cr1[i].radius < test_int[i].center + test_int[i].radius)
        {
            fprintf(fp2, "CR_c1 = [%.10e,\t%.10e]\t", CR_c1[i] - CR_cr1[i].radius, CR_c1[i] + CR_cr1[i].radius);
            fprintf(fp2, "Error\n");
            fprint_binary(fp2, test_int[i].center);
            fprint_binary(fp2, test_int[i].radius);
            fprint_binary(fp2, CR_c1[i]);
            fprint_binary(fp2, CR_cr1[i].center);
            fprint_binary(fp2, CR_cr1[i].radius);
        }
        
        if (CR_c2[i] - CR_cr2[i].radius > test_int[i].center - test_int[i].radius || 
            CR_c2[i] + CR_cr2[i].radius < test_int[i].center + test_int[i].radius)
        {
            fprintf(fp2, "CR_c2 = [%.10e,\t%.10e]\t", CR_c2[i] - CR_cr2[i].radius, CR_c2[i] + CR_cr2[i].radius);
            fprintf(fp2, "Error\n");
            fprint_binary(fp2, test_int[i].center);
            fprint_binary(fp2, test_int[i].radius);
            fprint_binary(fp2, CR_c2[i]);
            fprint_binary(fp2, CR_cr2[i].center);
            fprint_binary(fp2, CR_cr2[i].radius);
        }
        
        if (CR_c3[i] - CR_cr3[i].radius > test_int[i].center - test_int[i].radius || 
            CR_c3[i] + CR_cr3[i].radius < test_int[i].center + test_int[i].radius)
        {
            fprintf(fp2, "CR_c3 = [%.10e,\t%.10e]\t", CR_c3[i] - CR_cr3[i].radius, CR_c3[i] + CR_cr3[i].radius);
            fprintf(fp2, "Error\n");
            fprint_binary(fp2, test_int[i].center);
            fprint_binary(fp2, test_int[i].radius);
            fprint_binary(fp2, CR_c3[i]);
            fprint_binary(fp2, CR_cr3[i].center);
            fprint_binary(fp2, CR_cr3[i].radius);
        }
        
        if (CR_c4[i] - CR_cr4[i].radius > test_int[i].center - test_int[i].radius || 
            CR_c4[i] + CR_cr4[i].radius < test_int[i].center + test_int[i].radius)
        {
            fprintf(fp2, "CR_c4 = [%.10e,\t%.10e]\t", CR_c4[i] - CR_cr4[i].radius, CR_c4[i] + CR_cr4[i].radius);
            fprintf(fp2, "Error\n");
            fprint_binary(fp2, test_int[i].center);
            fprint_binary(fp2, test_int[i].radius);
            fprint_binary(fp2, CR_c4[i]);
            fprint_binary(fp2, CR_cr4[i].center);
            fprint_binary(fp2, CR_cr4[i].radius);
        }
        
        if (CR_c5[i] - CR_cr5[i].radius > test_int[i].center - test_int[i].radius || 
            CR_c5[i] + CR_cr5[i].radius < test_int[i].center + test_int[i].radius)
        {
            fprintf(fp2, "CR_c5 = [%.10e,\t%.10e]\t", CR_c5[i] - CR_cr5[i].radius, CR_c5[i] + CR_cr5[i].radius);
            fprintf(fp2, "Error\n");
            fprint_binary(fp2, test_int[i].center);
            fprint_binary(fp2, test_int[i].radius);
            fprint_binary(fp2, CR_c5[i]);
            fprint_binary(fp2, CR_cr5[i].center);
            fprint_binary(fp2, CR_cr5[i].radius);
        }
        
        if (CR_c6[i] - CR_cr6[i].radius > test_int[i].center - test_int[i].radius || 
            CR_c6[i] + CR_cr6[i].radius < test_int[i].center + test_int[i].radius)
        {
            fprintf(fp2, "CR_c6 = [%.10e,\t%.10e]\t", CR_c6[i] - CR_cr6[i].radius, CR_c6[i] + CR_cr6[i].radius);
            fprintf(fp2, "Error\n");
            fprint_binary(fp2, test_int[i].center);
            fprint_binary(fp2, test_int[i].radius);
            fprint_binary(fp2,CR_c6[i]);
            fprint_binary(fp2, CR_cr6[i].center);
            fprint_binary(fp2, CR_cr6[i].radius);
        }
    }
    fclose(fp2);

    FILE *fp3 = fopen("IS_FPINT_res.txt", "a");
    if (!fp3) 
    {
        perror("Failed to open IS_FPINT_res.txt");
        return 1;
    }

    for (int i = 0; i < 2 * N; i += 1000)
    {
        fprintf(fp3, "test_int2[%d] = [%.10e,\t%.10e]\n", i, test_int2[i].inf, test_int2[i].sup);
        fprint_binary(fp3, test_int2[i].inf);
        fprint_binary(fp3, test_int2[i].sup);
        fprintf(fp3, "IS1:\n");
        fprint_binary(fp3, IS_c1[i]);
        fprint_binary(fp3, IS_cr1[i].center);
        fprint_binary(fp3, IS_cr1[i].radius);
        fprintf(fp3, "IS2:\n");
        fprint_binary(fp3, IS_c2[i]);
        fprint_binary(fp3, IS_cr2[i].center);
        fprint_binary(fp3, IS_cr2[i].radius);
        fprintf(fp3, "IS3:\n");
        fprint_binary(fp3, IS_c3[i]);
        fprint_binary(fp3, IS_cr3[i].center);
        fprint_binary(fp3, IS_cr3[i].radius);
        fprintf(fp3, "IS4:\n");
        fprint_binary(fp3, IS_c4[i]);
        fprint_binary(fp3, IS_cr4[i].center);
        fprint_binary(fp3, IS_cr4[i].radius);
        fprintf(fp3, "IS5:\n");
        fprint_binary(fp3, IS_c5[i]);
        fprint_binary(fp3, IS_cr5[i].center);
        fprint_binary(fp3, IS_cr5[i].radius);
        fprintf(fp3, "IS6:\n");
        fprint_binary(fp3, IS_c6[i]);
        fprint_binary(fp3, IS_cr6[i].center);
        fprint_binary(fp3, IS_cr6[i].radius);
    }

    for (int i = 0; i < 2 * N; i ++)
    {     
        if (IS_c1[i] - IS_cr1[i].radius > test_int2[i].inf || 
            IS_c1[i] + IS_cr1[i].radius < test_int2[i].sup)
        {
            fprintf(fp3, "IS_c1 = [%.10e,\t%.10e]\t", IS_c1[i] - IS_cr1[i].radius, IS_c1[i] + IS_cr1[i].radius);
            fprintf(fp3, "Error\n");
            fprint_binary(fp3, IS_c1[i]);
            fprint_binary(fp3, IS_cr1[i].center);
            fprint_binary(fp3, IS_cr1[i].radius);
        }
        
        if (IS_c2[i] - IS_cr2[i].radius > test_int2[i].inf ||
            IS_c2[i] + IS_cr2[i].radius < test_int2[i].sup)
        {
            fprintf(fp3, "IS_c2 = [%.10e,\t%.10e]\t", IS_c2[i] - IS_cr2[i].radius, IS_c2[i] + IS_cr2[i].radius);
            fprintf(fp3, "Error\n");
            fprint_binary(fp3, IS_c2[i]);
            fprint_binary(fp3, IS_cr2[i].center);
            fprint_binary(fp3, IS_cr2[i].radius);
        }
        
        if (IS_c3[i] - IS_cr3[i].radius > test_int2[i].inf ||
            IS_c3[i] + IS_cr3[i].radius < test_int2[i].sup)
        {
            fprintf(fp3, "IS_c3 = [%.10e,\t%.10e]\t", IS_c3[i] - IS_cr3[i].radius, IS_c3[i] + IS_cr3[i].radius);
            fprintf(fp3, "Error\n");
            fprint_binary(fp3, IS_c3[i]);
            fprint_binary(fp3, IS_cr3[i].center);
            fprint_binary(fp3, IS_cr3[i].radius);
        }
        
        if (IS_c4[i] - IS_cr4[i].radius > test_int2[i].inf || 
            IS_c4[i] + IS_cr4[i].radius < test_int2[i].sup)
        {
            fprintf(fp3, "IS_c4 = [%.10e,\t%.10e]\t", IS_c4[i] - IS_cr4[i].radius, IS_c4[i] + IS_cr4[i].radius);
            fprintf(fp3, "Error\n");
            fprint_binary(fp3, IS_c4[i]);
            fprint_binary(fp3, IS_cr4[i].center);
            fprint_binary(fp3, IS_cr4[i].radius);
        }
        
        if (IS_c5[i] - IS_cr5[i].radius > test_int2[i].inf || 
            IS_c5[i] + IS_cr5[i].radius < test_int2[i].sup)
        {
            fprintf(fp3, "IS_c5 = [%.10e,\t%.10e]\t", IS_c5[i] - IS_cr5[i].radius, IS_c5[i] + IS_cr5[i].radius);
            fprintf(fp3, "Error\n");
            fprint_binary(fp3, IS_c5[i]);
            fprint_binary(fp3, IS_cr5[i].center);
            fprint_binary(fp3, IS_cr5[i].radius);
        }
        
        if (IS_c6[i] - IS_cr6[i].radius > test_int2[i].inf || 
            IS_c6[i] + IS_cr6[i].radius < test_int2[i].sup)
        {
            fprintf(fp3, "IS_c6 = [%.10e,\t%.10e]\t", IS_c6[i] - IS_cr6[i].radius, IS_c6[i] + IS_cr6[i].radius);
            fprintf(fp3, "Error\n");
            fprint_binary(fp3, IS_c6[i]);
            fprint_binary(fp3, IS_cr6[i].center);
            fprint_binary(fp3, IS_cr6[i].radius);
        }
    }
    fclose(fp3);

    free(test_int);

    free(CR_c1);
    free(CR_c2);
    free(CR_c3);
    free(CR_c4);
    free(CR_c5);

    free(CR_cr1);
    free(CR_cr2);
    free(CR_cr3);
    free(CR_cr4);
    free(CR_cr5);

    free(IS_c1);
    free(IS_c2);
    free(IS_c3);
    free(IS_c4);
    free(IS_c5);

    free(IS_cr1);
    free(IS_cr2);
    free(IS_cr3);
    free(IS_cr4);
    free(IS_cr5);

    return 0;
}

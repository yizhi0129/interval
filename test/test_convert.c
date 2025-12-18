#include "convert.h"

#include <stdlib.h>
#include <stdio.h>
#include <math.h>


#define N 20000

#define MAX_RATE 1e+300


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
    FP_INT *CR_c1b = malloc(2 * N * sizeof(FP_INT));
    FP_INT *CR_c1_adj = malloc(2 * N * sizeof(FP_INT));
    FP_INT *CR_c1_adjbis = malloc(2 * N * sizeof(FP_INT));
    FP_INT *CR_c1_adjter = malloc(2 * N * sizeof(FP_INT));

    C_R *CR_cr1 = malloc(2 * N * sizeof(C_R));
    C_R *CR_cr2 = malloc(2 * N * sizeof(C_R));
    C_R *CR_cr3 = malloc(2 * N * sizeof(C_R));
    C_R *CR_cr4 = malloc(2 * N * sizeof(C_R));
    C_R *CR_cr5 = malloc(2 * N * sizeof(C_R));
    C_R *CR_cr6 = malloc(2 * N * sizeof(C_R));
    C_R *CR_cr1b = malloc(2 * N * sizeof(C_R));
    C_R *CR_cr1_adj = malloc(2 * N * sizeof(C_R));
    C_R *CR_cr1_adjbis = malloc(2 * N * sizeof(C_R));
    C_R *CR_cr1_adjter = malloc(2 * N * sizeof(C_R));

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

    double e1_max = 0;
    double e1_min = MAX_RATE;
    double e2_max = 0;
    double e2_min = MAX_RATE;
    double e3_max = 0;
    double e3_min = MAX_RATE;

    double r1_max = 0;
    double r1_min = MAX_RATE;
    double r2_max = 0;
    double r2_min = MAX_RATE;
    double r3_max = 0;
    double r3_min = MAX_RATE;

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
        double c2 = - ((double)rand() / RAND_MAX) * pow(10, c_power2);   
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

    double start1b = get_time_ms();
    for (int i = 0; i < 2 * N; i ++)
    {
        CR_c1b[i] = CR_FP1b(test_int[i]);
    }
    double end1b = get_time_ms();

    double start1_adj = get_time_ms();
    for (int i = 0; i < 2 * N; i ++)
    {
        CR_c1_adj[i] = CR_FP1_adj(test_int[i]);
    }
    double end1_adj = get_time_ms();

    double start1_adjbis = get_time_ms();
    for (int i = 0; i < 2 * N; i ++)
    {
        CR_c1_adjbis[i] = CR_FP1_adjbis(test_int[i]);
    }
    double end1_adjbis = get_time_ms();

    double start1_adjter = get_time_ms();
    for (int i = 0; i < 2 * N; i ++)
    {
        CR_c1_adjter[i] = CR_FP1_adjter(test_int[i]);
    }
    double end1_adjter = get_time_ms();
  
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
        CR_cr1b[i] = FP_CR(CR_c1b[i]);
        CR_cr1_adj[i] = FP_CR(CR_c1_adj[i]);
        CR_cr1_adjbis[i] = FP_CR(CR_c1_adjbis[i]);
        CR_cr1_adjter[i] = FP_CR(CR_c1_adjter[i]);
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
    fprintf(fp, "%.17e\t%.17e\t%.17e\t%.17e\t%.17e\t%.17e\t%.17e\t%.17e\t%.17e\t%.17e\t%.17e\t%.17e\t%.17e\t%.17e\t%.17e\t%.17e\t%.17e\t(ms)\n", 
        end1 - start1, end2 - start2, end3 - start3, end4 - start4, (end0 - start0) / 16.0, 
        end111 - start111, end222 - start222, end333 - start333, end444 - start444, end5 - start5, end555 - start555, 
        end6 - start6, end666 - start666, end1b - start1b, end1_adj - start1_adj, 
        end1_adjbis - start1_adjbis, end1_adjter - start1_adjter); 
    fclose(fp);  

    FILE *fp2 = fopen("CR_FPINT_res.txt", "a");
    if (!fp2) 
    {
        perror("Failed to open convert_res.txt");
        return 1;
    }

    FILE *fpdil = fopen("dilat.txt", "a");
    for(int i = 0; i < 2 * N; i += 1000)
    {
        fprintf(fp2, "test_int[%d] = [%.17e,\t%.17e]\n", i, test_int[i].center - test_int[i].radius, test_int[i].center + test_int[i].radius);
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
        fprintf(fp2, "CR1b:\n");
        fprint_binary(fp2, CR_c1b[i]);
        fprint_binary(fp2, CR_cr1b[i].center);
        fprint_binary(fp2, CR_cr1b[i].radius);
        fprintf(fp2, "\n");
        fprintf(fp2, "CR1adj:\n");
        fprint_binary(fp2, CR_c1_adj[i]);
        fprint_binary(fp2, CR_cr1_adj[i].center);
        fprint_binary(fp2, CR_cr1_adj[i].radius);
        fprintf(fp2, "\n");
        fprintf(fp2, "CR1adjbis:\n");
        fprint_binary(fp2, CR_c1_adjbis[i]);
        fprint_binary(fp2, CR_cr1_adjbis[i].center);
        fprint_binary(fp2, CR_cr1_adjbis[i].radius);
        fprintf(fp2, "\n");
        fprintf(fp2, "CR1adjter:\n");
        fprint_binary(fp2, CR_c1_adjter[i]);
        fprint_binary(fp2, CR_cr1_adjter[i].center);
        fprint_binary(fp2, CR_cr1_adjter[i].radius);
        fprintf(fp2, "\n");

        double r1 = 2 * CR_cr1_adj[i].radius / test_int[i].radius;
        double r2 = 3 * CR_cr1_adjbis[i].radius / test_int[i].radius;
        double r3 = 4 * CR_cr1_adjter[i].radius / test_int[i].radius;

        double e1 = test_int[i].radius / test_int[i].center;
        double e1_rate = 2 * CR_cr1_adj[i].radius / CR_cr1_adj[i].center / e1;
        double e2_rate = 3 * CR_cr1_adjbis[i].radius / CR_cr1_adjbis[i].center / e1;
        double e3_rate = 4 * CR_cr1_adjter[i].radius / CR_cr1_adjter[i].center / e1;

        r1_max = fmax(r1_max, r1);
        r1_min = fmin(r1_min, r1);
        r2_max = fmax(r2_max, r2);
        r2_min = fmin(r2_min, r2);
        r3_max = fmax(r3_max, r3);
        r3_min = fmin(r3_min, r3);
        e1_max = fmax(e1_max, e1_rate);
        e1_min = fmin(e1_min, e1_rate);
        e2_max = fmax(e2_max, e2_rate);
        e2_min = fmin(e2_min, e2_rate);
        e3_max = fmax(e3_max, e3_rate);
        e3_min = fmin(e3_min, e3_rate);
    }

    fprintf(fpdil, "%.17e\t%.17e\t%.17e\t%.17e\t%.17e\t%.17e\t%.17e\t%.17e\t%.17e\t%.17e\t%.17e\t%.17e\n", 
        e1_max, e1_min, r1_max, r1_min, e2_max, e2_min, r2_max, r2_min, e3_max, e3_min, r3_max, r3_min);
    fclose(fpdil);

    for (int i = 0; i < 2 * N; i ++)
    {   
        if (CR_c1[i] - CR_cr1[i].radius > test_int[i].center - test_int[i].radius || 
            CR_c1[i] + CR_cr1[i].radius < test_int[i].center + test_int[i].radius)
        {
            fprintf(fp2, "CR_c1 = [%.17e,\t%.17e]\t", CR_c1[i] - CR_cr1[i].radius, CR_c1[i] + CR_cr1[i].radius);
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
            fprintf(fp2, "CR_c2 = [%.17e,\t%.17e]\t", CR_c2[i] - CR_cr2[i].radius, CR_c2[i] + CR_cr2[i].radius);
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
            fprintf(fp2, "CR_c3 = [%.17e,\t%.17e]\t", CR_c3[i] - CR_cr3[i].radius, CR_c3[i] + CR_cr3[i].radius);
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
            fprintf(fp2, "CR_c4 = [%.17e,\t%.17e]\t", CR_c4[i] - CR_cr4[i].radius, CR_c4[i] + CR_cr4[i].radius);
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
            fprintf(fp2, "CR_c5 = [%.17e,\t%.17e]\t", CR_c5[i] - CR_cr5[i].radius, CR_c5[i] + CR_cr5[i].radius);
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
            fprintf(fp2, "CR_c6 = [%.17e,\t%.17e]\t", CR_c6[i] - CR_cr6[i].radius, CR_c6[i] + CR_cr6[i].radius);
            fprintf(fp2, "Error\n");
            fprint_binary(fp2, test_int[i].center);
            fprint_binary(fp2, test_int[i].radius);
            fprint_binary(fp2,CR_c6[i]);
            fprint_binary(fp2, CR_cr6[i].center);
            fprint_binary(fp2, CR_cr6[i].radius);
        }

        if (CR_c1b[i] - CR_cr1b[i].radius > test_int[i].center - test_int[i].radius || 
            CR_c1b[i] + CR_cr1b[i].radius < test_int[i].center + test_int[i].radius)
        {
            fprintf(fp2, "CR_c1b = [%.17e,\t%.17e]\t", CR_c1b[i] - CR_cr1b[i].radius, CR_c1b[i] + CR_cr1b[i].radius);
            fprintf(fp2, "Error\n");
            fprint_binary(fp2, test_int[i].center);
            fprint_binary(fp2, test_int[i].radius);
            fprint_binary(fp2, CR_c1b[i]);
            fprint_binary(fp2, CR_cr1b[i].center);
            fprint_binary(fp2, CR_cr1b[i].radius);
        }

        if (CR_c1_adj[i] - 2 * CR_cr1_adj[i].radius > test_int[i].center - test_int[i].radius || 
            CR_c1_adj[i] + 2 * CR_cr1_adj[i].radius < test_int[i].center + test_int[i].radius)
        {
            fprintf(fp2, "CR_c1_adj = [%.17e,\t%.17e]\t", CR_c1_adj[i] - CR_cr1_adj[i].radius, CR_c1_adj[i] + CR_cr1_adj[i].radius);
            fprintf(fp2, "Error\n");
            fprint_binary(fp2, test_int[i].center);
            fprint_binary(fp2, test_int[i].radius);
            fprint_binary(fp2, CR_c1_adj[i]);
            fprint_binary(fp2, CR_cr1_adj[i].center);
            fprint_binary(fp2, CR_cr1_adj[i].radius);
        }

        if (CR_c1_adjbis[i] - 3 * CR_cr1_adjbis[i].radius > test_int[i].center - test_int[i].radius || 
            CR_c1_adjbis[i] + 3 * CR_cr1_adjbis[i].radius < test_int[i].center + test_int[i].radius)
        {
            fprintf(fp2, "CR_c1_adjbis = [%.17e,\t%.17e]\t", CR_c1_adjbis[i] - CR_cr1_adjbis[i].radius, CR_c1_adjbis[i] + CR_cr1_adjbis[i].radius);
            fprintf(fp2, "Error\n");
            fprint_binary(fp2, test_int[i].center);
            fprint_binary(fp2, test_int[i].radius);
            fprint_binary(fp2, CR_c1_adjbis[i]);
            fprint_binary(fp2, CR_cr1_adjbis[i].center);
            fprint_binary(fp2, CR_cr1_adjbis[i].radius);
        }

        if (CR_c1_adjter[i] - 4 * CR_cr1_adjter[i].radius > test_int[i].center - test_int[i].radius || 
            CR_c1_adjter[i] + 4 * CR_cr1_adjter[i].radius < test_int[i].center + test_int[i].radius)
        {
            fprintf(fp2, "CR_c1_adjter = [%.17e,\t%.17e]\t", CR_c1_adjter[i] - CR_cr1_adjter[i].radius, CR_c1_adjter[i] + CR_cr1_adjter[i].radius);
            fprintf(fp2, "Error\n");
            fprint_binary(fp2, test_int[i].center);
            fprint_binary(fp2, test_int[i].radius);
            fprint_binary(fp2, CR_c1_adjter[i]);
            fprint_binary(fp2, CR_cr1_adjter[i].center);
            fprint_binary(fp2, CR_cr1_adjter[i].radius);
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
        fprintf(fp3, "test_int2[%d] = [%.17e,\t%.17e]\n", i, test_int2[i].inf, test_int2[i].sup);
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
            fprintf(fp3, "IS_c1 = [%.17e,\t%.17e]\t", IS_c1[i] - IS_cr1[i].radius, IS_c1[i] + IS_cr1[i].radius);
            fprintf(fp3, "Error\n");
            fprint_binary(fp3, IS_c1[i]);
            fprint_binary(fp3, IS_cr1[i].center);
            fprint_binary(fp3, IS_cr1[i].radius);
        }
        
        if (IS_c2[i] - IS_cr2[i].radius > test_int2[i].inf ||
            IS_c2[i] + IS_cr2[i].radius < test_int2[i].sup)
        {
            fprintf(fp3, "IS_c2 = [%.17e,\t%.17e]\t", IS_c2[i] - IS_cr2[i].radius, IS_c2[i] + IS_cr2[i].radius);
            fprintf(fp3, "Error\n");
            fprint_binary(fp3, IS_c2[i]);
            fprint_binary(fp3, IS_cr2[i].center);
            fprint_binary(fp3, IS_cr2[i].radius);
        }
        
        if (IS_c3[i] - IS_cr3[i].radius > test_int2[i].inf ||
            IS_c3[i] + IS_cr3[i].radius < test_int2[i].sup)
        {
            fprintf(fp3, "IS_c3 = [%.17e,\t%.17e]\t", IS_c3[i] - IS_cr3[i].radius, IS_c3[i] + IS_cr3[i].radius);
            fprintf(fp3, "Error\n");
            fprint_binary(fp3, IS_c3[i]);
            fprint_binary(fp3, IS_cr3[i].center);
            fprint_binary(fp3, IS_cr3[i].radius);
        }
        
        if (IS_c4[i] - IS_cr4[i].radius > test_int2[i].inf || 
            IS_c4[i] + IS_cr4[i].radius < test_int2[i].sup)
        {
            fprintf(fp3, "IS_c4 = [%.17e,\t%.17e]\t", IS_c4[i] - IS_cr4[i].radius, IS_c4[i] + IS_cr4[i].radius);
            fprintf(fp3, "Error\n");
            fprint_binary(fp3, IS_c4[i]);
            fprint_binary(fp3, IS_cr4[i].center);
            fprint_binary(fp3, IS_cr4[i].radius);
        }
        
        if (IS_c5[i] - IS_cr5[i].radius > test_int2[i].inf || 
            IS_c5[i] + IS_cr5[i].radius < test_int2[i].sup)
        {
            fprintf(fp3, "IS_c5 = [%.17e,\t%.17e]\t", IS_c5[i] - IS_cr5[i].radius, IS_c5[i] + IS_cr5[i].radius);
            fprintf(fp3, "Error\n");
            fprint_binary(fp3, IS_c5[i]);
            fprint_binary(fp3, IS_cr5[i].center);
            fprint_binary(fp3, IS_cr5[i].radius);
        }
        
        if (IS_c6[i] - IS_cr6[i].radius > test_int2[i].inf || 
            IS_c6[i] + IS_cr6[i].radius < test_int2[i].sup)
        {
            fprintf(fp3, "IS_c6 = [%.17e,\t%.17e]\t", IS_c6[i] - IS_cr6[i].radius, IS_c6[i] + IS_cr6[i].radius);
            fprintf(fp3, "Error\n");
            fprint_binary(fp3, IS_c6[i]);
            fprint_binary(fp3, IS_cr6[i].center);
            fprint_binary(fp3, IS_cr6[i].radius);
        }
    }
    fclose(fp3);

    free(test_int);
    free(test_int2);

    free(CR_c1);
    free(CR_c2);
    free(CR_c3);
    free(CR_c4);
    free(CR_c5);
    free(CR_c6);
    free(CR_c1b);
    free(CR_c1_adj);
    free(CR_c1_adjbis);
    free(CR_c1_adjter);

    free(CR_cr1);
    free(CR_cr2);
    free(CR_cr3);
    free(CR_cr4);
    free(CR_cr5);
    free(CR_cr6);
    free(CR_cr1b);
    free(CR_cr1_adj);
    free(CR_cr1_adjbis);
    free(CR_cr1_adjter);

    free(IS_c1);
    free(IS_c2);
    free(IS_c3);
    free(IS_c4);
    free(IS_c5);
    free(IS_c6);

    free(IS_cr1);
    free(IS_cr2);
    free(IS_cr3);
    free(IS_cr4);
    free(IS_cr5);
    free(IS_cr6);

    return 0;
}

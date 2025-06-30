#include "functions.h"

#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <fenv.h>

#define MAX_RATE 1e+300

// old conversion

int main(int argc, char **argv)
{
    int n = atoi(argv[1]);
    
    int N = n * n;
    char file1[64];
    char file2[64];
    char res11[64];
    char res12[64];
    char res13[64];
    char res14[64];
    char res21[64];
    char res22[64];
    char res23[64];
    char res24[64];
    char dil1[64];
    char dil2[64];

    sprintf(file1, "matprod_time1_%d.txt", n);
    sprintf(file2, "matprod_time2_%d.txt", n);
    sprintf(res11, "C1_1_%d.txt", n);
    sprintf(res12, "C2_1_%d.txt", n);
    sprintf(res13, "C3_1_%d.txt", n);
    sprintf(res14, "C4_1_%d.txt", n);
    sprintf(res21, "C1_2_%d.txt", n);
    sprintf(res22, "C2_2_%d.txt", n);
    sprintf(res23, "C3_2_%d.txt", n);
    sprintf(res24, "C4_2_%d.txt", n);
    sprintf(dil1, "dil_1_%d.txt", n);
    sprintf(dil2, "dil_2_%d.txt", n);

    double *mAp = malloc(N * sizeof(double));
    double *rAp = malloc(N * sizeof(double));
    double *mBp = malloc(N * sizeof(double));
    double *rBp = malloc(N * sizeof(double));

    double *iAp = malloc(N * sizeof(double));
    double *iBp = malloc(N * sizeof(double));
    double *sAp = malloc(N * sizeof(double));
    double *sBp = malloc(N * sizeof(double));

    double *mAm = malloc(N * sizeof(double));
    double *rAm = malloc(N * sizeof(double));
    double *mBm = malloc(N * sizeof(double));
    double *rBm = malloc(N * sizeof(double));

    double *iAm = malloc(N * sizeof(double));
    double *iBm = malloc(N * sizeof(double));
    double *sAm = malloc(N * sizeof(double));
    double *sBm = malloc(N * sizeof(double));

    double *mC1 = malloc(N * sizeof(double));
    double *rC1 = malloc(N * sizeof(double));
    double *mC1_tilde = malloc(N * sizeof(double));
    double *rC1_tilde = malloc(N * sizeof(double));
    C_R *C1 = malloc(N * sizeof(C_R));
    double *iC1 = malloc(N * sizeof(double));
    double *sC1 = malloc(N * sizeof(double));

    double *mC2 = malloc(N * sizeof(double));
    double *rC2 = malloc(N * sizeof(double));
    double *mC2_tilde = malloc(N * sizeof(double));
    double *rC2_tilde = malloc(N * sizeof(double));
    C_R *C2 = malloc(N * sizeof(C_R));
    double *iC2 = malloc(N * sizeof(double));
    double *sC2 = malloc(N * sizeof(double));

    double *mC3 = malloc(N * sizeof(double));
    double *rC3 = malloc(N * sizeof(double));
    double *mC3_tilde = malloc(N * sizeof(double));
    double *rC3_tilde = malloc(N * sizeof(double));
    C_R *C3 = malloc(N * sizeof(C_R));
    double *iC3 = malloc(N * sizeof(double));
    double *sC3 = malloc(N * sizeof(double));

    double *mC4 = malloc(N * sizeof(double));
    double *rC4 = malloc(N * sizeof(double));
    double *mC4_tilde = malloc(N * sizeof(double));
    double *rC4_tilde = malloc(N * sizeof(double));
    C_R *C4 = malloc(N * sizeof(C_R));
    double *iC4 = malloc(N * sizeof(double));
    double *sC4 = malloc(N * sizeof(double));

    const double c_exp_min = -2.0;
    const double c_exp_max = 2.0;
    const double r_exp_min = -10.0;
    const double r_exp_max = -6.0;

    srand(get_time_ms());

    for (int i = 0; i < N; i ++)
    {
        double c_power1 = c_exp_min + ((double)rand() / RAND_MAX) * (c_exp_max - c_exp_min);
        double c1 = ((double)rand() / RAND_MAX) * pow(10, c_power1);   
        double r_power1 = r_exp_min + ((double)rand() / RAND_MAX) * (r_exp_max - r_exp_min);
        double r1 = pow(10, r_power1) * fabs(c1);  

        double c_power2 = c_exp_min + ((double)rand() / RAND_MAX) * (c_exp_max - c_exp_min);
        double c2 =  - ((double)rand() / RAND_MAX) * pow(10, c_power2);   
        double r_power2 = r_exp_min + ((double)rand() / RAND_MAX) * (r_exp_max - r_exp_min);
        double r2 = pow(10, r_power2) * fabs(c2);  

        double c_power3 = c_exp_min + ((double)rand() / RAND_MAX) * (c_exp_max - c_exp_min);
        double c3 = ((double)rand() / RAND_MAX) * pow(10, c_power3);   
        double r_power3 = r_exp_min + ((double)rand() / RAND_MAX) * (r_exp_max - r_exp_min);
        double r3 = pow(10, r_power3) * fabs(c3);  

        double c_power4 = c_exp_min + ((double)rand() / RAND_MAX) * (c_exp_max - c_exp_min);
        double c4 =  - ((double)rand() / RAND_MAX) * pow(10, c_power4);   
        double r_power4 = r_exp_min + ((double)rand() / RAND_MAX) * (r_exp_max - r_exp_min);
        double r4 = pow(10, r_power4) * fabs(c4);

        iAp[i] = c1;
        sAp[i] = c1 + r1;
        sAm[i] = c2;
        iAm[i] = c2 - r2;
        iBp[i] = c3;
        sBp[i] = c3 + r3;
        sBm[i] = c4;
        iBm[i] = c4 - r4;
    }

    ref1_pp(iAp, sAp, iBp, sBp, iC1, sC1, n);
    ref2_mm(iAm, sAm, iBm, sBm, iC2, sC2, n);
    ref3_pm(iAp, sAp, iBm, sBm, iC3, sC3, n);
    ref4_mp(iAm, sAm, iBp, sBp, iC4, sC4, n);

    mat_is_cr(iAp, sAp, iBp, sBp, mAp, rAp, mBp, rBp, n);
    mat_is_cr(iAm, sAm, iBm, sBm, mAm, rAm, mBm, rBm, n);

    // calculation
    double start2 = get_time_ms();
    int_mat_mult(mAp, rAp, mBp, rBp, mC1, rC1, n);
    double end2 = get_time_ms();

    double start3 = get_time_ms();
    int_mat_mult(mAm, rAm, mBm, rBm, mC2, rC2, n);
    double end3 = get_time_ms();
        
    double start4 = get_time_ms();
    int_mat_mult(mAp, rAp, mBm, rBm, mC3, rC3, n);
    double end4 = get_time_ms();

    double start5 = get_time_ms();
    int_mat_mult(mAm, rAm, mBp, rBp, mC4, rC4, n);
    double end5 = get_time_ms();

    // restore
    for (int i = 0; i < N; i ++)
    {
        C1[i] = (C_R){mC1[i], rC1[i]};
        C2[i] = (C_R){mC2[i], rC2[i]};
        C3[i] = (C_R){mC3[i], rC3[i]};
        C4[i] = (C_R){mC4[i], rC4[i]};
    }
    // conversion
    double start6 = get_time_ms();
    for (int i = 0; i < N; i ++)
    {
        mC1_tilde[i] = CR_FP1(C1[i]);
        mC2_tilde[i] = CR_FP1(C2[i]);
        mC3_tilde[i] = CR_FP1(C3[i]);
        mC4_tilde[i] = CR_FP1(C4[i]);             
    }
    double end6 = get_time_ms();

    // read r_tilde
    for (int i = 0; i < N; i ++)
    {
        rC1_tilde[i] = FP_CR(mC1_tilde[i]).radius;
        rC2_tilde[i] = FP_CR(mC2_tilde[i]).radius;
        rC3_tilde[i] = FP_CR(mC3_tilde[i]).radius;
        rC4_tilde[i] = FP_CR(mC4_tilde[i]).radius; 
    }

    FILE *result1[4];
    result1[0] = fopen(res11, "a");
    result1[1] = fopen(res12, "a");
    result1[2] = fopen(res13, "a");
    result1[3] = fopen(res14, "a");

    double e1_max = 0.0, e2_max = 0.0, e3_max = 0.0, e4_max = 0.0;
    double e1_min = MAX_RATE, e2_min = MAX_RATE, e3_min = MAX_RATE, e4_min = MAX_RATE;
    double r1_max = 0.0, r2_max = 0.0, r3_max = 0.0, r4_max = 0.0;
    double r1_min = MAX_RATE, r2_min = MAX_RATE, r3_min = MAX_RATE, r4_min = MAX_RATE;

    FILE *f = fopen("err1.txt", "a");
    FILE *freq1 = fopen("freq1.txt", "a");

    FILE *d1 = fopen(dil1, "a");

    for (int i = 0; i < N; i ++)
    {
        double e1 = rC1[i] / mC1[i];
        double e2 = rC2[i] / mC2[i];
        double e3 = - rC3[i] / mC3[i];
        double e4 = - rC4[i] / mC4[i];
    
        double e1_tilde = rC1_tilde[i] / mC1_tilde[i];
        double e2_tilde = rC2_tilde[i] / mC2_tilde[i];
        double e3_tilde = - rC3_tilde[i] / mC3_tilde[i];
        double e4_tilde = - rC4_tilde[i] / mC4_tilde[i];

        double e1_rate = e1_tilde / e1;
        double e2_rate = e2_tilde / e2;
        double e3_rate = e3_tilde / e3;
        double e4_rate = e4_tilde / e4;

        double r1_rate = rC1_tilde[i] / rC1[i];
        double r2_rate = rC2_tilde[i] / rC2[i];
        double r3_rate = rC3_tilde[i] / rC3[i];
        double r4_rate = rC4_tilde[i] / rC4[i];

        fprintf(d1, "%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%d\n", e1_rate, r1_rate, 
            e2_rate, r2_rate, e3_rate, r3_rate, e4_rate, r4_rate, i);

        if (r1_rate > 1000.0)
        {
            fprintf(freq1, "n = %d C1[%d]\n", n, i);
            fprint_binary(freq1, mC1[i]);
            fprint_binary(freq1, rC1[i]);
            fprint_binary(freq1, mC1_tilde[i]);
            fprint_binary(freq1, rC1_tilde[i]);
        } 

        if (r2_rate > 1000.0)
        {
            fprintf(freq1, "n = %d C2[%d]\n", n, i);
            fprint_binary(freq1, mC2[i]);
            fprint_binary(freq1, rC2[i]);
            fprint_binary(freq1, mC2_tilde[i]);
            fprint_binary(freq1, rC2_tilde[i]);
        }

        if (r3_rate > 1000.0)
        {
            fprintf(freq1, "n = %d C3[%d]\n", n, i);
            fprint_binary(freq1, mC3[i]);
            fprint_binary(freq1, rC3[i]);
            fprint_binary(freq1, mC3_tilde[i]);
            fprint_binary(freq1, rC3_tilde[i]);
        }

        if (r4_rate > 1000.0)
        {
            fprintf(freq1, "n = %d C4[%d]\n", n, i);
            fprint_binary(freq1, mC4[i]);
            fprint_binary(freq1, rC4[i]);
            fprint_binary(freq1, mC4_tilde[i]);
            fprint_binary(freq1, rC4_tilde[i]);
        }

        // check inclusion
        fesetround(FE_UPWARD);
        if (iC1[i] + rC1[i] < mC1[i] || sC1[i] > mC1[i] + rC1[i])
        {
            fprintf(f, "Error algo n = %d C1[%d]\niC1 = %.10f sC1 = %.10f\nmC1 = %.10f rC1 = %.10f\n",
                 n, i, iC1[i], sC1[i], mC1[i], rC1[i]);
            fprint_binary(f, iC1[i]);
            fprint_binary(f, sC1[i]);
            fprint_binary(f, mC1[i]);
            fprint_binary(f, rC1[i]);
        }
        if (iC2[i] + rC2[i] < mC2[i] || sC2[i] > mC2[i] + rC2[i])
        {
            fprintf(f, "Error algo n = %d C2[%d]\niC2 = %.10f sC2 = %.10f\nmC2 = %.10f rC2 = %.10f\n", 
                 n, i, iC2[i], sC2[i], mC2[i], rC2[i]);
            fprint_binary(f, iC2[i]);
            fprint_binary(f, sC2[i]);
            fprint_binary(f, mC2[i]);
            fprint_binary(f, rC2[i]);
        }
        if (iC3[i] + rC3[i] < mC3[i] || sC3[i] > mC3[i] + rC3[i] || mC3_tilde[i] > iC3[i] + rC3_tilde[i] || mC3_tilde[i] + rC3_tilde[i] < sC3[i])
        {
            fprintf(f, "Error algo n = %d C3[%d]\niC3 = %.10f sC3 = %.10f\nmC3 = %.10f rC3 = %.10f\n",
                 n, i, iC3[i], sC3[i], mC3[i], rC3[i]);
            fprint_binary(f, iC3[i]);
            fprint_binary(f, sC3[i]);
            fprint_binary(f, mC3[i]);
            fprint_binary(f, rC3[i]);
        }
        if (iC4[i] + rC4[i] < mC4[i] || sC4[i] > mC4[i] + rC4[i] || mC4_tilde[i] > iC4[i] + rC4_tilde[i] || mC4_tilde[i] + rC4_tilde[i] < sC4[i])
        {
            fprintf(f, "Error algo n = %d C4[%d]\niC4 = %.10f sC4 = %.10f\nmC4 = %.10f rC4 = %.10f\n",
                 n, i, iC4[i], sC4[i], mC4[i], rC4[i]);
            fprint_binary(f, iC4[i]);
            fprint_binary(f, sC4[i]);
            fprint_binary(f, mC4[i]);
            fprint_binary(f, rC4[i]);
        }
        fesetround(FE_TONEAREST);

        e1_max = fmax(e1_max, e1_rate);
        e1_min = fmin(e1_min, e1_rate);
        e2_max = fmax(e2_max, e2_rate);
        e2_min = fmin(e2_min, e2_rate);
        e3_max = fmax(e3_max, e3_rate);
        e3_min = fmin(e3_min, e3_rate);
        e4_max = fmax(e4_max, e4_rate);
        e4_min = fmin(e4_min, e4_rate);

        r1_max = fmax(r1_max, r1_rate);
        r1_min = fmin(r1_min, r1_rate);
        r2_max = fmax(r2_max, r2_rate);
        r2_min = fmin(r2_min, r2_rate);
        r3_max = fmax(r3_max, r3_rate);
        r3_min = fmin(r3_min, r3_rate);
        r4_max = fmax(r4_max, r4_rate);
        r4_min = fmin(r4_min, r4_rate);
    }

    fclose(d1);
    fclose(freq1);

    fprintf(result1[0], "%.6e\t%.6e\t%.6e\t%.6e\n", e1_max, e1_min, r1_max, r1_min);
    fprintf(result1[1], "%.6e\t%.6e\t%.6e\t%.6e\n", e2_max, e2_min, r2_max, r2_min);
    fprintf(result1[2], "%.6e\t%.6e\t%.6e\t%.6e\n", e3_max, e3_min, r3_max, r3_min);
    fprintf(result1[3], "%.6e\t%.6e\t%.6e\t%.6e\n", e4_max, e4_min, r4_max, r4_min);

    fclose(result1[0]);
    fclose(result1[1]);
    fclose(result1[2]);
    fclose(result1[3]);
    fclose(f);
        
    FILE *fp1 = fopen(file1, "a");
    fprintf(fp1, "%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t(ms)\n", 
        end2 - start2, end3 - start3, end4 - start4, end5 - start5, 0.25 * (end6 - start6));
    fclose(fp1);

    for (int i = 0; i < N; i ++)
    {
        double c_power1 = c_exp_min + ((double)rand() / RAND_MAX) * (c_exp_max - c_exp_min);
        double c1 = ((double)rand() / RAND_MAX) * pow(10, c_power1);   

        double c_power2 = c_exp_min + ((double)rand() / RAND_MAX) * (c_exp_max - c_exp_min);
        double c2 =  - ((double)rand() / RAND_MAX) * pow(10, c_power2);   

        double c_power3 = c_exp_min + ((double)rand() / RAND_MAX) * (c_exp_max - c_exp_min);
        double c3 = ((double)rand() / RAND_MAX) * pow(10, c_power3);   

        double c_power4 = c_exp_min + ((double)rand() / RAND_MAX) * (c_exp_max - c_exp_min);
        double c4 =  - ((double)rand() / RAND_MAX) * pow(10, c_power4);   

        mAp[i] = c1;
        rAp[i] = 2 * FP_CR(c1).radius;
        mAm[i] = c2;
        rAm[i] = 2 * FP_CR(c2).radius;
        mBp[i] = c3;
        rBp[i] = 2 * FP_CR(c3).radius;
        mBm[i] = c4;
        rBm[i] = 2 * FP_CR(c4).radius;
    }

    mat_cr_is(mAp, rAp, mBp, rBp, iAp, sAp, iBp, sBp, n);
    mat_cr_is(mAm, rAm, mBm, rBm, iAm, sAm, iBm, sBm, n);

    ref1_pp(iAp, sAp, iBp, sBp, iC1, sC1, n);
    ref2_mm(iAm, sAm, iBm, sBm, iC2, sC2, n);
    ref3_pm(iAp, sAp, iBm, sBm, iC3, sC3, n);
    ref4_mp(iAm, sAm, iBp, sBp, iC4, sC4, n);

    // calculation
    double start7 = get_time_ms();
    int_mat_mult(mAp, rAp, mBp, rBp, mC1, rC1, n);
    double end7 = get_time_ms();

    double start8 = get_time_ms();
    int_mat_mult(mAm, rAm, mBm, rBm, mC2, rC2, n);
    double end8 = get_time_ms();

    double start9 = get_time_ms();
    int_mat_mult(mAp, rAp, mBm, rBm, mC3, rC3, n);
    double end9 = get_time_ms();

    double start10 = get_time_ms();
    int_mat_mult(mAm, rAm, mBp, rBp, mC4, rC4, n);
    double end10 = get_time_ms();

    // restore
    for (int i = 0; i < N; i ++)
    {
        C1[i] = (C_R){mC1[i], rC1[i]};
        C2[i] = (C_R){mC2[i], rC2[i]};
        C3[i] = (C_R){mC3[i], rC3[i]};
        C4[i] = (C_R){mC4[i], rC4[i]};
    }
    // conversion
    double start11 = get_time_ms();
    for (int i = 0; i < N; i ++)
    {
        mC1_tilde[i] = CR_FP1(C1[i]);
        mC2_tilde[i] = CR_FP1(C2[i]);
        mC3_tilde[i] = CR_FP1(C3[i]);
        mC4_tilde[i] = CR_FP1(C4[i]);             
    }
    double end11 = get_time_ms();
    
    // read r_tilde
    for (int i = 0; i < N; i ++)
    {
        rC1_tilde[i] = FP_CR(mC1_tilde[i]).radius;
        rC2_tilde[i] = FP_CR(mC2_tilde[i]).radius;
        rC3_tilde[i] = FP_CR(mC3_tilde[i]).radius;
        rC4_tilde[i] = FP_CR(mC4_tilde[i]).radius; 
    }

    FILE *result2[4];
    result2[0] = fopen(res21, "a");
    result2[1] = fopen(res22, "a");
    result2[2] = fopen(res23, "a");
    result2[3] = fopen(res24, "a");

    e1_max = 0.0, e2_max = 0.0, e3_max = 0.0, e4_max = 0.0;
    e1_min = MAX_RATE, e2_min = MAX_RATE, e3_min = MAX_RATE, e4_min = MAX_RATE;
    r1_max = 0.0, r2_max = 0.0, r3_max = 0.0, r4_max = 0.0;
    r1_min = MAX_RATE, r2_min = MAX_RATE, r3_min = MAX_RATE, r4_min = MAX_RATE;

    FILE *f2 = fopen("err2.txt", "a");
    FILE *freq2 = fopen("freq2.txt", "a");

    FILE *d2 = fopen(dil2, "a");

    for (int i = 0; i < N; i ++)
    {
        double e1 = rC1[i] / mC1[i];
        double e2 = rC2[i] / mC2[i];
        double e3 = - rC3[i] / mC3[i];
        double e4 = - rC4[i] / mC4[i];
    
        double e1_tilde = rC1_tilde[i] / mC1_tilde[i];
        double e2_tilde = rC2_tilde[i] / mC2_tilde[i];
        double e3_tilde = - rC3_tilde[i] / mC3_tilde[i];
        double e4_tilde = - rC4_tilde[i] / mC4_tilde[i];

        double e1_rate = e1_tilde / e1;
        double e2_rate = e2_tilde / e2;
        double e3_rate = e3_tilde / e3;
        double e4_rate = e4_tilde / e4;

        double r1_rate = rC1_tilde[i] / rC1[i];
        double r2_rate = rC2_tilde[i] / rC2[i];
        double r3_rate = rC3_tilde[i] / rC3[i];
        double r4_rate = rC4_tilde[i] / rC4[i];

        fprintf(d2, "%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%d\n", e1_rate, r1_rate, 
            e2_rate, r2_rate, e3_rate, r3_rate, e4_rate, r4_rate, i);

        if (r1_rate > 1000.0)
        {
            fprintf(freq2, "n = %d C1[%d]\n", n, i);
            fprint_binary(freq2, mC1[i]);
            fprint_binary(freq2, rC1[i]);
            fprint_binary(freq2, mC1_tilde[i]);
            fprint_binary(freq2, rC1_tilde[i]);
        }

        if (r2_rate > 1000.0)
        {
            fprintf(freq2, "n = %d C2[%d]\n", n, i);
            fprint_binary(freq2, mC2[i]);
            fprint_binary(freq2, rC2[i]);
            fprint_binary(freq2, mC2_tilde[i]);
            fprint_binary(freq2, rC2_tilde[i]);
        }

        if (r3_rate > 1000.0)
        {
            fprintf(freq2, "n = %d C3[%d]\n", n, i);
            fprint_binary(freq2, mC3[i]);
            fprint_binary(freq2, rC3[i]);
            fprint_binary(freq2, mC3_tilde[i]);
            fprint_binary(freq2, rC3_tilde[i]);
        }

        if (r4_rate > 1000.0)
        {
            fprintf(freq2, "n = %d C4[%d]\n", n, i);
            fprint_binary(freq2, mC4[i]);
            fprint_binary(freq2, rC4[i]);
            fprint_binary(freq2, mC4_tilde[i]);
            fprint_binary(freq2, rC4_tilde[i]);
        }

        // check inclusion
        fesetround(FE_UPWARD);
        if (iC1[i] + rC1[i] < mC1[i] || sC1[i] > mC1[i] + rC1[i])
         {
            fprintf(f2, "Error algo n = %d C1[%d]\niC1 = %.10f sC1 = %.10f\nmC1 = %.10f rC1 = %.10f\n",
                 n, i, iC1[i], sC1[i], mC1[i], rC1[i]);
            fprint_binary(f2, iC1[i]);
            fprint_binary(f2, sC1[i]);
            fprint_binary(f2, mC1[i]);
            fprint_binary(f2, rC1[i]);
        }
        if (iC2[i] + rC2[i] < mC2[i] || sC2[i] > mC2[i] + rC2[i])
        {
            fprintf(f2, "Error algo n = %d C2[%d]\niC2 = %.10f sC2 = %.10f\nmC2 = %.10f rC2 = %.10f\n", 
                 n, i, iC2[i], sC2[i], mC2[i], rC2[i]);
            fprint_binary(f2, iC2[i]);
            fprint_binary(f2, sC2[i]);
            fprint_binary(f2, mC2[i]);
            fprint_binary(f2, rC2[i]);
        }
        if (iC3[i] + rC3[i] < mC3[i] || sC3[i] > mC3[i] + rC3[i])
        {
            fprintf(f2, "Error algo n = %d C3[%d]\niC3 = %.10f sC3 = %.10f\nmC3 = %.10f rC3 = %.10f\n",
                 n, i, iC3[i], sC3[i], mC3[i], rC3[i]);
            fprint_binary(f2, iC3[i]);
            fprint_binary(f2, sC3[i]);
            fprint_binary(f2, mC3[i]);
            fprint_binary(f2, rC3[i]);
        }
        if (iC4[i] + rC4[i] < mC4[i] || sC4[i] > mC4[i] + rC4[i])
        {
            fprintf(f2, "Error algo n = %d C4[%d]\niC4 = %.10f sC4 = %.10f\nmC4 = %.10f rC4 = %.10f\n",
                 n, i, iC4[i], sC4[i], mC4[i], rC4[i]);
            fprint_binary(f2, iC4[i]);
            fprint_binary(f2, sC4[i]);
            fprint_binary(f2, mC4[i]);
            fprint_binary(f2, rC4[i]);
        }
        fesetround(FE_TONEAREST);

        e1_max = fmax(e1_max, e1_rate);
        e1_min = fmin(e1_min, e1_rate);
        e2_max = fmax(e2_max, e2_rate);
        e2_min = fmin(e2_min, e2_rate);
        e3_max = fmax(e3_max, e3_rate);
        e3_min = fmin(e3_min, e3_rate);
        e4_max = fmax(e4_max, e4_rate);
        e4_min = fmin(e4_min, e4_rate);

        r1_max = fmax(r1_max, r1_rate);
        r1_min = fmin(r1_min, r1_rate);
        r2_max = fmax(r2_max, r2_rate);
        r2_min = fmin(r2_min, r2_rate);
        r3_max = fmax(r3_max, r3_rate);
        r3_min = fmin(r3_min, r3_rate);
        r4_max = fmax(r4_max, r4_rate);
        r4_min = fmin(r4_min, r4_rate);
    }

    fclose(d2);
    fclose(freq2);

    fprintf(result2[0], "%.6e\t%.6e\t%.6e\t%.6e\n", e1_max, e1_min, r1_max, r1_min);
    fprintf(result2[1], "%.6e\t%.6e\t%.6e\t%.6e\n", e2_max, e2_min, r2_max, r2_min);
    fprintf(result2[2], "%.6e\t%.6e\t%.6e\t%.6e\n", e3_max, e3_min, r3_max, r3_min);
    fprintf(result2[3], "%.6e\t%.6e\t%.6e\t%.6e\n", e4_max, e4_min, r4_max, r4_min);

    fclose(result2[0]);
    fclose(result2[1]);
    fclose(result2[2]);
    fclose(result2[3]);
    fclose(f2);
        
    FILE *fp2 = fopen(file2, "a");
    fprintf(fp2, "%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t(ms)\n", 
        end7 - start7, end8 - start8, end9 - start9, end10 - start10, 0.25 * (end11 - start11));
    fclose(fp2);

    free(mAp);
    free(rAp);
    free(mBp);
    free(rBp);

    free(mAm);
    free(rAm);
    free(mBm);
    free(rBm);

    free(mC1);
    free(rC1);
    free(mC1_tilde);
    free(rC1_tilde);
    free(mC2);
    free(rC2);
    free(mC2_tilde);
    free(rC2_tilde);
    free(mC3);
    free(rC3);
    free(mC3_tilde);
    free(rC3_tilde);
    free(mC4);
    free(rC4);
    free(mC4_tilde);
    free(rC4_tilde);

    free(C1);
    free(C2);
    free(C3);
    free(C4);

    free(iAp);
    free(iAm);
    free(iBp);
    free(iBm);
    free(sAp);
    free(sAm);
    free(sBp);
    free(sBm);
    free(iC1);
    free(iC2);
    free(iC3);
    free(iC4);
    free(sC1);
    free(sC2);
    free(sC3);
    free(sC4);

    return 0;
}
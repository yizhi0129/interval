#include "functions.h"

#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <fenv.h>

#define MAX_RATE 1e+300
//#define FIXED_RAD 1e-15

// new conversion 3 * uls

/*
double randn() 
{
    static int hasSpare = 0;
    static double spare;
    
    if (hasSpare) 
    {
        hasSpare = 0;
        return spare;
    }

    hasSpare = 1;

    double u, v, s;
    do {
        u = (double)rand() / RAND_MAX * 2.0 - 1.0; // [-1,1]
        v = (double)rand() / RAND_MAX * 2.0 - 1.0;
        s = u * u + v * v;
    } while (s >= 1.0 || s == 0.0);

    s = sqrt(-2.0 * log(s) / s);
    spare = v * s;
    return u * s;
}
*/

int main(int argc, char **argv)
{
    int n = atoi(argv[1]);
    
    int N = n * n;
    char file1[64];
    char file2[64];
    char res1[64];
    char res2[64];
    //char dil1[64];
    //char dil2[64];
    //char test[64];

    sprintf(file1, "matprod_time1_%d.txt", n);
    sprintf(file2, "matprod_time2_%d.txt", n);
    sprintf(res1, "C1_%d.txt", n);
    sprintf(res2, "C2_%d.txt", n);
    //sprintf(dil1, "dil_1_%d.txt", n);
    //sprintf(dil2, "dil_2_%d.txt", n);
    //sprintf(test, "Ct_%d.txt", n);

    double *mAp = malloc(N * sizeof(double));
    double *rAp = malloc(N * sizeof(double));
    double *mBp = malloc(N * sizeof(double));
    double *rBp = malloc(N * sizeof(double));

    double *iAp = malloc(N * sizeof(double));
    double *iBp = malloc(N * sizeof(double));
    double *sAp = malloc(N * sizeof(double));
    double *sBp = malloc(N * sizeof(double));

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

    double *ref_mC1 = malloc(N * sizeof(double));
    double *ref_rC1 = malloc(N * sizeof(double));
    double *ref_mC2 = malloc(N * sizeof(double));
    double *ref_rC2 = malloc(N * sizeof(double));
    
    /*
    double *mAt = malloc(N * sizeof(double));
    double *rAt = malloc(N * sizeof(double));
    double *mBt = malloc(N * sizeof(double));
    double *rBt = malloc(N * sizeof(double));
    double *mCt = malloc(N * sizeof(double));
    double *rCt = malloc(N * sizeof(double));
    */

    const double c_exp_min = -2.0;
    const double c_exp_max = 2.0;
    const double r_exp_min = -17.0;
    const double r_exp_max = -4.0;

    // 2 ^ (-52) ~ 2 * 10 ^ (-16)

    srand(get_time_ms());

    /*
    for (int i = 0; i < N; i ++)
    {
        rAt[i] = FIXED_RAD;
        rBt[i] = FIXED_RAD;
        mAt[i] = randn();
        mBt[i] = randn();
    }
    */

    for (int i = 0; i < N; i ++)
    {
        double c_power1 = c_exp_min + ((double)rand() / RAND_MAX) * (c_exp_max - c_exp_min);
        double c1 = ((double)rand() / RAND_MAX) * pow(10, c_power1);   
        double r_power1 = r_exp_min + ((double)rand() / RAND_MAX) * (r_exp_max - r_exp_min);
        double r1 = pow(10, r_power1) * fabs(c1);  

        /*
        double c_power2 = c_exp_min + ((double)rand() / RAND_MAX) * (c_exp_max - c_exp_min);
        double c2 =  - ((double)rand() / RAND_MAX) * pow(10, c_power2);   
        double r_power2 = r_exp_min + ((double)rand() / RAND_MAX) * (r_exp_max - r_exp_min);
        double r2 = pow(10, r_power2) * fabs(c2);  
        */
        
        double c_power3 = c_exp_min + ((double)rand() / RAND_MAX) * (c_exp_max - c_exp_min);
        double c3 = ((double)rand() / RAND_MAX) * pow(10, c_power3);   
        double r_power3 = r_exp_min + ((double)rand() / RAND_MAX) * (r_exp_max - r_exp_min);
        double r3 = pow(10, r_power3) * fabs(c3);  

        /*
        double c_power4 = c_exp_min + ((double)rand() / RAND_MAX) * (c_exp_max - c_exp_min);
        double c4 =  - ((double)rand() / RAND_MAX) * pow(10, c_power4);   
        double r_power4 = r_exp_min + ((double)rand() / RAND_MAX) * (r_exp_max - r_exp_min);
        double r4 = pow(10, r_power4) * fabs(c4);
        */

        iAp[i] = c1;
        sAp[i] = c1 + r1;
        
        iBp[i] = c3;
        sBp[i] = c3 + r3;

    }

    ref1_pp(iAp, sAp, iBp, sBp, iC1, sC1, n);

    mat_is_cr(iAp, sAp, mAp, rAp, n);
    mat_is_cr(iBp, sBp, mBp, rBp, n);

    mat_is_cr(iC1, sC1, ref_mC1, ref_rC1, n);

    //int_mat_mult(mAt, rAt, mBt, rBt, mCt, rCt, n);

    // calculation
    double start1 = get_time_ms();
    int_mat_mult(mAp, rAp, mBp, rBp, mC1, rC1, n);
    double end1 = get_time_ms();

    // restore
    for (int i = 0; i < N; i ++)
    {
        C1[i] = (C_R){mC1[i], rC1[i]};
    }

    // conversion
    double start2 = get_time_ms();
    for (int i = 0; i < N; i ++)
    {
        mC1_tilde[i] = CR_FP1_adjbis(C1[i]);            
    }
    double end2 = get_time_ms();

    // read r_tilde
    for (int i = 0; i < N; i ++)
    {
        rC1_tilde[i] = FP_CR3(mC1_tilde[i]).radius;
    }

    FILE *result1 = fopen(res1, "a");

    double e1_max = 0.0;
    double e1_min = MAX_RATE;
    double r1_max = 0.0;
    double r1_min = MAX_RATE;

    double b1_max = 0.0;

    //double dil_max = 0.0, dil_mean = 0.0;

    FILE *f = fopen("err1.txt", "a");
    FILE *freq1 = fopen("freq1.txt", "a");

    //FILE *d1 = fopen(dil1, "a");

    for (int i = 0; i < N; i ++)
    {
        double e1 = rC1[i] / mC1[i];
    
        double e1_tilde = rC1_tilde[i] / mC1_tilde[i];

        double e1_rate = e1_tilde / e1;

        double r1_rate = rC1_tilde[i] / rC1[i];

        double b1 = fabs((mC1_tilde[i] - mC1[i]) / mC1[i]);

        //double dil = fabs(rCt[i] / mCt[i]) / FIXED_RAD;

        //fprintf(d1, "%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%d\n", e1_rate, r1_rate, 
        //    e2_rate, r2_rate, e3_rate, r3_rate, e4_rate, r4_rate, i);

        if (r1_rate > 1000.0)
        {
            fprintf(freq1, "n = %d C1[%d]\n", n, i);
            fprint_binary(freq1, mC1[i]);
            fprint_binary(freq1, rC1[i]);
            fprint_binary(freq1, mC1_tilde[i]);
            fprint_binary(freq1, rC1_tilde[i]);
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
        fesetround(FE_TONEAREST);

        e1_max = fmax(e1_max, e1_rate);
        e1_min = fmin(e1_min, e1_rate);

        r1_max = fmax(r1_max, r1_rate);
        r1_min = fmin(r1_min, r1_rate);

        b1_max = fmax(b1, b1_max);

        //dil_max = fmax(dil, dil_max);
        //dil_mean += dil;
    }
    //dil_mean /= N;

    //fclose(d1);
    fclose(freq1);
    

    fprintf(result1, "%.17e\t%.17e\t%.17e\t%.17e\t%.17e\n", e1_max, e1_min, r1_max, r1_min, b1_max);

    fclose(result1); 
    fclose(f);
        
    FILE *fp1 = fopen(file1, "a");
    double t1 = end1 - start1;
    double t2 = end2 - start2;
    double percent1 = t2 / t1 * 100; 
    fprintf(fp1, "%.17e\t%.17e\t%.17e\n", t1, t2, percent1);
    fclose(fp1);

    for (int i = 0; i < N; i ++)
    {
        double c_power1 = c_exp_min + ((double)rand() / RAND_MAX) * (c_exp_max - c_exp_min);
        double c1 = ((double)rand() / RAND_MAX) * pow(10, c_power1);   

        /*
        double c_power2 = c_exp_min + ((double)rand() / RAND_MAX) * (c_exp_max - c_exp_min);
        double c2 =  - ((double)rand() / RAND_MAX) * pow(10, c_power2);   
        */

        double c_power3 = c_exp_min + ((double)rand() / RAND_MAX) * (c_exp_max - c_exp_min);
        double c3 = ((double)rand() / RAND_MAX) * pow(10, c_power3); 

        /*
        double c_power4 = c_exp_min + ((double)rand() / RAND_MAX) * (c_exp_max - c_exp_min);
        double c4 =  - ((double)rand() / RAND_MAX) * pow(10, c_power4);   
        */
       
        mAp[i] = c1;
        rAp[i] = FP_CR3(c1).radius;


        mBp[i] = c3;
        rBp[i] = FP_CR3(c3).radius;

    }

    mat_cr_is(mAp, rAp, iAp, sAp, n);
    mat_cr_is(mBp, rBp, iBp, sBp, n);

    ref1_pp(iAp, sAp, iBp, sBp, iC2, sC2, n);

    mat_is_cr(iC2, sC2, ref_mC2, ref_rC2, n);

    // calculation
    double start3 = get_time_ms();
    int_mat_mult(mAp, rAp, mBp, rBp, mC2, rC2, n);
    double end3 = get_time_ms();

    // restore
    for (int i = 0; i < N; i ++)
    {
        C2[i] = (C_R){mC2[i], rC2[i]};
    }

    // conversion
    double start4 = get_time_ms();
    for (int i = 0; i < N; i ++)
    {
        mC2_tilde[i] = CR_FP1_adjbis(C2[i]);       
    }
    double end4 = get_time_ms();
    
    // read r_tilde
    for (int i = 0; i < N; i ++)
    {
        rC2_tilde[i] = FP_CR3(mC2_tilde[i]).radius;
    }

    FILE *result2 = fopen(res2, "a");

    double e2_max = 0.0;
    double e2_min = MAX_RATE;
    double r2_max = 0.0;
    double r2_min = MAX_RATE;

    double b2_max = 0.0;

    FILE *f2 = fopen("err2.txt", "a");
    FILE *freq2 = fopen("freq2.txt", "a");

    //FILE *d2 = fopen(dil2, "a");

    for (int i = 0; i < N; i ++)
    {
        double e2 = rC2[i] / mC2[i];
    
        double e2_tilde = rC2_tilde[i] / mC2_tilde[i];

        double e2_rate = e2_tilde / e2;

        double r2_rate = rC2_tilde[i] / rC2[i];

        double b2 = fabs((mC2_tilde[i] - mC2[i]) / mC2[i]);

        //fprintf(d2, "%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%d\n", e1_rate, r1_rate, 
        //    e2_rate, r2_rate, e3_rate, r3_rate, e4_rate, r4_rate, i);

        if (r2_rate > 1000.0)
        {
            fprintf(freq2, "n = %d C2[%d]\n", n, i);
            fprint_binary(freq2, mC2[i]);
            fprint_binary(freq2, rC2[i]);
            fprint_binary(freq2, mC2_tilde[i]);
            fprint_binary(freq2, rC2_tilde[i]);
        }

        // check inclusion
        fesetround(FE_UPWARD);
        if (iC2[i] + rC2[i] < mC2[i] || sC2[i] > mC2[i] + rC2[i])
        {
            fprintf(f2, "Error algo n = %d C2[%d]\niC2 = %.10f sC2 = %.10f\nmC2 = %.10f rC2 = %.10f\n", 
                 n, i, iC2[i], sC2[i], mC2[i], rC2[i]);
            fprint_binary(f2, iC2[i]);
            fprint_binary(f2, sC2[i]);
            fprint_binary(f2, mC2[i]);
            fprint_binary(f2, rC2[i]);
        }
        fesetround(FE_TONEAREST);

        e2_max = fmax(e2_max, e2_rate);
        e2_min = fmin(e2_min, e2_rate);

        r2_max = fmax(r2_max, r2_rate);
        r2_min = fmin(r2_min, r2_rate);

        b2_max = fmax(b2, b2_max);
    }

    //fclose(d2);
    fclose(freq2);
    

    fprintf(result2, "%.17e\t%.17e\t%.17e\t%.17e\t%.17e\n", e2_max, e2_min, r2_max, r2_min, b2_max);

    fclose(result2);
    fclose(f2);
        
    FILE *fp2 = fopen(file2, "a");
    double t3 = end3 - start3;
    double t4 = end4 - start4;
    double percent2 = t4 / t3 * 100;
    fprintf(fp2, "%.17e\t%.17e\t%.17e\n", t3, t4, percent2);
    fclose(fp2);

    //FILE *fpt = fopen(test, "a");
    //fprintf(fpt, "%.17e\t%.17e\n", dil_max, dil_mean);
    //fclose(fpt);

    /*
    for (int i = 0; i < N; i += 400)
    {
        printf("\nmCt = ");
        print_binary(mCt[i]);
        printf("\nrCt = ");
        print_binary(rCt[i]);
    }
    */

    free(mAp);
    free(rAp);
    free(mBp);
    free(rBp);

    free(mC1);
    free(rC1);
    free(mC1_tilde);
    free(rC1_tilde);
    free(mC2);
    free(rC2);
    free(mC2_tilde);
    free(rC2_tilde);

    free(C1);
    free(C2);

    free(iAp);
    free(iBp);
    free(sAp);
    free(sBp);
    free(iC1);
    free(iC2);
    free(sC1);
    free(sC2);
    
    /*
    free(mAt);
    free(rAt);
    free(mBt);
    free(rBt);
    free(mCt);
    free(rCt);
    */

    return 0;
}
#include <cblas.h>

#include <fenv.h>
#include <math.h>
#include <stdlib.h>

#include "convert.h"
#include "functions.h"

#define EPS set_pow2(-53)
#define REALMIN set_pow2(-1022)

#define ITERMAX 15
#define TOLERANCE set_pow2(-52)


// an old version which does not work 
void mult_old(double *mA, double *mB, double *rA, double *rB, int k, double *mC, double *rC, double *Id, double *Ones)
{
	fesetround(FE_TONEAREST);
	cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans, k, k, k, 1, mA, k, mB, k, 0, mC, k); //mC
    double *abs_mA = (double *)malloc(k * k * sizeof(double));
    double *abs_mB = (double *)malloc(k * k * sizeof(double));
	for(int i = 0; i < k; i ++) 
    {
		for(int j = 0; j < k; j ++)
        {
			abs_mA[i * k + j] = fabs(mA[i * k + j]);
			abs_mB[i * k + j] = fabs(mB[i * k + j]);
		}
	}
	double *temp_Id = malloc(k * k * sizeof(double));
	double *temp_rB = malloc(k * k * sizeof(double));
	cblas_dcopy(k * k, Id, 1, temp_Id, 1);
	cblas_dcopy(k * k, rB, 1, temp_rB, 1);
	fesetround(FE_UPWARD);	
	cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans, k, k, k, (k + 2) * EPS, abs_mB, k, Id, k, 1, rB, k); //rB	
	cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans, k, k, k, 1, abs_mA, k, rB, k, REALMIN, Ones, k); //Ones
	cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans, k, k, k, 1, abs_mB, k, temp_Id, k, 1, temp_rB, k); //temp_rB
	cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans, k, k, k, 1, rA, k, temp_rB, k, 1, rC, k); //rC
	cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans, k, k, k, 1, Ones, k, Id, k, 1, rC, k); //ones
	free(temp_Id);
	free(temp_rB);
    free(abs_mA);
    free(abs_mB);
	fesetround(FE_TONEAREST);
}


// <mA, rA> * <mB, rB> = <mC, rC>
void int_mat_mult(double *mA, double *rA, double *mB, double *rB, double *mC, double *rC, int n)
{
    int size = n * n;
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans, n, n, n, 1, mA, n, mB, n, 0, mC, n); // mC = rnd_n(mA * mB)
    double *rBb = (double *)malloc(sizeof(double) * size);
    double *abs_mA = (double *)malloc(sizeof(double) * size);
    double *abs_mB = (double *)malloc(sizeof(double) * size);

    cblas_dcopy(size, mA, 1, abs_mA, 1);
    cblas_dcopy(size, mB, 1, abs_mB, 1);

    int signA = get_sign_bit(mA[0]);
    int signB = get_sign_bit(mB[0]);
    if (signA)
    {
        for (int i = 0; i < size; i ++)
        {
            abs_mA[i] = - mA[i];
        }
    }
    // abs_mA = |mA|
    if (signB)
    {
        for (int i = 0; i < size; i ++)
        {
            abs_mB[i] = - mB[i];
        }
    }
    // abs_mB = |mB|
    
    cblas_dcopy(size, rB, 1, rBb, 1);

    fesetround(FE_UPWARD);
    double a = (n + 2) * EPS;
    cblas_daxpy(size, a, abs_mB, 1, rBb, 1);   // rBb = rnd_up(rB + (n + 2) * EPS * |mB|)  
    cblas_daxpy(size, 1, abs_mB, 1, rB, 1);   // rB = rnd_up(rB + |mB|)
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans, n, n, n, 1, abs_mA, n, rBb, n, 0, rC, n); // rC = rnd_n(|mA| * rBb)
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans, n, n, n, 1, rA, n, rB, n, 1, rC, n); // rC = rnd_n(rA * rB + rC)
    for (int i = 0; i < size; i ++)
    {
        rC[i] += REALMIN; 
    }
    // rC = rnd_up(rC + REALMIN)
    fesetround(FE_TONEAREST);

    free(abs_mA);
    free(abs_mB);
    free(rBb);
}

void mat_is_cr(double *iA, double *sA, double *iB, double *sB, double *mA, double *rA, double *mB, double *rB, int n)
{
    int size = n * n;
    fesetround(FE_UPWARD);
    cblas_dcopy(size, iA, 1, mA, 1);
    cblas_dcopy(size, iB, 1, mB, 1);
    cblas_daxpy(size, 1.0, sA, 1, mA, 1);
    cblas_daxpy(size, 1.0, sB, 1, mB, 1);
    for (int i = 0; i < size; i ++)
    {
        mA[i] *= 0.5;
        mB[i] *= 0.5;
    }
    cblas_dcopy(size, mA, 1, rA, 1);
    cblas_dcopy(size, mB, 1, rB, 1);
    cblas_daxpy(size, -1.0, iA, 1, rA, 1);
    cblas_daxpy(size, -1.0, iB, 1, rB, 1);
    fesetround(FE_TONEAREST);
}

void mat_cr_is(double *mA, double *rA, double *mB, double *rB, double *iA, double *sA, double *iB, double *sB, int n)
{
    int size = n * n;
    cblas_dcopy(size, mA, 1, iA, 1);
    cblas_dcopy(size, mB, 1, iB, 1);
    cblas_dcopy(size, mA, 1, sA, 1);
    cblas_dcopy(size, mB, 1, sB, 1);
    fesetround(FE_UPWARD);
    cblas_daxpy(size, 1.0, rA, 1, sA, 1); // sA = rnd_up(mA + rA)
    cblas_daxpy(size, 1.0, rB, 1, sB, 1); // sB = rnd_up(mB + rB)
    fesetround(FE_DOWNWARD);
    cblas_daxpy(size, -1.0, rA, 1, iA, 1); // iA = rnd_down(mA - rA)
    cblas_daxpy(size, -1.0, rB, 1, iB, 1); // iB = rnd_down(mB - rB)
    fesetround(FE_TONEAREST);
}

// 0 < iA < sA, 0 < iB < sB
void ref1_pp(double *iA, double *sA, double *iB, double *sB, double *iC, double *sC, int n)
{
    fesetround(FE_DOWNWARD);
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans, n, n, n, 1, iA, n, iB, n, 0, iC, n); // iC = rnd_down(iA * iB)
    fesetround(FE_UPWARD);
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans, n, n, n, 1, sA, n, sB, n, 0, sC, n); // sC = rnd_up(sA * sB)
    fesetround(FE_TONEAREST);
}

// iA < sA < 0, iB < sB < 0
void ref2_mm(double *iA, double *sA, double *iB, double *sB, double *iC, double *sC, int n)
{
    fesetround(FE_DOWNWARD);
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans, n, n, n, 1, sA, n, sB, n, 0, iC, n); // iC = rnd_down(sA * sB)
    fesetround(FE_UPWARD);
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans, n, n, n, 1, iA, n, iB, n, 0, sC, n); // sC = rnd_up(iA * iB)
    fesetround(FE_TONEAREST);
}


// 0 < iA < sA, iB < sB < 0
void ref3_pm(double *iA, double *sA, double *iB, double *sB, double *iC, double *sC, int n)
{
    fesetround(FE_DOWNWARD);
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans, n, n, n, 1, sA, n, iB, n, 0, iC, n); // iC = rnd_down(sA * iB)
    fesetround(FE_UPWARD);
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans, n, n, n, 1, iA, n, sB, n, 0, sC, n); // sC = rnd_up(iA * sB)
    fesetround(FE_TONEAREST);
}

// iA < sA < 0, 0 < iB < sB
void ref4_mp(double *iA, double *sA, double *iB, double *sB, double *iC, double *sC, int n)
{
    fesetround(FE_DOWNWARD);
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans, n, n, n, 1, iA, n, sB, n, 0, iC, n); // iC = rnd_down(iA * sB)
    fesetround(FE_UPWARD);
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans, n, n, n, 1, sA, n, iB, n, 0, sC, n); // sC = rnd_up(sA * iB)
    fesetround(FE_TONEAREST);
}


C_R inverse(C_R cr_x)
{
    fesetround(FE_DOWNWARD);
	double c1 = ((-1) / (-fabs(cr_x.center) - cr_x.radius));
	fesetround(FE_UPWARD);
	double c2 = ((-1)/(-fabs(cr_x.center) + cr_x.radius));
	double c = (c1 + 0.5 * (c2 - c1));
	cr_x.radius = (c - c1);
	fesetround(FE_TONEAREST);
    int s = get_sign_bit(cr_x.center);
    cr_x.center = set_sign_bit(c, s);
    return cr_x;
}


C_R intersection(C_R x, C_R y)
{
    fesetround(FE_DOWNWARD);
    double inf = fmax(x.center - x.radius, y.center - y.radius);
    fesetround(FE_UPWARD);
    double sup = fmin(x.center + x.radius, y.center + y.radius);  
    if (inf > sup)
    {
        printf("Error: intersection is empty!\n");
        return (C_R){0, 0};
    }
    double center = (inf + sup) / 2;
    double radius = center - inf;
    fesetround(FE_TONEAREST);   
    C_R cr = {center, radius};
    return cr;
}



C_R newton(C_R cr_x)
{   
    printf("C_R Newton:\n");
    for (int i = 0; i < ITERMAX; i ++)
    {
        double c1 = cr_x.center;
        double r1 = cr_x.radius;
        C_R cr_inv = inverse((C_R){c1, r1});
        double c2 = c1 - (c1 * c1 - 2) / 2 * cr_inv.center;
        double r2 = fabs((c1 * c1 - 2) / 2 * cr_inv.radius);
        cr_x = intersection((C_R){c1, r1}, (C_R){c2, r2});
        printf("\niteration %d:\nc = ",i);
        print_binary(cr_x.center);
        printf("\nr = ");
        print_binary(cr_x.radius);
        printf("\n\n");
        if (fabs(cr_x.center - c1) <= EPS && cr_x.radius <= TOLERANCE)
        {
            break;
        }
    }
    return cr_x;
}
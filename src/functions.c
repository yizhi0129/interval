#include <cblas.h>

#include <fenv.h>
#include <math.h>
#include <stdlib.h>

#include "convert.h"
#include "functions.h"


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


    for (int i = 0; i < size; i ++)
    {
        abs_mA[i] = fabs(mA[i]);
        abs_mB[i] = fabs(mB[i]);
    }
    // abs_mA = |mA|, abs_mB = |mB|
    
    cblas_dcopy(size, rB, 1, rBb, 1);

    fesetround(FE_UPWARD);
    double a = (n + 2) * EPS;
    cblas_daxpy(size, a, abs_mB, 1, rBb, 1);   // rBb = rnd_up(rB + (n + 2) * EPS * |mB|)  
    cblas_daxpy(size, 1, abs_mB, 1, rB, 1);   // rB = rnd_up(rB + |mB|)
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans, n, n, n, 1, abs_mA, n, rBb, n, 0, rC, n); // rC = rnd_up(|mA| * rBb)
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans, n, n, n, 1, rA, n, rB, n, 1, rC, n); // rC = rnd_up(rA * rB + rC)
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

void mat_is_cr(double *iX, double *sX, double *mX, double *rX, int n)
{
    int size = n * n;
    fesetround(FE_UPWARD);
    cblas_dcopy(size, iX, 1, mX, 1);
    cblas_daxpy(size, 1.0, sX, 1, mX, 1);
    for (int i = 0; i < size; i ++)
    {
        mX[i] *= 0.5;
    }
    cblas_dcopy(size, mX, 1, rX, 1);
    cblas_daxpy(size, -1.0, iX, 1, rX, 1);
    fesetround(FE_TONEAREST);
}

void mat_cr_is(double *mX, double *rX, double *iX, double *sX, int n)
{
    int size = n * n;
    cblas_dcopy(size, mX, 1, iX, 1);
    cblas_dcopy(size, mX, 1, sX, 1);
    fesetround(FE_UPWARD);
    cblas_daxpy(size, 1.0, rX, 1, sX, 1); // sX = rnd_up(mX + rX)
    fesetround(FE_DOWNWARD);
    cblas_daxpy(size, -1.0, rX, 1, iX, 1); // iX = rnd_down(mX - rX)
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
	double c = (c1 + (c2 - c1) / 2);
	cr_x.radius = (c - c1);
	fesetround(FE_TONEAREST);
    int s = get_sign_bit(cr_x.center);
    cr_x.center = set_sign_bit(c, s);
    return cr_x;
}


C_R intersection(C_R x, C_R y)
{
    fesetround(FE_DOWNWARD);
    double inf_x = x.center - x.radius;
    double inf_y = y.center - y.radius;
    double inf = fmax(inf_x, inf_y);
    fesetround(FE_UPWARD);
    double sup_x = x.center + x.radius;
    double sup_y = y.center + y.radius;
    double sup = fmin(sup_x, sup_y);
    if (inf > sup)
    {
        printf("inf_x = %lf, sup_x = %lf, inf_y = %lf, sup_y = %lf\n", inf_x, sup_x, inf_y, sup_y);
        printf("Error: intersection is empty! back to previous interval\n");
        return y;
    }
    double center = (inf + sup) / 2;
    double radius = center - inf;
    fesetround(FE_TONEAREST);   
    C_R cr = {center, radius};
    return cr;
}



C_R newton_pr(C_R cr_x)
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
        printf("iteration %d:\nc = ", i + 1);
        print_binary(cr_x.center);
        printf("\nr = ");
        print_binary(cr_x.radius);
        printf("\n");
        if (fabs(cr_x.center - c1) <= EPS && cr_x.radius <= TOLERANCE)
        {
            break;
        }
    }
    return cr_x;
}

C_R newton_res(C_R cr_x)
{   
    for (int i = 0; i < ITERMAX; i ++)
    {
        double c1 = cr_x.center;
        double r1 = cr_x.radius;
        C_R cr_inv = inverse((C_R){c1, r1});
        double c2 = c1 - (c1 * c1 - 2) / 2 * cr_inv.center;
        double r2 = fabs((c1 * c1 - 2) / 2 * cr_inv.radius);
        cr_x = intersection((C_R){c1, r1}, (C_R){c2, r2});
        if (fabs(cr_x.center - c1) <= EPS && cr_x.radius <= TOLERANCE)
        {
            break;
        }
    }
    return cr_x;
}

C_R interval_add(C_R x, C_R y)
{
    double center = x.center + y.center;
    fesetround(FE_UPWARD);
    double radius = x.radius + y.radius + fabs(center) * EPS;
    fesetround(FE_TONEAREST);
    return (C_R){center, radius};
}

C_R interval_sub(C_R x, C_R y)
{
    double center = x.center - y.center;
    fesetround(FE_UPWARD);
    double radius = x.radius + y.radius + fabs(center) * EPS;
    fesetround(FE_TONEAREST);
    return (C_R){center, radius};
}

C_R interval_mult(C_R x, C_R y)
{
    double center = x.center * y.center;
    fesetround(FE_UPWARD);
    double radius = ETA + EPS * fabs(center) + (fabs(x.center) + x.radius) * y.radius + x.radius * fabs(y.center);
    fesetround(FE_TONEAREST);
    return (C_R){center, radius};
}


// A in band storage format diagonal n, subdiagonal n-1, superdiagonal n-1
// b is a vector of size n
void interval_GS_tridiag(C_R *A, C_R *b, C_R *x, int n)
{
    C_R *x_prev = malloc(n * sizeof(C_R));
    for (int i = 0; i < n; i ++)
    {
        x_prev[i] = x[i]; 
    }

    int count = 0;

    for (int i = 0; i < 1000; i ++) 
    {  
        for (int i = 0; i < n; i ++) 
        {
            C_R sum = b[i];

            if (i > 0) 
            {
                C_R prod = interval_mult(A[n + i - 1], x[i - 1]);
                sum = interval_sub(sum, prod);
            }

            if (i < n - 1) 
            {
                C_R prod = interval_mult(A[2 * n - 1 + i], x_prev[i + 1]);
                sum = interval_sub(sum, prod);
            }

            //printf("A[%d]: %lf, %lf\n", i, A[i].center, A[i].radius);
            C_R Aii_inv = inverse(A[i]);
            //printf("Inverse of A[%d]: %lf, %lf\n", i, Aii_inv.center, Aii_inv.radius);
            x[i] = interval_mult(sum, Aii_inv);  // x_i^{(k+1)} = sum / A_ii

            x[i] = intersection(x[i], x_prev[i]);
        }

        count ++;

        // Check convergence
        double max_diff1 = 0.0, max_diff2 = 0.0;
        for (int i = 0; i < n; i ++)
        {
            double diff1 = fabs((x[i].center - x_prev[i].center) / x[i].center);
            double diff2 = fabs((x[i].radius - x_prev[i].radius) / x[i].radius);
            if (diff1 > max_diff1)
            {
                max_diff1 = diff1;
            }
            if (diff2 > max_diff2)
            {
                max_diff2 = diff2;
            }
        }
        //printf("%d %lf %lf\n", count, max_diff1, max_diff2);

        if (max_diff1 < TOLERANCE && max_diff2 < TOLERANCE)
        {
            break;
        }

        for (int i = 0; i < n; i ++)
        {
            x_prev[i] = x[i];  // Update x_prev for the next iteration
        }
    }

    free(x_prev);
}

// A in CSR format
// idx is the index of the first non-zero element in each row
// col_id is the column index of each non-zero element
// b is a vector of size n
void interval_GS_CSR(C_R *A, int *idx, int *col_id, C_R *b, C_R *x, int n)
{
    C_R *x_prev = malloc(n * sizeof(C_R));
    for (int i = 0; i < n; i ++)
    {
        x_prev[i] = x[i]; 
    }

    for (int iter = 0; iter < 1000; iter ++) 
    {  
        //printf("Iteration %d:\n", iter + 1);
        for (int i = 0; i < n; i ++) 
        {
            C_R sum = b[i];
            C_R Aii_inv = { 0.0, 0.0 };

            for (int k = idx[i]; k < idx[i + 1]; k ++)
            {
                int j = col_id[k];

                if (j < i) 
                {
                    C_R prod = interval_mult(A[k], x[j]);
                    sum = interval_sub(sum, prod);
                }
                else if (j == i)
                {
                    Aii_inv = inverse(A[k]);
                }
                else if (j > i) 
                {
                    C_R prod = interval_mult(A[k], x_prev[j]);
                    sum = interval_sub(sum, prod);
                }              
            }
            x[i] = interval_mult(sum, Aii_inv);  // x_i^{(k+1)} = sum / A_ii

            x[i] = intersection(x[i], x_prev[i]);
        }

        // Check convergence
        double max_diff1 = 0.0, max_diff2 = 0.0;
        for (int i = 0; i < n; i ++)
        {
            double diff1 = fabs((x[i].center - x_prev[i].center) / x[i].center);
            double diff2 = fabs((x[i].radius - x_prev[i].radius) / x[i].radius);
            if (diff1 > max_diff1)
            {
                max_diff1 = diff1;
            }
            if (diff2 > max_diff2)
            {
                max_diff2 = diff2;
            }
        }

        if (max_diff1 < TOLERANCE && max_diff2 < TOLERANCE)
        {
            break;
        }

        for (int i = 0; i < n; i ++)
        {
            x_prev[i] = x[i];  // Update x_prev for the next iteration
        }
    }

    free(x_prev);
}


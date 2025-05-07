#ifndef _FUNCTIONS_H_
#define _FUNCTIONS_H_

#include "convert.h"

// old version of matrix multiplication: does not work
void mult_old(double *mA, double *mB, double *rA, double *rB, int k, double *mC, double *rC, double *Id, double *Ones);

void int_mat_mult(double *mA, double *rA, double *mB, double *rB, double *mC, double *rC, int n);

void mat_is_cr(double *iA, double *sA, double *iB, double *sB, double *mA, double *rA, double *mB, double *rB, int n);
void mat_cr_is(double *mA, double *rA, double *mB, double *rB, double *iA, double *sA, double *iB, double *sB, int n);

void ref1_pp(double *iA, double *sA, double *iB, double *sB, double *iC, double *sC, int n);
void ref2_mm(double *iA, double *sA, double *iB, double *sB, double *iC, double *sC, int n);
void ref3_pm(double *iA, double *sA, double *iB, double *sB, double *iC, double *sC, int n);
void ref4_mp(double *iA, double *sA, double *iB, double *sB, double *iC, double *sC, int n);

C_R inverse(C_R cr_x);

C_R intersection(C_R cr_x, C_R cr_y);

C_R newton(C_R cr_x);

#endif
#include "softposit.h"
#include "posit_interval.h"

#include "convert.h"

#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include <fenv.h>


int main()
{
    double a = 1.23456;
    double b = 5.9876; 
    double c = -3.1415926;
    double d = -2.7182818;

    posit8_t p8_a = convertDoubleToP8(a);
    posit8_t p8_b = convertDoubleToP8(b);
    posit8_t p8_c = convertDoubleToP8(c);
    posit8_t p8_d = convertDoubleToP8(d);

    posit16_t p16_a = convertDoubleToP16(a);
    posit16_t p16_b = convertDoubleToP16(b);
    posit16_t p16_c = convertDoubleToP16(c);
    posit16_t p16_d = convertDoubleToP16(d);

    posit32_t p32_a = convertDoubleToP32(a);
    posit32_t p32_b = convertDoubleToP32(b);
    posit32_t p32_c = convertDoubleToP32(c);
    posit32_t p32_d = convertDoubleToP32(d);

    printf("prec\ta\tb\tc\td\n");
    printf("8\t%.17e\t%.17e\t%.17e\t%.17e\n", convertP8ToDouble(p8_a), convertP8ToDouble(p8_b), convertP8ToDouble(p8_c), convertP8ToDouble(p8_d));
    printf("16\t%.17e\t%.17e\t%.17e\t%.17e\n", convertP16ToDouble(p16_a), convertP16ToDouble(p16_b), convertP16ToDouble(p16_c), convertP16ToDouble(p16_d));
    printf("32\t%.17e\t%.17e\t%.17e\t%.17e\n", convertP32ToDouble(p32_a), convertP32ToDouble(p32_b), convertP32ToDouble(p32_c), convertP32ToDouble(p32_d));
    printf("\n");

    posit16_t p8_a_16 = p8_to_p16(p8_a);
    posit16_t p8_b_16 = p8_to_p16(p8_b);
    posit16_t p8_c_16 = p8_to_p16(p8_c);
    posit16_t p8_d_16 = p8_to_p16(p8_d);

    posit32_t p16_a_32 = p16_to_p32(p16_a);
    posit32_t p16_b_32 = p16_to_p32(p16_b);
    posit32_t p16_c_32 = p16_to_p32(p16_c);
    posit32_t p16_d_32 = p16_to_p32(p16_d);

    // a + b
    posit16_t s1 = p16_add(p16_a, p16_b);
    posit16_t s11 = p16_add(p8_a_16, p8_b_16);
    posit32_t s2 = p32_add(p32_a, p32_b);
    posit32_t s22 = p32_add(p16_a_32, p16_b_32);

    posit16_t s1_up = p16_add_up(p16_a, p16_b);
    posit16_t s11_up = p16_add_up(p8_a_16, p8_b_16);
    posit32_t s2_up = p32_add_up(p32_a, p32_b);
    posit32_t s22_up = p32_add_up(p16_a_32, p16_b_32);

    posit16_t s1_down = p16_add_down(p16_a, p16_b);
    posit16_t s11_down = p16_add_down(p8_a_16, p8_b_16);
    posit32_t s2_down = p32_add_down(p32_a, p32_b);
    posit32_t s22_down = p32_add_down(p16_a_32, p16_b_32);

    double sum = a + b;
    fesetround(FE_UPWARD);
    double sum_up = a + b;
    fesetround(FE_DOWNWARD);
    double sum_down = a + b;
    fesetround(FE_TONEAREST);
    printf("a+b\t%.17e\t%.17e\t%.17e\n", sum, sum_up, sum_down);
    printf("p16\t%.17e\t%.17e\t%.17e\n", convertP16ToDouble(s1), convertP16ToDouble(s1_up), convertP16ToDouble(s1_down));
    printf("p8->p16\t%.17e\t%.17e\t%.17e\n", convertP16ToDouble(s11), convertP16ToDouble(s11_up), convertP16ToDouble(s11_down));    
    printf("p32\t%.17e\t%.17e\t%.17e\n", convertP32ToDouble(s2), convertP32ToDouble(s2_up), convertP32ToDouble(s2_down));
    printf("p16->p32\t%.17e\t%.17e\t%.17e\n", convertP32ToDouble(s22), convertP32ToDouble(s22_up), convertP32ToDouble(s22_down));
    printf("\n");

    // c + d
    posit16_t s3 = p16_add(p16_c, p16_d);
    posit16_t s33 = p16_add(p8_c_16, p8_d_16);
    posit32_t s4 = p32_add(p32_c, p32_d);
    posit32_t s44 = p32_add(p16_c_32, p16_d_32);

    posit16_t s3_up = p16_add_up(p16_c, p16_d);
    posit16_t s33_up = p16_add_up(p8_c_16, p8_d_16);
    posit32_t s4_up = p32_add_up(p32_c, p32_d);
    posit32_t s44_up = p32_add_up(p16_c_32, p16_d_32);

    posit16_t s3_down = p16_add_down(p16_c, p16_d);
    posit16_t s33_down = p16_add_down(p8_c_16, p8_d_16);
    posit32_t s4_down = p32_add_down(p32_c, p32_d);
    posit32_t s44_down = p32_add_down(p16_c_32, p16_d_32);

    double sum_c_d = c + d;
    fesetround(FE_UPWARD);
    double sum_c_d_up = c + d;
    fesetround(FE_DOWNWARD);
    double sum_c_d_down = c + d;
    fesetround(FE_TONEAREST);

    printf("c+d\t%.17e\t%.17e\t%.17e\n", sum_c_d, sum_c_d_up, sum_c_d_down);
    printf("p16\t%.17e\t%.17e\t%.17e\n", convertP16ToDouble(s3), convertP16ToDouble(s3_up), convertP16ToDouble(s3_down));
    printf("p8->p16\t%.17e\t%.17e\t%.17e\n", convertP16ToDouble(s33), convertP16ToDouble(s33_up), convertP16ToDouble(s33_down));
    printf("p32\t%.17e\t%.17e\t%.17e\n", convertP32ToDouble(s4), convertP32ToDouble(s4_up), convertP32ToDouble(s4_down));
    printf("p16->p32\t%.17e\t%.17e\t%.17e\n", convertP32ToDouble(s44), convertP32ToDouble(s44_up), convertP32ToDouble(s44_down));
    printf("\n");

    // a + c
    posit16_t s5 = p16_add(p16_a, p16_c);
    posit16_t s55 = p16_add(p8_a_16, p8_c_16);
    posit32_t s6 = p32_add(p32_a, p32_c);
    posit32_t s66 = p32_add(p16_a_32, p16_c_32);

    posit16_t s5_up = p16_add_up(p16_a, p16_c);
    posit16_t s55_up = p16_add_up(p8_a_16, p8_c_16);
    posit32_t s6_up = p32_add_up(p32_a, p32_c);
    posit32_t s66_up = p32_add_up(p16_a_32, p16_c_32);

    posit16_t s5_down = p16_add_down(p16_a, p16_c);
    posit16_t s55_down = p16_add_down(p8_a_16, p8_c_16);
    posit32_t s6_down = p32_add_down(p32_a, p32_c);
    posit32_t s66_down = p32_add_down(p16_a_32, p16_c_32);

    double sum_a_c = a + c;
    fesetround(FE_UPWARD);
    double sum_a_c_up = a + c;
    fesetround(FE_DOWNWARD);
    double sum_a_c_down = a + c;
    fesetround(FE_TONEAREST);

    printf("a+c\t%.17e\t%.17e\t%.17e\n", sum_a_c, sum_a_c_up, sum_a_c_down);
    printf("p16\t%.17e\t%.17e\t%.17e\n", convertP16ToDouble(s5), convertP16ToDouble(s5_up), convertP16ToDouble(s5_down));
    printf("p8->p16\t%.17e\t%.17e\t%.17e\n", convertP16ToDouble(s55), convertP16ToDouble(s55_up), convertP16ToDouble(s55_down));
    printf("p32\t%.17e\t%.17e\t%.17e\n", convertP32ToDouble(s6), convertP32ToDouble(s6_up), convertP32ToDouble(s6_down));
    printf("p16->p32\t%.17e\t%.17e\t%.17e\n", convertP32ToDouble(s66), convertP32ToDouble(s66_up), convertP32ToDouble(s66_down));
    printf("\n");

    // d + b
    posit16_t s7 = p16_add(p16_d, p16_b);
    posit16_t s77 = p16_add(p8_d_16, p8_b_16);
    posit32_t s8 = p32_add(p32_d, p32_b);
    posit32_t s88 = p32_add(p16_d_32, p16_b_32);

    posit16_t s7_up = p16_add_up(p16_d, p16_b);
    posit16_t s77_up = p16_add_up(p8_d_16, p8_b_16);
    posit32_t s8_up = p32_add_up(p32_d, p32_b);
    posit32_t s88_up = p32_add_up(p16_d_32, p16_b_32);

    posit16_t s7_down = p16_add_down(p16_d, p16_b);
    posit16_t s77_down = p16_add_down(p8_d_16, p8_b_16);
    posit32_t s8_down = p32_add_down(p32_d, p32_b);
    posit32_t s88_down = p32_add_down(p16_d_32, p16_b_32);

    double sum_d_b = d + b;
    fesetround(FE_UPWARD);
    double sum_d_b_up = d + b;
    fesetround(FE_DOWNWARD);
    double sum_d_b_down = d + b;
    fesetround(FE_TONEAREST);

    printf("d+b\t%.17e\t%.17e\t%.17e\n", sum_d_b, sum_d_b_up, sum_d_b_down);
    printf("p16\t%.17e\t%.17e\t%.17e\n", convertP16ToDouble(s7), convertP16ToDouble(s7_up), convertP16ToDouble(s7_down));
    printf("p8->p16\t%.17e\t%.17e\t%.17e\n", convertP16ToDouble(s77), convertP16ToDouble(s77_up), convertP16ToDouble(s77_down));
    printf("p32\t%.17e\t%.17e\t%.17e\n", convertP32ToDouble(s8), convertP32ToDouble(s8_up), convertP32ToDouble(s8_down));
    printf("p16->p32\t%.17e\t%.17e\t%.17e\n", convertP32ToDouble(s88), convertP32ToDouble(s88_up), convertP32ToDouble(s88_down));
    printf("\n");

    // b - a
    posit16_t s9 = p16_sub(p16_b, p16_a);
    posit16_t s99 = p16_sub(p8_b_16, p8_a_16);
    posit32_t s10 = p32_sub(p32_b, p32_a);
    posit32_t s110 = p32_sub(p16_b_32, p16_a_32);

    posit16_t s9_up = p16_sub_up(p16_b, p16_a);
    posit16_t s99_up = p16_sub_up(p8_b_16, p8_a_16);
    posit32_t s10_up = p32_sub_up(p32_b, p32_a);
    posit32_t s110_up = p32_sub_up(p16_b_32, p16_a_32);

    posit16_t s9_down = p16_sub_down(p16_b, p16_a);
    posit16_t s99_down = p16_sub_down(p8_b_16, p8_a_16);
    posit32_t s10_down = p32_sub_down(p32_b, p32_a);
    posit32_t s110_down = p32_sub_down(p16_b_32, p16_a_32);

    double sum_b_a = b - a;
    fesetround(FE_UPWARD);
    double sum_b_a_up = b - a;
    fesetround(FE_DOWNWARD);
    double sum_b_a_down = b - a;
    fesetround(FE_TONEAREST);

    printf("b-a\t%.17e\t%.17e\t%.17e\n", sum_b_a, sum_b_a_up, sum_b_a_down);
    printf("p16\t%.17e\t%.17e\t%.17e\n", convertP16ToDouble(s9), convertP16ToDouble(s9_up), convertP16ToDouble(s9_down));
    printf("p8->p16\t%.17e\t%.17e\t%.17e\n", convertP16ToDouble(s99), convertP16ToDouble(s99_up), convertP16ToDouble(s99_down));
    printf("p32\t%.17e\t%.17e\t%.17e\n", convertP32ToDouble(s10), convertP32ToDouble(s10_up), convertP32ToDouble(s10_down));
    printf("p16->p32\t%.17e\t%.17e\t%.17e\n", convertP32ToDouble(s110), convertP32ToDouble(s110_up), convertP32ToDouble(s110_down));
    printf("\n");

    // c - d
    s11 = p16_sub(p16_c, p16_d);
    posit16_t s111 = p16_sub(p8_c_16, p8_d_16);
    posit32_t s12 = p32_sub(p32_c, p32_d);
    posit32_t s122 = p32_sub(p16_c_32, p16_d_32);

    s11_up = p16_sub_up(p16_c, p16_d);
    posit16_t s111_up = p16_sub_up(p8_c_16, p8_d_16);
    posit32_t s12_up = p32_sub_up(p32_c, p32_d);
    posit32_t s122_up = p32_sub_up(p16_c_32, p16_d_32);

    s11_down = p16_sub_down(p16_c, p16_d);
    posit16_t s111_down = p16_sub_down(p8_c_16, p8_d_16);
    posit32_t s12_down = p32_sub_down(p32_c, p32_d);
    posit32_t s122_down = p32_sub_down(p16_c_32, p16_d_32);

    sum_c_d = c - d;
    fesetround(FE_UPWARD);
    sum_c_d_up = c - d;
    fesetround(FE_DOWNWARD);
    sum_c_d_down = c - d;
    fesetround(FE_TONEAREST);

    printf("c-d\t%.17e\t%.17e\t%.17e\n", sum_c_d, sum_c_d_up, sum_c_d_down);
    printf("p16\t%.17e\t%.17e\t%.17e\n", convertP16ToDouble(s11), convertP16ToDouble(s11_up), convertP16ToDouble(s11_down));
    printf("p8->p16\t%.17e\t%.17e\t%.17e\n", convertP16ToDouble(s111), convertP16ToDouble(s111_up), convertP16ToDouble(s111_down));
    printf("p32\t%.17e\t%.17e\t%.17e\n", convertP32ToDouble(s12), convertP32ToDouble(s12_up), convertP32ToDouble(s12_down));
    printf("p16->p32\t%.17e\t%.17e\t%.17e\n", convertP32ToDouble(s122), convertP32ToDouble(s122_up), convertP32ToDouble(s122_down));
    printf("\n");

    // b - c
    posit16_t s13 = p16_sub(p16_b, p16_c);
    posit16_t s131 = p16_sub(p8_b_16, p8_c_16);
    posit32_t s14 = p32_sub(p32_b, p32_c);
    posit32_t s142 = p32_sub(p16_b_32, p16_c_32);

    posit16_t s13_up = p16_sub_up(p16_b, p16_c);
    posit16_t s131_up = p16_sub_up(p8_b_16, p8_c_16);
    posit32_t s14_up = p32_sub_up(p32_b, p32_c);
    posit32_t s142_up = p32_sub_up(p16_b_32, p16_c_32);

    posit16_t s13_down = p16_sub_down(p16_b, p16_c);
    posit16_t s131_down = p16_sub_down(p8_b_16, p8_c_16);
    posit32_t s14_down = p32_sub_down(p32_b, p32_c);
    posit32_t s142_down = p32_sub_down(p16_b_32, p16_c_32);

    double sum_b_c = b - c;
    fesetround(FE_UPWARD);
    double sum_b_c_up = b - c;
    fesetround(FE_DOWNWARD);
    double sum_b_c_down = b - c;
    fesetround(FE_TONEAREST);

    printf("b-c\t%.17e\t%.17e\t%.17e\n", sum_b_c, sum_b_c_up, sum_b_c_down);
    printf("p16\t%.17e\t%.17e\t%.17e\n", convertP16ToDouble(s13), convertP16ToDouble(s13_up), convertP16ToDouble(s13_down));
    printf("p8->p16\t%.17e\t%.17e\t%.17e\n", convertP16ToDouble(s131), convertP16ToDouble(s131_up), convertP16ToDouble(s131_down));
    printf("p32\t%.17e\t%.17e\t%.17e\n", convertP32ToDouble(s14), convertP32ToDouble(s14_up), convertP32ToDouble(s14_down));
    printf("p16->p32\t%.17e\t%.17e\t%.17e\n", convertP32ToDouble(s142), convertP32ToDouble(s142_up), convertP32ToDouble(s142_down));
    printf("\n");

    // d - a
    posit16_t s15 = p16_sub(p16_d, p16_a);
    posit16_t s151 = p16_sub(p8_d_16, p8_a_16);
    posit32_t s16 = p32_sub(p32_d, p32_a);
    posit32_t s162 = p32_sub(p16_d_32, p16_a_32);

    posit16_t s15_up = p16_sub_up(p16_d, p16_a);
    posit16_t s151_up = p16_sub_up(p8_d_16, p8_a_16);
    posit32_t s16_up = p32_sub_up(p32_d, p32_a);
    posit32_t s162_up = p32_sub_up(p16_d_32, p16_a_32);

    posit16_t s15_down = p16_sub_down(p16_d, p16_a);
    posit16_t s151_down = p16_sub_down(p8_d_16, p8_a_16);
    posit32_t s16_down = p32_sub_down(p32_d, p32_a);
    posit32_t s162_down = p32_sub_down(p16_d_32, p16_a_32);

    double sum_d_a = d - a;
    fesetround(FE_UPWARD);
    double sum_d_a_up = d - a;
    fesetround(FE_DOWNWARD);
    double sum_d_a_down = d - a;
    fesetround(FE_TONEAREST);

    printf("d-a\t%.17e\t%.17e\t%.17e\n", sum_d_a, sum_d_a_up, sum_d_a_down);
    printf("p16\t%.17e\t%.17e\t%.17e\n", convertP16ToDouble(s15), convertP16ToDouble(s15_up), convertP16ToDouble(s15_down));
    printf("p8->p16\t%.17e\t%.17e\t%.17e\n", convertP16ToDouble(s151), convertP16ToDouble(s151_up), convertP16ToDouble(s151_down));
    printf("p32\t%.17e\t%.17e\t%.17e\n", convertP32ToDouble(s16), convertP32ToDouble(s16_up), convertP32ToDouble(s16_down));
    printf("p16->p32\t%.17e\t%.17e\t%.17e\n", convertP32ToDouble(s162), convertP32ToDouble(s162_up), convertP32ToDouble(s162_down));
    printf("\n");
    
    return 0;
}
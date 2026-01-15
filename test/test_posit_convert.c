#include "softposit.h"
#include "posit_interval.h"

#include "convert.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define N 20000

int main(int argc, char **argv)
{
    int total = 2 * N;
    C_R *test_int8 = malloc(total * sizeof(C_R));
    C_R *test_int16 = malloc(total * sizeof(C_R));

    const double c_exp_min8 = 0.4;
    const double c_exp_max8 = 1.8;
    const double r_exp_min8 = -2.0;
    const double r_exp_max8 = -1.0;

    const double c_exp_min16 = -2.0;
    const double c_exp_max16 = 7.0;
    const double r_exp_min16 = -5.0;
    const double r_exp_max16 = -1.0;

    srand(get_time_ms());    

    for (int i = 0; i < N; i ++) 
    {
        double c_power1 = c_exp_min8 + ((double)rand() / RAND_MAX) * (c_exp_max8 - c_exp_min8);
        double c1 = ((double)rand() / RAND_MAX) * pow(10, c_power1);   
        double r_power1 = r_exp_min8 + ((double)rand() / RAND_MAX) * (r_exp_max8 - r_exp_min8);
        double r1 = pow(10, r_power1) * fabs(c1);    
        test_int8[i].center = c1;
        test_int8[i].radius = r1;

        double c_power2 = c_exp_min8 + ((double)rand() / RAND_MAX) * (c_exp_max8 - c_exp_min8);
        double c2 = - ((double)rand() / RAND_MAX) * pow(10, c_power2);   
        double r_power2 = r_exp_min8 + ((double)rand() / RAND_MAX) * (r_exp_max8 - r_exp_min8);
        double r2 = pow(10, r_power2) * fabs(c2);     
        test_int8[i + N].center = c2;
        test_int8[i + N].radius = r2;

        double c_power11 = c_exp_min16 + ((double)rand() / RAND_MAX) * (c_exp_max16 - c_exp_min16);
        double c11 = ((double)rand() / RAND_MAX) * pow(10, c_power11);   
        double r_power11 = r_exp_min16 + ((double)rand() / RAND_MAX) * (r_exp_max16 - r_exp_min16);
        double r11 = pow(10, r_power11) * fabs(c11);    
        test_int16[i].center = c11;
        test_int16[i].radius = r11;

        double c_power22 = c_exp_min16 + ((double)rand() / RAND_MAX) * (c_exp_max16 - c_exp_min16);
        double c22 = - ((double)rand() / RAND_MAX) * pow(10, c_power22);   
        double r_power22 = r_exp_min16 + ((double)rand() / RAND_MAX) * (r_exp_max16 - r_exp_min8);
        double r22 = pow(10, r_power22) * fabs(c22);     
        test_int16[i + N].center = c22;
        test_int16[i + N].radius = r22;
    }

    posit8_t * p8_c = malloc(total * sizeof(posit8_t));
    posit8_t * p8_r = malloc(total * sizeof(posit8_t));
    posit16_t * p16_c = malloc(total * sizeof(posit16_t));
    posit16_t * p16_r = malloc(total * sizeof(posit16_t));

    posit32_t * ipb_8 = malloc(total * sizeof(posit32_t));
    posit32_t * ipb_16 = malloc(total * sizeof(posit32_t)); 

    p16_FP_INT_P * p16_interval = malloc(total * sizeof(p16_FP_INT_P));
    p32_FP_INT_P * p32_interval = malloc(total * sizeof(p32_FP_INT_P));

    posit32_t * ipb_8_convert = malloc(total * sizeof(posit32_t));
    posit32_t * ipb_16_convert = malloc(total * sizeof(posit32_t));

    double ipb_8_avg = 0.0, ipb_8_c_avg = 0.0, ipb_16_avg = 0.0, ipb_16_c_avg = 0.0, avg_rate8 = 0.0, avg_rate16 = 0.0;
    double ipb_8_max = 0.0, ipb_8_c_max = 0.0, ipb_16_max = 0.0, ipb_16_c_max = 0.0, max_rate8 = 0.0, max_rate16 = 0.0;
    double ipb_8_min = INFINITY, ipb_8_c_min = INFINITY, ipb_16_min = INFINITY, ipb_16_c_min = INFINITY, min_rate8 = INFINITY, min_rate16 = INFINITY;

    posit32_t tmp8 = convertDoubleToP32(CONST_IPB_8);
    posit32_t tmp16 = convertDoubleToP32(CONST_IPB_16);
    for (int i = 0; i < total; i ++)
    {
        p8_c[i] = convertDoubleToP8(test_int8[i].center);
        p8_r[i] = convertDoubleToP8(test_int8[i].radius);
        posit32_t tmp = p8_to_p32(p8_r[i]);
        ipb_8[i] = p32_div(tmp8, tmp);

        p16_c[i] = convertDoubleToP16(test_int16[i].center);
        p16_r[i] = convertDoubleToP16(test_int16[i].radius);
        posit32_t tmp2 = p16_to_p32(p16_r[i]);
        ipb_16[i] = p32_div(tmp16, tmp2);
    }


    double start1 = get_time_ms();
    for (int i = 0; i < total; i ++)
    {
        p16_interval[i] = p8_compression(p8_c[i], p8_r[i]);
    }
    double end1 = get_time_ms();

    double start2 = get_time_ms();
    for (int i = 0; i < total; i ++)
    {
        p32_interval[i] = p16_compression(p16_c[i], p16_r[i]);
    }
    double end2 = get_time_ms();

    posit16_t * p16_r_tilde = malloc(total * sizeof(posit16_t));
    posit32_t * p32_r_tilde = malloc(total * sizeof(posit32_t));
    
    double start3 = get_time_ms();
    for (int i = 0; i < total; i ++)
    {
        p16_r_tilde[i] = p16_read_3r(p16_interval[i].posit16);
    }
    double end3 = get_time_ms();


    double start4 = get_time_ms();
    for (int i = 0; i < total; i ++)
    {
        p32_r_tilde[i] = p32_read_3r(p32_interval[i].posit32);
    }
    double end4 = get_time_ms();

    for (int i = 0; i < total; i ++)
    {
        double factor16 = 0.5 / p16_interval[i].prec;
        posit32_t tmp1 = convertDoubleToP32(factor16);
        posit32_t tmp2 = p16_to_p32(p16_r_tilde[i]);
        ipb_8_convert[i] = p32_div(tmp1, tmp2); 

        double factor32 = 0.5 / p32_interval[i].prec;
        posit32_t tmp3 = convertDoubleToP32(factor32);
        ipb_16_convert[i] = p32_div(tmp3, p32_r_tilde[i]); 

        double ipb8_d = convertP32ToDouble(ipb_8[i]);
        double ipb8_c_d = convertP32ToDouble(ipb_8_convert[i]);
        double temp8 = ipb8_c_d / ipb8_d;
        ipb_8_avg += ipb8_d;
        ipb_8_c_avg += ipb8_c_d;
        avg_rate8 += temp8;
        ipb_8_max = fmax(ipb_8_max, ipb8_d);
        ipb_8_c_max = fmax(ipb_8_c_max, ipb8_c_d);
        max_rate8 = fmax(max_rate8, temp8);
        ipb_8_min = fmin(ipb_8_min, ipb8_d);
        ipb_8_c_min = fmin(ipb_8_c_min, ipb8_c_d);
        min_rate8 = fmin(min_rate8, temp8);

        double ipb16_d = convertP32ToDouble(ipb_16[i]);
        double ipb16_c_d = convertP32ToDouble(ipb_16_convert[i]);
        double temp16 = ipb16_c_d / ipb16_d;
        ipb_16_avg += ipb16_d;
        ipb_16_c_avg += ipb16_c_d;
        avg_rate16 += temp16;
        ipb_16_max = fmax(ipb_16_max, ipb16_d);
        ipb_16_c_max = fmax(ipb_16_c_max, ipb16_c_d);
        max_rate16 = fmax(max_rate16, temp16);
        ipb_16_min = fmin(ipb_16_min, ipb16_d);
        ipb_16_c_min = fmin(ipb_16_c_min, ipb16_c_d);
        min_rate16 = fmin(min_rate16, temp16);
    }
    ipb_8_avg /= total;
    ipb_8_c_avg /= total;
    avg_rate8 /= total;
    ipb_16_avg /= total;
    ipb_16_c_avg /= total;
    avg_rate16 /= total;

    /*
    for (int j = 0; j < total; j += 100)
    {
        printf("%d\n\n", j);

        printf("origin:\nc8 = ");
        print_binary(test_int8[j].center);
        printf("r8 = ");
        print_binary(test_int8[j].radius);
        printf("p8_c = ");
        p8_print_binary(p8_c[j]);
        printf("p8_r = ");
        p8_print_binary(p8_r[j]);
        double ipb8_d = convertP32ToDouble(ipb_8[j]);
        printf("ipb8: %.17e\n", ipb8_d);

        printf("convert:\nc16 = ");
        p16_print_binary(p16_interval[j].posit16);
        double c16 = convertP16ToDouble(p16_interval[j].posit16);
        print_binary(c16);
        printf("precision16: %d\n", p16_interval[j].prec);
        double ipb8_c_d = convertP32ToDouble(ipb_8_convert[j]);
        printf("ipb8 convert: %.17e\n", ipb8_c_d);

        printf("read radius:\nr16 = ");
        p16_print_binary(p16_r_tilde[j]);
        double r16 = convertP16ToDouble(p16_r_tilde[j]);
        print_binary(r16);
        printf("\n");

        printf("origin:\nc16 = ");
        print_binary(test_int16[j].center);
        printf("r16 = ");
        print_binary(test_int16[j].radius);
        printf("p16_c = ");
        p16_print_binary(p16_c[j]);
        printf("p16_r = ");
        p16_print_binary(p16_r[j]);
        double ipb16_d = convertP32ToDouble(ipb_16[j]);
        printf("ipb16: %.17e\n", ipb16_d);
        
        printf("convert:\nc32 = ");
        p32_print_binary(p32_interval[j].posit32);
        double c32 = convertP32ToDouble(p32_interval[j].posit32);
        print_binary(c32);
        printf("precision32: %d\n", p32_interval[j].prec);
        double ipb16_c_d = convertP32ToDouble(ipb_16_convert[j]);
        printf("ipb16 convert: %.17e\n", ipb16_c_d);
        
        printf("read radius:\nr32 = ");
        p32_print_binary(p32_r_tilde[j]);
        double r32 = convertP32ToDouble(p32_r_tilde[j]);
        print_binary(r32);
        printf("\n");
    }
    */
    

    FILE *fp = fopen("posit_convert_time.txt", "a");
    if (fp == NULL) 
    {
        perror("Error opening file: convert time");
        return 1;
    }
    fprintf(fp, "%.17e\t%.17e\t%.17e\t%.17e\t(ms)\n",
            end1 - start1, end2 - start2, end3 - start3, end4 - start4);

    fclose(fp);

    FILE *fp2 = fopen("posit_ipb.txt", "a");
    if (fp2 == NULL)
    {
        perror("Error opening file: ipb");
        return 1;
    } 
    fprintf(fp2, "%.17e\t%.17e\t%.17e\t%.17e\t%.17e\t%.17e\t%.17e\t%.17e\t%.17e\t%.17e\t%.17e\t%.17e\n",
            ipb_8_avg, ipb_8_min, ipb_8_max, ipb_8_c_avg, ipb_8_c_min, ipb_8_c_max, 
            ipb_16_avg, ipb_16_min, ipb_16_max, ipb_16_c_avg, ipb_16_c_min, ipb_16_c_max);
    fclose(fp2);

    FILE *fp3 = fopen("posit_ipb_rate.txt", "a");
    if (fp3 == NULL)
    {
        perror("Error opening file: ipb rate");
        return 1;
    } 
    fprintf(fp3, "%.17e\t%.17e\t%.17e\t%.17e\t%.17e\t%.17e\n",
        avg_rate8, min_rate8, max_rate8, 
        avg_rate16, min_rate16, max_rate16);
    fclose(fp3);
    
    free(test_int8);
    free(test_int16);
    free(p8_c);
    free(p8_r);
    free(p16_c);
    free(p16_r);
    free(p16_interval);
    free(p32_interval);
    free(p16_r_tilde);
    free(p32_r_tilde);
    free(ipb_8);
    free(ipb_8_convert);
    free(ipb_16);
    free(ipb_16_convert);

    return 0;
}
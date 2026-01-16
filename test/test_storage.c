#include "convert.h"

#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include <math.h>
#include <mpfr.h>


#define N 20000


int n_bits_u32(int precision)
{
    return (precision <= 20) ? 33 : 65;
}

int n_bits_fd(int precision)
{
    return (precision <= 23) ? 33 : 65;
}

int n_bits_mixed(int precision)
{
    return int_max(20, ((precision + 7 )/ 8) * 8 + 4);
}

double info_per_bit(double radius, int n_bits)
{
    return 0.5 / (radius * n_bits);
}

int main(int argc, char ** argv)
{
    int total = 2 * N;
    C_R *test_int2 = malloc(total * sizeof(C_R));

    int n_mask = (total + 7) >> 3; // ceil(total / 8)
    uint8_t *mask = malloc(n_mask * sizeof(uint8_t));
    int count_1 = 0;

    uint8_t *mask2 = malloc(n_mask * sizeof(uint8_t));
    int count_2 = 0;

    uint8_t *mask3 = malloc(N * sizeof(uint8_t));
    int n_double = 0, count_u32 = 0, count_u16 = 0, count_u8 = 0;

    FP_INT_PREC *cp_fp = malloc(total * sizeof(FP_INT_PREC));
    double *val1 = malloc(total * sizeof(double)); 
    double *val2 = malloc(total * sizeof(double));
    double *val3 = malloc(total * sizeof(double));

    FP_INT_PREC *cp_mpfr = malloc(total * sizeof(FP_INT_PREC));
    mpfr_t *mpfr_c = malloc(total * sizeof(mpfr_t));

    const double c_exp_min = -6.0;
    const double c_exp_max = 9.0;
    const double r_exp_min = -53.0;
    const double r_exp_max = -1.0;

    srand(get_time_ms());  

    for (int i = 0; i < N; i ++)
    {
        double c_power1 = c_exp_min + ((double)rand() / RAND_MAX) * (c_exp_max - c_exp_min);
        double c1 = ((double)rand() / RAND_MAX) * pow(10, c_power1);   
        double r_power1 = r_exp_min + ((double)rand() / RAND_MAX) * (r_exp_max - r_exp_min);
        double r1 = pow(2, r_power1) * fabs(c1);    
        test_int2[i].center = c1;
        test_int2[i].radius = r1;

        double c_power2 = c_exp_min + ((double)rand() / RAND_MAX) * (c_exp_max - c_exp_min);
        double c2 =  - ((double)rand() / RAND_MAX) * pow(10, c_power2);   
        double r_power2 = r_exp_min + ((double)rand() / RAND_MAX) * (r_exp_max - r_exp_min);
        double r2 = pow(2, r_power2) * fabs(c2);     
        test_int2[i + N].center = c2;
        test_int2[i + N].radius = r2;
    }

    // double to mpfr_t converssion
    double start1 = get_time_ms();
    for (int i = 0; i < total; i ++)
    {
        cp_mpfr[i] = CR_FP_mpfr(test_int2[i]);
        mpfr_init2(mpfr_c[i], cp_mpfr[i].precision);
        mpfr_set_d(mpfr_c[i], cp_mpfr[i].center, MPFR_RNDZ);
    }
    double end1 = get_time_ms();

    // double fp-int conversion
    double start2 = get_time_ms();
    for (int i = 0; i < total; i ++)
    {
        cp_fp[i] = CR_FP_p(test_int2[i]);
    }
    double end2 = get_time_ms();


    // information per bit (before and after fp-int conversion + improved storage, and rate)
    double max_ipb_i = 0.0, min_ipb_i = INFINITY, mean_ipb_i = 0.0;

    double max_ipb_u32 = 0.0, min_ipb_u32 = INFINITY, mean_ipb_u32 = 0.0;
    double max_ipb_fd = 0.0, min_ipb_fd = INFINITY, mean_ipb_fd = 0.0;
    double max_ipb_mixed = 0.0, min_ipb_mixed = INFINITY, mean_ipb_mixed = 0.0;

    double max_rate_u32 = 0.0, mean_rate_u32 = 0.0, min_rate_u32 = INFINITY, 
        max_rate_fd = 0.0, mean_rate_fd = 0.0, min_rate_fd = INFINITY,
        max_rate_mixed = 0.0, mean_rate_mixed = 0.0, min_rate_mixed = INFINITY;

    int n_bits_i = 128;

    for (int i = 0; i < total; i ++)
    {
        double ipb_i = info_per_bit(test_int2[i].radius, n_bits_i);
        C_R new = FP_CR3(cp_fp[i].center);

        int n_u32 = n_bits_u32(cp_fp[i].precision);
        int n_fd = n_bits_fd(cp_fp[i].precision);
        int n_mixed = n_bits_mixed(cp_fp[i].precision);

        double ipb_u32 = info_per_bit(new.radius, n_u32);
        double ipb_fd = info_per_bit(new.radius, n_fd);
        double ipb_mixed = info_per_bit(new.radius, n_mixed);

        double rate_u32 = ipb_u32 / ipb_i;
        double rate_fd = ipb_fd / ipb_i;
        double rate_mixed = ipb_mixed / ipb_i;

        mean_ipb_i += ipb_i;
        mean_ipb_u32 += ipb_u32;
        mean_ipb_fd += ipb_fd;
        mean_ipb_mixed += ipb_mixed;

        max_ipb_i = fmax(max_ipb_i, ipb_i);
        max_ipb_u32 = fmax(max_ipb_u32, ipb_u32);
        max_ipb_fd = fmax(max_ipb_fd, ipb_fd);
        max_ipb_mixed = fmax(max_ipb_mixed, ipb_mixed);

        min_ipb_i = fmin(min_ipb_i, ipb_i);
        min_ipb_u32 = fmin(min_ipb_u32, ipb_u32);
        min_ipb_fd = fmin(min_ipb_fd, ipb_fd);
        min_ipb_mixed = fmin(min_ipb_mixed, ipb_mixed);

        mean_rate_u32 += rate_u32;
        mean_rate_fd += rate_fd;
        mean_rate_mixed += rate_mixed;

        max_rate_u32 = fmax(max_rate_u32, rate_u32);
        max_rate_fd = fmax(max_rate_fd, rate_fd);
        max_rate_mixed = fmax(max_rate_mixed, rate_mixed);

        min_rate_u32 = fmin(min_rate_u32, rate_u32);
        min_rate_fd = fmin(min_rate_fd, rate_fd);
        min_rate_mixed = fmin(min_rate_mixed, rate_mixed);
    }
    mean_ipb_i /= total;
    mean_ipb_u32 /= total;
    mean_ipb_fd /= total;
    mean_ipb_mixed /= total;

    mean_rate_u32 /= total;
    mean_rate_fd /= total;
    mean_rate_mixed /= total;


    // set mask bits: 1 if precision <= 20, 0 otherwise
    double mask1start = get_time_ms();
    for (int i = 0; i < n_mask; i ++)
    {
        mask[i] = 0;
        for (int j = 0; j < 8; j ++)
        {
            int idx = i * 8 + j;
            if (idx >= total) break;

            if (cp_fp[i * 8 + j].precision <= 20)
            {
                mask[i] |= (1 << j);
            }
        }
    }
    double mask1end = get_time_ms();

    // set mask2 bits: 1 if precision <= 23, 0 otherwise
    double mask2start = get_time_ms();
    for (int i = 0; i < n_mask; i ++)
    {
        mask2[i] = 0;
        for (int j = 0; j < 8; j ++)
        {
            int idx = i * 8 + j;
            if (idx >= total) break;

            if (cp_fp[i * 8 + j].precision <= 23)
            {
                mask2[i] |= (1 << j);
            }
        }
    }
    double mask2end = get_time_ms();

    // set mask3 bits
    double mask3start = get_time_ms();
    for (int i = 0; i < N; i ++)
    {
        int idx = i * 2;
        int idx2 = i * 2 + 1;
        uint8_t m1 = precision_to_mask3(cp_fp[idx].precision);
        uint8_t m2 = precision_to_mask3(cp_fp[idx2].precision);
        mask3[i] = (m1 << 4) | m2;
    }
    double mask3end = get_time_ms();

    // count number of 1s in mask and mask2
    for (int i = 0; i < n_mask; i ++)
    {
        for (int j = 0; j < 8; j ++) 
        {
            int idx = i * 8 + j;
            if (idx >= total) break;
            count_1 += (mask[i] >> j) & 1;
            count_2 += (mask2[i] >> j) & 1;
        }
    }

    // count number of n_double, uint32_t, uint16_t and uint8_t in mask3
    for (int i = 0; i < N; i ++)
    {
        uint8_t m = mask3[i];
        count_u8 += m & 1; 
        count_u16 += (m >> 1) & 1;
        count_u32 += (m >> 2) & 1;
        n_double += (m >> 3) & 1;
        count_u8 += (m >> 4) & 1;
        count_u16 += (m >> 5) & 1;
        count_u32 += (m >> 6) & 1;
        n_double += (m >> 7) & 1;
    }

    // allocate necessary memory
    int n_32 = N * 4 - count_1;
    uint32_t *res32 = malloc(n_32 * sizeof(uint32_t));

    struct D_F res1;
    res1.d = malloc((total - count_2) * sizeof(double));
    res1.f = malloc(count_2 * sizeof(float));

    struct MP res_mp;
    res_mp.d = malloc(n_double * sizeof(double));
    res_mp.u32 = malloc(count_u32 * sizeof(uint32_t));
    res_mp.u16 = malloc(count_u16 * sizeof(uint16_t));
    res_mp.u8 = malloc(count_u8 * sizeof(uint8_t));

    // copy reduced results: double to uint32_t
    int k = 0;  // res32 index

    double start3 = get_time_ms();
    for (int i = 0; i < n_mask; i ++) 
    {
        for (int j = 7; j >= 0; j --) 
        {
            int idx = i * 8 + j;
            if (idx >= total) continue;

            cast dc;
            dc.d = cp_fp[idx].center;

            if ((mask[i] >> j) & 1) 
            {
                res32[k ++] = (uint32_t)(dc.u64 >> 32);
            } 
            else 
            {
                res32[k ++] = (uint32_t)(dc.u64 >> 32);       // high part
                res32[k ++] = (uint32_t)(dc.u64 & 0xFFFFFFFF); // low part
            }
        }
    }
    double end3 = get_time_ms();

    // copy reduced results: double to double and float
    int id_f = 0; // index for float array
    int id_d = 0; // index for double array

    double start3b = get_time_ms();
    for (int i = 0; i < n_mask; i ++) 
    {
        for (int j = 7; j >= 0; j --) 
        {
            int idx = i * 8 + j;
            if (idx >= total) continue;

            if ((mask2[i] >> j) & 1) 
            {
                res1.f[id_f ++] = (float)cp_fp[idx].center;
            } 
            else 
            {
                res1.d[id_d ++] = cp_fp[idx].center;
            }
        }
    }
    double end3b = get_time_ms();

    // copy reduced results: double to double, uint32_t, uint16_t and uint8_t
    int id_d2 = 0; // index for double array
    int id_u32 = 0; // index for uint32_t array
    int id_u16 = 0; // index for uint16_t array
    int id_u8 = 0; // index for uint8_t array

    double start3c = get_time_ms();
    for (int i = 0; i < N; i ++) 
    {
        uint8_t m2 = mask3[i] & 0x0F;
        uint8_t m1 = mask3[i] >> 4;
        int ind1 = i * 2;
        int ind2 = i * 2 + 1;
        cast dc;

        dc.d = cp_fp[ind1].center;

        if (m1 == 0x08) 
        {
            res_mp.d[id_d2 ++] = dc.d;
        }  
        else 
        {
            if (m1 & 0x04) 
            {  // u32
                res_mp.u32[id_u32 ++] = (uint32_t)(dc.u64 >> 32);
                if (m1 & 0x02) 
                {  // u16
                    res_mp.u16[id_u16 ++] = (uint16_t)(dc.u64 >> 16);
                    if (m1 & 0x01) 
                    {  // u8
                        res_mp.u8[id_u8 ++] = (uint8_t)(dc.u64 >> 8);
                    }
                } 
                else if (m1 & 0x01) 
                {  // u32 + u8
                    res_mp.u8[id_u8 ++] = (uint8_t)(dc.u64 >> 24);
                }
            } 
            else if (m1 & 0x02) 
            {  // u16 only
                res_mp.u16[id_u16 ++] = (uint16_t)(dc.u64 >> 48);
                if (m1 & 0x01) 
                {  // u16 + u8
                    res_mp.u8[id_u8 ++] = (uint8_t)(dc.u64 >> 40);
                }
            } 
            else if (m1 & 0x01) 
            {  // u8 only
                res_mp.u8[id_u8 ++] = (uint8_t)(dc.u64 >> 56);
            }
        }

        dc.d = cp_fp[ind2].center;

        if (m2 == 0x08) 
        {
            res_mp.d[id_d2 ++] = dc.d;
        } 
        else 
        {
            if (m2 & 0x04) 
            {
                res_mp.u32[id_u32 ++] = (uint32_t)(dc.u64 >> 32);
                if (m2 & 0x02) 
                {
                    res_mp.u16[id_u16 ++] = (uint16_t)(dc.u64 >> 16);
                    if (m2 & 0x01) 
                    {
                        res_mp.u8[id_u8 ++] = (uint8_t)(dc.u64 >> 8);
                    }
                } 
                else if (m2 & 0x01) 
                {
                    res_mp.u8[id_u8 ++] = (uint8_t)(dc.u64 >> 24);
                }
            } 
            else if (m2 & 0x02) 
            {
                res_mp.u16[id_u16 ++] = (uint16_t)(dc.u64 >> 48);
                if (m2 & 0x01) 
                {
                    res_mp.u8[id_u8 ++] = (uint8_t)(dc.u64 >> 40);
                }
            } 
            else if (m2 & 0x01) 
            {
                res_mp.u8[id_u8 ++] = (uint8_t)(dc.u64 >> 56);
            }
        }
    }
    double end3c = get_time_ms();
    

    // read uint32_t to double
    k = 0;

    double start4 = get_time_ms();
    for (int i = 0; i < n_mask; i ++) 
    {
        for (int j = 7; j >= 0; j --) 
        {
            int idx = i * 8 + j;
            if (idx >= total) continue;

            uint64_t high = (uint64_t)res32[k ++] << 32;
            uint64_t low = 0;

            if (!((mask[i] >> j) & 1)) 
            {
                low = (uint64_t)res32[k++];
            }

            uint64_t combined = high | low;

            cast cast;
            cast.u64 = combined;

            val1[idx] = cast.d;
        }
    }
    double end4 = get_time_ms();


    // read double and float to double
    id_f = 0;
    id_d = 0;

    double start4b = get_time_ms();
    for (int i = 0; i < n_mask; i ++) 
    {
        for (int j = 7; j >= 0; j --) 
        {
            int idx = i * 8 + j;
            if (idx >= total) continue;

            if ((mask2[i] >> j) & 1) 
            {
                val2[idx] = (double)res1.f[id_f ++];
            } 
            else 
            {
                val2[idx] = res1.d[id_d ++];
            }
        }
    }
    double end4b = get_time_ms();

    // read double, uint32_t, uint16_t and uint8_t to double
    id_d2 = 0;
    id_u32 = 0;
    id_u16 = 0;
    id_u8 = 0;

    double start4c = get_time_ms();
    for (int i = 0; i < N; i ++) 
    {
        uint8_t m2 = mask3[i] & 0x0F;
        uint8_t m1 = mask3[i] >> 4;
        int ind1 = i * 2;
        int ind2 = i * 2 + 1;

        cast dc;

        if (m1 == 0x08) 
        {
            val3[ind1] = res_mp.d[id_d2 ++];
        } 
        else 
        {
            dc.u64 = 0;
            if (m1 & 0x04) 
            {
                dc.u64 |= (uint64_t)res_mp.u32[id_u32 ++] << 32;
                if (m1 & 0x02) 
                {
                    dc.u64 |= (uint64_t)res_mp.u16[id_u16 ++] << 16;
                    if (m1 & 0x01)
                        dc.u64 |= (uint64_t)res_mp.u8[id_u8 ++] << 8;
                } 
                else if (m1 & 0x01) 
                {
                    dc.u64 |= (uint64_t)res_mp.u8[id_u8 ++] << 24;
                }
            } 
            else if (m1 & 0x02) 
            {
                dc.u64 |= (uint64_t)res_mp.u16[id_u16 ++] << 48;
                if (m1 & 0x01)
                    dc.u64 |= (uint64_t)res_mp.u8[id_u8 ++] << 40;
            } 
            else if (m1 & 0x01) 
            {
                dc.u64 |= (uint64_t)res_mp.u8[id_u8 ++] << 56;
            }
            val3[ind1] = dc.d;
        }

        if (m2 == 0x08) 
        {
            val3[ind2] = res_mp.d[id_d2 ++];
        } 
        else 
        {
            dc.u64 = 0;
            if (m2 & 0x04) 
            {
                dc.u64 |= (uint64_t)res_mp.u32[id_u32 ++] << 32;
                if (m2 & 0x02) 
                {
                    dc.u64 |= (uint64_t)res_mp.u16[id_u16 ++] << 16;
                    if (m2 & 0x01)
                        dc.u64 |= (uint64_t)res_mp.u8[id_u8 ++] << 8;
                } 
                else if (m2 & 0x01) 
                {
                    dc.u64 |= (uint64_t)res_mp.u8[id_u8 ++] << 24;
                }
            } 
            else if (m2 & 0x02) 
            {
                dc.u64 |= (uint64_t)res_mp.u16[id_u16 ++] << 48;
                if (m2 & 0x01)
                    dc.u64 |= (uint64_t)res_mp.u8[id_u8 ++] << 40;
            } 
            else if (m2 & 0x01) 
            {
                dc.u64 |= (uint64_t)res_mp.u8[id_u8 ++] << 56;
            }
            val3[ind2] = dc.d;
        }
    }
    double end4c = get_time_ms();


    // time
    FILE *fp = fopen("compression_time.txt", "a");
    fprintf(fp, "%.17e\t%.17e\t%.17e\t%.17e\t%.17e\t%.17e\t%.17e\t%.17e\t%.17e\t%.17e\t%.17e\n", 
            end1 - start1, end2 - start2, end3 - start3, end4 - start4, 
            end3b - start3b, end4b - start4b, 
            end3c - start3c, end4c - start4c, 
            mask1end - mask1start, mask2end - mask2start, mask3end - mask3start);
    fclose(fp);

    // information per bit
    FILE *fp2 = fopen("compression_ipb.txt", "a");
    fprintf(fp2, "%.17e\t%.17e\t%.17e\t%.17e\t%.17e\t%.17e\t%.17e\t%.17e\t%.17e\t%.17e\t%.17e\t%.17e\n", 
        mean_ipb_i, min_ipb_i, max_ipb_i, 
        mean_ipb_u32, min_ipb_u32,max_ipb_u32, 
        mean_ipb_fd, min_ipb_fd, max_ipb_fd, 
        mean_ipb_mixed, min_ipb_mixed, max_ipb_mixed);
    fclose(fp2);

    FILE *fp3 = fopen("compression_ipb_rate.txt", "a");
    fprintf(fp3, "%.17e\t%.17e\t%.17e\t%.17e\t%.17e\t%.17e\t%.17e\t%.17e\t%.17e\n", 
        mean_rate_u32, min_rate_u32, max_rate_u32, 
        mean_rate_fd, min_rate_fd, max_rate_fd, 
        mean_rate_mixed, min_rate_mixed, max_rate_mixed);
    fclose(fp3);

    for (int i = 0; i < total; i ++)
    {
        if (fabs(val1[i] - cp_fp[i].center) > 1e-10)  
        {
            printf("Error: %d\n", i);
            printf("Original: %.17e\n", cp_fp[i].center);
            printf("Recovered1: %.17e\n", val1[i]);
        }
        if (fabs(val2[i] - cp_fp[i].center) > 1e-10)
        {
            printf("Error: %d\n", i);
            printf("Original: %.17e\n", cp_fp[i].center);
            printf("Recovered2: %.17e\n", val2[i]);
        }
        if (fabs(val3[i] - cp_fp[i].center) > 1e-10)
        {
            printf("Error: %d\n", i);
            printf("Original: %.17e\n", cp_fp[i].center);
            printf("Recovered3: %.17e\n", val3[i]);
        }
    }

    printf("%d\t%d\t%d\t%d\t%d\t%d\t%d\n", n_mask, count_1, count_2, n_double, count_u32, count_u16, count_u8);

    for (int i = 0; i < total; i ++)
    {
        mpfr_clear(mpfr_c[i]);
    }

    free(cp_mpfr);
    free(test_int2);
    free(mpfr_c);
    free(mask);
    free(res32);
    free(cp_fp);
    free(mask2);
    free(res1.d);
    free(res1.f);
    free(val1);
    free(val2);
    free(val3);
    free(mask3);
    free(res_mp.d);
    free(res_mp.u32);
    free(res_mp.u16);
    free(res_mp.u8);

    mpfr_free_cache();

    return 0;
}
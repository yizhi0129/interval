#include "posit_interval.h"
#include "convert.h"

#include "internals.h"
#include "softposit.h"

void p8_print_binary(posit8_t x)
{
    uint8_t y = x.v;
    for (int i = 7; i >= 0; i --)
    {
        int b = (y >> i) & 1;
        printf("%d", b);
    }
    printf("\n");
}

void p16_print_binary(posit16_t x)
{
    uint16_t y = x.v;
    for (int i = 15; i >= 0; i --)
    {
        int b = (y >> i) & 1;
        printf("%d", b);
    }
    printf("\n");
}

void p32_print_binary(posit32_t x)
{
    uint32_t y = x.v;
    for (int i = 31; i >= 0; i --)
    {
        int b = (y >> i) & 1;
        printf("%d", b);
    }
    printf("\n");
}


posit8_t p8_read_3r(posit8_t c)
{
    union ui8_p8 uZ;
	double d8;
	uZ.p = c;

	if (uZ.ui == 0 || uZ.ui==0x7F || uZ.ui==0x81 || uZ.ui == 0x80)
    {
		printf("Error: not valid center!\n");
        return uZ.p;
	}

    bool regS, sign;
	uint_fast8_t shift = 2, reg, frac;
	int_fast8_t k = 0;

	sign = signP8UI(uZ.ui);
	if (sign) uZ.ui = -uZ.ui & 0xFF;
	regS = signregP8UI(uZ.ui);

	uint_fast8_t tmp = (uZ.ui<<2) & 0xFF;
	if (regS)
    {
		while (tmp>>7)
        {
			k ++;
			shift ++;
			tmp = (tmp<<1) & 0xFF;
		}
		reg = k + 1;
	}
	else
    {
		k = -1;
		while (!(tmp>>7))
        {
			k --;
			shift ++;
			tmp = (tmp<<1) & 0xFF;
		}
		reg = -k;
		tmp &= 0x7F;
	}
    frac = (tmp & 0x7F) >> shift;

    int frac_bits = 6 - reg;

    if (frac == 0) 
    {
        d8 = 3.0 * set_pow2(k);
        return convertDoubleToP8(d8);
    }

    int rightmost_pos = -1;
    for (int i = 0; i < frac_bits; i ++) 
	{
        if ((frac >> i) & 0x1) 
		{
            rightmost_pos = i;
            break;
        }
    }

	if (rightmost_pos == -1) 
	{
        double d8 = 3.0 * set_pow2(k);
        return convertDoubleToP8(d8);
    }

    int total_exp = k - frac_bits + rightmost_pos;

    d8 = 3.0 * set_pow2(total_exp);
    return convertDoubleToP8(d8);
}

posit16_t p16_read_3r(posit16_t c)
{
    union ui16_p16 uZ;
	double d16;
	uZ.p = c;

	if (uZ.ui == 0 || uZ.ui == 0x7FFF || uZ.ui==0x8001 || uZ.ui == 0x8000)
    {
        printf("Error: not valid center!\n");
		return uZ.p;
	}

	bool regS, sign;
	uint_fast16_t reg, shift = 2, frac;
	int_fast16_t k = 0;
	int_fast8_t exp;

	sign = signP16UI(uZ.ui);
	if (sign)
		uZ.ui = -uZ.ui & 0xFFFF;
	regS = signregP16UI(uZ.ui);

	uint_fast16_t tmp = (uZ.ui<<2) & 0xFFFF;
	if (regS)
    {
		while (tmp>>15)
        {
			k ++;
			shift ++;
			tmp= (tmp<<1) & 0xFFFF;
		}
		reg = k + 1;
	}
	else
    {
		k = -1;
		while (!(tmp>>15))
        {
			k --;
			shift ++;
			tmp = (tmp<<1) & 0xFFFF;
		}
		reg = -k;
		tmp &= 0x7FFF;
	}
	exp = tmp>>14;
	frac = (tmp & 0x3FFF) >> shift;

    int frac_bits = 13 - reg;

    if (frac == 0) 
    {
        d16 = 3.0 * set_pow2(k);
        return convertDoubleToP16(d16);
    }

    int rightmost_pos = -1;
    for (int i = 0; i < frac_bits; i ++) 
	{
        if ((frac >> i) & 0x1) 
		{
            rightmost_pos = i;
            break;
        }
    }

	if (rightmost_pos == -1) 
	{
        double d16 = 3.0 * set_pow2(k);
        return convertDoubleToP16(d16);
    }

    int total_exp = 2 * k + exp - frac_bits + rightmost_pos;
    d16 = 3.0 * set_pow2(total_exp);
    return convertDoubleToP16(d16);
}

posit32_t p32_read_3r(posit32_t c)
{
    union ui32_p32 uZ;
	double d32;
	uZ.p = c;

	if (uZ.ui == 0 || uZ.ui==0x7FFFFFFF || uZ.ui==0x80000001 || uZ.ui == 0x80000000)
    {
		printf("Error: not valid center!\n");
        return uZ.p;
	}

	bool regS, sign;
	uint_fast32_t reg, shift = 2, frac, tmp;
	int_fast32_t k = 0;
	int_fast8_t exp;

	sign = signP32UI(uZ.ui);
	if (sign)
		uZ.ui = -uZ.ui & 0xFFFFFFFF;
	regS = signregP32UI(uZ.ui);

	tmp = (uZ.ui<<2) & 0xFFFFFFFF;
	if (regS)
    {
		while (tmp>>31)
        {
			k ++;
			shift ++;
			tmp= (tmp<<1) & 0xFFFFFFFF;
		}
		reg = k + 1;
	}
	else
    {
		k = -1;
		while (!(tmp>>31))
        {
			k --;
			shift ++;
			tmp = (tmp<<1) & 0xFFFFFFFF;
		}
		reg = -k;
		tmp &= 0x7FFFFFFF;
	}
	exp = tmp>>29;

	frac = (tmp & 0x1FFFFFFF) >> shift;

    int frac_bits = 28 - reg;

    if (frac == 0) 
    {
        d32 = 3.0 * set_pow2(4 * k + exp);
        return convertDoubleToP32(d32);
    }

    int rightmost_pos = -1;
    for (int i = 0; i < frac_bits; i ++) 
	{
        if ((frac >> i) & 0x1) 
		{
            rightmost_pos = i;
            break;
        }
    }

	if (rightmost_pos == -1) 
	{
        double d32 = 3.0 * set_pow2(k);
        return convertDoubleToP32(d32);
    }

    int total_exp = 4 * k + exp - frac_bits + rightmost_pos;
    d32 = 3.0 * set_pow2(total_exp);
    return convertDoubleToP32(d32);
}


void p8_read_power(posit8_t x, bool * sign, int8_t * k, uint8_t * reg, uint8_t * n_frac)
{
	union ui8_p8 uZ;
	uZ.p = x;

	if (uZ.ui == 0x80 || uZ.ui == 0x00) 
	{
        *sign = false;
        *k = 0;
        *reg = 0;
        *n_frac = 0;
        return;
    }

	bool regS;
	*k = 0;

	*sign = signP8UI(uZ.ui);
	if (*sign) uZ.ui = -uZ.ui & 0xFF;
	regS = signregP8UI(uZ.ui);

	uint_fast8_t tmp = (uZ.ui<<2) & 0xFF;
	if (regS)
    {
		while (tmp>>7)
        {
			(*k) ++;
			tmp = (tmp<<1) & 0xFF;
		}
		*reg = (*k) + 1;
	}
	else
    {
		*k = -1;
		while (!(tmp>>7))
        {
			(*k) --;
			tmp = (tmp<<1) & 0xFF;
		}
		tmp &= 0x7F;
		*reg = -(*k);
	}
	*n_frac = 6 - (*reg);
}

void p16_read_power(posit16_t x, bool * sign, int8_t * k, uint8_t * reg, int8_t * exp, uint8_t * n_frac)
{
	union ui16_p16 uZ;
	uZ.p = x;

	if (uZ.ui == 0x8000 || uZ.ui == 0x0000) 
	{
        *sign = false;
        *k = 0;
        *reg = 0;
        *exp = 0;
        *n_frac = 0;
        return;
    }

	bool regS;
	*k = 0;

	*sign = signP16UI(uZ.ui);
	if (*sign)
		uZ.ui = -uZ.ui & 0xFFFF;
	regS = signregP16UI(uZ.ui);

	uint_fast16_t tmp = (uZ.ui<<2) & 0xFFFF;
	if (regS)
    {
		while (tmp>>15)
        {
			(*k) ++;
			tmp= (tmp<<1) & 0xFFFF;
		}
		*reg = (*k) + 1;
	}
	else
    {
		*k = -1;
		while (!(tmp>>15))
        {
			(*k) --;
			tmp = (tmp<<1) & 0xFFFF;
		}
		tmp &= 0x7FFF;
		*reg = -(*k);
	}
	*n_frac = 13 - (*reg);
	*exp = (tmp >> 14) & 0x1;
}

void p32_read_power(posit32_t x, bool * sign, int8_t * k, uint8_t * reg, int8_t * exp, uint8_t * n_frac)
{
	union ui32_p32 uZ;
	uZ.p = x;

	if (uZ.ui == 0x80000000 || uZ.ui == 0x00000000) 
	{
        *sign = false;
        *k = 0;
        *reg = 0;
        *exp = 0;
        *n_frac = 0;
        return;
    }

	bool regS;
	*k = 0;

	*sign = signP32UI(uZ.ui);
	if (*sign)
		uZ.ui = -uZ.ui & 0xFFFFFFFF;
	regS = signregP32UI(uZ.ui);

	uint_fast32_t tmp = (uZ.ui<<2) & 0xFFFFFFFF;
	if (regS)
	{
		while (tmp>>31)
		{
			(*k) ++;
			tmp= (tmp<<1) & 0xFFFFFFFF;
		}
		*reg = (*k) + 1;
	}
	else
	{
		*k = -1;
		while (!(tmp>>31))
		{
			(*k) --;
			tmp = (tmp<<1) & 0xFFFFFFFF;
		}
		tmp &= 0x7FFFFFFF;
		*reg = -(*k);
	}
	*n_frac = 28 - (*reg);
	*exp = (tmp >> 29) & 0x3;
}


posit8_t p8_truncate(posit8_t x, int index)
{
	union ui8_p8 uZ;
	uZ.p = x;

	if (index < 0 || index >= 8) return x;

	bool sign = signP8UI(uZ.ui);
	if (sign) uZ.ui = -uZ.ui & 0xFF;

	for (int i = 0; i < index; i ++)
	{
		uZ.ui &= ~(1 << i);
	}
	uZ.ui |= (1 << index);

	if (sign) uZ.ui = -uZ.ui & 0xFF;

	return uZ.p;
}

posit16_t p16_truncate(posit16_t x, int index)
{
	union ui16_p16 uZ;
	uZ.p = x;

	if (index < 0 || index >= 16) return x;

	bool sign = signP16UI(uZ.ui);
	if (sign) uZ.ui = -uZ.ui & 0xFFFF;

	for (int i = 0; i < index; i ++)
	{
		uZ.ui &= ~(1 << i);
	}
	uZ.ui |= (1 << index);

	if (sign) uZ.ui = -uZ.ui & 0xFFFF;

	return uZ.p;
}

posit32_t p32_truncate(posit32_t x, int index)
{
	union ui32_p32 uZ;
	uZ.p = x;

	if (index < 0 || index >= 32) return x;

	bool sign = signP32UI(uZ.ui);
	if (sign) uZ.ui = -uZ.ui & 0xFFFFFFFF;

	for (int i = 0; i < index; i ++)
	{
		uZ.ui &= ~(1 << i);
	}
	uZ.ui |= (1 << index);

	if (sign) uZ.ui = -uZ.ui & 0xFFFFFFFF;

	return uZ.p;
}


p16_FP_INT_P p8_compression(posit8_t c, posit8_t r)
{
	bool s_c, s_r;
	int8_t k_c, k_r;
	uint8_t reg_c, reg_r, n_frac_c, n_frac_r;
	p8_read_power(c, &s_c, &k_c, &reg_c, &n_frac_c);
	p8_read_power(r, &s_r, &k_r, &reg_r, &n_frac_r);
	int p = int_max(k_c - k_r, 0);
	if (p <= n_frac_c)
	{
		c = p8_truncate(c, n_frac_c - p);
		return (p16_FP_INT_P){.posit16 = p8_to_p16(c), .prec = 8};
	}
	else
	{
		posit16_t c_bis = p8_to_p16(c);
		int8_t exp;
		p16_read_power(c_bis, &s_c, &k_c, &reg_c, &exp, &n_frac_c);
		p = int_min(2 * k_c + exp - k_r, n_frac_c);
		c_bis = p16_truncate(c_bis, n_frac_c - p);
		return (p16_FP_INT_P){.posit16 = c_bis, .prec = 15 - n_frac_c + p};
	}
}

p32_FP_INT_P p16_compression(posit16_t c, posit16_t r)
{
	bool s_c, s_r;
	int8_t k_c, k_r, e_c, e_r;
	uint8_t reg_c, reg_r, n_frac_c, n_frac_r;
	p16_read_power(c, &s_c, &k_c, &reg_c, &e_c, &n_frac_c);
	p16_read_power(r, &s_r, &k_r, &reg_r, &e_r, &n_frac_r);
	int p = int_max(2 * k_c + e_c - 2 * k_r - e_r, 0);
	if (p <= n_frac_c)
	{
		c = p16_truncate(c, n_frac_c - p);
		return (p32_FP_INT_P){.posit32 = p16_to_p32(c), .prec = 16};
	}
	else
	{
		posit32_t c_bis = p16_to_p32(c);
		int8_t exp;
		p32_read_power(c_bis, &s_c, &k_c, &reg_c, &exp, &n_frac_c);
		p = int_min(4 * k_c + exp - 2 * k_r - e_r, n_frac_c);
		c_bis = p32_truncate(c_bis, n_frac_c - p);
		return (p32_FP_INT_P){.posit32 = c_bis, .prec = 31 - n_frac_c + p};
	}
}
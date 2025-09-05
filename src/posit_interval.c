#include "posit_interval.h"
#include "convert.h"

#include "internals.h"
#include "softposit.h"

#include <stdlib.h>
#include <math.h>

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

	if (uZ.ui == 0 || uZ.ui == MAXPOS16 || uZ.ui==0x8001 || uZ.ui == 0x8000)
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
		tmp &= MAXPOS16;
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

	if (uZ.ui == 0 || uZ.ui==MAXPOS32 || uZ.ui==0x80000001 || uZ.ui == NAR32)
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
		tmp &= MAXPOS32;
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
		tmp &= MAXPOS16;
		*reg = -(*k);
	}
	*n_frac = 13 - (*reg);
	*exp = (tmp >> 14) & 0x1;
}

void p32_read_power(posit32_t x, bool * sign, int8_t * k, uint8_t * reg, int8_t * exp, uint8_t * n_frac)
{
	union ui32_p32 uZ;
	uZ.p = x;

	if (uZ.ui == NAR32 || uZ.ui == 0x00000000) 
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
		tmp &= MAXPOS32;
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


// find next
posit8_t p8_nextup(posit8_t x)
{
	union ui8_p8 uZ;
    uZ.p = x;

    if (uZ.ui == NAR8 || uZ.ui == 0x00 || uZ.ui == MAXPOS8) return x;

    uZ.ui += 1;
    return uZ.p;
}

posit16_t p16_nextup(posit16_t x)
{
	union ui16_p16 uZ;
	uZ.p = x;

	if (uZ.ui == NAR16 || uZ.ui == 0x0000 || uZ.ui == MAXPOS16) return x;

	uZ.ui += 1;
	return uZ.p;
}

posit32_t p32_nextup(posit32_t x)
{
	union ui32_p32 uZ;
	uZ.p = x;

	if (uZ.ui == NAR32 || uZ.ui == 0x00000000 || uZ.ui == MAXPOS32) return x;

	uZ.ui += 1;
	return uZ.p;
}

posit8_t p8_nextdown(posit8_t x)
{
	union ui8_p8 uZ;
    uZ.p = x;

    if (uZ.ui == NAR8 || uZ.ui == 0x00 || uZ.ui == MINNEG8) return x;

    uZ.ui -= 1;
    return uZ.p;	
}

posit16_t p16_nextdown(posit16_t x)
{
	union ui16_p16 uZ;
    uZ.p = x;

    if (uZ.ui == NAR16 || uZ.ui == 0x0000 || uZ.ui == MINNEG16) return x;

    uZ.ui -= 1;
    return uZ.p;	
}

posit32_t p32_nextdown(posit32_t x)
{
	union ui32_p32 uZ;
    uZ.p = x;

    if (uZ.ui == NAR32 || uZ.ui == 0x00 || uZ.ui == MINNEG32) return x;

    uZ.ui -= 1;
    return uZ.p;	
}


// round up
posit16_t softposit_addMagsP16_up(uint_fast16_t uiA, uint_fast16_t uiB)
{
	uint_fast16_t regA, uiX, uiY;
	uint_fast32_t frac32A, frac32B;
	uint_fast16_t fracA=0,  regime, tmp;
	bool sign, regSA, regSB, rcarry=0, bitNPlusOne=0, bitsMore=0;
	int_fast8_t kA=0, expA;
	int_fast16_t shiftRight;
	union ui16_p16 uZ;

	sign = signP16UI( uiA ); //sign is always positive.. actually don't have to do this.
	if (sign){
		uiA = -uiA & 0xFFFF;
		uiB = -uiB & 0xFFFF;
	}

	if ((int_fast16_t)uiA < (int_fast16_t)uiB){
		uiX = uiA;
		uiY = uiB;
		uiA = uiY;
		uiB = uiX;
	}
	regSA = signregP16UI( uiA );
	regSB = signregP16UI( uiB );

	tmp = (uiA<<2) & 0xFFFF;
	if (regSA){
		while (tmp>>15){
			kA++;
			tmp= (tmp<<1) & 0xFFFF;
		}
	}
	else{
		kA=-1;
		while (!(tmp>>15)){
			kA--;
			tmp= (tmp<<1) & 0xFFFF;
		}
		tmp&=MAXPOS16;
	}
	expA = tmp>>14;
	frac32A = (0x4000 | tmp) << 16;
	shiftRight = kA;

	tmp = (uiB<<2) & 0xFFFF;
	if (regSB){
		while (tmp>>15){
			shiftRight--;
			tmp= (tmp<<1) & 0xFFFF;
		}
		frac32B = (0x4000 | tmp) <<16;
	}
	else{
		shiftRight++;
		while (!(tmp>>15)){
			shiftRight++;
			tmp= (tmp<<1) & 0xFFFF;
		}
		tmp&=MAXPOS16;
		frac32B = ( (0x4000 | tmp) <<16 ) & MAXPOS32;
	}

	//This is 2kZ + expZ; (where kZ=kA-kB and expZ=expA-expB)
	shiftRight = (shiftRight<<1) + expA - (tmp>>14);

	if (shiftRight==0){
		frac32A += frac32B;
		//rcarry is one
		if (expA) kA ++;
		expA^=1;
		frac32A>>=1;
	}
	else{
		//Manage CLANG (LLVM) compiler when shifting right more than number of bits
		(shiftRight>31) ? (frac32B=0): (frac32B >>= shiftRight); //frac32B >>= shiftRight

		frac32A += frac32B;
		rcarry = NAR32 & frac32A; //first left bit
		if(rcarry){
			if (expA) kA ++;
			expA^=1;
			frac32A>>=1;
		}
	}
	if(kA<0){
		regA = (-kA & 0xFFFF);
		regSA = 0;
		regime = 0x4000>>regA;
	}
	else{
		regA = kA+1;
		regSA=1;
		regime = MAXPOS16 - (MAXPOS16>>regA);
	}
	if(regA>14){
		//max or min pos. exp and frac does not matter.
		if (!sign) (regSA) ? (uZ.ui= MAXPOS16): (uZ.ui=0x1);
		else (regSA) ? (uZ.ui= MINNEG16): (uZ.ui=0x0);
	}
	else{
		//remove hidden bits
		frac32A = (frac32A & 0x3FFFFFFF) >>(regA + 1) ;
		fracA = frac32A>>16;
		if (regA!=14) bitNPlusOne = (frac32A>>15) & 0x1;
		else if (frac32A>0){
			fracA=0;
			bitsMore =1;
		}
		if (regA==14 && expA) bitNPlusOne = 1;
		uZ.ui = packToP16UI(regime, regA, expA, fracA);

		if ( frac32A&0x7FFF ) bitsMore=1;

		if (bitNPlusOne || bitsMore)
		{
			if (!sign) return p16_nextup(uZ.p);
		}
	}

	if (sign) uZ.ui = -uZ.ui & 0xFFFF;
	return uZ.p;
}

posit32_t softposit_addMagsP32_up(uint_fast32_t uiA, uint_fast32_t uiB) 
{
	uint_fast16_t regA, regB;
	uint_fast64_t frac64A=0, frac64B=0;
	uint_fast32_t fracA=0, regime, tmp;
	bool sign, regSA, regSB, rcarry=0, bitNPlusOne=0, bitsMore=0;
	int_fast8_t kA=0;
	int_fast32_t expA;
	int_fast16_t shiftRight;
	union ui32_p32 uZ;

	sign = signP32UI( uiA );
	if (sign){
		uiA = -uiA & 0xFFFFFFFF;
		uiB = -uiB & 0xFFFFFFFF;
	}

	if ((int_fast32_t)uiA < (int_fast32_t)uiB){
		uiA ^= uiB;
		uiB ^= uiA;
		uiA ^= uiB;
	}
	regSA = signregP32UI( uiA );
    regSB = signregP32UI( uiB );

    tmp = (uiA<<2)&0xFFFFFFFF;
	if (regSA){
		while (tmp>>31){
			kA++;
			tmp= (tmp<<1) & 0xFFFFFFFF;
		}
	}
	else{
		kA=-1;
		while (!(tmp>>31)){
			kA--;
			tmp= (tmp<<1) & 0xFFFFFFFF;
		}
		tmp&=MAXPOS32;
	}

	expA = tmp>>29; //to get 2 bits
	frac64A = ((0x40000000ULL | tmp<<1) & 0x7FFFFFFFULL) <<32;
	shiftRight = kA;

	tmp = (uiB<<2) & 0xFFFFFFFF;
	if (regSB){
		while (tmp>>31){
			shiftRight--;
			tmp= (tmp<<1) & 0xFFFFFFFF;
		}
	}
	else{
		shiftRight++;
		while (!(tmp>>31)){
			shiftRight++;
			tmp= (tmp<<1) & 0xFFFFFFFF;
		}
		tmp&=MAXPOS32;
	}
	frac64B = ((0x40000000ULL | tmp<<1) & 0x7FFFFFFFULL) <<32;
	//This is 4kZ + expZ; (where kZ=kA-kB and expZ=expA-expB)
	shiftRight = (shiftRight<<2) + expA - (tmp>>29);

	//Manage CLANG (LLVM) compiler when shifting right more than number of bits
	(shiftRight>63) ? (frac64B=0): (frac64B >>= shiftRight); //frac64B >>= shiftRight

	frac64A += frac64B;

	rcarry = NAR64 & frac64A; //first left bit
	if (rcarry){
		expA++;
		if (expA>3){
			kA ++;
			expA&=0x3;
		}
		frac64A>>=1;
	}
	if(kA<0){
		regA = -kA;
		regSA = 0;
		regime = 0x40000000>>regA;
	}
	else{
		regA = kA+1;
		regSA=1;
		regime = MAXPOS32 - (MAXPOS32>>regA);
	}

	if(regA>30){
		//max or min pos. exp and frac does not matter.
		if (!sign) (regSA) ? (uZ.ui= MAXPOS32): (uZ.ui=0x1);
		else (regSA) ? (uZ.ui= MAXPOS32): (uZ.ui=0x0);
	}
	else{
		//remove hidden bits
		frac64A = (frac64A & 0x3FFFFFFFFFFFFFFF) >>(regA + 2) ; // 2 bits exp

		fracA = frac64A>>32;

		if (regA<=28){
			bitNPlusOne |= (NAR32 & frac64A) ;
			expA <<= (28-regA);
		}
		else {
			if (regA==30){
				bitNPlusOne = expA&0x2;
				bitsMore = (expA&0x1);
				expA = 0;
			}
			else if (regA==29){
				bitNPlusOne = expA&0x1;
				expA>>=1;
			}
			if (fracA>0){
				fracA=0;
				bitsMore =1;
			}
		}

		uZ.ui = packToP32UI(regime, expA, fracA);

		if (0x7FFFFFFF & frac64A) bitsMore=1;
		
		if (bitNPlusOne || bitsMore)
		{
			if (!sign) return p32_nextup(uZ.p);
		}
	}
	if (sign) uZ.ui = -uZ.ui & 0xFFFFFFFF;
	return uZ.p;
}

posit16_t softposit_subMagsP16_up(uint_fast16_t uiA, uint_fast16_t uiB)
{
	uint_fast16_t regA;
	uint_fast32_t frac32A, frac32B;
	uint_fast16_t fracA=0, regime, tmp;
	bool sign=0, regSA, regSB, ecarry=0, bitNPlusOne=0, bitsMore=0;
	int_fast16_t shiftRight;
	int_fast8_t kA=0, expA;
    union ui16_p16 uZ;

    //Both uiA and uiB are actually the same signs if uiB inherits sign of sub
    //Make both positive
    sign = signP16UI( uiA );
    (sign)?(uiA = (-uiA & 0xFFFF)): (uiB = (-uiB & 0xFFFF));

    if (uiA==uiB){ //essential, if not need special handling
		uZ.ui = 0;
		return uZ.p;
	}
    if(uiA<uiB){
		uiA ^= uiB;
		uiB ^= uiA;
		uiA ^= uiB;
		(sign) ? (sign = 0 ) : (sign=1); //A becomes B
	}

    regSA = signregP16UI( uiA );
    regSB = signregP16UI( uiB );

    tmp = (uiA<<2) & 0xFFFF;
	if (regSA){
		while (tmp>>15){
			kA++;
			tmp= (tmp<<1) & 0xFFFF;
		}
	}
	else{
		kA=-1;
		while (!(tmp>>15)){
			kA--;
			tmp= (tmp<<1) & 0xFFFF;
		}
		tmp&=MAXPOS16;
	}
	expA = tmp>>14;
	frac32A = (0x4000 | tmp) << 16;
	shiftRight = kA;

	tmp = (uiB<<2) & 0xFFFF;
	if (regSB){
		while (tmp>>15){
			shiftRight--;
			tmp= (tmp<<1) & 0xFFFF;
		}
	}
	else{
		shiftRight++;
		while (!(tmp>>15)){
			shiftRight++;
			tmp= (tmp<<1) & 0xFFFF;
		}
		tmp&=MAXPOS16;
	}
	frac32B = (0x4000 | tmp) <<16;
	//This is 2kZ + expZ; (where kZ=kA-kB and expZ=expA-expB)

	shiftRight = (shiftRight<<1) + expA - (tmp>>14);

	if (shiftRight!=0){
		if (shiftRight>=29){
			uZ.ui = uiA;
			if (!sign) return uZ.p;
			else 
			{
				uZ.ui = -uZ.ui & 0xFFFF;
				return p16_nextup(uZ.p);
			}
		}
		else
			frac32B >>= shiftRight;
	}

	frac32A -= frac32B;

	while((frac32A>>29)==0){
		kA--;
		frac32A<<=2;
	}
	ecarry = (0x40000000 & frac32A)>>30;
	if(!ecarry){
		if (expA==0) kA--;
		expA^=1;
		frac32A<<=1;
	}

	if(kA<0){
		regA = (-kA & 0xFFFF);
		regSA = 0;
		regime = 0x4000>>regA;
	}
	else{
		regA = kA+1;
		regSA=1;
		regime = MAXPOS16 - (MAXPOS16>>regA);
	}

	if(regA>14){
		//max or min pos. exp and frac does not matter.
		if (!sign) (regSA) ? (uZ.ui= MAXPOS16): (uZ.ui=0x1);
		else (regSA) ? (uZ.ui= MAXPOS16): (uZ.ui=0x0);
	}
	else{
		//remove hidden bits
		frac32A = (frac32A & 0x3FFFFFFF) >>(regA + 1) ;
		fracA = frac32A>>16;
		if (regA!=14) bitNPlusOne = (frac32A>>15) & 0x1;
		else if (frac32A>0){
			fracA=0;
			bitsMore =1;
		}
		if (regA==14 && expA) bitNPlusOne = 1;
		uZ.ui = packToP16UI(regime, regA, expA, fracA);
		if ( frac32A&0x7FFF ) bitsMore=1;
		if (bitNPlusOne || bitsMore)
		{
			if (!sign) return p16_nextup(uZ.p);
		}
	}
	if (sign) uZ.ui = -uZ.ui & 0xFFFF;
	return uZ.p;
}

posit32_t softposit_subMagsP32_up(uint_fast32_t uiA, uint_fast32_t uiB) 
{
	uint_fast16_t regA, regB;
	uint_fast64_t frac64A=0, frac64B=0;
	uint_fast32_t fracA=0, regime, tmp;
	bool sign, regSA, regSB, ecarry=0, bitNPlusOne=0, bitsMore=0;
	int_fast8_t kA=0;
	int_fast32_t expA=0;
	int_fast16_t shiftRight;
	union ui32_p32 uZ;

	sign = signP32UI( uiA );
	if (sign)
		uiA = -uiA & 0xFFFFFFFF;
	else
		uiB = -uiB & 0xFFFFFFFF;

	if (uiA==uiB){ //essential, if not need special handling
		uZ.ui = 0;
		return uZ.p;
	}
	if ((int_fast32_t)uiA < (int_fast32_t)uiB){
		uiA ^= uiB;
		uiB ^= uiA;
		uiA ^= uiB;
		(sign) ? (sign = 0 ) : (sign=1); //A becomes B
	}
	regSA = signregP32UI( uiA );
	regSB = signregP32UI( uiB );

	tmp = (uiA<<2)&0xFFFFFFFF;
	if (regSA){
		while (tmp>>31){
			kA++;
			tmp= (tmp<<1) & 0xFFFFFFFF;
		}
	}
	else{
		kA=-1;
		while (!(tmp>>31)){
			kA--;
			tmp= (tmp<<1) & 0xFFFFFFFF;
		}
		tmp&=MAXPOS32;
	}

	expA = tmp>>29; //to get 2 bits
	frac64A = ((0x40000000ULL | tmp<<1) & 0x7FFFFFFFULL) <<32;
	shiftRight = kA;


	tmp = (uiB<<2) & 0xFFFFFFFF;
	if (regSB){
		while (tmp>>31){
			shiftRight--;
			tmp= (tmp<<1) & 0xFFFFFFFF;
		}

	}
	else{
		shiftRight++;
		while (!(tmp>>31)){
			shiftRight++;
			tmp= (tmp<<1) & 0xFFFFFFFF;
		}
		tmp&=MAXPOS32;

	}
	frac64B = ((0x40000000ULL | tmp<<1) & 0x7FFFFFFFULL) <<32;

	//This is 4kZ + expZ; (where kZ=kA-kB and expZ=expA-expB)
	shiftRight = (shiftRight<<2) + expA - (tmp>>29);
	if (shiftRight>63){
		uZ.ui = uiA;
		if (!sign) return uZ.p;
		else 
		{
			uZ.ui = -uZ.ui & 0xFFFFFFFF;
		 	return p32_nextup(uZ.p);
		}
	}
	else
		(frac64B >>= shiftRight);

	frac64A -= frac64B;

	while((frac64A>>59)==0){
		kA--;
		frac64A<<=4;
	}
	ecarry = (0x4000000000000000 & frac64A);//(0x4000000000000000 & frac64A)>>62;
	while (!ecarry){
		if (expA==0){
			kA--;
			expA=3;
		}
		else
			expA--;
		frac64A<<=1;
		ecarry = (0x4000000000000000 & frac64A);
	}

	if(kA<0){
		regA = -kA;
		regSA = 0;
		regime = 0x40000000>>regA;
	}
	else{
		regA = kA+1;
		regSA=1;
		regime = MAXPOS32 - (MAXPOS32>>regA);
	}
	if(regA>30){
		//max or min pos. exp and frac does not matter.
		if (!sign) (regSA) ? (uZ.ui= MAXPOS32): (uZ.ui=0x1);
		else (regSA) ? (uZ.ui= MAXPOS32): (uZ.ui=0x0);
	}
	else{
		//remove hidden bits
		frac64A = (frac64A & 0x3FFFFFFFFFFFFFFF) >>(regA + 2) ; // 2 bits exp

		fracA = frac64A>>32;

		if (regA<=28){
			bitNPlusOne |= (NAR32 & frac64A) ;
			expA <<= (28-regA);
		}
		else {
			if (regA==30){
				bitNPlusOne = expA&0x2;
				bitsMore = (expA&0x1);
				expA = 0;
			}
			else if (regA==29){
				bitNPlusOne = expA&0x1;
				expA>>=1;
			}
			if (fracA>0){
				fracA=0;
				bitsMore =1;
			}

		}

		uZ.ui = packToP32UI(regime, expA, fracA);
		if (0x7FFFFFFF & frac64A) bitsMore=1;
		if (bitNPlusOne || bitsMore)
		{
			if (!sign) return p32_nextup(uZ.p);
		}
	}
	if (sign) uZ.ui = -uZ.ui & 0xFFFFFFFF;
	return uZ.p;
}


posit16_t p16_add_up(posit16_t a, posit16_t b)
{
	union ui16_p16 uA, uB;
    uint_fast16_t uiA, uiB;
    union ui16_p16 uZ;

    uA.p = a;
	uiA = uA.ui;
	uB.p = b;
	uiB = uB.ui;


    //Zero or infinity
	if (uiA==0 || uiB==0){ // Not required but put here for speed
		uZ.ui = uiA | uiB;
		return uZ.p;
	}
	else if ( uiA==0x8000 || uiB==0x8000 ){
		uZ.ui = 0x8000;
		return uZ.p;
	}

	//different signs
	if ((uiA^uiB)>>15)
		return softposit_subMagsP16_up(uiA, uiB);
	else
		 return softposit_addMagsP16_up(uiA, uiB);
}

posit32_t p32_add_up(posit32_t a, posit32_t b)
{
 	union ui32_p32 uA, uB, uZ;
    uint_fast32_t uiA, uiB;

    uA.p = a;
	uiA = uA.ui;
	uB.p = b;
	uiB = uB.ui;

#ifdef SOFTPOSIT_EXACT
		uZ.ui.exact = (uiA.ui.exact & uiB.ui.exact);
#endif

    //Zero or infinity
	if (uiA==0 || uiB==0){ // Not required but put here for speed
#ifdef SOFTPOSIT_EXACT
		uZ.ui.v = uiA | uiB;
		uZ.ui.exact = (uiA.ui.exact & uiB.ui.exact);
#else
		uZ.ui = uiA | uiB;
#endif
		return uZ.p;
	}
	else if ( uiA==NAR32 || uiB==NAR32 ){
		//printf("in infinity\n");
#ifdef SOFTPOSIT_EXACT
		uZ.ui.v = NAR32;
		uZ.ui.exact = 0;
#else
		uZ.ui = NAR32;
#endif
		return uZ.p;
	}

	//different signs
	if ((uiA^uiB)>>31)
		return softposit_subMagsP32_up(uiA, uiB);
	else
		return softposit_addMagsP32_up(uiA, uiB);
}

posit16_t p16_sub_up(posit16_t a, posit16_t b)
{
	union ui16_p16 uA, uB;
    uint_fast16_t uiA, uiB;
    union ui16_p16 uZ;

    uA.p = a;
	uiA = uA.ui;
	uB.p = b;
	uiB = uB.ui;

#ifdef SOFTPOSIT_EXACT
		uZ.ui.exact = (uiA.ui.exact & uiB.ui.exact);
#endif

    //infinity
	if ( uiA==0x8000 || uiB==0x8000 ){
#ifdef SOFTPOSIT_EXACT
		uZ.ui.v = 0x8000;
		uZ.ui.exact = 0;
#else
		uZ.ui = 0x8000;
#endif
		return uZ.p;
	}
    //Zero
	else if ( uiA==0 || uiB==0 ){
#ifdef SOFTPOSIT_EXACT
		uZ.ui.v = (uiA | -uiB);
		uZ.ui.exact = 0;
#else
		uZ.ui = (uiA | -uiB);
#endif
		return uZ.p;
	}

	//different signs
	if ((uiA^uiB)>>15)
		return softposit_addMagsP16_up(uiA, (-uiB & 0xFFFF));
	else
		return softposit_subMagsP16_up(uiA, (-uiB & 0xFFFF));

}

posit32_t p32_sub_up(posit32_t a, posit32_t b)
{
	union ui32_p32 uA, uB, uZ;
	uint_fast32_t uiA, uiB;

	uA.p = a;
	uiA = uA.ui;
	uB.p = b;
	uiB = uB.ui;

#ifdef SOFTPOSIT_EXACT
		uZ.ui.exact = (uiA.ui.exact & uiB.ui.exact);
#endif

	//infinity
	if ( uiA==NAR32 || uiB==NAR32 ){
#ifdef SOFTPOSIT_EXACT
		uZ.ui.v = NAR32;
		uZ.ui.exact = 0;
#else
		uZ.ui = NAR32;
#endif
		return uZ.p;
	}
	//Zero
	else if ( uiA==0 || uiB==0 ){
#ifdef SOFTPOSIT_EXACT
		uZ.ui.v = (uiA | -uiB);
		uZ.ui.exact = 0;
#else
		uZ.ui = (uiA | -uiB);
#endif
		return uZ.p;
	}

	//different signs
	if ((uiA^uiB)>>31)
		return softposit_addMagsP32_up(uiA, (-uiB & 0xFFFFFFFF));
	else
		return softposit_subMagsP32_up(uiA, (-uiB & 0xFFFFFFFF));
}

posit16_t p16_mul_up(posit16_t a, posit16_t b)
{
	union ui16_p16 uA, uB, uZ;
	uint_fast16_t uiA, uiB;
	uint_fast16_t regA, fracA, regime, tmp;
	bool signA, signB, signZ, regSA, regSB, bitNPlusOne=0, bitsMore=0, rcarry;
	int_fast8_t expA;
	int_fast8_t kA=0;
	uint_fast32_t frac32Z;

	uA.p = a;
	uiA = uA.ui;
	uB.p = b;
	uiB = uB.ui;


	//NaR or Zero
	if ( uiA==0x8000 || uiB==0x8000 ){
		uZ.ui = 0x8000;
		return uZ.p;
	}
	else if (uiA==0 || uiB==0){
		uZ.ui = 0;
		return uZ.p;
	}

	signA = signP16UI( uiA );
	signB = signP16UI( uiB );
	signZ = signA ^ signB;

	if(signA) uiA = (-uiA & 0xFFFF);
	if(signB) uiB = (-uiB & 0xFFFF);

	regSA = signregP16UI(uiA);
	regSB = signregP16UI(uiB);

	tmp = (uiA<<2) & 0xFFFF;
	if (regSA){
		while (tmp>>15){
			kA++;
			tmp= (tmp<<1) & 0xFFFF;
		}
	}
	else{
		kA=-1;
		while (!(tmp>>15)){
			kA--;
			tmp= (tmp<<1) & 0xFFFF;
		}
		tmp&=0x7FFF;
	}
	expA = tmp>>14;
	fracA = (0x4000 | tmp);

	tmp = (uiB<<2) & 0xFFFF;
	if (regSB){
		while (tmp>>15){
			kA++;
			tmp= (tmp<<1) & 0xFFFF;
		}
	}
	else{
		kA--;
		while (!(tmp>>15)){
			kA--;
			tmp= (tmp<<1) & 0xFFFF;
		}
		tmp&=0x7FFF;
	}
	expA += tmp>>14;
	frac32Z = (uint_fast32_t) fracA * (0x4000 | tmp);

	if (expA>1){
		kA++;
		expA ^=0x2;
	}

	rcarry = frac32Z>>29;//3rd bit of frac32Z
	if (rcarry){
		if (expA) kA ++;
		expA^=1;
		frac32Z>>=1;
	}

	if(kA<0){
		regA = (-kA & 0xFFFF);
		regSA = 0;
		regime = 0x4000>>regA;
	}
	else{
		regA = kA+1;
		regSA=1;
		regime = 0x7FFF - (0x7FFF>>regA);
	}

	if(regA>14){
		//max or min pos. exp and frac does not matter.
		if (!signZ) (regSA) ? (uZ.ui= 0x7FFF): (uZ.ui=0x1);
		else (regSA) ? (uZ.ui= 0x7FFF): (uZ.ui=0x0);
	}
	else{
		//remove carry and rcarry bits and shift to correct position
		frac32Z = (frac32Z&0xFFFFFFF) >> (regA-1);
		fracA = (uint_fast16_t) (frac32Z>>16);

		if (regA!=14) bitNPlusOne |= (0x8000 & frac32Z) ;
		else if (fracA>0){
			fracA=0;
			bitsMore =1;
		}
		if (regA==14 && expA) bitNPlusOne = 1;

		//sign is always zero
		uZ.ui = packToP16UI(regime, regA, expA, fracA);

		if (bitNPlusOne || bitsMore)
		{
			if (!signZ) return p16_nextup(uZ.p);
		}
	}

	if (signZ) uZ.ui = -uZ.ui & 0xFFFF;
	return uZ.p;
}

posit32_t p32_mul_up(posit32_t a, posit32_t b)
{
	union ui32_p32 uA, uB, uZ;
	uint_fast32_t uiA, uiB;
	uint_fast32_t regA, fracA, regime, tmp;
	bool signA, signB, signZ, regSA, regSB, bitNPlusOne=0, bitsMore=0, rcarry;
	int_fast32_t expA;
	int_fast8_t kA=0;
	uint_fast64_t frac64Z;

	uA.p = a;
	uiA = uA.ui;
	uB.p = b;
	uiB = uB.ui;

#ifdef SOFTPOSIT_EXACT
		uZ.ui.exact = (uiA.ui.exact & uiB.ui.exact);
#endif
	//NaR or Zero
	if ( uiA==0x80000000 || uiB==0x80000000 ){

#ifdef SOFTPOSIT_EXACT
		uZ.ui.v = 0x80000000;
		uZ.ui.exact = 0;
#else
		uZ.ui = 0x80000000;
#endif
		return uZ.p;
	}
	else if (uiA==0 || uiB==0){
#ifdef SOFTPOSIT_EXACT

		uZ.ui.v = 0;
		if ( (uiA==0 && uiA.ui.exact) || (uiB==0 && uiB.ui.exact) )
			uZ.ui.exact = 1;
		else
			uZ.ui.exact = 0;
#else
		uZ.ui = 0;
#endif
		return uZ.p;
	}

	signA = signP32UI( uiA );
	signB = signP32UI( uiB );
	signZ = signA ^ signB;

	if(signA) uiA = (-uiA & 0xFFFFFFFF);
	if(signB) uiB = (-uiB & 0xFFFFFFFF);

	regSA = signregP32UI(uiA);
	regSB = signregP32UI(uiB);

	tmp = (uiA<<2)&0xFFFFFFFF;
	if (regSA){

		while (tmp>>31){

			kA++;
			tmp= (tmp<<1) & 0xFFFFFFFF;
		}
	}
	else{
		kA=-1;
		while (!(tmp>>31)){
			kA--;
			tmp= (tmp<<1) & 0xFFFFFFFF;
		}
		tmp&=0x7FFFFFFF;
	}
	expA = tmp>>29; //to get 2 bits
	fracA = ((tmp<<1) | 0x40000000) & 0x7FFFFFFF;

	tmp = (uiB<<2)&0xFFFFFFFF;
	if (regSB){
		while (tmp>>31){
			kA++;
			tmp= (tmp<<1) & 0xFFFFFFFF;
		}
	}
	else{
		kA--;
		while (!(tmp>>31)){
			kA--;
			tmp= (tmp<<1) & 0xFFFFFFFF;
		}
		tmp&=0x7FFFFFFF;
	}
	expA += tmp>>29;
	frac64Z = (uint_fast64_t) fracA * (((tmp<<1) | 0x40000000) & 0x7FFFFFFF);

	if (expA>3){
		kA++;
		expA&=0x3; // -=4
	}

	rcarry = frac64Z>>61;//3rd bit of frac64Z
	if (rcarry){
		expA++;
		if (expA>3){
			kA ++;
			expA&=0x3;
		}
		frac64Z>>=1;
	}

	if(kA<0){
		regA = -kA;
		regSA = 0;
		regime = 0x40000000>>regA;
	}
	else{
		regA = kA+1;
		regSA=1;
		regime = 0x7FFFFFFF - (0x7FFFFFFF>>regA);
	}


	if(regA>30){
		//max or min pos. exp and frac does not matter.
		if (!signZ) (regSA) ? (uZ.ui= 0x7FFFFFFF): (uZ.ui=0x1);
		else (regSA) ? (uZ.ui= 0x7FFFFFFF): (uZ.ui=0x0);
	}
	else{
		//remove carry and rcarry bits and shift to correct position (2 bits exp, so + 1 than 16 bits)
		frac64Z = (frac64Z&0xFFFFFFFFFFFFFFF) >> regA;
		fracA = (uint_fast32_t) (frac64Z>>32);
		if (regA<=28){
			bitNPlusOne |= (0x80000000 & frac64Z);
			expA<<= (28-regA);
		}
		else {
			if (regA==30){
				bitNPlusOne = expA&0x2;
				bitsMore = (expA&0x1);
				expA = 0;
			}
			else if (regA==29){
				bitNPlusOne = expA&0x1;
				expA>>=1; //taken care of by the pack algo
			}
			if (fracA>0){
				fracA=0;
				bitsMore =1;
			}

		}
		//sign is always zero
		uZ.ui = packToP32UI(regime, expA, fracA);

		if (bitNPlusOne || bitsMore) 
		{
			if (!signZ) return p32_nextup(uZ.p);	
		}
	}

	if (signZ) uZ.ui = -uZ.ui & 0xFFFFFFFF;
	return uZ.p;
}

posit16_t p16_div_up(posit16_t a, posit16_t b) 
{
	union ui16_p16 uA, uB, uZ;
	uint_fast16_t uiA, uiB, fracA, fracB, regA, regime, tmp;
	bool signA, signB, signZ, regSA, regSB, bitNPlusOne=0, bitsMore=0, rcarry;
	int_fast8_t expA, kA=0;
	uint_fast32_t frac32A, frac32Z, rem;
	div_t divresult;

	uA.p = a;
	uiA = uA.ui;
	uB.p = b;
	uiB = uB.ui;

	//Zero or infinity
	if ( uiA==0x8000 || uiB==0x8000 || uiB==0){
#ifdef SOFTPOSIT_EXACT
		uZ.ui.v = 0x8000;
		uZ.ui.exact = 0;
#else
		uZ.ui = 0x8000;
#endif
		return uZ.p;
	}
	else if (uiA==0){
#ifdef SOFTPOSIT_EXACT
		uZ.ui.v = 0;
		if ( (uiA==0 && uiA.ui.exact) || (uiB==0 && uiB.ui.exact) )
			uZ.ui.exact = 1;
		else
			uZ.ui.exact = 0;
#else
		uZ.ui = 0;
#endif
		return uZ.p;
	}

	signA = signP16UI( uiA );
	signB = signP16UI( uiB );
	signZ = signA ^ signB;
	if(signA) uiA = (-uiA & 0xFFFF);
	if(signB) uiB = (-uiB & 0xFFFF);
	regSA = signregP16UI(uiA);
	regSB = signregP16UI(uiB);

	tmp = (uiA<<2) & 0xFFFF;
	if (regSA){
		while (tmp>>15){
			kA++;
			tmp= (tmp<<1) & 0xFFFF;
		}
	}
	else{
		kA=-1;
		while (!(tmp>>15)){
			kA--;
			tmp= (tmp<<1) & 0xFFFF;
		}
		tmp&=0x7FFF;
	}
	expA = tmp>>14;
	fracA = (0x4000 | tmp);
	frac32A = fracA<<14;

	tmp = (uiB<<2) & 0xFFFF;
	if (regSB){
		while (tmp>>15){
			kA--;
			tmp= (tmp<<1) & 0xFFFF;
		}
		fracB = (0x4000 | tmp);
	}
	else{
		kA++;
		while (!(tmp>>15)){
			kA++;
			tmp= (tmp<<1) & 0xFFFF;
		}
		tmp&=0x7FFF;
		fracB = (0x4000 | (0x7FFF & tmp));
	}
	expA -= tmp>>14;

	divresult = div (frac32A,fracB);
	frac32Z = divresult.quot;
	rem = divresult.rem;

	if (expA<0){
		expA=1;
		kA--;
	}
	if (frac32Z!=0){
		rcarry = frac32Z >> 14; // this is the hidden bit (14th bit) , extreme right bit is bit 0
		if (!rcarry){
			if (expA==0) kA --;
			expA^=1;
			frac32Z<<=1;
		}
	}
	if(kA<0){
		regA = (-kA & 0xFFFF);
		regSA = 0;
		regime = 0x4000>>regA;
	}
	else{
		regA = kA+1;
		regSA=1;
		regime = 0x7FFF - (0x7FFF>>regA);
	}

	if(regA>14){
		//max or min pos. exp and frac does not matter.
		if (!signZ) (regSA) ? (uZ.ui= 0x7FFF): (uZ.ui=0x1);
		else (regSA) ? (uZ.ui= 0x7FFF): (uZ.ui=0x0);
	}
	else{
		//remove carry and rcarry bits and shift to correct position
		frac32Z &= 0x3FFF;
		fracA = (uint_fast16_t)frac32Z >> (regA+1);

		if (regA!=14) bitNPlusOne = (frac32Z >> regA) & 0x1;
		else if (fracA>0){
			fracA=0;
			bitsMore =1;
		}
		if (regA==14 && expA) bitNPlusOne = 1;

		//sign is always zero
		uZ.ui = packToP16UI(regime, regA, expA, fracA);

		if (bitNPlusOne){
			( ((1<<regA)-1) & frac32Z ) ? (bitsMore=1) : (bitsMore=0);
			if (rem) bitsMore =1;
		}

		if (bitNPlusOne || bitsMore)
		{
			if (!signZ) uZ.p = p16_nextup(uZ.p);
		}
	}
	if (signZ) uZ.ui = -uZ.ui & 0xFFFF;

	return uZ.p;
}

posit32_t p32_div_up(posit32_t a, posit32_t b)
{
    union ui32_p32 uA, uB, uZ;
    uint_fast32_t uiA, uiB, fracA, fracB, regA, regime, regB, tmp;
    bool signA, signB, signZ, regSA, regSB, bitNPlusOne=0, bitsMore=0, rcarry;
	int_fast8_t kA=0;
	int_fast32_t expA;
	uint_fast64_t frac64A, frac64Z, rem;
	lldiv_t divresult;

	uA.p = a;
	uiA = uA.ui;
	uB.p = b;
	uiB = uB.ui;

	//Zero or infinity
	if ( uiA==0x80000000 || uiB==0x80000000 || uiB==0){
#ifdef SOFTPOSIT_EXACT
		uZ.ui.v = 0x80000000;
		uZ.ui.exact = 0;
#else
		uZ.ui = 0x80000000;
#endif
		return uZ.p;
	}
	else if (uiA==0){
#ifdef SOFTPOSIT_EXACT

		uZ.ui.v = 0;
		if ( (uiA==0 && uiA.ui.exact) || (uiB==0 && uiB.ui.exact) )
			uZ.ui.exact = 1;
		else
			uZ.ui.exact = 0;
#else
		uZ.ui = 0;
#endif
		return uZ.p;
	}

	signA = signP32UI( uiA );
	signB = signP32UI( uiB );
	signZ = signA ^ signB;
	if(signA) uiA = (-uiA & 0xFFFFFFFF);
	if(signB) uiB = (-uiB & 0xFFFFFFFF);
	regSA = signregP32UI(uiA);
	regSB = signregP32UI(uiB);

	tmp = (uiA<<2)&0xFFFFFFFF;
	if (regSA){
		while (tmp>>31){
			kA++;
			tmp= (tmp<<1) & 0xFFFFFFFF;
		}
	}
	else{
		kA=-1;
		while (!(tmp>>31)){
			kA--;
			tmp= (tmp<<1) & 0xFFFFFFFF;
		}
		tmp&=0x7FFFFFFF;
	}
	expA = tmp>>29; //to get 2 bits
	fracA = ((tmp<<1) | 0x40000000) & 0x7FFFFFFF;
	frac64A = (uint64_t) fracA << 30;

	tmp = (uiB<<2)&0xFFFFFFFF;
	if (regSB){
		while (tmp>>31){
			kA--;
			tmp= (tmp<<1) & 0xFFFFFFFF;
		}
	}
	else{
		kA++;
		while (!(tmp>>31)){
			kA++;
			tmp= (tmp<<1) & 0xFFFFFFFF;
		}
		tmp&=0x7FFFFFFF;
	}
	expA -= tmp>>29;
	fracB = ((tmp<<1) | 0x40000000) & 0x7FFFFFFF;

	divresult = lldiv (frac64A,(uint_fast64_t)fracB);
	frac64Z = divresult.quot;
	rem = divresult.rem;

	if (expA<0){
		expA+=4;
		kA--;
	}
	if (frac64Z!=0){
		rcarry = frac64Z >> 30; // this is the hidden bit (14th bit) , extreme right bit is bit 0
		if (!rcarry){
			if (expA==0){
				kA--;
				expA=3;
			}
			else
				expA--;
			frac64Z<<=1;
		}
	}

	if(kA<0){
		regA = -kA;
		regSA = 0;
		regime = 0x40000000>>regA;
	}
	else{
		regA = kA+1;
		regSA=1;
		regime = 0x7FFFFFFF - (0x7FFFFFFF>>regA);
	}
	if(regA>30){
		//max or min pos. exp and frac does not matter.
		if (!signZ) (regSA) ? (uZ.ui= 0x7FFFFFFF): (uZ.ui=0x1);
		else (regSA) ? (uZ.ui= 0x7FFFFFFF): (uZ.ui=0x0);
	}
	else{
		//remove carry and rcarry bits and shift to correct position
		frac64Z &= 0x3FFFFFFF;

		fracA = (uint_fast32_t)frac64Z >> (regA+2);

		if (regA<=28){
			bitNPlusOne = (frac64Z >> (regA +1)) & 0x1;
			expA<<= (28-regA);
			if (bitNPlusOne) ( ((1<<(regA+1))-1) & frac64Z ) ? (bitsMore=1) : (bitsMore=0);
		}
		else {
			if (regA==30){
				bitNPlusOne = expA&0x2;
				bitsMore = (expA&0x1);
				expA = 0;
			}
			else if (regA==29){
				bitNPlusOne = expA&0x1;
				expA>>=1; //taken care of by the pack algo
			}
			if (frac64Z>0){
				fracA=0;
				bitsMore =1;
			}

		}
		if (rem) bitsMore =1;

		uZ.ui = packToP32UI(regime, expA, fracA);
		if (bitNPlusOne || bitsMore) 
		{
			if (!signZ) uZ.p = p32_nextup(uZ.p);
		}
	}

	if (signZ) uZ.ui = -uZ.ui & 0xFFFFFFFF;
	return uZ.p;
}



// round down
posit16_t softposit_addMagsP16_down(uint_fast16_t uiA, uint_fast16_t uiB)
{
	uint_fast16_t regA, uiX, uiY;
	uint_fast32_t frac32A, frac32B;
	uint_fast16_t fracA=0,  regime, tmp;
	bool sign, regSA, regSB, rcarry=0, bitNPlusOne=0, bitsMore=0;
	int_fast8_t kA=0, expA;
	int_fast16_t shiftRight;
	union ui16_p16 uZ;

	sign = signP16UI( uiA ); //sign is always positive.. actually don't have to do this.
	if (sign){
		uiA = -uiA & 0xFFFF;
		uiB = -uiB & 0xFFFF;
	}

	if ((int_fast16_t)uiA < (int_fast16_t)uiB){
		uiX = uiA;
		uiY = uiB;
		uiA = uiY;
		uiB = uiX;
	}
	regSA = signregP16UI( uiA );
	regSB = signregP16UI( uiB );

	tmp = (uiA<<2) & 0xFFFF;
	if (regSA){
		while (tmp>>15){
			kA++;
			tmp= (tmp<<1) & 0xFFFF;
		}
	}
	else{
		kA=-1;
		while (!(tmp>>15)){
			kA--;
			tmp= (tmp<<1) & 0xFFFF;
		}
		tmp&=MAXPOS16;
	}
	expA = tmp>>14;
	frac32A = (0x4000 | tmp) << 16;
	shiftRight = kA;

	tmp = (uiB<<2) & 0xFFFF;
	if (regSB){
		while (tmp>>15){
			shiftRight--;
			tmp= (tmp<<1) & 0xFFFF;
		}
		frac32B = (0x4000 | tmp) <<16;
	}
	else{
		shiftRight++;
		while (!(tmp>>15)){
			shiftRight++;
			tmp= (tmp<<1) & 0xFFFF;
		}
		tmp&=MAXPOS16;
		frac32B = ( (0x4000 | tmp) <<16 ) & MAXPOS32;
	}

	//This is 2kZ + expZ; (where kZ=kA-kB and expZ=expA-expB)
	shiftRight = (shiftRight<<1) + expA - (tmp>>14);

	if (shiftRight==0){
		frac32A += frac32B;
		//rcarry is one
		if (expA) kA ++;
		expA^=1;
		frac32A>>=1;
	}
	else{
		//Manage CLANG (LLVM) compiler when shifting right more than number of bits
		(shiftRight>31) ? (frac32B=0): (frac32B >>= shiftRight); //frac32B >>= shiftRight

		frac32A += frac32B;
		rcarry = NAR32 & frac32A; //first left bit
		if(rcarry){
			if (expA) kA ++;
			expA^=1;
			frac32A>>=1;
		}
	}
	if(kA<0){
		regA = (-kA & 0xFFFF);
		regSA = 0;
		regime = 0x4000>>regA;
	}
	else{
		regA = kA+1;
		regSA=1;
		regime = MAXPOS16 - (MAXPOS16>>regA);
	}
	if(regA>14){
		//max or min pos. exp and frac does not matter.
		if (!sign) (regSA) ? (uZ.ui= MAXPOS16): (uZ.ui=0x0);
		else (regSA) ? (uZ.ui= MINNEG16): (uZ.ui=0x1);
	}
	else{
		//remove hidden bits
		frac32A = (frac32A & 0x3FFFFFFF) >>(regA + 1) ;
		fracA = frac32A>>16;
		if (regA!=14) bitNPlusOne = (frac32A>>15) & 0x1;
		else if (frac32A>0){
			fracA=0;
			bitsMore =1;
		}
		if (regA==14 && expA) bitNPlusOne = 1;
		uZ.ui = packToP16UI(regime, regA, expA, fracA);

		if ( frac32A&0x7FFF ) bitsMore=1;

		if (bitNPlusOne || bitsMore)
		{
			if (sign) uZ.p = p16_nextup(uZ.p);
		}
	}

	if (sign) uZ.ui = -uZ.ui & 0xFFFF;
	return uZ.p;
}

posit32_t softposit_addMagsP32_down(uint_fast32_t uiA, uint_fast32_t uiB) 
{
	uint_fast16_t regA, regB;
	uint_fast64_t frac64A=0, frac64B=0;
	uint_fast32_t fracA=0, regime, tmp;
	bool sign, regSA, regSB, rcarry=0, bitNPlusOne=0, bitsMore=0;
	int_fast8_t kA=0;
	int_fast32_t expA;
	int_fast16_t shiftRight;
	union ui32_p32 uZ;

	sign = signP32UI( uiA );
	if (sign){
		uiA = -uiA & 0xFFFFFFFF;
		uiB = -uiB & 0xFFFFFFFF;
	}

	if ((int_fast32_t)uiA < (int_fast32_t)uiB){
		uiA ^= uiB;
		uiB ^= uiA;
		uiA ^= uiB;
	}
	regSA = signregP32UI( uiA );
    regSB = signregP32UI( uiB );

    tmp = (uiA<<2)&0xFFFFFFFF;
	if (regSA){
		while (tmp>>31){
			kA++;
			tmp= (tmp<<1) & 0xFFFFFFFF;
		}
	}
	else{
		kA=-1;
		while (!(tmp>>31)){
			kA--;
			tmp= (tmp<<1) & 0xFFFFFFFF;
		}
		tmp&=MAXPOS32;
	}

	expA = tmp>>29; //to get 2 bits
	frac64A = ((0x40000000ULL | tmp<<1) & 0x7FFFFFFFULL) <<32;
	shiftRight = kA;

	tmp = (uiB<<2) & 0xFFFFFFFF;
	if (regSB){
		while (tmp>>31){
			shiftRight--;
			tmp= (tmp<<1) & 0xFFFFFFFF;
		}
	}
	else{
		shiftRight++;
		while (!(tmp>>31)){
			shiftRight++;
			tmp= (tmp<<1) & 0xFFFFFFFF;
		}
		tmp&=MAXPOS32;
	}
	frac64B = ((0x40000000ULL | tmp<<1) & 0x7FFFFFFFULL) <<32;
	//This is 4kZ + expZ; (where kZ=kA-kB and expZ=expA-expB)
	shiftRight = (shiftRight<<2) + expA - (tmp>>29);

	//Manage CLANG (LLVM) compiler when shifting right more than number of bits
	(shiftRight>63) ? (frac64B=0): (frac64B >>= shiftRight); //frac64B >>= shiftRight

	frac64A += frac64B;

	rcarry = NAR64 & frac64A; //first left bit
	if (rcarry){
		expA++;
		if (expA>3){
			kA ++;
			expA&=0x3;
		}
		frac64A>>=1;
	}
	if(kA<0){
		regA = -kA;
		regSA = 0;
		regime = 0x40000000>>regA;
	}
	else{
		regA = kA+1;
		regSA=1;
		regime = MAXPOS32 - (MAXPOS32>>regA);
	}

	if(regA>30){
		//max or min pos. exp and frac does not matter.
		if (sign) (regSA) ? (uZ.ui= MAXPOS32): (uZ.ui=0x1);
		else (regSA) ? (uZ.ui= MAXPOS32): (uZ.ui=0x0);
	}
	else{
		//remove hidden bits
		frac64A = (frac64A & 0x3FFFFFFFFFFFFFFF) >>(regA + 2) ; // 2 bits exp

		fracA = frac64A>>32;

		if (regA<=28){
			bitNPlusOne |= (NAR32 & frac64A) ;
			expA <<= (28-regA);
		}
		else {
			if (regA==30){
				bitNPlusOne = expA&0x2;
				bitsMore = (expA&0x1);
				expA = 0;
			}
			else if (regA==29){
				bitNPlusOne = expA&0x1;
				expA>>=1;
			}
			if (fracA>0){
				fracA=0;
				bitsMore =1;
			}
		}

		uZ.ui = packToP32UI(regime, expA, fracA);

		if (0x7FFFFFFF & frac64A) bitsMore=1;
		
		if (bitNPlusOne || bitsMore)
		{
			if (sign) uZ.p = p32_nextup(uZ.p);
		}
	}
	if (sign) uZ.ui = -uZ.ui & 0xFFFFFFFF;
	return uZ.p;
}

posit16_t softposit_subMagsP16_down(uint_fast16_t uiA, uint_fast16_t uiB)
{
	uint_fast16_t regA;
	uint_fast32_t frac32A, frac32B;
	uint_fast16_t fracA=0, regime, tmp;
	bool sign=0, regSA, regSB, ecarry=0, bitNPlusOne=0, bitsMore=0;
	int_fast16_t shiftRight;
	int_fast8_t kA=0, expA;
    union ui16_p16 uZ;

    //Both uiA and uiB are actually the same signs if uiB inherits sign of sub
    //Make both positive
    sign = signP16UI( uiA );
    (sign)?(uiA = (-uiA & 0xFFFF)): (uiB = (-uiB & 0xFFFF));

    if (uiA==uiB){ //essential, if not need special handling
		uZ.ui = 0;
		return uZ.p;
	}
    if(uiA<uiB){
		uiA ^= uiB;
		uiB ^= uiA;
		uiA ^= uiB;
		(sign) ? (sign = 0 ) : (sign=1); //A becomes B
	}

    regSA = signregP16UI( uiA );
    regSB = signregP16UI( uiB );

    tmp = (uiA<<2) & 0xFFFF;
	if (regSA){
		while (tmp>>15){
			kA++;
			tmp= (tmp<<1) & 0xFFFF;
		}
	}
	else{
		kA=-1;
		while (!(tmp>>15)){
			kA--;
			tmp= (tmp<<1) & 0xFFFF;
		}
		tmp&=MAXPOS16;
	}
	expA = tmp>>14;
	frac32A = (0x4000 | tmp) << 16;
	shiftRight = kA;

	tmp = (uiB<<2) & 0xFFFF;
	if (regSB){
		while (tmp>>15){
			shiftRight--;
			tmp= (tmp<<1) & 0xFFFF;
		}
	}
	else{
		shiftRight++;
		while (!(tmp>>15)){
			shiftRight++;
			tmp= (tmp<<1) & 0xFFFF;
		}
		tmp&=MAXPOS16;
	}
	frac32B = (0x4000 | tmp) <<16;
	//This is 2kZ + expZ; (where kZ=kA-kB and expZ=expA-expB)

	shiftRight = (shiftRight<<1) + expA - (tmp>>14);

	if (shiftRight!=0){
		if (shiftRight>=29){
			uZ.ui = uiA;
			if (!sign) return p16_nextdown(uZ.p);
			else 
			{
				uZ.ui = -uZ.ui & 0xFFFF;
				return uZ.p;
			}
		}
		else
			frac32B >>= shiftRight;
	}

	frac32A -= frac32B;

	while((frac32A>>29)==0){
		kA--;
		frac32A<<=2;
	}
	ecarry = (0x40000000 & frac32A)>>30;
	if(!ecarry){
		if (expA==0) kA--;
		expA^=1;
		frac32A<<=1;
	}

	if(kA<0){
		regA = (-kA & 0xFFFF);
		regSA = 0;
		regime = 0x4000>>regA;
	}
	else{
		regA = kA+1;
		regSA=1;
		regime = MAXPOS16 - (MAXPOS16>>regA);
	}

	if(regA>14){
		//max or min pos. exp and frac does not matter.
		if (sign) (regSA) ? (uZ.ui= MAXPOS16): (uZ.ui=0x1);
		else (regSA) ? (uZ.ui= MAXPOS16): (uZ.ui=0x0);
	}
	else{
		//remove hidden bits
		frac32A = (frac32A & 0x3FFFFFFF) >>(regA + 1) ;
		fracA = frac32A>>16;
		if (regA!=14) bitNPlusOne = (frac32A>>15) & 0x1;
		else if (frac32A>0){
			fracA=0;
			bitsMore =1;
		}
		if (regA==14 && expA) bitNPlusOne = 1;
		uZ.ui = packToP16UI(regime, regA, expA, fracA);

		if ( frac32A&0x7FFF ) bitsMore=1;

		if (bitNPlusOne || bitsMore)
		{
			if (sign) uZ.p = p16_nextup(uZ.p);
		}
	}
	if (sign) uZ.ui = -uZ.ui & 0xFFFF;
	return uZ.p;
}

posit32_t softposit_subMagsP32_down(uint_fast32_t uiA, uint_fast32_t uiB) 
{
	uint_fast16_t regA, regB;
	uint_fast64_t frac64A=0, frac64B=0;
	uint_fast32_t fracA=0, regime, tmp;
	bool sign, regSA, regSB, ecarry=0, bitNPlusOne=0, bitsMore=0;
	int_fast8_t kA=0;
	int_fast32_t expA=0;
	int_fast16_t shiftRight;
	union ui32_p32 uZ;

	sign = signP32UI( uiA );
	if (sign)
		uiA = -uiA & 0xFFFFFFFF;
	else
		uiB = -uiB & 0xFFFFFFFF;

	if (uiA==uiB){ //essential, if not need special handling
		uZ.ui = 0;
		return uZ.p;
	}
	if ((int_fast32_t)uiA < (int_fast32_t)uiB){
		uiA ^= uiB;
		uiB ^= uiA;
		uiA ^= uiB;
		(sign) ? (sign = 0 ) : (sign=1); //A becomes B
	}
	regSA = signregP32UI( uiA );
	regSB = signregP32UI( uiB );

	tmp = (uiA<<2)&0xFFFFFFFF;
	if (regSA){
		while (tmp>>31){
			kA++;
			tmp= (tmp<<1) & 0xFFFFFFFF;
		}
	}
	else{
		kA=-1;
		while (!(tmp>>31)){
			kA--;
			tmp= (tmp<<1) & 0xFFFFFFFF;
		}
		tmp&=MAXPOS32;
	}

	expA = tmp>>29; //to get 2 bits
	frac64A = ((0x40000000ULL | tmp<<1) & 0x7FFFFFFFULL) <<32;
	shiftRight = kA;


	tmp = (uiB<<2) & 0xFFFFFFFF;
	if (regSB){
		while (tmp>>31){
			shiftRight--;
			tmp= (tmp<<1) & 0xFFFFFFFF;
		}

	}
	else{
		shiftRight++;
		while (!(tmp>>31)){
			shiftRight++;
			tmp= (tmp<<1) & 0xFFFFFFFF;
		}
		tmp&=MAXPOS32;

	}
	frac64B = ((0x40000000ULL | tmp<<1) & 0x7FFFFFFFULL) <<32;

	//This is 4kZ + expZ; (where kZ=kA-kB and expZ=expA-expB)
	shiftRight = (shiftRight<<2) + expA - (tmp>>29);
	if (shiftRight>63){
		uZ.ui = uiA;
		if (!sign) return p32_nextdown(uZ.p);
		else 
		{
			uZ.ui = -uZ.ui & 0xFFFFFFFF;
		 	return uZ.p;
		}
	}
	else
		(frac64B >>= shiftRight);

	frac64A -= frac64B;

	while((frac64A>>59)==0){
		kA--;
		frac64A<<=4;
	}
	ecarry = (0x4000000000000000 & frac64A);//(0x4000000000000000 & frac64A)>>62;
	while (!ecarry){
		if (expA==0){
			kA--;
			expA=3;
		}
		else
			expA--;
		frac64A<<=1;
		ecarry = (0x4000000000000000 & frac64A);
	}

	if(kA<0){
		regA = -kA;
		regSA = 0;
		regime = 0x40000000>>regA;
	}
	else{
		regA = kA+1;
		regSA=1;
		regime = MAXPOS32 - (MAXPOS32>>regA);
	}
	if(regA>30){
		//max or min pos. exp and frac does not matter.
		if (!sign) (regSA) ? (uZ.ui= MAXPOS32): (uZ.ui=0x1);
		else (regSA) ? (uZ.ui= MAXPOS32): (uZ.ui=0x0);
	}
	else{
		//remove hidden bits
		frac64A = (frac64A & 0x3FFFFFFFFFFFFFFF) >>(regA + 2) ; // 2 bits exp

		fracA = frac64A>>32;

		if (regA<=28){
			bitNPlusOne |= (NAR32 & frac64A) ;
			expA <<= (28-regA);
		}
		else {
			if (regA==30){
				bitNPlusOne = expA&0x2;
				bitsMore = (expA&0x1);
				expA = 0;
			}
			else if (regA==29){
				bitNPlusOne = expA&0x1;
				expA>>=1;
			}
			if (fracA>0){
				fracA=0;
				bitsMore =1;
			}

		}

		uZ.ui = packToP32UI(regime, expA, fracA);
		if (0x7FFFFFFF & frac64A) bitsMore=1;
		if (bitNPlusOne || bitsMore)
		{
			if (sign) uZ.p = p32_nextup(uZ.p);
		}
	}
	if (sign) uZ.ui = -uZ.ui & 0xFFFFFFFF;
	return uZ.p;
}


posit16_t p16_add_down(posit16_t a, posit16_t b)
{
	union ui16_p16 uA, uB;
    uint_fast16_t uiA, uiB;
    union ui16_p16 uZ;

    uA.p = a;
	uiA = uA.ui;
	uB.p = b;
	uiB = uB.ui;


    //Zero or infinity
	if (uiA==0 || uiB==0){ // Not required but put here for speed
		uZ.ui = uiA | uiB;
		return uZ.p;
	}
	else if ( uiA==0x8000 || uiB==0x8000 ){
		uZ.ui = 0x8000;
		return uZ.p;
	}

	//different signs
	if ((uiA^uiB)>>15)
		return softposit_subMagsP16_down(uiA, uiB);
	else
		 return softposit_addMagsP16_down(uiA, uiB);
}

posit32_t p32_add_down(posit32_t a, posit32_t b)
{
 	union ui32_p32 uA, uB, uZ;
    uint_fast32_t uiA, uiB;

    uA.p = a;
	uiA = uA.ui;
	uB.p = b;
	uiB = uB.ui;

#ifdef SOFTPOSIT_EXACT
		uZ.ui.exact = (uiA.ui.exact & uiB.ui.exact);
#endif

    //Zero or infinity
	if (uiA==0 || uiB==0){ // Not required but put here for speed
#ifdef SOFTPOSIT_EXACT
		uZ.ui.v = uiA | uiB;
		uZ.ui.exact = (uiA.ui.exact & uiB.ui.exact);
#else
		uZ.ui = uiA | uiB;
#endif
		return uZ.p;
	}
	else if ( uiA==NAR32 || uiB==NAR32 ){
		//printf("in infinity\n");
#ifdef SOFTPOSIT_EXACT
		uZ.ui.v = NAR32;
		uZ.ui.exact = 0;
#else
		uZ.ui = NAR32;
#endif
		return uZ.p;
	}

	//different signs
	if ((uiA^uiB)>>31)
		return softposit_subMagsP32_down(uiA, uiB);
	else
		return softposit_addMagsP32_down(uiA, uiB);
}

posit16_t p16_sub_down(posit16_t a, posit16_t b)
{
	union ui16_p16 uA, uB;
    uint_fast16_t uiA, uiB;
    union ui16_p16 uZ;

    uA.p = a;
	uiA = uA.ui;
	uB.p = b;
	uiB = uB.ui;

#ifdef SOFTPOSIT_EXACT
		uZ.ui.exact = (uiA.ui.exact & uiB.ui.exact);
#endif

    //infinity
	if ( uiA==0x8000 || uiB==0x8000 ){
#ifdef SOFTPOSIT_EXACT
		uZ.ui.v = 0x8000;
		uZ.ui.exact = 0;
#else
		uZ.ui = 0x8000;
#endif
		return uZ.p;
	}
    //Zero
	else if ( uiA==0 || uiB==0 ){
#ifdef SOFTPOSIT_EXACT
		uZ.ui.v = (uiA | -uiB);
		uZ.ui.exact = 0;
#else
		uZ.ui = (uiA | -uiB);
#endif
		return uZ.p;
	}

	//different signs
	if ((uiA^uiB)>>15)
		return softposit_addMagsP16_down(uiA, (-uiB & 0xFFFF));
	else
		return softposit_subMagsP16_down(uiA, (-uiB & 0xFFFF));

}

posit32_t p32_sub_down(posit32_t a, posit32_t b)
{
	union ui32_p32 uA, uB, uZ;
	uint_fast32_t uiA, uiB;

	uA.p = a;
	uiA = uA.ui;
	uB.p = b;
	uiB = uB.ui;

#ifdef SOFTPOSIT_EXACT
		uZ.ui.exact = (uiA.ui.exact & uiB.ui.exact);
#endif

	//infinity
	if ( uiA==NAR32 || uiB==NAR32 ){
#ifdef SOFTPOSIT_EXACT
		uZ.ui.v = NAR32;
		uZ.ui.exact = 0;
#else
		uZ.ui = NAR32;
#endif
		return uZ.p;
	}
	//Zero
	else if ( uiA==0 || uiB==0 ){
#ifdef SOFTPOSIT_EXACT
		uZ.ui.v = (uiA | -uiB);
		uZ.ui.exact = 0;
#else
		uZ.ui = (uiA | -uiB);
#endif
		return uZ.p;
	}

	//different signs
	if ((uiA^uiB)>>31)
		return softposit_addMagsP32_down(uiA, (-uiB & 0xFFFFFFFF));
	else
		return softposit_subMagsP32_down(uiA, (-uiB & 0xFFFFFFFF));
}

posit16_t p16_mul_down(posit16_t a, posit16_t b)
{
	union ui16_p16 uA, uB, uZ;
	uint_fast16_t uiA, uiB;
	uint_fast16_t regA, fracA, regime, tmp;
	bool signA, signB, signZ, regSA, regSB, bitNPlusOne=0, bitsMore=0, rcarry;
	int_fast8_t expA;
	int_fast8_t kA=0;
	uint_fast32_t frac32Z;

	uA.p = a;
	uiA = uA.ui;
	uB.p = b;
	uiB = uB.ui;


	//NaR or Zero
	if ( uiA==0x8000 || uiB==0x8000 ){
		uZ.ui = 0x8000;
		return uZ.p;
	}
	else if (uiA==0 || uiB==0){
		uZ.ui = 0;
		return uZ.p;
	}

	signA = signP16UI( uiA );
	signB = signP16UI( uiB );
	signZ = signA ^ signB;

	if(signA) uiA = (-uiA & 0xFFFF);
	if(signB) uiB = (-uiB & 0xFFFF);

	regSA = signregP16UI(uiA);
	regSB = signregP16UI(uiB);

	tmp = (uiA<<2) & 0xFFFF;
	if (regSA){
		while (tmp>>15){
			kA++;
			tmp= (tmp<<1) & 0xFFFF;
		}
	}
	else{
		kA=-1;
		while (!(tmp>>15)){
			kA--;
			tmp= (tmp<<1) & 0xFFFF;
		}
		tmp&=0x7FFF;
	}
	expA = tmp>>14;
	fracA = (0x4000 | tmp);

	tmp = (uiB<<2) & 0xFFFF;
	if (regSB){
		while (tmp>>15){
			kA++;
			tmp= (tmp<<1) & 0xFFFF;
		}
	}
	else{
		kA--;
		while (!(tmp>>15)){
			kA--;
			tmp= (tmp<<1) & 0xFFFF;
		}
		tmp&=0x7FFF;
	}
	expA += tmp>>14;
	frac32Z = (uint_fast32_t) fracA * (0x4000 | tmp);

	if (expA>1){
		kA++;
		expA ^=0x2;
	}

	rcarry = frac32Z>>29;//3rd bit of frac32Z
	if (rcarry){
		if (expA) kA ++;
		expA^=1;
		frac32Z>>=1;
	}

	if(kA<0){
		regA = (-kA & 0xFFFF);
		regSA = 0;
		regime = 0x4000>>regA;
	}
	else{
		regA = kA+1;
		regSA=1;
		regime = 0x7FFF - (0x7FFF>>regA);
	}

	if(regA>14){
		//max or min pos. exp and frac does not matter.
		if (signZ) (regSA) ? (uZ.ui= 0x7FFF): (uZ.ui=0x1);
		else (regSA) ? (uZ.ui= 0x7FFF): (uZ.ui=0x0);
	}
	else{
		//remove carry and rcarry bits and shift to correct position
		frac32Z = (frac32Z&0xFFFFFFF) >> (regA-1);
		fracA = (uint_fast16_t) (frac32Z>>16);

		if (regA!=14) bitNPlusOne |= (0x8000 & frac32Z) ;
		else if (fracA>0){
			fracA=0;
			bitsMore =1;
		}
		if (regA==14 && expA) bitNPlusOne = 1;

		//sign is always zero
		uZ.ui = packToP16UI(regime, regA, expA, fracA);
		//n+1 frac bit is 1. Need to check if another bit is 1 too if not round to even
		if (bitNPlusOne || bitsMore)
		{
			if (signZ) uZ.p = p16_nextup(uZ.p);
		}
	}

	if (signZ) uZ.ui = -uZ.ui & 0xFFFF;
	return uZ.p;
}

posit32_t p32_mul_down(posit32_t a, posit32_t b)
{
	union ui32_p32 uA, uB, uZ;
	uint_fast32_t uiA, uiB;
	uint_fast32_t regA, fracA, regime, tmp;
	bool signA, signB, signZ, regSA, regSB, bitNPlusOne=0, bitsMore=0, rcarry;
	int_fast32_t expA;
	int_fast8_t kA=0;
	uint_fast64_t frac64Z;

	uA.p = a;
	uiA = uA.ui;
	uB.p = b;
	uiB = uB.ui;

#ifdef SOFTPOSIT_EXACT
		uZ.ui.exact = (uiA.ui.exact & uiB.ui.exact);
#endif
	//NaR or Zero
	if ( uiA==0x80000000 || uiB==0x80000000 ){

#ifdef SOFTPOSIT_EXACT
		uZ.ui.v = 0x80000000;
		uZ.ui.exact = 0;
#else
		uZ.ui = 0x80000000;
#endif
		return uZ.p;
	}
	else if (uiA==0 || uiB==0){
#ifdef SOFTPOSIT_EXACT

		uZ.ui.v = 0;
		if ( (uiA==0 && uiA.ui.exact) || (uiB==0 && uiB.ui.exact) )
			uZ.ui.exact = 1;
		else
			uZ.ui.exact = 0;
#else
		uZ.ui = 0;
#endif
		return uZ.p;
	}

	signA = signP32UI( uiA );
	signB = signP32UI( uiB );
	signZ = signA ^ signB;

	if(signA) uiA = (-uiA & 0xFFFFFFFF);
	if(signB) uiB = (-uiB & 0xFFFFFFFF);

	regSA = signregP32UI(uiA);
	regSB = signregP32UI(uiB);

	tmp = (uiA<<2)&0xFFFFFFFF;
	if (regSA){

		while (tmp>>31){

			kA++;
			tmp= (tmp<<1) & 0xFFFFFFFF;
		}
	}
	else{
		kA=-1;
		while (!(tmp>>31)){
			kA--;
			tmp= (tmp<<1) & 0xFFFFFFFF;
		}
		tmp&=0x7FFFFFFF;
	}
	expA = tmp>>29; //to get 2 bits
	fracA = ((tmp<<1) | 0x40000000) & 0x7FFFFFFF;

	tmp = (uiB<<2)&0xFFFFFFFF;
	if (regSB){
		while (tmp>>31){
			kA++;
			tmp= (tmp<<1) & 0xFFFFFFFF;
		}
	}
	else{
		kA--;
		while (!(tmp>>31)){
			kA--;
			tmp= (tmp<<1) & 0xFFFFFFFF;
		}
		tmp&=0x7FFFFFFF;
	}
	expA += tmp>>29;
	frac64Z = (uint_fast64_t) fracA * (((tmp<<1) | 0x40000000) & 0x7FFFFFFF);

	if (expA>3){
		kA++;
		expA&=0x3; // -=4
	}

	rcarry = frac64Z>>61;//3rd bit of frac64Z
	if (rcarry){
		expA++;
		if (expA>3){
			kA ++;
			expA&=0x3;
		}
		frac64Z>>=1;
	}

	if(kA<0){
		regA = -kA;
		regSA = 0;
		regime = 0x40000000>>regA;
	}
	else{
		regA = kA+1;
		regSA=1;
		regime = 0x7FFFFFFF - (0x7FFFFFFF>>regA);
	}


	if(regA>30){
		//max or min pos. exp and frac does not matter.
		if (!signZ) (regSA) ? (uZ.ui= 0x7FFFFFFF): (uZ.ui=0x1);
		else (regSA) ? (uZ.ui= 0x7FFFFFFF): (uZ.ui=0x0);
	}
	else{
		//remove carry and rcarry bits and shift to correct position (2 bits exp, so + 1 than 16 bits)
		frac64Z = (frac64Z&0xFFFFFFFFFFFFFFF) >> regA;
		fracA = (uint_fast32_t) (frac64Z>>32);
		if (regA<=28){
			bitNPlusOne |= (0x80000000 & frac64Z);
			expA<<= (28-regA);
		}
		else {
			if (regA==30){
				bitNPlusOne = expA&0x2;
				bitsMore = (expA&0x1);
				expA = 0;
			}
			else if (regA==29){
				bitNPlusOne = expA&0x1;
				expA>>=1; //taken care of by the pack algo
			}
			if (fracA>0){
				fracA=0;
				bitsMore =1;
			}

		}
		//sign is always zero
		uZ.ui = packToP32UI(regime, expA, fracA);
		//n+1 frac bit is 1. Need to check if another bit is 1 too if not round to even
		if (bitNPlusOne || bitsMore)
		{
			if (signZ) uZ.p = p32_nextup(uZ.p);
		}
	}

	if (signZ) uZ.ui = -uZ.ui & 0xFFFFFFFF;
	return uZ.p;
}


posit16_t p16_div_down(posit16_t a, posit16_t b) 
{
	union ui16_p16 uA, uB, uZ;
	uint_fast16_t uiA, uiB, fracA, fracB, regA, regime, tmp;
	bool signA, signB, signZ, regSA, regSB, bitNPlusOne=0, bitsMore=0, rcarry;
	int_fast8_t expA, kA=0;
	uint_fast32_t frac32A, frac32Z, rem;
	div_t divresult;

	uA.p = a;
	uiA = uA.ui;
	uB.p = b;
	uiB = uB.ui;

	//Zero or infinity
	if ( uiA==0x8000 || uiB==0x8000 || uiB==0){
#ifdef SOFTPOSIT_EXACT
		uZ.ui.v = 0x8000;
		uZ.ui.exact = 0;
#else
		uZ.ui = 0x8000;
#endif
		return uZ.p;
	}
	else if (uiA==0){
#ifdef SOFTPOSIT_EXACT
		uZ.ui.v = 0;
		if ( (uiA==0 && uiA.ui.exact) || (uiB==0 && uiB.ui.exact) )
			uZ.ui.exact = 1;
		else
			uZ.ui.exact = 0;
#else
		uZ.ui = 0;
#endif
		return uZ.p;
	}

	signA = signP16UI( uiA );
	signB = signP16UI( uiB );
	signZ = signA ^ signB;
	if(signA) uiA = (-uiA & 0xFFFF);
	if(signB) uiB = (-uiB & 0xFFFF);
	regSA = signregP16UI(uiA);
	regSB = signregP16UI(uiB);

	tmp = (uiA<<2) & 0xFFFF;
	if (regSA){
		while (tmp>>15){
			kA++;
			tmp= (tmp<<1) & 0xFFFF;
		}
	}
	else{
		kA=-1;
		while (!(tmp>>15)){
			kA--;
			tmp= (tmp<<1) & 0xFFFF;
		}
		tmp&=0x7FFF;
	}
	expA = tmp>>14;
	fracA = (0x4000 | tmp);
	frac32A = fracA<<14;

	tmp = (uiB<<2) & 0xFFFF;
	if (regSB){
		while (tmp>>15){
			kA--;
			tmp= (tmp<<1) & 0xFFFF;
		}
		fracB = (0x4000 | tmp);
	}
	else{
		kA++;
		while (!(tmp>>15)){
			kA++;
			tmp= (tmp<<1) & 0xFFFF;
		}
		tmp&=0x7FFF;
		fracB = (0x4000 | (0x7FFF & tmp));
	}
	expA -= tmp>>14;

	divresult = div (frac32A,fracB);
	frac32Z = divresult.quot;
	rem = divresult.rem;

	if (expA<0){
		expA=1;
		kA--;
	}
	if (frac32Z!=0){
		rcarry = frac32Z >> 14; // this is the hidden bit (14th bit) , extreme right bit is bit 0
		if (!rcarry){
			if (expA==0) kA --;
			expA^=1;
			frac32Z<<=1;
		}
	}
	if(kA<0){
		regA = (-kA & 0xFFFF);
		regSA = 0;
		regime = 0x4000>>regA;
	}
	else{
		regA = kA+1;
		regSA=1;
		regime = 0x7FFF - (0x7FFF>>regA);
	}

	if(regA>14){
		//max or min pos. exp and frac does not matter.
		if (signZ) (regSA) ? (uZ.ui= 0x7FFF): (uZ.ui=0x1);
		else (regSA) ? (uZ.ui= 0x7FFF): (uZ.ui=0x0);
	}
	else{
		//remove carry and rcarry bits and shift to correct position
		frac32Z &= 0x3FFF;
		fracA = (uint_fast16_t)frac32Z >> (regA+1);

		if (regA!=14) bitNPlusOne = (frac32Z >> regA) & 0x1;
		else if (fracA>0){
			fracA=0;
			bitsMore =1;
		}
		if (regA==14 && expA) bitNPlusOne = 1;

		//sign is always zero
		uZ.ui = packToP16UI(regime, regA, expA, fracA);

		if (bitNPlusOne){
			( ((1<<regA)-1) & frac32Z ) ? (bitsMore=1) : (bitsMore=0);
			if (rem) bitsMore =1;
		}

		if (bitNPlusOne || bitsMore)
		{
			if (signZ) uZ.p = p16_nextup(uZ.p);
		}
	}
	if (signZ) uZ.ui = -uZ.ui & 0xFFFF;

	return uZ.p;
}


posit32_t p32_div_down(posit32_t a, posit32_t b)
{
    union ui32_p32 uA, uB, uZ;
    uint_fast32_t uiA, uiB, fracA, fracB, regA, regime, regB, tmp;
    bool signA, signB, signZ, regSA, regSB, bitNPlusOne=0, bitsMore=0, rcarry;
	int_fast8_t kA=0;
	int_fast32_t expA;
	uint_fast64_t frac64A, frac64Z, rem;
	lldiv_t divresult;

	uA.p = a;
	uiA = uA.ui;
	uB.p = b;
	uiB = uB.ui;

	//Zero or infinity
	if ( uiA==0x80000000 || uiB==0x80000000 || uiB==0){
#ifdef SOFTPOSIT_EXACT
		uZ.ui.v = 0x80000000;
		uZ.ui.exact = 0;
#else
		uZ.ui = 0x80000000;
#endif
		return uZ.p;
	}
	else if (uiA==0){
#ifdef SOFTPOSIT_EXACT

		uZ.ui.v = 0;
		if ( (uiA==0 && uiA.ui.exact) || (uiB==0 && uiB.ui.exact) )
			uZ.ui.exact = 1;
		else
			uZ.ui.exact = 0;
#else
		uZ.ui = 0;
#endif
		return uZ.p;
	}

	signA = signP32UI( uiA );
	signB = signP32UI( uiB );
	signZ = signA ^ signB;
	if(signA) uiA = (-uiA & 0xFFFFFFFF);
	if(signB) uiB = (-uiB & 0xFFFFFFFF);
	regSA = signregP32UI(uiA);
	regSB = signregP32UI(uiB);

	tmp = (uiA<<2)&0xFFFFFFFF;
	if (regSA){
		while (tmp>>31){
			kA++;
			tmp= (tmp<<1) & 0xFFFFFFFF;
		}
	}
	else{
		kA=-1;
		while (!(tmp>>31)){
			kA--;
			tmp= (tmp<<1) & 0xFFFFFFFF;
		}
		tmp&=0x7FFFFFFF;
	}
	expA = tmp>>29; //to get 2 bits
	fracA = ((tmp<<1) | 0x40000000) & 0x7FFFFFFF;
	frac64A = (uint64_t) fracA << 30;

	tmp = (uiB<<2)&0xFFFFFFFF;
	if (regSB){
		while (tmp>>31){
			kA--;
			tmp= (tmp<<1) & 0xFFFFFFFF;
		}
	}
	else{
		kA++;
		while (!(tmp>>31)){
			kA++;
			tmp= (tmp<<1) & 0xFFFFFFFF;
		}
		tmp&=0x7FFFFFFF;
	}
	expA -= tmp>>29;
	fracB = ((tmp<<1) | 0x40000000) & 0x7FFFFFFF;

	divresult = lldiv (frac64A,(uint_fast64_t)fracB);
	frac64Z = divresult.quot;
	rem = divresult.rem;

	if (expA<0){
		expA+=4;
		kA--;
	}
	if (frac64Z!=0){
		rcarry = frac64Z >> 30; // this is the hidden bit (14th bit) , extreme right bit is bit 0
		if (!rcarry){
			if (expA==0){
				kA--;
				expA=3;
			}
			else
				expA--;
			frac64Z<<=1;
		}
	}

	if(kA<0){
		regA = -kA;
		regSA = 0;
		regime = 0x40000000>>regA;
	}
	else{
		regA = kA+1;
		regSA=1;
		regime = 0x7FFFFFFF - (0x7FFFFFFF>>regA);
	}
	if(regA>30){
		//max or min pos. exp and frac does not matter.
		if (signZ) (regSA) ? (uZ.ui= 0x7FFFFFFF): (uZ.ui=0x1);
		else (regSA) ? (uZ.ui= 0x7FFFFFFF): (uZ.ui=0x0);
	}
	else{
		//remove carry and rcarry bits and shift to correct position
		frac64Z &= 0x3FFFFFFF;

		fracA = (uint_fast32_t)frac64Z >> (regA+2);

		if (regA<=28){
			bitNPlusOne = (frac64Z >> (regA +1)) & 0x1;
			expA<<= (28-regA);
			if (bitNPlusOne) ( ((1<<(regA+1))-1) & frac64Z ) ? (bitsMore=1) : (bitsMore=0);
		}
		else {
			if (regA==30){
				bitNPlusOne = expA&0x2;
				bitsMore = (expA&0x1);
				expA = 0;
			}
			else if (regA==29){
				bitNPlusOne = expA&0x1;
				expA>>=1; //taken care of by the pack algo
			}
			if (frac64Z>0){
				fracA=0;
				bitsMore =1;
			}

		}
		if (rem) bitsMore =1;

		uZ.ui = packToP32UI(regime, expA, fracA);
		if (bitNPlusOne || bitsMore)
		{
			if (signZ) uZ.p = p32_nextup(uZ.p);
		} 
	}

	if (signZ) uZ.ui = -uZ.ui & 0xFFFFFFFF;
	return uZ.p;
}




// interval tools

C_R_p16 reciprocal_p16(C_R_p16 x)
{
	posit16_t c1, c2, c, tmp;
	tmp = absP16(x.center);
	c1 = p16_sub_down(F0_16, tmp);
	c1 = p16_sub_down(c1, x.radius);
	c1 = p16_div_down(Fm1_16, c1);
	c2 = p16_sub_up(F0_16, tmp);
	c2 = p16_add_up(c2, x.radius);
	c2 = p16_div_up(Fm1_16, c2);
	c = p16_sub_up(c2, c1);
	c = p16_div_up(c, F2_16);
	c = p16_add_up(c, c1);
	x.radius = p16_sub_up(c, c1);
	bool signC = signP16UI(x.center.v);
	if (signC) 
	{
		union ui16_p16 uZ;
		uZ.p = c;
		uZ.ui = -uZ.ui & 0xFFFF;
		x.center = uZ.p;
	}
	return x;
}

C_R_p32 reciprocal_p32(C_R_p32 x)
{
	posit32_t c1, c2, c, tmp;
	tmp = absP32(x.center);
	c1 = p32_sub_down(F0_32, tmp);
	c1 = p32_sub_down(c1, x.radius);
	c1 = p32_div_down(Fm1_32, c1);
	c2 = p32_sub_up(F0_32, tmp);
	c2 = p32_add_up(c2, x.radius);
	c2 = p32_div_up(Fm1_32, c2);
	c = p32_sub_up(c2, c1);
	c = p32_div_up(c, F2_32);
	c = p32_add_up(c, c1);
	x.radius = p32_sub_up(c, c1);
	bool signC = signP32UI(x.center.v);
	if (signC) 
	{
		union ui32_p32 uZ;
		uZ.p = c;
		uZ.ui = -uZ.ui & 0xFFFFFFFF;
		x.center = uZ.p;
	}
	return x;
}

C_R_p16 intersect_p16(C_R_p16 x,C_R_p16 y)
{
	posit16_t inf_x, inf_y, inf, sup_x, sup_y, sup, c, r;
	bool le1, le2, lt;
	inf_x = p16_sub_down(x.center, x.radius);
	inf_y = p16_sub_down(y.center, y.radius);
	le1 = p16_le(inf_x, inf_y); 
	if (le1)
	{
		inf = inf_y;
	}
	else
	{
		inf = inf_x;
	}
	sup_x = p16_add_up(x.center, x.radius);
	sup_y = p16_add_up(y.center, y.radius);
	le2 = p16_le(sup_x, sup_y);
	if (le2)
	{
		sup = sup_x;
	}
	else
	{
		sup = sup_y;
	}
	lt = p16_lt(sup, inf);
	if (lt)
	{
		double t1 = convertP16ToDouble(inf_x);
		double t2 = convertP16ToDouble(sup_x);
		double t3 = convertP16ToDouble(inf_y);
		double t4 = convertP16ToDouble(sup_y);
		printf("inf_x = %.17e\t", t1);
		printf("sup_x = %.17e\n", t2);
		printf("inf_y = %.17e\t", t3);
		printf("sup_y = %.17e\n", t4);
        printf("Error: intersection is empty! back to previous interval\n");
        return y;
	}
	c = p16_add_up(inf, sup);
	c = p16_div_up(c, F2_16);
	r = p16_sub_up(c, inf);
	C_R_p16 res = {c, r};
	return res;
}

C_R_p32 intersect_p32(C_R_p32 x, C_R_p32 y)
{
	posit32_t inf_x, inf_y, inf, sup_x, sup_y, sup, c, r;
	bool le1, le2, lt;
	inf_x = p32_sub_down(x.center, x.radius);
	inf_y = p32_sub_down(y.center, y.radius);
	le1 = p32_le(inf_x, inf_y); 
	if (le1)
	{
		inf = inf_y;
	}
	else
	{
		inf = inf_x;
	}
	sup_x = p32_add_up(x.center, x.radius);
	sup_y = p32_add_up(y.center, y.radius);
	le2 = p32_le(sup_x, sup_y);
	if (le2)
	{
		sup = sup_x;
	}
	else
	{
		sup = sup_y;
	}
	lt = p32_lt(sup, inf);
	if (lt)
	{
		double t1 = convertP32ToDouble(inf_x);
		double t2 = convertP32ToDouble(sup_x);
		double t3 = convertP32ToDouble(inf_y);
		double t4 = convertP32ToDouble(sup_y);
		printf("inf_x = %.17e\t", t1);
		printf("sup_x = %.17e\n", t2);
		printf("inf_y = %.17e\t", t3);
		printf("sup_y = %.17e\n", t4);
        printf("Error: intersection is empty! back to previous interval\n");
        return y;
	}
	c = p32_add_up(inf, sup);
	c = p32_div_up(c, F2_32);
	r = p32_sub_up(c, inf);
	C_R_p32 res = {c, r};
	return res;
}

C_R_p16 interval_add_p16(C_R_p16 x, C_R_p16 y)
{
	posit16_t c, r;
	posit16_t c_up, c_down, eps;
	bool lt;
	c_up = p16_add_up(x.center, y.center);
	c_down = p16_add_down(x.center, y.center);
	c = p16_add(x.center, y.center);
	c_up = p16_sub_up(c_up, c);
	c_down = p16_sub_up(c, c_down);
	lt = p16_lt(c_down, c_up);
	if (lt)
	{
		eps = p16_div_up(c_up, F2_16);
	}
	else
	{
		eps = p16_div_up(c_down, F2_16);
	}
	r = p16_add_up(x.radius, y.radius);
	r = p16_add_up(r, eps);
	C_R_p16 res = {c, r};
	return res;
}

C_R_p32 interval_add_p32(C_R_p32 x, C_R_p32 y)
{
	posit32_t c, r;
	posit32_t c_up, c_down, eps;
	bool lt;
	c_up = p32_add_up(x.center, y.center);
	c_down = p32_add_down(x.center, y.center);
	c = p32_add(x.center, y.center);
	c_up = p32_sub_up(c_up, c);
	c_down = p32_sub_up(c, c_down);
	lt = p32_lt(c_down, c_up);
	if (lt)
	{
		eps = p32_div_up(c_up, F2_32);
	}
	else
	{
		eps = p32_div_up(c_down, F2_32);
	}
	r = p32_add_up(x.radius, y.radius);
	r = p32_add_up(r, eps);
	C_R_p32 res = {c, r};
	return res;
}

C_R_p16 interval_sub_p16(C_R_p16 x, C_R_p16 y)
{
	posit16_t c, r;
	posit16_t c_up, c_down, eps;
	bool lt;
	c_up = p16_sub_up(x.center, y.center);
	c_down = p16_sub_down(x.center, y.center);
	c = p16_sub(x.center, y.center);
	c_up = p16_sub_up(c_up, c);
	c_down = p16_sub_up(c, c_down);
	lt = p16_lt(c_down, c_up);
	if (lt)
	{
		eps = p16_div_up(c_up, F2_16);
	}
	else
	{
		eps = p16_div_up(c_down, F2_16);
	}
	r = p16_add_up(x.radius, y.radius);
	r = p16_add_up(r, eps);
	C_R_p16 res = {c, r};
	return res;
}

C_R_p32 interval_sub_p32(C_R_p32 x, C_R_p32 y)
{
	posit32_t c, r;
	posit32_t c_up, c_down, eps;
	bool lt;
	c_up = p32_sub_up(x.center, y.center);
	c_down = p32_sub_down(x.center, y.center);
	c = p32_sub(x.center, y.center);
	c_up = p32_sub_up(c_up, c);
	c_down = p32_sub_up(c, c_down);
	lt = p32_lt(c_down, c_up);
	if (lt)
	{
		eps = p32_div_up(c_up, F2_32);
	}
	else
	{
		eps = p32_div_up(c_down, F2_32);
	}
	r = p32_add_up(x.radius, y.radius);
	r = p32_add_up(r, eps);
	C_R_p32 res = {c, r};
	return res;
}

C_R_p16 interval_mul_p16(C_R_p16 x, C_R_p16 y)
{
	posit16_t c, r;
	posit16_t c_up, c_down, eps, tmp1, tmp2;
	posit16_t abs1, abs2;
	bool lt;
	c_up = p16_mul_up(x.center, y.center);
	c_down = p16_mul_down(x.center, y.center);
	c = p16_mul(x.center, y.center);
	c_up = p16_sub_up(c_up, c);
	c_down = p16_sub_up(c, c_down);
	lt = p16_lt(c_down, c_up);
	if (lt)
	{
		eps = p16_div_up(c_up, F2_16);
	}
	else
	{
		eps = p16_div_up(c_down, F2_16);
	}	
	abs1 = absP16(x.center);
	abs2 = absP16(y.center);
	tmp1 = p16_add_up(abs1, x.radius);
	tmp1 = p16_mul_up(tmp1, y.radius);
	tmp2 = p16_mul_up(x.radius, abs2);
	r = p16_add_up(tmp1, tmp2);
	r = p16_add_up(eps, r);
}

C_R_p32 interval_mul_p32(C_R_p32 x, C_R_p32 y)
{
	posit32_t c, r;
	posit32_t c_up, c_down, eps, tmp1, tmp2;
	posit32_t abs1, abs2;
	bool lt;
	c_up = p32_mul_up(x.center, y.center);
	c_down = p32_mul_down(x.center, y.center);
	c = p32_mul(x.center, y.center);
	c_up = p32_sub_up(c_up, c);
	c_down = p32_sub_up(c, c_down);
	lt = p32_lt(c_down, c_up);
	if (lt)
	{
		eps = p32_div_up(c_up, F2_32);
	}
	else
	{
		eps = p32_div_up(c_down, F2_32);
	}	
	abs1 = absP32(x.center);
	abs2 = absP32(y.center);
	tmp1 = p32_add_up(abs1, x.radius);
	tmp1 = p32_mul_up(tmp1, y.radius);
	tmp2 = p32_mul_up(x.radius, abs2);
	r = p32_add_up(tmp1, tmp2);
	r = p32_add_up(eps, r);
}


// A in band storage format diagonal n, subdiagonal n-1, superdiagonal n-1
// b is a vector of size n
void GS_tridiag_p16(int n, C_R_p16 *A, C_R_p16 *b, C_R_p16 *x)
{
    C_R_p16 *x_prev = malloc(n * sizeof(C_R_p16));
    for (int i = 0; i < n; i ++)
    {
        x_prev[i] = x[i]; 
    }

    int count = 0;

    for (int i = 0; i < 1000; i ++) 
    {  
        for (int i = 0; i < n; i ++) 
        {
            C_R_p16 sum = b[i];

            if (i > 0) 
            {
                C_R_p16 prod = interval_mul_p16(A[n + i - 1], x[i - 1]);
                sum = interval_sub_p16(sum, prod);
            }

            if (i < n - 1) 
            {
                C_R_p16 prod = interval_mul_p16(A[2 * n - 1 + i], x_prev[i + 1]);
                sum = interval_sub_p16(sum, prod);
            }

            C_R_p16 Aii_inv = reciprocal_p16(A[i]);

            x[i] = interval_mul_p16(sum, Aii_inv); 
            x[i] = intersect_p16(x[i], x_prev[i]);
        }

        count ++;
		printf("iteration %d\n", count);

        // Check convergence
        bool cvg = true;
        for (int i = 0; i < n; i ++)
        {
			union ui16_p16 c1, c2, r1, r2;
			c1.p = x[i].center;
			c2.p = x_prev[i].center;
			r1.p = x[i].radius;
			r2.p = x_prev[i].radius;
			if (c1.ui != c2.ui || r1.ui != r2.ui)
			{
				cvg = false;
				break;
			}
        }

        if (cvg)
        {
            break;
        }

        for (int i = 0; i < n; i ++)
        {
            x_prev[i] = x[i];
        }
    }

    free(x_prev);
}

void GS_tridiag_p32(int n, C_R_p32 *A, C_R_p32 *b, C_R_p32 *x)
{
    C_R_p32 *x_prev = malloc(n * sizeof(C_R_p32));
    for (int i = 0; i < n; i ++)
    {
        x_prev[i] = x[i]; 
    }

    int count = 0;

    for (int i = 0; i < 1000; i ++) 
    {  
        for (int i = 0; i < n; i ++) 
        {
            C_R_p32 sum = b[i];

            if (i > 0) 
            {
                C_R_p32 prod = interval_mul_p32(A[n + i - 1], x[i - 1]);
                sum = interval_sub_p32(sum, prod);
            }

            if (i < n - 1) 
            {
                C_R_p32 prod = interval_mul_p32(A[2 * n - 1 + i], x_prev[i + 1]);
                sum = interval_sub_p32(sum, prod);
            }

            C_R_p32 Aii_inv = reciprocal_p32(A[i]);

            x[i] = interval_mul_p32(sum, Aii_inv); 
            x[i] = intersect_p32(x[i], x_prev[i]);
        }

        count ++;

        // Check convergence
        bool cvg = true;
        for (int i = 0; i < n; i ++)
        {
			union ui32_p32 c1, c2, r1, r2;
			c1.p = x[i].center;
			c2.p = x_prev[i].center;
			r1.p = x[i].radius;
			r2.p = x_prev[i].radius;
			if (c1.ui != c2.ui || r1.ui != r2.ui)
			{
				cvg = false;
				break;
			}
        }

        if (cvg)
        {
            break;
        }

        for (int i = 0; i < n; i ++)
        {
            x_prev[i] = x[i];
        }
    }
	free(x_prev);
}


// A in CSR format
// idx is the index of the first non-zero element in each row
// col_id is the column index of each non-zero element
// b is a vector of size n
void GS_csr_p16(int n, C_R_p16 *A, int *idx, int *col_id, C_R_p16 *b, C_R_p16 *x)
{
    C_R_p16 *x_prev = malloc(n * sizeof(C_R_p16));
    for (int i = 0; i < n; i ++)
    {
        x_prev[i] = x[i]; 
    }

    for (int iter = 0; iter < 1000; iter ++) 
    {  
        for (int i = 0; i < n; i ++) 
        {
            C_R_p16 sum = b[i];
            C_R_p16 Aii_inv = { {0.0}, {0.0} };

            for (int k = idx[i]; k < idx[i + 1]; k ++)
            {
                int j = col_id[k];

                if (j < i) 
                {
                    C_R_p16 prod = interval_mul_p16(A[k], x[j]);
                    sum = interval_sub_p16(sum, prod);
                }
                else if (j == i)
                {
                    Aii_inv = reciprocal_p16(A[k]);
                }
                else if (j > i) 
                {
                    C_R_p16 prod = interval_mul_p16(A[k], x_prev[j]);
                    sum = interval_sub_p16(sum, prod);
                }              
            }
            x[i] = interval_mul_p16(sum, Aii_inv);  // x_i^{(k+1)} = sum / A_ii

            x[i] = intersect_p16(x[i], x_prev[i]);
        }

        // Check convergence
        bool cvg = true;
        for (int i = 0; i < n; i ++)
        {
			union ui16_p16 c1, c2, r1, r2;
			c1.p = x[i].center;
			c2.p = x_prev[i].center;
			r1.p = x[i].radius;
			r2.p = x_prev[i].radius;
			if (c1.ui != c2.ui || r1.ui != r2.ui)
			{
				cvg = false;
				break;
			}
        }

        if (cvg)
        {
            break;
        }

        for (int i = 0; i < n; i ++)
        {
            x_prev[i] = x[i]; 
        }
    }

    free(x_prev);
}

void GS_csr_p32(int n, C_R_p32 *A, int *idx, int *col_id, C_R_p32 *b, C_R_p32 *x)
{
    C_R_p32 *x_prev = malloc(n * sizeof(C_R_p32));
    for (int i = 0; i < n; i ++)
    {
        x_prev[i] = x[i]; 
    }

    for (int iter = 0; iter < 1000; iter ++) 
    {  
        for (int i = 0; i < n; i ++) 
        {
            C_R_p32 sum = b[i];
            C_R_p32 Aii_inv = { {0.0}, {0.0} };

            for (int k = idx[i]; k < idx[i + 1]; k ++)
            {
                int j = col_id[k];

                if (j < i) 
                {
                    C_R_p32 prod = interval_mul_p32(A[k], x[j]);
                    sum = interval_sub_p32(sum, prod);
                }
                else if (j == i)
                {
                    Aii_inv = reciprocal_p32(A[k]);
                }
                else if (j > i) 
                {
                    C_R_p32 prod = interval_mul_p32(A[k], x_prev[j]);
                    sum = interval_sub_p32(sum, prod);
                }              
            }
            x[i] = interval_mul_p32(sum, Aii_inv);  // x_i^{(k+1)} = sum / A_ii

            x[i] = intersect_p32(x[i], x_prev[i]);
        }

        // Check convergence
        bool cvg = true;
        for (int i = 0; i < n; i ++)
        {
			union ui32_p32 c1, c2, r1, r2;
			c1.p = x[i].center;
			c2.p = x_prev[i].center;
			r1.p = x[i].radius;
			r2.p = x_prev[i].radius;
			if (c1.ui != c2.ui || r1.ui != r2.ui)
			{
				cvg = false;
				break;
			}
        }

        if (cvg)
        {
            break;
        }

        for (int i = 0; i < n; i ++)
        {
            x_prev[i] = x[i]; 
        }
    }

    free(x_prev);
}

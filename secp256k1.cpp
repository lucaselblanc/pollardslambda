/***********************************************************************************************************
* This file is part of the Pollard's Lambda distribution: (https://github.com/lucaselblanc/pollardslambda) *
* Copyright (c) 2024, 2026 Lucas Leblanc.                                                                  *
* Distributed under the MIT software license, see the accompanying.                                        *
* file COPYING or https://www.opensource.org/licenses/mit-license.php.                                     *
************************************************************************************************************/

/*******************************************
* Pollard's Lambda Algorithm for SECP256K1 *
* Written by Lucas Leblanc                 *
********************************************/

#include "secp256k1.h"

#ifndef _addcarry_u64
#define _addcarry_u64(c_in, a, b, out) __builtin_ia32_addcarryx_u64(c_in, a, b, (unsigned long long*)(out))
#endif

#ifndef _subborrow_u64
#define _subborrow_u64(b_in, a, b, out) __builtin_ia32_sbb_u64(b_in, a, b, (unsigned long long*)(out))
#endif

static inline int64_t _mul128(int64_t a, int64_t b, int64_t *h) {
    __int128 p = (__int128)a * b;
    *h = (int64_t)(p >> 64);
    return (int64_t)p;
}

#define __shiftright128(a,b,n) (((a)>>(n))|((b)<<(64-(n))))
#define __shiftleft128(a,b,n) (((b)<<(n))|((a)>>(64-(n))))
#define SWAP_STEP(tmp,x,y) do { tmp = x; x = y; y = tmp; } while(0)
#define MSK62 0x3FFFFFFFFFFFFFFFULL

static inline void secp256k1_reduce_final(uint64_t *r, unsigned char carry) {
    uint64_t diff[4];
    unsigned char b = 0;
    b = _subborrow_u64(0, r[0], 0xFFFFFFFEFFFFFC2FULL, &diff[0]);
    b = _subborrow_u64(b, r[1], 0xFFFFFFFFFFFFFFFFULL, &diff[1]);
    b = _subborrow_u64(b, r[2], 0xFFFFFFFFFFFFFFFFULL, &diff[2]);
    b = _subborrow_u64(b, r[3], 0xFFFFFFFFFFFFFFFFULL, &diff[3]);
    if (carry || b == 0) {
        r[0] = diff[0];
        r[1] = diff[1];
        r[2] = diff[2];
        r[3] = diff[3];
    }
}

constexpr uint64_t P_CONST[4] = { 0xFFFFFFFEFFFFFC2FULL, 0xFFFFFFFFFFFFFFFFULL, 0xFFFFFFFFFFFFFFFFULL, 0xFFFFFFFFFFFFFFFFULL };
constexpr uint64_t N_CONST[4] = { 0xBFD25E8CD0364141ULL, 0xBAAEDCE6AF48A03BULL, 0xFFFFFFFFFFFFFFFEULL, 0xFFFFFFFFFFFFFFFFULL };
constexpr uint64_t GX_CONST[4] = { 0x59F2815B16F81798ULL, 0x029BFCDB2DCE28D9ULL, 0x55A06295CE870B07ULL, 0x79BE667EF9DCBBACULL };
constexpr uint64_t GY_CONST[4] = { 0x9C47D08FFB10D4B8ULL, 0xFD17B448A6855419ULL, 0x5DA4FBFC0E1108A8ULL, 0x483ADA7726A3C465ULL };
constexpr uint64_t R2_MOD_P[4] = { 0x000007A2000E90A1, 0x0000000000000001, 0x0ULL, 0x0ULL };
constexpr uint64_t ZERO_MONT[4] = { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL };
constexpr uint64_t ONE_MONT[4] = { 0x00000001000003D1ULL, 0x0ULL, 0x0ULL, 0x0ULL };
constexpr uint64_t SEVEN_MONT[4] = {0x700001AB7ULL, 0x0ULL, 0x0ULL, 0x0ULL};
constexpr uint64_t P_CONST_MINUS_2[4] = { 0xFFFFFFFEFFFFFC2DULL, 0xFFFFFFFFFFFFFFFFULL, 0xFFFFFFFFFFFFFFFFULL, 0xFFFFFFFFFFFFFFFFULL };
constexpr uint64_t LAMBDA_N[4] = { 0xDF02967C1B23BD72ULL, 0xA5261C028812645AULL, 0xA5261C028812645AULL, 0x5363AD4CC05C30E0ULL };
constexpr uint64_t BETA_P[4] = { 0xB315ECECBB640683ULL, 0x9CF0497512F58995ULL, 0x6E64479EAC3434E9ULL, 0x7AE96A2B657C0710ULL };
constexpr uint64_t MINUS_B1[4] = { 0x6F547FA90ABFE4C3ULL, 0xE4437ED6010E8828ULL, 0x0ULL, 0x0ULL };
constexpr uint64_t MINUS_B2[4] = { 0x8A280AC50774346DULL, 0xD765CDA83DB1562CULL, 0xCFFFFFFFFFFFFFFEULL, 0xFFFFFFFFFFFFFFFFULL };
constexpr uint64_t G1[4] = { 0xE893209A45DBB031ULL, 0x3DAA8A1471E8CA7FULL, 0xE86C90E49284EB15ULL, 0x3086D221A7D46BCDULL };
constexpr uint64_t G2[4] = { 0x1571B4AE8AC47F71ULL, 0x221208AC9DF506C6ULL, 0x6F547FA90ABFE4C4ULL, 0xE4437ED6010E8828ULL };
constexpr uint64_t MU_P = 0xD838091DD2253531ULL;
ECPointJacobian* preCompG = nullptr;
ECPointJacobian* preCompGphi = nullptr;
ECPointJacobian* preCompH = nullptr;
ECPointJacobian* preCompHphi = nullptr;
ECPointJacobian* jacNorm = nullptr;
ECPointJacobian* jacEndo = nullptr;
ECPointJacobian* jacNormH = nullptr;
ECPointJacobian* jacEndoH = nullptr;

void montgomeryReduceP(uint64_t *result, const uint64_t *inputHigh, const uint64_t *inputLow) {
    uint64_t temp[8];
    for (int i = 0; i < 4; i++) {
        temp[i]     =  inputLow[i];
        temp[i + 4] = inputHigh[i];
    }

    uint64_t extra = 0;

    for (int i = 0; i < 4; i++) {
        uint64_t ui = (uint64_t)((uint128_t)temp[i] * (uint128_t)MU_P);
        uint128_t carry = 0;
        for (int j = 0; j < 4; j++) {
            uint128_t prod = (uint128_t)ui * (uint128_t)P_CONST[j] + (uint128_t)temp[i + j] + carry;
            temp[i + j] = (uint64_t)prod;
            carry = prod >> 64;
        }

        uint128_t s = (uint128_t)temp[i + 4] + carry;
        temp[i + 4] = (uint64_t)s;
        carry = s >> 64;

        for (int j = i + 5; j < 8; ++j) {
            uint128_t sum = (uint128_t)temp[j] + carry;
            temp[j] = (uint64_t)sum;
            carry = sum >> 64;
        }

        extra += (uint64_t)carry;
    }

    for (int i = 0; i < 4; i++) result[i] = temp[i + 4];

    uint64_t diff[4];
    uint128_t borrow = 0;
    for (int i = 0; i < 4; i++) {
        uint128_t sub = (uint128_t)result[i] - (uint128_t)P_CONST[i] - borrow;
        diff[i] = (uint64_t)sub;
        borrow = (sub >> 127) & 1;
    }

    if (extra != 0 || borrow == 0) {
        for (int i = 0; i < 4; i++) result[i] = diff[i];
    }
}

static inline void imm_umul(const uint64_t *x, uint64_t y, uint64_t *dst) {
    uint64_t h, carry;
    unsigned char c = 0;
    dst[0] = _umul128(x[0], y, &h); carry = h;
    c = _addcarry_u64(c, _umul128(x[1], y, &h), carry, &dst[1]); carry = h;
    c = _addcarry_u64(c, _umul128(x[2], y, &h), carry, &dst[2]); carry = h;
    c = _addcarry_u64(c, _umul128(x[3], y, &h), carry, &dst[3]); carry = h;
    _addcarry_u64(c, 0ULL, carry, &dst[4]);
}

void toMontgomeryP(uint64_t *result, const uint64_t *a) {
    uint64_t aLocal[4];
    for (int i = 0; i < 4; i++) aLocal[i] = a[i];

    uint64_t temp[8] = {0};

    for (int i = 0; i < 4; i++) {
        uint64_t prod[5];
        imm_umul(R2_MOD_P, aLocal[i], prod);

        unsigned char c = 0;
        c = _addcarry_u64(c, temp[i + 0], prod[0], &temp[i + 0]);
        c = _addcarry_u64(c, temp[i + 1], prod[1], &temp[i + 1]);
        c = _addcarry_u64(c, temp[i + 2], prod[2], &temp[i + 2]);
        c = _addcarry_u64(c, temp[i + 3], prod[3], &temp[i + 3]);
        c = _addcarry_u64(c, temp[i + 4], prod[4], &temp[i + 4]);

        for (int k = i + 5; k < 8 && c; k++) {
            c = _addcarry_u64(c, temp[k], 0, &temp[k]);
        }
    }

    uint64_t low[4], high[4];
    for (int i = 0; i < 4; i++) {
        low[i] = temp[i];
        high[i] = temp[i + 4];
    }

    montgomeryReduceP(result, high, low);
}

void fromMontgomeryP(uint64_t *result, const uint64_t *a) {
    uint64_t zero[4] = {0, 0, 0, 0};
    montgomeryReduceP(result, zero, a);
}

void modMulK1(uint64_t *result, const uint64_t *a, const uint64_t *b) {
    asm volatile (
        "movq 0(%[b]), %%rdx\n\t"
        "mulx 0(%[a]), %%r8, %%r9\n\t"
        "mulx 8(%[a]), %%rax, %%r10\n\t"
        "addq %%rax, %%r9\n\t"
        "mulx 16(%[a]), %%rax, %%r11\n\t"
        "adcq %%rax, %%r10\n\t"
        "mulx 24(%[a]), %%rax, %%r12\n\t"
        "adcq %%rax, %%r11\n\t"
        "adcq $0, %%r12\n\t"

        "movq 8(%[b]), %%rdx\n\t"
        "xorq %%rax, %%rax\n\t"
        "mulx 0(%[a]), %%rax, %%rbx\n\t"
        "adcxq %%rax, %%r9\n\t"
        "adoxq %%rbx, %%r10\n\t"
        "mulx 8(%[a]), %%rax, %%rbx\n\t"
        "adcxq %%rax, %%r10\n\t"
        "adoxq %%rbx, %%r11\n\t"
        "mulx 16(%[a]), %%rax, %%rbx\n\t"
        "adcxq %%rax, %%r11\n\t"
        "adoxq %%rbx, %%r12\n\t"
        "mulx 24(%[a]), %%rax, %%r13\n\t"
        "adcxq %%rax, %%r12\n\t"
        "movq $0, %%rax\n\t"
        "adcxq %%rax, %%r13\n\t"
        "adoxq %%rax, %%r13\n\t"

        "movq 16(%[b]), %%rdx\n\t"
        "xorq %%rax, %%rax\n\t"
        "mulx 0(%[a]), %%rax, %%rbx\n\t"
        "adcxq %%rax, %%r10\n\t"
        "adoxq %%rbx, %%r11\n\t"
        "mulx 8(%[a]), %%rax, %%rbx\n\t"
        "adcxq %%rax, %%r11\n\t"
        "adoxq %%rbx, %%r12\n\t"
        "mulx 16(%[a]), %%rax, %%rbx\n\t"
        "adcxq %%rax, %%r12\n\t"
        "adoxq %%rbx, %%r13\n\t"
        "mulx 24(%[a]), %%rax, %%r14\n\t"
        "adcxq %%rax, %%r13\n\t"
        "movq $0, %%rax\n\t"
        "adcxq %%rax, %%r14\n\t"
        "adoxq %%rax, %%r14\n\t"

        "movq 24(%[b]), %%rdx\n\t"
        "xorq %%rax, %%rax\n\t"
        "mulx 0(%[a]), %%rax, %%rbx\n\t"
        "adcxq %%rax, %%r11\n\t"
        "adoxq %%rbx, %%r12\n\t"
        "mulx 8(%[a]), %%rax, %%rbx\n\t"
        "adcxq %%rax, %%r12\n\t"
        "adoxq %%rbx, %%r13\n\t"
        "mulx 16(%[a]), %%rax, %%rbx\n\t"
        "adcxq %%rax, %%r13\n\t"
        "adoxq %%rbx, %%r14\n\t"
        "mulx 24(%[a]), %%rax, %%r15\n\t"
        "adcxq %%rax, %%r14\n\t"
        "movq $0, %%rax\n\t"
        "adcxq %%rax, %%r15\n\t"
        "adoxq %%rax, %%r15\n\t"

        "movabsq $0x1000003D1, %%rdx\n\t"
        "mulx %%r12, %%rax, %%rbx\n\t"
        "addq %%rax, %%r8\n\t"
        "adcq %%rbx, %%r9\n\t"

        "mulx %%r13, %%rax, %%rbx\n\t"
        "adcq %%rax, %%r9\n\t"
        "adcq %%rbx, %%r10\n\t"

        "mulx %%r14, %%rax, %%rbx\n\t"
        "adcq %%rax, %%r10\n\t"
        "adcq %%rbx, %%r11\n\t"

        "mulx %%r15, %%rax, %%r12\n\t"
        "adcq %%rax, %%r11\n\t"
        "adcq $0, %%r12\n\t"

        "mulx %%r12, %%rax, %%rbx\n\t"
        "addq %%rax, %%r8\n\t"
        "adcq %%rbx, %%r9\n\t"
        "adcq $0, %%r10\n\t"
        "adcq $0, %%r11\n\t"

        "movq $0, %%rax\n\t"
        "movabsq $0x1000003D1, %%rcx\n\t"
        "cmovc %%rcx, %%rax\n\t"
        "addq %%rax, %%r8\n\t"
        "adcq $0, %%r9\n\t"
        "adcq $0, %%r10\n\t"
        "adcq $0, %%r11\n\t"

        "movq %%r8, %%rax\n\t"
        "addq %%rcx, %%rax\n\t"
        "movq %%r9, %%rbx\n\t"
        "adcq $0, %%rbx\n\t"
        "movq %%r10, %%r12\n\t"
        "adcq $0, %%r12\n\t"
        "movq %%r11, %%r13\n\t"
        "adcq $0, %%r13\n\t"

        "cmovc %%rax, %%r8\n\t"
        "cmovc %%rbx, %%r9\n\t"
        "cmovc %%r12, %%r10\n\t"
        "cmovc %%r13, %%r11\n\t"

        "movq %%r8, 0(%[result])\n\t"
        "movq %%r9, 8(%[result])\n\t"
        "movq %%r10, 16(%[result])\n\t"
        "movq %%r11, 24(%[result])\n\t"
        :
        : [result] "r" (result), [a] "r" (a), [b] "r" (b)
        : "rax", "rbx", "rcx", "rdx", "r8", "r9", "r10", "r11", "r12", "r13", "r14", "r15", "memory", "cc"
    );
}

void modSquareK1(uint64_t *r, const uint64_t *a) {
    asm volatile (
        "movq 0(%[a]), %%rdx\n\t"
        "mulx 8(%[a]), %%r9, %%r10\n\t"
        "mulx 16(%[a]), %%rax, %%r11\n\t"
        "addq %%rax, %%r10\n\t"
        "mulx 24(%[a]), %%rax, %%r12\n\t"
        "adcq %%rax, %%r11\n\t"
        "adcq $0, %%r12\n\t"

        "movq 8(%[a]), %%rdx\n\t"
        "xorq %%rax, %%rax\n\t"
        "mulx 16(%[a]), %%rax, %%rbx\n\t"
        "adcxq %%rax, %%r11\n\t"
        "adoxq %%rbx, %%r12\n\t"
        "mulx 24(%[a]), %%rax, %%r13\n\t"
        "adcxq %%rax, %%r12\n\t"
        "movq $0, %%rax\n\t"
        "adcxq %%rax, %%r13\n\t"
        "adoxq %%rax, %%r13\n\t"

        "movq 16(%[a]), %%rdx\n\t"
        "mulx 24(%[a]), %%rax, %%r14\n\t"
        "addq %%rax, %%r13\n\t"
        "adcq $0, %%r14\n\t"

        "xorq %%r15, %%r15\n\t"
        "addq %%r9, %%r9\n\t"
        "adcq %%r10, %%r10\n\t"
        "adcq %%r11, %%r11\n\t"
        "adcq %%r12, %%r12\n\t"
        "adcq %%r13, %%r13\n\t"
        "adcq %%r14, %%r14\n\t"
        "adcq $0, %%r15\n\t"

        "movq 0(%[a]), %%rdx\n\t"
        "mulx %%rdx, %%r8, %%rax\n\t"
        "addq %%rax, %%r9\n\t"

        "movq 8(%[a]), %%rdx\n\t"
        "mulx %%rdx, %%rax, %%rbx\n\t"
        "adcq %%rax, %%r10\n\t"
        "adcq %%rbx, %%r11\n\t"

        "movq 16(%[a]), %%rdx\n\t"
        "mulx %%rdx, %%rax, %%rbx\n\t"
        "adcq %%rax, %%r12\n\t"
        "adcq %%rbx, %%r13\n\t"

        "movq 24(%[a]), %%rdx\n\t"
        "mulx %%rdx, %%rax, %%rbx\n\t"
        "adcq %%rax, %%r14\n\t"
        "adcq %%rbx, %%r15\n\t"

        "movabsq $0x1000003D1, %%rdx\n\t"
        "mulx %%r12, %%rax, %%rbx\n\t"
        "addq %%rax, %%r8\n\t"
        "adcq %%rbx, %%r9\n\t"

        "mulx %%r13, %%rax, %%rbx\n\t"
        "adcq %%rax, %%r9\n\t"
        "adcq %%rbx, %%r10\n\t"

        "mulx %%r14, %%rax, %%rbx\n\t"
        "adcq %%rax, %%r10\n\t"
        "adcq %%rbx, %%r11\n\t"

        "mulx %%r15, %%rax, %%r12\n\t"
        "adcq %%rax, %%r11\n\t"
        "adcq $0, %%r12\n\t"

        "mulx %%r12, %%rax, %%rbx\n\t"
        "addq %%rax, %%r8\n\t"
        "adcq %%rbx, %%r9\n\t"
        "adcq $0, %%r10\n\t"
        "adcq $0, %%r11\n\t"

        "movq $0, %%rax\n\t"
        "movabsq $0x1000003D1, %%rcx\n\t"
        "cmovc %%rcx, %%rax\n\t"
        "addq %%rax, %%r8\n\t"
        "adcq $0, %%r9\n\t"
        "adcq $0, %%r10\n\t"
        "adcq $0, %%r11\n\t"

        "movq %%r8, %%rax\n\t"
        "addq %%rcx, %%rax\n\t"
        "movq %%r9, %%rbx\n\t"
        "adcq $0, %%rbx\n\t"
        "movq %%r10, %%r12\n\t"
        "adcq $0, %%r12\n\t"
        "movq %%r11, %%r13\n\t"
        "adcq $0, %%r13\n\t"

        "cmovc %%rax, %%r8\n\t"
        "cmovc %%rbx, %%r9\n\t"
        "cmovc %%r12, %%r10\n\t"
        "cmovc %%r13, %%r11\n\t"

        "movq %%r8, 0(%[r])\n\t"
        "movq %%r9, 8(%[r])\n\t"
        "movq %%r10, 16(%[r])\n\t"
        "movq %%r11, 24(%[r])\n\t"
        :
        : [r] "r" (r), [a] "r" (a)
        : "rax", "rbx", "rcx", "rdx", "r8", "r9", "r10", "r11", "r12", "r13", "r14", "r15", "memory", "cc"
    );
}

struct Int5 {
    uint64_t d[5];

    inline void clear() { memset(d, 0, sizeof(d)); }
    inline bool isZero() const { return (d[0] | d[1] | d[2] | d[3] | d[4]) == 0; }
    inline bool isOne() const { return d[0] == 1 && (d[1] | d[2] | d[3] | d[4]) == 0; }
    inline bool isNegative() const { return ((int64_t)d[4]) < 0; }

    inline void neg() {
        unsigned char c = 0;
        c = _subborrow_u64(c, 0, d[0], &d[0]);
        c = _subborrow_u64(c, 0, d[1], &d[1]);
        c = _subborrow_u64(c, 0, d[2], &d[2]);
        c = _subborrow_u64(c, 0, d[3], &d[3]);
        _subborrow_u64(c, 0, d[4], &d[4]);
    }

    inline void add(const Int5 *a, const Int5 *b) {
        unsigned char c = 0;
        c = _addcarry_u64(c, a->d[0], b->d[0], &d[0]);
        c = _addcarry_u64(c, a->d[1], b->d[1], &d[1]);
        c = _addcarry_u64(c, a->d[2], b->d[2], &d[2]);
        c = _addcarry_u64(c, a->d[3], b->d[3], &d[3]);
        _addcarry_u64(c, a->d[4], b->d[4], &d[4]);
    }

    inline uint64_t addCh(const Int5 *a, uint64_t ca, const Int5 *b, uint64_t cb) {
        uint64_t carry;
        unsigned char c = 0;
        c = _addcarry_u64(c, a->d[0], b->d[0], &d[0]);
        c = _addcarry_u64(c, a->d[1], b->d[1], &d[1]);
        c = _addcarry_u64(c, a->d[2], b->d[2], &d[2]);
        c = _addcarry_u64(c, a->d[3], b->d[3], &d[3]);
        c = _addcarry_u64(c, a->d[4], b->d[4], &d[4]);
        _addcarry_u64(c, ca, cb, &carry);
        return carry;
    }

    inline uint64_t addCh(const Int5 *a, uint64_t ca) {
        uint64_t carry;
        unsigned char c = 0;
        c = _addcarry_u64(c, d[0], a->d[0], &d[0]);
        c = _addcarry_u64(c, d[1], a->d[1], &d[1]);
        c = _addcarry_u64(c, d[2], a->d[2], &d[2]);
        c = _addcarry_u64(c, d[3], a->d[3], &d[3]);
        c = _addcarry_u64(c, d[4], a->d[4], &d[4]);
        _addcarry_u64(c, ca, 0, &carry);
        return carry;
    }

    inline uint64_t mult(const Int5 *a, uint64_t b) {
        unsigned char c = 0;
        uint64_t h, carry;
        d[0] = _umul128(a->d[0], b, &h); carry = h;
        c = _addcarry_u64(c, _umul128(a->d[1], b, &h), carry, &d[1]); carry = h;
        c = _addcarry_u64(c, _umul128(a->d[2], b, &h), carry, &d[2]); carry = h;
        c = _addcarry_u64(c, _umul128(a->d[3], b, &h), carry, &d[3]); carry = h;
        c = _addcarry_u64(c, _umul128(a->d[4], b, &h), carry, &d[4]); carry = h;
        return carry;
    }

    inline uint64_t iMult(const Int5 *a, int64_t b) {
        uint64_t carry;
        Int5 tmp;
        if (b < 0) {
            unsigned char c = 0;
            c = _subborrow_u64(c, 0, a->d[0], &tmp.d[0]);
            c = _subborrow_u64(c, 0, a->d[1], &tmp.d[1]);
            c = _subborrow_u64(c, 0, a->d[2], &tmp.d[2]);
            c = _subborrow_u64(c, 0, a->d[3], &tmp.d[3]);
            c = _subborrow_u64(c, 0, a->d[4], &tmp.d[4]);
            b = -b;
        } else {
            tmp = *a;
        }
        unsigned char c = 0;
        uint64_t h;
        d[0] = _umul128(tmp.d[0], (uint64_t)b, &h); carry = h;
        c = _addcarry_u64(c, _umul128(tmp.d[1], (uint64_t)b, &h), carry, &d[1]); carry = h;
        c = _addcarry_u64(c, _umul128(tmp.d[2], (uint64_t)b, &h), carry, &d[2]); carry = h;
        c = _addcarry_u64(c, _umul128(tmp.d[3], (uint64_t)b, &h), carry, &d[3]); carry = h;
        int64_t h_signed;
        c = _addcarry_u64(c, (uint64_t)_mul128((int64_t)tmp.d[4], b, &h_signed), carry, &d[4]);
        carry = (uint64_t)h_signed;
        return carry;
    }
};

static inline void shiftR62(uint64_t *d) {
    d[0] = __shiftright128(d[0], d[1], 62);
    d[1] = __shiftright128(d[1], d[2], 62);
    d[2] = __shiftright128(d[2], d[3], 62);
    d[3] = __shiftright128(d[3], d[4], 62);
    d[4] = ((int64_t)d[4]) >> 62;
}

static inline void shiftR62(uint64_t *d, uint64_t h) {
    d[0] = __shiftright128(d[0], d[1], 62);
    d[1] = __shiftright128(d[1], d[2], 62);
    d[2] = __shiftright128(d[2], d[3], 62);
    d[3] = __shiftright128(d[3], d[4], 62);
    d[4] = __shiftright128(d[4], h, 62);
}

static inline void matrixVecMul(Int5 *u, Int5 *v, int64_t _11, int64_t _12, int64_t _21, int64_t _22) {
    Int5 t1, t2, t3, t4;
    t1.iMult(u, _11);
    t2.iMult(v, _12);
    t3.iMult(u, _21);
    t4.iMult(v, _22);
    u->add(&t1, &t2);
    v->add(&t3, &t4);
}

static inline void matrixVecMul(Int5 *u, Int5 *v, int64_t _11, int64_t _12, int64_t _21, int64_t _22, uint64_t *cu, uint64_t *cv) {
    Int5 t1, t2, t3, t4;
    uint64_t c1 = t1.iMult(u, _11);
    uint64_t c2 = t2.iMult(v, _12);
    uint64_t c3 = t3.iMult(u, _21);
    uint64_t c4 = t4.iMult(v, _22);
    *cu = u->addCh(&t1, c1, &t2, c2);
    *cv = v->addCh(&t3, c3, &t4, c4);
}

static inline void divStep62(Int5 *u, Int5 *v, int *pos, int64_t *uu, int64_t *uv, int64_t *vu, int64_t *vv) {
    uint64_t u0 = u->d[0];
    uint64_t v0 = v->d[0];
    uint64_t uh, vh, w, x;

    while (*pos >= 1 && (u->d[*pos] | v->d[*pos]) == 0) (*pos)--;
    if (*pos == 0) {
        uh = u->d[0];
        vh = v->d[0];
    } else {
        uint64_t s = __builtin_clzll(u->d[*pos] | v->d[*pos]);
        if (s == 0) {
            uh = u->d[*pos];
            vh = v->d[*pos];
        } else {
            uh = __shiftleft128(u->d[*pos - 1], u->d[*pos], (uint8_t)s);
            vh = __shiftleft128(v->d[*pos - 1], v->d[*pos], (uint8_t)s);
        }
    }

    int bitCount = 62;
    __m128i _u, _v, _t;
    ((int64_t*)&_u)[0] = 1; ((int64_t*)&_u)[1] = 0;
    ((int64_t*)&_v)[0] = 0; ((int64_t*)&_v)[1] = 1;

    while (true) {
        uint64_t zeros = __builtin_ctzll(v0 | (1ULL << bitCount));
        vh >>= zeros;
        v0 >>= zeros;
        _u = _mm_slli_epi64(_u, (int)zeros);
        bitCount -= (int)zeros;
        if (bitCount <= 0) break;

        if (vh < uh) {
            SWAP_STEP(w, uh, vh);
            SWAP_STEP(x, u0, v0);
            SWAP_STEP(_t, _u, _v);
        }
        vh -= uh;
        v0 -= u0;
        _v = _mm_sub_epi64(_v, _u);
    }

    *uu = ((int64_t*)&_u)[0];
    *uv = ((int64_t*)&_u)[1];
    *vu = ((int64_t*)&_v)[0];
    *vv = ((int64_t*)&_v)[1];
}

void modinvP(uint64_t *result, const uint64_t *a) {
    Int5 mod;
    mod.d[0] = 0xFFFFFFFEFFFFFC2FULL;
    mod.d[1] = 0xFFFFFFFFFFFFFFFFULL;
    mod.d[2] = 0xFFFFFFFFFFFFFFFFULL;
    mod.d[3] = 0xFFFFFFFFFFFFFFFFULL;
    mod.d[4] = 0;

    Int5 u = mod;
    Int5 v;
    v.d[0] = a[0]; v.d[1] = a[1]; v.d[2] = a[2]; v.d[3] = a[3]; v.d[4] = 0;
    Int5 r; r.clear();
    Int5 s; s.clear(); s.d[0] = 1;

    int64_t uu, uv, vu, vv;
    uint64_t carryR, carryS;
    int pos = 4;
    while (pos >= 1 && (u.d[pos] | v.d[pos]) == 0) pos--;

    Int5 r0_P, s0_P;

    while (!v.isZero()) {
        divStep62(&u, &v, &pos, &uu, &uv, &vu, &vv);

        matrixVecMul(&u, &v, uu, uv, vu, vv);
        if (u.isNegative()) { u.neg(); uu = -uu; uv = -uv; }
        if (v.isNegative()) { v.neg(); vu = -vu; vv = -vv; }

        matrixVecMul(&r, &s, uu, uv, vu, vv, &carryR, &carryS);

        uint64_t r0 = (r.d[0] * MU_P) & MSK62;
        uint64_t s0 = (s.d[0] * MU_P) & MSK62;
        r0_P.mult(&mod, r0);
        s0_P.mult(&mod, s0);
        carryR = r.addCh(&r0_P, carryR);
        carryS = s.addCh(&s0_P, carryS);

        shiftR62(u.d);
        shiftR62(v.d);
        shiftR62(r.d, carryR);
        shiftR62(s.d, carryS);
    }

    if (u.isNegative()) {
        u.neg();
        r.neg();
    }

    if (!u.isOne()) {
        memset(result, 0, 32);
        return;
    }

    while (r.isNegative()) {
        r.add(&r, &mod);
    }

    while (true) {
        uint64_t diff[5];
        unsigned char b = 0;
        b = _subborrow_u64(b, r.d[0], mod.d[0], &diff[0]);
        b = _subborrow_u64(b, r.d[1], mod.d[1], &diff[1]);
        b = _subborrow_u64(b, r.d[2], mod.d[2], &diff[2]);
        b = _subborrow_u64(b, r.d[3], mod.d[3], &diff[3]);
        b = _subborrow_u64(b, r.d[4], mod.d[4], &diff[4]);
        if (!b) {
            r.d[0] = diff[0]; r.d[1] = diff[1]; r.d[2] = diff[2]; r.d[3] = diff[3]; r.d[4] = diff[4];
        } else {
            break;
        }
    }

    result[0] = r.d[0];
    result[1] = r.d[1];
    result[2] = r.d[2];
    result[3] = r.d[3];
}

void modMulMontP(uint64_t *result, const uint64_t *a, const uint64_t *b) {
    uint64_t temp[8] = {0};

    for (int i = 0; i < 4; i++) {
        uint64_t carry = 0;
        for (int j = 0; j < 4; j++) {
            uint128_t prod =
                (uint128_t)a[i] * (uint128_t)b[j]
              + (uint128_t)temp[i + j]
              + (uint128_t)carry;

            temp[i + j] = (uint64_t)prod;
            carry = (uint64_t)(prod >> 64);
        }
        temp[i + 4] += carry;
    }

    montgomeryReduceP(result, temp + 4, temp);
}

void modSqrMontP(uint64_t *out, const uint64_t *in) {
    modMulMontP(out, in, in);
}

void scalarReduceN(uint64_t *r, const uint64_t *k) {
    bool ge = true;
    for (int i = 3; i >= 0; i--) {
        if (k[i] > N_CONST[i]) { ge = true; break; }
        if (k[i] < N_CONST[i]) { ge = false; break; }
    }

    if (ge) {
        uint128_t borrow = 0;
        for (int i = 0; i < 4; i++) {
            uint128_t res = (uint128_t)k[i] - N_CONST[i] - borrow;
            r[i] = (uint64_t)res;
            borrow = (res >> 127) & 1;
        }
    } else {
        for (int i = 0; i < 4; i++) r[i] = k[i];
    }
}

void modExpMontP(uint64_t *res, const uint64_t *base, const uint64_t *exp)
{
    uint64_t one[4] = {1,0,0,0};

    toMontgomeryP(res, one);

    uint64_t acc[4];
    for (int i = 0; i < 4; i++)
        acc[i] = base[i];

    for (int word = 3; word >= 0; word--) {
        for (int bit = 63; bit >= 0; bit--) {
            modSqrMontP(res, res);
            if ((exp[word] >> bit) & 1ULL) {
                modMulMontP(res, res, acc);
            }
        }
    }
}

void sqrtModP(uint64_t y[4], const uint64_t v[4]) {
    uint64_t exp[4];

    exp[0] = 0xFFFFFFFFBFFFFF0CULL;
    exp[1] = 0xFFFFFFFFFFFFFFFFULL;
    exp[2] = 0xFFFFFFFFFFFFFFFFULL;
    exp[3] = 0x3FFFFFFFFFFFFFFFULL;

    modExpMontP(y, v, exp);
}

void jacobianSetInfinity(ECPointJacobian *point) {
    for (int i = 0; i < 4; i++) {
        point->X[i] = ONE_MONT[i];
        point->Y[i] = ONE_MONT[i];
        point->Z[i] = 0;
    }
    point->infinity = 1;
}

bool jacobianIsInfinity(const ECPointJacobian *point) {
    uint64_t z_zero = 0;
    for (int i = 0; i < 4; i++) {
        z_zero |= point->Z[i];
    }
    return point->infinity || (z_zero == 0);
}

void jacobianToAffine(ECPointAffine *aff, const ECPointJacobian *jac) {
    if (jacobianIsInfinity(jac)) {
        for (int i = 0; i < 4; i++) aff->x[i] = aff->y[i] = 0;
        aff->infinity = 1;
        return;
    }

    uint64_t z_std[4], zInv[4], zInv2[4], zInv3[4];
    uint64_t x_std[4], y_std[4];

    fromMontgomeryP(z_std, jac->Z);
    modinvP(zInv, z_std);
    modSquareK1(zInv2, zInv);
    modMulK1(zInv3, zInv2, zInv);
    fromMontgomeryP(x_std, jac->X);
    fromMontgomeryP(y_std, jac->Y);
    modMulK1(aff->x, x_std, zInv2);
    modMulK1(aff->y, y_std, zInv3);

    aff->infinity = 0;
}

void affineToJacobian(ECPointJacobian *jac, const ECPointAffine *aff) {

    if (aff->infinity) {
        jacobianSetInfinity(jac);
        return;
    }

    toMontgomeryP(jac->X, aff->x);
    toMontgomeryP(jac->Y, aff->y);

    for (int i = 0; i < 4; i++) { jac->Z[i] = ONE_MONT[i]; }

    jac->infinity = 0;
}

void jacobianDouble(ECPointJacobian *R, const ECPointJacobian *P) {
    if (jacobianIsInfinity(P) || ((P->Y[0] | P->Y[1] | P->Y[2] | P->Y[3]) == 0)) {
        jacobianSetInfinity(R);
        return;
    }

    uint64_t A[4], B[4], C[4], D[4], E[4], F[4], newZ[4], tmp[4];

    bool z_is_one = (P->Z[0] == ONE_MONT[0] && P->Z[1] == ONE_MONT[1] && P->Z[2] == ONE_MONT[2] && P->Z[3] == ONE_MONT[3]);

    if (z_is_one) {
        for(int i = 0; i < 4; ++i) newZ[i] = P->Y[i];
    } else {
        modMulMontP(newZ, P->Y, P->Z);
    }
    modAddP(newZ, newZ, newZ);
    
    modMulMontP(A, P->X, P->X);
    modMulMontP(B, P->Y, P->Y);
    modMulMontP(C, B, B);
    modMulMontP(D, P->X, B);
    modAddP(D, D, D);
    modAddP(D, D, D);
    modAddP(E, A, A);
    modAddP(E, E, A);
    modMulMontP(F, E, E);
    modSubP(R->X, F, D);
    modSubP(R->X, R->X, D);
    modSubP(tmp, D, R->X);
    modMulMontP(R->Y, E, tmp);
    modAddP(C, C, C);
    modAddP(C, C, C);
    modAddP(C, C, C);
    modSubP(R->Y, R->Y, C);

    for(int i=0; i<4; i++) R->Z[i] = newZ[i];
    R->infinity = 0;
}

void jacobianAdd(ECPointJacobian *R, const ECPointJacobian *P, const ECPointJacobian *Q) {
    if (jacobianIsInfinity(P)) { *R = *Q; return; }
    if (jacobianIsInfinity(Q)) { *R = *P; return; }

    uint64_t Z1Z1[4], Z2Z2[4], U1[4], U2[4], S1[4], S2[4];
    uint64_t H[4], I[4], J[4], r[4], V[4], SJ[4], newZ[4], newX[4], newY[4];

    bool z1_is_one = (P->Z[0] == ONE_MONT[0] && P->Z[1] == ONE_MONT[1] && P->Z[2] == ONE_MONT[2] && P->Z[3] == ONE_MONT[3]);
    bool z2_is_one = (Q->Z[0] == ONE_MONT[0] && Q->Z[1] == ONE_MONT[1] && Q->Z[2] == ONE_MONT[2] && Q->Z[3] == ONE_MONT[3]);

    if (z1_is_one) {
        for(int i = 0; i < 4; ++i) { Z1Z1[i] = ONE_MONT[i]; U2[i] = Q->X[i]; S2[i] = Q->Y[i]; }
    } else {
        modMulMontP(Z1Z1, P->Z, P->Z);
        modMulMontP(U2, Q->X, Z1Z1);
        modMulMontP(S2, Q->Y, Z1Z1);
        modMulMontP(S2, S2, P->Z);
    }
    
    if (z2_is_one) {
        for(int i = 0; i < 4; ++i) { Z2Z2[i] = ONE_MONT[i]; U1[i] = P->X[i]; S1[i] = P->Y[i]; }
    } else {
        modMulMontP(Z2Z2, Q->Z, Q->Z);
        modMulMontP(U1, P->X, Z2Z2);
        modMulMontP(S1, P->Y, Z2Z2);
        modMulMontP(S1, S1, Q->Z);
    }

    modSubP(H, U2, U1);
    modSubP(r, S2, S1);

    if ((H[0] | H[1] | H[2] | H[3]) == 0) {
        if ((r[0] | r[1] | r[2] | r[3]) == 0) {
            jacobianDouble(R, P);
        } else {
            jacobianSetInfinity(R);
        }
        return;
    }

    modMulMontP(I, H, H);
    modMulMontP(J, I, H);
    modMulMontP(V, U1, I);
    modMulMontP(newX, r, r);
    modSubP(newX, newX, J);
    modSubP(newX, newX, V);
    modSubP(newX, newX, V);
    modSubP(newY, V, newX);
    modMulMontP(newY, newY, r);
    modMulMontP(SJ, S1, J);
    modSubP(newY, newY, SJ);

    if (z1_is_one && z2_is_one) {
        for(int i = 0; i < 4; ++i) newZ[i] = H[i];
    } else if (z1_is_one) {
        modMulMontP(newZ, Q->Z, H);
    } else if (z2_is_one) {
        modMulMontP(newZ, P->Z, H);
    } else {
        modMulMontP(newZ, P->Z, Q->Z);
        modMulMontP(newZ, newZ, H);
    }

    for(int i=0; i<4; i++) {
        R->X[i] = newX[i];
        R->Y[i] = newY[i];
        R->Z[i] = newZ[i];
    }
    R->infinity = 0;
}

void endomorphismMap(ECPointJacobian *R, const ECPointJacobian *P) {
    uint64_t beta_mont[4];
    toMontgomeryP(beta_mont, BETA_P);
    modMulMontP(R->X, P->X, beta_mont);

    if (R != P) {
        for (int i = 0; i < 4; i++) {
            R->Y[i] = P->Y[i];
            R->Z[i] = P->Z[i];
        }
        R->infinity = P->infinity;
    }
}

void initPreCompG(int windowSize) {
    int dnorm = (128 + windowSize - 1) / windowSize;
    int dphi = (128 + windowSize - 1) / windowSize;
    int tableSize = (1 << windowSize) - 1;

    uint64_t gxMont[4], gyMont[4];
    toMontgomeryP(gxMont, GX_CONST);
    toMontgomeryP(gyMont, GY_CONST);

    ECPointJacobian g, ge;
    for (int i = 0; i < 4; i++) {
        g.X[i]  = gxMont[i];
        g.Y[i]  = gyMont[i];
        g.Z[i]  = ONE_MONT[i];
        ge.X[i] = gxMont[i];
        ge.Y[i] = gyMont[i];
        ge.Z[i] = ONE_MONT[i];
    }
    g.infinity  = 0;
    ge.infinity = 0;

    endomorphismMap(&ge, &ge);

    jacNorm[0]  = g;
    jacEndo[0] = ge;

    for (int j = 1; j < windowSize; j++) {
        jacNorm[j]  = jacNorm[j - 1];
        jacEndo[j]  = jacEndo[j - 1];

        for (int i = 0; i < dnorm; i++)
            jacobianDouble(&jacNorm[j], &jacNorm[j]);

        for (int i = 0; i < dphi; i++)
            jacobianDouble(&jacEndo[j], &jacEndo[j]);
    }

    for (int i = 1; i <= tableSize; i++) {
        jacobianSetInfinity(&preCompG[i - 1]);

        for (int j = 0; j < windowSize; j++) {
            if ((i >> j) & 1) {
                ECPointJacobian tmp;
                jacobianAdd(&tmp, &preCompG[i - 1], &jacNorm[j]);
                preCompG[i - 1] = tmp;
            }
        }

        jacobianSetInfinity(&preCompGphi[i - 1]);

        for (int j = 0; j < windowSize; j++) {
            if ((i >> j) & 1) {
                ECPointJacobian tmp;
                jacobianAdd(&tmp, &preCompGphi[i - 1], &jacEndo[j]);
                preCompGphi[i - 1] = tmp;
            }
        }
    }
}

void initPreCompH(const ECPointJacobian *h, int windowSize) {
    int dnorm = (128 + windowSize - 1) / windowSize;
    int dphi = (128 + windowSize - 1) / windowSize;
    int tableSize = (1 << windowSize) - 1;

    ECPointJacobian H, HE;
    for (int i = 0; i < 4; i++) {
        H.X[i]  = h->X[i];
        H.Y[i]  = h->Y[i];
        H.Z[i]  = h->Z[i];
        HE.X[i] = h->X[i];
        HE.Y[i] = h->Y[i];
        HE.Z[i] = h->Z[i];
    }
    H.infinity  = h->infinity;
    HE.infinity = h->infinity;

    endomorphismMap(&HE, &HE);

    jacNormH[0] = H;
    jacEndoH[0] = HE;

    for (int j = 1; j < windowSize; j++) {
        jacNormH[j]  = jacNormH[j - 1];
        jacEndoH[j]  = jacEndoH[j - 1];

        for (int i = 0; i < dnorm; i++)
            jacobianDouble(&jacNormH[j], &jacNormH[j]);

        for (int i = 0; i < dphi; i++)
            jacobianDouble(&jacEndoH[j], &jacEndoH[j]);
    }

    for (int i = 1; i <= tableSize; i++) {
        jacobianSetInfinity(&preCompH[i - 1]);

        for (int j = 0; j < windowSize; j++) {
            if ((i >> j) & 1) {
                ECPointJacobian tmp;
                jacobianAdd(&tmp, &preCompH[i - 1], &jacNormH[j]);
                preCompH[i - 1] = tmp;
            }
        }

        jacobianSetInfinity(&preCompHphi[i - 1]);

        for (int j = 0; j < windowSize; j++) {
            if ((i >> j) & 1) {
                ECPointJacobian tmp;
                jacobianAdd(&tmp, &preCompHphi[i - 1], &jacEndoH[j]);
                preCompHphi[i - 1] = tmp;
            }
        }
    }
}

void scalarMulShiftVar(uint64_t r[4], const uint64_t a[4], const uint64_t b[4], int shift) {
    uint64_t t[8] = {0};

    for (int i = 0; i < 4; i++) {
        uint64_t prod[5];
        imm_umul(b, a[i], prod);

        unsigned char c = 0;
        c = _addcarry_u64(c, t[i + 0], prod[0], &t[i + 0]);
        c = _addcarry_u64(c, t[i + 1], prod[1], &t[i + 1]);
        c = _addcarry_u64(c, t[i + 2], prod[2], &t[i + 2]);
        c = _addcarry_u64(c, t[i + 3], prod[3], &t[i + 3]);
        c = _addcarry_u64(c, t[i + 4], prod[4], &t[i + 4]);

        for (int k = i + 5; k < 8 && c; k++) {
            c = _addcarry_u64(c, t[k], 0, &t[k]);
        }
    }

    int word_shift = shift / 64;
    int bit_shift  = shift % 64;
    for (int i = 0; i < 4; i++) {
        uint128_t v = 0;
        if (i + word_shift < 8) v = uint128_t(t[i + word_shift]);
        if (bit_shift && i + word_shift + 1 < 8) {
            v += uint128_t(t[i + word_shift + 1]) << 64;
        }
        r[i] = (uint64_t)(v >> bit_shift);
    }
}

void scalarMul(uint64_t r[4], const uint64_t a[4], const uint64_t b[4]) {
    uint64_t t[8] = {0};

    for (int i = 0; i < 4; i++) {
        uint64_t prod[5];
        imm_umul(b, a[i], prod);

        unsigned char c = 0;
        c = _addcarry_u64(c, t[i + 0], prod[0], &t[i + 0]);
        c = _addcarry_u64(c, t[i + 1], prod[1], &t[i + 1]);
        c = _addcarry_u64(c, t[i + 2], prod[2], &t[i + 2]);
        c = _addcarry_u64(c, t[i + 3], prod[3], &t[i + 3]);
        c = _addcarry_u64(c, t[i + 4], prod[4], &t[i + 4]);

        for (int k = i + 5; k < 8 && c; k++) {
            c = _addcarry_u64(c, t[k], 0, &t[k]);
        }
    }

    static const uint64_t N_C[4] = {
        0x402DA1732FC9BEBFULL,
        0x4551231950B75FC4ULL,
        0x0000000000000001ULL,
        0x0000000000000000ULL
    };

    uint64_t S[8] = {t[0], t[1], t[2], t[3], 0, 0, 0, 0};
    for (int i = 0; i < 4; i++) {
        uint64_t prod[5];
        imm_umul(N_C, t[i + 4], prod);
        unsigned char c = 0;
        c = _addcarry_u64(c, S[i + 0], prod[0], &S[i + 0]);
        c = _addcarry_u64(c, S[i + 1], prod[1], &S[i + 1]);
        c = _addcarry_u64(c, S[i + 2], prod[2], &S[i + 2]);
        c = _addcarry_u64(c, S[i + 3], prod[3], &S[i + 3]);
        c = _addcarry_u64(c, S[i + 4], prod[4], &S[i + 4]);
        for (int k = i + 5; k < 8 && c; k++) {
            c = _addcarry_u64(c, S[k], 0, &S[k]);
        }
    }

    uint64_t S_high[3] = {S[4], S[5], S[6]};
    uint64_t S2[8] = {S[0], S[1], S[2], S[3], 0, 0, 0, 0};
    for (int i = 0; i < 3; i++) {
        uint64_t prod[5];
        imm_umul(N_C, S_high[i], prod);
        unsigned char c = 0;
        c = _addcarry_u64(c, S2[i + 0], prod[0], &S2[i + 0]);
        c = _addcarry_u64(c, S2[i + 1], prod[1], &S2[i + 1]);
        c = _addcarry_u64(c, S2[i + 2], prod[2], &S2[i + 2]);
        c = _addcarry_u64(c, S2[i + 3], prod[3], &S2[i + 3]);
        c = _addcarry_u64(c, S2[i + 4], prod[4], &S2[i + 4]);
        for (int k = i + 5; k < 8 && c; k++) {
            c = _addcarry_u64(c, S2[k], 0, &S2[k]);
        }
    }

    uint64_t h2 = S2[4];
    uint64_t res[4] = {S2[0], S2[1], S2[2], S2[3]};
    unsigned char top_carry = 0;
    if (h2 != 0) {
        uint64_t prod[5];
        imm_umul(N_C, h2, prod);
        unsigned char c = 0;
        c = _addcarry_u64(c, res[0], prod[0], &res[0]);
        c = _addcarry_u64(c, res[1], prod[1], &res[1]);
        c = _addcarry_u64(c, res[2], prod[2], &res[2]);
        c = _addcarry_u64(c, res[3], prod[3], &res[3]);
        top_carry = c;
    }

    if (top_carry) {
        unsigned char c = 0;
        c = _addcarry_u64(c, res[0], N_C[0], &res[0]);
        c = _addcarry_u64(c, res[1], N_C[1], &res[1]);
        c = _addcarry_u64(c, res[2], N_C[2], &res[2]);
        _addcarry_u64(c, res[3], N_C[3], &res[3]);
    }

    while (true) {
        uint64_t diff[4];
        unsigned char b = 0;
        b = _subborrow_u64(b, res[0], N_CONST[0], &diff[0]);
        b = _subborrow_u64(b, res[1], N_CONST[1], &diff[1]);
        b = _subborrow_u64(b, res[2], N_CONST[2], &diff[2]);
        b = _subborrow_u64(b, res[3], N_CONST[3], &diff[3]);
        if (!b) {
            res[0] = diff[0]; res[1] = diff[1]; res[2] = diff[2]; res[3] = diff[3];
        } else {
            break;
        }
    }

    for (int i = 0; i < 4; i++) r[i] = res[i];
}

void scalarAdd(uint64_t r[4], const uint64_t a[4], const uint64_t b[4]) {
    uint64_t sum[4];
    unsigned char c = 0;
    c = _addcarry_u64(c, a[0], b[0], &sum[0]);
    c = _addcarry_u64(c, a[1], b[1], &sum[1]);
    c = _addcarry_u64(c, a[2], b[2], &sum[2]);
    c = _addcarry_u64(c, a[3], b[3], &sum[3]);

    uint64_t diff[4];
    unsigned char borrow = 0;
    borrow = _subborrow_u64(borrow, sum[0], N_CONST[0], &diff[0]);
    borrow = _subborrow_u64(borrow, sum[1], N_CONST[1], &diff[1]);
    borrow = _subborrow_u64(borrow, sum[2], N_CONST[2], &diff[2]);
    borrow = _subborrow_u64(borrow, sum[3], N_CONST[3], &diff[3]);

    if (c || borrow == 0) {
        r[0] = diff[0]; r[1] = diff[1]; r[2] = diff[2]; r[3] = diff[3];
    } else {
        r[0] = sum[0]; r[1] = sum[1]; r[2] = sum[2]; r[3] = sum[3];
    }
}

void scalarSub(uint64_t r[4], const uint64_t a[4], const uint64_t b[4]) {
    uint64_t diff[4];
    unsigned char borrow = 0;
    borrow = _subborrow_u64(borrow, a[0], b[0], &diff[0]);
    borrow = _subborrow_u64(borrow, a[1], b[1], &diff[1]);
    borrow = _subborrow_u64(borrow, a[2], b[2], &diff[2]);
    borrow = _subborrow_u64(borrow, a[3], b[3], &diff[3]);

    if (borrow) {
        unsigned char c = 0;
        c = _addcarry_u64(c, diff[0], N_CONST[0], &r[0]);
        c = _addcarry_u64(c, diff[1], N_CONST[1], &r[1]);
        c = _addcarry_u64(c, diff[2], N_CONST[2], &r[2]);
        c = _addcarry_u64(c, diff[3], N_CONST[3], &r[3]);
    } else {
        r[0] = diff[0]; r[1] = diff[1]; r[2] = diff[2]; r[3] = diff[3];
    }
}

int scalarIsZero(const uint64_t a[4]) {
    return ((a[0] | a[1] | a[2] | a[3]) == 0ULL);
}

void scalarNeg(uint64_t r[4], const uint64_t a[4]) {
    if (scalarIsZero(a)) {
        r[0] = r[1] = r[2] = r[3] = 0;
        return;
    }

    unsigned char borrow = 0;
    borrow = _subborrow_u64(borrow, N_CONST[0], a[0], &r[0]);
    borrow = _subborrow_u64(borrow, N_CONST[1], a[1], &r[1]);
    borrow = _subborrow_u64(borrow, N_CONST[2], a[2], &r[2]);
    borrow = _subborrow_u64(borrow, N_CONST[3], a[3], &r[3]);
}

void scalarSplitLambda(uint64_t r1[4], uint64_t r2[4], const uint64_t k[4]) {
    uint64_t c1[4], c2[4], t1[4], t2[4];

    scalarMulShiftVar(c1, k, G1, 384);
    scalarMulShiftVar(c2, k, G2, 384);
    scalarMul(t1, c1, MINUS_B1);
    scalarMul(t2, c2, MINUS_B2);
    scalarAdd(r2, t1, t2);
    scalarMul(t1, r2, LAMBDA_N);
    scalarNeg(t1, t1);
    scalarAdd(r1, t1, k);
}

void jacobianScalarMultPhi(ECPointJacobian *result, ECPointJacobian *preCompTable, ECPointJacobian *preCompTablePhi, const uint64_t *scalar, int windowSize) {
    uint64_t r1[4], r2[4];
    int d = (128 + windowSize - 1) / windowSize;

    scalarSplitLambda(r1, r2, scalar);
    jacobianSetInfinity(result);

    for (int col = d - 1; col >= 0; col--) {
        if (!jacobianIsInfinity(result)) {
            jacobianDouble(result, result);
        }

        int idx1 = 0;
        int idx2 = 0;

        for (int row = 0; row < windowSize; row++) {
            int bitIndex = row * d + col;
            if (bitIndex >= 128) continue;
            int limb = bitIndex / 64;
            int shift = bitIndex % 64;
            uint64_t bit1 = (r1[limb] >> shift) & 1ULL;
            uint64_t bit2 = (r2[limb] >> shift) & 1ULL;
            idx1 |= (bit1 << row);
            idx2 |= (bit2 << row);
        }

        if (idx1 != 0) {
            ECPointJacobian tmp;
            jacobianAdd(&tmp, result, &preCompTable[idx1 - 1]);
            *result = tmp;
        }

        if (idx2 != 0) {
            ECPointJacobian tmp;
            jacobianAdd(&tmp, result, &preCompTablePhi[idx2 - 1]);
            *result = tmp;
        }
    }
}

void pointAddJacobian(ECPointJacobian *R, const ECPointJacobian *P, const ECPointJacobian *Q) { jacobianAdd(R, P, Q); }
void pointDoubleJacobian(ECPointJacobian *R, const ECPointJacobian *P) { jacobianDouble(R, P); }
void serializePublicKey(unsigned char *out, const ECPointAffine *publicKey) {
    unsigned char prefix = (publicKey->y[0] & 1ULL) ? 0x03 : 0x02;
    out[0] = prefix;

    for (int i = 0; i < 4; i++) {
        uint64_t word = publicKey->x[3 - i];
        out[1 + i*8 + 0] = (word >> 56) & 0xFF;
        out[1 + i*8 + 1] = (word >> 48) & 0xFF;
        out[1 + i*8 + 2] = (word >> 40) & 0xFF;
        out[1 + i*8 + 3] = (word >> 32) & 0xFF;
        out[1 + i*8 + 4] = (word >> 24) & 0xFF;
        out[1 + i*8 + 5] = (word >> 16) & 0xFF;
        out[1 + i*8 + 6] = (word >> 8) & 0xFF;
        out[1 + i*8 + 7] = word & 0xFF;
    }
}

void generatePublicKey(ECPointJacobian *preCompTable, ECPointJacobian *preCompTablePhi, unsigned char *out, const uint64_t *PRIV_KEY, int windowSize) {
    ECPointAffine pub;
    ECPointJacobian pub_jac;

    jacobianScalarMultPhi(&pub_jac, preCompTable, preCompTablePhi, PRIV_KEY, windowSize);
    jacobianToAffine(&pub, &pub_jac);
    serializePublicKey(out, &pub);
}

bool decompressAndValidatePublicKey(ECPointAffine* out, const unsigned char compressed[33]) {
    unsigned char prefix = compressed[0];
    if (prefix != 0x02 && prefix != 0x03) {
        return false;
    }

    for (int i = 0; i < 4; i++) {
        uint64_t word = 0;
        for (int j = 0; j < 8; j++) {
            word = (word << 8) | compressed[1 + i * 8 + j];
        }
        out->x[3 - i] = word;
    }

    bool is_zero = true;
    for (int i = 0; i < 4; i++) {
        if (out->x[i] != 0) { is_zero = false; break; }
    }
    if (is_zero) return false;

    for (int i = 3; i >= 0; i--) {
        if (out->x[i] > P_CONST[i]) return false;
        if (out->x[i] < P_CONST[i]) break;
        if (i == 0) return false;
    }

    uint64_t x_mont[4], x2[4], x3[4], rhs[4];

    toMontgomeryP(x_mont, out->x);
    modMulMontP(x2, x_mont, x_mont);
    modMulMontP(x3, x2, x_mont);
    modAddP(rhs, x3, SEVEN_MONT);

    uint64_t y_mont[4];
    sqrtModP(y_mont, rhs);

    uint64_t y2_test[4];
    modMulMontP(y2_test, y_mont, y_mont);
    for (int i = 0; i < 4; i++) {
        if (y2_test[i] != rhs[i]) {
            return false;
        }
    }

    fromMontgomeryP(out->y, y_mont);

    uint64_t is_odd   = out->y[0] & 1ULL;
    uint64_t want_odd = (prefix == 0x03);

    if (is_odd != want_odd) {
        modSubP(out->y, P_CONST, out->y);
    }

    out->infinity = 0;
    return true;
}

void decompressPublicKey(ECPointAffine* out, const unsigned char compressed[33]) {
    decompressAndValidatePublicKey(out, compressed);
}

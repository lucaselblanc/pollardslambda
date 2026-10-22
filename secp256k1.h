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

#ifndef EC_SECP256K1_H
#define EC_SECP256K1_H

#include <stdint.h>
#include <x86intrin.h>

#ifndef _addcarry_u64
#define _addcarry_u64(c_in, a, b, out) __builtin_ia32_addcarryx_u64(c_in, a, b, (unsigned long long*)(out))
#endif

#ifndef _subborrow_u64
#define _subborrow_u64(b_in, a, b, out) __builtin_ia32_sbb_u64(b_in, a, b, (unsigned long long*)(out))
#endif

static inline uint64_t _umul128(uint64_t a, uint64_t b, uint64_t *h) {
    unsigned __int128 p = (unsigned __int128)a * b;
    *h = (uint64_t)(p >> 64);
    return (uint64_t)p;
}

#include "parallel_hashmap/phmap.h"
#include <condition_variable>
#include <openssl/sha.h>
#include <fstream>
#include <unistd.h>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <random>
#include <thread>
#include <atomic>
#include <mutex>
#include <chrono>
#include <limits>
#include <climits>
#include <ctime>
#include <cmath>
#include <cstring>
#include <cstdint>
#include <tuple>
#include <algorithm>
#include <cerrno>
#include <cstdio>
#include <fcntl.h>
#include <sys/stat.h>
#include <sys/types.h>

struct uint256_t {
    uint64_t limbs[4];
};

typedef struct {
    uint64_t x[4];
    uint64_t y[4];
    int infinity;
} ECPointAffine;

typedef struct {
    uint64_t X[4];
    uint64_t Y[4];
    uint64_t Z[4];
    int infinity;
} ECPointJacobian;

using uint128_t = unsigned __int128;

extern ECPointJacobian* preCompG;
extern ECPointJacobian* preCompGphi;
extern ECPointJacobian* preCompH;
extern ECPointJacobian* preCompHphi;
extern ECPointJacobian* jacNorm;
extern ECPointJacobian* jacEndo;
extern ECPointJacobian* jacNormH;
extern ECPointJacobian* jacEndoH;

uint256_t modinv(uint256_t base, uint256_t mod);

void affineToJacobian(ECPointJacobian *jac, const ECPointAffine *aff);
void decompressPublicKey(ECPointAffine* out, const unsigned char compressed[33]);
void endomorphismMap(ECPointJacobian *R, const ECPointJacobian *P);
void fromMontgomeryP(uint64_t *result, const uint64_t *a);
void generatePublicKey(ECPointJacobian *preCompTable, ECPointJacobian *preCompTablePhi, unsigned char *out, const uint64_t *PRIV_KEY, int windowSize);
void initPreCompG(int windowSize);
void initPreCompH(const ECPointJacobian *h, int windowSize);
void jacobianScalarMultPhi(ECPointJacobian *result, ECPointJacobian *preCompTable, ECPointJacobian *preCompTablePhi, const uint64_t *scalar, int windowSize);
void jacobianDouble(ECPointJacobian *R, const ECPointJacobian *P);
void jacobianAdd(ECPointJacobian *R, const ECPointJacobian *P, const ECPointJacobian *Q);
void jacobianToAffine(ECPointAffine *aff, const ECPointJacobian *jac);
void jacobianSetInfinity(ECPointJacobian *point);
bool jacobianIsInfinity(const ECPointJacobian *P);
void modMulMontP(uint64_t *result, const uint64_t *a, const uint64_t *b);
void modMulK1(uint64_t *result, const uint64_t *a, const uint64_t *b);
void modSquareK1(uint64_t *result, const uint64_t *a);
void modinvP(uint64_t *result, const uint64_t *a);
static inline void modAddP(uint64_t *result, const uint64_t *a, const uint64_t *b) {
    uint64_t sum[4];
    unsigned char c = 0;
    c = _addcarry_u64(c, a[0], b[0], &sum[0]);
    c = _addcarry_u64(c, a[1], b[1], &sum[1]);
    c = _addcarry_u64(c, a[2], b[2], &sum[2]);
    c = _addcarry_u64(c, a[3], b[3], &sum[3]);

    uint64_t diff[4];
    unsigned char borrow = 0;
    borrow = _subborrow_u64(borrow, sum[0], 0xFFFFFFFEFFFFFC2FULL, &diff[0]);
    borrow = _subborrow_u64(borrow, sum[1], 0xFFFFFFFFFFFFFFFFULL, &diff[1]);
    borrow = _subborrow_u64(borrow, sum[2], 0xFFFFFFFFFFFFFFFFULL, &diff[2]);
    borrow = _subborrow_u64(borrow, sum[3], 0xFFFFFFFFFFFFFFFFULL, &diff[3]);

    uint64_t mask = (c || borrow == 0);
    result[0] = mask ? diff[0] : sum[0];
    result[1] = mask ? diff[1] : sum[1];
    result[2] = mask ? diff[2] : sum[2];
    result[3] = mask ? diff[3] : sum[3];
}

static inline void modSubP(uint64_t *result, const uint64_t *a, const uint64_t *b) {
    uint64_t diff[4];
    unsigned char borrow = 0;
    borrow = _subborrow_u64(borrow, a[0], b[0], &diff[0]);
    borrow = _subborrow_u64(borrow, a[1], b[1], &diff[1]);
    borrow = _subborrow_u64(borrow, a[2], b[2], &diff[2]);
    borrow = _subborrow_u64(borrow, a[3], b[3], &diff[3]);

    uint64_t mask = (0ULL - (uint64_t)borrow) & 0x1000003D1ULL;
    borrow = _subborrow_u64(0, diff[0], mask, &result[0]);
    borrow = _subborrow_u64(borrow, diff[1], 0, &result[1]);
    borrow = _subborrow_u64(borrow, diff[2], 0, &result[2]);
    _subborrow_u64(borrow, diff[3], 0, &result[3]);
}
void modExpMontP(uint64_t *res, const uint64_t *base, const uint64_t *exp);
void pointInitJacobian(ECPointJacobian *P);
void pointAddJacobian(ECPointJacobian *R, const ECPointJacobian *P, const ECPointJacobian *Q);
void pointDoubleJacobian(ECPointJacobian *R, const ECPointJacobian *P);
void scalarReduceN(uint64_t *r, const uint64_t *k);
void scalarMul(uint64_t r[4], const uint64_t a[4], const uint64_t b[4]);
void scalarAdd(uint64_t r[4], const uint64_t a[4], const uint64_t b[4]);
void scalarSub(uint64_t r[4], const uint64_t a[4], const uint64_t b[4]);
void scalarNeg(uint64_t r[4], const uint64_t a[4]);
int scalarIsZero(const uint64_t a[4]);
void serializePublicKey(unsigned char *out, const ECPointAffine *publicKey);
void toMontgomeryP(uint64_t *result, const uint64_t *a);

#endif /* EC_SECP256K1_H */
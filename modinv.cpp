/*********************************************************************************************************
 * This file is part of the Almost Inverse distribution: (https://github.com/lucaselblanc/almostinverse) *
 * Copyright (c) 2024, 2026 Lucas Leblanc.                                                               *
 * Distributed under the MIT software license, see the accompanying.                                     *
 * file COPYING or https://www.opensource.org/licenses/mit-license.php.                                  *
 *********************************************************************************************************/

/*******************************************************
 * High-Performance 256-bit Modular Inversion (Pornin) *
 * Based on Thomas Pornin's DivStep62 variable time    *
 * Optimized x86-64 SIMD Carry-Chain Implementation    *
 *******************************************************/

#include "secp256k1.h"

#ifndef _addcarry_u64
#define _addcarry_u64(c_in, a, b, out) __builtin_ia32_addcarryx_u64(c_in, a, b, (unsigned long long*)(out))
#endif

#ifndef _subborrow_u64
#define _subborrow_u64(b_in, a, b, out) __builtin_ia32_sbb_u64(b_in, a, b, (unsigned long long*)(out))
#endif

static inline uint64_t _umul128_inv(uint64_t a, uint64_t b, uint64_t *h) {
    unsigned __int128 p = (unsigned __int128)a * b;
    *h = (uint64_t)(p >> 64);
    return (uint64_t)p;
}

static inline int64_t _mul128_inv(int64_t a, int64_t b, int64_t *h) {
    __int128 p = (__int128)a * b;
    *h = (int64_t)(p >> 64);
    return (int64_t)p;
}

#define __shiftright128_inv(a,b,n) (((a)>>(n))|((b)<<(64-(n))))
#define __shiftleft128_inv(a,b,n) (((b)<<(n))|((a)>>(64-(n))))
#define SWAP_INV(tmp,x,y) do { tmp = x; x = y; y = tmp; } while(0)
#define MSK62_INV 0x3FFFFFFFFFFFFFFFULL

struct Int5Inv {
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

    inline void add(const Int5Inv *a, const Int5Inv *b) {
        unsigned char c = 0;
        c = _addcarry_u64(c, a->d[0], b->d[0], &d[0]);
        c = _addcarry_u64(c, a->d[1], b->d[1], &d[1]);
        c = _addcarry_u64(c, a->d[2], b->d[2], &d[2]);
        c = _addcarry_u64(c, a->d[3], b->d[3], &d[3]);
        _addcarry_u64(c, a->d[4], b->d[4], &d[4]);
    }

    inline uint64_t addCh(const Int5Inv *a, uint64_t ca, const Int5Inv *b, uint64_t cb) {
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

    inline uint64_t addCh(const Int5Inv *a, uint64_t ca) {
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

    inline uint64_t mult(const Int5Inv *a, uint64_t b) {
        unsigned char c = 0;
        uint64_t h, carry;
        d[0] = _umul128_inv(a->d[0], b, &h); carry = h;
        c = _addcarry_u64(c, _umul128_inv(a->d[1], b, &h), carry, &d[1]); carry = h;
        c = _addcarry_u64(c, _umul128_inv(a->d[2], b, &h), carry, &d[2]); carry = h;
        c = _addcarry_u64(c, _umul128_inv(a->d[3], b, &h), carry, &d[3]); carry = h;
        c = _addcarry_u64(c, _umul128_inv(a->d[4], b, &h), carry, &d[4]); carry = h;
        return carry;
    }

    inline uint64_t iMult(const Int5Inv *a, int64_t b) {
        uint64_t carry;
        Int5Inv tmp;
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
        d[0] = _umul128_inv(tmp.d[0], (uint64_t)b, &h); carry = h;
        c = _addcarry_u64(c, _umul128_inv(tmp.d[1], (uint64_t)b, &h), carry, &d[1]); carry = h;
        c = _addcarry_u64(c, _umul128_inv(tmp.d[2], (uint64_t)b, &h), carry, &d[2]); carry = h;
        c = _addcarry_u64(c, _umul128_inv(tmp.d[3], (uint64_t)b, &h), carry, &d[3]); carry = h;
        int64_t h_signed;
        c = _addcarry_u64(c, (uint64_t)_mul128_inv((int64_t)tmp.d[4], b, &h_signed), carry, &d[4]);
        carry = (uint64_t)h_signed;
        return carry;
    }
};

static inline void shiftR62_inv(uint64_t *d) {
    d[0] = __shiftright128_inv(d[0], d[1], 62);
    d[1] = __shiftright128_inv(d[1], d[2], 62);
    d[2] = __shiftright128_inv(d[2], d[3], 62);
    d[3] = __shiftright128_inv(d[3], d[4], 62);
    d[4] = ((int64_t)d[4]) >> 62;
}

static inline void shiftR62_inv(uint64_t *d, uint64_t h) {
    d[0] = __shiftright128_inv(d[0], d[1], 62);
    d[1] = __shiftright128_inv(d[1], d[2], 62);
    d[2] = __shiftright128_inv(d[2], d[3], 62);
    d[3] = __shiftright128_inv(d[3], d[4], 62);
    d[4] = __shiftright128_inv(d[4], h, 62);
}

static inline void matrixVecMul_inv(Int5Inv *u, Int5Inv *v, int64_t _11, int64_t _12, int64_t _21, int64_t _22) {
    Int5Inv t1, t2, t3, t4;
    t1.iMult(u, _11);
    t2.iMult(v, _12);
    t3.iMult(u, _21);
    t4.iMult(v, _22);
    u->add(&t1, &t2);
    v->add(&t3, &t4);
}

static inline void matrixVecMul_inv(Int5Inv *u, Int5Inv *v, int64_t _11, int64_t _12, int64_t _21, int64_t _22, uint64_t *cu, uint64_t *cv) {
    Int5Inv t1, t2, t3, t4;
    uint64_t c1 = t1.iMult(u, _11);
    uint64_t c2 = t2.iMult(v, _12);
    uint64_t c3 = t3.iMult(u, _21);
    uint64_t c4 = t4.iMult(v, _22);
    *cu = u->addCh(&t1, c1, &t2, c2);
    *cv = v->addCh(&t3, c3, &t4, c4);
}

static inline void divStep62_inv(Int5Inv *u, Int5Inv *v, int *pos, int64_t *uu, int64_t *uv, int64_t *vu, int64_t *vv) {
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
            uh = __shiftleft128_inv(u->d[*pos - 1], u->d[*pos], (uint8_t)s);
            vh = __shiftleft128_inv(v->d[*pos - 1], v->d[*pos], (uint8_t)s);
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
            SWAP_INV(w, uh, vh);
            SWAP_INV(x, u0, v0);
            SWAP_INV(_t, _u, _v);
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

static inline uint64_t compute_mu64(uint64_t m0) {
    uint64_t x = 1;
    for (int i = 0; i < 6; i++) {
        x = x * (2ULL - m0 * x);
    }
    return -x;
}

uint256_t modinv(uint256_t base, uint256_t mod) {
    if ((base.limbs[0] | base.limbs[1] | base.limbs[2] | base.limbs[3]) == 0) return uint256_t{};
    if ((mod.limbs[0] | mod.limbs[1] | mod.limbs[2] | mod.limbs[3]) == 0) return uint256_t{};
    if ((mod.limbs[0] & 1ULL) == 0) return uint256_t{};

    uint64_t mu = compute_mu64(mod.limbs[0]);

    Int5Inv mod5;
    mod5.d[0] = mod.limbs[0]; mod5.d[1] = mod.limbs[1];
    mod5.d[2] = mod.limbs[2]; mod5.d[3] = mod.limbs[3];
    mod5.d[4] = 0;

    Int5Inv u = mod5;
    Int5Inv v;
    v.d[0] = base.limbs[0]; v.d[1] = base.limbs[1];
    v.d[2] = base.limbs[2]; v.d[3] = base.limbs[3];
    v.d[4] = 0;

    Int5Inv r; r.clear();
    Int5Inv s; s.clear(); s.d[0] = 1;

    int64_t uu, uv, vu, vv;
    uint64_t carryR, carryS;
    int pos = 4;
    while (pos >= 1 && (u.d[pos] | v.d[pos]) == 0) pos--;

    Int5Inv r0_P, s0_P;

    while (!v.isZero()) {
        divStep62_inv(&u, &v, &pos, &uu, &uv, &vu, &vv);

        matrixVecMul_inv(&u, &v, uu, uv, vu, vv);
        if (u.isNegative()) { u.neg(); uu = -uu; uv = -uv; }
        if (v.isNegative()) { v.neg(); vu = -vu; vv = -vv; }

        matrixVecMul_inv(&r, &s, uu, uv, vu, vv, &carryR, &carryS);

        uint64_t r0 = (r.d[0] * mu) & MSK62_INV;
        uint64_t s0 = (s.d[0] * mu) & MSK62_INV;
        r0_P.mult(&mod5, r0);
        s0_P.mult(&mod5, s0);
        carryR = r.addCh(&r0_P, carryR);
        carryS = s.addCh(&s0_P, carryS);

        shiftR62_inv(u.d);
        shiftR62_inv(v.d);
        shiftR62_inv(r.d, carryR);
        shiftR62_inv(s.d, carryS);
    }

    if (u.isNegative()) {
        u.neg();
        r.neg();
    }

    if (!u.isOne()) {
        return uint256_t{};
    }

    while (r.isNegative()) {
        r.add(&r, &mod5);
    }

    while (true) {
        uint64_t diff[5];
        unsigned char b = 0;
        b = _subborrow_u64(b, r.d[0], mod5.d[0], &diff[0]);
        b = _subborrow_u64(b, r.d[1], mod5.d[1], &diff[1]);
        b = _subborrow_u64(b, r.d[2], mod5.d[2], &diff[2]);
        b = _subborrow_u64(b, r.d[3], mod5.d[3], &diff[3]);
        b = _subborrow_u64(b, r.d[4], mod5.d[4], &diff[4]);
        if (!b) {
            r.d[0] = diff[0]; r.d[1] = diff[1]; r.d[2] = diff[2]; r.d[3] = diff[3]; r.d[4] = diff[4];
        } else {
            break;
        }
    }

    uint256_t res;
    res.limbs[0] = r.d[0];
    res.limbs[1] = r.d[1];
    res.limbs[2] = r.d[2];
    res.limbs[3] = r.d[3];
    return res;
}
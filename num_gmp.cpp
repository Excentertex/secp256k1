#include <assert.h>
#include <string>
#include <string.h>
#include <stdlib.h>
#include <gmp.h>

#include "num.h"

typedef struct {
    int initialized;
    gmp_randstate_t rng;
} secp256k1_num_state_t;

static secp256k1_num_state_t secp256k1_num_state = {};

void secp256k1_num_start(void) {
    if (secp256k1_num_state.initialized)
        return;
    secp256k1_num_state.initialized = 1;
    gmp_randinit_default(rng);
}

void secp256k1_num_init(secp256k1_num_t *r) {
    mpz_init(r->bn);
}

void secp256k1_num_free(secp256k1_num_t *r) {
    mpz_clear(r->bn);
}

void secp256k1_num_copy(secp256k1_num_t *r, const secp256k1_num_t *a) {
    mpz_set(r->bn, a->bn);
}

void secp256k1_num_get_bin(unsigned char *r, unsigned int rlen, const secp256k1_num_t *a);
void secp256k1_num_set_bin(secp256k1_num_t *r, unsigned char *a, unsigned int alen);
void secp256k1_num_set_int(secp256k1_num_t *r, int a);
void secp256k1_num_mod_inverse(secp256k1_num_t *r, const secp256k1_num_t *a, const secp256k1_num_t *m);
void secp256k1_num_mod_mul(secp256k1_num_t *r, const secp256k1_num_t *a, const secp256k1_num_t *b, const secp256k1_num_t *m);
int  secp256k1_num_cmp(const secp256k1_num_t *a, const secp256k1_num_t *b);
void secp256k1_num_add(secp256k1_num_t *r, const secp256k1_t *a, const secp256k1_t *b);
void secp256k1_num_sub(secp256k1_num_t *r, const secp256k1_t *a, const secp256k1_t *b);
void secp256k1_num_mul(secp256k1_num_t *r, const secp256k1_t *a, const secp256k1_t *b);
void secp256k1_num_div(secp256k1_num_t *r, const secp256k1_t *a, const secp256k1_t *b);
void secp256k1_num_mod(secp256k1_num_t *r, const secp256k1_t *a, const secp256k1_t *b);
int  secp256k1_num_bits(const secp256k1_num_t *a);
int  secp256k1_num_shift(secp256k1_num_t *r, int bits);
int  secp256k1_num_is_zero(const secp256k1_num_t *a);
int  secp256k1_num_is_odd(const secp256k1_num_t *a);
int  secp256k1_num_is_neg(const secp256k1_num_t *a);
int  secp256k1_num_get_bit(const secp256k1_num_t *a, int pos);
void secp256k1_num_inc(secp256k1_num_t *r);
void secp256k1_num_set_hex(secp256k1_num_t *r, const unsigned char *a, int alen);
void secp256k1_num_get_hex(unsigned char *r, int *rlen, const secp256k1_num_t *a);
void secp256k1_num_split(secp256k1_num_t *rl, secp256k1_num_t *rh, const secp256k1_num_t *a);

Number &Number::operator=(const Number &x) {
    mpz_set(bn, x.bn);
    return *this;
}

void Number::SetNumber(const Number &x) {
    mpz_set(bn, x.bn);
}

Number::Number(const unsigned char *bin, int len) {
    mpz_init(bn);
    SetBytes(bin,len);
}

void Number::SetBytes(const unsigned char *bin, unsigned int len) {
    mpz_import(bn, len, 1, 1, 1, 0, bin);
}

bool Number::CheckBit(int pos) const {
    return mpz_tstbit(bn, pos);
}

void Number::GetBytes(unsigned char *bin, unsigned int len) {
    unsigned int size = (mpz_sizeinbase(bn,2)+7)/8;
    assert(size <= len);
    memset(bin,0,len);
    size_t count = 0;
    mpz_export(bin + len - size, &count, 1, 1, 1, 0, bn);
    assert(count == 0 || size == count);
}

void Number::SetInt(int x) {
    mpz_set_si(bn, x);
}

void Number::SetModInverse(const Number &x, const Number &m) {
    mpz_invert(bn, x.bn, m.bn);
}

void Number::SetModMul(const Number &a, const Number &b, const Number &m) {
    mpz_mul(bn, a.bn, b.bn);
    mpz_mod(bn, bn, m.bn);
}

void Number::SetAdd(const Number &a1, const Number &a2) {
    mpz_add(bn, a1.bn, a2.bn);
}

void Number::SetSub(const Number &a1, const Number &a2) {
    mpz_sub(bn, a1.bn, a2.bn);
}

void Number::SetMult(const Number &a1, const Number &a2) {
    mpz_mul(bn, a1.bn, a2.bn);
}

void Number::SetDiv(const Number &a1, const Number &a2) {
    mpz_tdiv_q(bn, a1.bn, a2.bn);
}

void Number::SetMod(const Number &a, const Number &m) {
    mpz_mod(bn, a.bn, m.bn);
}

int Number::Compare(const Number &a) const {
    return mpz_cmp(bn, a.bn);
}

int Number::GetBits() const {
    return mpz_sizeinbase(bn,2);
}

int Number::ShiftLowBits(int bits) {
    int ret = mpz_get_ui(bn) & ((1 << bits) - 1);
    mpz_fdiv_q_2exp(bn, bn, bits);
    return ret;
}

bool Number::IsZero() const {
    return mpz_size(bn) == 0;
}

bool Number::IsOdd() const {
    return mpz_get_ui(bn) & 1;
}

bool Number::IsNeg() const {
    return mpz_sgn(bn) < 0;
}

void Number::Negate() {
    mpz_neg(bn, bn);
}

void Number::Shift1() {
    mpz_fdiv_q_2exp(bn, bn, 1);
}

void Number::Inc() {
    mpz_add_ui(bn, bn, 1);
}

void Number::SetHex(const std::string &str) {
    mpz_set_str(bn, str.c_str(), 16);
}

void Number::SetPseudoRand(const Number &max) {
    number_state.gen(bn, max.bn);
}

void Number::SplitInto(int bits, Number &low, Number &high) const {
    mpz_t tmp;
    mpz_init_set_ui(tmp,1);
    mpz_mul_2exp(tmp,tmp,bits);
    mpz_sub_ui(tmp,tmp,1);
    mpz_and(low.bn, bn, tmp);
    mpz_clear(tmp);
    mpz_fdiv_q_2exp(high.bn, bn, bits);
}

std::string Number::ToString() const {
    char *str = (char*)malloc(mpz_sizeinbase(bn,16) + 2);
    mpz_get_str(str, 16, bn);
    std::string ret(str);
    free(str);
    return ret;
}

}

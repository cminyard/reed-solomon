/**
 * @file galois_field.h
 * @brief Finite Field (GF(2^m)) routines for Reed–Solomon encoding/decoding.
 *
 * This header declares:
 *   - galois_field structure
 *   - Basic GF arithmetic functions
 *   - Initialization routine for generating all tables
 *
 * All implementation details are in galois_field.c.
 */

#ifndef GALOIS_FIELD_H
#define GALOIS_FIELD_H

#include <stdint.h>
#include <assert.h>

/* Define for debugging with valgrind. */
#define GF_DYN_ALLOC 0

#if GF_DYN_ALLOC
#include <stdlib.h>
#include <stdio.h>
#endif

/* -------------------------------------------------------------------------
 * Configuration
 * ------------------------------------------------------------------------- */
/* Maximum GF size = 2^8 = 256 */
#define GF_EXP_MAX 8
#define GF_MAX (1 << GF_EXP_MAX)

typedef uint8_t gf_sym;

/* -------------------------------------------------------------------------
 * GF tables and polynomial data
 * ------------------------------------------------------------------------- */
struct galois_field {
    unsigned int Np; /* parent field length */
#if GF_DYN_ALLOC
    gf_sym *exp;
    gf_sym *log;
#else
    /* This is one bigger than necessary, see galois_field.c. */
    gf_sym exp[GF_MAX];    /* Exponential table */

    gf_sym log[GF_MAX];    /* Logarithm table */
#endif
};

/* -------------------------------------------------------------------------
 * GF(2^m) arithmetic primitives
 * ------------------------------------------------------------------------- */

/**
 * @brief alpha ^ a, a must be a valid gf number.
 */
static inline gf_sym
gf_exp(struct galois_field *gf, gf_sym a)
{
    return gf->exp[a];
}

/**
 * @brief alpha ^ a, a may be an arbitrary number.
 */
static inline gf_sym
gf_exp_o(struct galois_field *gf, unsigned int a)
{
    return gf->exp[a % gf->Np];
}

/**
 * @brief GF addition (same as XOR).
 */
static inline gf_sym
gf_add(gf_sym a, gf_sym b)
{
    return a ^ b;
}

/**
 * @brief GF multiplication using exp/log tables.
 */
static inline gf_sym
gf_mul(struct galois_field *gf, gf_sym a, gf_sym b)
{
    int idx;

    if (a == 0 || b == 0)
	return 0;

    idx = gf->log[a] + gf->log[b];
    if (idx >= gf->Np)
	idx -= gf->Np;

    return gf->exp[idx];
}

/**
 * @brief GF division using exp/log tables.
 */
static inline gf_sym
gf_div(struct galois_field *gf, gf_sym a, gf_sym b)
{
    int idx;

    /* FIXME - handle this. */
    assert(b > 0); /* Division by zero. */

    if (a == 0)
	return 0;

    idx = (int) gf->log[a] - (int) gf->log[b];
    if (idx < 0)
	idx += gf->Np;

    return gf->exp[idx];
}

/**
 * @brief Raise base to an integer power (base^power) in GF.
 */
static inline
gf_sym gf_pow(struct galois_field *gf, gf_sym base, int power)
{
    int logv;
    int x;

    if (base == 0)
	return 0;

    logv = gf->log[base];
    x = (logv * power) % gf->Np;
    if (x < 0)
	x += gf->Np;
    return gf->exp[x];
}

/**
 * @brief like gf_pow(), but base is already in log format.
 */
static inline
gf_sym gf_pow_nl(struct galois_field *gf, gf_sym base, int power)
{
    int x;

    if (base == gf->Np)
	return 0;

    x = (base * power) % gf->Np;
    if (x < 0)
	x += gf->Np;
    return gf->exp[x];
}

/**
 * @brief Multiplicative inverse in GF.
 */
static inline gf_sym
gf_inv(struct galois_field *gf, gf_sym a)
{
    if (a == 0)
	return 0;

    /*
     * gf->exp has a value for gf->Np to handle when gf->log[a] == 0.
     * See galois_field.c for details.
     */
    return gf->exp[gf->Np - gf->log[a]];
}

/* -------------------------------------------------------------------------
 * Initialization
 * ------------------------------------------------------------------------- */

/**
 * @brief Initialize GF(2^m).
 *
 * @param m  GF size parameter (1–8), GF size = 2^m
 *
 * @return 0 on success, non-zero on failure.
 */
int galois_field_init(struct galois_field *gf,
		      unsigned int m, unsigned int poly);

#endif /* GALOIS_FIELD_H */

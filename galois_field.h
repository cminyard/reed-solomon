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

/* -------------------------------------------------------------------------
 * Configuration
 * ------------------------------------------------------------------------- */
/* Maximum GF size = 2^8 = 256 */
#define GALOIS_FIELD_EXP_MAX 8
#define GALOIS_FIELD_MAX (1 << GALOIS_FIELD_EXP_MAX)

typedef uint16_t galois_field_val;

/* -------------------------------------------------------------------------
 * GF tables and polynomial data
 * ------------------------------------------------------------------------- */
struct galois_field {
    unsigned int Np; /* parent field length */
    galois_field_val exp[2 * GALOIS_FIELD_MAX];/* Exponential table */
    galois_field_val log[GALOIS_FIELD_MAX];    /* Logarithm table */
};

/* -------------------------------------------------------------------------
 * GF(2^m) arithmetic primitives
 * ------------------------------------------------------------------------- */

/**
 * @brief GF addition (same as XOR).
 */
static inline galois_field_val
galois_field_add(galois_field_val a, galois_field_val b)
{
    return a ^ b;
}

/**
 * @brief GF multiplication using exp/log tables.
 */
static inline galois_field_val
galois_field_mul(struct galois_field *gf,
		 galois_field_val a, galois_field_val b)
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
static inline galois_field_val
galois_field_div(struct galois_field *gf,
		 galois_field_val a, galois_field_val b)
{
    int idx;

    if (b == 0) /* Division by zero. */
	return 0;

    if (a == 0)
	return 0;

    idx = gf->log[a] - gf->log[b];
    if (idx < 0)
	idx += gf->Np;

    return gf->exp[idx];
}

/**
 * @brief Raise base to an integer power (base^power) in GF.
 */
static inline
galois_field_val galois_field_pow(struct galois_field *gf,
				  galois_field_val base, int power)
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
 * @brief Multiplicative inverse in GF.
 */
static inline galois_field_val
galois_field_inv(struct galois_field *gf, galois_field_val a)
{
    if (a == 0)
	return 0;

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
int galois_field_init(struct galois_field *gf, int m);

#endif /* GALOIS_FIELD_H */

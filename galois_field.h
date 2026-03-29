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
galois_field_val galois_field_add(galois_field_val a, galois_field_val b);

/**
 * @brief GF multiplication using exp/log tables.
 */
galois_field_val galois_field_mul(struct galois_field *gf,
				  galois_field_val a, galois_field_val b);

/**
 * @brief GF division using exp/log tables.
 */
galois_field_val galois_field_div(struct galois_field *gf,
				  galois_field_val a, galois_field_val b);

/**
 * @brief Raise base to an integer power (base^power) in GF.
 */
galois_field_val galois_field_pow(struct galois_field *gf,
				  galois_field_val base, int power);

/**
 * @brief Multiplicative inverse in GF.
 */
galois_field_val galois_field_inv(struct galois_field *gf,
				  galois_field_val a);

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

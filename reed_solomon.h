/**
 * @file reed_solomon.h
 * @brief Reed–Solomon encoding/decoding.
 *
 * This header declares:
 *   - reed_solomon structure
 *   - reed_solomon encoder structure
 *   - reed_solomon decoder structure
 *   - Encoder, decoder, and initialization functions
 *
 * All implementation details are in reed_solomon.c.
 */

#ifndef REED_SOLOMON_H
#define REED_SOLOMON_H

#include <stdint.h>
#include "galois_field.h"

/* Maximum value T may be. */
#define REED_SOLOMON_MAX_T 32

struct reed_solomon {
    unsigned int m;  /* GF size parameter m → GF(2^m) */
    unsigned int T;  /* Number of parity symbols (generator degree) */

    struct galois_field gf;

    /* FIXME - allocate these based on m. */

    /* Generator polynomial g(x) */
    galois_field_val generator[GALOIS_FIELD_MAX];

    /* Bit representation table */
    int symbol_bits[GALOIS_FIELD_MAX][GALOIS_FIELD_EXP_MAX];
};

/**
 * @brief Initialize GF(2^m) and construct RS generator polynomial.
 *
 * @param rs The structure to fill in with data
 * @param m  GF size parameter (1–8), GF size = 2^m
 * @param T  Number of parity symbols
 *
 * The length of the actual data may vary up to (2^m - T) bytes long,
 * and that is passed into the encode and decode functions.  The
 * arrays and lengths passed into those function include T, see them
 * for details.
 *
 * @return 0 on success, non-zero on failure.
 */
int reed_solomon_init(struct reed_solomon *rs, unsigned int m, unsigned int T);

struct reed_solomon_encoder {
    struct reed_solomon *rs;

    unsigned int N;  /* Codeword length (shortened) */
    unsigned int S;  /* Shortening amount = gf.Np - N */
    unsigned int K;  /* Number of information symbols */

    /* Needs to be 2^M - T elements; */
    galois_field_val u[GALOIS_FIELD_MAX];
    /* Needs to be T elements. */
    galois_field_val parity[REED_SOLOMON_MAX_T];
};

/**
 * @brief Systematic Reed–Solomon encoding.
 *
 * @param rse A encoder structure with rs set.
 * @param buf Buffer of len + T bytes, the last T bytes are replaced with parity.
 * @param len Length of buf, including the T bytes.
 *
 * The length of the buffer may be up to 2^M symbols long, meaning
 * that the actual data may be (2^M - T) symbols long.
 *
 * @return 0 on success, non-zero on error.
 */
int reed_solomon_encode(struct reed_solomon_encoder *rse,
			 uint8_t *buf, unsigned int len);

/* We can process up to T/2 errors.  More than that we ignore. */
#define REED_SOLOMON_MAX_ERR (REED_SOLOMON_MAX_T / 2)

struct reed_solomon_decoder {
    struct reed_solomon *rs;

    unsigned int N;  /* Codeword length (shortened) */
    unsigned int S;  /* Shortening amount = gf.Np - N */
    unsigned int K;  /* Number of information symbols */

    galois_field_val C[GALOIS_FIELD_MAX]; /* current polynomial */
    galois_field_val B[GALOIS_FIELD_MAX]; /* previous polynomial */
    galois_field_val Temp[GALOIS_FIELD_MAX];

    galois_field_val A[REED_SOLOMON_MAX_ERR][REED_SOLOMON_MAX_ERR];
    galois_field_val e[REED_SOLOMON_MAX_ERR];

    /* galois_field_val recv_sym_p[Np]; */
    galois_field_val recv_sym_p[GALOIS_FIELD_MAX];

    galois_field_val synd[REED_SOLOMON_MAX_T];
    galois_field_val sigma[REED_SOLOMON_MAX_ERR + 1];
    unsigned int error_pos[REED_SOLOMON_MAX_ERR];
};

/**
 * @brief Decode a shortened systematic Reed–Solomon codeword.
 *
 * @param buf Buffer, including the parity.
 * @param lan Length of the buffer.
 *
 * The first K bytes of buf are replaced with the corrected data.
 *
 * The length of the buffer may be up to 2^M symbols long, meaning
 * that the actual data may be (2^M - T), or K symbols long.
 *
 * @return 0 on success, non-zero on error.
 */
int reed_solomon_decode(struct reed_solomon_decoder *rsd,
			uint8_t *buf, unsigned int len,
			unsigned int *err_count);

#endif /* REED_SOLOMON_H */

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

struct reed_solomon {
    unsigned int m;  /* GF size parameter m → GF(2^m) */
    unsigned int N;  /* Codeword length (shortened) */
    unsigned int S;  /* Shortening amount = gf.Np - N */
    unsigned int K;  /* Number of information symbols */
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
 * @param N  Codeword length (shortened)
 * @param K  Information symbol length
 *
 * @return 0 on success, negative on failure.
 */
int reed_solomon_init(struct reed_solomon *rs, int m, int N, int K);

struct reed_solomon_encoder {
    /* FIXME - allocate these based on K and T. */

    struct reed_solomon *rs;

    /*
     * FIXME - do these need to be galois_field_val?  If they are uint8_t,
     * they can be used directly from the user buffer provided.
     */
    /* galois_field_val u[K]; */
    galois_field_val u[GALOIS_FIELD_MAX];
    /* galois_field_val parity[T]; */
    galois_field_val parity[GALOIS_FIELD_MAX];
};

/**
 * @brief Systematic Reed–Solomon encoding.
 *
 * @param rse A encoder structure with rs set.
 * @param buf Buffer of N bytes, the last (N - K) bits are replaced with parity.
 */
void reed_solomon_encode(struct reed_solomon_encoder *rse, uint8_t *buf);

struct reed_solomon_decoder {
    struct reed_solomon *rs;
    /*
     * FIXME - allocate these based on parms. Make sure to fix the
     * memset in reed_solomon_decode().
     */

    galois_field_val C[GALOIS_FIELD_MAX]; /* current polynomial */
    galois_field_val B[GALOIS_FIELD_MAX]; /* previous polynomial */
    galois_field_val Temp[GALOIS_FIELD_MAX];

    /* galois_field_val A[error_count][error_count]; */
    galois_field_val A[GALOIS_FIELD_MAX][GALOIS_FIELD_MAX];
    /* galois_field_val B[error_count]; */
    /* galois_field_val B[GALOIS_FIELD_MAX]; - Reused the "B" above */
    /* galois_field_val e[error_count]; */
    galois_field_val e[GALOIS_FIELD_MAX];

    /* galois_field_val recv_sym_p[Np]; */
    galois_field_val recv_sym_p[GALOIS_FIELD_MAX];
    /* galois_field_val synd[T]; */
    galois_field_val synd[GALOIS_FIELD_MAX];
    /* galois_field_val sigma[rs->T / 2 + 1]; */
    galois_field_val sigma[GALOIS_FIELD_MAX];
    /* unsigned int error_pos[rs->T / 2]; */
    unsigned int error_pos[GALOIS_FIELD_MAX];
};

/**
 * @brief Decode a shortened systematic Reed–Solomon codeword.
 *
 * @param buf Buffer, N bytes long.
 *
 * The first K bytes of buf are replaced with the corrected data.
 */
unsigned int reed_solomon_decode(struct reed_solomon_decoder *rsd, uint8_t *buf);

#endif /* REED_SOLOMON_H */

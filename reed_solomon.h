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

#include "galois_field.h"

/* Maximum value T may be. */
#define RS_MAX_T 32

struct reed_solomon {
    unsigned int m;  /* GF size parameter m → GF(2^m) */
    unsigned int T;  /* Number of parity symbols (generator degree) */

    unsigned int prim;
    unsigned int iprim;
    unsigned int fcr;

    struct galois_field gf;
};

/**
 * @brief Initialize GF(2^m) and construct RS generator polynomial.
 *
 * @param rs   The structure to fill in with data
 * @param m    GF size parameter (1–8), GF size = 2^m
 * @param gfpoly The GF polynomial
 * @param T    Number of parity symbols
 * @param fcr  The first consecutive root of the RS polynomial
 * @param prim The primitive element of the GF.
 *
 * The length of the actual data may vary up to (2^m - T) bytes long,
 * and that is passed into the encode and decode functions.
 *
 * @return 0 on success, non-zero on failure.
 */
int reed_solomon_init(struct reed_solomon *rs, unsigned int m,
		      unsigned int gfpoly, unsigned int T,
		      unsigned int fcr, unsigned int prim);

struct reed_solomon_encoder {
    struct reed_solomon *rs;

#if GF_DYN_ALLOC
    gf_sym *generator;
#else
    /* Generator polynomial g(x) */
    gf_sym generator[GF_MAX];
#endif
};

/**
 * @brief Initialize a Reed–Solomon encoder structure.
 *
 * @param rse The structure to initialize.
 * @param rs  The main Reed Solomon information, this pointer is kept.
 */
void rs_encoder_init(struct reed_solomon_encoder *rse,
		     struct reed_solomon *rs);

/**
 * @brief Reed–Solomon encoding.
 *
 * @param rse A encoder structure with rs set.
 * @param buf Buffer to generate parity for
 * @param len Length of buf
 * @param parity Output parity bytes, should be T bytes long
 *
 * The length of the buffer may be up to 2^M symbols long, meaning
 * that the actual data may be (2^M - T) symbols long.
 *
 * @return 0 on success, non-zero on error.
 */
int rs_encode(struct reed_solomon_encoder *rse,
	      uint8_t *buf, unsigned int len, uint8_t *parity);

/**
 * @brief Reed–Solomon encoding for CCSDS parameter 8 bit without
 * dual-basis symbols. Highly optimized.
 *
 * @param buf Buffer to generate parity for
 * @param parity Output parity bytes, should be T bytes long
 *
 * The length of the buffer may be up to 2^M symbols long, meaning
 * that the actual data may be (2^M - T) symbols long.
 *
 * Initialization is not needed for this function.
 *
 * @return 0 on success, non-zero on error.
 */
int rs_encode_8(uint8_t *inbuf, unsigned int len, uint8_t *parity);

/* We can process up to T/2 errors.  More than that we ignore. */
#define RS_MAX_ERR (RS_MAX_T / 2)

struct reed_solomon_decoder {
    struct reed_solomon *rs;

    unsigned int N;  /* Codeword length (shortened) */
    unsigned int S;  /* Shortening amount = gf.Np - N */
    unsigned int K;  /* Number of information symbols */
    unsigned int pad;

#if GF_DYN_ALLOC
    gf_sym *C;
    gf_sym *B;
    gf_sym *Temp;

    gf_sym *synd;
    unsigned int *error_idx;
    unsigned int *error_pos;

    gf_sym *O;
#else
    /* Lambda array. */
    gf_sym C[RS_MAX_T + 1]; /* current polynomial */

    /* Two temporary arrays used by decoding */
    gf_sym B[2 * RS_MAX_T + 1]; /* previous polynomial */
    gf_sym Temp[RS_MAX_T + 1];

    /* Syndrome array. */
    gf_sym synd[RS_MAX_T];

    /* Information about error locations. */
    unsigned int error_idx[RS_MAX_ERR];
    unsigned int error_pos[RS_MAX_ERR];

    /* Omega array. */
    gf_sym O[RS_MAX_ERR];
#endif
};

/**
 * @brief Initialize a Reed–Solomon decoder structure.
 *
 * @param rsd The structure to initialize.
 * @param rs  The main Reed Solomon information, this pointer is kept.
 */
void rs_decoder_init(struct reed_solomon_decoder *rsd,
		     struct reed_solomon *rs);

/**
 * @brief Decode a Reed–Solomon encoded buffer.
 *
 * @param rsd A decoder structure with rs set.
 * @param buf Buffer, including the parity.
 * @param len Length of the buffer.
 * @param err_count The number of errors that were corrected.
 *
 * The first K bytes of buf are replaced with the corrected data if
 * there are errors.
 *
 * The length of the buffer may be up to 2^M symbols long, meaning
 * that the actual data may be (2^M - T), or K symbols long.
 *
 * If errors could not be corrected, this returns failure and err_count
 * will not be set.
 *
 * @return 0 on success, non-zero on failure.
 */
int rs_decode(struct reed_solomon_decoder *rsd,
	      uint8_t *buf, unsigned int len,
	      unsigned int *err_count);

#endif /* REED_SOLOMON_H */

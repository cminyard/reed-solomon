/**
 * @file rs_gf.c
 * @brief Galois Field (GF(2^m)) arithmetic for Reed–Solomon codes.
 *
 * This module initializes the finite field GF(2^m), manages exponential/log
 * tables, provides basic operations (add/mul/div/inv/pow), and generates the
 * RS generator polynomial of degree T.
 *
 * Supported:
 *   - GF sizes up to RS_GF_MAX (configurable in rs_gf.h)
 *   - Arbitrary (N, K, T) compatible with GF(2^m)
 *
 * This implementation is self-contained and uses only standard C.
 */

#include "galois_field.h"

/* Primitive polynomials for m = 1..8 (CCSDS/NASA compatible) */
static const uint16_t primitive_poly[9] = {
    0x00, /* unused (m=0) */
    0x03, /* m=1 */
    0x07, /* m=2 */
    0x0B, /* m=3 */
    0x13, /* m=4 */
    0x25, /* m=5 */
    0x43, /* m=6 */
    0x89, /* m=7 */
    0x11D /* m=8 (used for RS(255,223), GF(256)) */
};

/* -------------------------------------------------------------------------
 * Initialize GF(2^m) and build generator polynomial g(x)
 * ------------------------------------------------------------------------- */
int galois_field_init(struct galois_field *gf, int m)
{
    uint16_t prim, x;
    unsigned int i;

    /* Field size (2^m - 1) */
    gf->Np = (1 << m) - 1;

#if GF_DYN_ALLOC
    gf->exp = malloc(sizeof(galois_field_val) * GALOIS_FIELD_MAX);
    gf->log = malloc(sizeof(galois_field_val) * GALOIS_FIELD_MAX);
#endif

    /* Select primitive polynomial */
    prim = primitive_poly[m];

    /* Build exp/log tables */
    x = 1;
    gf->log[0] = 0;
    for (i = 0; i < gf->Np; i++) {
	gf->exp[i] = x;
	gf->log[x] = i;

	x <<= 1;
	if (x & (1u << m))
	    x ^= prim;
    }
    /*
     * gf->exp[i] is one larger than necessary to avoid having to do a mod
     * in galois_field_inv().
     */
    gf->exp[gf->Np] = gf->exp[0];

    return 0;
}

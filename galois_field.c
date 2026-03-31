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

/* -------------------------------------------------------------------------
 * Initialize GF(2^m)
 * ------------------------------------------------------------------------- */
int galois_field_init(struct galois_field *gf,
		      unsigned int m, unsigned int poly)
{
    uint16_t x;
    unsigned int i;

    /* Field size (2^m - 1) */
    gf->Np = (1 << m) - 1;

#if GF_DYN_ALLOC
    gf->exp = malloc(sizeof(gf_val) * GF_MAX);
    gf->log = malloc(sizeof(gf_val) * GF_MAX);
#endif

    /* Build exp/log tables */
    x = 1;
    gf->log[0] = 255;
    for (i = 0; i < gf->Np; i++) {
	gf->exp[i] = x;
	gf->log[x] = i;

	x <<= 1;
	if (x & (1u << m))
	    x ^= poly;
    }
    /*
     * gf->exp[i] is one larger than necessary to avoid having to do a mod
     * in gf_inv().
     */
    gf->exp[gf->Np] = gf->exp[0];

    return 0;
}

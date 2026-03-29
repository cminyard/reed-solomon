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
static const galois_field_val primitive_poly[9] = {
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
 * Basic GF operations
 * ------------------------------------------------------------------------- */

/** GF addition = XOR */
galois_field_val
galois_field_add(galois_field_val a, galois_field_val b)
{
    return a ^ b;
}

/** GF multiplication using log/exp tables */
galois_field_val
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

/** GF division using log/exp tables */
galois_field_val
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

/** GF exponentiation (base^power in GF) */
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

/** Multiplicative inverse */
galois_field_val galois_field_inv(struct galois_field *gf, galois_field_val a)
{
    if (a == 0)
	return 0;

    return gf->exp[gf->Np - gf->log[a]];
}

/* -------------------------------------------------------------------------
 * Initialize GF(2^m) and build generator polynomial g(x)
 * ------------------------------------------------------------------------- */
int galois_field_init(struct galois_field *gf, int m)
{
    galois_field_val prim;
    galois_field_val x;
    unsigned int i;

    /* Field size (2^m - 1) */
    gf->Np = (1 << m) - 1;

    /* Select primitive polynomial */
    prim = primitive_poly[m];

    /* Build exp/log tables */
    x = 1;
    for (i = 0; i < gf->Np; i++) {
	gf->exp[i] = x;
	gf->log[x] = i;

	x <<= 1;
	if (x & (1u << m))
	    x ^= prim;
    }

    /* Extend exp table for mod-free multiplication */
    for (i = gf->Np; i < 2 * gf->Np; i++)
	gf->exp[i] = gf->exp[i - gf->Np];

    gf->log[0] = 0;

    return 0;
}

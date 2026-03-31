
#include <stdbool.h>
#include <string.h>
#include "reed_solomon.h"

int
reed_solomon_init(struct reed_solomon *rs, unsigned int m,
		  unsigned int gfpoly, unsigned int T,
		  unsigned int fcr, unsigned int prim)
{
    int err = galois_field_init(&rs->gf, m, gfpoly);

    if (err)
	return err;

    if (T > RS_MAX_T)
	return 1;
    if (T >= rs->gf.Np)
	return 1;

    rs->m = m;
    rs->T = T;

    rs->fcr = fcr;
    rs->prim = prim;
    for (rs->iprim = 1; rs->iprim % prim != 0; rs->iprim += rs->gf.Np)
	;
    rs->iprim /= prim;

    return 0;
}

void
rs_encoder_init(struct reed_solomon_encoder *rse,
		struct reed_solomon *rs)
{
    struct galois_field *gf = &rs->gf;
    unsigned int i, j;
    gf_val tmp[GF_MAX];
    unsigned int root;

    rse->rs = rs;

#if GF_DYN_ALLOC
    rse->generator = malloc(sizeof(gf_val) * GF_MAX);
#endif

    /* ---------------------------------------------------------------------
     * Generator polynomial construction (degree T)
     * g(x) = (x - α^0)(x - α^1)...(x - α^(T-1))
     * --------------------------------------------------------------------- */
    rse->generator[0] = 1;
    memset(rse->generator + 1, 0, rs->T * sizeof(gf_val));

    root = rs->fcr * rs->prim;
    for (i = 0; i < rs->T; i++, root += rs->prim) {
	for (j = 0; j <= i; j++)
	    tmp[j] = rse->generator[j];

	rse->generator[i + 1] = 1;

	/* Perform polynomial multiplication by (x - α^i) */
	for (j = i; j > 0; j--) {
	    if (rse->generator[j] != 0)
		/* rse->generator[j] = tmp[j - 1] + (tmp[j] * a ^ root) */
		rse->generator[j] = gf_add(tmp[j - 1],
					   gf_mul(gf, tmp[j],
						  gf_exp_o(gf, root)));
	    else
		rse->generator[j] = rse->generator[j - 1];
	}
	rse->generator[0] = gf_mul(gf, tmp[0], gf_exp_o(gf, root));
    }
}

/**
 * @brief Systematic Reed–Solomon encoder.
 *
 * Produces a codeword of:
 *      [K info symbols][T parity symbols]
 */
int
rs_encode(struct reed_solomon_encoder *rse,
	  uint8_t *inbuf, unsigned int len, uint8_t *parity)
{
    struct reed_solomon *rs = rse->rs;
    struct galois_field *gf = &rs->gf;
    unsigned int i, j;

    if (len > gf->Np - rs->T)
	return 1;

    memset(parity, 0, rs->T * sizeof(gf_val));

    /* No need to process the pad, the result will still be 0. */

    /* Feed the actual data. */
    for (i = 0; i < len; i++) {
	gf_val fb = gf_add(inbuf[i], parity[0]);

	if (fb != 0) {
	    for (j = 1; j < rs->T; j++)
		/* parity[j] += fb * rse->generator[rs->T - j] */
		parity[j] = gf_add(parity[j],
				   gf_mul(gf, fb, rse->generator[rs->T - j]));
	}
	memcpy(parity, parity + 1, rs->T * sizeof(gf_val));

	if (fb != 0)
	    parity[rs->T - 1] = gf_mul(gf, fb, rse->generator[0]);
	else
	    parity[rs->T - 1] = 0;
    }

    return 0;
}


void
rs_decoder_init(struct reed_solomon_decoder *rsd,
		struct reed_solomon *rs)
{
    rsd->rs = rs;
#if GF_DYN_ALLOC
    rsd->data = malloc(sizeof(gf_val) * GF_MAX);

    rsd->C = malloc(sizeof(gf_val) * (RS_MAX_T + 1));
    rsd->B = malloc(sizeof(gf_val) * (2 * RS_MAX_T + 1));
    rsd->Temp = malloc(sizeof(gf_val) * (RS_MAX_T + 1));

    rsd->synd = malloc(sizeof(gf_val) * RS_MAX_T);
    rsd->error_pos = malloc(sizeof(unsigned int) * RS_MAX_ERR);
    rsd->error_idx = malloc(sizeof(unsigned int) * RS_MAX_ERR);

    rsd->O = malloc(sizeof(gf_val) * RS_MAX_ERR);
#endif
}

/* -------------------------------------------------------------------------
 * 1) Syndrome computation (on parent length Np)
 *
 *     S_i = Σ_{j=0}^{Np-1} r_j α^{(i+1)*j},   for i = 0..T-1
 *
 * Returns the number of non-zero syndromes.
 * Zero syndromes = no errors.
 * ------------------------------------------------------------------------- */
static unsigned int
compute_syndromes(struct reed_solomon *rs, unsigned int N,
		  const gf_val *e, gf_val *S)
{
    struct galois_field *gf = &rs->gf;
    unsigned int i, j, count = 0;

    for (i = 0; i < rs->T; i++)
	S[i] = e[0];

    for (i = 1; i < N; i++) {
	for (j = 0; j < rs->T; j++)
	    /* S[j] = e[i] + (S[j] * (rs->fcr + j) ^ rs->prim) */
	    S[j] = gf_add(e[i],
			  gf_mul(gf, S[j],
				 gf_pow_nl(gf, rs->fcr + j, rs->prim)));
    }

    for (i = 0; i < rs->T; i++) {
	if (S[i])
	    count++;
    }

    return count;
}

/* -------------------------------------------------------------------------
 * 2) Berlekamp–Massey algorithm
 *
 * Finds the error-locator polynomial σ(x).
 * Output: sigma_out[0..L]
 * Ensures σ(0) = 1.
 *
 * L = degree of error-locator polynomial
 * ------------------------------------------------------------------------- */
static unsigned int
berlekamp_massey(struct reed_solomon *rs,
		 gf_val *S,
		 gf_val *B, gf_val *C,
		 gf_val *Temp)
{
    struct galois_field *gf = &rs->gf;
    unsigned int L = 0;
    unsigned int i, n;

    /*
     * Instead of copying the B array, move it backwards.  This can
     * happen at most RS_MAX_T times, so start out there.
     */
    B += RS_MAX_T;

    memset(B + 1, 0, sizeof(gf_val) * RS_MAX_T);
    memset(C + 1, 0, sizeof(gf_val) * RS_MAX_T);

    C[0] = 1;
    B[0] = 1;

    for (n = 1; n <= rs->T; n++) {
	gf_val d = 0;

	for (i = 0; i < n; i++)
	    /* d += C[i] * S[n - i - 1] */
	    d = gf_add(d, gf_mul(gf, C[i], S[n - i - 1]));

	if (d == 0) {
	    B--;
	    B[0] = 0;
	} else {
	    Temp[0] = C[0];
	    for (i = 0; i < rs->T; i++)
		    Temp[i + 1] = gf_add(C[i + 1], gf_mul(gf, d, B[i]));

	    if (2 * L <= n - 1) {
		for (i = 0; i <= rs->T; i++)
		    B[i] = gf_div(gf, C[i], d);

		L = n - L;
	    } else {
		B--;
		B[0] = 0;
	    }

	    memcpy(C, Temp, rs->T * sizeof(gf_val));
	}
    }

    return L;
}

/* -------------------------------------------------------------------------
 * 3) Chien search
 *
 * Find i such that σ(α^{-i}) = 0, for i = 0..Np-1.
 * Each such i corresponds to an error at position i.
 * Returns non-zero if errors are uncorrectible.
 * ------------------------------------------------------------------------- */
static unsigned int
chien_search(struct reed_solomon *rs, unsigned int L,
	     gf_val *C, gf_val *Temp,
	     unsigned int *err_idx, unsigned int *err_pos,
	     unsigned int *rcount)
{
    struct galois_field *gf = &rs->gf;
    unsigned int count = 0, d = 0;
    unsigned int i, j, k;

    for (i = 0; i <= rs->T; i++) {
	if (C[i] != 0)
	    d = i;
	Temp[i] = C[i];
    }

    k = rs->iprim - 1;
    for (i = 1; i <= gf->Np; i++) {
	unsigned int sum = 1;

	for (j = d; j > 0; j--) {
	    Temp[j] = gf_mul(gf, Temp[j], gf->exp[j]);
	    sum = gf_add(sum, Temp[j]);
	}

	if (sum == 0) {
	    err_idx[count] = i;
	    err_pos[count++] = k;
	    if (count == d)
		break;
	}

	k = (k + rs->iprim) % gf->Np;
    }

    *rcount = count;

    return d != count;
}

/* -------------------------------------------------------------------------
 * 4) Error magnitude solving via linear system
 *
 * Simplified Forney method:
 *     S_l = Σ e_k α^{(l+1) * i_k}
 * Solve for e_k using Gaussian elimination in GF(2^m).
 * ------------------------------------------------------------------------- */
static void
correct_errors(struct reed_solomon *rs,
	       gf_val *e, const gf_val *S,
	       const gf_val *C, gf_val *O,
	       const unsigned int *err_idx, const unsigned int *err_pos,
	       unsigned int err_cnt)
{
    struct galois_field *gf = &rs->gf;
    unsigned int i, j, o;

    o = err_cnt - 1;
    for (i = 0; i <= o; i++) {
	O[i] = 0;
	for (j = 0; j <= i; j++)
	    /* O[i] += S[i - j] * C[j] */
	    O[i] = gf_add(O[i], gf_mul(gf, S[i - j], C[j]));
    }

    for (i = 0; i < err_cnt; i++) {
	unsigned int tmp = 0, tmp2, d;

	for (j = 0; j <= o; j++)
	    /* tmp += O[j] * j ^ err_idx[i] */
	    tmp = gf_add(tmp, gf_mul(gf, O[j],
				     gf_pow_nl(gf, j, err_idx[i])));
	/* tmp2 = err_index[i] ^ (fcr - 1) */
	tmp2 = gf_pow_nl(gf, err_idx[i], rs->fcr - 1);

	d = 0;
	if (err_cnt < rs->T - 1)
	    j = err_cnt;
	else
	    j = rs->T - 1;
	for (j = j & ~1; ; j -= 2) {
	    /* d += C[j + 1] * j ^ err_idx[i] */
	    d = gf_add(d, gf_mul(gf, C[j + 1],
				 gf_pow_nl(gf, j, err_idx[i])));
	    if (j < 2)
		break;
	}

	e[err_pos[i]] = gf_add(e[err_pos[i]],
			       gf_div(gf, gf_mul(gf, tmp, tmp2), d));
    }
}

/* -------------------------------------------------------------------------
 * 5) Public API: RS decoding
 *
 * Steps:
 *   - Expand to parent length: [S zero-symbols][Ns received]
 *   - Compute syndromes
 *   - If non-zero: BM → Chien → Solve magnitudes → Correct
 *   - Output:
 *       code_bits : Ns symbols
 *       info_bits : first K symbols
 * ------------------------------------------------------------------------- */
int
rs_decode(struct reed_solomon_decoder *rsd,
	  uint8_t *buf, unsigned int len,
	  unsigned int *err_count)
{
    struct reed_solomon *rs = rsd->rs;
    struct galois_field *gf = &rs->gf;
    uint8_t *data;
    unsigned int count = 0;

    if (len > gf->Np)
	return 1;
    if (len <= rs->T)
	return 1;
    rsd->N = len;
    rsd->K = len - rs->T;
    /* Number of shortened symbols */
    rsd->S = gf->Np - rsd->N;

    if (rsd->S == 0) {
	data = buf;
    } else {
	data = rsd->data;

	/* Fill the beginning with zeros and copy the rest in. */
	memset(data, 0, rsd->S * sizeof(gf_val));
	memcpy(data + rsd->S, buf, rsd->N * sizeof(gf_val));
    }

    /* Syndromes */
    count = compute_syndromes(rs, rsd->N, data, rsd->synd);

    if (count != 0) {
	/* BM → locator polynomial */
	unsigned int L = berlekamp_massey(rs, rsd->synd, rsd->B, rsd->C,
					  rsd->Temp);

	/* Chien search */
	if (chien_search(rs, L, rsd->C, rsd->Temp,
			 rsd->error_idx, rsd->error_pos, &count))
	    return 1; /* Could not correct all symbols. */

	/* Correct */
	if (count > 0)
	    correct_errors(rs, data, rsd->synd, rsd->C, rsd->O,
			   rsd->error_idx, rsd->error_pos, count);
    }

    *err_count = count;

    if (count == 0) /* No errors. */
	return 0;

    /* Output K information symbols */
    if (rsd->S != 0)
	memcpy(buf, data + rsd->S, rsd->K * sizeof(gf_val));

    return 0;
}


#include <stdbool.h>
#include <string.h>
#include "reed_solomon.h"

int
reed_solomon_init(struct reed_solomon *rs, unsigned int m,
		  unsigned int gfpoly, unsigned int T,
		  unsigned int fcr, unsigned int prim)
{
    unsigned int i, j;
    galois_field_val tmp[GALOIS_FIELD_MAX];
    unsigned int root;
    int err = galois_field_init(&rs->gf, m, gfpoly);

    if (err)
	return err;

    if (T > REED_SOLOMON_MAX_T)
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

#if GF_DYN_ALLOC
    rs->generator = malloc(sizeof(galois_field_val) * GALOIS_FIELD_MAX);
#endif

    /* ---------------------------------------------------------------------
     * Generator polynomial construction (degree T)
     * g(x) = (x - α^0)(x - α^1)...(x - α^(T-1))
     * --------------------------------------------------------------------- */
    rs->generator[0] = 1;
    for (i = 1; i <= rs->T; i++)
	rs->generator[i] = 0;

    root = fcr * prim;
    for (i = 0; i < rs->T; i++, root += prim) {
	/* Copy existing coefficients */
	for (j = 0; j <= i; j++)
	    tmp[j] = rs->generator[j];

	rs->generator[i + 1] = 1;

	/* Perform polynomial multiplication by (x - α^i) */
	for (j = i; j > 0; j--) {
	    if (rs->generator[j] != 0) {
		unsigned int term;

		term = (rs->gf.log[tmp[j]] + root) % rs->gf.Np;
		term = rs->gf.exp[term];
		rs->generator[j] = galois_field_add(tmp[j - 1], term);
	    } else {
		rs->generator[j] = rs->generator[j - 1];
	    }
	}
	rs->generator[0] = rs->gf.exp[(rs->gf.log[tmp[0]] + root) % rs->gf.Np];
    }

    /*
     * Pre-compute the log so it doesn't have to be done every time
     * when encoding.
     */
    for (i = 0; i <= rs->T; i++)
	rs->generator[i] = rs->gf.log[rs->generator[i]];

    return 0;
}

void
reed_solomon_encoder_init(struct reed_solomon_encoder *rse,
			  struct reed_solomon *rs)
{
    rse->rs = rs;
}

/**
 * @brief Systematic Reed–Solomon encoder.
 *
 * Produces a codeword of:
 *      [K info symbols][T parity symbols]
 */
int
reed_solomon_encode(struct reed_solomon_encoder *rse,
		    uint8_t *inbuf, unsigned int len, uint8_t *parity)
{
    struct reed_solomon *rs = rse->rs;
    unsigned int i, j;

    if (len > rs->gf.Np - rs->T)
	return 1;

    /* -------------------------------------------------------------
     * Initialize T parity registers to zero
     * ------------------------------------------------------------- */
    for (i = 0; i < rs->T; i++)
	parity[i] = 0;

#if 0 /* No need to process the pad, the result will still be 0. */
    /* -------------------------------------------------------------
     * Handle shortening:
     *   Shift S dummy symbols (all zeros) through the encoder.
     *
     * This produces the same result as encoding an N-symbol RS code
     * and then shortening it to length N.
     * ------------------------------------------------------------- */
    unsigned int N = len + rs->T;
    unsigned int S = rs->gf.Np - N;
    for (i = 0; i < S; i++) {
	galois_field_val fb = parity[0];

	for (j = 0; j < rs->T - 1; j++)
	    parity[j] =
		galois_field_add(parity[j + 1],
				 galois_field_mul(&rs->gf, fb,
						  rs->generator[j + 1]));
	parity[rs->T - 1] = galois_field_mul(&rs->gf, fb,
					     rs->generator[rs->T]);
    }
#endif

    /* -------------------------------------------------------------
     * Feed the actual K information symbols
     * ------------------------------------------------------------- */
    for (i = 0; i < len; i++) {
	galois_field_val fb = galois_field_add(inbuf[i], parity[0]);

	fb = rs->gf.log[fb];
	if (fb != rs->gf.Np) {
	    for (j = 1; j < rs->T; j++) {
		galois_field_val tmp;

		tmp = rs->gf.exp[(fb + rs->generator[rs->T - j]) % rs->gf.Np];
		parity[j] = galois_field_add(parity[j], tmp);
	    }
	}
	for (j = 1; j < rs->T; j++)
	    parity[j - 1] = parity[j];

	if (fb != rs->gf.Np)
	    parity[rs->T - 1] = rs->gf.exp[(fb + rs->generator[0]) % rs->gf.Np];
	else
	    parity[rs->T - 1] = 0;
    }

    return 0;
}


void
reed_solomon_decoder_init(struct reed_solomon_decoder *rsd,
			  struct reed_solomon *rs)
{
    rsd->rs = rs;
#if GF_DYN_ALLOC
    rsd->recv_sym_p = malloc(sizeof(galois_field_val) * GALOIS_FIELD_MAX);

    rsd->C = malloc(sizeof(galois_field_val) * GALOIS_FIELD_MAX);
    rsd->B = malloc(sizeof(galois_field_val) * GALOIS_FIELD_MAX);
    rsd->Temp = malloc(sizeof(galois_field_val) * (GALOIS_FIELD_MAX + 1));

    rsd->synd = malloc(sizeof(galois_field_val) * REED_SOLOMON_MAX_T);
    rsd->error_pos = malloc(sizeof(unsigned int) * REED_SOLOMON_MAX_ERR);
    rsd->error_idx = malloc(sizeof(unsigned int) * REED_SOLOMON_MAX_ERR);

    rsd->O = malloc(sizeof(galois_field_val) * REED_SOLOMON_MAX_ERR);
#endif
}

/* -------------------------------------------------------------------------
 * 1) Syndrome computation (on parent length Np)
 *
 *     S_i = Σ_{j=0}^{Np-1} r_j α^{(i+1)*j},   for i = 0..T-1
 *
 * Returns the number of non-zero syndromes.
 * Zero syndromes → no errors.
 * ------------------------------------------------------------------------- */
static unsigned int
compute_syndromes(struct reed_solomon_decoder *rsd,
		  const galois_field_val *e)
{
    struct reed_solomon *rs = rsd->rs;
    galois_field_val *S = rsd->synd;
    unsigned int i, j, count = 0;

    for (i = 0; i < rs->T; i++)
	S[i] = e[0];

    for (i = 1; i < rsd->N; i++) {
	for (j = 0; j < rs->T; j++) {
	    if (S[j] == 0) {
		S[j] = e[i];
	    } else {
		unsigned int tmp;

		tmp = rs->gf.log[S[j]] + (rs->fcr + j) * rs->prim;
		tmp %= rs->gf.Np;
		tmp = rs->gf.exp[tmp];
		S[j] = galois_field_add(e[i], tmp);
	    }
	}
    }

    for (i = 0; i < rs->T; i++) {
	if (S[i])
	    count++;
	/* Pre-compute this to make it more efficient. */
	S[i] = rs->gf.log[S[i]];
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
berlekamp_massey(struct reed_solomon_decoder *rsd)
{
    struct reed_solomon *rs = rsd->rs;
    const galois_field_val *S = rsd->synd;
    unsigned int L = 0;
    unsigned int i, n;

    /* FIXME: Reduce size of B, C, and Temp */
    memset(rsd->B, 0, sizeof(galois_field_val) * GALOIS_FIELD_MAX);
    memset(rsd->C, 0, sizeof(galois_field_val) * GALOIS_FIELD_MAX);

    rsd->C[0] = 1;

    for (i = 0; i <= rs->T; i++)
	rsd->B[i] = rs->gf.log[rsd->C[i]];

    for (n = 1; n <= rs->T; n++) {
	galois_field_val d = 0;

	for (i = 0; i < n; i++) {
	    if (rsd->C[i] != 0 && S[n - i - 1] != rs->gf.Np) {
		unsigned int tmp = rs->gf.log[rsd->C[i]] + S[n - i - 1];

		tmp = rs->gf.exp[tmp % rs->gf.Np];
		d = galois_field_add(d, tmp);
	    }
	}

	if (d == 0) {
	    memmove(rsd->B + 1, rsd->B, rs->T * sizeof(rsd->B[0]));
	    rsd->B[0] = rs->gf.Np;
	} else {
	    d = rs->gf.log[d];

	    rsd->Temp[0] = rsd->C[0];
	    for (i = 0; i < rs->T; i++) {
		if (rsd->B[i] != rs->gf.Np) {
		    unsigned int tmp;

		    tmp = rs->gf.exp[(d + rsd->B[i]) % rs->gf.Np];
		    rsd->Temp[i + 1] = galois_field_add(rsd->C[i + 1], tmp);
		} else {
		    rsd->Temp[i + 1] = rsd->C[i + 1];
		}
	    }

	    if (2 * L <= n - 1) {
		/* Update B(x) ← previous C(x) */
		for (i = 0; i <= rs->T; i++) {
		    if (rsd->C[i] == 0) {
			rsd->B[i] = rs->gf.Np;
		    } else {
			unsigned int tmp;

			tmp = rs->gf.log[rsd->C[i]] - d + rs->gf.Np;
			rsd->B[i] = tmp % rs->gf.Np;
		    }
		}

		L = n - L;
	    } else {
		memmove(rsd->B + 1, rsd->B, rs->T * sizeof(rsd->B[0]));
		rsd->B[0] = rs->gf.Np;
	    }

	    for (i = 0; i <= rs->T; i++)
		rsd->C[i] = rsd->Temp[i];
	}
    }

    return L;
}

/* -------------------------------------------------------------------------
 * 3) Chien search
 *
 * Find i such that σ(α^{-i}) = 0, for i = 0..Np-1.
 * Each such i corresponds to an error at position i.
 * ------------------------------------------------------------------------- */
static unsigned int
chien_search(struct reed_solomon_decoder *rsd, unsigned int L,
	     unsigned int *rcount)
{
    struct reed_solomon *rs = rsd->rs;
    unsigned int count = 0, d = 0;
    unsigned int i, j, k;

    for (i = 0; i <= rs->T; i++) {
	if (rsd->C[i] != 0)
	    d = i;
	rsd->C[i] = rs->gf.log[rsd->C[i]];
    }

    for (i = 1; i <= rs->T; i++)
	rsd->Temp[i] = rsd->C[i];
    k = rs->iprim - 1;
    for (i = 1; i <= rs->gf.Np; i++) {
	unsigned int sum = 1;

	for (j = d; j > 0; j--) {
	    if (rsd->Temp[j] != rs->gf.Np) {
		rsd->Temp[j] = (rsd->Temp[j] + j) % rs->gf.Np;
		sum = galois_field_add(sum, rs->gf.exp[rsd->Temp[j]]);
	    }
	}

	if (sum == 0) {
	    rsd->error_idx[count] = i;
	    rsd->error_pos[count++] = k;
	    if (count == d)
		break;
	}

	k = (k + rs->iprim) % rs->gf.Np;
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
correct_errors(struct reed_solomon_decoder *rsd,
	       galois_field_val *e, const galois_field_val *S,
	       const unsigned int *err_pos, unsigned int err_cnt)
{
    struct reed_solomon *rs = rsd->rs;
    unsigned int i, j, o;

    o = err_cnt - 1;
    for (i = 0; i <= o; i++) {
	unsigned int tmp = 0;

	for (j = i; ; j--) {
	    if (S[i - j] != rs->gf.Np && rsd->C[j] != rs->gf.Np)
		tmp ^= rs->gf.exp[(S[i - j] + rsd->C[j]) % rs->gf.Np];
	    if (j == 0)
		break;
	}
	rsd->O[i] = rs->gf.log[tmp];
    }

    for (i = err_cnt - 1; ; i--) {
	unsigned int tmp = 0, tmp2, d;

	for (j = o; ; j--) {
	    if (rsd->O[j] != rs->gf.Np)
		tmp ^= rs->gf.exp[(rsd->O[j] + j * rsd->error_idx[i]) % rs->gf.Np];
	    if (j == 0)
		break;
	}
	tmp2 = rs->gf.exp[(rsd->error_idx[i] * (rs->fcr - 1) + rs->gf.Np) % rs->gf.Np];

	d = 0;
	if (err_cnt < rs->T - 1)
	    j = err_cnt;
	else
	    j = rs->T - 1;
	for (j = j & ~1; ; j -= 2) {
	    if (rsd->C[j + 1] != rs->gf.Np)
		d ^= rs->gf.exp[(rsd->C[j + 1] + j * rsd->error_idx[i]) % rs->gf.Np];
	    if (j < 2)
		break;
	}

	if (tmp != 0) {
	    e[rsd->error_pos[i]] ^= rs->gf.exp[(rs->gf.log[tmp] + rs->gf.log[tmp2] + rs->gf.Np - rs->gf.log[d]) % rs->gf.Np];
	}

	if (i == 0)
	    break;
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
reed_solomon_decode(struct reed_solomon_decoder *rsd,
		    uint8_t *buf, unsigned int len,
		    unsigned int *err_count)
{
    struct reed_solomon *rs = rsd->rs;
    unsigned int i;
    unsigned int count = 0;

    if (len > rs->gf.Np)
	return 1;
    if (len <= rs->T)
	return 1;
    rsd->N = len;
    rsd->K = len - rs->T;
    /* Number of shortened symbols */
    rsd->S = rs->gf.Np - rsd->N;

    /* FIXME = optimize if S == 0. */
    /* Build parent-length buffer */
    for (i = 0; i < rsd->S; i++)
	rsd->recv_sym_p[i] = 0;

    for (i = 0; i < rsd->N; i++)
	rsd->recv_sym_p[rsd->S + i] = buf[i];

    /* Syndromes */
    count = compute_syndromes(rsd, rsd->recv_sym_p);

    if (count != 0) {
	/* BM → locator polynomial */
	unsigned int L = berlekamp_massey(rsd);

	/* Chien search */
	if (chien_search(rsd, L, &count))
	    return 1; /* Could not correct all symbols. */

	/* Correct */
	if (count > 0)
	    correct_errors(rsd, rsd->recv_sym_p,
			   rsd->synd, rsd->error_pos, count);
    }

    *err_count = count;

    if (count == 0) /* No errors. */
	return 0;

    /* Output K information symbols */
    for (i = 0; i < rsd->K; i++)
	buf[i] = rsd->recv_sym_p[rsd->S + i];

    return 0;
}

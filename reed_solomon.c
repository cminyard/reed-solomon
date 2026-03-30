
#include <stdbool.h>
#include <string.h>
#include "reed_solomon.h"

int
reed_solomon_init(struct reed_solomon *rs, unsigned int m, unsigned int T)
{
    unsigned int i, j;
    galois_field_val tmp[GALOIS_FIELD_MAX];
    galois_field_val g0;
    galois_field_val inv_g0;
    int err = galois_field_init(&rs->gf, m);

    if (err)
	return err;

    if (T > REED_SOLOMON_MAX_T)
	return 1;
    if (T >= rs->gf.Np)
	return 1;

    rs->m = m;
    rs->T = T;

#if GF_DYN_ALLOC
    rs->generator = malloc(sizeof(galois_field_val) * GALOIS_FIELD_MAX);
#endif

    /* ---------------------------------------------------------------------
     * Generator polynomial construction (degree T)
     * g(x) = (x - α^0)(x - α^1)...(x - α^(T-1))
     * --------------------------------------------------------------------- */
    for (i = 0; i <= rs->T; i++)
	rs->generator[i] = 0;
    rs->generator[0] = 1;

    for (i = 0; i < rs->T; i++) {
	/* Copy existing coefficients */
	for (j = 0; j <= i; j++)
	    tmp[j] = rs->generator[j];

	rs->generator[i + 1] = 0;

	/* Perform polynomial multiplication by (x - α^i) */
	for (j = i + 1; j >= 1; j--) {
	    galois_field_val term;

	    if (j <= i)
		term = galois_field_mul(&rs->gf, tmp[j], rs->gf.exp[i]);
	    else
		term =0;
	    rs->generator[j] = galois_field_add(tmp[j - 1], term);
	}
	rs->generator[0] = galois_field_mul(&rs->gf,
					    tmp[0], rs->gf.exp[i]);
    }

    /* Normalize g(x) so that g[0] = 1 */
    g0 = rs->generator[0];
    inv_g0 = galois_field_inv(&rs->gf, g0);
    for (j = 0; j <= rs->T; j++)
	rs->generator[j] = galois_field_mul(&rs->gf, rs->generator[j], inv_g0);

    return 0;
}

void
reed_solomon_encoder_init(struct reed_solomon_encoder *rse,
			  struct reed_solomon *rs)
{
    rse->rs = rs;
#if GF_DYN_ALLOC
    rse->u = malloc(sizeof(galois_field_val) * GALOIS_FIELD_MAX);
    rse->parity = malloc(sizeof(galois_field_val) * REED_SOLOMON_MAX_T);
#endif
}

/**
 * @brief Systematic Reed–Solomon encoder.
 *
 * Produces a codeword of:
 *      [K info symbols][T parity symbols]
 */
int
reed_solomon_encode(struct reed_solomon_encoder *rse,
		    uint8_t *buf, unsigned int len)
{
    struct reed_solomon *rs = rse->rs;
    unsigned int i, j;

    if (len > rs->gf.Np)
	return 1;
    if (len <= rs->T)
	return 1;
    rse->N = len;
    rse->K = len - rs->T;
    /* Number of shortened symbols */
    rse->S = rs->gf.Np - rse->N;

    /* -------------------------------------------------------------
     * Initialize T parity registers to zero
     * ------------------------------------------------------------- */
    for (i = 0; i < rse->K; i++)
	rse->u[i] = buf[i];
    for (i = 0; i < rs->T; i++)
	rse->parity[i] = 0;

    /* -------------------------------------------------------------
     * Handle shortening:
     *   Shift S dummy symbols (all zeros) through the encoder.
     *
     * This produces the same result as encoding an N-symbol RS code
     * and then shortening it to length N.
     * ------------------------------------------------------------- */
    for (i = 0; i < rse->S; i++) {
	galois_field_val fb = rse->parity[0];

	for (j = 0; j < rs->T - 1; j++)
	    rse->parity[j] =
		galois_field_add(rse->parity[j + 1],
				 galois_field_mul(&rs->gf, fb,
						  rs->generator[j + 1]));
	rse->parity[rs->T - 1] = galois_field_mul(&rs->gf, fb,
						  rs->generator[rs->T]);
    }

    /* -------------------------------------------------------------
     * Feed the actual K information symbols
     * ------------------------------------------------------------- */
    for (i = 0; i < rse->K; i++) {
	galois_field_val fb = galois_field_add(buf[i], rse->parity[0]);

	for (j = 0; j < rs->T - 1; j++)
	    rse->parity[j] =
		galois_field_add(rse->parity[j + 1],
				 galois_field_mul(&rs->gf, fb,
						  rs->generator[j + 1]));
	rse->parity[rs->T - 1] = galois_field_mul(&rs->gf, fb,
						  rs->generator[rs->T]);
    }

    for (i = 0; i < rs->T; i++)
	buf[i + rse->K] = rse->parity[i];

    return 0;
}


void
reed_solomon_decoder_init(struct reed_solomon_decoder *rsd,
			  struct reed_solomon *rs)
{
    rsd->rs = rs;
#if GF_DYN_ALLOC
    unsigned int i;
    rsd->C = malloc(sizeof(galois_field_val) * GALOIS_FIELD_MAX);
    rsd->B = malloc(sizeof(galois_field_val) * GALOIS_FIELD_MAX);
    rsd->Temp = malloc(sizeof(galois_field_val) * GALOIS_FIELD_MAX);

    rsd->A = malloc(sizeof(galois_field_val *) * REED_SOLOMON_MAX_ERR);
    for (i = 0; i < REED_SOLOMON_MAX_ERR; i++)
	rsd->A[i] = malloc(sizeof(galois_field_val) * REED_SOLOMON_MAX_ERR);

    rsd->e = malloc(sizeof(galois_field_val) * REED_SOLOMON_MAX_ERR);
    rsd->recv_sym_p = malloc(sizeof(galois_field_val) * GALOIS_FIELD_MAX);

    rsd->synd = malloc(sizeof(galois_field_val) * REED_SOLOMON_MAX_T);
    rsd->sigma = malloc(sizeof(galois_field_val) * (REED_SOLOMON_MAX_ERR + 1));
    rsd->error_pos = malloc(sizeof(unsigned int) * REED_SOLOMON_MAX_ERR);
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
compute_syndromes(struct reed_solomon *rs,
		  const galois_field_val *recv_sym_p, galois_field_val *S)
{
    unsigned int i, j, count = 0;

    for (i = 0; i < rs->T; i++) {
	galois_field_val sum = 0;
	unsigned int si = i + 1; /* Evaluate at α^(i+1) */

	for (j = 0; j < rs->gf.Np; j++) {
	    galois_field_val k = (si * j) % rs->gf.Np;

	    sum = galois_field_add(sum,
				   galois_field_mul(&rs->gf,
						    recv_sym_p[j],
						    rs->gf.exp[k]));
	}
	S[i] = sum;
	if (sum)
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
berlekamp_massey(struct reed_solomon_decoder *rsd,
		 const galois_field_val *S, galois_field_val *sigma_out)
{
    struct reed_solomon *rs = rsd->rs;
    unsigned int t = rs->T / 2;
    unsigned int L = 0;
    unsigned int m_shift = 1;
    galois_field_val bbb = 1;
    unsigned int i, n;

    memset(rsd->B, 0, sizeof(galois_field_val) * GALOIS_FIELD_MAX);
    memset(rsd->C, 0, sizeof(galois_field_val) * GALOIS_FIELD_MAX);

    rsd->C[0] = 1;
    rsd->B[0] = 1;

    for (n = 0; n < rs->T - 1; n++) {
	galois_field_val d = S[n];

	for (i = 1; i <= L; i++)
	    d = galois_field_add(d,
				 galois_field_mul(&rs->gf,
						  rsd->C[i], S[n - i]));
	if (d != 0) {
	    galois_field_val coef;

	    for (i = 0; i <= rs->T; i++)
		rsd->Temp[i] = rsd->C[i];

	    coef = galois_field_div(&rs->gf, d, bbb);

	    /* C(x) ← C(x) - coef * x^m_shift * B(x) */
	    for (i = 0; i <= rs->T; i++) {
		int idx = i + m_shift;

		if (idx <= rs->T)
		    rsd->C[idx] = galois_field_add(rsd->C[idx],
						   galois_field_mul(&rs->gf,
							    coef, rsd->B[i]));
	    }

	    if (2 * L <= n) {
		/* Update B(x) ← previous C(x) */
		for (i = 0; i <= rs->T; i++)
		    rsd->B[i] = rsd->Temp[i];
		L = n + 1 - L;
		bbb = d;
		m_shift = 1;
	    } else {
		m_shift++;
	    }
	} else {
	    m_shift++;
	}
    }

    /* Copy result */
    for (i = 0; i <= t && i <= L; i++)
	sigma_out[i] = rsd->C[i];

    /* Ensure σ(0) = 1 */
    if (sigma_out[0] == 0)
	sigma_out[0] = 1;

    return L;
}

/* -------------------------------------------------------------------------
 * 3) Chien search
 *
 * Find i such that σ(α^{-i}) = 0, for i = 0..Np-1.
 * Each such i corresponds to an error at position i.
 * ------------------------------------------------------------------------- */
static unsigned int
chien_search(struct reed_solomon *rs,
	     const galois_field_val *sigma, unsigned int L,
	     unsigned int *error_pos)
{
    unsigned int count = 0;
    unsigned int i, j;

    for (i = 0; i < rs->gf.Np; i++) {
	galois_field_val x_inv;
	galois_field_val sum = 0;
	galois_field_val power = 1;

	if (i == 0)
	    x_inv = 1;
	else
	    x_inv = rs->gf.exp[rs->gf.Np - i];

	for (j = 0; j <= L; j++) {
	    if (sigma[j] != 0)
		sum = galois_field_add(sum,
				       galois_field_mul(&rs->gf,
							sigma[j], power));
	    power = galois_field_mul(&rs->gf, power, x_inv);
	}

	if (sum == 0)
	    error_pos[count++] = i;

	if (count > L)
	    break;
    }

    return count;
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
	       galois_field_val *recv_sym_p, const galois_field_val *S,
	       const unsigned int *error_pos, unsigned int error_count)
{
    struct reed_solomon *rs = rsd->rs;
    unsigned int i, r, c;

    if (error_count == 0)
	return;

    /* Construct linear system */
    for (r = 0; r < error_count; r++) {
	rsd->B[r] = S[r];
	for (c = 0; c < error_count; c++) {
	    int pos = error_pos[c];
	    int exp = ((r + 1) * pos) % rs->gf.Np;

	    rsd->A[r][c] = rs->gf.exp[exp];
	}
    }

    /* Gaussian elimination over GF(2^m) */
    for (i = 0; i < error_count; i++) {
	galois_field_val piv = rsd->A[i][i];
	galois_field_val inv;

	if (piv == 0) {
	    int swap_r = -1;
	    for (int r = i + 1; r < error_count; r++)
		if (rsd->A[r][i] != 0) {
		    swap_r = r;
		    break;
		}

	    if (swap_r >= 0) {
		galois_field_val tmp;

		for (c = 0; c < error_count; c++) {
		    tmp = rsd->A[i][c];
		    rsd->A[i][c] = rsd->A[swap_r][c];
		    rsd->A[swap_r][c] = tmp;
		}

		tmp = rsd->B[i];
		rsd->B[i] = rsd->B[swap_r];
		rsd->B[swap_r] = tmp;
		piv = rsd->A[i][i];
	    }
	}

	if (piv == 0)
	    continue;

	inv = galois_field_inv(&rs->gf, piv);
	for (c = 0; c < error_count; c++)
	    rsd->A[i][c] = galois_field_mul(&rs->gf, rsd->A[i][c], inv);
	rsd->B[i] = galois_field_mul(&rs->gf, rsd->B[i], inv);

	for (r = 0; r < error_count; r++) {
	    galois_field_val factor;

	    if (r == i)
		continue;
	    factor = rsd->A[r][i];
	    if (factor == 0)
		continue;

	    for (c = 0; c < error_count; c++)
		rsd->A[r][c] =
		    galois_field_add(rsd->A[r][c],
				     galois_field_mul(&rs->gf, factor,
						      rsd->A[i][c]));
	    rsd->B[r] = galois_field_add(rsd->B[r],
					 galois_field_mul(&rs->gf, factor,
							  rsd->B[i]));
	}
    }

    for (i = 0; i < error_count; i++)
	rsd->e[i] = rsd->B[i];

    /* Apply error corrections */
    for (i = 0; i < error_count; i++) {
	int pos = error_pos[i];

	recv_sym_p[pos] ^= rsd->e[i];
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
    unsigned int t = rs->T / 2;
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

    /* Build parent-length buffer */
    for (i = 0; i < rsd->S; i++)
	rsd->recv_sym_p[i] = 0;

    for (i = 0; i < rsd->N; i++)
	rsd->recv_sym_p[rsd->S + i] = buf[i];

    /* Syndromes */
    count = compute_syndromes(rs, rsd->recv_sym_p, rsd->synd);

    if (count != 0) {
	/* BM → locator polynomial */
	unsigned int L = berlekamp_massey(rsd, rsd->synd, rsd->sigma);

	if (L > t) {
	    count = L;
	} else {
	    /* Chien search */
	    count = chien_search(rs, rsd->sigma, L, rsd->error_pos);

	    /* Correct */
	    if (count > 0 && count <= t)
		correct_errors(rsd, rsd->recv_sym_p,
			       rsd->synd, rsd->error_pos, count);
	}
    }

    *err_count = count;

    if (count == 0) /* No errors. */
	return 0;

    /* Output K information symbols */
    for (i = 0; i < rsd->K; i++)
	buf[i] = rsd->recv_sym_p[rsd->S + i];

    return 0;
}

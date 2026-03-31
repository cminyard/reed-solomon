
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

    return 0;
}

void
reed_solomon_encoder_init(struct reed_solomon_encoder *rse,
			  struct reed_solomon *rs)
{
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
    for (i = 1; i <= rs->T; i++)
	rse->generator[i] = 0;

    root = rs->fcr * rs->prim;
    for (i = 0; i < rs->T; i++, root += rs->prim) {
	/* Copy existing coefficients */
	for (j = 0; j <= i; j++)
	    tmp[j] = rse->generator[j];

	rse->generator[i + 1] = 1;

	/* Perform polynomial multiplication by (x - α^i) */
	for (j = i; j > 0; j--) {
	    if (rse->generator[j] != 0) {
		unsigned int term;

		term = (rs->gf.log[tmp[j]] + root) % rs->gf.Np;
		term = rs->gf.exp[term];
		rse->generator[j] = gf_add(tmp[j - 1], term);
	    } else {
		rse->generator[j] = rse->generator[j - 1];
	    }
	}
	rse->generator[0] = rs->gf.exp[(rs->gf.log[tmp[0]] + root) % rs->gf.Np];
    }

    /*
     * Pre-compute the log so it doesn't have to be done every time
     * when encoding.
     */
    for (i = 0; i <= rs->T; i++)
	rse->generator[i] = rs->gf.log[rse->generator[i]];
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
	gf_val fb = parity[0];

	for (j = 0; j < rs->T - 1; j++)
	    parity[j] =
		gf_add(parity[j + 1],
		       gf_mul(&rs->gf, fb, se->generator[j + 1]));
	parity[rs->T - 1] = gf_mul(&rs->gf, fb,
					     rse->generator[rs->T]);
    }
#endif

    /* -------------------------------------------------------------
     * Feed the actual K information symbols
     * ------------------------------------------------------------- */
    for (i = 0; i < len; i++) {
	gf_val fb = gf_add(inbuf[i], parity[0]);

	fb = rs->gf.log[fb];
	if (fb != rs->gf.Np) {
	    for (j = 1; j < rs->T; j++) {
		gf_val tmp;

		tmp = rs->gf.exp[(fb + rse->generator[rs->T - j]) % rs->gf.Np];
		parity[j] = gf_add(parity[j], tmp);
	    }
	}
	for (j = 1; j < rs->T; j++)
	    parity[j - 1] = parity[j];

	if (fb != rs->gf.Np)
	    parity[rs->T - 1] = rs->gf.exp[(fb + rse->generator[0]) % rs->gf.Np];
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
    rsd->recv_sym_p = malloc(sizeof(gf_val) * GF_MAX);

    rsd->C = malloc(sizeof(gf_val) * (REED_SOLOMON_MAX_T + 1));
    rsd->B = malloc(sizeof(gf_val) * (REED_SOLOMON_MAX_T + 1));
    rsd->Temp = malloc(sizeof(gf_val) * (REED_SOLOMON_MAX_T + 1));

    rsd->synd = malloc(sizeof(gf_val) * REED_SOLOMON_MAX_T);
    rsd->error_pos = malloc(sizeof(unsigned int) * REED_SOLOMON_MAX_ERR);
    rsd->error_idx = malloc(sizeof(unsigned int) * REED_SOLOMON_MAX_ERR);

    rsd->O = malloc(sizeof(gf_val) * REED_SOLOMON_MAX_ERR);
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
	for (j = 0; j < rs->T; j++) {
	    if (S[j] == 0) {
		S[j] = e[i];
	    } else {
		unsigned int tmp;

		tmp = gf->log[S[j]] + (rs->fcr + j) * rs->prim;
		tmp %= gf->Np;
		tmp = gf->exp[tmp];
		S[j] = gf_add(e[i], tmp);
	    }
	}
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

    memset(B, 0, sizeof(gf_val) * (REED_SOLOMON_MAX_T + 1));
    memset(C, 0, sizeof(gf_val) * (REED_SOLOMON_MAX_T + 1));

    C[0] = 1;

    for (i = 0; i <= rs->T; i++)
	B[i] = gf->log[C[i]];

    for (n = 1; n <= rs->T; n++) {
	gf_val d = 0;

	for (i = 0; i < n; i++) {
	    if (C[i] != 0 && S[n - i - 1] != 0) {
		gf_val t = gf_mul(gf, C[i], S[n - i - 1]);

		d = gf_add(d, t);
	    }
	}

	if (d == 0) {
	    memmove(B + 1, B, rs->T * sizeof(B[0]));
	    B[0] = gf->Np;
	} else {
	    d = gf->log[d];

	    Temp[0] = C[0];
	    for (i = 0; i < rs->T; i++) {
		if (B[i] != gf->Np) {
		    unsigned int tmp;

		    tmp = gf->exp[(d + B[i]) % gf->Np];
		    Temp[i + 1] = gf_add(C[i + 1], tmp);
		} else {
		    Temp[i + 1] = C[i + 1];
		}
	    }

	    if (2 * L <= n - 1) {
		/* Update B(x) ← previous C(x) */
		for (i = 0; i <= rs->T; i++) {
		    if (C[i] == 0) {
			B[i] = gf->Np;
		    } else {
			B[i] = (gf->log[C[i]] - d + gf->Np) % gf->Np;
		    }
		}

		L = n - L;
	    } else {
		memmove(B + 1, B, rs->T * sizeof(B[0]));
		B[0] = gf->Np;
	    }

	    for (i = 0; i <= rs->T; i++)
		C[i] = Temp[i];
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
	C[i] = gf->log[C[i]];
    }

    for (i = 1; i <= rs->T; i++)
	Temp[i] = C[i];
    k = rs->iprim - 1;
    for (i = 1; i <= gf->Np; i++) {
	unsigned int sum = 1;

	for (j = d; j > 0; j--) {
	    if (Temp[j] != gf->Np) {
		Temp[j] = (Temp[j] + j) % gf->Np;
		sum = gf_add(sum, gf->exp[Temp[j]]);
	    }
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
	unsigned int tmp = 0;

	for (j = i; ; j--) {
	    if (S[i - j] != 0 && C[j] != gf->Np)
		tmp ^= gf->exp[(gf->log[S[i - j]] + C[j]) % gf->Np];
	    if (j == 0)
		break;
	}
	O[i] = gf->log[tmp];
    }

    for (i = err_cnt - 1; ; i--) {
	unsigned int tmp = 0, tmp2, d;

	for (j = o; ; j--) {
	    if (O[j] != gf->Np)
		tmp ^= gf->exp[(O[j] + j * err_idx[i]) % gf->Np];
	    if (j == 0)
		break;
	}
	tmp2 = gf->exp[(err_idx[i] * (rs->fcr - 1) + gf->Np) % gf->Np];

	d = 0;
	if (err_cnt < rs->T - 1)
	    j = err_cnt;
	else
	    j = rs->T - 1;
	for (j = j & ~1; ; j -= 2) {
	    if (C[j + 1] != gf->Np)
		d ^= gf->exp[(C[j + 1] + j * err_idx[i]) % gf->Np];
	    if (j < 2)
		break;
	}

	if (tmp != 0) {
	    e[err_pos[i]] ^= gf->exp[(gf->log[tmp] + gf->log[tmp2] + gf->Np - gf->log[d]) % gf->Np];
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
    count = compute_syndromes(rs, rsd->N, rsd->recv_sym_p, rsd->synd);

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
	    correct_errors(rs, rsd->recv_sym_p, rsd->synd, rsd->C, rsd->O,
			   rsd->error_idx, rsd->error_pos, count);
    }

    *err_count = count;

    if (count == 0) /* No errors. */
	return 0;

    /* Output K information symbols */
    for (i = 0; i < rsd->K; i++)
	buf[i] = rsd->recv_sym_p[rsd->S + i];

    return 0;
}

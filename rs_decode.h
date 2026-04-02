/* -------------------------------------------------------------------------
 * 1) Syndrome computation (on parent length Np)
 *
 *     S_i = Σ_{j=0}^{Np-1} r_j α^{(i+1)*j},   for i = 0..T-1
 *
 * Returns the number of non-zero syndromes.
 * Zero syndromes = no errors.
 * ------------------------------------------------------------------------- */
static GF_FORCE_INLINE unsigned int
RS_NAME(compute_syndromes)(struct reed_solomon *rs, unsigned int pad,
			   const gf_sym *e, gf_sym *S)
{
    unsigned int i, j, count = 0;

    /* No need to process the pad, the result will be all zeros. */

    for (i = 0; i < RS_T; i++)
	S[i] = e[0];

    for (i = pad + 1; i < GF_NP; i++) {
	for (j = 0; j < RS_T; j++) {
	    if (S[j] == 0) {
		S[j] = e[i - pad];
	    } else {
		/* S[j] = e[i] + (S[j] * (rs->fcr + j) ^ rs->prim) */
		S[j] = e[i - pad] ^ GF_EXP((GF_LOG(S[j])
					+ (rs->fcr + j) * rs->prim) % GF_NP);
		/*
		 * Direct calculation above is much faster than the below
		 * because it doesn't do all the error checks.
		 *
		 * S[j] = gf_add(e[i = pad],
		 *	      gf_mul_el(gf, S[j],
		 *			gf_pow_l_l(gf, rs->fcr + j, rs->prim)));
		 */
	    }
	}
    }

    for (i = 0; i < RS_T; i++) {
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
static GF_FORCE_INLINE unsigned int
RS_NAME(berlekamp_massey)(struct reed_solomon *rs,
			  gf_sym *S,
			  gf_sym *B, gf_sym *C,
			  gf_sym *Temp)
{
    unsigned int L = 0;
    unsigned int i, n;

    for (i = 0; i < RS_T; i++)
	/*
	 * Put S into log format so it doesn't have to be recomputed
	 * when used.  It doesn't matter if S[i] is zero upon entry.
	 */
	S[i] = GF_LOG(S[i]);

    /*
     * Instead of copying the B array, move it backwards.  This can
     * happen at most RS_MAX_T times, so start out there.
     */
    B += RS_MAX_T;

    memset(C + 1, 0, sizeof(gf_sym) * RS_T);
    C[0] = 1;

    /* B array is in log format. */
    for (i = 0; i <= RS_T; i++)
	B[i] = GF_NP;
    B[0] = GF_LOG(1);

    for (n = 1; n <= RS_T; n++) {
	gf_sym d = 0;

	for (i = 0; i < n; i++)
	    /* d += C[i] * S[n - i - 1] */
	    d = gf_add(d, GF_MUL_EL(C[i], S[n - i - 1]));

	if (d == 0) {
	    B--;
	    B[0] = GF_NP;
	} else {
	    d = GF_LOG(d);
	    Temp[0] = C[0];
	    for (i = 0; i < RS_T; i++)
		Temp[i + 1] = gf_add(C[i + 1], GF_MUL_LL(d, B[i]));

	    if (2 * L <= n - 1) {
		for (i = 0; i <= RS_T; i++)
		    B[i] = GF_DIV_EL_L(C[i], d);

		L = n - L;
	    } else {
		B--;
		B[0] = GF_NP;
	    }

	    memcpy(C, Temp, RS_T * sizeof(gf_sym));
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
static GF_FORCE_INLINE unsigned int
RS_NAME(chien_search)(struct reed_solomon *rs, unsigned int L,
		      gf_sym *C, gf_sym *Temp,
		      unsigned int *err_idx, unsigned int *err_pos,
		      unsigned int *rcount)
{
    unsigned int count = 0, d = 0;
    unsigned int i, j, k;

    for (i = 0; i <= RS_T; i++) {
	if (C[i] != 0)
	    d = i;
	/* Convert C[i] to log format for faster processing. */
	C[i] = GF_LOG(C[i]);
	Temp[i] = C[i];
    }

    k = rs->iprim - 1;
    for (i = 1; i <= GF_NP; i++) {
	unsigned int sum = 1;

	for (j = d; j > 0; j--) {
	    Temp[j] = GF_MUL_LL_L(Temp[j], j);
	    sum = gf_add(sum, GF_EXP(Temp[j]));
	}

	if (sum == 0) {
	    err_idx[count] = i;
	    err_pos[count++] = k;
	    if (count == d)
		break;
	}

	k = (k + rs->iprim) % GF_NP;
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
static GF_FORCE_INLINE void
RS_NAME(correct_errors)(struct reed_solomon *rs,
			gf_sym *e, unsigned int pad, const gf_sym *S,
			const gf_sym *C, gf_sym *O,
			const unsigned int *err_idx,
			const unsigned int *err_pos,
			unsigned int err_cnt)
{
    unsigned int i, j, o;

    o = err_cnt - 1;
    for (i = 0; i <= o; i++) {
	O[i] = 0;
	for (j = 0; j <= i; j++)
	    /* O[i] += S[i - j] * C[j] */
	    O[i] = gf_add(O[i], GF_MUL_LL(S[i - j], C[j]));
	/* Convert O[i] to log format for faster processing. */
	O[i] = GF_LOG(O[i]);
    }

    for (i = 0; i < err_cnt; i++) {
	unsigned int tmp = 0, tmp2, d;

	for (j = 0; j <= o; j++)
	    /* tmp += O[j] * j ^ err_idx[i] */
	    tmp = gf_add(tmp, GF_MUL_LL(O[j], GF_POW_L_L(j, err_idx[i])));
	/* tmp2 = err_index[i] ^ (fcr - 1) */
	tmp2 = GF_POW_L(err_idx[i], rs->fcr - 1);

	d = 0;
	if (err_cnt < RS_T - 1)
	    j = err_cnt;
	else
	    j = RS_T - 1;
	for (j = j & ~1; ; j -= 2) {
	    /* d += C[j + 1] * j ^ err_idx[i] */
	    d = gf_add(d, GF_MUL_LL(C[j + 1],
				    GF_POW_L_L(j, err_idx[i])));
	    if (j < 2)
		break;
	}

	if (err_pos[i] >= pad)
	    e[err_pos[i] - pad] = gf_add(e[err_pos[i] - pad],
					 GF_DIV(GF_MUL(tmp, tmp2), d));
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
RS_DEC_START()
{
    unsigned int count = 0, pad;

    if (len > GF_NP)
	return 1;
    if (len <= RS_T)
	return 1;

    pad = GF_NP - len;

    /* Syndromes */
    count = RS_NAME(compute_syndromes)(rs, pad, data, rsd->synd);

    if (count == 0) {
	*err_count = 0;
	return 0;
    }

    /* BM → locator polynomial */
    unsigned int L = RS_NAME(berlekamp_massey)(rs, rsd->synd, rsd->B, rsd->C,
					       rsd->Temp);

    /* Chien search */
    if (RS_NAME(chien_search)(rs, L, rsd->C, rsd->Temp,
			      rsd->error_idx, rsd->error_pos, &count))
	return 1; /* Could not correct all symbols. */

    /* Correct */
    if (count > 0)
	RS_NAME(correct_errors)(rs, data, pad, rsd->synd, rsd->C, rsd->O,
				rsd->error_idx, rsd->error_pos, count);

    *err_count = count;

    return 0;
}
RS_DEC_END

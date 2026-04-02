/**
 * @brief Systematic Reed–Solomon encoder.
 *
 * Produces a codeword of:
 *      [K info symbols][T parity symbols]
 */

int
rs_encode_8(uint8_t *inbuf, unsigned int len, uint8_t *parity)
{
    unsigned int i, j, k;
    static gf_v16ss zero = { 0 };
    static gf_v16ss shift_mask = { 1, 2, 3, 4, 5, 6, 7, 8 };
    static gf_v16ss vecnp = {
	GF_NP, GF_NP, GF_NP, GF_NP,
	GF_NP, GF_NP, GF_NP, GF_NP
    };
    gf_v16ss p[4] = { 0, 0, 0, 0 };

    if (len > GF_NP - RS_T)
	return 1;

    /* No need to process the pad, the result will still be 0. */

    /* Feed the actual data. */
    for (i = 0; i < len; i++) {
	uint8_t fb1 = GF_LOG(gf_add(inbuf[i], p[0][0]));
 
        if (fb1 != GF_NP) {
	    gf_v16ss fb = zero + fb1, tmp, cmp;

	    for (j = 0; j < 4; j++) {
		/* Do the galois field multiply. */
		tmp = fb + SIMD_GEN[j];

		/*
		 * Now do the modulo. cmp will be -1 if they are not
		 * the same.  So multiply that by Np and add it to
		 * subtract Np from every element that is >= Np.
		 */
		cmp = tmp >= vecnp;
		cmp *= vecnp;
		tmp += cmp;

		/* Convert from log to exp domain. */
		for (k = 0; k < 8; k++)
		    tmp[k] = GF_EXP(tmp[k]);

		/* Do the galois field add. */
		p[j] ^= tmp;
	    }
	}

	/*
	 * This shifts all the data to the left and adds the new
	 * element ot the end.
	 */
	p[0] = __builtin_shuffle(p[0], p[1], shift_mask);
	p[1] = __builtin_shuffle(p[1], p[2], shift_mask);
	p[2] = __builtin_shuffle(p[2], p[3], shift_mask);
	p[3] = __builtin_shuffle(p[3], p[0], shift_mask);
	p[3][7] = GF_MUL_LL(fb1, RS_GENERATOR[0]);
    }

    for (i = 0; i < RS_T; ) {
	for (j = 0; j < 4; j++) {
	    for (k = 0; k < 8; k++)
		parity[i++] = p[j][k];
	}
    }
	
    return 0;
}

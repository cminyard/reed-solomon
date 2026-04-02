/**
 * @brief Systematic Reed–Solomon encoder.
 *
 * Produces a codeword of:
 *      [K info symbols][T parity symbols]
 */

RS_ENC_START()
{
    unsigned int i, j, k;
    /*
     * Does a shift to the left with __builtin_shuffle.  See gcc
     * docs for details.
     */
    static gf_v16ss shift_mask = { 1, 2, 3, 4, 5, 6, 7, 7 };
    static gf_v16ss zerov = { 0 };
    /* Parity array, as 4 vectors of 8 elements. */
    gf_v16ss p[SIMD_LEN];

    if (len > GF_NP - RS_T)
	return 1;

    if (!CAN_DO_SIMD)
	return SIMD_FALLBACK;

    for (i = 0; i < SIMD_LEN; i++)
	p[i] = zerov;

    /* No need to process the pad, the result will still be 0. */

    /* Feed the actual data. */
    for (i = 0; i < len; i++) {
	uint8_t fb = GF_LOG(gf_add(inbuf[i], p[0][0]));
 
        if (fb != GF_NP) {
	    gf_v16ss fbv = { fb, fb, fb, fb, fb, fb, fb, fb };
	    gf_v16ss tmp, cmp;

	    for (j = 0; j < SIMD_LEN; j++) {
		/* Do the galois field multiply. */
		tmp = fbv + SIMD_GEN[j];

		/*
		 * Now do the modulo.  After the add above elements
		 * may be >= Np, so subtract Np from those elements to
		 * do the modulo.
		 */
		cmp = tmp >= VECNP;
		/* cmp will be -1 for each element >= Np, 0 for others. */
		cmp *= VECNP;
		/* cmp will be -Np for each element >= Np, 0 for others. */
		tmp += cmp;
		/* Each element >= Np will have Np subtracted from it now. */

		/* Convert from log to exp domain. */
		for (k = 0; k < 8; k++)
		    tmp[k] = GF_EXP(tmp[k]);

		/* Do the galois field add. */
		p[j] ^= tmp;
	    }
	}

	/*
	 * This shifts all the data to the left and adds the new
	 * element at the end.
	 */
	p[0] = __builtin_shuffle(p[0], shift_mask);
	p[0][7] = p[1][0];
	p[1] = __builtin_shuffle(p[1], shift_mask);
	p[1][7] = p[2][0];
	p[2] = __builtin_shuffle(p[2], shift_mask);
	p[2][7] = p[3][0];
	p[3] = __builtin_shuffle(p[3], shift_mask);
	p[3][7] = GF_MUL_LL(fb, RS_GENERATOR[0]);
    }

    for (i = 0, j = 0; j < SIMD_LEN; j++) {
	for (k = 0; k < 8; k++)
	    parity[i++] = p[j][k];
    }
	
    return 0;
}
RS_ENC_END

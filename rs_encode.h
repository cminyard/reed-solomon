/**
 * @brief Systematic Reed–Solomon encoder.
 *
 * Produces a codeword of:
 *      [K info symbols][T parity symbols]
 */
RS_ENC_START()
{
    unsigned int i, j;

    if (RS_LEN > GF_NP - RS_T)
	return 1;

    memset(parity, 0, RS_T * sizeof(gf_sym));

    /* No need to process the pad, the result will still be 0. */

    /* Feed the actual data. */
    for (i = 0; i < RS_LEN; i++) {
	gf_sym fb = GF_LOG[gf_add(inbuf[i], parity[0])];

	if (fb != GF_NP) {
	    for (j = 1; j < RS_T; j++)
		/* parity[j] += fb * rse->generator[RS_T - j] */
#if RS_USE_GF
		parity[j] = gf_add(parity[j],
				   gf_mul_ll(gf, fb, RS_GENERATOR[RS_T - j]));
#else
		parity[j] = gf_add(parity[j],
				   GF_EXP[do_mod(fb + RS_GENERATOR[RS_T - j])]);
#endif
	}
	memmove(parity, parity + 1, (RS_T - 1) * sizeof(gf_sym));

	if (fb != GF_NP) {
#if RS_USE_GF
	    parity[RS_T - 1] = gf_mul_ll(gf, fb, RS_GENERATOR[0]);
#else
	    parity[RS_T - 1] = GF_EXP[do_mod(fb + RS_GENERATOR[0])];
#endif
	} else {
	    parity[RS_T - 1] = 0;
	}
    }

    return 0;
}
RS_ENC_END

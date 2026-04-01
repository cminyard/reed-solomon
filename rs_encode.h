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
	gf_sym fb = GF_LOG(gf_add(inbuf[i], parity[0]));

	if (fb != GF_NP) {
	    for (j = 1; j < RS_T; j++)
		/* parity[j] += fb * rse->generator[RS_T - j] */
		parity[j] = gf_add(parity[j],
				   GF_MUL_LL(gf, fb, RS_GENERATOR[RS_T - j]));
	}
	memmove(parity, parity + 1, (RS_T - 1) * sizeof(gf_sym));

	if (fb != GF_NP) {
	    parity[RS_T - 1] = GF_MUL_LL(gf, fb, RS_GENERATOR[0]);
	} else {
	    parity[RS_T - 1] = 0;
	}
    }

    return 0;
}
RS_ENC_END

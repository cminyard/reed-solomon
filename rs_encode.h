/**
 * @brief Systematic Reed–Solomon encoder.
 *
 * Produces a codeword of:
 *      [K info symbols][T parity symbols]
 */
RS_ENC_START()
{
    unsigned int i, j, p = 0, c;

    /*
     * This uses the parity array as a circular array of symbols, "p"
     * is the current starting position of the array.  On a x86 64-bit
     * processor with a RS(223,255) coder this is not really any
     * faster than just copying, but it's got to be faster on less
     * capable processors.
     */
#define NEXTP(x) ((x) == (RS_T - 1) ? 0 : (x) + 1)

    if (RS_LEN > GF_NP - RS_T)
	return 1;

    memset(parity, 0, RS_T * sizeof(gf_sym));

    /* No need to process the pad, the result will still be 0. */

    /* Feed the actual data. */
    for (i = 0; i < RS_LEN; i++) {
	gf_sym fb = GF_LOG(gf_add(inbuf[i], parity[p]));

	if (fb != GF_NP) {
	    for (c = 1, j = NEXTP(p); c < RS_T; c++, j = NEXTP(j))
		/* parity[j] += fb * rse->generator[RS_T - c] */
		parity[j] = gf_add(parity[j],
				   GF_MUL_LL(fb, RS_GENERATOR[RS_T - c]));
	}

	parity[p] = GF_MUL_LL(fb, RS_GENERATOR[0]);
	p = NEXTP(p);
    }

    /* Shift the array around to the right place. */
    gf_shift(parity, RS_T, p);

    return 0;
}
RS_ENC_END

#undef NEXTP

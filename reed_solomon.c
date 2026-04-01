
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
    gf_sym tmp[GF_MAX];
    unsigned int root;

    rse->rs = rs;

#if GF_DYN_ALLOC
    rse->generator = malloc(sizeof(gf_sym) * GF_MAX);
#endif

    /* ---------------------------------------------------------------------
     * Generator polynomial construction (degree T)
     * g(x) = (x - α^0)(x - α^1)...(x - α^(T-1))
     * --------------------------------------------------------------------- */
    rse->generator[0] = 1;
    memset(rse->generator + 1, 0, rs->T * sizeof(gf_sym));

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

    /* Put the generators into log format to speed calculations. */
    for (i = 0; i < rs->T; i++)
	rse->generator[i] = gf_log_nc(gf, rse->generator[i]);
}


#define RS_T rs->T
#define GF_NP gf->Np
#define GF_EXP gf->exp
#define GF_LOG gf->log
#define RS_LEN len
#define RS_USE_GF 1
#define RS_GENERATOR rse->generator
#define RS_ENC_START() \
int \
rs_encode(struct reed_solomon_encoder *rse, \
	  uint8_t *inbuf, unsigned int len, uint8_t *parity) \
{ \
    struct reed_solomon *rs = rse->rs; \
    struct galois_field *gf = &rs->gf;
#define RS_ENC_END \
}
#include "rs_encode.h"

#undef RS_T
#undef GF_NP
#undef GF_EXP
#undef GF_LOG
#undef RS_LEN
#undef RS_USE_GF
#undef RS_GENERATOR
#undef RS_ENC_START
#undef RS_ENC_END

#define RS_T 32
#define GF_NP 255
#define GF_EXP CCSDS_exp
#define GF_LOG CCSDS_log
#define RS_GENERATOR CCSDS_gen
#define RS_LEN len
#define RS_USE_GF 0
#define RS_ENC_START() \
int \
rs_encode_8(uint8_t *inbuf, unsigned int len, uint8_t *parity)	\
{
#define RS_ENC_END \
}
const uint8_t CCSDS_exp[] = {
      1,   2,   4,   8,  16,  32,  64, 128,
    135, 137, 149, 173, 221,  61, 122, 244,
    111, 222,  59, 118, 236,  95, 190, 251,
    113, 226,  67, 134, 139, 145, 165, 205,
     29,  58, 116, 232,  87, 174, 219,  49,
     98, 196,  15,  30,  60, 120, 240, 103,
    206,  27,  54, 108, 216,  55, 110, 220,
     63, 126, 252, 127, 254, 123, 246, 107,
    214,  43,  86, 172, 223,  57, 114, 228,
     79, 158, 187, 241, 101, 202,  19,  38,
     76, 152, 183, 233,  85, 170, 211,  33,
     66, 132, 143, 153, 181, 237,  93, 186,
    243,  97, 194,   3,   6,  12,  24,  48,
     96, 192,   7,  14,  28,  56, 112, 224,
     71, 142, 155, 177, 229,  77, 154, 179,
    225,  69, 138, 147, 161, 197,  13,  26,
     52, 104, 208,  39,  78, 156, 191, 249,
    117, 234,  83, 166, 203,  17,  34,  68,
    136, 151, 169, 213,  45,  90, 180, 239,
     89, 178, 227,  65, 130, 131, 129, 133,
    141, 157, 189, 253, 125, 250, 115, 230,
     75, 150, 171, 209,  37,  74, 148, 175,
    217,  53, 106, 212,  47,  94, 188, 255,
    121, 242,  99, 198,  11,  22,  44,  88,
    176, 231,  73, 146, 163, 193,   5,  10,
     20,  40,  80, 160, 199,   9,  18,  36,
     72, 144, 167, 201,  21,  42,  84, 168,
    215,  41,  82, 164, 207,  25,  50, 100,
    200,  23,  46,  92, 184, 247, 105, 210,
     35,  70, 140, 159, 185, 245, 109, 218,
     51, 102, 204,  31,  62, 124, 248, 119,
    238,  91, 182, 235,  81, 162, 195,   0,
};

const uint8_t CCSDS_log[] = {
    255,   0,   1,  99,   2, 198, 100, 106,
      3, 205, 199, 188, 101, 126, 107,  42,
      4, 141, 206,  78, 200, 212, 189, 225,
    102, 221, 127,  49, 108,  32,  43, 243,
      5,  87, 142, 232, 207, 172,  79, 131,
    201, 217, 213,  65, 190, 148, 226, 180,
    103,  39, 222, 240, 128, 177,  50,  53,
    109,  69,  33,  18,  44,  13, 244,  56,
      6, 155,  88,  26, 143, 121, 233, 112,
    208, 194, 173, 168,  80, 117, 132,  72,
    202, 252, 218, 138, 214,  84,  66,  36,
    191, 152, 149, 249, 227,  94, 181,  21,
    104,  97,  40, 186, 223,  76, 241,  47,
    129, 230, 178,  63,  51, 238,  54,  16,
    110,  24,  70, 166,  34, 136,  19, 247,
     45, 184,  14,  61, 245, 164,  57,  59,
      7, 158, 156, 157,  89, 159,  27,   8,
    144,   9, 122,  28, 234, 160, 113,  90,
    209,  29, 195, 123, 174,  10, 169, 145,
     81,  91, 118, 114, 133, 161,  73, 235,
    203, 124, 253, 196, 219,  30, 139, 210,
    215, 146,  85, 170,  67,  11,  37, 175,
    192, 115, 153, 119, 150,  92, 250,  82,
    228, 236,  95,  74, 182, 162,  22, 134,
    105, 197,  98, 254,  41, 125, 187, 204,
    224, 211,  77, 140, 242,  31,  48, 220,
    130, 171, 231,  86, 179, 147,  64, 216,
     52, 176, 239,  38,  55,  12,  17,  68,
    111, 120,  25, 154,  71, 116, 167, 193,
     35,  83, 137, 251,  20,  93, 248, 151,
     46,  75, 185,  96,  15, 237,  62, 229,
    246, 135, 165,  23,  58, 163,  60, 183,
};

const uint8_t CCSDS_gen[] = {
      0, 249,  59,  66,   4,  43, 126, 251,
     97,  30,   3, 213,  50,  66, 170,   5,
     24,   5, 170,  66,  50, 213,   3,  30,
     97, 251, 126,  43,   4,  66,  59, 249,
      0,
};

static GF_FORCE_INLINE gf_sym
do_mod(unsigned int val)
{
    while (val >= 255) {
	val -= 255;
	val = (val >> 8) + (val & 255);
    }
    return val;
}
#include "rs_encode.h"
#undef RS_T
#undef GF_NP
#undef GF_EXP
#undef GF_LOG
#undef RS_LEN
#undef RS_USE_GF
#undef RS_GENERATOR
#undef RS_ENC_START
#undef RS_ENC_END

void
rs_decoder_init(struct reed_solomon_decoder *rsd,
		struct reed_solomon *rs)
{
    rsd->rs = rs;
#if GF_DYN_ALLOC
    rsd->C = malloc(sizeof(gf_sym) * (RS_MAX_T + 1));
    rsd->B = malloc(sizeof(gf_sym) * (2 * RS_MAX_T + 1));
    rsd->Temp = malloc(sizeof(gf_sym) * (RS_MAX_T + 1));

    rsd->synd = malloc(sizeof(gf_sym) * RS_MAX_T);
    rsd->error_pos = malloc(sizeof(unsigned int) * RS_MAX_ERR);
    rsd->error_idx = malloc(sizeof(unsigned int) * RS_MAX_ERR);

    rsd->O = malloc(sizeof(gf_sym) * RS_MAX_ERR);
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
static GF_FORCE_INLINE unsigned int
compute_syndromes(struct reed_solomon *rs, unsigned int pad,
		  const gf_sym *e, gf_sym *S)
{
    struct galois_field *gf = &rs->gf;
    unsigned int i, j, count = 0;

    if (pad > 0) {
	memset(S, 0, rs->T * sizeof(gf_sym));
	for (i = 1; i < pad; i++) {
	    for (j = 0; j < rs->T; j++) {
		if (S[j] == 0) {
		    S[j] = 0;
		} else {
		    /* S[j] = (S[j] * (rs->fcr + j) ^ rs->prim) */
		    S[j] = gf->exp[(gf->log[S[j]]
				    + (rs->fcr + j) * rs->prim) % gf->Np];
		    /*
		     * Direct calculation above is much faster than the below
		     * because it doesn't do all the error checks.
		     *
		     * S[j] = gf_mul_el(gf, S[j],
		     *		    gf_pow_l_l(gf, rs->fcr + j, rs->prim)));
		     */
		}
	    }
	}
    } else {
	for (i = 0; i < rs->T; i++)
	    S[i] = e[0];
	i = 1;
    }

    for (; i < gf->Np; i++) {
	for (j = 0; j < rs->T; j++) {
	    if (S[j] == 0) {
		S[j] = e[i - pad];
	    } else {
		/* S[j] = e[i] + (S[j] * (rs->fcr + j) ^ rs->prim) */
		S[j] = e[i - pad] ^ gf->exp[(gf->log[S[j]]
				       + (rs->fcr + j) * rs->prim) % gf->Np];
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
static GF_FORCE_INLINE unsigned int
berlekamp_massey(struct reed_solomon *rs,
		 gf_sym *S,
		 gf_sym *B, gf_sym *C,
		 gf_sym *Temp)
{
    struct galois_field *gf = &rs->gf;
    unsigned int L = 0;
    unsigned int i, n;

    for (i = 0; i < rs->T; i++)
	/*
	 * Put S into log format so it doesn't have to be recomputed
	 * when used.  It doesn't matter if S[i] is zero upon entry.
	 */
	S[i] = gf_log_nc(gf, S[i]);

    /*
     * Instead of copying the B array, move it backwards.  This can
     * happen at most RS_MAX_T times, so start out there.
     */
    B += RS_MAX_T;

    memset(C + 1, 0, sizeof(gf_sym) * rs->T);
    C[0] = 1;

    /* B array is in log format. */
    for (i = 0; i <= rs->T; i++)
	B[i] = gf->Np;
    B[0] = gf_log_nc(gf, 1);

    for (n = 1; n <= rs->T; n++) {
	gf_sym d = 0;

	for (i = 0; i < n; i++)
	    /* d += C[i] * S[n - i - 1] */
	    d = gf_add(d, gf_mul_el(gf, C[i], S[n - i - 1]));

	if (d == 0) {
	    B--;
	    B[0] = gf->Np;
	} else {
	    d = gf_log_nc(gf, d);
	    Temp[0] = C[0];
	    for (i = 0; i < rs->T; i++)
		Temp[i + 1] = gf_add(C[i + 1], gf_mul_ll(gf, d, B[i]));

	    if (2 * L <= n - 1) {
		for (i = 0; i <= rs->T; i++)
		    B[i] = gf_div_el_l(gf, C[i], d);

		L = n - L;
	    } else {
		B--;
		B[0] = gf->Np;
	    }

	    memcpy(C, Temp, rs->T * sizeof(gf_sym));
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
chien_search(struct reed_solomon *rs, unsigned int L,
	     gf_sym *C, gf_sym *Temp,
	     unsigned int *err_idx, unsigned int *err_pos,
	     unsigned int *rcount)
{
    struct galois_field *gf = &rs->gf;
    unsigned int count = 0, d = 0;
    unsigned int i, j, k;

    for (i = 0; i <= rs->T; i++) {
	if (C[i] != 0)
	    d = i;
	/* Convert C[i] to log format for faster processing. */
	C[i] = gf_log_nc(gf, C[i]);
	Temp[i] = C[i];
    }

    k = rs->iprim - 1;
    for (i = 1; i <= gf->Np; i++) {
	unsigned int sum = 1;

	for (j = d; j > 0; j--) {
	    Temp[j] = gf_mul_ll_l(gf, Temp[j], j);
	    sum = gf_add(sum, gf_exp(gf, Temp[j]));
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
static GF_FORCE_INLINE void
correct_errors(struct reed_solomon *rs,
	       gf_sym *e, unsigned int pad, const gf_sym *S,
	       const gf_sym *C, gf_sym *O,
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
	    O[i] = gf_add(O[i], gf_mul_ll(gf, S[i - j], C[j]));
	/* Convert O[i] to log format for faster processing. */
	O[i] = gf_log(gf, O[i]);
    }

    for (i = 0; i < err_cnt; i++) {
	unsigned int tmp = 0, tmp2, d;

	for (j = 0; j <= o; j++)
	    /* tmp += O[j] * j ^ err_idx[i] */
	    tmp = gf_add(tmp, gf_mul_ll(gf, O[j],
					gf_pow_l_l(gf, j, err_idx[i])));
	/* tmp2 = err_index[i] ^ (fcr - 1) */
	tmp2 = gf_pow_l(gf, err_idx[i], rs->fcr - 1);

	d = 0;
	if (err_cnt < rs->T - 1)
	    j = err_cnt;
	else
	    j = rs->T - 1;
	for (j = j & ~1; ; j -= 2) {
	    /* d += C[j + 1] * j ^ err_idx[i] */
	    d = gf_add(d, gf_mul_ll(gf, C[j + 1],
				    gf_pow_l_l(gf, j, err_idx[i])));
	    if (j < 2)
		break;
	}

	if (err_pos[i] >= pad)
	    e[err_pos[i] - pad] = gf_add(e[err_pos[i] - pad],
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
	  uint8_t *data, unsigned int len,
	  unsigned int *err_count)
{
    struct reed_solomon *rs = rsd->rs;
    struct galois_field *gf = &rs->gf;
    unsigned int count = 0;

    if (len > gf->Np)
	return 1;
    if (len <= rs->T)
	return 1;
    rsd->N = len;
    rsd->K = len - rs->T;
    /* Number of shortened symbols */
    rsd->S = gf->Np - rsd->N;

    /* Syndromes */
    count = compute_syndromes(rs, rsd->S, data, rsd->synd);

    if (count == 0) {
	*err_count = 0;
	return 0;
    }

    /* BM → locator polynomial */
    unsigned int L = berlekamp_massey(rs, rsd->synd, rsd->B, rsd->C,
				      rsd->Temp);

    /* Chien search */
    if (chien_search(rs, L, rsd->C, rsd->Temp,
		     rsd->error_idx, rsd->error_pos, &count))
	return 1; /* Could not correct all symbols. */

    /* Correct */
    if (count > 0)
	correct_errors(rs, data, rsd->S, rsd->synd, rsd->C, rsd->O,
		       rsd->error_idx, rsd->error_pos, count);

    *err_count = count;

    return 0;
}

/*
 * The decoder doesn't seem to benefit much from the optimizations
 * done for the encoder, so just call it for now.
 */
int
rs_decode_8(uint8_t *data, unsigned int len,
	    unsigned int *err_count)
{
    struct reed_solomon rs;
    struct reed_solomon_decoder rsd;

    reed_solomon_init(&rs, 8, 0x187, 32, 112, 11);
    rs_decoder_init(&rsd, &rs);

    return rs_decode(&rsd, data, len, err_count);
}


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

/* Pre-computed parameters for CCSDS version. */
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
CCSDS_do_mod(unsigned int val)
{
    while (val >= 255) {
	val -= 255;
	val = (val >> 8) + (val & 255);
    }
    return val;
}

#define GF_SWAP(type, x, y) do { type tmp = (x); (x) = (y); (y) = tmp; } while (0)

/* Reverse the order of the array. */
static GF_FORCE_INLINE void gf_reverse(gf_sym *x, int start, int end)
{
    while (start < end) {
	GF_SWAP(gf_sym, x[start], x[end]);
        start++;
        end--;
    }
}

/* Rotate the array x that is n bytes long left by k elements. */
static GF_FORCE_INLINE void gf_shift(gf_sym *x, unsigned int n, unsigned int k)
{
    gf_reverse(x, 0, n - 1);
    gf_reverse(x, 0, n - k - 1);
    gf_reverse(x, n - k, n - 1);
}

/* Define the general-purpose Reed-Solomon coder. */
#define RS_T rs->T
#define GF_NP rs->gf.Np
#define GF_EXP(x) gf_exp((&rs->gf), x)
#define GF_LOG(x) gf_log_nc((&rs->gf), x)
#define GF_MUL(x, y) gf_mul((&rs->gf), x, y)
#define GF_MUL_LL(x, y) gf_mul_ll((&rs->gf), x, y)
#define GF_MUL_EL(x, y) gf_mul_el((&rs->gf), x, y)
#define GF_MUL_LL_L(x, y) gf_mul_ll_l((&rs->gf), x, y)
#define GF_DIV(x, y) gf_div((&rs->gf), x, y)
#define GF_DIV_EL_L(x, y) gf_div_el_l((&rs->gf), x, y)
#define GF_POW_L(x, y) gf_pow_l((&rs->gf), x, y)
#define GF_POW_L_L(x, y) gf_pow_l_l((&rs->gf), x, y)
#define RS_LEN len
#define RS_USE_GF 1
#define RS_GENERATOR rse->generator
#define RS_ENC_START() \
int \
rs_encode(struct reed_solomon_encoder *rse, \
	  uint8_t *inbuf, unsigned int len, uint8_t *parity) \
{ \
    struct reed_solomon *rs = rse->rs;
#define RS_ENC_END \
}
#define RS_DEC_START() \
int \
rs_decode(struct reed_solomon_decoder *rsd, \
	  uint8_t *data, unsigned int len, \
	  unsigned int *err_count) \
{ \
    struct reed_solomon *rs = rsd->rs;
#define RS_DEC_END \
}
#define RS_NAME(x) x

#include "rs_encode.h"
#include "rs_decode.h"

#undef RS_T
#undef GF_NP
#undef GF_EXP
#undef GF_LOG
#undef GF_MUL
#undef GF_MUL_LL
#undef GF_MUL_EL
#undef GF_MUL_LE
#undef GF_MUL_LL_L
#undef GF_DIV
#undef GF_DIV_EL_L
#undef GF_POW_L_L
#undef GF_POW_L
#undef RS_LEN
#undef RS_USE_GF
#undef RS_GENERATOR
#undef RS_ENC_START
#undef RS_ENC_END
#undef RS_DEC_START
#undef RS_DEC_END
#undef RS_NAME

/* Define the high-performance CCSDS Reed-Solomon coder. */
#define RS_T 32
#define GF_NP 255
#define do_mod(x) CCSDS_do_mod(x)
#define GF_EXP(x) (((x) == GF_NP) ? 0 : CCSDS_exp[x])
#define GF_LOG(x) CCSDS_log[x]
#define GF_MUL(x, y) (((x) == 0 || (y) == 0) ? 0 : GF_EXP(do_mod(GF_LOG(x) + GF_LOG(y))))
#define GF_MUL_LL(x, y) (((x) == GF_NP || (y) == GF_NP) ? 0 : GF_EXP(do_mod(x + y)))
#define GF_MUL_EL(x, y) (((x) == 0 || (y) == GF_NP) ? 0 : GF_EXP(do_mod(GF_LOG(x) + y)))
#define GF_MUL_LL_L(x, y) (((x) == GF_NP || (y) == GF_NP) ? GF_NP : do_mod(x + y))

#define GF_DIV(x, y) (((x) == 0) ? 0 : GF_EXP(do_mod(GF_LOG(x) - GF_LOG(y) + GF_NP)))
#define GF_DIV_EL_L(x, y) (((x) == 0) ? GF_NP : do_mod(GF_LOG(x) - (y) + GF_NP))
#define GF_POW_L(x, y) ((x) == GF_NP ? 0 : GF_EXP(do_mod(((x) * (y)) % GF_NP + GF_NP)))
#define GF_POW_L_L(x, y) ((x) == GF_NP ? GF_NP : do_mod(((x) * (y)) % GF_NP + GF_NP))
#define RS_GENERATOR CCSDS_gen
#define RS_LEN len
#define RS_USE_GF 0
#define RS_ENC_START() \
int \
rs_encode_8(uint8_t *inbuf, unsigned int len, uint8_t *parity)	\
{
#define RS_ENC_END \
}
#define RS_DEC_START() \
int \
rs_decode_8(uint8_t *data, unsigned int len, \
	    unsigned int *err_count)	     \
{ \
    struct reed_solomon rs_data; \
    struct reed_solomon *rs = &rs_data; \
    struct reed_solomon_decoder rsd_data; \
    struct reed_solomon_decoder *rsd = &rsd_data; \
    reed_solomon_init(rs, 8, 0x187, 32, 112, 11); \
    rs_decoder_init(rsd, rs);
#define RS_DEC_END \
}
#define RS_NAME(x) CCSDS_ ## x

#if DO_SIMD && 1
typedef int16_t gf_v16ss __attribute__ ((vector_size (16)));
static gf_v16ss CCSDS_simd_generator[4] = {
    {   0, 249,  59,  66,   4,  43, 126, 251 },
    {  97,  30,   3, 213,  50,  66, 170,   5 },
    {  24,   5, 170,  66,  50, 213,   3,  30 },
    {  97, 251, 126,  43,   4,  66,  59, 249 }
};
#define SIMD_GEN CCSDS_simd_generator
#include "rs_encode_simd.h"
#undef SIMD_GEN
#else
#include "rs_encode.h"
#endif
#include "rs_decode.h"

#undef RS_T
#undef GF_NP
#undef GF_EXP
#undef GF_LOG
#undef GF_MUL
#undef GF_MUL_LL
#undef GF_MUL_EL
#undef GF_MUL_LE
#undef GF_MUL_LL_L
#undef GF_DIV
#undef GF_DIV_EL_L
#undef GF_POW_L
#undef GF_POW_L_L
#undef RS_LEN
#undef RS_USE_GF
#undef RS_GENERATOR
#undef RS_ENC_START
#undef RS_ENC_END

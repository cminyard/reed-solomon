
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <time.h>

#define DO_RS_CHECK 1
#define DO_LIBFEC_CHECK 1

#include "reed_solomon.h"

#if DO_LIBFEC_CHECK
#include <fec.h>
#endif

static unsigned int
test_one(unsigned int num_errs,
	 unsigned int enc_len,
	 struct reed_solomon_encoder *rse,
	 struct reed_solomon_decoder *rsd)
{
    unsigned int i;
    uint8_t origbuf[255];
#if DO_RS_CHECK
    uint8_t buf[255];
#endif
#if DO_LIBFEC_CHECK
    uint8_t buf2[255];
#endif
    bool errpos[255] = { false };
    unsigned int errcount = 0;

    if (enc_len == 0)
	/* If not specified, pick a random length between 1 and 223. */
	enc_len = rand() % 223 + 1;

    for (i = 0; i < enc_len; i++) {
	origbuf[i] = rand();
#if DO_RS_CHECK
	buf[i] = origbuf[i];
#endif
#if DO_LIBFEC_CHECK
	buf2[i] = origbuf[i];
#endif
    }

#if DO_RS_CHECK
    /* Add the parity bytes. */
    rs_encode(rse, buf, enc_len, buf + enc_len);
#endif

#if DO_LIBFEC_CHECK
    /* Do it with libfec to compare. */
    encode_rs_8(buf2, buf2 + enc_len, 223 - enc_len);

#if DO_RS_CHECK
    for (i = 0; i < enc_len + 32; i++) {
	if (buf[i] != buf2[i]) {
	    printf("Diff1(%d): %2.2x %2.2x\n", i, buf[i], buf2[i]);
	    return 1;
	}
    }
#endif
#endif

#if 0
    printf("Encoded:");
    for (i = 0; i < 255; i++) {
	if (i % 16 == 0)
	    printf("\n%2.2x:", i);
	printf(" %2.2x", buf[i]);
    }
    printf("\n");
#endif

    /* Inject some errors. */
    for (i = 0; i < num_errs; ) {
	unsigned int pos = rand() % (enc_len + 32) * 8;

	if (errpos[pos / 8])
	    continue;

#if DO_RS_CHECK
	buf[pos / 8] ^= 1 << (pos % 8);
#endif
#if DO_LIBFEC_CHECK
	buf2[pos / 8] ^= 1 << (pos % 8);
#endif
#if 0
	printf("Injecting error at byte %u\n", pos / 8);
#endif
	errpos[pos / 8] = true;
	i++;
    }
#if DO_RS_CHECK
    if (rs_decode(rsd, buf, enc_len + 32, &errcount)) {
	if (num_errs <= 16)
	    printf("Decode error\n");
	return 1;
    } else {
	if (errcount != num_errs) {
	    printf("Error count mismatch: %u %u\n", num_errs, errcount);
	    return 1;
	}
    }

    for (i = 0; i < enc_len; i++) {
	if (buf[i] != origbuf[i]) {
	    printf("Data mismatch\n");
	    return 1;
	}
    }
#endif

#if DO_LIBFEC_CHECK
    errcount = decode_rs_8(buf2, NULL, 0, 223 - enc_len);
    if (errcount != num_errs) {
	printf("Bad err count 1\n");
    }

    for (i = 0; i < enc_len; i++) {
	if (buf2[i] != origbuf[i]) {
	    printf("Diff2(%d): %2.2x %2.2x\n", i, origbuf[i], buf2[i]);
	    return 1;
	}
    }
#endif

    return 0;
}

static unsigned int
test_loop(unsigned int num_loops,
	  unsigned int num_errs,
	  unsigned int enc_len,
	  struct reed_solomon_encoder *rse,
	  struct reed_solomon_decoder *rsd)
{
    unsigned int i, errcount = 0;

    for (i = 0; i < num_loops; i++)
	errcount += test_one(num_errs, enc_len, rse, rsd);
    printf("Testing with %u errors: %u output errors\n", num_errs, errcount);
    return errcount;
}

int
main(int argc, char *argv[])
{
    struct reed_solomon rs;
    struct reed_solomon_encoder rse = { .rs = &rs };
    struct reed_solomon_decoder rsd = { .rs = &rs };
    unsigned int loops = 100;
    bool do_cpu_usage = false;
    unsigned int num_errs = 8;
    unsigned int enc_len = 0;
    unsigned int i;
    int arg, err = 0;

    for (arg = 1; arg < argc; arg++) {
	if (argv[arg][0] != '-')
	    break;
	if (strcmp(argv[arg], "-c") == 0) {
	    do_cpu_usage = true;
	} else if (strcmp(argv[arg], "-l") == 0) {
	    arg++;
	    if (arg >= argc) {
		fprintf(stderr, "No data supplied for -l\n");
		return 1;
	    }
	    loops = strtoul(argv[arg], NULL, 0);
	} else if (strcmp(argv[arg], "-n") == 0) {
	    arg++;
	    if (arg >= argc) {
		fprintf(stderr, "No data supplied for -n\n");
		return 1;
	    }
	    num_errs = strtoul(argv[arg], NULL, 0);
	} else if (strcmp(argv[arg], "-e") == 0) {
	    arg++;
	    if (arg >= argc) {
		fprintf(stderr, "No data supplied for -e\n");
		return 1;
	    }
	    enc_len = strtoul(argv[arg], NULL, 0);
	} else {
	    fprintf(stderr, "unknown option: %s\n", argv[arg]);
	    return 1;
	}
    }
    srand(time(NULL));

    reed_solomon_init(&rs, 8, 0x187, 32, 112, 11);
    rs_encoder_init(&rse, &rs);
    rs_decoder_init(&rsd, &rs);

    if (do_cpu_usage) {
	test_loop(loops, num_errs, enc_len, &rse, &rsd);
	return 0;
    }

    for (i = 0; i < 32; i++) {
	unsigned int errs = test_loop(loops, i, enc_len, &rse, &rsd);

	if (errs == loops && i <= 16) {
	    err = 1;
	    break;
	}
	if (errs != loops && i > 16) {
	    err = 1;
	    break;
	}
    }

    return err;
}

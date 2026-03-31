
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <time.h>
#include "reed_solomon.h"
#include <fec.h>

#define ERR_COUNT_CHECK 1

static unsigned int
test_one(unsigned int num_errs,
	 struct reed_solomon_encoder *rse,
	 struct reed_solomon_decoder *rsd)
{
    unsigned int i;
    uint8_t origbuf[255], buf[255];
    bool errpos[255] = { false };
    unsigned int err = 0;
    unsigned int errcount = 0;

    for (i = 0; i < 223; i++) {
	buf[i] = rand();
	origbuf[i] = buf[i];
    }

    /* Add the parity bytes. */
    reed_solomon_encode(rse, buf, 223, buf + 223);
#if 0
    uint8_t buf2[255];

    for (i = 0; i < 223; i++)
	buf2[i] = origbuf[i];

    encode_rs_8(buf2, buf2 + 223, 0);

    for (i = 0; i < 255; i++) {
	if (buf[i] != buf2[i])
	    printf("Diff(%d): %2.2x %2.2x\n", i, buf[i], buf2[i]);
    }
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
	unsigned int pos = rand() % 255 * 8;

	if (errpos[pos / 8])
	    continue;

	buf[pos / 8] ^= 1 << (pos % 8);
#if 0
	buf2[pos / 8] ^= 1 << (pos % 8);
#endif
#if 0
	printf("Injecting error at byte %u\n", pos / 8);
#endif
	errpos[pos / 8] = true;
	i++;
    }
    if (reed_solomon_decode(rsd, buf, 255, &errcount)) {
	return 1;
    } else {
#if ERR_COUNT_CHECK
	if (errcount != num_errs) {
	    printf("Error count mismatch: %u %u\n", num_errs, errcount);
	    return 1;
	}
#endif
    }
#if 0
    printf("errcount = %u\n", errcount);
    errcount = decode_rs_8(buf2, NULL, 0, 0);
    printf("errcount(2) = %u\n", errcount);

    for (i = 0; i < 223; i++) {
	if (buf[i] != origbuf[i])
	    printf("Diff(%d): %2.2x %2.2x\n", i, buf[i], buf2[i]);
    }
#endif

    for (i = 0; i < 223; i++) {
	if (buf[i] != origbuf[i]) {
	    printf("Data mismatch\n");
	    err = 1;
	    break;
	}
    }

    return err;
}

static unsigned int
test_loop(unsigned int num_loops,
	  unsigned int num_errs,
	  struct reed_solomon_encoder *rse,
	  struct reed_solomon_decoder *rsd)
{
    unsigned int i, errcount = 0;

    for (i = 0; i < num_loops; i++)
	errcount += test_one(num_errs, rse, rsd);
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
	} else {
	    fprintf(stderr, "unknown option: %s\n", argv[arg]);
	    return 1;
	}
    }
    srand(time(NULL));

    reed_solomon_init(&rs, 8, 0x187, 32, 112, 11);
    reed_solomon_encoder_init(&rse, &rs);
    reed_solomon_decoder_init(&rsd, &rs);

    for (i = 0; i < 128; i++) {
	unsigned int errs = test_loop(loops, i, &rse, &rsd);

	if (do_cpu_usage)
	    break;

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

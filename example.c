
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <time.h>
#include "reed_solomon.h"

#define ERR_COUNT_CHECK 0

static unsigned int
test_one(unsigned int num_errs,
	 struct reed_solomon_encoder *rse,
	 struct reed_solomon_decoder *rsd)
{
    unsigned int i;
    uint8_t origbuf[255], buf[255];
    bool errpos[255] = { false };
    unsigned int err = 0;
#if ERR_COUNT_CHECK
    unsigned int errcount = 0;
#endif

    for (i = 0; i < 223; i++) {
	buf[i] = rand();
	origbuf[i] = buf[i];
    }
    for (; i < 255; i++)
	buf[i] = 0;

    /* Add the parity bytes. */
    reed_solomon_encode(rse, buf);

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
	unsigned int pos = rand() % rse->rs->N * 8;

	if (errpos[pos / 8])
	    continue;

	buf[pos / 8] ^= 1 << (pos % 8);
	errpos[pos / 8] = true;
	i++;
    }
#if ERR_COUNT_CHECK
    errcount = reed_solomon_decode(rsd, buf);
    if (errcount != num_errs)
	printf("Error count mismatch: %u %u\n", num_errs, errcount);
#else
    reed_solomon_decode(rsd, buf);
#endif

    for (i = 0; i < 223; i++) {
	if (buf[i] != origbuf[i]) {
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
    unsigned int i;
    int err = 0;

    srand(time(NULL));

    reed_solomon_init(&rs, 8, 255, 223);

    for (i = 0; ; i++) {
	unsigned int errs = test_loop(100, i, &rse, &rsd);

	if (errs == 100)
	    break;
    }

    return err;
}

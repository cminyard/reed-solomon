
all: tester

# There is SIMD implemented for rs_encoder_8, but it doesn't seem to
# improve the performance at all, at least on the AMD processor I'm
# using.
DO_SIMD = 1

# X86 support
# popcnt was added in the Nehalem CPU
ARCH = -march=native
# If you set DO_SIMD, adding -msse2, or whatever for your processor,
# will be required.
ARCH += -msse2

CFLAGS = -g -Wall -O2 -DCONVCODE_TESTS -DDO_SIMD=$(DO_SIMD) $(ARCH)

tester: tester.o reed_solomon.o galois_field.o
	$(CC) $(CFLAGS) -o $@ $^ -lfec

tester.o: tester.c reed_solomon.h galois_field.h

reed_solomon.o: reed_solomon.c reed_solomon.h galois_field.h \
	rs_encode.h rs_decode.h rs_encode_simd.h

galois_field.o: galois_field.c galois_field.h

check:	tester
	./tester

clean:
	rm -f *.o tester

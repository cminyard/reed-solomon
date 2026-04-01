
all: tester

# No SIMD support yet, just a placeholder.
DO_SIMD = 0

# X86 support
# popcnt was added in the Nehalem CPU
ARCH = -march=native
# If you set DO_SIMD, adding -msse2 will be required.
#ARCH += -msse2

CFLAGS = -g -Wall -O2 -DCONVCODE_TESTS -DDO_SIMD=$(DO_SIMD) $(ARCH)

tester: tester.o reed_solomon.o galois_field.o
	gcc $(CFLAGS) -o $@ $^ -lfec

tester.o: tester.c reed_solomon.h galois_field.h

reed_solomon.o: reed_solomon.c reed_solomon.h galois_field.h \
	rs_encode.h rs_decode.h

galois_field.o: galois_field.c galois_field.h

check:	tester
	./tester

clean:
	rm -f *.o tester

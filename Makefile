
all: example

# X86 support
# popcnt was added in the Nehalem CPU
ARCH = -march=native
# If you set DO_SIMD, adding -msse2 will be required.
#ARCH += -msse2

CFLAGS = -g -Wall -DCONVCODE_TESTS -DDO_SIMD=$(DO_SIMD) $(ARCH)

example: example.o reed_solomon.o galois_field.o
	gcc $(CFLAGS) -o $@ $^

example.o: example.c reed_solomon.h galois_field.h

reed_solomon.o: reed_solomon.c reed_solomon.h galois_field.h

galois_field.o: galois_field.c galois_field.h

clean:
	rm -f *.o example

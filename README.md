Reed-Solomon Coder
==================

This is a production-quality RS coder implementation as a couple of C
files.  It's taken from
https://github.com/fujiyama-kota-comm/fec-rs-codec which I found easy
to understand and well written, and perhaps good for playing with, but
not good for a production environment.  And, it turns out, it was
wrong and I had to rewrite all the algorithms pretty much from
scratch.  Not much is left of the original except the overall
structure.

See reed_solomon.h and tester.c for how to use.

The Galois field implementation is in galois\_field.c and the API is in
galois\_field.h

There are higher-performance versions of the functions specifically
for CCSDS coding (without dual basis at the moment).  They don't
actually have that much better performance, about 10% better.  But
they are a little more convenient if that's what you are doing.

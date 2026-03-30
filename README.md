Reed-Solomon Coder
==================

This is a production-quality RS coder implementation as a couple of C
files.  It's taken from
https://github.com/fujiyama-kota-comm/fec-rs-codec which I found easy
to understand and well written, and perhaps good for playing with, but
not good for a production environment.

See reed_solomon.h and example.c for how to use.

The Galois field implementation is in galois\_field.c and the API is in
galois\_field.h

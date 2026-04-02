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

## CCSDS version

There are higher-performance versions of the functions specifically
for CCSDS coding (without dual basis at the moment).  They don't
actually have that much better performance, about 10% better in my
testing.  But they are a little more convenient if that's what you are
doing.

## SIMD

There is SIMD code in here using the gcc vector operations.  It should
be portable to systems with 16-byte vector operations.  You can enable
it by setting it in the makefile and setting the ARCH per the comments.

In my tests it doesn't improve performance at all on encode and slows
it down a lot on decode.

It looks like the compiler produces decent code for the vector
operations, for the most part, at least in my looking at it.

It appears the primary issue is having to do the log() and exp()
operations, since that data has to be pulled out a symbol at a time
and processed, then put back into the vector.  This is done once in
the encoder and twice in the decoder.  If those could be done with
vector operations then SIMD would be a lot faster, but I haven't
figured out a way to do it.

Also, on x86, it appears getting data into and out of the mmx
registers requires going through memory.  There doesn't appear to be a
way to say "set each element to this value" or "pull this element out
and put it in a processor register.

I'm not an expert in these matters, so I'm leaving it here and maybe
someone can improve it.

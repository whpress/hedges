# HEDGES

A package for encoding and decoding arbitrary byte data to and
from strands of DNA using a robust an error-correcting code (ECC)
that accounts for substitutions, insertions, and deletions.

### HEDGES Error-Correcting Code for DNA Storage Corrects Indels and Allows Sequence Constraints

**William H. Press, John A. Hawkins, Stephen Knox Jones Jr, Jeffrey M. Schaub, and Ilya J. Finkelstein**

*Proc Natl Acad Sci*. 117 (31) 18489-18496 (July 16, 2020)

### History

HEDGES was originally implemented as C++ modules callable from Python
via the PyObject* interface. However, over years, this proved fragile, because
Python updates, and updates of the C++ runtime libraries,
tended to break the interface.

This repository gives a new, pure C++ implementation
that ought to be less susceptible to bitrot. The identical code compiles
and executes (as of November, 2024) in both Linux
and Windows with no differences necessary in the source files.

Here supplied are the new source files, a Linux script `build_linux_bash.sh`
that compiles them under Linux, a Windows batch file
`build_windows_command_line.bat` that compiles them in the Windows C++ CL
command-line environment, and a directory `Windows_VisualStudio_IDE` that compiles
them as a Visual Studio project (Visual Studio Community 2022).

### Numerical Recipes conventions

The source files make use of a version of the Numerical Recipes class
library, here supplied as `nr3b.h`. Principally, this typedefs the integer
and float types with unique names for cross-platform compatibility, and defines
some classes for vectors and matrices.

For a casual reading of the source code, all you really need to know
is
```
typedef int32_t Int; // 32 bit integer
typedef uint32_t Uint; // unsigned 32 bit integer
typedef uint64_t Ullong; // unsigned 64 bit integer
typedef char Char; // 8 bit integer
typedef unsigned char Uchar; // unsigned byte
typedef double Doub; // default double
typedef bool Bool; // boolean
```
and that vector and matrix declarations look like this
```
VecUchar foo(1000); // declare vector of bytes, length 1000
MatUchar bah(100,400); // declare matrix of bytes, 100 rows x 400 cols
```
The vector and matrix classes have a bunch of overloaded constructors,
which you can find in `nr3b.h`, but which should be obvious for
casual code-reading. Similarly, there are some obvious vector and
matrix manipulations, like assignment, row-of-matrix-as-a-vector, etc.

### Usage

What is supplied are functions for mix-and-match use, not a single
software package.  The program `DNAcode_demo.cpp` illustrates the
use of the most important ones, and here is a quick guide to the
most important.

Prerequisite: Read the PNAS paper to understand how any very long
string of bytes gets divided up into packets, each packet containing
strands. The strands contain their packet number and sequence number,
so the idea is that bytes are encoded into an unordered
pool of DNA strands, which can then be sequenced, ordered into
packets, decoded (with error correction, including on the packet
and sequence number), and then re-assembled to the original string
of bytes.

```class GetWizBytes;```
is a helper program that streams bytes from a
file for tests, defaulting to the text of the Wizard of Oz
in Esperanto if available, or else random bytes.

```Createmesspacket_out createmesspacket(Int packno);```\
streams the next bytes from `GetWizBytes` into a numbered packet, returning
the packet and the bytes used in a structure,\
```struct Createmesspacket_out { MatUchar packet; VecUchar plaintext; };```\

`MatUchar protectmesspacket(const MatUchar& packetin);`\
fills the outer-code Reed-Solomon check bytes into the packet,
returning a protected packet.

`MatUchar messtodna(const MatUchar& mpacket);`\
performs the inner-code encoding of the packet, returning
a matrix whose rows are DNA strands, that is, vectors of bytes whose numerical values
(0, 1, 2, 3) to be interpreted as (A, C, G, T).

For testing purposes,\
```MatUchar createerrors(const MatUchar& dnabag, Doub srate, Doub drate, Doub irate);```
takes a matrix whose rows are DNA strands and creates errors with
specified rates for substitutions, deletions, and insertions.

That's the encoding. Now the decoding:

```Dnatomess_out dnatomess(const MatUchar& dnapacket);```\
decodes the inner code on the strands of an assumed packet.
This is for testing. In real life, one would decode the strands of
an arbitrary "bag" of strands, then use their decoded packet
and serial numbers to arrange them into packets before applying
the next step.

The real meat of HEDGES, called within `dnatomess()`, is the decoding function\
```Decode_out decode(VecUchar &codetext)```
which takes as input a vector of DNA (values 0..3) and decodes it, with error
correction, to a vector of byte plaintext. The structure returned contains
also diagnostic information:\
```struct Decode_out {Int errcode; VecUchar plaintext; Int nhypo; Doub finalscore; Int finaloffset; Int finalseq;};```


The struct returned by the higher-level `dnatomess` has the whole
error-corrected packet, a matrix marking where erasures are
detected, and some statistics:\
```struct Dnatomess_out { MatUchar mpacket; MatUchar epacket; Int baddecodes; Int erasures; };```

Then,\
```Correctmesspacket_out correctmesspacket(const MatUchar& packetin, const MatUchar& epacket);```\
decodes the outer code Reed-Solomon code (using the erasure information), returning a\
```struct Correctmesspacket_out { MatUchar packet; Int tot_detect; Int tot_uncorrect; Int max_detect; Int max_uncorrect; Int toterrcodes; };```\
with a corrected packet and some statistics.

Finally,\
```VecUchar extractplaintext(const MatUchar& cpacket);```\
returns the string of plaintext contained in the packet.

There are also a bunch of lower level functions for getting
and setting parameters, helper functions, etc. These have one-line
descriptions in the header file `DNAcode.h`. Best read the
code to understand these in detail.

The program `main_test_all()` shows how the various parameters
are set, how the above functions are invoked, and whether
(for a given error rate) perfect packet decodes are obtained.
See the PNAS paper for statistics applicable to very large
plaintext corpora.

### What You Can Do Next?

Try compiling and running the program `DNAcode_demo.cpp` with the parameters
as now set. It should print diagnostic information for 20 packets.
The output will show a large number of errors detected and corrected,
but there should show no errors after the inner and outer decodes are
applied.

Because different random errors are intentionally introduced each time
the test is run, the results will not be exactly identical from run to
run. With the default program settings for errors introduced and
coderate, there can be, rarely, a bad decode. See PNAS paper for
what coderates and error rates can give exponentially small decode errors.

Then, for your own applications, study how `DNAcode_demo.cpp` makes use
of the lower-level functions provided.

Good luck!

### LICENSE (MIT License)

Copyright 2020 by William H. Press, John A. Hawkins, Stephen K. Jones jr, and
Ilya J. Finkelstein

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights to use,
copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software,
and to permit persons to whom the Software is furnished to do so, subject to the
following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,
INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE
OR OTHER DEALINGS IN THE SOFTWARE.


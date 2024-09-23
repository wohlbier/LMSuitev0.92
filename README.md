Tested on Ubuntu 20.04 x86 and aarch64 with gfortran 9.3.0.

Build dependencies
- Fortran compiler
- cmake

The sections in the [UsersManual](./Documentation/UsersManual.pdf) that
describe the building and the running of the code and paths to output are
no longer accurate. `CMake` and out of tree builds were not a thing in 2003.
Use the steps below to build and run, and use your thinking cap to figure out
discrepancies between your experience and the manual. I will entertain pull
requests as warranted.

```
mkdir build
cd build
../arch/gfortran.sh
make
./bin/lmsuite ../inputs/
ls ./outputs
```

===============================================================================

The original author of lmsuite is

John G. Wohlbier

Carnegie Mellon University

Software Engineering Institute

jgwohlbier sei cmu edu

lmsuite was written while doing graduate work at the University of Wisconsin
under Prof. John Booske.

The author has not worked on lmsuite or on traveling wave tubes in over
~~a~~ two decades.

As of 2024-09-23 the code builds on both aarch64 and x86 using
cmake 3.30.3 and gfortran 14.2.0. Older versions of the tools probably will
also work going pretty far back.

As of 2014-11-14 the code was compiled with Intel fortran 14.0.4 on
a linux system and run on the default problem. Quick inspection of the outputs
indicate that the code is doing the right thing, but there is no guarantee that
results are correct or that problems will not be discovered. Any additional
support for lmsuite from the author, including support for other compilers
and/or platforms, would require funding.

If you use lmsuite in published research, please cite the following article:
J.G. Wohlbier, J.H. Booske, and I. Dobson.
The Multifrequency Spectral Eulerian (MUSE) Model of a Traveling Wave
Tube. IEEE Trans. Plasma Sci. Vol. 30, no. 3, June 2002.

jgw
2014-11-14

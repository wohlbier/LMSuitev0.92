Build dependencies
- Fortran compiler
- cmake

```
mkdir build
cd build
../arch/gfortran.sh
./bin/lmsuite
```

The original author of lmsuite is
John G. Wohlbier
Carnegie Mellon University
Software Engineering Institute
jgwohlbier@sei.cmu.edu

lmsuite was written while doing graduate work at the University of Wisconsin
under Prof. John Booske.

The author has not worked on lmsuite or on traveling wave tubes in over a
decade. As of 2014-11-14 the code was compiled with Intel fortran 14.0.4 on
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

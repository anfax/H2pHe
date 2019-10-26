# H2pHe
Potential energy surface for H<sub>2</sub><sup>+</sup>--He

**Requirements**

(1) RKHS toolkit : Download from https://github.com/MeuwlyGroup/RKHS

**Compile a particular PES**

A test program file (pes_test.f90) is given and it can be compiled for the FCI PES as

`gfortran RKHS.f90 pes_test.f90 He--H2+_FCI.f90`

**Running the executable**

Before running the executable make sure that the pes1.csv and or pes1.kernel file for that PES present at the same directory.

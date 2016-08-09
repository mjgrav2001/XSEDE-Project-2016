README for Cntor

Use:
Command line flags:
  -h [ --help ]             Produce help message
  -b [ --begin-energy ] arg Energy at which to begin calculating occupancies
  -s [ --energy-step ] arg  Step size through energy values
  -e [ --end-energy ] arg   Ending energy
  -n [ --n-atoms ] arg      Number of atoms in system
  -g [ --group-size ] arg   Size of MPI process group inverting 1 matrix
  -t [ --n-threads ] arg    Number of threads for each process
  -d [ --dense ]            Solve with MKL dense solver rather than PETSc
  -z [ --hermitian ]        Solve with methods specifically for hermitian matrices

Programming changes:
    New functions:
    In general.h: applyCmdLine() takes Parameters object of defaults and command line flags, returns new parameters object
    In general.h: ReadIVvalues() reads the DETEfEoutput.txt file and fills an array with the modulus of each T(E) value



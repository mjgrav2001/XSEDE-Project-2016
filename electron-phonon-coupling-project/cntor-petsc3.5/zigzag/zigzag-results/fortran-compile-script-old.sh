gfortran cntor-current-b0-0d32-blas-10-0-3600.f -o a.out -Wall -O3 -openmp -I$(MKLROOT)/include $(MKLROOT)/lib/intel64/libmkl_scalapack_lp64.a -Wl,--start-group $(MKLROOT)/lib/intel64/libmkl_intel_lp64.a $(MKLROOT)/lib/intel64/libmkl_core.a $(MKLROOT)/lib/intel64/libmkl_intel_thread.a -Wl,--end-group $(MKLROOT)/lib/intel64/libmkl_blacs_intelmpi_lp64.a -lpthread -lm


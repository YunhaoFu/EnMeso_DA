FC=mpiifort -O0 -g -DDEBUG -traceback -check bounds #-check pointers -check uninit -debug all -warn all
FCFLAGS=-cpp -I$(NETCDF_DIR)/include
#LDFLAGS=-L$(NETCDF_DIR)/lib -lnetcdff -lnetcdf -lhdf5_hl -lhdf5 -lz -L$(MKLROOT)/lib/intel64 -lmkl_rt -lpthread -ldl
LDFLAGS=-L$(NETCDF_DIR)/lib -lnetcdff -lnetcdf -L$(HDF5_DIR)/lib -lhdf5_hl -lhdf5 -lz -L$(MKLROOT)/lib/intel64 -lmkl_rt -lpthread -ldl

EnKF.exe : letkf_core.o kdtree.o precision.o enkf_obs.o enkf_mpi.o EnKF_MPI_IF.o EnKF_DM.o EnKF_IO.o letkf_driver.o EnKF_Meso.o

	$(FC) $^ -o $@ $(LDFLAGS)

%.o : %.f90
	$(FC) $^ -c $(FCFLAGS)

clean:
	rm -f *.o *.mod *.exe


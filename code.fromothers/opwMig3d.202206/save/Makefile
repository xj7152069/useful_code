include $(ICLUSTERHOME)/sys_configure

FFTW_INC_DIR=-I/usr/local/include
FFTW_LIB_DIR=-L/usr/local/lib

MPI_INC_DIR=-I$(MPI_PATH)/include
FB_IO_DIR=.

INC_DIR=-I$(SRC_ROOT)/include 
BIN_DIR=$(INSTALL_DIR)/bin
LIB_DIR=$(INSTALL_DIR)/lib

MPI_LIB= -L$(MPI_LIB_DIR) -lfmpich -lmpich
PRIV_LIB=-L$(LIB_DIR) -lmigT
#MKL_PATH=/mnt/data5/apps/intel
#2015.01.27, FengBo modified MKL path.
MKL_PATH=/mnt/data5/apps/intel/mkl

ifeq ($(IA) , ia32)
	MKL_LIB = -L$(MKL_PATH)/lib/32 -lmkl_ia32 -lguide -lpthread
endif

ifeq ($(IA) ,em64t)
	MKL_LIB = -L$(MKL_PATH)/lib/em64t -lmkl_em64t -lguide -lpthread
endif

ifeq ($(IA) , ia64)
	MKL_LIB = -L$(MKL_PATH)/lib/64 -lmkl_ipf -lguide -lpthread
endif

F90=$(MPI_BIN_DIR)/mpif90
FC=$(MPI_BIN_DIR)/mpif90
MC=$(MPI_BIN_DIR)/mpicc
CC=icc
FFLAGS=-O3

###########################################################################################################
#		########   For 3d Processing. ########                                                    #
#		########   Author : Feng Bo.  ########                                                    #
#		########   Date   : 2009.08.  ########                                                    #
#		########   Version: v1.       ########                                                    #
#	step 1. 3d PlaneWave Decomposition.                                                               #
#		! cmpdata(ns,ntr,nmx,nmy)  => pwdat3d(ns,nmx,nmy,nphr)                                    #
#	step 2. 3d PlaneWave Migration.                                                                   #
#		! pwdat3d(ns,nmx,nmy,nphr) => pwmig3d(ns,nmx,nmy,nphr)                                    #
###########################################################################################################

#	step 1. 3d PlaneWave Decomposition.
#pwd3dsrc=./3d_Phx_planewave_decomposition_for_iCluster_v16.omp.c
pwd3dsrc=./3d_Phx_planewave_decomposition_for_iCluster_v16.omp.forTianjin.c
pwd3dexec=./pwd3domp.exe
pw3d:
	$(MC) -o $(pwd3dexec) $(pwd3dsrc) -std=c99 -openmp $(INC_DIR) $(PRIV_LIB) $(MPI_LIB) $(MKL_LIB) -lm

#	step 2. 3d PlaneWave Migration.
pwmig3dsrc=./3d_Phx_planewave_acoustic_pstm_for_iCluster_v19f.f
pwmig3dexec=./pwmig3d.exe
fbio.o: gdbfc_fbrw_io.c
	$(CC) -c  gdbfc_fbrw_io.c -o fbio.o
mig3d: fbio.o
	$(F90) $(FFLAGS) -openmp -o $(pwmig3dexec) $(pwmig3dsrc) $(FB_IO_DIR)/fbio.o $(PRIV_LIB) $(MPI_LIB) $(MKL_LIB)

clean:
	-rm  -f *.o $(pwd3dexec) $(pwmig3dexec)

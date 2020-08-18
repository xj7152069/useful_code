
  FC      = mpifort
  FFLAGS  = -g -O

  CC      = mpicc
  CFLAGS  = -g -O
  

  files =	lsqrDataModule.o	lsqrblas.o			lsqrblasInterface.o \
  			lsqrModule.o		lsqrCheckModule.o	lsqrTestModule.o   	\
			smooth2.o       	\
			tomo.o

all: lsqrDataModule		lsqrblas 		lsqrblasInterface lsqrModule 	\
	 lsqrCheckModule	lsqrTestModule 	smooth2 \
	 tomo tomography_modify

lsqrDataModule:
	$(FC) ${FFLAGS} -c lsqrDataModule.f90

lsqrblas:
	$(FC) ${FFLAGS} -c lsqrblas.f90

lsqrblasInterface:
	$(FC) ${FFLAGS} -c lsqrblasInterface.f90

lsqrModule:
	$(FC) ${FFLAGS} -c lsqrModule.f90

lsqrCheckModule:
	$(FC) ${FFLAGS} -c lsqrCheckModule.f90

lsqrTestModule:
	$(FC) ${FFLAGS} -c lsqrTestModule.f90

smooth2:
	$(CC) $(CFLAGS) -c smooth2.c

tomo:
	${FC} ${FFLAGS} -c tomography_lsqr_0313_modify.f90 -o tomo.o

tomography_modify: ${files}
	${FC} ${FFLAGS} ${files} -o Tomography_modify


clean:
	\rm -f *.o *.mod

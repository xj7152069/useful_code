  FC      = mpifort
#  FFLAGS  =-O3
  FFLAGS  = -g -O0
  #-check bounds

  CC      = mpicc
#  CFLAGS  = -O3
  CFLAGS  = -g -O0
  #-check-pointers=rw
  
  LIB     = 

  files =	lsqrDataModule.o	lsqrblas.o	global.o	lsqrblasInterface.o \
  			lsqrModule.o	AnalysisGLI.o	lsqrCheckModule.o	lsqrTestModule.o   	\
			smooth2.o  	illumination_balance.o		\
			tomo.o

all: lsqrDataModule		lsqrblas 	global	AnalysisGLI lsqrblasInterface lsqrModule 	\
	 lsqrCheckModule	lsqrTestModule 	smooth2  illumination_balance 	\
	 tomo tomography

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

global:
	$(FC) ${FFLAGS} -c global.f90

AnalysisGLI:
	$(FC) ${FFLAGS} -c AnalysisGLI.f90

smooth2:
	$(CC) $(CFLAGS) -c smooth2.c

illumination_balance:
	$(CC) $(CFLAGS) -c illumination_balance.c	

tomo:
	${FC} ${FFLAGS} -c QTomo_lsqr_wpy.f90 -o tomo.o
#	${FC} ${FFLAGS} -c RayTomo_lsqr_170930.f90 -o tomo.o
#	${FC} ${FFLAGS} -c Forward.f90 -o tomo.o

tomography: ${files}
	${FC} ${FFLAGS} ${files} -o Q_Ray_Inv.exe 
#	${FC} ${FFLAGS} ${files} -o RayInv.exe 
#	${FC} ${FFLAGS} ${files} -o RayInv.exe -llfFunC


clean:
	\rm -f *.o *.mod

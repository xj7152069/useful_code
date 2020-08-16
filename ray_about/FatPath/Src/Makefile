  COMPLIER=
  LIBPATH=
  
  FC      = ${COMPLIER}mpiifort
  FFLAGS  = -g -O0

  CC      = ${COMPLIER}mpiicc
  CFLAGS  = -g -O0
  
  LIB     = ${LIBPATH} -lmkl_intel_ilp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm -ldl

  files =	lsqrDataModule.o	lsqrblas.o	global.o	lsqrblasInterface.o \
  			lsqrModule.o	AnalysisGLI.o	lsqrCheckModule.o	lsqrTestModule.o \
			smooth2.o    openCfile.o   	illumination_balance.o	calculate_fatray_kernels.o \
			tomo.o LinearInter.o RayingTracingFat.o

all: lsqrDataModule	lsqrblas global	AnalysisGLI lsqrblasInterface lsqrModule 	\
	 lsqrCheckModule lsqrTestModule smooth2 openCfile illumination_balance calculate_fatray_kernels	\
	 LinearInter RayingTracingFat tomo tomography

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

openCfile:
	$(CC) $(CFLAGS) -c ./include/openCfile.c

illumination_balance:
	$(CC) $(CFLAGS) -c illumination_balance.c	

tomo:
	${FC} ${FFLAGS} -c FatTomo_lsqr_0102.f90 -o tomo.o

calculate_fatray_kernels:
	$(CC) $(CFLAGS) -c calculate_fatray_kernels.c	

LinearInter:
	$(FC) $(CFLAGS) -c LinearInter.f90

RayingTracingFat:
	$(FC) $(CFLAGS) -c RayingTracingFat.f90


tomography: ${files}
	${FC} ${FFLAGS} ${files} -o BeamInv.exe

clean:
	\rm -f *.o *.mod

# before including this file,
# the ProjRoot var needs to be set
BinDir     = ${ProjRoot}/bin
LibDir     = ${ProjRoot}/lib ${HOME}/pengdev/lib
IncDirList = ${ProjRoot}/include ${inc_peng}
ObjDir     = ${ProjRoot}/obj

CXX=icpc

ifeq ($(CXX),g++)
	flag_warn = -Wno-unused-result
endif

ifeq ($(CXX),icpc)
	flag_warn = -Wno-unused -wd2218,381
endif

FLAG_DEBUG =
flag_profiler=
#FLAG_DEBUG = -DDEBUG
#flag_profiler= -pg

CPP_CMPL=mpicxx -cxx=${CXX}
C_CMPL  =mpicc -cc=${CC}
F_CMPL  =gfortran
LINKER  =$(CPP_CMPL)

flag_arma 	  = -DARMA_ALLOW_FAKE_GCC -DARMA_DONT_USE_WRAPPER
C_CMPL_FLAG   = -fopenmp -O3 ${FLAG_DEBUG}
CPP_CMPL_FLAG = -std=c++11 -fopenmp -fmax-errors=2 -O3 \
				$(flag_profiler) $(flag_warn) ${flag_arma} ${FLAG_DEBUG}
F_CMPL_FLAG   = -cpp -O3 ${FLAG_DEBUG}

# list of libs in use (with no -l)
LIB_LIST	  = pengpack m mkl_rt
# list of libs in use (with -l)
LIB_FLAG      = $(addprefix -L,${LibDir}) $(addprefix -l,${LIB_LIST})
# list of inc
INC_FLAG      = $(addprefix -I,${IncDirList})
# LINKER_FLAG!
LINKER_FLAG   = -std=c++11 -fopenmp -fmax-errors=2 -O3 \
				$(flag_profiler) ${LIB_FLAG} $(flag_warn)

ifeq ($(JOBS),)
  JOBS := $(shell grep -c ^processor /proc/cpuinfo 2>/dev/null)
  ifeq ($(JOBS),)
    JOBS := 1
  endif
endif


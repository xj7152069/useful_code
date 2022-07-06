# before including this file,
# the ProjRoot var needs to be set
BinDir     = ${ProjRoot}/bin
LibDir     =
IncDirList = ${ProjRoot}/include
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

#701 server, 192.168.1.200
mpich2intel=/home/apps/intel2018/compilers_and_libraries_2018.1.148/linux/mpi/intel64
#222 server, 192.168.1.222
#mpich2intel=/home/apps/mpich214intel2017
F_CMPL  = $(mpich2intel)/bin/mpif90
#C_CMPL  = $(mpich2intel)/bin/mpicc
C_CMPL  = $(mpich2intel)/bin/mpiicc
CPP_CMPL= $(mpich2intel)/bin/mpicxx
LibDir	= $(mpich2intel)/lib
XJ_ADD_LibDir   = /data3/xj/opwMig3d.202206/lib/
# list of libs in use (with no -l)
LIB_LIST	= m
XJ_ADD_LIB_LIST	= timespacelslrt2d

#Ocean-313 server, 100.64.164.255
#mpich2intel=/home/apps/intelOneAPI/mpi/2021.5.0
#F_CMPL  = $(mpich2intel)/bin/mpiifort
#C_CMPL  = $(mpich2intel)/bin/mpiicc
#CPP_CMPL= $(mpich2intel)/bin/mpicpc
#LibDir	= $(mpich2intel)/lib
# list of libs in use (with no -l)
#LIB_LIST	  =

#F_CMPL  =ifort
#LINKER  =$(CPP_CMPL)
LINKER  =$(C_CMPL)

flag_arma 	  = -DARMA_ALLOW_FAKE_GCC -DARMA_DONT_USE_WRAPPER
#C_CMPL_FLAG   = -std=c99 -qopenmp -check-pointers=rw ${FLAG_DEBUG}
C_CMPL_FLAG   = -std=c99 -qopenmp ${FLAG_DEBUG}
#C_CMPL_FLAG   = -std=c99 -qopenmp -O3 ${FLAG_DEBUG}
CPP_CMPL_FLAG = -std=c++11 -fopenmp -fmax-errors=2 -O3 \
				$(flag_profiler) $(flag_warn) ${flag_arma} ${FLAG_DEBUG}
F_CMPL_FLAG   = -O3 ${FLAG_DEBUG}


# list of libs in use (with -l)
LIB_FLAG      = $(addprefix -L,${LibDir}) $(addprefix -l,${LIB_LIST})
XJ_ADD_LIB_FLAG	= $(addprefix -L,${XJ_ADD_LibDir}) $(addprefix -l,${XJ_ADD_LIB_LIST}) -lstdc++ -lirc

# list of inc
INC_FLAG      = $(addprefix -I,${IncDirList})
# LINKER_FLAG!
#LINKER_FLAG   = -std=c++11 -fopenmp -fmax-errors=2 -O3 \
#				$(flag_profiler) ${LIB_FLAG} $(flag_warn)
LINKER_FLAG   = -std=c99 -xhost -qopenmp -O3 \
				$(flag_profiler) ${LIB_FLAG} \
				$(flag_profiler) ${XJ_ADD_LIB_FLAG}

ifeq ($(JOBS),)
  JOBS := $(shell grep -c ^processor /proc/cpuinfo 2>/dev/null)
  ifeq ($(JOBS),)
    JOBS := 1
  endif
endif


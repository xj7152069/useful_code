#makefile for c c++ fotran
#BinDir =$(bin_peng)
divider="--------------------------------------------------------------------------------"
BinDir =../bin
LibDir =../lib
IncDir =../lib
ObjDir =../obj
BackDir=~/
Date=$(shell date +%F_time%H)
#--------------------------------------#
Projn   =ELASTIC_WAVE_FD
#--------------------------------------#
CPP_CMPL_0=icpc
CPP_CMPL=$(CPP_CMPL_0) -O3
#CPP_CMPL=icc -std=c++11 -O2
C_CMPL  =icc -O3
F_CMPL  =gfortran	-O3 -cpp
LINKER  =$(CPP_CMPL)

C_CMPL_FLAG   = -Wno-unused-result
#CPP_CMPL_FLAG = -Wno-unused-result -fmax-errors=5
CPP_CMPL_FLAG = -fmax-errors=5
F_CMPL_FLAG   =
LINKER_FLAG   =
LINKLIBS	  = -lm

CSrc   =$(wildcard *.c)
CppSrc =$(wildcard *.cpp)
F90Src =$(wildcard *.f90)
F77Src =$(wildcard *.f)

CObj0  =$(patsubst %.c,%.o,${CSrc})
CppObj0=$(patsubst %.cpp,%.o,${CppSrc})
F90Obj0=$(patsubst %.f90,%.o,${F90Src})
F77Obj0=$(patsubst %.f,%.o,${F77Src})

CObj   =$(addprefix $(ObjDir)/,$(CObj0))
CppObj =$(addprefix $(ObjDir)/,$(CppObj0))
F90Obj =$(addprefix $(ObjDir)/,$(F90Obj0))
F77Obj =$(addprefix $(ObjDir)/,$(F77Obj0))

OBJ    =$(CObj) $(CppObj) $(F90Obj) $(F77Obj)
EXE    =$(BinDir)/$(Projn).e
BACK   =$(BackDir)/$(Projn)_back_$(Date).tgz
Back_Content=$(shell basename $(shell pwd))
ex_tar ="--exclude=*.sw? --exclude=*.log --exclude=*.[oe]"

all:
	mkdir -p $(ObjDir)
	make -j 10 $(EXE)

$(EXE):$(OBJ)
	@echo $(divider)
	@echo '			Linking ---> ' $(EXE)
	@echo $(divider)
	$(LINKER) $^ -o $(EXE)  $(LINKLIBS) $(LINKER_FLAG) -fopenmp
$(ObjDir):
	mkdir -p $(ObjDir)
$(CppObj):$(ObjDir)/%.o:%.cpp
	@echo $(divider)
	@echo '			Making '$<'  --> *.o'
	@echo $(divider)
	$(CPP_CMPL) -c  -I$(IncDir) -I./inc $< -o $@ $(CPP_CMPL_FLAG) -fopenmp
$(CObj):$(ObjDir)/%.o:%.c
	@echo $(divider)
	@echo '			Making '$<'  --> *.o'
	@echo $(divider)
	$(C_CMPL)   -c  -I$(IncDir) -I./inc  $< -o $@ $(C_CMPL_FLAG)
$(F90Obj):$(ObjDir)/%.o:%.f90
	@echo $(divider)
	@echo '			Making *.f90  --> *.o'
	@echo $(divider)
	$(F_CMPL)  -c  $< -o $@
$(F77Obj):$(ObjDir)/%.o:%.f
	@echo $(divider)
	@echo '			Making *.f  --> *.o'
	@echo $(divider)
	$(F_CMPL)  -c  $< -o $@

backup:
	@#make clean;
	mkdir -p $(BackDir)
	cd .. && tar $(ex_tar) -czvf $(BACK) $(Back_Content)
	ls -hl $(BACK)
b:backup

clean:
	@#CREAT EMPTY files
	@touch $(EXE)
	@touch $(OBJ)
	rm  $(EXE) $(OBJ)
	@pwd; echo $(Date)

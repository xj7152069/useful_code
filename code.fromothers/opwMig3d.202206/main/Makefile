#makefile for   [projname]
#--------------------------------------#
#ProjRoot = $(shell bash get_proj_root_by_gitdir.sh)
ProjRoot = ..
#--------------------------------------#
include proj_config.mk

#======= Obj  =======
CppSrc       = $(wildcard *.cpp)
CSrc         = $(wildcard *.c) $(wildcard ${ProjRoot}/include/*.c)
F90Src       = $(wildcard *.f90)
CppObj0      = $(patsubst %.cpp,%.o,${CppSrc})
CObj0        = $(patsubst %.c,%.o,${CSrc})
F90Obj0      = $(patsubst %.f90,%.o,${F90Src})

CppObj       = $(addprefix $(ObjDir)/,$(CppObj0))
CObj         = $(addprefix $(ObjDir)/,$(CObj0))
F90Obj       = $(addprefix $(ObjDir)/,$(F90Obj0))

OBJ          = $(CppObj) $(F90Obj) ${CObj}
# the non-main OBJ
OBJ_NonMain  = $(filter-out $(ObjDir)/main%, ${OBJ})
OBJ_Main     = $(filter $(ObjDir)/main%, ${OBJ})
#=== multi-EXEs =====
EXE          = $(filter-out $(BinDir)/_%, \
			   $(patsubst $(ObjDir)/main_%.o,$(BinDir)/%.e,${OBJ_Main}))

all:
	#@mkdir -p obj bin
	make $(EXE) -j ${JOBS};
	#mktag

$(EXE):$(BinDir)/%.e:$(ObjDir)/main_%.o ${OBJ_NonMain}
	@echo '			Linking ---> ' $@
	$(LINKER) $^ -o $@  $(LINKER_FLAG)

$(CppObj):$(ObjDir)/%.o:%.cpp
	@echo '			Compiling '$<'  --> *.o'
	$(CPP_CMPL) $(INC_FLAG) -c $< -o $@ ${CPP_CMPL_FLAG}

$(CObj):$(ObjDir)/%.o:%.c
	@echo '			Compiling '$<'  --> *.o'
	$(C_CMPL) -c  $(INC_FLAG) $< -o $@ $(C_CMPL_FLAG)

$(F90Obj):$(ObjDir)/%.o:%.f90
	@echo '			Compiling '$<'  --> *.o'
	$(F_CMPL)  -c  $< -o $@ $(F_CMPL_FLAG)

lib:
	#Xiangjian's C++ code.
	icpc -O3 -c -qopenmp ./cpp_code/src/slantstack3d.omp.cpp -o ../lib/timespacelslrt2d.o -xhost
	#g++ -O3 -c ./cpp_code/src/slantstack3d.cpp -o ../lib/timespacelslrt2d.o
	#g++ -O3 -c -fopenmp ./cpp_code/src/slantstack3d.omp.cpp -o ../lib/timespacelslrt2d.o
	ar crs ../lib/libtimespacelslrt2d.a ../lib/timespacelslrt2d.o

clean:
	@#CREAT EMPTY files
	#@mkdir -p obj bin
	@touch $(EXE) $(OBJ) tags
	rm  $(EXE) $(OBJ) tags
	rm ../lib/*.a ../lib/*.o
	@pwd; echo $(Date)

c:clean

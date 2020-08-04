
#include ../config/Makefile.inc


CUR_DIR		= $(shell pwd)

HDS_DIR		= $(CUR_DIR)/hds
BIN_DIR		= $(CUR_DIR)/bin
SRC_DIR		= $(CUR_DIR)/src
OBJ_DIR		= $(CUR_DIR)/obj
INC_DIR		= $(CUR_DIR)/include

#vpath	%.c $(INC_DIR)

#compiler
#FC			= $(MPICH)/bin/mpif90
#CC			= $(MPICH)/bin/mpicc
#FC			=  /home/apps/mpich214intel2017/bin/mpif90
FC			= mpif90
CC			= mpicc

FFLAG		= -f90=ifort  -O3
CFLAG       = -std=c99 -O3  -lm
#CFLAG		= -std=c99 -O3  -lm  -g   -pg
FFLAGS		= $(FFLAG)  -qopenmp 
#FFLAGS		= $(FFLAG)  -qopenmp -check bounds -traceback  -g
CFLAGS		= $(CFLAG)  -qopenmp 
#LIBF		= -f90=ifort -std=c99 -lm -lfftw3f -O3 -w -Wall -g
	

# for main program (f90).
#fb_Seiscope_Toolbox_Dir=$(PackageDir)/fb_Seiscope_Toolbox
#INC_Seiscope = -I $(fb_Seiscope_Toolbox_Dir)/src/COMMON/include
#FLIB_Seiscope	=  -L $(fb_Seiscope_Toolbox_Dir)/lib -lSEISCOPE_OPTIM


#include
#INCLUDE	    = -I$(HDS_DIR) -I$(INC_DIR) -I$(OBJ_DIR)  $(INC_Seiscope) -I$(MPICH_INC) -I$(FFTW3_INC) 
INCLUDE	    = -I$(HDS_DIR) -I$(INC_DIR) -I$(OBJ_DIR)  


#lib
LIBS		= -lfftw3f_threads -lfftw3
LIBF	= -lifcore
LIB		=  -lm
#LIB		= $(LFFTW3F_LIB) $(LFFTW3_LIB) $(LIBF) -lm

CLIB	= -L /home/apps/intel/compilers_and_libraries/linux/lib/intel64 -liomp5
#FLIB	= $(FLIB_Seiscope)  $(LIB)
FLIB	=  $(LIB)


#code
#INCLUDE 
inc_src =	$(INC_DIR)/global.f90  \
			$(INC_DIR)/rtm_2d_single_shot_extrapolation_module.f90 \
			$(INC_DIR)/rtm_2d_single_shot_extrapolation_allocate.f90 \
			$(INC_DIR)/Extrapolation_2D_Without_PML_One_Step.f90 \
			$(INC_DIR)/Extrapolation_2D_ISO_RTM_One_Step_v1.f90 \
			$(INC_DIR)/Extrapolation_2D_ISO_RTM_One_Step_v2.f90 \
			$(INC_DIR)/Form_Coe_extrapolation_2d.f90 \
			$(INC_DIR)/Form_Dynamic_Filename.f90 \
			$(INC_DIR)/Get_Current_Shot_Velocity_2d.f90 \
			$(INC_DIR)/index_shotgathers_and_write_info.f90 \
			$(INC_DIR)/Input_Shot_Gather_2d.f90 \
			$(INC_DIR)/local_to_all_image_result_2d.f90 \
			$(INC_DIR)/merge_shot_files.f90 \
			$(INC_DIR)/Read_shotgathers_info.f90 \
			$(INC_DIR)/Source_and_Receiver_Position_2d.f90 \
			$(INC_DIR)/Trace_Sinc_Interp.f90 \
			$(INC_DIR)/Wavelet_Forming.f90 \
			$(INC_DIR)/Cor_Imaging_Condition_2d.f90 \
			$(INC_DIR)/Laplace_time_2d.f90 \
			$(INC_DIR)/maxinum_1d.f90  \
			$(INC_DIR)/coefficient_1st.f90  \
			$(INC_DIR)/coefficient_2nd.f90 \
			$(INC_DIR)/process_directwave_2d.f90


finc_src   = $(filter %.f90,$(inc_src))
finc_src_o = $(finc_src:.f90=.o)
finc_obj_name = $(notdir $(finc_src_o))
finc_obj = $(addprefix $(OBJ_DIR)/,$(finc_obj_name))
$(finc_obj): $(finc_src)
	$(FC) $(FFLAGS) $(INCLUDE) -c  $(filter %$(*F).f90,$(finc_src))  -o $@   -module $(OBJ_DIR) 
	
cinc_src   = $(filter %.c,$(inc_src))
cinc_src_o = $(cinc_src:.c=.o)
cinc_obj_name = $(notdir $(cinc_src_o))
cinc_obj = $(addprefix $(OBJ_DIR)/,$(cinc_obj_name))
$(cinc_obj): $(cinc_src)
	$(CC) $(INCLUDE) -c  $(filter %$(*F).c,$(cinc_src))  $(CFLAGS) -o $@ 

inc_obj = $(finc_obj) $(cinc_obj)


######################################################################
#model2d
model2d_main_src  = \
			  $(SRC_DIR)/model2d/model_2d_single_shot.f90 \
			  $(SRC_DIR)/model2d/model_2d_multi_shot.f90 \
			  $(SRC_DIR)/model2d/model_2d_main.f90 

model2d_main_src_o = $(model2d_main_src:.f90=.o)
model2d_main_obj_name = $(notdir $(model2d_main_src_o))
model2d_main_obj = $(addprefix $(OBJ_DIR)/,$(model2d_main_obj_name))
$(model2d_main_obj): $(model2d_main_src)
	$(FC) $(FFLAGS) $(INCLUDE) -c  $(filter %$(*F).f90,$(model2d_main_src))  -o $@   -module $(OBJ_DIR)  $(FLIB)


model2d_obj =  $(inc_obj)  $(model2d_main_obj)
model2d_bin	= $(BIN_DIR)/model2d

model2d: $(model2d_obj)
	$(FC)  $(FFLAGS)  $(INCLUDE) -o  $(model2d_bin) $(model2d_obj)  $(FLIB)


######################################################################
#rtm2d
rtm2d_main_src  =  \
			$(SRC_DIR)/rtm2d/rtm_2d_multi_shot.f90 \
			$(SRC_DIR)/rtm2d/rtm_2d_single_shot.f90 \
			$(SRC_DIR)/rtm2d/rtm_2d_main.f90 

rtm2d_main_src_o = $(rtm2d_main_src:.f90=.o)
rtm2d_main_obj_name = $(notdir $(rtm2d_main_src_o))
rtm2d_main_obj = $(addprefix $(OBJ_DIR)/,$(rtm2d_main_obj_name))
$(rtm2d_main_obj): $(rtm2d_main_src)
	$(FC) $(FFLAGS) $(INCLUDE) -c  $(filter %$(*F).f90,$(rtm2d_main_src))  -o $@   -module $(OBJ_DIR)  $(FLIB)

rtm2d_obj =  $(inc_obj)  $(rtm2d_main_obj)
rtm2d_bin = $(BIN_DIR)/rtm2d

rtm2d: $(rtm2d_obj)
	$(FC)  $(FFLAGS)  $(INCLUDE) -o  $(rtm2d_bin) $(rtm2d_obj)  $(FLIB)



######################################################################




######################################################################
.PHONY:clean
clean:
	rm $(OBJ_DIR)/*.o  $(OBJ_DIR)/*.mod



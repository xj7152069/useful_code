#!/bin/bash

set -x
#lenovo-222.
MPI_F90=/home/apps/mpich214intel2017/bin/mpif90
#FC=ifort
FC=gfortran
CompileFlag_F="-fopenmp -O3 -c "  # for gfortran

#wpi-200.
#MPI_F90=/home/fb/apps/mpich214intel2019/bin/mpif90

#Ocean-313.
#MPI_F90=/home/apps/mpich214IntelOneAPI/bin/mpif90
#MPI_F90=/home/apps/intelOneAPI/mpi/2021.5.0/bin/mpiifort

#static compling.
#CompileFlag_F="-c -check bounds -traceback"
#CompileFlag_F="-qopenmp -O3 -c "  # for ifort
#CompileFlag_F="-qopenmp -c "

$FC $CompileFlag_F sub_FFT1D.f -o FFT1D.o
$FC $CompileFlag_F sub_FFT2D_OMP.f -o FFT2D_OMP.o
#$FC $CompileFlag_F sub_subroutines.f -o subroutines.o
#$FC $CompileFlag_F sub_extrapolation_WXFD_visco.f -o extrapolation_WXFD_visco.o
$FC $CompileFlag_F sub_extrapolation_WXFD_acoustic.f -o extrapolation_WXFD_acoustic.o
$FC $CompileFlag_F sub_extrapolation_SSF_acoustic.f -o extrapolation_SSF_acoustic.o
$FC $CompileFlag_F sub_readOPW_Data.f90 -o readOPW_Data.o
$FC $CompileFlag_F sub_readParMig3d.f90 -o readParMig3d.o
exit
main=main_opwPreStackTimeMigration3D.v23.f90
$MPI_F90 $CompileFlag_F $main -o main_opwpstm3d.o

#linking
EXEC=../../bin/opwmig3d.e
#$MPI_F90 -check bounds -traceback *.o -o $EXEC
$MPI_F90 -qopenmp -O3 *.o -o $EXEC
#$MPI_F90 -qopenmp *.o -o $EXEC

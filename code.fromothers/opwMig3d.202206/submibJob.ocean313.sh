#!/bin/bash
set -x

dir=`pwd`
nproc=33
MPIEXEC=/home/apps/intelOneAPI/mpi/2021.5.0/bin/mpiexec
OPWD_par=${dir}/par.ocean313/makeCMPgatherIndex.Luojia3D.allfiles.par
PSTM_par=${dir}/par.ocean313/opwPSTM3D.Luojia3D.par

job=cmpIndex
bin_index=${dir}/bin/makeCMPgatherIndex_v2.e
par=${OPWD_par}
#$bin_index $par
#time $bin_index $par 1>${job}.log.txt 2>${job}.err.txt &

job=opwd
bin_opwd=${dir}/bin/3d_Phr_OPWD_mpi_v4.e
par=${OPWD_par}
export OMP_NUM_THREADS=32
nproc=9
#time $MPIEXEC -n $nproc $bin_opwd $par
#(time $MPIEXEC -n $nproc $bin_opwd $par) 1>log.${job}.txt 2>err.${job}.txt &

job=pstm3d
bin_mig3d=${dir}/bin/opwmig3d.e
data_dir=/home/fb/data/syntheticData
#impulse reponse.
PSTM_par=${data_dir}/par.ocean313/opwPSTM3D.syntheticData3D.impulseResponse.par
#true OPW data migration.
PSTM_par=${data_dir}/par.ocean313/opwPSTM3D.syntheticData3D.trueOPWdata.par
par=${PSTM_par}
nproc=9
(time $MPIEXEC -n $nproc $bin_mig3d $par) 1>log.${job}.txt 2>err.${job}.txt &

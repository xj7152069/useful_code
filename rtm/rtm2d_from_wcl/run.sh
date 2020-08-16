#!/bin/bash
set -x

export OMP_NUM_THREADS=45
#mpiexec=/home/apps/mpich214intel2017/bin/mpiexec
mpiexec=mpiexec
dir=.
hostfile=$dir/nodefile

bin=$dir/bin/model2d
par=./par/model2d/parameter_model_2d
log=./log/log_model2d
err=./log/err_model2d

#bin=$dir/bin/rtm2d
#par=./par/rtm2d/parameter_rtm_2d
#log=./log/log_rtm2d
#err=./log/err_rtm2d


time $mpiexec  -f  $hostfile  $bin $par  1>$log  2>$err  &

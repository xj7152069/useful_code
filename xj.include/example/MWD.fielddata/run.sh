rm -f ./mpi.demultiple3d.fielddata.out
#compile:
time mpic++ -O3 demultiple3d.fielddata.continue.mpi.cpp -o mpi.demultiple3d.fielddata.out -llapack -lblas -lgfortran -std=c++11 -lpthread -fopenmp

#run:
nohup time mpirun  -machinefile ./nodefile ./mpi.demultiple3d.fielddata.out ../data/r2413_pzproc.su ./result/r2413.demultiple.su 323 128 4097 50 25 0.002 10 10 0.01 32&

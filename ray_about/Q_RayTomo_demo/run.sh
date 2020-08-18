#mpiexec -np 5 -machinefile nodefile ./Q_Ray_Inv.exe >log&
mpirun -np 5 ./Q_Ray_Inv.exe >log&

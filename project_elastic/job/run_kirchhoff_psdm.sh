gfortran -O2 -c ../include/xj.ray.f90/caltime.sigle.f90 -o ../lib/caltime.sigle.o #生成fortranCode.o
g++ -O2 -c main.cpp -o ../lib/main.o #生成main.o
g++ -O2 -o ../obj/main_kirchhoff_psdm.exe ../lib/*.o -lgfortran #生成program test
#time ./main.exe  #运行程序




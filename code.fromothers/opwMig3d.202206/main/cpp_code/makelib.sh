#make lib for LS-LRT
g++ -O3 -c ./code/slantstack3d.cpp -o ./timespacelslrt2d.o
ar crs ./libtimespacelslrt2d.a ./timespacelslrt2d.o
rm -f ./timespacelslrt2d.o

 c++ -O2 main.cpp -o main.out\
 -L ./include/lib -llapack -lblas\
 -lgfortran -lpthread
  #c++ -O2 obn_get_com_recvgather.cpp -o getcr.out -llapack -lblas -lgfortran -lpthread

 #c++ -O2 main_radon3d_orig.cpp -o radon.out -llapack -lblas -lgfortran -std=c++11 -lpthread
 #c++ -O2 demultiple.cpp -o demultiple.out -llapack -lblas -lgfortran -std=c++11 -lpthread
 #c++ -O2 main.coding.cpp -o main.coding.out -llapack -lblas -lgfortran -std=c++11 -lpthread
 #
 #c++ -O2 model_complex.cpp -o model.out -llapack -lblas -lgfortran -std=c++11 -lpthread

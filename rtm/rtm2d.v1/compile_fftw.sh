#2017-05-19
#mpif90 2d_acoustic_rtm_pml_cpu_v1.f90   -o  2d_acoustic_rtm_pml_cpu_v1.e  -O3  -fopenmp -I/usr/include  -lfftw3f_threads -lfftw3f -lm
#mpif90 2d_acoustic_rtm_pml_cpu.f90   -o  2d_acoustic_rtm_pml_cpu.e  -O3  -fopenmp -I/usr/include  -lfftw3f_threads -lfftw3f -lm

#mpif90 2d_acoustic_rtm_pml_cpu_v1_test.f90   -o  2d_acoustic_rtm_pml_cpu_v1_test.e  -O3  -fopenmp -I/usr/include  -lfftw3f_threads -lfftw3f -lm

#2017.12-05
mpif90 -f90=ifort ./src/2d_acoustic_rtm_pml_cpu_v1.f90   -o  ./bin/2d_acoustic_rtm_pml_cpu_v1.e  -O3  -qopenmp -I/usr/include  -lfftw3f_threads -lfftw3f -lm
#mpif90 -f90=ifort ./src/2d_acoustic_rtm_pml_cpu_v1.f90   -o  ./bin/2d_acoustic_rtm_pml_cpu_v1.e  -O3  -qopenmp -I/usr/include  -lfftw3f_threads -lfftw3f -lm

#2019-01-13

#mpif90 -f90=ifort  ./src/2d_acoustic_rtm_pml_cpu_optimal_stack.f90   -o  \
#				   ./bin/2d_acoustic_rtm_pml_cpu_optimal_stack.e \
#					-O3  -qopenmp -I/usr/include  -lfftw3f_threads -lfftw3f -lm

!*=======================================================================*
!*                                                                       *
!*                            2D  RTM                                    * 
!*                                                                       *
!*=======================================================================*
!*                                                                       *
!*                      ALL RIGHTS RESERVED                              *
!*           WAVE PHENOMENA AND INVERSION IMAGING GROUP (WPI)            *
!*      BY SCHOOL OF OCEAN AND EARTH SCIENCE TONGJI UNIVERSITY           *
!*                                                                       *
!*=======================================================================*
!*            REMARKS:                                                   *
!*            Written By : Wu Chengliang                                 *
!*                    05.01, 2020 -- 05.05, 2020                         *
!*            Modify  By : Wu Chengliang                                 *
!*                    05.05, 2020 -- 05.05, 2020                         *
!*			 Version  :	NO.1										     *
!*																		 *
!*=======================================================================*
!*Technical Reference:                                                   *
!*版本说明
!*
!*
!*
!*************************************************************************


!*************************************************************************
!*						The main program begin	                         *
!*************************************************************************
	program rtm_2d

	use global
	use omp_lib
	implicit none
	include 'mpif.h'

	!Data dictionary: declare variable types, definitons, & units
	integer	::	mpi_np, myid
	integer     ::  para_int(para_int_fnum)
	integer*8   ::  para_long(para_long_fnum)
	real        ::  para_float(para_float_fnum)
	real*8      ::  para_double(para_double_fnum)
	character(len=para_char_flen)	::	para_char(para_char_fnum)


	!for mig
	character(len=para_char_flen)	::	fn_par
	character(len=para_char_flen)	::	fn_vel
	character(len=para_char_flen)	::	fn_cs, fn_grad
	character(len=para_char_flen)  ::  fn_shot_info
	integer	::	nvx, nvz
	real	::	dvx, dvz
	integer	::	iflag_shot_info
	integer	::	ns_start, ns_interval, ns_end
	integer	::	lt,	lt_real
	real	::	dt, dt_real
	real	::	depth
	real	::	dx, dz
	real	::	aperture_x_l, aperture_x_r
	real	::	aperture_z_u, aperture_z_d
	real	::	boundary_x_l, boundary_x_r
	real	::	boundary_z_u, boundary_z_d
	real	::	cvx_initial, cvz_initial
	real	::	fmain
	integer	::	ntrace_max
	integer	::	iflag_mig, iflag_mig_iter
	integer	::	ierr
	integer ::	master
	integer	::	threads
	character(MPI_MAX_PROCESSOR_NAME)   ::  processorname
	integer ::  processorlen
	integer ::  processorpid

	!===================== Initializing the mpi =========================
	CALL MPI_INIT( IERR )
	CALL MPI_COMM_SIZE( MPI_COMM_WORLD, MPI_NP, IERR )
    CALL MPI_COMM_RANK( MPI_COMM_WORLD, MYID, IERR )

	!********************************************************************
	!*---------Inputing the parameters for running the code ------------*
	master=0
	if(myid .eq. master )then
		CALl Getarg(1,fn_par)
	endif
	CALL MPI_BCAST(FN_PAR, para_char_flen, MPI_CHARACTER, MASTER, MPI_COMM_WORLD, IERR)
    CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)

	!********************************************************************
	threads=omp_get_max_threads()
	CALL omp_set_num_threads(threads)

	if(myid .eq. master) then
		!info
		!printf parameter
		print*, 'The OPENMP threads is:', threads
	endif

	CALL MPI_Get_processor_name(processorname, processorlen, IERR)
	print*, 'myid', myid, ' processorname : ', trim(processorname)


	!***********************************************************************
    CALL rtm2d_ReadPar(fn_par, &
			fn_vel, fn_cs, fn_grad, &
			iflag_shot_info, fn_shot_info, &
			nvx, nvz, & 
			dvx, dvz, &
			ns_start, ns_interval, ns_end, &
			lt, lt_real, &
			dt, dt_real, &
			depth, &
			dx, dz, &
			aperture_x_l, aperture_x_r, &
			aperture_z_u, aperture_z_d, &
			boundary_x_l, boundary_x_r, &
			boundary_z_u, boundary_z_d, &
			cvx_initial, cvz_initial, &
			fmain, ntrace_max)



	para_char(1)    =	fn_vel
	para_char(2)	=	fn_cs
	para_char(3)	=	fn_grad
	para_char(4)	=	fn_shot_info

	para_int(1)		=	nvx	
	para_int(2)		=	nvz
	para_int(3)		=	lt		
	para_int(4)		=	lt_real	
	para_int(5)		=	ns_start	
	para_int(6)		=	ns_interval	
	para_int(7)		=	ns_end	
	para_int(8)		=	iflag_shot_info
	para_int(9)		=	ntrace_max


	para_float(1)	=	fmain		
	para_float(2)	=	dt		
	para_float(3)	=	dt_real		
	para_float(4)	=	aperture_x_l	
	para_float(5)	=	aperture_x_r
	para_float(6)	=	aperture_z_u	
	para_float(7)	=	aperture_z_d
	para_float(8)	=	boundary_x_l
	para_float(9)	=	boundary_x_r
	para_float(10)	=	boundary_z_u
	para_float(11)	=	boundary_z_d
	para_float(12)	=	depth		
	para_float(13)	=	dx		
	para_float(14)	=	dz	
	para_float(15)	=	dvx			
	para_float(16)	=	dvz			
	para_float(17)	=	cvx_initial
	para_float(18)	=	cvz_initial


	!********************************************************************
    if( myid .eq. master ) then
		write(*,*) '********************************************'
		write(*,*)'        Run the parameter file is :'
		write(*,*)'     ',trim(fn_par)
		write(*,*)'		The parameter variables in this file is :'
		write(*,*)'		input shot gathers filename ', trim(fn_cs)
		write(*,*)'		grad filename ', trim(fn_grad)
		write(*,*)'		shot info filename ', trim(fn_shot_info)
		write(*,*)'		nvx, nvz',nvx, nvz
		write(*,*)'		dvx, dvz',dvx, dvz
		write(*,*)'		lt, lt_real',lt, lt_real
		write(*,*)'		dt, dt_real',dt, dt_real
		write(*,*)'		ns_start, ns_interval, ns_end', ns_start, ns_interval, ns_end
		write(*,*)'		aperture_x_l, aperture_x_r, '
		write(*,*) 		aperture_x_l, aperture_x_r
		write(*,*)'		aperture_z_u, aperture_z_d', aperture_z_u, aperture_z_d
		write(*,*)'		boundary_x_l, boundary_x_r, '
		write(*,*) 		boundary_x_l, boundary_x_r
		write(*,*)'		boundary_z_u, boundary_z_d', boundary_z_u, boundary_z_d
		write(*,*)'		cvx_initial, cvz_initial', cvx_initial, cvz_initial
		write(*,*)'		dvx, dvz', dvx, dvz
		write(*,*)'		depth', depth
		write(*,*)'		fmain', fmain
        write(*,*) '********************************************'
    endif
	!===================== Begin the main program =======================
    if( myid .eq. master ) then
		write(*,*) '********************************************'
		write(*,*) '         2D rtm                             '
        write(*,*) '                 begins                     '
        write(*,*) '********************************************'
    endif

    CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)

    CALL rtm_2d_multi_shot_mpi(para_int, para_long, para_float, para_double, para_char, &
				mpi_np, myid)


	!=================== End the main program ===========================
    if( myid .eq. master ) then
		write(*,*) '********************************************'
		write(*,*) '         2D rtm                             '
		write(*,*) '                 end                        '
        write(*,*) '********************************************'
    endif

	!*===============  Closing the MPI working environment ==============
    CALL MPI_FINALIZE(IERR)
	!********************************************************************
	!*                   The Main Program end !                         *
	!********************************************************************
	stop
    end	program
!***********************************************************************

    subroutine rtm2d_ReadPar(fn_par, &
			fn_vel, fn_cs, fn_grad, &
			iflag_shot_info, fn_shot_info, &
			nvx, nvz, & 
			dvx, dvz, &
			ns_start, ns_interval, ns_end, &
			lt, lt_real, &
			dt, dt_real, &
			depth, &
			dx, dz, &
			aperture_x_l, aperture_x_r, &
			aperture_z_u, aperture_z_d, &
			boundary_x_l, boundary_x_r, &
			boundary_z_u, boundary_z_d, &
			cvx_initial, cvz_initial, &
			fmain, ntrace_max)

	use global
	implicit none
	!Data dictionary: declare variable types, definitons, & units
	!Dummy Variables
	character(len=para_char_flen)	::	fn_par
	character(len=para_char_flen)	::	fn_vel
	character(len=para_char_flen)	::	fn_cs, fn_grad
	character(len=para_char_flen)   ::	fn_shot_info
	integer ::	ierr
	integer ::	nvx, nvz
	real	::	dvx, dvz
	integer ::	ns_start, ns_interval, ns_end
	integer ::	lt, lt_real
	real	::	dt, dt_real
	real	::	depth, dz
	real	::	dx
	real	::	aperture_x_l, aperture_x_r
	real	::	aperture_z_u, aperture_z_d
	real	::	boundary_x_l, boundary_x_r
	real	::	boundary_z_u, boundary_z_d
	real	::	cvx_initial, cvz_initial
	real	::	fmain
	integer	::	ntrace_max
	integer	::	iflag_shot_info

	character(len=para_char_flen)	::	par_name

	open(33, file=fn_par, status='old', iostat=ierr)
	if(ierr /= 0)then
		write(*,*)	'Parameter file cannot open right, please check it !!!'
		write(*,*)	'Shutting down the program'
		stop
	endif

    !*初始的速度模型
	read(33, '(a)')	par_name
    read(33, '(a)') fn_vel

	!*梯度结果
	read(33, '(a)')	par_name
    read(33, '(a)') fn_grad

    !*fn1	the shot gather filename
	read(33, '(a)')	par_name
    read(33, '(a)') fn_cs

	!*iflag_shot_info
	read(33, '(a)')	par_name
	read(33, *) iflag_shot_info

	!*有关炮检点坐标信息的文件
	read(33, '(a)')	par_name
	read(33, '(a)') fn_shot_info

	!*	the dimension of the velocity model
	read(33, '(a)')	par_name
    read(33, *)	nvx, nvz

	!*	the sampling rate of the velocity model(m)
	read(33, '(a)')	par_name
    read(33, *)	dvx, dvz

	!*	the first, interval and last shot
	read(33, '(a)')	par_name
    read(33, *)	ns_start, ns_interval, ns_end

	!*	the first, interval and last shot
	read(33, '(a)')	par_name
    read(33, *)	ntrace_max

	!*	the imaging grid interval
	read(33, '(a)')	par_name
    read(33, *)	dx

	!*	the migration aperture
	read(33, '(a)')	par_name
    read(33, *)	aperture_x_l, aperture_x_r, aperture_z_u, aperture_z_d

	!*	'the trace length (sample point number) and dt(ms)
	read(33, '(a)')	par_name
    read(33, *)	lt, dt, lt_real, dt_real

	!*	the imaging depth(m) and step_length(m)
	read(33, '(a)')	par_name
    read(33, *)	depth, dz

	!*	the velocity starting coordinate (x)(m)
	read(33, '(a)')	par_name
    read(33, *)	cvx_initial, cvz_initial

	!*	the width of absorbing boundary of x axis
	read(33, '(a)')	par_name
    read(33, *)	boundary_x_l, boundary_x_r, boundary_z_u, boundary_z_d

	!*	the main frequency(hz)
	read(33, '(a)')	par_name
    read(33, *)	fmain

    close(33)

    return
    end	subroutine

!***********************************************************************

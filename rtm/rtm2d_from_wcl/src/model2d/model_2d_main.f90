!*=======================================================================*
!*                                                                       *
!*                            2D  model                                  * 
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
	program	model_2D

	use global
	use omp_lib
	implicit none
	include 'mpif.h'

	!Data dictionary: declare variable types, definitons, & units
	integer     ::  para_int(para_int_fnum)
	integer*8   ::  para_long(para_long_fnum)
	real        ::  para_float(para_float_fnum)
	real*8      ::  para_double(para_double_fnum)
	character(len=para_char_flen)	::	para_char(para_char_fnum)


	character(len=para_char_flen)	::	fn_par
	!mpi parameters
	integer	::	ierr, mpi_np, myid
	integer	::	threads
	integer ::	master
	character(MPI_MAX_PROCESSOR_NAME)   ::  processorname
	integer ::  processorlen
	integer ::  processorpid

	
	! variables from the parameter file.
	character(len=para_char_flen)	::	fn_cs, fn_vel
	character(len=para_char_flen)	::	fn_shot_info
	character(len=para_char_flen)  ::  fn_dynamic
	integer	::	flag_shot_info
	integer	::	flag_output_sig
	integer	::	nvx, nvz
	real	::	dvx, dvz
	integer	::	ns_total
	integer	::	ntrace_max

	integer	::	lt
	real	::	dt
	real	::	dx, dz
	real	::	aperture_x_l, aperture_x_r
	real	::	aperture_z_u, aperture_z_d
	real	::	boundary_x_l, boundary_x_r
	real	::	boundary_z_u, boundary_z_d
	real	::	cvx_initial, cvz_initial
	real	::	fmain

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
	CALL MPI_BCAST(FN_PAR, PARA_CHAR_FLEN, MPI_CHARACTER, MASTER, MPI_COMM_WORLD, IERR)
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
    CALL ReadPar(fn_par, &
			fn_cs, fn_vel, fn_dynamic, &
			fn_shot_info, &
			flag_shot_info, &
			flag_output_sig, &
			nvx, nvz, & 
			dvx, dvz, &
			ns_total, ntrace_max, &
			lt, &
			dt, &
			dx, dz, &
			aperture_x_l, aperture_x_r, &
			aperture_z_u, aperture_z_d, &
			boundary_x_l, boundary_x_r, &
			boundary_z_u, boundary_z_d, &
			cvx_initial, cvz_initial, &
			fmain )

	para_char(1)	=	fn_vel
	para_char(2)	=	fn_cs
	para_char(3)	=	fn_dynamic
	para_char(4)	=	fn_shot_info

	para_int(1)		=	nvx	
	para_int(2)		=	nvz
	para_int(3)		=	lt		
	para_int(4)		=	ns_total	
	para_int(5)		=	ntrace_max
	para_int(6)		=	flag_shot_info
	para_int(7)		=	flag_output_sig


	para_float(1)	=	fmain		
	para_float(2)	=	dt		
	para_float(3)	=	aperture_x_l	
	para_float(4)	=	aperture_x_r
	para_float(5)	=	aperture_z_u	
	para_float(6)	=	aperture_z_d
	para_float(7)	=	boundary_x_l
	para_float(8)	=	boundary_x_r
	para_float(9)	=	boundary_z_u
	para_float(10)	=	boundary_z_d
	para_float(11)	=	dx		
	para_float(12)	=	dz	
	para_float(13)	=	dvx			
	para_float(14)	=	dvz			
	para_float(15)	=	cvx_initial
	para_float(16)	=	cvz_initial

	!********************************************************************
	!===================== Begin the main program =======================
    if( myid .eq. master ) then
		write(*,*) '********************************************'
		write(*,*) '         2D acoustic wave modelling         '
        write(*,*) '                 begins                     '
        write(*,*) '********************************************'
    endif

    CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)

    CALL model_2d_mpi(para_int, para_long, para_float, para_double, para_char, &
				mpi_np, myid)

	!*===============  Closing the MPI working environment ==============
    CALL MPI_FINALIZE(IERR)

	!=================== End the main program ===========================
    if( myid .eq. master ) then
		write(*,*) '********************************************'
		write(*,*) '         2D acoustic wave modelling         '
		write(*,*) '                 end                        '
        write(*,*) '********************************************'
    endif

	!********************************************************************
	!*                   The Main Program end !                         *
	!********************************************************************
	stop
    end	program
!***********************************************************************

    subroutine ReadPar(fn_par, &
			fn_cs, fn_vel, fn_dynamic, &
			fn_shot_info, &
			flag_shot_info, &
			flag_output_sig, &
			nvx, nvz, & 
			dvx, dvz, &
			ns_total, ntrace_max, &
			lt, &
			dt, &
			dx, dz, &
			aperture_x_l, aperture_x_r, &
			aperture_z_u, aperture_z_d, &
			boundary_x_l, boundary_x_r, &
			boundary_z_u, boundary_z_d, &
			cvx_initial, cvz_initial, &
			fmain )

	use global
	implicit none
	!Data dictionary: declare variable types, definitons, & units
	!Dummy Variables
	character(len=para_char_flen)	::	fn_par
	character(len=para_char_flen)	::	fn_cs, fn_vel
	character(len=para_char_flen)   ::	fn_dynamic
	character(len=para_char_flen)   ::	fn_shot_info
	integer	::	flag_shot_info
	integer	::	flag_output_sig
	integer ::	ierr
	integer ::	nvx, nvz
	real	::	dvx, dvz
	integer	::	ns_total
	integer	::	ntrace_max
	integer ::	lt
	real	::	dt
	real	::	dz
	real	::	dx
	real	::	aperture_x_l, aperture_x_r
	real	::	aperture_z_u, aperture_z_d
	real	::	boundary_x_l, boundary_x_r
	real	::	boundary_z_u, boundary_z_d
	real	::	cvx_initial, cvz_initial
	real	::	fmain

	character(len=para_char_flen)	::	par_name

	open(33, file=fn_par, status='old', iostat=ierr)
	if(ierr /= 0)then
		write(*,*)	'Parameter file cannot open right, please check it !!!'
		write(*,*)	'Shutting down the program'
		stop
	endif

    !*fn2	the velocity model filename'
	read(33, '(a)')	par_name
    read(33, '(a)') fn_vel

    !*fn1	the shot gather filename
	read(33, '(a)')	par_name
    read(33, '(a)') fn_cs

	!*	flag_shot_info
	read(33, '(a)')	par_name
    read(33, *)	flag_shot_info

	!*fn3: the image profile
	read(33, '(a)')	par_name
    read(33, '(a)') fn_shot_info

	!*	flag_output_sig
	read(33, '(a)')	par_name
    read(33, *)	flag_output_sig

	!*fn_dynamic	the dynamic temp file for  storage gradient result
	read(33, '(a)')	par_name
	read(33, '(a)') fn_dynamic

	!*	'the trace length (sample point number) and dt(ms)
	read(33, '(a)')	par_name
    read(33, *)	lt, dt

	!*	the first, interval and last shot
	read(33, '(a)')	par_name
    read(33, *)	ns_total,	ntrace_max

	!*	the dimension of the velocity model
	read(33, '(a)')	par_name
    read(33, *)	nvx, nvz

	!*	the sampling rate of the velocity model(m)
	read(33, '(a)')	par_name
    read(33, *)	dvx, dvz

	!*	the main frequency(hz)
	read(33, '(a)')	par_name
    read(33, *)	fmain

	!*	the imaging grid interval
	read(33, '(a)')	par_name
    read(33, *)	dx, dz

	!*	the velocity starting coordinate (x)(m)
	read(33, '(a)')	par_name
    read(33, *)	cvx_initial, cvz_initial

	!*	the migration aperture
	read(33, '(a)')	par_name
    read(33, *)	aperture_x_l, aperture_x_r, aperture_z_u, aperture_z_d

	!*	the width of absorbing boundary of x axis
	read(33, '(a)')	par_name
    read(33, *)	boundary_x_l, boundary_x_r, boundary_z_u, boundary_z_d


    close(33)

    return
    end	subroutine

!***********************************************************************

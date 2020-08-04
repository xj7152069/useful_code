!*=======================================================================*
!*                                                                       *
!*             2D REVERSE TIME PRESTACK DEPTH MIGRATION	      			 *
!*                                                                       *
!*=======================================================================*
!*			THE PROGRAM IS RUNNING ON 2D SINGLE SHOT GATHER              * 
!*=======================================================================*
!*                                                                       *
!*                      ALL RIGHTS RESERVED                              *
!*           WAVE PHENOMENA AND INVERSION IMAGING GROUP (WPI)            *
!*      BY SCHOOL OF OCEAN AND EARTH SCIENCE TONGJI UNIVERSITY           *
!*                                                                       *
!*=======================================================================*
!*                                                                       *
!*           Authers  : Wu Chengliang				                     *
!*           Date     : 2016, 05-29                                      *
!*           Modified : Wu Chengliang                                    *
!*           Data     : 2016, 05-29                                      *
!*           End data : 2016, 06-01                                      *
!*			 Version  :	NO.0										     *
!*																		 *
!*=======================================================================*
!*版本说明
!*2016-05-29: NO.0
!*	2维逆时偏移程序:
!*采用时间2阶、空间10阶的有限差分计算，边界条件采用PML吸收边界条件.
!*震源波场反传采用波场重构的方式实现
!*成像条件：采用互相关成像条件
!*2017-05-19: NO.1
!*	测试波场反传策略对重构波场的影响
!*边界存储的点数为差分的节点数
!*
!*************************************************************************
	include '../include/global_rtm.f90'
	include '../include/fft_fortran.f90'

!*************************************************************************
!*						The main program begin	                         *
!*************************************************************************
	program	 Reverse_Time_Migration_Acousitc_2D_PML

	implicit none
	include 'mpif.h'
!*======================================================================
!*  Variable demonstrations:
!*
!*		fn1: the shot gather file name
!*		fn2: the velocity file name
!*		fn3: the image profile
!*		fn_single: the single gradient profile file
!*		fn_dynamic:the dynamic file storage image result 
!*
!*		ns_start:	the start shot number for 2d rtm
!*		ns_end:		the final shot number for 2d rtm
!*		ns_interval:the step length of shot number for 2d reverse time shot migration
!*    
!*		lt:			the extrapolation sample point number of the trace length
!*		lt_real:	the really sample point number of the trace length
!*		dt:			the extrapolation sample rate(ms)
!*		dt_real:	the really sample rate(ms)

!*		depth:		the exptropalation/imaging depth(m) 
!*		dz:			the exptropalation step length(m)
!*		dx:			the interval of imaging grid(m)

!*		aperture_x_l:the migration aperture of left of x_direction(m)
!*		aperture_x_r:the migration aperture of right of x_direction(m)
!*
!*		boundary_x_l:the width of the absorbing boundary zone along left of x_direction(m)
!*		boundary_x_r:the width of the absorbing boundary zone along right of x_direction(m)
!*		boundary_z_u:the width of the absorbing boundary zone along up of z_direction(m)
!*		boundary_z_d:the width of the absorbing boundary zone along dowm of z_direction(m)

!*		cvx_initial:the initial coodinate of the velocity model
!*		nvx, nvz:	the dimension of the velocity model
!*		dvx, dvz:	the sample rate of the velocity model    

!*		fmain:		the main frequency of the source function 
!*
!*======================================================================

	!Data dictionary: declare variable types, definitons, & units
	integer	::	ierr, mpi_np, myid
	integer ::	master

	character(len=256)	::	fn_par, fn1, fn2, fn3
	character(len=256)  ::  fn_single
	character(len=128)  ::  fn_dynamic
	integer	::	ns_start, ns_interval, ns_end
	integer	::	lt,	lt_real
	real	::	dt, dt_real
	real	::	depth
	real	::	dx, dz
	real	::	aperture_x_l, aperture_x_r
	real	::	boundary_x_l, boundary_x_r
	real	::	boundary_z_u, boundary_z_d
	real	::	cvx_initial
	integer	::	nvx, nvz
	real	::	dvx, dvz
	real	::	fmain
	real	::	max_offset_x
	real	::  min_offset_x
	integer ::  shot_direct
	integer	::	order

	integer ::	ns_total
	integer ::	nz

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
	CALL MPI_BCAST(FN_PAR, 256, MPI_CHARACTER, MASTER, MPI_COMM_WORLD, IERR)

    CALL ReadPar(fn_par, fn1, fn2, fn3, &
			fn_single, fn_dynamic, &
			ns_start, ns_interval, ns_end, &
			lt, lt_real, &
			dt, dt_real, &
			depth, &
			dx, dz, &
			aperture_x_l, aperture_x_r, &
			boundary_x_l, boundary_x_r, &
			boundary_z_u, boundary_z_d, &
			Cvx_initial, &
			nvx, nvz, & 
			dvx, dvz, &
			fmain, &
			max_offset_x, min_offset_x, &
			shot_direct, order)

	!********************************************************************
	if(myid .eq. master) then
		write(*,*)	'the shot gather filename',					trim(fn1)
		write(*,*)	'the velocity model filename',				trim(fn2)
		write(*,*)	'the migrated and stacked profile',			trim(fn3)
		write(*,*)  'the single gradient profile file',         trim(fn_single)
		write(*,*)  'the dynamic file storage image result',	trim(fn_dynamic)
		write(*,*)	'the first, interval and last shot',		ns_start, ns_interval, ns_end
		write(*,*)	'the trace length (sample point number)',	lt, lt_real
		write(*,*)	'the trace length (sample point number)',	dt, dt_real
		write(*,*)	'the imaging depth(m)',						depth
		write(*,*)	'the imaging grid interval',				dx, dz
		write(*,*)	'the migration aperture',					aperture_x_l, aperture_x_r
		write(*,*)	'the width of absorbing boundary',			boundary_x_l, boundary_x_r
		write(*,*)	'the width of absorbing boundary',			boundary_z_u, boundary_z_d
		write(*,*)	'the velocity starting coordinate (x)(m)',	Cvx_initial
		write(*,*)	'the dimension of the velocity model',		nvx, nvz
		write(*,*)	'the sampling rate of the velocity model',	dvx, dvz
		write(*,*)	'the main frequency(hz)',					fmain
		write(*,*)  'the maximum offset of x direction',        max_offset_x
		write(*,*)  'the minimum offset of x direction',        min_offset_x
		write(*,*)  'judge direct wave exist',                  shot_direct
		write(*,*) "                                       "
		write(*,*) "********* Check parameters ! **********"
		write(*,*) "                                       "
	endif

	!===================== Calulate The Image Depth ======================
    nz=nint(depth/dz)
    if(nz.gt.nvz) then
		write(*,*) 'the extrapolation depth is too large and '
		write(*,*) 'is beyond the depth of the defined velocity'
		write(*,*) 'field, please modify it'
		stop
    end if

	!==================== Calulate The Image Shot Num ====================
    ns_total=(ns_end-ns_start)/ns_interval+1

	!********************************************************************
	!	decide whether the parameters  are right

	!********************************************************************
	!===================== Begin the main program =======================
    if( myid .eq. master ) then
		write(*,*) '********************************************'
		write(*,*) '     Shot gather reverse time migrating     '
        write(*,*) '                 begins                     '
        write(*,*) '********************************************'
    endif
    
    CALL Multiple_Shot_RTM(fn1, fn2, fn3, &
			fn_single, fn_dynamic, &
			ns_start, ns_interval, ns_end, ns_total, &
			lt, lt_real, &
			dt, dt_real, &
			depth, nz, &
			dx, dz, &
			aperture_x_l, aperture_x_r, &
			boundary_x_l, boundary_x_r, &
			boundary_z_u, boundary_z_d, &
			Cvx_initial, &
			nvx, nvz, & 
			dvx, dvz, &
			fmain, master, &
			max_offset_x, min_offset_x, &
			shot_direct, order, &
			myid, mpi_np)

	!*===============  Closing the MPI working environment ==============
    CALL MPI_FINALIZE(IERR)

	!=================== End the main program ===========================
    if( myid .eq. master ) then
        write(*,*) '********************************************'
        write(*,*) '                                            '
        write(*,*) '     All shots migrating finished           '
        write(*,*) '                                            '
        write(*,*) '********************************************'
    endif

	!********************************************************************
	!*                   The Main Program end !                         *
	!********************************************************************
	stop
    end	program

!***********************************************************************
    subroutine Multiple_Shot_RTM(fn1, fn2, fn3, &
			fn_single, fn_dynamic, &
			ns_start, ns_interval, ns_end, ns_total, &
			lt, lt_real, &
			dt, dt_real, &
			depth, nz, &
			dx, dz, &
			aperture_x_l, aperture_x_r, &
			boundary_x_l, boundary_x_r, &
			boundary_z_u, boundary_z_d, &
			Cvx_initial, &
			nvx, nvz, & 
			dvx, dvz, &
			fmain, master, &
			max_offset_x, min_offset_x, &
			shot_direct, order, &
			myid, mpi_np)

	use	global
	use header_module
	implicit none
	include 'mpif.h'

	!Data dictionary: declare variable types, definitons, & units
	!Dummy variables
    integer status(mpi_status_size)
	character(len=256)	::	fn1, fn2, fn3
	character(len=256)  ::  fn_single
	character(len=128)  ::  fn_dynamic
    integer	::	ns_start, ns_interval, ns_end, ns_total
	integer	::	lt, lt_real
	real	::	dt, dt_real
	real	::	depth
	integer	::	nz
	real	::	dx, dz
	real	::	aperture_x_l, aperture_x_r
	real	::	boundary_x_l, boundary_x_r
	real	::	boundary_z_u, boundary_z_d
	real	::	Cvx_initial
	integer	::	nvx, nvz
	real	::	dvx, dvz
	real	::	fmain
	integer ::	master
	real	::	max_offset_x
	real	::	min_offset_x
	integer	::	order
	integer ::	shot_direct
	integer	::	myid, mpi_np

	!Local variables
    character(len=256)	::	fn_dynamic_all
	integer	::	i
	integer ::	isfirst
	integer ::	isend
	integer ::	ix, iz
	integer ::	ierr
	integer ::	iter
	integer ::	iss
	integer ::	ntrace
	integer ::	krec_save
	integer	::	xs
	integer	::	cx_min
	integer	::	cx_max
	integer ::	nx_with_apert
	integer ::	nz_with_apert
	integer ::	nnx, nnz
	integer ::	nvxx_shift
	integer	::	nshot
	integer ::	nx
	integer ::	nxz
	integer ::	nx_apert_l, nx_apert_r
	integer ::	nx_bound_l, nx_bound_r
	integer ::	nz_bound_u, nz_bound_d
	integer ::	iflag

	real	::	dt2
	real	::	cvx_final


    integer,allocatable	::  index_table(:,:)
	real,allocatable	::	buf(:)
	

	integer	::	r_table(9)
    integer*8	::	l1, l2, l3, l4, l5, l6, l7, l8, l9, l10
	integer*8	::	l11, l12, l13, l14, l15, l16, l17, l18, l19, l20
	integer*8	::	l21, l22, l23, l24, l25, l26, l27, l28, l29, l30
	integer*8	::	l31, l32, l33, l34
	integer*8	::	l_total
	integer*8	::	lmax_remain

!*=========================================================================
!*  Variable demonstrations:
!*
!*		cx_min, cx_max: the maximum and minimum coordinates
!*                          of the current shot gather
!*		cvx_final:	the maximum coordinates of the velocity model
!* 
!*		nx_apert: the sample point number of the added migration aperture
!*		nx_bound, nz_bound: the sample point number of the added 
!*                                    absorbing boundary range
!*
!*		nx: the sample point number of the shot and receiver defined range
!*		nx_witg_apert: the sample point number after addiing migration 
!*                                   aperture 
!*		nnx, nnz: the sample point number of the wavefield extrapolation range
!*=========================================================================

	!=========== allocate velocity and stack index working buff =========
	allocate(index_table(9,ns_total),stat=ierr)
	if(ierr.ne.0)then
		write(*,*)"Can not allocate working memory about index_table, stop!!!"
		stop
	endif
	index_table=0.0

	!*=====================================================================*
	iflag=1		!the flag of controlling the print  on-off

	!*change the dt ms to s 
	dt = dt*0.001
	dt_real = dt_real*0.001
	dt2=dt*dt

	!*comoute some parameters for calculation
    cvx_final =  cvx_initial +(nvx-1)*dvx

    nx_apert_l = aperture_x_l/dx +0.5
    nx_apert_r = aperture_x_r/dx +0.5

    nx_bound_l = boundary_x_l/dx +0.5
    nx_bound_r = boundary_x_r/dx +0.5

    nz_bound_u = boundary_z_u/dz +0.5
    nz_bound_d = boundary_z_d/dz +0.5

	if(iflag .eq. 1)then
		if(myid .eq.master)then
			write(*,*)'dt=',dt, 'dt_real=',dt_real,'dt2=',dt2
			write(*,*)'nx_apert_l=',nx_apert_l,'nx_apert_r=',nx_apert_r
			write(*,*)'nx_bound_l=',nx_bound_l,'nx_bound_r=',nx_bound_r
			write(*,*)'master=',master
		endif
	endif

	!========= The Master Node Open The Seismic Data And Do The Index ========
	!*打开炮集文件
	open(11, file=fn1, access='direct', recl=lbyte*(60+lt_real), status='old', iostat=ierr)
	if(ierr /= 0)then
		write(*,*)	'Shot gather file cannot open right, please check it !!!'
		write(*,*)	'Shutting down the program'
		stop
	endif

	!打开速度文件
	open(12, file=fn2, access='direct', status='old', recl=lbyte*nvz, iostat=ierr)
	if(ierr /= 0)then
		write(*,*)  'Velocity file cannot open right, please check it!!!'
		write(*,*)  'Velocity file ',fn2
		stop
	endif

	!========================================================================
    if( myid .eq. master ) then
		print*, 'Form assignment card started...'
        CALL Index_Data(ns_start, ns_total, ns_interval,&
              cx_min, cx_max, nshot, index_table, lt_real, &
			  max_offset_x, min_offset_x)
        print*,'nshot/ns=',nshot,'/',ns_total
		write(*,*)'Index end! wait for MPI_BCAST'
    else
		!注意fn_dynamic文件路径是否正确(包括节点是否有储存能力、是否有足够储存空间等)
		CALL Form_Dynamic_Filename(myid, fn_dynamic_all, fn_dynamic, ns_start)   
		open(14, file=trim(adjustl(fn_dynamic_all)), access='direct', status='replace', &
			recl=lbyte*nz, iostat=ierr)
		if(ierr /= 0)then
			write(*,*)  'Dynamic file cannot open right, please check it!!!'
			write(*,*)  'fn_dynamic_all=',fn_dynamic_all
			stop
		endif
    endif

	!========================= MPI Enviroment Initialed ======================
    CALL MPI_BCAST(NSHOT,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
    CALL MPI_BCAST(NX,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
    CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)

	!*=======================================================================*
	isfirst=1

    if( myid .eq. master ) then
		do i=1, nshot
			CALL MPI_RECV(R_TABLE, 9, MPI_INTEGER, MPI_ANY_SOURCE,&
					MPI_ANY_TAG, MPI_COMM_WORLD, STATUS, IERR)
			isend=status(mpi_source)
			CALL MPI_SEND(INDEX_TABLE(1,I), 9, MPI_INTEGER, ISEND, I,&
					MPI_COMM_WORLD, IERR)
        enddo
        do i = 1, mpi_np-1
			CALL MPI_RECV(R_TABLE, 9, MPI_INTEGER, MPI_ANY_SOURCE,&
					MPI_ANY_TAG, MPI_COMM_WORLD, STATUS, IERR)
			isend=status(mpi_source)
			r_table(1) = 0
			r_table(2) = 0
			r_table(3) = 0
			r_table(4) = 0
			r_table(5) = 0
			r_table(6) = 0
			r_table(7) = 0
			r_table(8) = 0
			r_table(9) = 0
			CALL MPI_SEND(R_TABLE, 9, MPI_INTEGER, ISEND, I,&
					MPI_COMM_WORLD, IERR)
        enddo
    else
		!*--------------------worker-process------------------------------------*
		r_table(1) = 0
        r_table(2) = 0
        r_table(3) = 0
        r_table(4) = 0
        r_table(5) = 0
        r_table(6) = 0
        r_table(7) = 0
        r_table(8) = 0
        r_table(9) = 0
        CALL MPI_SEND(R_TABLE, 9, MPI_INTEGER,0, 0, MPI_COMM_WORLD, IERR)

90      continue
        CALL MPI_RECV(INDEX_TABLE, 9, MPI_INTEGER, 0, MPI_ANY_TAG, MPI_COMM_WORLD, STATUS, IERR)

        if( index_table(1,1) .ne. 0 ) then
			iss      = index_table(1,1)
			ntrace   = index_table(2,1)
			xs       = index_table(3,1)

			cx_min   = index_table(5,1)
			cx_max   = index_table(6,1)


			krec_save= index_table(9,1)

			!*the imaging range after adding the migration aperture
			nx = (cx_max-cx_min)/dx + 1.5
			nx_with_apert = aperture_x_l/dx + nx + aperture_x_r/dx
			nz_with_apert = nz

			!*the wave extrapolation range after adding the absorbing boundary ---*
			nnx = boundary_x_l/dx +  nx_with_apert + boundary_x_r/dx
			nnz = boundary_z_u/dz +  nz_with_apert + boundary_z_d/dz

			!*the shift between the current shot imaging range and the total velocity field range
			nvxx_shift=((cx_min-aperture_x_l)-cvx_initial)/dvx+1

			if(iflag .eq. 1)then
				if(myid .eq. 1)then
					write(*,*)	'iss',	iss
					write(*,*)	'ntrace',	ntrace
					write(*,*)	'xs',	xs
					write(*,*)	'cx_min',	cx_min
					write(*,*)	'cx_max',	cx_max
					write(*,*)	'krec_save',	krec_save
					write(*,*)	'nx=',	nx
					write(*,*)	'nx_with_apert',	nx_with_apert
					write(*,*)	'nz_with_apert',	nz_with_apert
					write(*,*)	'nnx',	nnx
					write(*,*)	'nnz',	nnz
					write(*,*)	'nvxx_shift',	nvxx_shift
				endif
			endif

			!*assigning the inner memory volume
			l1  = 1										!
			l2  = l1 + lt								!wavelet
			l3  = l2 + lt								!trace
			l4  = l3 + lt_real							!trace_real
			l5  = l4 + nx*lt							!shot_gather
			l6  = l5 + nvz								!vv_tmp
			l7  = l6 + nnx*nnz							!vv
			l8  = l7 + nx_with_apert					!vv_surface
			l9  = l8 + nx_with_apert*nz					!image
			l10 = l9 + nvx*nz							!image_tmp

			l11 = l10 + (5+nnx+5)*(5+nnz+5)				!wavefieldforward1_wf1
			l12 = l11 + (5+nnx+5)*(5+nnz+5)				!wavefieldforward2_wf2
			l13 = l12 + (5+nnx+5)*(5+nnz+5)				!wavefieldreverse1_wr1
			l14 = l13 + (5+nnx+5)*(5+nnz+5)				!wavefieldreverse2_wr2

			!pml_for x direction
			l15 = l14 + nnz*(nx_bound_l+nx_bound_r)		!u1_x 
			l16 = l15 + nnz*(nx_bound_l+nx_bound_r)		!u1_x_last
			l17 = l16 + nnz*(nx_bound_l+nx_bound_r)		!u2_x
			l18 = l17 + nnz*(nx_bound_l+nx_bound_r)		!u2_x_last
			l19 = l18 + nnz*(nx_bound_l+nx_bound_r)		!u2_x_temp
			l20 = l19 + nnz*(nx_bound_l+nx_bound_r)		!u2_x_temp_last
			l21 = l20 + nnz*(nx_bound_l+nx_bound_r)		!u3_x
			l22 = l21 + nnz*(nx_bound_l+nx_bound_r)		!u3_x_last
	
			!pml_for z direction
			l23 = l22 + nnx*(nz_bound_u+nz_bound_d)		!u1_z
			l24 = l23 + nnx*(nz_bound_u+nz_bound_d)		!u1_z_last
			l25 = l24 + nnx*(nz_bound_u+nz_bound_d)		!u2_z
			l26 = l25 + nnx*(nz_bound_u+nz_bound_d)		!u2_z_last
			l27 = l26 + nnx*(nz_bound_u+nz_bound_d)		!u2_z_temp
			l28 = l27 + nnx*(nz_bound_u+nz_bound_d)		!u2_z_temp_last
			l29 = l28 + nnx*(nz_bound_u+nz_bound_d)		!u3_z
			l30 = l29 + nnx*(nz_bound_u+nz_bound_d)		!u3_z_last

			l31 = l30 + nx_with_apert*lt*order				!top
			l32 = l31 + nx_with_apert*lt*order				!bot
			l33 = l32 + nz_with_apert*lt*order				!lef
			l34 = l33 + nz_with_apert*lt*order				!rig
			!*=====================================================================*

			l_total=l34
			if(myid .eq. 1)then
				write(*,*) 'memory_needed=',l_total*4*1.0/1024/1024,'m' 
				write(*,*) 'wavefield_slice',(l11-l10)*1.0*4/1024/1024,'m' 
			endif

			!*checking whether the total inner memory application enough
			if(lmax_remain.lt.0) then
				write(*,*) 'lmax_remain=',lmax_remain
				write(*,*) 'the acquired inner memory is not enough'
				write(*,*) 'please enlarge it'
			!	stop
			end if

			allocate( buf(l_total),stat=ierr)
			if(ierr.ne.0.0)then
				write(*,*) 'Can not allocate single shot working memory ' 
				write(*,*)'Please check it and shut down the program......'
				stop
			endif

			buf=0.0

			!*=====================================================================*
			!* single shot reverse time extrapolation and getting image
			CALL Single_Shot_RTM_2D( &
					buf(l1), buf(l2), buf(l3), &
					buf(l4), &
					buf(l5), buf(l6), buf(l7), &
					buf(l8), buf(l9), &
					buf(l10), buf(l11), buf(l12), buf(l13), &
					buf(l14), buf(l15), &
					buf(l16), buf(l17), &
					buf(l18), buf(l19), &
					buf(l20), buf(l21), &
					buf(l22), buf(l23), &
					buf(l24), buf(l25), &
					buf(l26), buf(l27), &
					buf(l28), buf(l29), &
					buf(l30), buf(l31), buf(l32), buf(l33), &
					fn1, iss, krec_save, &
					nx, nz, nnx, nnz, &
					nx_with_apert, nz_with_apert, &
					lt, lt_real, &
					nx_bound_l, nz_bound_u, &
					nx_bound_r, nz_bound_d, &
					nx_apert_l, nx_apert_r, &
					nvx, nvz, nvxx_shift, &
					max_offset_x, min_offset_x, &
					shot_direct, cx_min, &
					dx, dz, dvx, dvz, &
					dt, dt2, dt_real, &
					fmain, cvx_initial, &
					order, myid)
          
			isfirst=0
			deallocate (buf)

			CALL MPI_SEND(INDEX_TABLE, 9, MPI_INTEGER, 0, 0, MPI_COMM_WORLD, IERR)
			goto 90
		endif
    endif
	!*=====================================================================*

	!归约结果
	CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )
	CALL Reduce_Result_Gradient(fn3, myid, master, nvx, nz)

	!删除临时文件
	if(myid .ne. master)then
!		close(14, status='delete',err=333)
		close(14)
333		continue
	endif

	deallocate(index_table)

	close(11)
	close(12)

	return
	end subroutine

!***********************************************************************
    subroutine Single_Shot_RTM_2D( &
					wavelet, trace, trace_real, &
					shot_gather, &
					vv_tmp, vv, vv_surface, &
					image, image_tmp, &
					wf1, wf2, wr1, wr2, &
					u1_x, u1_x_last, &
					u2_x, u2_x_last, &
					u2_x_temp, u2_x_temp_last, &
					u3_x, u3_x_last, &
					u1_z, u1_z_last, &
					u2_z, u2_z_last, &
					u2_z_temp, u2_z_temp_last, &
					u3_z, u3_z_last, &
					top, bot, lef, rig, &
					fn1, iss, krec_save, &
					nx, nz, nnx, nnz, &
					nx_with_apert, nz_with_apert, &
					lt, lt_real, &
					nx_bound_l, nz_bound_u, &
					nx_bound_r, nz_bound_d, &
					nx_apert_l, nx_apert_r, &
					nvx, nvz, nvxx_shift, &
					max_offset_x, min_offset_x, &
					shot_direct, cx_min, &
					dx, dz, dvx, dvz, &
					dt, dt2, dt_real, &
					fmain, cvx_initial, &
					order, myid)

!*=====================================================================*
!*	Variable demonstration:
!*   
!*		wavelet:		the working buffer for the source function
!*		trace:			the working buffer for the shot gather inputting
!*		trace_real:		the working buffer for the actually single shot gather 
!*		shot_gather:	the working buffer for the current shot gather
!*
!*		vv:				the working buffer for the current shot velocty field
!*		vv_surface:		the working buffer for the current surface velocty field
!*
!*		wf1, wf2 :		the working buffer for forward time extrapoaltion
!*		wr1, wr2 :		the working buffer for reverse time extrapolation
!*
!*		fn1:			the input shot gather file
!*		iss:			the shot number
!*		krec_save:		the current shot beginning position in fn1
!*		vv_all:			the working buffer for inputting the one_layer velocity field
!*
!*		nx, nz:			the dimension of the current shot gather
!*		nnx, nnz:		the wave filed extrapolation range after adding 
!*					the random boundary 
!*		nx_with_apert: the dimension of the current shot gather after
!*                  adding the migration aperture, they define 
!*                  the imaging range of the current shot gather
!*
!*		lt:				the extrapolation sample point number
!*		lt_real:		the actually sample point number
!*		nx_bound, nz_bound : the sample point number of the added 
!*					absorbing boundary range
!*		nx_apert, ny_apert: the sample point number of the added migration aperture
!*		nvx, nvy:		the dimension of the total velocity field
!*		nvxx_shift:		the shift between the current shot imaging
!*					range and the total velocity field range
!*
!*		max_offset_x:	the maxinum of offset of x direction
!*		min_offset_x:	the mininun of offset of x direction
!*
!*		shot_direct:	the flag decide whether direct wave exist
!*		cx_min:			the mininum  coordinate of current shot gather of x direction
!*
!*		dx, dz:			the imaging grid interval(m)
!*		dvx, dvz:		the sampling rate of the velocity model(m)
!*		dt:				the extrapolation sampling rate of time (s)
!*		dt_real:		the actually sampling rate of time (s)
!*		dt2:			dt*dt

!*		fmain:			the main frequency(hz)
!*		cvx_initial:	the velocity starting coordinate (x)(m)
!*
!*		Local variables:
!*		ns_x, ns_y:		the current shot position in the current velocity  
!*					field and the current extrapolation range 
!*
!*=====================================================================*

	use global
	implicit none
	!Data dictionary: declare variable types, definitons, & units
	!Dummy Variables
	character(len=256)	::	fn1
	integer ::	myid

	integer ::	iss
	integer	::	krec_save
    integer ::	nx, nz
    integer ::	nnx, nnz
    integer ::	nx_with_apert, nz_with_apert
    integer ::	lt
	integer ::	lt_real
    integer ::	nx_bound_l, nz_bound_u
    integer ::	nx_bound_r, nz_bound_d
    integer ::	nx_apert_l, nx_apert_r
	integer ::	nvx, nvz
	integer ::	nvxx_shift
	real	::  max_offset_x
	real	::  min_offset_x
	integer ::  shot_direct
	integer	::	cx_min

	real	::	dx, dz
	real	::	dvx, dvz
	real	::	dt
	real	::	dt2
	real	::	dt_real
	real	::	fmain
	real	::	cvx_initial

	integer	::	order

	!declare buf variables for single shot migration
	real	::	wavelet(lt)
	real	::	trace(lt)
	real	::	trace_real(lt_real)
	real	::	shot_gather(lt,nx)
	real	::	vv_tmp(nvz)
	real	::	vv(nnz, nnx)
	real	::	vv_surface(nx_with_apert)
	real	::	image(nz, nx_with_apert)					
	real	::	image_tmp(nz, nvx)	

    real	::	wf1(-4:nnz+5, -4:nnx+5)
    real	::	wf2(-4:nnz+5, -4:nnx+5)
    real	::	wr1(-4:nnz+5, -4:nnx+5)
    real	::	wr2(-4:nnz+5, -4:nnx+5)

    !wavefield for pml field
    !forward
    !x direction
    real	::	u1_x(nnz,nx_bound_l+nx_bound_r)
    real	::	u1_x_last(nnz,nx_bound_l+nx_bound_r)
    real	::	u2_x(nnz,nx_bound_l+nx_bound_r)
    real	::	u2_x_last(nnz,nx_bound_l+nx_bound_r)
    real	::	u2_x_temp(nnz,nx_bound_l+nx_bound_r)
    real	::	u2_x_temp_last(nnz,nx_bound_l+nx_bound_r)
    real	::	u3_x(nnz,nx_bound_l+nx_bound_r)
    real	::	u3_x_last(nnz,nx_bound_l+nx_bound_r)

    !z direction
    real	::	u1_z(nz_bound_u+nz_bound_d,nnx)
    real	::	u1_z_last(nz_bound_u+nz_bound_d,nnx)
    real	::	u2_z(nz_bound_u+nz_bound_d,nnx)
    real	::	u2_z_last(nz_bound_u+nz_bound_d,nnx)
    real	::	u2_z_temp(nz_bound_u+nz_bound_d,nnx)
    real	::	u2_z_temp_last(nz_bound_u+nz_bound_d,nnx)
    real	::	u3_z(nz_bound_u+nz_bound_d,nnx)
    real	::	u3_z_last(nz_bound_u+nz_bound_d,nnx)
	
    real	::	top(order, nx_with_apert, lt)	!v1
    real	::	bot(order, nx_with_apert, lt)	!v1
    real	::	lef(order, nz_with_apert, lt)	!v1
    real	::	rig(order, nz_with_apert, lt)	!v1

	!Local Variables
	integer	::	ix, iz
	integer ::	it, iit
	integer ::	iix, inx
	integer ::	iiz, inz
	integer ::	ixx
	integer ::	nx_shift
    integer ::	ns_x
	integer ::	ns_z
	integer ::	nr_z
	integer ::	istep
	integer ::	nwt
	integer ::	ierr
	integer ::	iflag
	integer	::	iii

	real	::	r			!attenuation coefficient of pml field
	real	::	coe(13)

	real	::	v1
	real	::	v2
	real	::	tmp

	!character(len=256)  ::  currt, currtsnap
	!character(len=256)  ::  snapx_fn= &
	!	'/data3/wcl/rtm/2d/snap1/snap_x_v1_1_'

	!*======================================================================= 
	iflag=1		!the flag of controlling the print  on-off
      
    !Initiallization  
    r=1.0e-3

	!*----------- zero the work buffer -------------------------------------*
    wavelet=0.0
    trace=0.0
	trace_real=0.0
    shot_gather=0.0
    vv=0.0
    vv_surface=0.0

    u1_x=0.0
    u1_x_last=0.0
    u2_x=0.0
    u2_x_last=0.0
    u2_x_temp=0.0
    u2_x_temp_last=0.0
    u3_x=0.0
    u3_x_last=0.0

    u1_z=0.0
    u1_z_last=0.0
    u2_z=0.0
    u2_z_last=0.0
    u2_z_temp=0.0
    u2_z_temp_last=0.0
    u3_z=0.0
    u3_z_last=0.0

    wf1=0.0
    wf2=0.0
    wr1=0.0
    wr2=0.0

	!*-------- input the current shot gather data into the working bufer====*
	CALL Input_Shot_Gather(trace, trace_real, shot_gather, &
				nx, dx, lt, dt, lt_real, dt_real, &
				ns_x, cx_min, krec_save, &
				max_offset_x, min_offset_x)

	if(iflag .eq. 1)then
		print*, 'nx, dx,  lt, dt, lt_real, dt_real,ns_x, cx_min, krec_save'
		print*, nx, dx,  lt, dt, lt_real, dt_real,ns_x, cx_min, krec_save
	endif

	!*======= forming the source function for the forward extrpolation ------*
    CALL Wavelet_Forming(wavelet, dt, lt, fmain, nwt)

	!*------ input the velocity field for the current shot extrapolation-----*
    CALL Get_Current_Shot_Velocity(vv_tmp, vv, vv_surface, &
				nvx, nvz, dvx, dvz,&
                nx_with_apert, nnx, nnz,&
                nx_bound_l, nz_bound_u,&
                nx_bound_r, nz_bound_d,&
                nx, nz,&
                nvxx_shift)

	!*-------------------------get coefficient ----------------------------*!
    CALL Form_Coe(coe, dx, dz, dt)

	!********** set up the scoure  parameter *******************************
    nx_shift = nx_apert_l + nx_bound_l

	ns_z=1
	nr_z=1
	

	!**********************************************************************!
	!cut direct wave
	!若shot_direct=1,则对原始炮集中直达波进行切除.
	if(shot_direct .eq. 1 )then
		CALL Cut_Mute_Direct_Wave(shot_gather, vv_surface, &
					ns_x, nx, nx_with_apert, &
					nz, lt, nx_apert_l, nwt, &
					max_offset_x, dx, dt )
	endif

	!*======================================================================= 
	!*--------- extrapolating source wavefield forward time ----------------*      
	!*=======================================================================
    istep=-1
    do it=1, lt+nwt				!!!! loop for forward source wavefield

		if(mod(it, 1000).eq. 1) then
			write(*,*)'shot=',iss,'myid=',myid,'  step=forward  ','it=',it
        endif

		!*======================================================================= 
		!储存波场快照，用于检查是否正确  lyh
		!if(modulo(it+500-nwt, 500).eq.1 .and. myid .eq. 999999 ) then
			!write(currt,'(I4)')it-nwt

			!write(currtsnap,*)trim(adjustl(snapx_fn)),trim(adjustl(currt))
			!open(133, file=trim(adjustl(currtsnap))//'.dat', access='direct',&
				!form='unformatted', status='replace', recl=lbyte*nx_with_apert*nz_with_apert, iostat=ierr)
			!write(133,rec=1)((wf2(iz, ix),iz=1+nz_bound_u,nz_with_apert+nz_bound_u), &
							!ix=1+nx_bound_l,nx_with_apert+nx_bound_l)
			!close(133)

		!endif  lyh

		!*======================================================================= 
        istep=istep+1
        if(mod(istep,2).eq.0)then	!*single time step extrapolation

			!*assignning the source function
			if(it.lt.nwt*2+100)then
				wf2(ns_z+nz_bound_u, ns_x+nx_shift) &     
					=wf2(ns_z+nz_bound_u, ns_x+nx_shift)+ wavelet(it)
			endif

			CALL  Extrapolation_Add_PML(vv, wf1, wf2, &
					u1_x, u1_x_last, &
					u2_x, u2_x_last, &
					u2_x_temp, u2_x_temp_last, &
					u3_x, u3_x_last, &
					u1_z, u1_z_last, &
					u2_z, u2_z_last, &
					u2_z_temp, u2_z_temp_last, &
					u3_z, u3_z_last, &
					coe, nnx, nnz, nx, nz, &
					nx_bound_l, nx_bound_r, &
					nz_bound_u, nz_bound_d, &
					dx, dz, dt, r, it)

			iit=it-nwt
			if(iit.ge.1) then
				do ix=1, nx_with_apert
					iix=ix
					inx=iix+nx_bound_l
					do iii=1, order
						top(iii,iix,iit)=wf1(nz_bound_u-order/2+iii, inx)
						bot(iii,iix,iit)=wf1(nz_bound_u+nz_with_apert+iii-1-order/2, inx)
					enddo
				enddo

				do iz=1, nz_with_apert
					iiz=iz
					inz=iiz+nz_bound_u
					do iii=1, order
						lef(iii,iiz,iit)=wf1(inz,nx_bound_l-order/2+iii)
						rig(iii,iiz,iit)=wf1(inz,nx_bound_l+nx_with_apert+iii-1-order/2)
					enddo
				enddo
			endif

			!检查波场外推过程是否出现nan情况
			if(isnan(wf1(int(nnz/2),int(nnx/2)))) then
				print*,'it',it, 'the snap appear isnan!!!!!!'
				stop
			endif

		else						!*single time step extrapolation

			!*assignning the source function
			if(it.lt.nwt*2+100)then
				wf1(ns_z+nz_bound_u, ns_x+nx_shift)&  
					=wf1(ns_z+nz_bound_u, ns_x+nx_shift)+ wavelet(it)
			endif

			!*====== single time step extrapolation
			CALL  Extrapolation_Add_PML(vv, wf2, wf1, &
					u1_x_last, u1_x, &
					u2_x_last, u2_x, &
					u2_x_temp_last, u2_x_temp, &
					u3_x_last, u3_x, &
					u1_z_last, u1_z, &
					u2_z_last, u2_z, &
					u2_z_temp_last, u2_z_temp, &
					u3_z_last, u3_z, &
					coe, nnx, nnz, nx, nz, &
					nx_bound_l, nx_bound_r, &
					nz_bound_u, nz_bound_d, &
					dx, dz, dt, r, it)
     
			iit=it-nwt
			if(iit.ge.1) then
				do ix=1, nx_with_apert
					iix=ix
					inx=iix+nx_bound_l
					do iii=1, order
						top(iii,iix,iit)=wf2(nz_bound_u-order/2+iii, inx)
						bot(iii,iix,iit)=wf2(nz_bound_u+nz_with_apert+iii-1-order/2,inx)
					enddo
				enddo
				do iz=1, nz_with_apert
					iiz=iz
					inz=iiz+nz_bound_u
					do iii=1, order
						lef(iii,iiz,iit)=wf2(inz,nx_bound_l-order/2+iii)
						rig(iii,iiz,iit)=wf2(inz,nx_bound_l+nx_with_apert+iii-1-order/2)
					enddo
				enddo
			endif

		endif						!*end single time step extrapolation

	enddo						!*loop end for forward source wavefield 
	!*======================================================================*


	!***********************************************************************!
	!*	Extrapolating the extrapolated source wavfield and the receiver 
	!*wavefield along the reverse time and extracting the imaging value 
	!*with the correlation of the two wavefield at each time step
	!***********************************************************************!

    !Reinitiallization for back propagation
    u1_x=0.0
    u1_x_last=0.0
    u2_x=0.0
    u2_x_last=0.0
    u2_x_temp=0.0
    u2_x_temp_last=0.0
    u3_x=0.0
    u3_x_last=0.0

    u1_z=0.0
    u1_z_last=0.0
    u2_z=0.0
    u2_z_last=0.0
    u2_z_temp=0.0
    u2_z_temp_last=0.0
    u3_z=0.0
    u3_z_last=0.0
!	wf1=0.0
!	wf2=0.0

	!***************************************************************************!
	istep=0
    do it=lt, 1, -1			!!!! loop for backrward source and receiver wavefield
		istep=istep+1

        if(mod(it, 1000).eq.1) then
			write(*,*)'shot=',iss,'myid=',myid,'  step=reverse  ','it=',it
        endif

		!if(modulo(it+500, 500).eq.1 .and. myid .eq. 999999 ) then
		!	write(currt,'(I4)')it

		!	write(currtsnap,*)trim(adjustl(snapx_fn)),trim(adjustl(currt))
		!	open(133, file=trim(adjustl(currtsnap))//'reverse.dat', access='direct',&
		!		form='unformatted', status='replace', recl=lbyte*nx_with_apert*nz_with_apert, iostat=ierr)
		!	write(133,rec=1)((wf2(iz, ix),iz=1+nz_bound_u,nz_with_apert+nz_bound_u), &
		!				ix=1+nx_bound_l,nx_with_apert+nx_bound_l)
		!	close(133)

		!endif


        if(mod(istep,2).eq.1)then	!!!!!single time step extrapolation	

			!*single time step -- extrapolation source wavefield reverse time --------*      
			do ix=1, nx_with_apert
				iix=ix
				inx=iix+nx_bound_l
				do iii=1, order
					wf2(nz_bound_u-order/2+iii, inx)=top(iii,iix,it)
					wf2(nz_bound_u+nz_with_apert+iii-1-order/2,inx)=bot(iii,iix,it)
				enddo
            enddo
            do iz=1, nz_with_apert
                iiz=iz
                inz=iiz+nz_bound_u
				do iii=1, order
					wf2(inz,nx_bound_l-order/2+iii)=lef(iii,iiz,it)
					wf2(inz,nx_bound_l+nx_with_apert+iii-1-order/2)=rig(iii,iiz,it)
				enddo
            enddo
             
			CALL Extrapolation_Without_PML(vv, wf1, wf2, coe, nnx, nnz, &
					nx_bound_l,nx_bound_r, nz_bound_u,nz_bound_d,&
					dx, dz, dt)

			!*single time step -- extrapolation recesiver wavefield reverse time  --*   
			do ix = 1, nx
				ixx = ix + nx_shift
				wr2(nr_z+nz_bound_u, ixx) = &
					wr2(nr_z+nz_bound_u, ixx) + shot_gather(it, ix)
			end do

			CALL  Extrapolation_Add_PML(vv, wr1, wr2, &
					u1_x, u1_x_last, &
					u2_x, u2_x_last, &
					u2_x_temp, u2_x_temp_last, &
					u3_x, u3_x_last, &
					u1_z, u1_z_last, &
					u2_z, u2_z_last, &
					u2_z_temp, u2_z_temp_last, &
					u3_z, u3_z_last, &
					coe, nnx, nnz, nx, nz, &
					nx_bound_l, nx_bound_r, &
					nz_bound_u, nz_bound_d, &
					dx, dz, dt, r, it)

			!*cross correlation image condition at every time slice
			if(it .gt. 2*nwt+1)then
				CALL Cor_Imaging_Condition(vv, image, &
						wf1, wr1, &
						nnx, nnz, nz, &
						nx_with_apert, &
						nx_bound_l, nz_bound_u, &
						dt, myid)
			endif

        else					!!!!!single time step extrapolation

			!*ingle time step -- extrapolation source wavefield reverse time -------*      
			do ix=1, nx_with_apert
				iix=ix
				inx=iix+nx_bound_l
				do iii=1, order
					wf1(nz_bound_u-order/2+iii, inx)=top(iii,iix,it)
					wf1(nz_bound_u+nz_with_apert+iii-1-order/2,inx)=bot(iii,iix,it)
				enddo
            enddo
            do iz=1, nz_with_apert
                iiz=iz
                inz=iiz+nz_bound_u
				do iii=1, order
					wf1(inz,nx_bound_l-order/2+iii)=lef(iii,iiz,it)
					wf1(inz,nx_bound_l+nx_with_apert+iii-1-order/2)=rig(iii,iiz,it)
				enddo
            enddo

			CALL Extrapolation_Without_PML(vv, wf2, wf1, coe, nnx, nnz, &
					nx_bound_l, nx_bound_r, nz_bound_u, nz_bound_d,&
					dx, dz, dt)

			!*single time step -- extrapolation receiver wavefield reverse time  ---*   
			do ix = 1, nx
				ixx = ix + nx_shift
				wr1(nr_z+nz_bound_u, ixx) = &
					wr1(nr_z+nz_bound_u, ixx) + shot_gather(it, ix)
			end do

			CALL  Extrapolation_Add_PML(vv, wr2, wr1, &
					u1_x_last, u1_x, &
					u2_x_last, u2_x, &
					u2_x_temp_last, u2_x_temp, &
					u3_x_last, u3_x, &
					u1_z_last, u1_z, &
					u2_z_last, u2_z, &
					u2_z_temp_last, u2_z_temp, &
					u3_z_last, u3_z, &
					coe, nnx, nnz, nx, nz, &
					nx_bound_l, nx_bound_r, &
					nz_bound_u, nz_bound_d, &
					dx, dz, dt, r, it)

			!*cross correlation image condition at every time slice
			if(it .gt. 2*nwt+1)then
				CALL Cor_Imaging_Condition(vv, image, &
						wf2, wr2, &
						nnx, nnz, nz, &
						nx_with_apert, &
						nx_bound_l, nz_bound_u, &
						dt, myid)
			endif

        endif				!!!!end single time step extrapolation

	enddo				!!! loop end for backward source and receiver wavefield
	!*===================================================================

	!*Writing the current shot imaging result onto hard disk
    CALL Write_Disk(image, image_tmp, &
			nx_with_apert, nvx, nvxx_shift, nz, &
			krec_save, myid, iss)

	!*======================================================================*
	!*     The current shot gather reverse time migration finished          *
	!*======================================================================*

    return
    end	subroutine

!***********************************************************************
	subroutine Reduce_Result_Gradient(fn3, myid, master, nvx, nz)

	use global
	implicit none
	include 'mpif.h'
    integer status(mpi_status_size)
	!Data dictionary: declare variable types, definitons, & units
	!Dummy Variables
	character(len=256)	::	fn3
	integer ::	myid, master
	integer ::	nz
	integer ::	nvx
	!Local Variables
	integer ::	ix, iz
	integer ::	nxz
	integer ::	ierr
	integer ::	irec

	real,allocatable ::	image_all(:,:)
	real,allocatable ::	image_single(:,:)


    if(myid .eq. master) then
		open(13, file=fn3, access='direct', status='replace', recl=lbyte*nz, iostat=ierr)
		if(ierr /= 0)then
			write(*,*)  'Image function file cannot open right, please check it !!!'
			write(*,*)  'fn3=',fn3
		endif
    endif

    nxz=nvx*nz

    allocate (image_all(nz, nvx), stat=ierr)
	if(ierr.ne.0)then
		write(*,*)"Can not allocate working memory about image_all, stop!!!"
		stop
	endif

    allocate (image_single(nz, nvx),stat=ierr)
	if(ierr.ne.0)then
		write(*,*)"Can not allocate working memory about image_single, stop!!!"
		stop
	endif

	do ix=1, nvx
		do iz=1, nz
			image_single(iz, ix)=0.0
		enddo
    enddo

    if(myid .ne. master)then
		do ix=1, nvx
			read(14, rec=ix, err=555)(image_single(iz, ix),iz=1, nz)
		enddo
    endif
555 continue

	CALL MPI_REDUCE(IMAGE_SINGLE, IMAGE_ALL, NXZ, MPI_REAL, MPI_SUM, &
			MASTER, MPI_COMM_WORLD, IERR)

	if(myid .eq. master) then
		do ix=1, nvx
			write(13, rec=ix)(image_all(iz, ix),iz=1, nz)
		enddo
    endif       

    deallocate(image_all)
    deallocate(image_single) 

    if(myid .eq. master) then
		close(13)
    endif

	return

	end subroutine

!*********************************************************************** 
	subroutine Cut_Mute_Direct_Wave(shot_gather, vv_surface, &
					ns_x, nx, nx_with_apert, &
					nz, lt, nx_apert_l, nwt, &
					max_offset_x, dx, dt )

!*======================================================================
!*  Variable demonstrations:
!*		shot_gather:		shot gather
!*		vv_surface:			surface velocity
!*		ns_x:				the number of source point location of x axis
!*		nx:					the dimension of current migration profile of x axis
!*		nx_with_apert:		the dimension of x axis of the current shot gather after
!*								adding the migration aperture
!*		nz:					the depth of current migration profile
!*		lt:					the length of record
!*		nx_apert_l:			the sample point number of the added migration aperture at left of x axis
!*		nwt:				the half dimension of wavelet width
!*		max_offset_x:		the maxinum offset of x axis  (m)
!*		dx:					the sample interval of x axis (m)
!*		dt:					the sample interval of record (s)
!*
!*======================================================================

	implicit none
	!Data dictionary: declare variable types, definitons, & units
	!Dummy Variables
	integer ::	ns_x	
	integer ::	nx		
	integer ::	nx_with_apert
	integer ::	nz		
	integer ::	lt		
    integer ::	nx_apert_l
	integer ::	nwt
	real	::	max_offset_x
	real	::	dx
	real	::	dt
	real	::	vv_surface(nx_with_apert)
	real	::	shot_gather(lt, nx)
	!Local Variables
	integer ::	ix,iz
	integer ::	iix
	integer ::	it
	integer ::	itmp
	integer ::	mute(nx_with_apert)
	real	::	temp

	if(ns_x .ge.1 .and. ns_x .le. nx_with_apert)then
		do ix=1, nx_with_apert
			if(abs(vv_surface(ix)).lt.1.0e-4)then
				write(*,*)	'The velocity of input euqal zero, please modifty it'
				stop
			endif
			temp=0.0
			if(ix .le. ns_x)then
				do iix=ix, ns_x
					temp=temp+dx/vv_surface(iix)
				enddo
			else
				do iix=ix, ns_x, -1
					temp=temp+dx/vv_surface(iix)
				enddo
			endif
			mute(ix)=nint(temp/dt)+0.5
			if( abs(dx*(ix-ns_x)) .lt. max_offset_x/4 )then
				mute(ix)=mute(ix)+2*nwt+101
			else
				mute(ix)=mute(ix)+2*nwt+1 
			endif
			if(mute(ix).lt.0)	mute(ix)=0
			if(mute(ix).gt.lt)	mute(ix)=lt
		enddo
	else
		write(*,*)	'The ns_x coordinate is error!!!'
		stop
	endif		
	
	do ix=1, nx
		iix=ix+nx_apert_l
		do it=1, lt
			if(it .lt. mute(iix) )then
				shot_gather(it,ix)=0.0
			endif
		enddo
	enddo
		
	return
	
	end subroutine

!***********************************************************************
    subroutine Write_Disk(image, image_tmp, &
					nx_with_apert, nvx, nvxx_shift, nz, &
					krec_save, myid, iss)

!*======================================================================
!*  Variable demonstrations:
!*		image:			the current image result
!*		iamge_tmp:		the tmp stack image result
!*		nx_with_apert:	the dimension of x axis of the current shot gather after
!*							adding the migration aperture
!*		nvx:			the dimension of total velocity of x axis
!*		nvxx_shift:		the shift about current image range and total image result
!*		nz:					the depth of current migration profile
!*
!*======================================================================
	implicit none
	!Data dictionary: declare variable types, definitons, & units
	!Dummy Variables
	integer ::	nx_with_apert
	integer ::	nvx
	integer ::	nz
	integer ::	nvxx_shift
	integer ::	krec_save
	integer ::	myid
	integer ::	iss

	real	::	image(nz, nx_with_apert)
	real	::	image_tmp(nz, nvx)
	!Local Variables
	integer ::	ix, iy, iz
	integer ::	iix, iiy
    integer	::	irec

	write(*,*)'shot',iss,'is writing disk!','myid=',myid

    do ix=1, nvx
        do iz=1, nz
            image_tmp(iz, ix)=0.0
        enddo
	enddo

    do ix=1, nvx
        irec=ix
        read(14, rec=irec, err=111)(image_tmp(iz, ix), iz=1, nz)
    enddo

111 continue

    do ix=1, nx_with_apert
        iix= ix + nvxx_shift
        if(iix .ge. 1 .and. iix .le. nvx) then
            do iz=1, nz
				image_tmp(iz, iix)=image_tmp(iz, iix) + image(iz, ix)
            enddo
        endif
    enddo

    do ix=1,nvx
        irec=ix
        write(14, rec=irec)(image_tmp(iz, ix), iz=1, nz)
    enddo

    write(*,*)'shot',iss,'is done!','myid=',myid

    return

    end	subroutine

!***********************************************************************
	subroutine Cor_Imaging_Condition(vv, image, &
					wf1, wr1, &
					nnx, nnz, nz, &
					nx_with_apert, &
					nx_bound_l, nz_bound_u, &
					dt, myid)

!*======================================================================
!*  Variable demonstrations:
!*		vv:				the working buffer for the current shot velocty field
!*		image:			the current image result
!*
!*		wf1, wr1:		the wavefiled of scoure and receiver
!*
!*		nnx, nnz:		the wave filed extrapolation range after adding the
!							absorbing boundary
!*		nz:				the depth of current migration profile
!*
!*		nx_with_apert:	the dimension of x axis of the current shot gather after
!*							adding the migration aperture
!*
!*		nx_bound_l:		the sample point number of the added absorbing boundary range
!*		nz_bound_u:		the sample point number of the added absorbing boundary range
!*
!*		dt:				the sample rate of  extrapolation time  
!*
!*======================================================================
	implicit none
	!Data dictionary: declare variable types, definitons, & units
	!Dummy Variables
	integer ::	nnx, nnz
	integer ::	nz
	integer ::	nx_with_apert
	integer ::	nx_bound_l
	integer ::	nz_bound_u
	integer ::	myid
	real	::	dt
	real	::	vv(nnz, nnx)
	real	::	image(nz, nx_with_apert)
	real	::	wf1(-4:nnz+5,-4:nnx+5)
	real	::	wr1(-4:nnz+5,-4:nnx+5)

	!Local Variables
	integer	::	ix, ixx
	integer ::	iz, izz
	real	::	tmp


        
	do ix=1, nx_with_apert
		ixx = ix + nx_bound_l
		do iz=1, nz
			izz=iz+nz_bound_u
			
			!**********************************************************!
			!互相关成像条件
			tmp=wf1(izz, ixx)*wr1(izz, ixx)
			image(iz, ix)=image(iz, ix) + tmp
			!**********************************************************!
        enddo
	enddo


    return
    end	subroutine

!***********************************************************************
	subroutine Extrapolation_Without_PML(vv, x1, x2, coe, nnx, nnz,&
				nx_bound_l,nx_bound_r, nz_bound_u,nz_bound_d,&
				dx, dz, dt)

    implicit none
	!Data dictionary: declare variable types, definitons, & units
    !dummy variables
    integer	::	nnx, nnz
	integer ::	nx_bound_l, nx_bound_r
    integer ::	nz_bound_u, nz_bound_d
    real	::	dt, dx, dz
    real	::	vv(nnz, nnx)
	real	::	coe(13)
	real	::	x1(-4:nnz+5, -4:nnx+5)
    real	::	x2(-4:nnz+5, -4:nnx+5)

    !local variables
    integer	::	ix, iz
	real	::	ddt, vv2

	!*=================================================================
	!*          now do the wavefield extrapolation
	!*=================================================================
    ddt=dt
    
!$omp parallel private(vv2) 
!$omp do

 	do ix=nx_bound_l+1, nnx-nx_bound_r
        do iz=nz_bound_u+1, nnz-nz_bound_d

			vv2=vv(iz, ix)*vv(iz, ix)

			x1(iz,ix) = coe(1)*x2(iz, ix)&
                      + coe(2)*x1(iz, ix)&
                      + coe(3)*vv2*x2(iz, ix)&
                      + coe(4)*vv2*(x2(iz, ix+1) + x2(iz, ix-1))&
                      + coe(5)*vv2*(x2(iz, ix+2) + x2(iz, ix-2))&
                      + coe(6)*vv2*(x2(iz, ix+3) + x2(iz, ix-3))&
                      + coe(7)*vv2*(x2(iz, ix+4) + x2(iz, ix-4))&
                      + coe(8)*vv2*(x2(iz, ix+5) + x2(iz, ix-5))&
                      + coe(9)*vv2*(x2(iz+1, ix) + x2(iz-1, ix))&
                      + coe(10)*vv2*(x2(iz+2, ix) + x2(iz-2, ix))&
                      + coe(11)*vv2*(x2(iz+3, ix) + x2(iz-3, ix))&
                      + coe(12)*vv2*(x2(iz+4, ix) + x2(iz-4, ix))&
                      + coe(13)*vv2*(x2(iz+5, ix) + x2(iz-5, ix))

		enddo
    enddo

!$omp end do
!$omp end parallel
    return

    end subroutine

!***********************************************************************
    subroutine  Extrapolation_Add_PML(vv, x1, x2, &
			u1_x, u1_x_present, &
			u2_x, u2_x_present, &
			u2_x_temp, u2_x_temp_present, &
			u3_x, u3_x_present, &
			u1_z, u1_z_present, &
			u2_z, u2_z_present, &
			u2_z_temp, u2_z_temp_present, &
			u3_z, u3_z_present, &
			coe, nnx, nnz, nx, nz, &
			nx_bound_l, nx_bound_r, &
			nz_bound_u, nz_bound_d, &
			dx, dz, dt, r, it)

    implicit none
	!Data dictionary: declare variable types, definitons, & units
    !Dummy variables
    integer ::	nnx, nnz
    integer ::	nx, nz
    integer ::	nx_bound_l, nx_bound_r
	integer ::	nz_bound_u, nz_bound_d
    real    ::	dx, dz, dt
    real    ::	r
	integer ::	it
  
    real	::	vv(nnz, nnx)
    real	::	x1(-4:nnz+5, -4:nnx+5)
    real	::	x2(-4:nnz+5, -4:nnx+5)
    real	::	coe(13)

    !wavefield within pml field  

    !x direction
    real	::	u1_x(nnz,nx_bound_l+nx_bound_r)
    real	::	u1_x_present(nnz,nx_bound_l+nx_bound_r)
    real	::	u2_x(nnz,nx_bound_l+nx_bound_r)
    real	::	u2_x_present(nnz,nx_bound_l+nx_bound_r)
    real	::	u2_x_temp(nnz,nx_bound_l+nx_bound_r)
    real	::	u2_x_temp_present(nnz,nx_bound_l+nx_bound_r)
    real	::	u3_x(nnz,nx_bound_l+nx_bound_r)
    real	::	u3_x_present(nnz,nx_bound_l+nx_bound_r)

    !z direction
    real	::	u1_z(nz_bound_u+nz_bound_d,nnx)
    real	::	u1_z_present(nz_bound_u+nz_bound_d,nnx)
    real	::	u2_z(nz_bound_u+nz_bound_d,nnx)
    real	::	u2_z_present(nz_bound_u+nz_bound_d,nnx)
    real	::	u2_z_temp(nz_bound_u+nz_bound_d,nnx)
    real	::	u2_z_temp_present(nz_bound_u+nz_bound_d,nnx)
    real	::	u3_z(nz_bound_u+nz_bound_d,nnx)
    real	::	u3_z_present(nz_bound_u+nz_bound_d,nnx)
 
    !Local variables
    integer ::	ix,iz
    real    ::	vv2
    real	::	ddt

    integer ::	grid_x_positive_pml,grid_x_negative_pml
    integer ::	grid_z_positive_pml,grid_z_negative_pml

    real	::	axis_x_positive_pml,axis_x_thickness_pml,axis_x_negative_pml
    real	::	axis_z_positive_pml,axis_z_thickness_pml,axis_z_negative_pml

    real	::	d_attenuation_x_present,deri_d_attenuation_x
    real	::	d_attenuation_z_present,deri_d_attenuation_z

	integer ::	iflag
      
	ddt=dt

	iflag=0

	!*=================================================================
	!*          now do the wavefield extrapolation
	!*=================================================================

!$omp parallel private(vv2,&
!$omp         d_attenuation_x_present,deri_d_attenuation_x,&
!$omp         d_attenuation_z_present,deri_d_attenuation_z,&
!$omp         axis_x_negative_pml,axis_x_thickness_pml,&
!$omp         grid_x_negative_pml,axis_x_positive_pml,&
!$omp         grid_x_positive_pml,&
!$omp         axis_z_negative_pml,axis_z_thickness_pml,&
!$omp         grid_z_negative_pml,axis_z_positive_pml,&
!$omp         grid_z_positive_pml)

!$omp do 

    do ix=1, nnx
        do iz=1, nnz
			if(iz>nz_bound_u.and.iz<nnz-nz_bound_d+1.and.&
				ix>nx_bound_l.and.ix<nnx-nx_bound_r+1)then

			vv2=vv(iz, ix)*vv(iz, ix)
			x1(iz,ix) = coe(1)*x2(iz, ix)&
					  + coe(2)*x1(iz, ix)&
                      + coe(3)*vv2*x2(iz, ix)&
                      + coe(4)*vv2*(x2(iz, ix+1) + x2(iz, ix-1))&
                      + coe(5)*vv2*(x2(iz, ix+2) + x2(iz, ix-2))&
                      + coe(6)*vv2*(x2(iz, ix+3) + x2(iz, ix-3))&
                      + coe(7)*vv2*(x2(iz, ix+4) + x2(iz, ix-4))&
                      + coe(8)*vv2*(x2(iz, ix+5) + x2(iz, ix-5))&
                      + coe(9)*vv2*(x2(iz+1, ix) + x2(iz-1, ix))&
                      + coe(10)*vv2*(x2(iz+2, ix) + x2(iz-2, ix))&
                      + coe(11)*vv2*(x2(iz+3, ix) + x2(iz-3, ix))&
                      + coe(12)*vv2*(x2(iz+4, ix) + x2(iz-4, ix))&
                      + coe(13)*vv2*(x2(iz+5, ix) + x2(iz-5, ix))
             

			!absorbing boundary condition
			!x_negative
			else if(ix<=nx_bound_l)then
				axis_x_negative_pml=real(ix-(nx_bound_l+1))*dx
				axis_x_thickness_pml=real(nx_bound_l)*dx
				grid_x_negative_pml=ix
              
				d_attenuation_x_present=(3.0*vv(iz,ix)/(2.0*&
					axis_x_thickness_pml))*(axis_x_negative_pml/&
					axis_x_thickness_pml)**2*log(1/r)

				deri_d_attenuation_x=(3.0*vv(iz,ix)/(2.0*&
					(axis_x_thickness_pml)))*(2.0*axis_x_negative_pml/&
					axis_x_thickness_pml**2)*log(1/r)


				!pay attention to the subsripts of x2 etc.
				u1_x(iz,ix)=(u1_x_present(iz,ix)*(2*dx**2+2*&
					d_attenuation_x_present*ddt*dx**2-&
					d_attenuation_x_present**2*dx**2*ddt**2)-dx**2*&
					u1_x(iz,ix)+vv(iz,ix)**2*ddt**2*&
					(x2(iz,ix+1)+&
					x2(iz,ix-1)-2*x2(iz,ix)))/&
					(dx**2+2*d_attenuation_x_present*ddt*dx**2)


				u2_x_temp(iz,ix)=(u2_x_temp_present(iz,ix)*(2*dx+2*&
					d_attenuation_x_present*ddt*dx-ddt**2*dx*&
					d_attenuation_x_present**2)-dx*u2_x_temp(iz,ix)&
					-vv(iz,ix)**2*ddt**2*deri_d_attenuation_x*&
					(x2(iz,ix+1)-x2(iz,ix)))/&
					(dx+2*d_attenuation_x_present*ddt*dx)

              
				u2_x(iz,ix)=u2_x_present(iz,ix)*&
					(1-d_attenuation_x_present*ddt)&
					+u2_x_temp(iz,ix)*ddt
              
				u3_x(iz,ix)=2*u3_x_present(iz,ix)-&
					u3_x(iz,ix)+vv(iz,ix)**2*&
					ddt**2*(x2(iz+1,ix)+x2(iz-1,ix)-2*&
					x2(iz,ix))/dz**2
              
				x1(iz,ix)=u1_x(iz,ix)+u2_x(iz,ix)+u3_x(iz,ix)

			!x_positive  
			else if(ix>=nnx-nx_bound_r+1)then
				axis_x_positive_pml=real(ix-(nnx-nx_bound_r))*dx
				axis_x_thickness_pml=real(nx_bound_r)*dx
				grid_x_positive_pml=ix-(nnx-nx_bound_r)+nx_bound_l
              
				d_attenuation_x_present=(3.0*vv(iz,ix)/(2.0*&
					axis_x_thickness_pml))*(axis_x_positive_pml/&
					axis_x_thickness_pml)**2*log(1/r)

				deri_d_attenuation_x=(3.0*vv(iz,ix)/(2.0*&
					(axis_x_thickness_pml)))*(2.0*axis_x_positive_pml/&
					axis_x_thickness_pml**2)*log(1/r)

				!pay attention to the subsripts of x2 etc.
				u1_x(iz,grid_x_positive_pml)=&
					(u1_x_present(iz,grid_x_positive_pml)*(2*dx**2+2*&
					d_attenuation_x_present*ddt*dx**2-&
					d_attenuation_x_present**2*dx**2*ddt**2)-dx**2*&
					u1_x(iz,grid_x_positive_pml)+vv(iz,ix)**2*&
					ddt**2*(x2(iz,ix+1)+&
					x2(iz,ix-1)-2*x2(iz,ix)))/&
					(dx**2+2*d_attenuation_x_present*ddt*dx**2)
              
				u2_x_temp(iz,grid_x_positive_pml)=&
					(u2_x_temp_present(iz,grid_x_positive_pml)*&
					(2*dx+2*d_attenuation_x_present*ddt*dx-ddt**2*dx*&
					d_attenuation_x_present**2)-dx*&
					u2_x_temp(iz,grid_x_positive_pml)&
					-vv(iz,ix)**2*ddt**2*deri_d_attenuation_x*&
					(x2(iz,ix+1)-x2(iz,ix)))/&
					(dx+2*d_attenuation_x_present*ddt*dx)
              
				u2_x(iz,grid_x_positive_pml)=&
					u2_x_present(iz,grid_x_positive_pml)*&
					(1-d_attenuation_x_present*ddt)&
					+u2_x_temp(iz,grid_x_positive_pml)*ddt
              
				u3_x(iz,grid_x_positive_pml)=&
					2*u3_x_present(iz,grid_x_positive_pml)-&
					u3_x(iz,grid_x_positive_pml)+vv(iz,ix)**2*&
					ddt**2*(x2(iz+1,ix)+x2(iz-1,ix)-2*&
					x2(iz,ix))/dz**2
              
				x1(iz,ix)=u1_x(iz,grid_x_positive_pml)+&
					u2_x(iz,grid_x_positive_pml)+&
					u3_x(iz,grid_x_positive_pml)

        
			!z_negative
			else if(iz<=nz_bound_u)then
				axis_z_negative_pml=real(iz-(nz_bound_u+1))*dz
				axis_z_thickness_pml=real(nz_bound_u)*dz
				grid_z_negative_pml=iz
              
				d_attenuation_z_present=(3.0*vv(iz,ix)/(2.0*&
					axis_z_thickness_pml))*(axis_z_negative_pml/&
					axis_z_thickness_pml)**2*log(1/r)

				deri_d_attenuation_z=(3.0*vv(iz,ix)/(2.0*&
					(axis_z_thickness_pml)))*(2.0*axis_z_negative_pml/&
					axis_z_thickness_pml**2)*log(1/r)

				!pay attention to the subsripts of x2 etc.
				u1_z(iz,ix)=(u1_z_present(iz,ix)*(2*dz**2+2*&
					d_attenuation_z_present*ddt*dz**2-&
					d_attenuation_z_present**2*dz**2*ddt**2)-dz**2*&
					u1_z(iz,ix)+vv(iz,ix)**2*ddt**2*&
					(x2(iz+1,ix)+&
					x2(iz-1,ix)-2*x2(iz,ix)))/&
					(dz**2+2*d_attenuation_z_present*ddt*dz**2)

				u2_z_temp(iz,ix)=(u2_z_temp_present(iz,ix)*(2*dz+2*&
					d_attenuation_z_present*ddt*dz-ddt**2*dz*&
					d_attenuation_z_present**2)-dz*u2_z_temp(iz,ix)&
					-vv(iz,ix)**2*ddt**2*deri_d_attenuation_z*&
					(x2(iz+1,ix)-x2(iz,ix)))/&
					(dz+2*d_attenuation_z_present*ddt*dz)
              
				u2_z(iz,ix)=u2_z_present(iz,ix)*&
					(1-d_attenuation_z_present*ddt)&
					+u2_z_temp(iz,ix)*ddt
              
				u3_z(iz,ix)=2*u3_z_present(iz,ix)-&
					u3_z(iz,ix)+vv(iz,ix)**2*&
					ddt**2*(x2(iz,ix+1)+x2(iz,ix-1)-2*&
					x2(iz,ix))/dx**2
              
				x1(iz,ix)=u1_z(iz,ix)+u2_z(iz,ix)+u3_z(iz,ix)

			!z_positive
			else if(iz>=nnz-nz_bound_d+1)then
				axis_z_positive_pml=real(iz-(nnz-nz_bound_d))*dz
				axis_z_thickness_pml=real(nz_bound_d)*dz
				grid_z_positive_pml=iz-(nnz-nz_bound_d)+nz_bound_u

				d_attenuation_z_present=(3.0*vv(iz,ix)/(2.0*&
					axis_z_thickness_pml))*(axis_z_positive_pml/&
					axis_z_thickness_pml)**2*log(1/r)

				deri_d_attenuation_z=(3.0*vv(iz,ix)/(2.0*&
					(axis_z_thickness_pml)))*(2.0*axis_z_positive_pml/&
					axis_z_thickness_pml**2)*log(1/r)

				!pay attention to the subsripts of x2 etc.
				u1_z(grid_z_positive_pml,ix)=&
					(u1_z_present(grid_z_positive_pml,ix)*(2*dz**2+2*&
					d_attenuation_z_present*ddt*dz**2-&
					d_attenuation_z_present**2*dz**2*ddt**2)-dz**2*&
					u1_z(grid_z_positive_pml,ix)+vv(iz,ix)**2*&
					ddt**2*(x2(iz+1,ix)+&
					x2(iz-1,ix)-2*x2(iz,ix)))/&
					(dz**2+2*d_attenuation_z_present*ddt*dz**2)

				u2_z_temp(grid_z_positive_pml,ix)=&
					(u2_z_temp_present(grid_z_positive_pml,ix)*(2*dz+2*&
					d_attenuation_z_present*ddt*dz-ddt**2*dz*&
					d_attenuation_z_present**2)-dz*&
					u2_z_temp(grid_z_positive_pml,ix)&
					-vv(iz,ix)**2*ddt**2*deri_d_attenuation_z*&
					(x2(iz+1,ix)-x2(iz,ix)))/&
					(dz+2*d_attenuation_z_present*ddt*dz)
              
				u2_z(grid_z_positive_pml,ix)=&
					u2_z_present(grid_z_positive_pml,ix)*&
					(1-d_attenuation_z_present*ddt)&
					+u2_z_temp(grid_z_positive_pml,ix)*ddt
              
				u3_z(grid_z_positive_pml,ix)=&
					2*u3_z_present(grid_z_positive_pml,ix)-&
					u3_z(grid_z_positive_pml,ix)+vv(iz,ix)**2*&
					ddt**2*(x2(iz,ix+1)+x2(iz,ix-1)-2*&
					x2(iz,ix))/dx**2
              
				x1(iz,ix)=u1_z(grid_z_positive_pml,ix)+&
					u2_z(grid_z_positive_pml,ix)+&
					u3_z(grid_z_positive_pml,ix)

			endif
		enddo
    enddo

!$omp end do
!$omp end parallel

    return
    end subroutine

!***********************************************************************
    subroutine Get_Current_Shot_Velocity(vv_tmp, vv, vv_surface, &
				nvx, nvz, dvx, dvz, &
                nx_with_apert, nnx, nnz, &
                nx_bound_l, nz_bound_u, &
                nx_bound_r, nz_bound_d, &
                nx, nz, &
                nvxx_shift)

!*=======================================================================
!*	Variable demonstrations:
!*
!*		vv_all: the total velocity field
!*		vv: the working buffer for the current shot velocty field
!*		vv_surface: the working buffer of surface for the current shot velocty field
!*
!*		nvx, nvz: the dimension of the total velocity field
!*		dvxm dvz: the interval of  the velocoty field
!*
!*		nx_with_apert: the dimension of the current shot gather after
!*							adding the migration aperture, they define 
!*							the imaging range of the current shot gather
!*
!*		nnx,nnz: the wave filed extrapolation range after adding 
!*                   the absorbing boundary 
!*    
!*		nx_bound, nz_bound : the sample point number of the added 
!*                              absorbing boundary range
!*
!*		nx, nz:	the dimension of the current shot gather after
!*
!*		nvxx_shift: the shift between the current shot imaging
!*                      range and the total velocity field range
!* 
!*=======================================================================
	implicit none
	!Data dictionary: declare variable types, definitons, & units
	!Dummy Variables
    integer ::	nnx, nnz
    integer ::	nvx, nvz
    integer ::	nx_with_apert
    integer ::	nvxx_shift
    integer ::	nx_bound_l, nz_bound_u
    integer ::	nx_bound_r, nz_bound_d
	integer	::	nx, nz
	real	::	dvx, dvz
	real	::	vv_tmp(nvz)
    real	::	vv(nnz, nnx)
    real	::	vv_surface(nx_with_apert)
	!Local Variables
	integer ::	ix, ixx, iix 
	integer ::	iz, iiz
	integer ::	ierr

    do ix=1, nx_with_apert
		ixx=ix+nvxx_shift
        iix=ix+nx_bound_l
        if(ixx.lt.1) ixx = 1
        if(ixx.gt.nvx) ixx = nvx

		vv_tmp=0.0
		read(12, rec=ixx)(vv_tmp(iz),iz=1, nvz)

		do iz=1, nz
			iiz=iz+nz_bound_u
			vv(iiz, iix)=vv_tmp(iz)
		enddo
    enddo

    !padding the velocity 
    do iiz=1,nz_bound_u
		do iix=1,nnx
			vv(iiz,iix)=vv(nz_bound_u+1,iix)
        enddo
    enddo
    do iiz=nnz-nz_bound_d+1,nnz
        do iix=1,nnx
			vv(iiz,iix)=vv(nnz-nz_bound_d,iix)
        enddo
    enddo   

    do iix=1,nx_bound_l
        do iiz=1,nnz
			vv(iiz,iix)=vv(iiz,nx_bound_l+1)
        enddo
    enddo
    do iix=nnx-nx_bound_r+1,nnx
        do iiz=1,nnz
			vv(iiz,iix)=vv(iiz,nnx-nx_bound_r)
        enddo
    enddo
	
	!vv_surface
	do ix=1, nx_with_apert
		vv_surface(ix)=vv(nz_bound_u+1,nx_bound_l+ix)
	enddo

    return
    end subroutine

!***********************************************************************
    subroutine	Index_Data(ns_start, ns_total, ns_interval,&
              cx_min, cx_max, nshot, index_table, lt_real, &
			  max_offset_x, min_offset_x)

	use global
	use header_module
	implicit none
	!Dummy Variables
    integer	::	ns_start, ns_total, ns_interval
	integer	::	cx_min, cx_max
	integer	::	nshot
	integer ::  lt_real
	real	::  max_offset_x
	real	::  min_offset_x
	integer	::	index_table(9, ns_total)

	!Local	Variables
	type(segy)  ::  head
	integer	::	ns_end
	integer ::	nss
	integer ::	irec
	integer ::	ierr
	integer ::	icount
	integer	::	is
	integer ::	ntrace
	integer	::	krec_save
    integer ::	csx, cgx
	integer ::	cxmin, cxmax
	integer ::	iflag

	iflag=1		!the flag of controlling the print  on-off


    ns_end=ns_start+(ns_total-1)*ns_interval
    irec=0
    icount=1

    do is=ns_start, ns_end, ns_interval
		cx_min=9999999
        cx_max=-9999999
        ntrace=0
        krec_save=0
1111    continue
        irec=irec+1
		
        read(11, rec=irec, err=5555)head
		nss=head%fldr
        if(nss.lt.is)  then
			goto 1111
        elseif(nss.eq.is) then
			ntrace=ntrace+1
			if(krec_save.eq.0) krec_save=irec

			csx=head%sx
			cgx=head%gx

			!判断偏移距大小是否在预设的范围内
			if(abs(csx-cgx) .gt. max_offset_x .or. abs(csx-cgx) .lt. min_offset_x)then
				goto 1111
			endif

			cxmin=min(csx, cgx)
			if(cxmin.lt.cx_min) cx_min=cxmin
			cxmax=max(csx, cgx)
			if(cxmax.gt.cx_max) cx_max=cxmax

			goto 1111
        elseif(nss.gt.is) then
			if (ntrace.eq.0) then
				irec=irec-1
				goto 2222
			endif
			if(ntrace/=0)then
				irec=irec-1
				index_table(1, icount)=is
				index_table(2, icount)=ntrace
				index_table(3, icount)=csx
				index_table(4, icount)=0
				index_table(5, icount)=cx_min
				index_table(6, icount)=cx_max
				index_table(7, icount)=0
				index_table(8, icount)=0
				index_table(9, icount)=krec_save
				nshot=icount
				icount=icount+1
			else
				write(*,*)'shot',is,'is absent'
			endif
        endif
2222	continue
    enddo
    go to 6666
5555	continue
    index_table(1, icount)=is
    index_table(2, icount)=ntrace
    index_table(3, icount)=csx
    index_table(4, icount)=0
    index_table(5, icount)=cx_min
    index_table(6, icount)=cx_max
    index_table(7, icount)=0
    index_table(8, icount)=0
    index_table(9, icount)=krec_save
    nshot=icount
6666	continue

	cx_min=999999999
	cx_max=-999999999

	do is=1, nshot
		cxmin=index_table(5, is)
		cxmax=index_table(6, is)
		if(cxmin.lt.cx_min) cx_min=cxmin
		if(cxmax.gt.cx_max) cx_max=cxmax
	enddo
	
	if(iflag .eq. 1)then
		do is=ns_start, nshot,	1
			write(*,*)	'index_table(1,', is,	'is)=',	index_table(1, is)
			write(*,*)	'index_table(2,', is,	'ntrace)=',		index_table(2, is)
			write(*,*)	'index_table(3,', is,	'csx)=',	index_table(3, is)
			write(*,*)	'index_table(4,', is,	'0)=',	index_table(4, is)
			write(*,*)	'index_table(5,', is,	'cx_min)=',	index_table(5, is)
			write(*,*)	'index_table(6,', is,	'cx_max)=',	index_table(6, is)
			write(*,*)	'index_table(7,', is,	'0)=',	index_table(7, is)
			write(*,*)	'index_table(8,', is,	'0)=',	index_table(8, is)
			write(*,*)	'index_table(9,', is,	'krec_save)=',	index_table(9, is)
		enddo
	endif

    return

    end	subroutine

!***********************************************************************
	subroutine Form_Dynamic_Filename(myid, fn4, fn4_1, ns_start)

	implicit none
	!Data dictionary: declare variable types, definitons, & units
	!Dummy Variables
	character(len=256)	::	fn4
	character(len=128)  ::	fn4_1
	integer	::	myid
	integer ::	ns_start
	!Local Variables
	character(len=5)	::	fn4_2
	character(len=3)	::	fn4_3
	character(len=5)	::	fn4_4
	character(len=1)	::	ns1
	character(len=2)	::	ns2
	character(len=3)	::	ns3
	character(len=4)	::	ns4
	character(len=5)	::	ns5
	character(len=1)	::	no1
	character(len=2)	::	no2
	character(len=3)	::	no3
	character(len=4)	::	no4
	character(len=5)	::	no5

	fn4_3='_id'
	if(myid.ge.1.and.myid.le.9)then
		write(no1,'(i1)') myid
		fn4_4='000'//no1
		if(ns_start.ge.1.and.ns_start.le.9)then
			write(ns1,'(i1)')ns_start
			fn4_2='0000'//ns1
			fn4=trim(fn4_1)//fn4_2//fn4_3//fn4_4
		!	write(*,*)'filename=',fn4,'   myid=',myid
		elseif(ns_start.ge.10.and.ns_start.le.99)then
			write(ns2,'(i2)')ns_start
			fn4_2='000'//ns2
			fn4=trim(fn4_1)//fn4_2//fn4_3//fn4_4
		!	write(*,*)'filename=',fn4,'   myid=',myid
        elseif(ns_start.ge.100.and.ns_start.le.999)then
			write(ns3,'(i3)')ns_start
			fn4_2='00'//ns3
			fn4=trim(fn4_1)//fn4_2//fn4_3//fn4_4
		!	write(*,*)'filename=',fn4,'   myid=',myid
        elseif(ns_start.ge.1000.and.ns_start.le.9999)then
			write(ns4,'(i4)')ns_start
			fn4_2='0'//ns4
			fn4=trim(fn4_1)//fn4_2//fn4_3//fn4_4
		!	write(*,*)'filename=',fn4,'   myid=',myid
        elseif(ns_start.ge.10000.and.ns_start.le.99999)then
			write(ns5,'(i5)')ns_start
			fn4_2=ns5
			fn4=trim(fn4_1)//fn4_2//fn4_3//fn4_4
		!	write(*,*)'filename=',fn4,'   myid=',myid
        else
			write(*,*)'form_dynamic_filename error'
			stop
        endif
	elseif(myid.ge.10.and.myid.le.99)then
        write(no2,'(i2)') myid
        fn4_4='00'//no2
        if(ns_start.ge.1.and.ns_start.le.9)then
			write(ns1,'(i1)')ns_start
			fn4_2='0000'//ns1
			fn4=trim(fn4_1)//fn4_2//fn4_3//fn4_4
		!	write(*,*)'filename=',fn4,'   myid=',myid
        elseif(ns_start.ge.10.and.ns_start.le.99)then
			write(ns2,'(i2)')ns_start
			fn4_2='000'//ns2
			fn4=trim(fn4_1)//fn4_2//fn4_3//fn4_4
		!	write(*,*)'filename=',fn4,'   myid=',myid
        elseif(ns_start.ge.100.and.ns_start.le.999)then
			write(ns3,'(i3)')ns_start
			fn4_2='00'//ns3
			fn4=trim(fn4_1)//fn4_2//fn4_3//fn4_4
		!	write(*,*)'filename=',fn4,'   myid=',myid
        elseif(ns_start.ge.1000.and.ns_start.le.9999)then
			write(ns4,'(i4)')ns_start
			fn4_2='0'//ns4
			fn4=trim(fn4_1)//fn4_2//fn4_3//fn4_4
		!	write(*,*)'filename=',fn4,'   myid=',myid
        elseif(ns_start.ge.10000.and.ns_start.le.99999)then
			write(ns5,'(i5)')ns_start
			fn4_2=ns5
			fn4=trim(fn4_1)//fn4_2//fn4_3//fn4_4
		!	write(*,*)'filename=',fn4,'   myid=',myid
        else
			write(*,*)'form_dynamic_filename error'
			stop
        endif
    elseif(myid.ge.100.and.myid.le.999)then
		write(no3,'(i3)') myid
		fn4_4='0'//no3
        if(ns_start.ge.1.and.ns_start.le.9)then
			write(ns1,'(i1)')ns_start
			fn4_2='0000'//ns1
			fn4=trim(fn4_1)//fn4_2//fn4_3//fn4_4
		!	write(*,*)'filename=',fn4,'   myid=',myid
        elseif(ns_start.ge.10.and.ns_start.le.99)then
			write(ns2,'(i2)')ns_start
			fn4_2='000'//ns2
			fn4=trim(fn4_1)//fn4_2//fn4_3//fn4_4
		!	write(*,*)'filename=',fn4,'   myid=',myid
        elseif(ns_start.ge.100.and.ns_start.le.999)then
			write(ns3,'(i3)')ns_start
			fn4_2='00'//ns3
			fn4=trim(fn4_1)//fn4_2//fn4_3//fn4_4
		!	write(*,*)'filename=',fn4,'   myid=',myid
        elseif(ns_start.ge.1000.and.ns_start.le.9999)then
			write(ns4,'(i4)')ns_start
			fn4_2='0'//ns4
			fn4=trim(fn4_1)//fn4_2//fn4_3//fn4_4
		!	write(*,*)'filename=',fn4,'   myid=',myid
        elseif(ns_start.ge.10000.and.ns_start.le.99999)then
			write(ns5,'(i5)')ns_start
			fn4_2=ns5
			fn4=trim(fn4_1)//fn4_2//fn4_3//fn4_4
		!	write(*,*)'filename=',fn4,'   myid=',myid
        else
			write(*,*)'form_dynamic_filename error'
			stop
        endif
	elseif(myid.ge.1000.and.myid.le.9999)then
        write(no4,'(i4)') myid
        fn4_4=no4
        if(ns_start.ge.1.and.ns_start.le.9)then
			write(ns1,'(i1)')ns_start
			fn4_2='0000'//ns1
			fn4=trim(fn4_1)//fn4_2//fn4_3//fn4_4
		!	write(*,*)'filename=',fn4,'   myid=',myid
        elseif(ns_start.ge.10.and.ns_start.le.99)then
			write(ns2,'(i2)')ns_start
			fn4_2='000'//ns2
			fn4=trim(fn4_1)//fn4_2//fn4_3//fn4_4
		!	write(*,*)'filename=',fn4,'   myid=',myid
        elseif(ns_start.ge.100.and.ns_start.le.999)then
			write(ns3,'(i3)')ns_start
			fn4_2='00'//ns3
			fn4=trim(fn4_1)//fn4_2//fn4_3//fn4_4
		!	write(*,*)'filename=',fn4,'   myid=',myid
        elseif(ns_start.ge.1000.and.ns_start.le.9999)then
			write(ns4,'(i4)')ns_start
			fn4_2='0'//ns4
			fn4=trim(fn4_1)//fn4_2//fn4_3//fn4_4
		!	write(*,*)'filename=',fn4,'   myid=',myid
        elseif(ns_start.ge.10000.and.ns_start.le.99999)then
			write(ns5,'(i5)')ns_start
			fn4_2=ns5
			fn4=trim(fn4_1)//fn4_2//fn4_3//fn4_4
		!	write(*,*)'filename=',fn4,'   myid=',myid
        else
			write(*,*)'form_dynamic_filename error'
			stop
        endif
	else
        write(*,*)'form_dynamic_filename error'
        stop
    endif

    return
    end	subroutine

!***********************************************************************
	subroutine Input_Shot_Gather(trace, trace_real, shot_gather, &
					nx, dx, lt, dt, lt_real, dt_real, &
					ns_x, cx_min, krec_save, &
					max_offset_x, min_offset_x)

!*==========================================================================
!*
!*     ns_x: the shot point position corresponding to the current shot and 
!*                 receiver range
!*     lt_real, dt_real: the sample point number and the sample rate from the 
!*                       inputed gather
!*
!*==========================================================================
	use global
	use header_module	
	implicit none
	!Data dictionary: declare variable types, definitons, & units
	!Dummy Variables
    integer ::	nx
	integer	::	lt, lt_real
	integer ::	ns_x
	integer ::	krec_save
	integer ::	cx_min
	real	::  max_offset_x
	real	::  min_offset_x
    real    ::	dx
	real	::	dt, dt_real

    real    ::	trace(lt)
	real	::	trace_real(lt_real)
    real    ::	shot_gather(lt, nx)
	!Local Variables
	type(segy)  ::  head
	integer	::	irec
	integer ::	ierr
	integer ::	it
	integer	::	lflag
    integer ::	nshot
	integer ::	ntrace_count_shot
	integer ::	npshot
	integer ::	interp_flag
	integer ::	mgx
    integer ::	sx, gx
	integer ::	num_trace
	integer ::	iflag
	

	iflag=0		!the flag of controlling the print  on-off
	!*----------------------------------------------------------------------*
    irec=krec_save
    lflag=1

1001  continue
    read(11, rec=irec, err=5555)head
	nshot=head%fldr

	!*-lflag=1 indicates that the first trace of the current shot is reading
	!into the working buffer
    if(lflag.eq.1) then
		ntrace_count_shot=1       ! trace number counter in a shot gather
        npshot=nshot
        if(modulo(dt_real, dt).eq.0.0.and.lt.gt.lt_real) then  !yzhou
			interp_flag = 2
        else
			interp_flag = 1
        end if
        lflag=2
    end if

    if(npshot.ne.nshot) then
	!*------ the current shot gather inputting finished !!
		goto 5555
    end if

    read(11, rec=irec)head,trace_real
	sx=head%sx
	gx=head%gx

	!*********************************************************************!
	!判断偏移距大小是否在预设的范围内
	if(abs(sx-gx) .gt. max_offset_x .or. abs(sx-gx) .lt. min_offset_x)then
		irec=irec+1
		goto 1001
	endif

	!*********************************************************************!

    if(interp_flag.eq.2) then
		CALL Trace_Sinc_Interp(trace, lt, dt, trace_real, lt_real, dt_real)
    else if(interp_flag.eq.1) then
		do it=1, lt
			trace(it) = trace_real(it)
		end do
    end if

	!*---------- the shot position of current shot -------------------------*
    mgx=(gx-cx_min)/(dx)+1.5
    ns_x=(sx-cx_min)/(dx)+1.5

    do it=1, lt
		shot_gather(it, mgx)=trace(it)
    enddo

	!*----------------------------------------------------------------------*
    irec=irec+1
    ntrace_count_shot=ntrace_count_shot+1
    goto 1001

5555  continue

    krec_save=irec
    num_trace=ntrace_count_shot-1

	return
    end	subroutine

!***********************************************************************
!*======================================================================*
!*     for the calculation stability, the temporal sample rate must be 
!*     dense enough. sometimes, the interpolation is necessary.
!*     here, the sinc interpolation is used.
!*     2010.2.11 modified by whz 
!*======================================================================*

!***********************************************************************
    subroutine Trace_Sinc_Interp(trace, lt, dt, trace_real, lt_real, dt_real)

	implicit none
	!Data dictionary: declare variable types, definitons, & units
	!Dummy Variables
	integer ::	lt, lt_real
	real	::	dt, dt_real
	real	::	trace(lt), trace_real(lt_real)

	!Local Variables
	!Data dictionary: declare constants
	integer,parameter	::	lw=4
	real,parameter	::	pi=3.14159265359
	!Data dictionary: declare variable types, definitons, & units
	integer ::	it, iit
	integer ::	iw
	integer	::	n_interped
	integer ::	it_start, it_final
	real	::	rdt
	real	::	tmp_hf
	real	::	tdif
	real	::	x
	real	::	sinc

    n_interped =nint(dt_real/dt)
    rdt=1.0/dt_real

    do it = 1, lt

		!iit=int(it/n_interped)
        !write(*,*)'iit-int'
        !write(*,*)iit
        iit = it/n_interped + 0.5
        !write(*,*)'itt-+0.5'
        !write(*,*)iit

        it_start = iit - lw
        if(it_start.lt.1) it_start = 1

        it_final = iit + lw
        if(it_final.gt.lt_real) it_final = lt_real

        tmp_hf = 0.0
        tdif = it * dt - it_start * dt_real

        do iw = it_start, it_final
			if(tdif.eq.0.0) then
				tmp_hf = tmp_hf + trace_real(iw)
			else
				x = pi * tdif * rdt
				sinc = sin(x)/x
				tmp_hf = tmp_hf + trace_real(iw)*sinc
			endif
			tdif = tdif - dt_real

		enddo
		trace(it)  = tmp_hf
	enddo		

    return
    end subroutine

!***********************************************************************
	subroutine Wavelet_Forming(wavelet, dt, lt, fmain, nwt)
!*==========================================================================
!*  assigning the source function for the shot wavefield extrapolation
!*  fmain: the main frequency of the source function
!*  here the ricker wavelet is used.
!*  2010,2,11, modified by herb.wang
!*  richer wavlet function: f(t) = (1.0-2.0*(pai*fmain*t)**2)exp(-(pai*fmain*t)**2)
!*==========================================================================

	implicit none
	!Dummy Variables
	!Data dictionary: declare variable types, definitons, & units
	integer ::	lt
	integer ::	nwt
	real	::	dt
	real	::	fmain
	real	::	wavelet(lt)

	!Local Variables
	!Data dictionary: declare constants
	real,parameter	::	pi=3.14159265359
	!Data dictionary: declare variable types, definitons, & units
	integer it
	real	tmain
	real	tp1, tp2
	real	time

	!tmain=1.2/fmain
    tmain=1.0/fmain
    nwt=tmain/dt-1
		
    do it=1, lt
		!time = it*sdt-100*sdt  !zhaolei_20100720
		time = (it-1)*dt-tmain
		!time = (it-1)*sdt!zero phase wavelet hephaestus_9_25

		tp1 = pi*fmain*time
		tp2 = tp1*tp1
		wavelet(it) = (1.0-2.0*tp2)*exp(-tp2)!*0.1e7
    enddo

    return
    end subroutine

!***********************************************************************
    subroutine Form_Coe(coe, dx, dz, dt)

	implicit none
	!Data dictionary: declare variable types, definitons, & units
	!Dummy Variables
    real	::	coe(13)
    real	::	dx, dz, dt
	!Local Variables
	real	::	ddt
	real	::	dtx2
	real	::	dtz2
	real	::	omga0, omga1, omga2, omga3, omga4, omga5

	ddt=dt
      
    dtx2=ddt*ddt/dx/dx/2.0
    dtz2=ddt*ddt/dz/dz/2.0

    omga0=-5.8544444444
    omga1=+3.3333333333
    omga2=-0.4761904762
    omga3=+0.0793650794
    omga4=-0.0099206349
    omga5=+0.0006349206

    coe(1) =+2
    coe(2) =-1
    coe(3) =omga0*(dtx2+dtz2)
    coe(4) =omga1*dtx2
    coe(5) =omga2*dtx2
    coe(6) =omga3*dtx2
    coe(7) =omga4*dtx2
    coe(8) =omga5*dtx2
    coe(9) =omga1*dtz2
    coe(10)=omga2*dtz2
    coe(11)=omga3*dtz2
    coe(12)=omga4*dtz2
    coe(13)=omga5*dtz2

    return
    end	subroutine

!***********************************************************************
    subroutine ReadPar(fn_par, fn1, fn2, fn3, &
					fn_single, fn_dynamic, &
					ns_start, ns_interval, ns_end, &
					lt, lt_real, &
					dt, dt_real, &
					depth, &
					dx, dz, &
					aperture_x_l, aperture_x_r, &
					boundary_x_l, boundary_x_r, &
					boundary_z_u, boundary_z_d, &
					Cvx_initial, &
					nvx, nvz, & 
					dvx, dvz, &
					fmain, &
					max_offset_x, min_offset_x, &
					shot_direct, order)

!*======================================================================
!*  Variable demonstrations:
!*
!*		fn1:			the shot gather file name
!*		fn2:			the velocity file name
!*		fn3:			the image profile
!*		fn_single:		the single gradient profile file
!*		fn_dynamic:		the dynamic temp file for  storage gradient result 
!*
!*		ns_start   :	the start shot number 
!*		ns_end:			the final shot number 
!*		ns_interval:	the step length of shot number for reverse time shot migration
!*    
!*		lt:				the exptropalation sample point number of the trace length
!*		lt_real:		the actually sample point number of the trace length
!*		dt:				the exptropalation sample rate
!*		dt_real:		the actually sample rate
!*
!*		depth:			the exptropalation/imaging depth 
!*		dz:				the exptropalation step length
!*		dx:				the interval of imaging grid
!*
!*		aperture_x:		the cross_line and in_line migration aperture
!*		boundary_x :	the width of the absorbing boundary zone along x_direction
!*		boundary_z :	the width of the absorbing boundary zone along z_direction
!*
!*		cvx_initial:	the initial coodinate of the velocity model
!*		nvx, nvz:		the dimension of the velocity model
!*		dvx, dvz:		the sample rate of the velocity model    

!*		fmain:			the main frequency of the source function 
!*		max_offset_x:	the maxinum of offset
!*		min_offset_x:	the mininun of offset
!*		shot_direct:	whether direct wave exist(shot_direct=1 represent
!					existance and 0 vice versa
!*
!*======================================================================
	implicit none
	!Data dictionary: declare variable types, definitons, & units
	!Dummy Variables
	character(len=256)	::	fn_par, fn1, fn2, fn3, fn_single
	character(len=128)  ::	fn_dynamic
	integer ::	ierr
	integer ::	ns_start, ns_interval, ns_end
	integer ::	lt, lt_real
	real	::	dt, dt_real
	real	::	depth, dz
	real	::	dx
	real	::	aperture_x_l, aperture_x_r
	real	::	boundary_x_l, boundary_x_r
	real	::	boundary_z_u, boundary_z_d
	real	::	cvx_initial
	integer ::	nvx, nvz
	real	::	dvx, dvz
	real	::	fmain
	real	::	max_offset_x
	real	::	min_offset_x
	integer ::	shot_direct
	integer	::	order

	character(len=256)	::	par_name

	open(33, file=fn_par, status='old', iostat=ierr)
	if(ierr /= 0)then
		write(*,*)	'Parameter file cannot open right, please check it !!!'
		write(*,*)	'Shutting down the program'
		stop
	endif

    !*fn1	the shot gather filename
	read(33, '(a)')	par_name
    read(33, '(a)') fn1

    !*fn2	the velocity model filename'
	read(33, '(a)')	par_name
    read(33, '(a)') fn2

	!*fn3: the image profile
	read(33, '(a)')	par_name
    read(33, '(a)') fn3

	!*fn_single:   the single gradient profile file
	read(33, '(a)')	par_name
    read(33, '(a)') fn_single

	!*fn_dynamic	the dynamic temp file for  storage gradient result
	read(33, '(a)')	par_name
	read(33, '(a)') fn_dynamic

	!*	the first, interval and last shot
	read(33, '(a)')	par_name
    read(33, *)	ns_start, ns_interval, ns_end

	!*	the imaging grid interval
	read(33, '(a)')	par_name
    read(33, *)	dx

	!*	the migration aperture
	read(33, '(a)')	par_name
    read(33, *)	aperture_x_l, aperture_x_r

	!*	'the trace length (sample point number) and dt(ms)
	read(33, '(a)')	par_name
    read(33, *)	lt, dt, lt_real, dt_real

	!*	the imaging depth(m) and step_length(m)
	read(33, '(a)')	par_name
    read(33, *)	depth, dz

	!*	the velocity starting coordinate (x)(m)
	read(33, '(a)')	par_name
    read(33, *)	cvx_initial

	!*	the width of absorbing boundary of x axis
	read(33, '(a)')	par_name
    read(33, *)	boundary_x_l, boundary_x_r

	!*	the width of absorbing boundary of z axis
	read(33, '(a)')	par_name
    read(33, *)	boundary_z_u, boundary_z_d

	!*	the dimension of the velocity model
	read(33, '(a)')	par_name
    read(33, *)	nvx, nvz

	!*	the sampling rate of the velocity model(m)
	read(33, '(a)')	par_name
    read(33, *)	dvx, dvz

	!*	the main frequency(hz)
	read(33, '(a)')	par_name
    read(33, *)	fmain

	!*	the maxinum of offset
	read(33, '(a)')	par_name
	read(33, *) max_offset_x

	!*	the mininun of offset
	read(33, '(a)')	par_name
	read(33, *)	min_offset_x

	!*	whether direct wave exist(shot_direct=1 represent existance
	!*and 0 vice versa)
	read(33, '(a)')	par_name
	read(33, *)	shot_direct
	
	!*反传波场记录的宽度
	read(33, '(a)')	par_name
	read(33, *)	order

    close(33)

    return
    end	subroutine

!***********************************************************************

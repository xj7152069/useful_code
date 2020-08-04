

!***********************************************************************
	  module	model_2D_Multi_SHOT_Parameter
	  use	global	
	  implicit none

	  	!==========================================================
	  	!requisite variable for rtm
	  	!read parameter
	  	character(len=para_char_flen)   :: fn_shot_info

		character(len=para_char_flen)	::	fn_cs, fn_vel, fn_mig
		character(len=para_char_flen)  	::  fn_single
		character(len=para_char_flen)  	::  fn_dynamic

		integer	::	ns_start, ns_interval, ns_end
		integer	::	lt, lt_real
		real	::	dt, dt_real
		real	::	depth
		real	::	dx, dz
		real	::	aperture_x_l, aperture_x_r
		real	::	aperture_z_u, aperture_z_d
		real	::	boundary_x_l, boundary_x_r
		real	::	boundary_z_u, boundary_z_d
		real	::	cvx_initial, cvz_initial
		integer	::	nvx, nvz
		real	::	dvx, dvz
		real	::	fmain
	  	integer	::  flag_shot_info
	  	integer	::	flag_output_sig

	  	!==========================================================
	  	!requisite variable for multi shot rtm
		character(len=para_char_flen)	::	fn_dynamic_all
		integer	::	ns_total
	  	integer	::	nshot
	  	integer	::	ntrace_max
	  	integer	::	ntrace


		!index_table info
		integer*8 ::	krec_save
	  	integer	::	ishot
		real	::	cx_min
		real	::	cx_max
		real	::	cz_min

		integer	::	cx_total_min
	  	integer	::	cx_total_max


		integer	::	nz, nx
		integer	::	nwt
	  	integer	::	ns
		integer ::	nx_apert_l, nx_apert_r
		integer ::	nz_apert_u, nz_apert_d
		integer ::	nx_bound_l, nx_bound_r
		integer ::	nz_bound_u, nz_bound_d
	  	integer	::	nx_with_apert, nz_with_apert
	  	integer	::	nnx, nnz
		integer ::	nvxx_shift, nvzz_shift
		integer	::	nx_shift, nz_shift
		integer	::	ns_x, ns_z, nr_x, nr_z

	  	integer	::	it_step

		real	::	vmax
	  	real	::	vmin
	  	real	::	mem_need

	  	!==========================================================
	  	!memory parameter
	  	real,parameter	::	mem_g	=	1024.0*1024.0*1024.0
	  	real,parameter	::	memmax	=	50.0
	  	integer,parameter	::	npad	=	5
	  	integer,parameter	::	order	=	5


	  	!==========================================================
	  	!info 
	  	integer,parameter	::	iflag_debug=1
	  		!=0		no debug  info
	  		!=1		print debug info
	  		!=2		some temp debug info
		character(len=para_char_flen)  ::  currt, currtsnap
		character(len=para_char_flen),parameter	::	fn_debug= './log/model2d_debug_'
	  	integer,parameter	::	debug_ishot = 1

	  	integer,parameter   ::  iflag_src = 1


		integer ::	ierr
		integer	::	isok

	  end module	model_2D_Multi_SHOT_Parameter

!***********************************************************************



!***********************************************************************
	subroutine	model_2d_mpi(para_int, para_long, &
					para_float, para_double, para_char, &
					mpi_np, myid )

	use global
	use shotgather_info
	use header_module
	use model_2D_Multi_SHOT_Parameter 
	implicit none
	include 'mpif.h'

	!Data dictionary: declare variable types, definitons, & units
	!Dummy Variables
    integer status(mpi_status_size)
	integer     ::  para_int(para_int_fnum)
	integer*8   ::  para_long(para_long_fnum)
	real        ::  para_float(para_float_fnum)
	real*8      ::  para_double(para_double_fnum)
	character(len=para_char_flen)	::	para_char(para_char_fnum)
	integer		::	mpi_np
	integer		::	myid

	!Local variables
	integer     ::  sig_para_int(para_int_fnum)
	integer*8   ::  sig_para_long(para_long_fnum)
	real        ::  sig_para_float(para_float_fnum)
	real*8      ::  sig_para_double(para_double_fnum)
	character(len=para_char_flen)	::	sig_para_char(para_char_fnum)

    integer*8	::	l1, l2, l3, l4, l5, l6, l7, l8, l9, l10
	integer*8	::	l11, l12, l13, l14, l15, l16, l17, l18, l19, l20
	integer*8	::	l21, l22, l23, l24, l25, l26, l27, l28, l29, l30
	integer*8	::	l31, l32, l33, l34, l35, l36, l37, l38, l39, l40
	integer*8	::	l41, l42, l43, l44, l45, l46, l47, l48, l49, l50
	integer*8	::	l51, l52, l53, l54, l55, l56, l57, l58, l59, l60
	integer*8	::	l61, l62, l63, l64, l65, l66, l67, l68, l69, l70
	integer*8	::	ltotal
	

	type(shotinfo),allocatable	::	shotreceiver(:)
	real,allocatable	::	vel(:,:), vv(:,:)
	real,allocatable	::	shot_gather(:,:)
	real,allocatable	::	wavelet(:)
	real,allocatable	::	gx(:), gz(:)
	integer,allocatable	::	index_gx(:), index_gz(:)
	real	::	sx, sz

	integer	::	master
	integer	::	i, ii
	integer ::	isend
	integer	::	ix, iz
	integer	::	iix, iiz
	integer	::	ivx, ivz
	integer	::	it, itr, itrace
	integer	::	shotnum
	real	::	xx, yy, zz

	real	::	time1, time2, time

	!处理接口函数参数传递
	master=0.0
	!===========================================================================
	nvx = para_int(1)
	nvz = para_int(2)
	allocate(vel(nvz, nvx), stat=ierr)
	vel=0.0
		!打开速度文件
		!采用初始的速度
		fn_vel = para_char(1)
		open(21, file=fn_vel, access='direct', status='old', recl=lbyte*nvz*nvx, iostat=ierr)
		if(ierr /= 0)then
			write(*,*)  'Velocity file cannot open right, please check it!!!'
			write(*,*)  'Velocity file ',fn_vel
			stop
		endif
		read(21, rec=1)((vel(ivz, ivx), ivz=1, nvz), ivx=1, nvx)
		close(21)

!	vel=2000.00

	vmax=maxval(vel)
	vmin=minval(vel)
	para_float(34)	= vmin
	para_float(35)	= vmax

	!===========================================================================
	CALL model_2D_MULTI_SHOT_Initialize_Checkset(para_int, para_long, para_float, &
			para_double, para_char, myid)

	!===========================================================================
	!炮检点坐标信息
	allocate(shotreceiver(ns_total))
    do ishot =1, ns_total
      		allocate(shotreceiver(ishot)%gx(ntrace_max))
      		allocate(shotreceiver(ishot)%gy(ntrace_max))
      		allocate(shotreceiver(ishot)%gz(ntrace_max))
	enddo


	if(flag_shot_info .eq. 1)then

	 	open(21, file=fn_shot_info, status='old')
	 	ishot=0
	 	do while (.TRUE.)
			read(21, *, iostat=ierr)ii, xx, yy, zz
         	if(ierr /= 0)then
        		exit
			endif
            if(ii .ne. -1)then
            	ishot = ishot +1
				if(ishot .gt. 1)then
					shotreceiver(ishot-1)%ntr = itrace
				endif
				if(ishot .gt. ns_total)	exit
				if(ishot .gt. 1)then
					shotreceiver(ishot)%krec = shotreceiver(ishot-1)%krec + itrace
				endif
				if(ishot  .eq. 1)then
					shotreceiver(ishot)%krec = 1
				endif
       			shotreceiver(ishot)%fldr = ii
      			shotreceiver(ishot)%sx = xx
      			shotreceiver(ishot)%sy = yy
      			shotreceiver(ishot)%sz = zz
       			itrace = 0

			else if(ii  .eq. -1)then
       			itrace = itrace +1
				if(itrace .gt. ntrace_max)then
					print*, 'Read the trace error!!!!'
					stop
				endif
      			shotreceiver(ishot)%gx(itrace)	=	xx
      			shotreceiver(ishot)%gy(itrace)	=	yy
      			shotreceiver(ishot)%gz(itrace)	=	zz

			endif
		enddo
		shotreceiver(ishot).ntr = itrace
		close(21)

	endif

	if(myid .eq. 0)then
		print*, 'ns_total, ishot', ns_total, ishot
		if(ishot .lt. ns_total)then
			print*, 'The read shot info error!!!!'
			stop
		endif

		print*, 'The shot and receiver info :      '
		do ishot=1, ns_total
			print*, 'fldr', ishot, shotreceiver(ishot)%fldr
			print*, 'krec', ishot, shotreceiver(ishot)%krec
			print*, 'sx  ', ishot, shotreceiver(ishot)%sx
			print*, 'sx  ', ishot, shotreceiver(ishot)%sy
			print*, 'sz  ', ishot, shotreceiver(ishot)%sz
			print*, 'ntr ', ishot, shotreceiver(ishot)%ntr
			print*, '      '
			print*, '      '
		enddo

	endif

	if(flag_shot_info  .eq. 0)then

	endif
	!========================================================================

	!========================= MPI Enviroment Initialed ======================
    CALL MPI_BCAST(ns_total, 1, MPI_INTEGER, MASTER, MPI_COMM_WORLD, IERR)
    CALL MPI_BARRIER(MPI_COMM_WORLD, IERR)

	!*=======================================================================*

    if( myid .eq. master ) then
		do i=1, ns_total
			CALL MPI_RECV(shotnum, 1, MPI_INTEGER, MPI_ANY_SOURCE,&
					MPI_ANY_TAG, MPI_COMM_WORLD, STATUS, IERR)
			isend=status(mpi_source)
			CALL MPI_SEND(i, 1, MPI_INTEGER, ISEND, I,&
					MPI_COMM_WORLD, IERR)
        enddo
        do i = 1, mpi_np-1
			CALL MPI_RECV(shotnum, 1, MPI_INTEGER, MPI_ANY_SOURCE,&
					MPI_ANY_TAG, MPI_COMM_WORLD, STATUS, IERR)
			isend=status(mpi_source)
			shotnum=0
			CALL MPI_SEND(shotnum, 1, MPI_INTEGER, ISEND, I,&
					MPI_COMM_WORLD, IERR)
        enddo
    else
		!*--------------------worker-process------------------------------------*
		shotnum=0
        CALL MPI_SEND(shotnum, 1, MPI_INTEGER,0, 0, MPI_COMM_WORLD, IERR)

90      continue
        CALL MPI_RECV(shotnum, 1, MPI_INTEGER, 0, MPI_ANY_TAG, MPI_COMM_WORLD, STATUS, IERR)

        if( shotnum .ne. 0 ) then
			!index_table info
			ishot    =  shotnum
			ntrace	 = shotreceiver(ishot)%ntr

			if(ntrace .lt. 1)then
				write(*,*)'The INput NTRACE IS ERROR'
			endif

			!===============================================================
			!读取当前炮数据，并记录当前炮的炮检点坐标
			allocate(gx(ntrace))
			allocate(gz(ntrace))
			allocate(index_gx(ntrace))
			allocate(index_gz(ntrace))

			cx_min = 999999.9
			cx_max = -99999.9
			cz_min = 999999.9

			sx = shotreceiver(ishot)%sx
			sz = shotreceiver(ishot)%sz

			do itrace=1, ntrace
				gx(itrace) = shotreceiver(ishot)%gx(itrace)
				gz(itrace) = shotreceiver(ishot)%gz(itrace)

				if(shotreceiver(ishot)%gx(itrace) .gt. cx_max)then
					cx_max	= shotreceiver(ishot)%gx(itrace)
				endif
				if(shotreceiver(ishot)%gx(itrace) .lt. cx_min)then
					cx_min	= shotreceiver(ishot)%gx(itrace)
				endif
				if(shotreceiver(ishot)%gz(itrace) .lt. cz_min)then
					cz_min	= shotreceiver(ishot)%gz(itrace)
				endif
			enddo

			cx_min = min(sx, cx_min)
			cx_max = max(sx, cx_max)

			print*, 'cx_min, cx_max, ', cx_min, cx_max

			!*-------- input the current shot gather data into the working bufer====*
			nx=(cx_max-cx_min)/dx+1
			nz=nvz

			!*the imaging range after adding the migration aperture
			nx_with_apert = nx_apert_l + nx + nx_apert_r
			nz_with_apert = nz_apert_u + nz + nz_apert_d

			!*the wave extrapolation range after adding the absorbing boundary ---*
			nnx = nx_bound_l +  nx_with_apert + nx_bound_r
			nnz = nz_bound_u +  nz_with_apert + nz_bound_d

			!*the shift between the current shot imaging range and the total velocity field range
			nvxx_shift = (cx_min-cvx_initial)/dvx - nx_apert_l
			nvzz_shift = (cz_min-cvz_initial)/dvz - nz_apert_u
			nx_shift = nx_apert_l + nx_bound_l
			nz_shift = nz_apert_u + nz_bound_u

			!========================================================
			!allocate some variables for single shot rtm 
			allocate(shot_gather(lt, ntrace))
			allocate(vv(nnz, nnx))
			allocate(wavelet(lt))

			!*input the velocity field for the current shot extrapolation-----*
    		call Get_Current_Shot_Velocity_2d(vel, vv, &
					nvz, nvx, nnz, nnx, &
					nz_bound_u, nx_bound_l, &
					nvzz_shift, nvxx_shift, &
					isok)


			call Source_and_Receiver_Position_2d(sx, sz, gx, gz, &
					ns_x, ns_z, index_gx, index_gz, &
					ntrace, cx_min, cz_min, &
					dx, dz, &
					nx_shift, nz_shift, &
					isok)

			!===============================================================
			!*======= forming the source function for the forward extrpolation ------*
    		CALL Wavelet_Forming(wavelet, dt, lt, fmain, nwt)


			!================================================================
			!单炮RTM需要的参数
			sig_para_int(1)		=	ishot
			sig_para_int(2)		=	myid
			sig_para_int(3)		=	lt
			sig_para_int(4)		=	ntrace
	!		sig_para_int(5)		=	ns
	!		sig_para_int(6)		=	order
			sig_para_int(7)		=	nwt
			sig_para_int(8)		=	nx
			sig_para_int(9)		=	nz
			sig_para_int(10)	=	nx_apert_l
			sig_para_int(11)	=	nx_apert_r
			sig_para_int(12)	=	nz_apert_u
			sig_para_int(13)	=	nz_apert_d
			sig_para_int(14)	=	nx_bound_l
			sig_para_int(15)	=	nx_bound_r
			sig_para_int(16)	=	nz_bound_u
			sig_para_int(17)	=	nz_bound_d
			sig_para_int(18)	=	nx_with_apert	
			sig_para_int(19)	=	nz_with_apert
			sig_para_int(20)	=	nnx
			sig_para_int(21)	=	nnz
			sig_para_int(22)	=	nx_shift
			sig_para_int(23)	=	nz_shift
	!		sig_para_int(24)	=	it_step
			sig_para_int(25)	=	npad

			sig_para_float(1)	=	dt
			sig_para_float(2)	=	dx
			sig_para_float(3)	=	dz


			if(iflag_debug .eq. 2 .and. ishot .eq. debug_ishot)then
				write(currt, '(I8)')nnz
				open(122, file=trim(adjustl(fn_debug))//'vv_nnz'//trim(adjustl(currt))//'.dat', access='direct', &
							recl=lbyte*nnz*nnx, status='replace')
				write(122, rec=1)((vv(iiz, iix), iiz=1, nnz), iix=1, nnx)
				close(122)
			endif

	
			if(iflag_debug .eq. 1)then
				open(111, file=trim(fn_debug)//'sig_parameters_output.txt', status='replace')	
				write(111, *) 'The current shot ', ishot, '  Parameter is shown as: '
				write(111, *)'ishot', ishot
				write(111, *)'it_step = ', it_step
				write(111, *)'ntr = ', ntrace
				write(111, *)'lt = ', lt, ' ns = ', ns, 'order ', order
				write(111, *)'nwt = ', nwt
				write(111, *)'dx = ', dx, ' dz = ', dz
				write(111, *)'dt = ', dt
				write(111, *)'nx_apert_l = ', nx_apert_l, ' nx_apert_r = ', nx_apert_r
				write(111, *)'nz_apert_u = ', nz_apert_u, ' nz_apert_d = ', nz_apert_d
				write(111, *)'nx_bound_l = ', nx_bound_l, ' nx_bound_r = ', nx_bound_r
				write(111, *)'nz_bound_u = ', nz_bound_u, ' nz_bound_d = ', nz_bound_d
				write(111, *)'nx = ', nx,  ' nz = ', nz
				write(111, *)'nx_with_apert = ', nx_with_apert
				write(111, *)'nnx = ', nnx, 'nnz = ', nnz
				write(111, *)'nvxx_shift = ', nvxx_shift,  'nvzz_shift = ', nvzz_shift
				write(111, *)'nx_shift = ', nx_shift, ' nz_shift = ', nz_shift
				write(111, *)'sx, sz', sx ,sz
				write(111, *)'ns_x, ns_z', ns_x, ns_z
				close(111)
			endif

			if(iflag_debug  .lt. 2)then
				print*, '    RUN Single Shot RTM: Begin : ', ishot
			!	time1=MPI_WTIME()
			endif
		!	goto 1212
			!*=====================================================================*
			!* single shot reverse time extrapolation and getting image
    		call model_2D_single_shot(sig_para_int, sig_para_long, &
					sig_para_float, sig_para_double, sig_para_char, &
					wavelet, shot_gather, vv, &
					ns_x, ns_z, &
					index_gx, index_gz, &
					lt, ntrace, nnx, nnz, nx, nz )
!1212		continue

    		call model_2d_Write_Disk(fn_dynamic, ishot, shot_gather, lt, ntrace, dt, &
					sx, sz, gx, gz, &
					isok)


			if(iflag_debug  .lt. 2)then
				print*, '    RUN Single Shot RTM: End : ', ishot
			!	time2=MPI_WTIME()
			!	print*, 'Time is : ', (time2-time1)/60.00
			endif

			deallocate(wavelet)
			deallocate(shot_gather)
			deallocate(vv)
			deallocate(gx)
			deallocate(gz)
			deallocate(index_gx)
			deallocate(index_gz)

			CALL MPI_SEND(shotnum, 1, MPI_INTEGER, 0, 0, MPI_COMM_WORLD, IERR)
			goto 90
		endif
    endif
	!*=====================================================================*

	!归约结果
	CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )

	if(myid .eq. 0)then
		call merge_shot_files(fn_cs, fn_dynamic, ns_total, lt, shotreceiver, flag_output_sig)
	endif

	CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )


    do ishot =1, ns_total
      	deallocate(shotreceiver(ishot)%gx)
      	deallocate(shotreceiver(ishot)%gy)
      	deallocate(shotreceiver(ishot)%gz)
	enddo
	deallocate(shotreceiver)
	deallocate(vel)

	return
	end subroutine

!*************************************************************************
	subroutine	model_2D_MULTI_SHOT_Initialize_Checkset(para_int, para_long, para_float, para_double, para_char, myid)
	
	use global
	use model_2D_Multi_SHOT_Parameter 
	implicit none

	!Data dictionary: declare variable types, definitons, & units
	!Dummy Variables
	integer     ::  para_int(para_int_fnum)
	integer*8   ::  para_long(para_long_fnum)
	real        ::  para_float(para_float_fnum)
	real*8      ::  para_double(para_double_fnum)
	character(len=para_char_flen)	::	para_char(para_char_fnum)
	integer		::	myid

	!Local Variables
	integer	::	lt_extra
	real	::	dt_step, dt_extra
	real	::	v1, v2, vt
	!高阶差分稳定性
	real	::	sc(5)

	fn_vel			=	para_char(1)
	fn_cs			=	para_char(2)
	fn_dynamic		=	para_char(3)
	fn_shot_info	=	para_char(4)


	nvx				=	para_int(1)
	nvz				=	para_int(2)
	lt				=	para_int(3)
	ns_total		=	para_int(4)
	ntrace_max		=	para_int(5)
	flag_shot_info	=	para_int(6)
	flag_output_sig	=	para_int(7)



	fmain		=	para_float(1)
	dt			=	para_float(2)
	aperture_x_l=	para_float(3)	
	aperture_x_r=	para_float(4)
	aperture_z_u=	para_float(5)	
	aperture_z_d=	para_float(6)
	boundary_x_l=	para_float(7)
	boundary_x_r=	para_float(8)
	boundary_z_u=	para_float(9)
	boundary_z_d=	para_float(10)
	dx			=	para_float(11)
	dz			=	para_float(12)
	dvx			=	para_float(13)
	dvz			=	para_float(14)
	cvx_initial =	para_float(15)
	cvz_initial =	para_float(16)


	!*change the dt ms to s 
	dt = dt*0.001

    nx_apert_l = aperture_x_l/dx + 1
    nx_apert_r = aperture_x_r/dx + 1
    nz_apert_u = aperture_z_u/dz + 1
    nz_apert_d = aperture_z_d/dz + 1

    nx_bound_l = boundary_x_l/dx + 1
    nx_bound_r = boundary_x_r/dx + 1
    nz_bound_u = boundary_z_u/dz + 1
    nz_bound_d = boundary_z_d/dz + 1

	!========================================================
	!PML and apert  check 
	if(nx_apert_l .lt. npad )then
		print*, 'The Aperture x left is too small, increase up ', npad
		nx_apert_l = npad
	endif

	if(nx_apert_r .lt. npad )then
		print*, 'The Aperture x right is too small, increase up ', npad
		nx_apert_r = npad
	endif

	if(nz_apert_u .lt. npad )then
		print*, 'The Aperture z up is too small, increase up ', npad
		nz_apert_u = npad
	endif

	if(nz_apert_d .lt. npad )then
		print*, 'The Aperture z down is too small, increase up ', npad
		nz_apert_d = npad
	endif
	
	if(nx_bound_l .lt. npad )then
		print*, 'The PML Boundary x left is too small, increase up ', npad
		nx_bound_l = npad
	endif

	if(nx_bound_r .lt. npad )then
		print*, 'The PML Boundary x right is too small, increase up ', npad
		nx_bound_r = npad
	endif

	if(nz_bound_u .lt. npad )then
		print*, 'The PML Boundary z up is too small, increase up ', npad
		nz_bound_u = npad
	endif

	if(nz_bound_d .lt. npad )then
		print*, 'The PML Boundary z down is too small, increase up ', npad
		nz_bound_d = npad
	endif

	!========================================================
	vmin	=	para_float(34)
	vmax	=	para_float(35)
	
	if(vmin .lt. 100.0 .or. vmax .gt. 20000.0)then
		print*, 'THe velocity model is not right !!!! '
		print*, 'vmin = ', vmin, '  vmax =  ', vmax 
		stop
	endif

	!========================================================
	sc(1)	=	1.0
	sc(2)	=	0.866
	sc(3)	=	0.813
	sc(4)	=	0.784
	sc(5)	=	0.765

	!========================================================
	!高阶差分稳定性
	dt_extra=dt
	v1=vmax*sqrt(1.0/(dx*dx) + 1.0/(dz*dz))
	vt=dt*v1
	dt_step=1.0
	do while  (vt .gt. sc(5) )
		dt_step=dt_step*2.0
		dt_extra=dt_extra/dt_step
		vt=dt_extra*v1

		if(dt_step .gt. 10.0)then
			print*, 'The interpolation inverval too large. Please check the parameter.'	
			print*, 'Shut down...'
			stop
		endif

	enddo
	lt_extra = lt*dt/dt_extra+0.5

	lt = lt_extra
	dt = dt_extra

	!========================================================
	if(ns_total .lt. 1)then
		write(*,*)	'THe shot numbers is not right, please check it !!!'
		write(*,*)	'Shutting down the program'
		stop
	endif


	if(iflag_debug .eq. 1)then
		if(myid .eq. 0)then
			open(111, file=trim(fn_debug)//'model_2d_multi_shot_parameters_output.txt', status='replace')	
			write(111, *)"fn_cs : ", trim(fn_cs)
			write(111, *)"fn_dynamic : ", trim(fn_dynamic)
			write(111, *)'fmain  = ', fmain
			write(111, *)'lt = ', lt
			write(111, *)'dt  = ',  dt
			write(111, *)'nvx, nvz = ', nvx, nvz
			write(111, *)'dvx, dvz = ', dvx, dvz
			write(111, *)'aperture_x_l = ', aperture_x_l 
			write(111, *)'aperture_x_r = ', aperture_x_r
			write(111, *)'aperture_z_u = ', aperture_z_u
			write(111, *)'aperture_z_d = ', aperture_z_d
			write(111, *)'boundary_x_l = ', boundary_x_l
			write(111, *)'boundary_x_r = ', boundary_x_r
			write(111, *)'boundary_z_u = ', boundary_z_u
			write(111, *)'boundary_z_d = ', boundary_z_d
			write(111, *)'dx, dz = ', dx, dz
			write(111, *)'cvx_initial = ', cvx_initial
			write(111, *)'cvz_initial = ', cvz_initial
			write(111, *)'ns_total = ', ns_total
			write(111, *)'nx_apert_l, nx_apert_r ', nx_apert_l, nx_apert_r
			write(111, *)'nz_apert_u, nz_apert_d ', nz_apert_u, nz_apert_d
			write(111, *)'nx_bound_l, nx_bound_r ', nx_bound_l, nx_bound_r
			write(111, *)'nz_bound_u, nz_bound_d ', nz_bound_u, nz_bound_d
			write(111, *)'vmin, vmax', vmin, vmax

			close(111)
		endif
	endif


	isok=0

	end subroutine

!***********************************************************************

    subroutine model_2d_Write_Disk(fn1, ishot, shot_gather, lt, ntrace, dt, &
					sx, sz, gx, gz, &
					isok)

	use global
	use header_module	
	implicit none
	!Data dictionary: declare variable types, definitons, & units
	!Dummy Variables
	character(len=para_char_flen)	::	fn1
	real	::	shot_gather(lt, ntrace)
	real	::	gx(ntrace), gz(ntrace)
	real	::	sx, sz
	integer	::	ishot
	integer	::	lt, ntrace
	real	::	dt
	integer	::	isok

	!Local Variables
	character(len=para_char_flen)	::	currt
	type(segy)::head
	integer	::	it, itr


	!Initiallization of the su_head
	head=segy(0,0,0,0,0,0,0,0,0,0,0,&
								0,0,0,0,0,0,0,0,0,0,&
								0,0,0,0,0,0,0,0,0,0,&
								0,0,0,0,0,0,0,0,0,0,&
								0,0,0,0,0,0,0,0,0,0,&
								0,0,0,0,0,0,0,0,0,0,&
								0,0,0,0,0,0,0,0,0,0,&
								0,0,0,0,0,0,0,0,0,0)

	head%fldr	=	ishot
	head%ns		=	lt
	head%dt		=	nint(dt*1.0E6)
	head%sx		=	sx

	write(currt, '(I8)')ishot
	open(113, file=trim(adjustl(fn1))//'_ishot'//trim(adjustl(currt))//'.su', access='direct', &
			recl=lbyte*(lt+60), status='replace')
	
	do itr=1, ntrace
			head%gx	=gx(itr)
			write(113, rec=itr)head, (shot_gather(it, itr), it=1, lt)
	enddo
	close(113)

	
	isok=0
    return

    end	subroutine

!***********************************************************************



















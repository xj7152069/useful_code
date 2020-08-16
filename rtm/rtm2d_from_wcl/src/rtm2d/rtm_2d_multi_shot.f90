

!***********************************************************************
	  module	rtm_2d_multi_shot_parameter
	  use	global	
	  implicit none

	  	!==========================================================
	  	!requisite variable for multi shot rtm
	  	!read parameter

		character(len=para_char_flen)	::	fn_cs, fn_vel, fn_grad
		character(len=para_char_flen)  	::  fn_shot_info
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

		integer	::	ns_total
	  	integer	::	nshot
	  	integer	::	ntrace_max
		integer ::	ntrace
	  	integer	::	iflag_shot_info
	  	integer	::	iflag_mig
		integer	::	iflag_mig_iter

		!index_table info
		integer*8 ::	krec_save
	  	integer	::	ishot
		real	::	cx_min
		real	::	cx_max
		real	::	cz_min


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
	  	integer	::	shotfldr

		real	::	vmax
	  	real	::	vmin
	  	real	::	mem_need

	  	!==========================================================
	  	!memory parameter
	  	real,parameter	::	mem_g	=	1024.0*1024.0*1024.0
	  	real,parameter	::	memmax	=	50.0
	  	integer,parameter	::	npad	=	5
	  	integer,parameter	::	order	=	10


	  	!==========================================================
	  	!info 
	  	integer,parameter	::	iflag_debug=1
	  		!=0		no debug  info
	  		!=1		print debug info
	  		!=2		some temp debug info
		character(len=para_char_flen)  ::  currt, currtsnap
		character(len=para_char_flen),parameter	::	fn_debug= './log/rtm2d_debug_'
	  	integer,parameter	::	debug_ishot = 1

	  	integer,parameter   ::  iflag_src = 0


		integer ::	ierr
		integer	::	isok

	  end module	rtm_2d_multi_shot_parameter

!*************************************************************************

	subroutine	rtm_2d_multi_shot_mpi(para_int, para_long, &
					para_float, para_double, para_char, &
					mpi_np, myid)

	use global
	use shotgather_info
	use header_module
	use rtm_2d_multi_shot_parameter
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
	real,allocatable	::	shotobs(:,:)
	real,allocatable	::	wavelet(:)
	real,allocatable	::	image(:,:)
	real,allocatable	::	image_vel(:,:)
	real,allocatable	::	grad(:,:)
	real,allocatable	::	gx(:), gz(:)
	integer,allocatable	::	index_gx(:), index_gz(:)
	real	::	sx, sz
	real	::	imagemin, imagemax

	integer	::	master
	integer	::	i, ii
	integer ::	isend
	integer	::	ix, iz
	integer	::	iix, iiz
	integer	::	ivx, ivz
	integer	::	it, itr, itrace
	integer	::	shotnum
	integer	::	iflag_single_shot
	integer	::	num_grad = 1
	real	::	misfit_loc, misfit_sum


	real	::	time1, time2, time


	!===========================================================================
	nvx = para_int(1)
	nvz = para_int(2)
	master=0

	!===========================================================================
	allocate(vel(nvz, nvx), stat=ierr)
	vel=0.0
	allocate(image_vel(nvz, nvx), stat=ierr)
	image_vel=0.0
	allocate(grad(nvz, nvx), stat=ierr)
	grad=0.0

		!打开速度文件
		fn_vel	=	para_char(1)
		open(21, file=trim(adjustl(fn_vel)), access='direct', status='old', recl=lbyte*nvz*nvx, iostat=ierr)
		if(ierr /= 0)then
			write(*,*)  'Velocity file cannot open right, please check it!!!'
			write(*,*)  'Velocity file ',fn_vel
			stop
		endif
		read(21, rec=1)((vel(ivz, ivx), ivz=1, nvz), ivx=1, nvx)
		close(21)

	vmax=maxval(vel)
	vmin=minval(vel)
	para_float(34)	= vmin
	para_float(35)	= vmax

	!===========================================================================
	CALL rtm_2d_multi_shot_initialize_checkset(para_int, para_long, para_float, &
			para_double, para_char, myid)

    CALL MPI_BARRIER(MPI_COMM_WORLD, IERR)
	!===========================================================================
	!炮检点坐标信息
	 allocate(shotreceiver(ns_total))
     do ishot =1, ns_total
      		allocate(shotreceiver(ishot)%gx(ntrace_max))
      		allocate(shotreceiver(ishot)%gy(ntrace_max))
      		allocate(shotreceiver(ishot)%gz(ntrace_max))
	 enddo

	!===========================================================================
	!从文件中读取需要的炮检点信息 
	if(myid .eq. 0)then
    	call index_shotgathers_and_write_info(fn_cs, fn_shot_info, &
					ns_start, ns_end, ns_interval, ns_total, nshot, &
					shotreceiver, ntrace_max, lt_real, &
					myid, isok)
	endif
    CALL MPI_BARRIER(MPI_COMM_WORLD, IERR)
  	CALL MPI_BCAST(nshot,  1, MPI_integer, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BARRIER(MPI_COMM_WORLD, IERR)

	iflag_shot_info = 1
	if(iflag_shot_info .eq. 1 )then
    	call Read_shotgathers_info(fn_shot_info, &
					ns_start, ns_end, ns_interval, ns_total, nshot, &
					shotreceiver, ntrace_max, lt_real, &
					num_grad, myid, isok)
	endif

	!*======================================================================*
    CALL MPI_BARRIER(MPI_COMM_WORLD, IERR)

	!*=======================================================================*
    if( myid .eq. master ) then
		do i=1, nshot
			CALL MPI_RECV(shotnum, 1, MPI_INTEGER, MPI_ANY_SOURCE,&
					MPI_ANY_TAG, MPI_COMM_WORLD, STATUS, IERR)
			isend=status(mpi_source)
			shotnum=i
			CALL MPI_SEND(shotnum, 1, MPI_INTEGER, ISEND, I,&
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
			ishot    = shotnum
			ntrace   = shotreceiver(ishot)%ntr
			shotfldr = shotreceiver(ishot)%fldr

			!===============================================================
			if(iflag_debug  .eq. 0)then
				print*, 'Run The Current Shot Info : ', ishot, shotfldr, ntrace
			endif

			if(ntrace .lt. 1)then
				write(*,*)'The INput NTRACE IS ERROR'
			endif

			!===============================================================
			!读取当前炮数据，并记录当前炮的炮检点坐标
			allocate(gx(ntrace))
			allocate(gz(ntrace))
			allocate(index_gx(ntrace))
			allocate(index_gz(ntrace))
			gx = 0.0
			gz = 0.0
			index_gx = 0
			index_gz = 0

			cx_min = 999999.9
			cx_max = -999999.9
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

			!===============================================================
			!*-------- input the current shot gather data into the working bufer====*
			nx=(cx_max-cx_min)/dx+1.5

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

			it_step = 1
			ns = lt
			!========================================================
			!
			l1  = 1
			l2  = l1 + nnx*nnz							!vv
			l3  = l2 + lt								!wavelet
			l4  = l3 + ntrace*3							!gx,gy,gz
			l5  = l4 + ntrace*3							!index_gx,index_gy,index_gz
			l6  = l5 + nz*nx							!image
			l7  = l6 + lt*ntrace*3						!shotobs,shotcal,shotres

			l10 = l7
			l11 = l10 + (npad+nnx+npad)*(npad+nnz+npad)				!wavefieldforward1_us1
			l12 = l11 + (npad+nnx+npad)*(npad+nnz+npad)				!wavefieldforward1_us2
			l13 = l12 + (npad+nnx+npad)*(npad+nnz+npad)				!wavefieldforward1_us3
			l14 = l13 + (npad+nnx+npad)*(npad+nnz+npad)				!wavefieldforward1_ur1
			l15 = l14 + (npad+nnx+npad)*(npad+nnz+npad)				!wavefieldforward1_ur2
			l16 = l15 + (npad+nnx+npad)*(npad+nnz+npad)				!wavefieldforward1_ur3

			!pml_for x direction
			l17 = l16 + nnz*(nx_bound_l+nx_bound_r)		!ux_a1 
			l18 = l17 + nnz*(nx_bound_l+nx_bound_r)		!ux_a2 
			l19 = l18 + nnz*(nx_bound_l+nx_bound_r)		!ux_a3 
			l20 = l19 + nnz*(nx_bound_l+nx_bound_r)		!ux_b1 
			l21 = l20 + nnz*(nx_bound_l+nx_bound_r)		!ux_b2 
			l22 = l21 + nnz*(nx_bound_l+nx_bound_r)		!ux_bp1 
			l23 = l22 + nnz*(nx_bound_l+nx_bound_r)		!ux_bp2
			l24 = l23 + nnz*(nx_bound_l+nx_bound_r)		!ux_bp3
			l25 = l24 + nnz*(nx_bound_l+nx_bound_r)		!ux_c1
			l26 = l25 + nnz*(nx_bound_l+nx_bound_r)		!ux_c2 
			l27 = l26 + nnz*(nx_bound_l+nx_bound_r)		!ux_c3 

			!pml_for z direction
			l28 = l27 + (nz_bound_u+nz_bound_d)*nnx		!uz_a1
			l29 = l28 + (nz_bound_u+nz_bound_d)*nnx		!uz_a2
			l30 = l29 + (nz_bound_u+nz_bound_d)*nnx		!uz_a3
			l31 = l30 + (nz_bound_u+nz_bound_d)*nnx		!uz_b1
			l32 = l31 + (nz_bound_u+nz_bound_d)*nnx		!uz_b2
			l33 = l32 + (nz_bound_u+nz_bound_d)*nnx		!uz_bp1
			l34 = l33 + (nz_bound_u+nz_bound_d)*nnx		!uz_bp2
			l35 = l34 + (nz_bound_u+nz_bound_d)*nnx		!uz_bp3
			l36 = l35 + (nz_bound_u+nz_bound_d)*nnx		!uz_c1
			l37 = l36 + (nz_bound_u+nz_bound_d)*nnx		!uz_c2
			l38 = l37 + (nz_bound_u+nz_bound_d)*nnx		!uz_c3


			if(iflag_src .eq. 1)then

				l50 = l49 + nx_with_apert*ns*order				!top
				l51 = l50 + nx_with_apert*ns*order				!bot
				l54 = l53 + nz_with_apert*ns*order				!lef
				l55 = l54 + nz_with_apert*ns*order				!rig

				ltotal=l55
				mem_need=ltotal*4.0/mem_g

				if(mem_need .gt. memmax)then
					print*, 'WARNING!!!'
					print*, 'The needed memory exceed the set threshold :', memmax, ' Gb. '	
					print*, 'Now, adjust the boundary store to reduce the memory...'
					do while  (mem_need .gt. memmax )
						print*, 'The boundary Value ', lt, ns, it_step
						print*, 'Current memory ', mem_need, 'G, Max Memory ', memmax
						ns=ns/2
						it_step=it_step*2
						l50 = l49 + nx_with_apert*ns*order				!top
						l51 = l50 + nx_with_apert*ns*order				!bot
						l54 = l53 + nz_with_apert*ns*order				!lef
						l55 = l54 + nz_with_apert*ns*order				!rig
						ltotal=l55
						mem_need=ltotal*4.0/mem_g

						if(it_step .gt. 5)then
							print*, 'The inverval too large. Please check the parameter.'	
							print*, 'Shut down...'
							stop
						endif
					enddo
				endif

			else if(iflag_src .eq. 0)then

			endif


			!allocate some variables for single shot rtm 
			allocate(image(nz_with_apert, nx_with_apert))
			allocate(shotobs(lt, ntrace))
			allocate(vv(nnz, nnx))
			allocate(wavelet(lt))
			image = 0.0
			shotobs = 0.0
			vv =0.0
			wavelet = 0.0


			!*input the velocity field for the current shot extrapolation-----*
    		call Get_Current_Shot_Velocity_2d(vel, vv, &
					nvz, nvx, nnz, nnx, &
					nz_bound_u, nx_bound_l, &
					nvzz_shift, nvxx_shift, &
					isok)

			!===============================================================
			!读取当前炮数据，并记录当前炮的炮检点坐标
			CALL Input_Shot_Gather_2d(shotobs, shotreceiver, &
						sx, sz, gx, gz, &
						ishot, shotfldr, ntrace, ns_total, &
						lt, lt_real, dt, dt_real, isok)

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
		!	sig_para_int(1)		=	ishot
			sig_para_int(1)		=	shotfldr
			sig_para_int(2)		=	myid
			sig_para_int(3)		=	lt
			sig_para_int(4)		=	ntrace
			sig_para_int(5)		=	ns
			sig_para_int(6)		=	order
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
			sig_para_int(24)	=	it_step
			sig_para_int(25)	=	num_grad
			sig_para_int(26)	=	npad

			sig_para_float(1)	=	dt
			sig_para_float(2)	=	dx
			sig_para_float(3)	=	dz


			if(iflag_debug .eq. 1 .and. ishot .eq. debug_ishot  .and. num_grad .eq. 1)then
				write(currt, '(I8)')nnz
				open(113, file=trim(fn_debug)//'nnz'//trim(adjustl(currt))//'_vv.dat', access='direct', &
						recl=lbyte*nnz*nnx, status='replace')
				write(113, rec=1)((vv(iz, ix), iz=1, nnz), ix=1, nnx)
				close(113)

				write(currt, '(I8)')nvz
				open(113, file=trim(fn_debug)//'nvz'//trim(adjustl(currt))//'_vel.dat', access='direct', &
						recl=lbyte*nvz*nvx, status='replace')
				write(113, rec=1)((vel(iz, ix), iz=1, nvz), ix=1, nvx)
				close(113)
			endif
	

			if(iflag_debug .eq. 1 .and. num_grad .eq. 1)then
				open(111, file=trim(fn_debug)//'lsmig_3d_sig_parameters_output.txt', status='replace')	
				write(111, *) 'The current shot ', shotfldr, '  Parameter is shown as: '
				write(111, *)'ishot', shotfldr
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
				write(111, *)'nx_with_apert = ', nx_with_apert, 'nz_with_apert = ', nz_with_apert
				write(111, *)'nnx = ', nnx, 'nnz = ', nnz
				write(111, *)'nvxx_shift = ', nvxx_shift,   'nvzz_shift = ', nvzz_shift
				write(111, *)'nx_shift = ', nx_shift,  ' nz_shift = ', nz_shift
				write(111, *)'sx, sz', sx ,sz
				write(111, *)'ns_x, ns_z', ns_x, ns_z
				write(111, *)'mem_need = ', mem_need
				write(111, *)'The memory usage : '
				write(111, *)'    tmp       : ', 4.0*(l6-l1)/mem_g
				write(111, *)'    wavefield : ', 4.0*(l11-l10)/mem_g, 4.0*(l16-l10)/mem_g
				write(111, *)'    pml x     : ', 4.0*(l27-l16)/mem_g
				write(111, *)'    pml z     : ', 4.0*(l49-l38)/mem_g
				write(111, *)'    boudary   : ', 4.0*(l55-l49)/mem_g
				close(111)
			endif

			if(iflag_debug  .eq. 0)then
				print*, '    RUN Single Shot RTM: Begin : ', ishot
			!	time1=MPI_WTIME()
			endif
			!*=====================================================================*
			!* single shot reverse time extrapolation and getting image

		!		if(ishot .eq. 16)then
    		call rtm_2d_single_shot(sig_para_int, sig_para_long, &
						sig_para_float, sig_para_double, sig_para_char, &
						wavelet, shotobs, vv, image, &
						sx, sz, gx, gz, &
						ns_x, ns_z, &
						index_gx, index_gz, &
						lt, ntrace, nnx, nnz, nx_with_apert, nz_with_apert, &
						misfit_loc)
		
		!		endif

    		call local_to_all_image_result_2d(image_vel, image, nz_with_apert, nx_with_apert, nvz, nvx, &
						nx_apert_l, nz_apert_u, &
						nvxx_shift, nvzz_shift, &
						isok)


				
			imagemax=maxval(image_vel)
			imagemin=minval(image_vel)
			print*, 'num_grad, ishot, misfit_loc', num_grad, ishot, misfit_loc
			print*, 'num_grad, ishot, max, min', num_grad, ishot, imagemax, imagemin


			if(iflag_debug .eq. 1 .and. ishot .eq. debug_ishot .and. num_grad .eq. 1)then
				write(currt, '(I8)')nz_with_apert
				open(113, file=trim(fn_debug)//'nz'//trim(adjustl(currt))//'_image_v2.dat', access='direct', &
						recl=lbyte*nz_with_apert*nx_with_apert, status='replace')
				write(113, rec=1)((image(iz, ix), iz=1, nz_with_apert), ix=1, nx_with_apert)
				close(113)

			endif

			if(iflag_debug  .eq. 0)then
				print*, '    RUN Single Shot RTM: End : ', ishot
			endif

			deallocate(wavelet)
			deallocate(shotobs)
			deallocate(vv)
			deallocate(image)
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
	CALL MPI_REDUCE(image_vel, grad, nvx*nvz, mpi_real, mpi_sum, master,mpi_comm_world,ierr)

	!*======================================================================*
    CALL MPI_BARRIER(MPI_COMM_WORLD, IERR)
	if(myid .eq. 0)then
		open(21, file=fn_grad, access='direct', status='replace', recl=lbyte*nvz*nvx, iostat=ierr)
		if(ierr /= 0)then
			write(*,*)  'Grad file cannot open right, please check it!!!'
			write(*,*)  'Grad file ',fn_grad
			stop
		endif
		write(21, rec=1)((grad(ivz, ivx), ivz=1, nvz), ivx=1, nvx)
		close(21)
	endif

	deallocate(vel)
	deallocate(image_vel)
	deallocate(grad)
    do ishot =1, ns_total
      	deallocate(shotreceiver(ishot)%gx)
      	deallocate(shotreceiver(ishot)%gy)
      	deallocate(shotreceiver(ishot)%gz)
	enddo
	deallocate(shotreceiver)


	return
	end subroutine

!***********************************************************************
	subroutine	rtm_2d_multi_shot_initialize_checkset(para_int, para_long, para_float, para_double, para_char, myid)
	
	use global
	use rtm_2d_multi_shot_parameter
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
	fn_grad			=	para_char(3)
	fn_shot_info	=	para_char(4)


	nvx				=	para_int(1)
	nvz				=	para_int(2)
	lt				=	para_int(3)
	lt_real			=	para_int(4)
	ns_start		=	para_int(5)
	ns_interval		=	para_int(6)
	ns_end			=	para_int(7)
	iflag_shot_info =	para_int(8)
	ntrace_max		=	para_int(9)



	fmain		=	para_float(1)
	dt			=	para_float(2)
	dt_real		=	para_float(3)
	aperture_x_l=	para_float(4)	
	aperture_x_r=	para_float(5)
	aperture_z_u=	para_float(6)	
	aperture_z_d=	para_float(7)
	boundary_x_l=	para_float(8)
	boundary_x_r=	para_float(9)
	boundary_z_u=	para_float(10)
	boundary_z_d=	para_float(11)
	depth		=	para_float(12)
	dx			=	para_float(13)
	dz			=	para_float(14)
	dvx			=	para_float(15)
	dvz			=	para_float(16)
	cvx_initial =	para_float(17)
	cvz_initial =	para_float(18)


	!*change the dt ms to s 
	dt = dt*0.001
	dt_real = dt_real*0.001

	ns_total =	(ns_end - ns_start)/ns_interval +1
	nz=depth/dz+1

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
	v1=vmax*sqrt(1.0/(dx*dx)  + 1.0/(dz*dz))
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

	open(11, file=fn_cs, access='direct', recl=lbyte*(60+lt_real), status='old', iostat=ierr)
	if(ierr /= 0)then
		write(*,*)	'Shot gather file cannot open right, please check it !!!'
		write(*,*)	'Shutting down the program'
		stop
	endif
	close(11)


	if(iflag_debug .eq. 1)then
		if(myid .eq. 0)then
			open(111, file=trim(fn_debug)//'multi_shot_parameters_output.txt', status='replace')	
			write(111, *)"fn_cs : ", trim(fn_cs)
			write(111, *)'fmain  = ', fmain
			write(111, *)'lt = ', lt
			write(111, *)'lt_real = ', lt_real
			write(111, *)'dt  = ',  dt
			write(111, *)'dt_real = ', dt_real
			write(111, *)'ns_start = ', ns_start
			write(111, *)'ns_interval = ', ns_interval
			write(111, *)'ns_end = ', ns_end
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
			write(111, *)'depth = ', depth
			write(111, *)'dx, dz = ', dx, dz
			write(111, *)'cvx_initial = ', cvx_initial
			write(111, *)'cvz_initial = ', cvz_initial
			write(111, *)'ns_total = ', ns_total
			write(111, *)'nz = ', nz
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






















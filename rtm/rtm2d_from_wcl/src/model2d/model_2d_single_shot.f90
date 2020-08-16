

!***********************************************************************
	  module	model_2d_single_shot_parameter
	  use	global	
	  implicit none

	  	!==========================================================
	  	!==========================================================
	  	!requisite variable for single shot rtm
	  	!requisite variable for rtm
	  	!read parameter
		integer	::	ishot
	  	integer	::	lt
	  	integer	::	ns
	  	integer	::	order
	  	integer	::	ntr
	  	integer	::	nwt
	  	integer	::	nvx, nvz
	  	integer	::	nx, nz
	  	integer	::	nnx, nnz
	  	integer	::	nx_with_apert, nz_with_apert
		integer ::	nx_apert_l, nx_apert_r
		integer ::	nz_apert_u, nz_apert_d
		integer ::	nx_bound_l, nx_bound_r
		integer ::	nz_bound_u, nz_bound_d
		integer ::	nvxx_shift, nvzz_shift
		integer	::	nx_shift, nz_shift
		integer	::	nr_x, nr_z

	  	real	::	dt
	  	real	::	dvx, dvz
		real	::	dx, dz
		real	::	cvx_initial, cvz_initial
		real	::	cx_min, cz_min

	  	integer	::	myid
	  	integer	::	it_step
		integer	::	npad



	  	!==========================================================
	  	!info 
		integer,parameter	::	iflag_src=1		
			!=1     source wavefiled reconstruct
	  		!=0		store the source waveifled 
	  	integer,parameter	::	iflag_debug=1
	  		!=0		no debug  info
	  		!=1		print debug info
	  		!=2		some temp debug info
	  	integer,parameter	::	debug_ishot = 1
	  	character(len=para_char_flen)	::	fn_debug='./log/model2d_debug_'

	  	integer,parameter	::	iflag_snap	=	0
	  		!=0		no snapshot output
	  		!=1		snapshot output
	  	integer,parameter	::	it_snap	=	500
		character(len=para_char_flen)  ::  currt, currtsnap
		character(len=para_char_flen)  ::  snapx_fn='./log/snap/model2d_snap_'

	  	!==========================================================
	  	!extrapolation coefficient
	  	real	::	coe(18)
	  	!pml parameter
	  	real,parameter	::	r=1.0e-3


	  	!==========================================================
		integer	::	ierr
		integer	::	isok

	  end module	model_2d_single_shot_parameter

!***********************************************************************
	  module	model_2d_single_shot_extrapolation

	  implicit none

    	!wavefield 
	  	type wavefield
	  		real,allocatable	::	u(:,:)
	  	end	type

		type(wavefield), target	::	us1, us2, us3

	  	type(wavefield), target	::	ux_a1, ux_a2, ux_a3
	  	type(wavefield), target	::	ux_b1, ux_b2
	  	type(wavefield), target	::	ux_bp1, ux_bp2, ux_bp3
	  	type(wavefield), target	::	ux_c1, ux_c2, ux_c3

	  	type(wavefield), target	::	uz_a1, uz_a2, uz_a3
	  	type(wavefield), target	::	uz_b1, uz_b2
	  	type(wavefield), target	::	uz_bp1, uz_bp2, uz_bp3
	  	type(wavefield), target	::	uz_c1, uz_c2, uz_c3


	  	type(wavefield), pointer	::	pus1, pus2, pus3

	  	type(wavefield), pointer	::	pux_a1, pux_a2, pux_a3
	  	type(wavefield), pointer	::	pux_b1, pux_b2
	  	type(wavefield), pointer	::	pux_bp1, pux_bp2, pux_bp3
	  	type(wavefield), pointer	::	pux_c1, pux_c2, pux_c3

	  	type(wavefield), pointer	::	puz_a1, puz_a2, puz_a3
	  	type(wavefield), pointer	::	puz_b1, puz_b2
	  	type(wavefield), pointer	::	puz_bp1, puz_bp2, puz_bp3
	  	type(wavefield), pointer	::	puz_c1, puz_c2, puz_c3

	  	type(wavefield), pointer	::	pt
	  	
	  	!pml
		real,allocatable	::	coe_1st(:), coe_2nd(:)
		real,allocatable	::	coe_1st_dx(:), coe_1st_dz(:)
		real,allocatable	::	coe_2nd_dx2(:), coe_2nd_dz2(:)
		real,allocatable	::	funa(:,:), dfuna(:,:)
		integer,allocatable	::	line_xl(:,:), line_xr(:,:)
		integer,allocatable	::	line_zu(:,:), line_zd(:,:)

	  	

	  end module model_2d_single_shot_extrapolation
!***********************************************************************

!***********************************************************************
    subroutine model_2D_single_shot(para_int, para_long, &
					para_float, para_double, para_char, &
					wavelet, shot_gather, vv, &
					ns_x, ns_z, &
					index_recgx, index_recgz, &
					ltd, ntrd, nnxd, nnzd, nxd, nzd )

	use	global
	use model_2d_single_shot_parameter
	use model_2d_single_shot_extrapolation
	implicit none
	include 'mpif.h'
	!Data dictionary: declare variable types, definitons, & units
	!Dummy Variables
	integer     ::  para_int(para_int_fnum)
	integer*8   ::  para_long(para_long_fnum)
	real        ::  para_float(para_float_fnum)
	real*8      ::  para_double(para_double_fnum)
	character(len=para_char_flen)	::	para_char(para_char_fnum)

	real	::	wavelet(ltd)
	real	::	shot_gather(ltd, ntrd)
	real	::	vv(nnzd, nnxd)
	integer	::	index_recgx(ntrd), index_recgz(ntrd)

	integer	::	ltd, ntrd
	integer	::	nnxd, nnzd
	integer	::	nxd, nzd
	integer	::	ns_x, ns_z

	!Local Variables
	integer	::	ix, iz
	integer ::	it, iit, iii, itt
	integer ::	iix, inx
	integer ::	iiz, inz
	integer	::	itr, ii
	real	::	time1, time2


	CALL model_2d_single_shot_initialize_checkset(para_int, para_long, para_float, para_double, para_char)

	CALL Form_Coe_extrapolation_2d(coe, dx, dz, dt)

	CALL model_2D_SINGLE_SHOT_ALLOCATE()

	call coefficient_2nd(npad, coe_2nd)
	call coefficient_1st(npad, coe_1st)
	do ii=1, npad
		coe_2nd_dx2(ii)=coe_2nd(ii)/(dx*dx)
		coe_2nd_dz2(ii)=coe_2nd(ii)/(dz*dz)
		coe_1st_dx(ii)=coe_1st(ii)/dx
		coe_1st_dz(ii)=coe_1st(ii)/dz
	enddo

	call  absorbing_function_2d(funa, dfuna, vv, nnz, nnx, dx, dz, &
				line_xl, line_xr, line_zu, line_zd, &
				nx_bound_l, nx_bound_r, &
				nz_bound_u, nz_bound_d, &
				r)

	!*======================================================================= 
	!*--------- extrapolating source wavefield forward time ----------------*      
	!*=======================================================================
	pus1%u(ns_z, ns_x)=wavelet(1)
	do it=2, lt + nwt
	
		itt=it-nwt
		if(iflag_debug .le. 2)then
			if(mod(itt, it_snap).eq. 1) then
				write(*,*)'shot=',ishot, 'myid=',myid,'  step=forward  ','it=',itt
        	endif
		endif

!		goto 1211

		CALL  Extrapolation_2D_Add_PML_One_Step(vv, pus1%u, pus2%u, pus3%u, &
					pux_a1%u, pux_a2%u, pux_a3%u, &
					pux_b1%u, pux_b2%u, pux_bp1%u, pux_bp2%u, pux_bp3%u, &
					pux_c1%u, pux_c2%u, pux_c3%u, &
					puz_a1%u, puz_a2%u, puz_a3%u, &
					puz_b1%u, puz_b2%u, puz_bp1%u, puz_bp2%u, puz_bp3%u, &
					puz_c1%u, puz_c2%u, puz_c3%u, &
					funa, dfuna, &
					line_xl, line_xr, line_zu, line_zd, &
					coe_2nd_dx2, coe_2nd_dz2, &
					coe_1st_dx, coe_1st_dz, &
					npad, nnx, nnz, nx, nz, &
					nx_bound_l, nx_bound_r, &
					nz_bound_u, nz_bound_d, &
					dx, dz, dt, r, itt)

!1211	continue
		goto 1212				
		CALL  Extrapolation_2D_Add_PML_One_Step_old(vv, pus1%u, pus2%u, pus3%u, &
					pux_a1%u, pux_a2%u, pux_a3%u, &
					pux_b1%u, pux_b2%u, pux_bp1%u, pux_bp2%u, pux_bp3%u, &
					pux_c1%u, pux_c2%u, pux_c3%u, &
					puz_a1%u, puz_a2%u, puz_a3%u, &
					puz_b1%u, puz_b2%u, puz_bp1%u, puz_bp2%u, puz_bp3%u, &
					puz_c1%u, puz_c2%u, puz_c3%u, &
					coe, nnx, nnz, nx, nz, &
					nx_bound_l, nx_bound_r, &
					nz_bound_u, nz_bound_d, &
					dx, dz, dt, r, itt)
1212	continue
		!===========================================================
		!add source function
		if(it .le. lt)then
			pus2%u(ns_z, ns_x)=pus2%u(ns_z, ns_x)+ wavelet(it-1)
		endif


		if(iflag_snap  .eq. 1 .and. ishot .eq. debug_ishot)then
			if(mod(itt, it_snap) .eq. 1)then
				write(currt, '(I8)')itt
				write(currtsnap,*)trim(adjustl(snapx_fn)),trim(adjustl(currt))
				write(currt, '(I8)')nnz
				open(113, file=trim(currtsnap)//'t_nnz'//trim(adjustl(currt))//'_us3.dat', access='direct', &
					recl=lbyte*nnx*nnz, status='replace')
				write(113, rec=1)((pus3%u(iiz, iix), iiz=1, nnz), iix=1, nnx)
				close(113)
			endif
		endif


		!===========================================================
		!store the shotcal 
		if(itt .gt. 0)then
			do itr=1, ntr
				nr_x=index_recgx(itr)
				nr_z=index_recgz(itr)
				if(nr_x .ge. 1 .and. nr_z .ge. 1 .and. nr_x .le. nnx .and. nr_z .le. nnz)then
					shot_gather(itt, itr)=pus3%u(nr_z, nr_x)
				endif
			enddo
		endif

		!===========================================================
		!update wavefield 
		pt   => pus1
		pus1 => pus2
		pus2 => pus3
		pus3 => pt

		!x direction
		pt     => pux_a1
		pux_a1 => pux_a2
		pux_a2 => pux_a3
		pux_a3 => pt

		pt		=> pux_bp1
		pux_bp1 => pux_bp2
		pux_bp2 => pux_bp3
		pux_bp3 => pt

		pt	   => pux_b1
		pux_b1 => pux_b2
		pux_b2 => pt

		pt	   => pux_c1
		pux_c1 => pux_c2
		pux_c2 => pux_c3
		pux_c3 => pt

		!z direction
		pt     => puz_a1
		puz_a1 => puz_a2
		puz_a2 => puz_a3
		puz_a3 => pt

		pt		=> puz_bp1
		puz_bp1 => puz_bp2
		puz_bp2 => puz_bp3
		puz_bp3 => pt

		pt	   => puz_b1
		puz_b1 => puz_b2
		puz_b2 => pt

		pt	   => puz_c1
		puz_c1 => puz_c2
		puz_c2 => puz_c3
		puz_c3 => pt

	enddo


	!***************************************************************************!

	CALL model_2D_SINGLE_SHOT_DEALLOCATE() 


	!*======================================================================*
	!*     The current shot gather reverse time migration finished          *
	!*======================================================================*
	isok=0

    return
    end	subroutine

!***********************************************************************
	subroutine	model_2d_single_shot_initialize_checkset(para_int, para_long, para_float, para_double, para_char)
	
	use global
	use model_2d_single_shot_parameter 
	implicit none

	!Data dictionary: declare variable types, definitons, & units
	!Dummy Variables
	integer     ::  para_int(para_int_fnum)
	integer*8   ::  para_long(para_long_fnum)
	real        ::  para_float(para_float_fnum)
	real*8      ::  para_double(para_double_fnum)
	character(len=para_char_flen)	::	para_char(para_char_fnum)

	ishot		=	para_int(1)
	myid		=	para_int(2)
	lt			=	para_int(3)
	ntr			=	para_int(4)
!	ns			=	para_int(5)
!	order		=	para_int(6)
	nwt			=	para_int(7)
	nx			=	para_int(8)
	nz			=	para_int(9)
	nx_apert_l	=	para_int(10)
	nx_apert_r	=	para_int(11)
	nz_apert_u	=	para_int(12)
	nz_apert_d	=	para_int(13)
	nx_bound_l	=	para_int(14)
	nx_bound_r	=	para_int(15)
	nz_bound_u	=	para_int(16)
	nz_bound_d	=	para_int(17)
	nx_with_apert	=	para_int(18)	
	nz_with_apert	=	para_int(19)
	nnx			=	para_int(20)
	nnz			=	para_int(21)
	nx_shift	=	para_int(22)
	nz_shift	=	para_int(23)
!	it_step		=	para_int(24)
	npad		=	para_int(25)


	dt			=	para_float(1)
	dx			=	para_float(2)
	dz			=	para_float(3)


	if(iflag_debug .eq. 1 .and. ishot .eq. debug_ishot)then
		open(111, file=trim(fn_debug)//'sig_parameters_checkset.txt', status='replace')	
		write(111, *)'ishot ', ishot
		write(111, *)'myid ', myid
		write(111, *)'lt ', lt
		write(111, *)'ntr ', ntr
		write(111, *)'nwt ', nwt
		write(111, *)'nx_apert_l = ', nx_apert_l, ' nx_apert_r = ', nx_apert_r
		write(111, *)'nz_apert_u = ', nz_apert_u, ' nz_apert_d = ', nz_apert_d
		write(111, *)'nx_bound_l = ', nx_bound_l, ' nx_bound_r = ', nx_bound_r
		write(111, *)'nz_bound_u = ', nz_bound_u, ' nz_bound_d = ', nz_bound_d
		write(111, *)'nx = ', nx,  ' nz = ', nz
		write(111, *)'nx_with_apert = ', nx_with_apert, 'nz_with_apert = ', nz_with_apert
		write(111, *)'nnx = ', nnx, 'nnz = ', nnz
		write(111, *)'nvxx_shift = ', nvxx_shift,  'nvzz_shift = ', nvzz_shift
		write(111, *)'nx_shift = ', nx_shift, ' nz_shift = ', nz_shift
		write(111, *)'dx = ', dx, ' dz = ', dz
		write(111, *)'dt = ', dt
		
		close(111)
	endif

	isok=0

	end subroutine

!***********************************************************************
	subroutine model_2D_SINGLE_SHOT_ALLOCATE()

	use model_2d_single_shot_parameter
	use model_2d_single_shot_extrapolation
	implicit none
	

	!allocate wavefiled
	allocate(us1%u(-4:nnz+5, -4:nnx+5), stat=ierr)
	if(ierr.ne.0)then
		write(*,*)"Can not allocate working memory about us1, stop!!!"
		stop
	endif

	allocate(us2%u(-4:nnz+5, -4:nnx+5), stat=ierr)
	if(ierr.ne.0)then
		write(*,*)"Can not allocate working memory about us2, stop!!!"
		stop
	endif

	allocate(us3%u(-4:nnz+5, -4:nnx+5), stat=ierr)
	if(ierr.ne.0)then
		write(*,*)"Can not allocate working memory about us3, stop!!!"
		stop
	endif


	!allocate pml
    allocate(ux_a1%u(nnz,nx_bound_l+nx_bound_r), stat=ierr)
	if(ierr.ne.0)then
		write(*,*)"Can not allocate working memory about ux_a1, stop!!!"
		stop
	endif

    allocate(ux_a2%u(nnz,nx_bound_l+nx_bound_r), stat=ierr)
	if(ierr.ne.0)then
		write(*,*)"Can not allocate working memory about ux_a2, stop!!!"
		stop
	endif

    allocate(ux_a3%u(nnz,nx_bound_l+nx_bound_r), stat=ierr)
	if(ierr.ne.0)then
		write(*,*)"Can not allocate working memory about ux_a3, stop!!!"
		stop
	endif

    allocate(ux_b1%u(nnz,nx_bound_l+nx_bound_r), stat=ierr)
	if(ierr.ne.0)then
		write(*,*)"Can not allocate working memory about ux_b1, stop!!!"
		stop
	endif

    allocate(ux_b2%u(nnz,nx_bound_l+nx_bound_r), stat=ierr)
	if(ierr.ne.0)then
		write(*,*)"Can not allocate working memory about ux_b2, stop!!!"
		stop
	endif

    allocate(ux_bp1%u(nnz,nx_bound_l+nx_bound_r), stat=ierr)
	if(ierr.ne.0)then
		write(*,*)"Can not allocate working memory about ux_bp1, stop!!!"
		stop
	endif

    allocate(ux_bp2%u(nnz,nx_bound_l+nx_bound_r), stat=ierr)
	if(ierr.ne.0)then
		write(*,*)"Can not allocate working memory about ux_bp2, stop!!!"
		stop
	endif

    allocate(ux_bp3%u(nnz,nx_bound_l+nx_bound_r), stat=ierr)
	if(ierr.ne.0)then
		write(*,*)"Can not allocate working memory about ux_bp3, stop!!!"
		stop
	endif

    allocate(ux_c1%u(nnz,nx_bound_l+nx_bound_r), stat=ierr)
	if(ierr.ne.0)then
		write(*,*)"Can not allocate working memory about ux_c1, stop!!!"
		stop
	endif

    allocate(ux_c2%u(nnz,nx_bound_l+nx_bound_r), stat=ierr)
	if(ierr.ne.0)then
		write(*,*)"Can not allocate working memory about ux_c2, stop!!!"
		stop
	endif

    allocate(ux_c3%u(nnz,nx_bound_l+nx_bound_r), stat=ierr)
	if(ierr.ne.0)then
		write(*,*)"Can not allocate working memory about ux_c3, stop!!!"
		stop
	endif


    !z direction
    allocate(uz_a1%u(nz_bound_u+nz_bound_d,nnx), stat=ierr)
	if(ierr.ne.0)then
		write(*,*)"Can not allocate working memory about uz_a1, stop!!!"
		stop
	endif

    allocate(uz_a2%u(nz_bound_u+nz_bound_d,nnx), stat=ierr)
	if(ierr.ne.0)then
		write(*,*)"Can not allocate working memory about uz_a2, stop!!!"
		stop
	endif

    allocate(uz_a3%u(nz_bound_u+nz_bound_d,nnx), stat=ierr)
	if(ierr.ne.0)then
		write(*,*)"Can not allocate working memory about uz_a3, stop!!!"
		stop
	endif

    allocate(uz_b1%u(nz_bound_u+nz_bound_d,nnx), stat=ierr)
	if(ierr.ne.0)then
		write(*,*)"Can not allocate working memory about uz_b1, stop!!!"
		stop
	endif

    allocate(uz_b2%u(nz_bound_u+nz_bound_d,nnx), stat=ierr)
	if(ierr.ne.0)then
		write(*,*)"Can not allocate working memory about uz_b2, stop!!!"
		stop
	endif

    allocate(uz_bp1%u(nz_bound_u+nz_bound_d,nnx), stat=ierr)
	if(ierr.ne.0)then
		write(*,*)"Can not allocate working memory about uz_bp1, stop!!!"
		stop
	endif

    allocate(uz_bp2%u(nz_bound_u+nz_bound_d,nnx), stat=ierr)
	if(ierr.ne.0)then
		write(*,*)"Can not allocate working memory about uz_bp2, stop!!!"
		stop
	endif

    allocate(uz_bp3%u(nz_bound_u+nz_bound_d,nnx), stat=ierr)
	if(ierr.ne.0)then
		write(*,*)"Can not allocate working memory about uz_bp3, stop!!!"
		stop
	endif

    allocate(uz_c1%u(nz_bound_u+nz_bound_d,nnx), stat=ierr)
	if(ierr.ne.0)then
		write(*,*)"Can not allocate working memory about uz_c1, stop!!!"
		stop
	endif

    allocate(uz_c2%u(nz_bound_u+nz_bound_d,nnx), stat=ierr)
	if(ierr.ne.0)then
		write(*,*)"Can not allocate working memory about uz_c2, stop!!!"
		stop
	endif

    allocate(uz_c3%u(nz_bound_u+nz_bound_d,nnx), stat=ierr)
	if(ierr.ne.0)then
		write(*,*)"Can not allocate working memory about uz_c3, stop!!!"
		stop
	endif

	allocate(coe_1st(npad), stat=ierr)
	if(ierr.ne.0)then
		write(*,*)"Can not allocate working memory about coe_1st, stop!!!"
		stop
	endif
	allocate(coe_2nd(npad), stat=ierr)
	if(ierr.ne.0)then
		write(*,*)"Can not allocate working memory about coe_2nd, stop!!!"
		stop
	endif
	allocate(coe_1st_dx(npad), stat=ierr)
	if(ierr.ne.0)then
		write(*,*)"Can not allocate working memory about coe_1st_dx, stop!!!"
		stop
	endif
	allocate(coe_1st_dz(npad), stat=ierr)
	if(ierr.ne.0)then
		write(*,*)"Can not allocate working memory about coe_1st_dz, stop!!!"
		stop
	endif
	allocate(coe_2nd_dx2(npad), stat=ierr)
	if(ierr.ne.0)then
		write(*,*)"Can not allocate working memory about coe_2nd_dx2, stop!!!"
		stop
	endif
	allocate(coe_2nd_dz2(npad), stat=ierr)
	if(ierr.ne.0)then
		write(*,*)"Can not allocate working memory about coe_2nd_dz2, stop!!!"
		stop
	endif
	allocate(funa(nnz, nnx), stat=ierr)
	if(ierr.ne.0)then
		write(*,*)"Can not allocate working memory about funa, stop!!!"
		stop
	endif
	allocate(dfuna(nnz, nnx), stat=ierr)
	if(ierr.ne.0)then
		write(*,*)"Can not allocate working memory about dfuna, stop!!!"
		stop
	endif
	allocate(line_xl(nx_bound_l, 2), stat=ierr)
	if(ierr.ne.0)then
		write(*,*)"Can not allocate working memory about line_xl, stop!!!"
		stop
	endif
	allocate(line_xr(nx_bound_r, 2), stat=ierr)
	if(ierr.ne.0)then
		write(*,*)"Can not allocate working memory about line_xr, stop!!!"
		stop
	endif
	allocate(line_zu(nz_bound_u, 2), stat=ierr)
	if(ierr.ne.0)then
		write(*,*)"Can not allocate working memory about line_zu, stop!!!"
		stop
	endif
	allocate(line_zd(nz_bound_d, 2), stat=ierr)
	if(ierr.ne.0)then
		write(*,*)"Can not allocate working memory about line_zd, stop!!!"
		stop
	endif



	!wavefields
	us1%u=0.0
	us2%u=0.0
	us3%u=0.0

	!x direction
    ux_a1%u=0.0
    ux_a2%u=0.0
    ux_a3%u=0.0
    ux_b1%u=0.0
    ux_b2%u=0.0
    ux_bp1%u=0.0
    ux_bp2%u=0.0
    ux_bp3%u=0.0
    ux_c1%u=0.0
    ux_c2%u=0.0
    ux_c3%u=0.0

    !z direction
    uz_a1%u=0.0
    uz_a2%u=0.0
    uz_a3%u=0.0
    uz_b1%u=0.0
    uz_b2%u=0.0
    uz_bp1%u=0.0
    uz_bp2%u=0.0
    uz_bp3%u=0.0
    uz_c1%u=0.0
    uz_c2%u=0.0
    uz_c3%u=0.0

	coe_1st=0.0
	coe_2nd=0.0
	coe_1st_dx=0.0
	coe_1st_dz=0.0
	coe_2nd_dx2=0.0
	coe_2nd_dz2=0.0
	funa=0.0
	dfuna=0.0
	line_xl=0
	line_xr=0
	line_zu=0
	line_zd=0



	!wavefields
	pus1 =>	us1
	pus2 => us2
	pus3 => us3

	!x direction
    pux_a1 => ux_a1
    pux_a2 => ux_a2
    pux_a3 => ux_a3
    pux_b1 => ux_b1
    pux_b2 => ux_b2
    pux_bp1 => ux_bp1
    pux_bp2 => ux_bp2
    pux_bp3 => ux_bp3
    pux_c1 => ux_c1
    pux_c2 => ux_c2
    pux_c3 => ux_c3

    !z direction
    puz_a1 => uz_a1
    puz_a2 => uz_a2
    puz_a3 => uz_a3
    puz_b1 => uz_b1
    puz_b2 => uz_b2
    puz_bp1 => uz_bp1
    puz_bp2 => uz_bp2
    puz_bp3 => uz_bp3
    puz_c1 => uz_c1
    puz_c2 => uz_c2
    puz_c3 => uz_c3


	isok=0

	end subroutine

!***********************************************************************

	subroutine model_2D_SINGLE_SHOT_DEALLOCATE()

	use model_2d_single_shot_parameter
	use model_2d_single_shot_extrapolation
	implicit none


	!deallocate wavefiled
	deallocate(us1%u)
	deallocate(us2%u)
	deallocate(us3%u)

	!deallocate pml
    !x direction
    deallocate(ux_a1%u)
    deallocate(ux_a2%u)
    deallocate(ux_a3%u)
    deallocate(ux_b1%u)
    deallocate(ux_b2%u)
    deallocate(ux_bp1%u)
    deallocate(ux_bp2%u)
    deallocate(ux_bp3%u)
    deallocate(ux_c1%u)
    deallocate(ux_c2%u)
    deallocate(ux_c3%u)

    !z direction
    deallocate(uz_a1%u)
    deallocate(uz_a2%u)
    deallocate(uz_a3%u)
    deallocate(uz_b1%u)
    deallocate(uz_b2%u)
    deallocate(uz_bp1%u)
    deallocate(uz_bp2%u)
    deallocate(uz_bp3%u)
    deallocate(uz_c1%u)
    deallocate(uz_c2%u)
    deallocate(uz_c3%u)

	deallocate(coe_1st)
	deallocate(coe_2nd)
	deallocate(coe_1st_dx)
	deallocate(coe_1st_dz)
	deallocate(coe_2nd_dx2)
	deallocate(coe_2nd_dz2)
	deallocate(funa)
	deallocate(dfuna)
	deallocate(line_xl)
	deallocate(line_xr)
	deallocate(line_zu)
	deallocate(line_zd)



	isok=0

	end subroutine

!***********************************************************************










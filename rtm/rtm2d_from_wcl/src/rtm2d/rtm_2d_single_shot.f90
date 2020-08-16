


!*************************************************************************
	  module	rtm_2d_single_shot_parameter
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
	  	integer	::	num_grad
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
	  	integer,parameter	::	debug_ishot = 16
	  	character(len=para_char_flen)	::	fn_debug='./log/rtm_sig_debug_'

	  	integer,parameter	::	iflag_snap	=	0
	  		!=0		no snapshot output
	  		!=1		snapshot output
	  	integer,parameter	::	it_snap	=	500
		character(len=para_char_flen)  ::  currt, currtsnap
		character(len=para_char_flen)  ::  snapx_fn=  &
			'/ssd/ssd_3.2T/wcl/result/rtm2d/snap/v1_rtm_snap_'

	  	!==========================================================
	  	!extrapolation coefficient
	  	real	::	coe(18)
	  	!pml parameter
	  	real,parameter	::	r=1.0e-3

	  	!==========================================================
	  	!angle gather
	 	real	::	poy_coe(5)
		integer,parameter	::	niter	=	10
		real,parameter		::	alpha	=	1.0

		integer	::	flag_sr
	  	integer	::	wflag

	  	!==========================================================
		integer	::	ierr
		integer	::	isok

	  end module	rtm_2d_single_shot_parameter

!***********************************************************************

    subroutine rtm_2d_single_shot(para_int, para_long, &
					para_float, para_double, para_char, &
					wavelet, shotobs, vv, image, &
					sx, sz, gx, gz, &
					ns_x, ns_z, &
					index_recgx, index_recgz, &
					ltd, ntrd, nnxd, nnzd, nx_with_apertd, nz_with_apertd, &
					misfit_loc)

	use	global
	use rtm_2d_single_shot_parameter 
	use module_rtm_extrapolation_2d
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
	real	::	shotobs(ltd, ntrd)
	real	::	vv(nnzd, nnxd)
	real	::	image(nz_with_apertd, nx_with_apertd)	
	real	::	sx, sz
	real	::	gx(ntrd), gz(ntrd)
	integer	::	index_recgx(ntrd), index_recgz(ntrd)

	integer	::	ltd, ntrd
	integer	::	nnxd, nnzd
	integer	::	nx_with_apertd, nz_with_apertd
	integer	::	ns_x, ns_z
	real	::	misfit_loc

	!Local Variables
	integer	::	ix, iz
	integer ::	it, iit, iii, itt
	integer ::	iix, inx
	integer ::	iiz, inz
	integer	::	itr
	real	::	vsurf
	real	::	v2
	real	::	misfit_sum
	real	::	time1, time2

	integer	::	ii



	call rtm_2d_single_shot_initialize_checkset(para_int, para_long, para_float, para_double, para_char)

	call Form_Coe_extrapolation_2d(coe, dx, dz, dt)

	call rtm_2d_single_shot_extrapolation_allocate(nnx, nnz, &
					nx_with_apert, nz_with_apert, &
					nx_bound_l, nx_bound_r, &
					nz_bound_u, nz_bound_d, &
					npad, order, ns, iflag_src)


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

	if(iflag_debug .eq. 1 .and. ishot .eq. debug_ishot)then
		open(122, file=trim(fn_debug)//'_funa.dat', &
				access='direct', recl=lbyte*nnz*nnx, status='replace')
		write(122, rec=1)((funa(inz, inx), inz=1, nnz), inx=1, nnx)
		close(122)
		open(122, file=trim(fn_debug)//'_dfuna.dat', &
				access='direct', recl=lbyte*nnz*nnx, status='replace')
		write(122, rec=1)((dfuna(inz, inx), inz=1, nnz), inx=1, nnx)
		close(122)
	endif


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

		!===========================================================
		!add source function
		if(it .le. lt)then
			pus2%u(ns_z, ns_x)=pus2%u(ns_z, ns_x)+ wavelet(it-1)
		endif


		if(iflag_snap  .eq. 1 .and. ishot .eq. debug_ishot)then
			if(mod(itt+1, it_snap) .eq. 1)then
				write(currt, '(I8)')itt
				write(currtsnap,*)trim(adjustl(snapx_fn)),trim(adjustl(currt))
				write(currt, '(I8)')nnz
				open(113, file=trim(currtsnap)//'t_nnz'//trim(adjustl(currt))//'_us3.dat', access='direct', &
					recl=lbyte*nnx*nnz, status='replace')
				write(113, rec=1)((pus3%u(iiz, iix), iiz=1, nnz), iix=1, nnx)
				close(113)

				write(currt, '(I8)')itt
				write(currtsnap,*)trim(adjustl(snapx_fn)),trim(adjustl(currt))
				write(currt, '(I8)')nz_with_apert
				open(113, file=trim(currtsnap)//'t_nz_with_apert'//trim(adjustl(currt))//'us3.dat', access='direct', &
					recl=lbyte*nz_with_apert*nx_with_apert, status='replace')
				do ix=1, nx_with_apert
					iix = ix +nx_bound_l
					do iz=1, nz_with_apert
						iiz = iz + nz_bound_u
						image_tmp(iz, ix) = pus3%u(iiz, iix)
					enddo
				enddo
				write(113, rec=1)((image_tmp(iz, ix), iz=1, nz_with_apert), ix=1, nx_with_apert)
				close(113)
			endif
		endif

		if(iflag_src  .eq. 1)then
			!*存储波场边界用于波场重构
			if(mod(itt, it_step) .eq. 0 .and. itt .gt. 0) then
				iit=itt/it_step
			!	!$omp parallel private(inx, inz) 
			!	!$omp do
					do ix=1, nx_with_apert
						inx=ix+nx_bound_l
						do iii=1, order
							top%u(iii, ix, iit)=pus3%u(nz_bound_u+iii, inx)
							bot%u(iii, ix, iit)=pus3%u(nz_bound_u+nz_with_apert-iii+1, inx)
						enddo
					enddo
			!	!$omp enddo

			!	!$omp do
					do iz=1, nz_with_apert
						inz=iz+nz_bound_u
						do iii=1, order
							lef%u(iii, iz, iit)=pus3%u(inz, nx_bound_l+iii)
							rig%u(iii, iz, iit)=pus3%u(inz, nx_bound_l+nx_with_apert-iii+1)
						enddo
					enddo
			!	!$omp enddo
			!	!$omp end parallel
			endif
		else if( iflag_src .eq. 0)then
			if(mod(itt, it_step) .eq. 0 .and. itt .gt. 0) then
				iit=itt/it_step
				ust%u(:, :, iit)=pus3%u
			endif
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


	!***********************************************************************!
	!*	Extrapolating the extrapolated source wavfield and the receiver 
	!*wavefield along the reverse time and extracting the imaging value 
	!*with the correlation of the two wavefield at each time step
	!***********************************************************************!
	!store the last two wavefields
	pt   => pus2
	pus2 => pus1
	pus1 => pt

	!============================================================================
	!*计算检波器附近的速度值
	vsurf=0.0
	do ix=1, nnx
		vsurf = vsurf+vv(ns_z, ix)/1000.0
	enddo
	vsurf=vsurf/(nnx)*1000.0

!	!debug wcl
!	vsurf=2000.0
!
	print*, 'num_grad, ishot, vsurf', num_grad, ishot, vsurf
	if(iflag_debug .eq. 1 .and. ishot .eq. debug_ishot)then
		write(currt, '(I8)')num_grad
		open(122, file=trim(fn_debug)//'num_grad'//trim(adjustl(currt))//'_shotobs_be.dat', access='direct', recl=lbyte*lt*ntr, status='replace')
		write(122, rec=1)((shotobs(it, itr), it=1, lt), itr=1, ntr)
		close(122)

	endif
	call process_directwave_2d(shotobs, misfit_sum, &
					sx, sz, gx, gz, index_recgx, index_recgz, &
					vsurf, lt, ntr, dt, nwt)

	if(iflag_debug .eq. 1 .and. ishot .eq. debug_ishot)then
		write(currt, '(I8)')num_grad
		open(122, file=trim(fn_debug)//'num_grad'//trim(adjustl(currt))//'_shotobs.dat', access='direct', recl=lbyte*lt*ntr, status='replace')
		write(122, rec=1)((shotobs(it, itr), it=1, lt), itr=1, ntr)
		close(122)
	endif


    !Reinitiallization for back propagation
	!x direction
    pux_a1%u=0.0
    pux_a2%u=0.0
    pux_a3%u=0.0
    pux_b1%u=0.0
    pux_b2%u=0.0
    pux_bp1%u=0.0
    pux_bp2%u=0.0
    pux_bp3%u=0.0
    pux_c1%u=0.0
    pux_c2%u=0.0
    pux_c3%u=0.0

    !z direction
    puz_a1%u=0.0
    puz_a2%u=0.0
    puz_a3%u=0.0
    puz_b1%u=0.0
    puz_b2%u=0.0
    puz_bp1%u=0.0
    puz_bp2%u=0.0
    puz_bp3%u=0.0
    puz_c1%u=0.0
    puz_c2%u=0.0
    puz_c3%u=0.0


	!***************************************************************************!
    do it=lt-2, 1, -1			!!!! loop for backrward source and receiver wavefield

		if(iflag_debug .le. 2)then
			if(mod(it, it_snap).eq. 1) then
				write(*,*)'shot=',ishot,'myid=',myid,'  step=reverse  ','it=',it
        	endif
		endif

		!============================================================================
		!震源端波场处理
		if(iflag_src  .eq. 1)then
			!波场重构，重新计算震源波场
			!source wavefiled extrapolation with reverse time
			if(mod(it, it_step) .eq. 0) then
				iit=it/it_step
			!	!$omp parallel private(inx, inz) 
			!	!$omp do
					do ix=1, nx_with_apert
						inx=ix+nx_bound_l
						do iii=1, order
							pus2%u(nz_bound_u+iii, inx) = top%u(iii, ix, iit+1)
							pus2%u(nz_bound_u+nz_with_apert-iii+1, inx) = bot%u(iii, ix, iit+1)
						enddo
					enddo
			!	!$omp enddo
			!	!$omp do
					do iz=1, nz_with_apert
						inz=iz+nz_bound_u
						do iii=1, order
							pus2%u(inz, nx_bound_l+iii) = lef%u(iii, iz, iit+1)
							pus2%u(inz, nx_bound_l+nx_with_apert-iii+1) = rig%u(iii, iz, iit+1)
						enddo
					enddo
			!	!$omp enddo
			!	!$omp end parallel
			endif

			!===========================================================
			!时逆过程，处理反传到二次源（震源），继续传播问题
			if(it + nwt +1 .le. lt)then
				pus2%u(ns_z, ns_x)=pus2%u(ns_z, ns_x) - wavelet(it+nwt+1)
			endif

			CALL Extrapolation_2D_Without_PML_One_Step(vv, pus1%u, pus2%u, pus3%u, &
					coe, nnx, nnz, &
					nx_bound_l, nx_bound_r, &
					nz_bound_u, nz_bound_d, &
					dx, dz, dt)

		else if( iflag_src .eq. 0)then
			!读取存储的波场
			if(mod(it, it_step) .eq. 0 ) then
				iit=it/it_step
				pus1%u = ust%u(:, :, iit+2)
				pus2%u = ust%u(:, :, iit+1)
				pus3%u = ust%u(:, :, iit)
			endif
		endif

		!============================================================================
		!reciever wavefiled extrapolation with reverse time
		do itr=1, ntr
			nr_x=index_recgx(itr)
			nr_z=index_recgz(itr)
			if( nr_x .ge. 1   .and. nr_z .ge. 1   .and. &
				nr_x .le. nnx .and. nr_z .le. nnz )then
				v2=vv(nr_z, nr_x)*vv(nr_z, nr_x)
				if(it .eq. lt-2)then
					pur1%u(nr_z, nr_x) = pur1%u(nr_z, nr_x) + shotobs(it+2, itr)*v2*dt*dt
				endif
					pur2%u(nr_z, nr_x) = pur2%u(nr_z, nr_x) + shotobs(it+1, itr)*v2*dt*dt
			endif
		enddo

		CALL  Extrapolation_2D_Add_PML_One_Step(vv, pur1%u, pur2%u, pur3%u, &
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

		if(iflag_snap  .eq. 1  .and. ishot  .eq. debug_ishot)then
			if(mod(it+1, it_snap) .eq. 1)then
				write(currt, '(I8)')it
				write(currtsnap,*)trim(adjustl(snapx_fn)),trim(adjustl(currt))
				write(currt, '(I8)')nnz
				open(113, file=trim(currtsnap)//'t_nnz'//trim(adjustl(currt))//'_ur3.dat', access='direct', &
					recl=lbyte*nnx*nnz, status='replace')
				write(113, rec=1)((pur3%u(iiz, iix), iiz=1, nnz), iix=1, nnx)
				close(113)


					write(currt, '(I8)')it
					write(currtsnap,*)trim(adjustl(snapx_fn)),trim(adjustl(currt))
					write(currt, '(I8)')nz_with_apert
					open(113, file=trim(currtsnap)//'t_nz_with_apert'//trim(adjustl(currt))//'ur3.dat', access='direct', &
						recl=lbyte*nz_with_apert*nx_with_apert, status='replace')
					do ix=1, nx_with_apert
						iix = ix +nx_bound_l
						do iz=1, nz_with_apert
							iiz = iz + nz_bound_u
							image_tmp(iz, ix) = pur3%u(iiz, iix)
						enddo
					enddo
					write(113, rec=1)((image_tmp(iz, ix), iz=1, nz_with_apert), ix=1, nx_with_apert)
					close(113)

			endif
		endif

		!============================================================================
		!*cross correlation image condition at every time slice
		if(it .gt. 2*nwt+1 .and. it .lt. lt)then
		CALL Cor_Imaging_Condition_2d(image, vv, pus1%u, pus2%u, pus3%u, pur1%u, pur2%u, pur3%u, &
							nnz, nnx, nz_with_apert, nx_with_apert, &
							nz_bound_u, nx_bound_l ,dt)
		endif


			if(iflag_snap  .eq. 1 .and. ishot .eq. debug_ishot)then
				if(mod(it, it_snap) .eq. 1)then
					write(currt, '(I8)')it
					write(currtsnap,*)trim(adjustl(snapx_fn)),trim(adjustl(currt))
					write(currt, '(I8)')nz_with_apert
					open(113, file=trim(currtsnap)//'t_nz_with_apert'//trim(adjustl(currt))//'_image.dat', access='direct', &
						recl=lbyte*nz_with_apert*nx_with_apert, status='replace')

					do ix=1, nx_with_apert
						iix = ix +nx_bound_l
						do iz=1, nz_with_apert
							iiz = iz + nz_bound_u
							image_tmp(iz, ix) = pus3%u(iiz, iix)*pur3%u(iiz, iix)
						enddo
					enddo
					write(113, rec=1)((image_tmp(iz, ix), iz=1, nz_with_apert), ix=1, nx_with_apert)
					close(113)
				endif
			endif


		!============================================================================
		!update wavefield 
		pt   => pur1
		pur1 => pur2
		pur2 => pur3
		pur3 => pt

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


		if(iflag_src .eq. 1)then
			!update source wavefield
			pt   => pus1
			pus1 => pus2
			pus2 => pus3
			pus3 => pt
		endif


	enddo				!!! loop end for backward source and receiver wavefield
	!*===================================================================


	call rtm_2d_single_shot_extrapolation_deallocate(iflag_src)

	!*======================================================================*
	!*     The current shot gather reverse time migration finished          *
	!*======================================================================*
	isok=0

    return
    end	subroutine

!***********************************************************************
	subroutine	rtm_2d_single_shot_initialize_checkset(para_int, para_long, para_float, para_double, para_char)
	
	use global
	use rtm_2d_single_shot_parameter 
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
	ns			=	para_int(5)
	order		=	para_int(6)
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
	it_step		=	para_int(24)
	num_grad		=	para_int(25)
	npad		=	para_int(26)


	dt			=	para_float(1)
	dx			=	para_float(2)
	dz			=	para_float(3)


	if(iflag_debug .eq. 1 .and. ishot .eq. debug_ishot)then
		open(111, file=trim(fn_debug)//'sig_parameters_checkset.txt', status='replace')	
		write(111, *)'ishot ', ishot
		write(111, *)'myid ', myid
		write(111, *)'lt ', lt
		write(111, *)'ntr ', ntr
		write(111, *)'ns ', ns
		write(111, *)'order ', order
		write(111, *)'npad ', npad
		write(111, *)'nwt ', nwt
		write(111, *)'it_step', it_step
		write(111, *)'nx_apert_l = ', nx_apert_l, ' nx_apert_r = ', nx_apert_r
		write(111, *)'nz_apert_u = ', nz_apert_u, ' nz_apert_d = ', nz_apert_d
		write(111, *)'nx_bound_l = ', nx_bound_l, ' nx_bound_r = ', nx_bound_r
		write(111, *)'nz_bound_u = ', nz_bound_u, ' nz_bound_d = ', nz_bound_d
		write(111, *)'nx = ', nx,  ' nz = ', nz
		write(111, *)'nx_with_apert = ', nx_with_apert, 'nz_with_apert = ', nz_with_apert
		write(111, *)'nnx = ', nnx, 'nnz = ', nnz
		write(111, *)'nx_shift = ', nx_shift, ' nz_shift = ', nz_shift
		write(111, *)'dx = ', dx, ' dz = ', dz
		write(111, *)'dt = ', dt
		
		close(111)
	endif

	isok=0

	end subroutine

!***********************************************************************



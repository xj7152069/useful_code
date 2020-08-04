
!***********************************************************************
	subroutine rtm_2d_single_shot_extrapolation_allocate(nnx, nnz, &
					nx_with_apert, nz_with_apert, &
					nx_bound_l, nx_bound_r, &
					nz_bound_u, nz_bound_d, &
					npad, order, ns, iflag_src)

	use module_rtm_extrapolation_2d
	implicit none
	integer		::	nnx, nnz
	integer		::	nx_with_apert, nz_with_apert
	integer		::	nx_bound_l, nx_bound_r
	integer		::	nz_bound_u, nz_bound_d
	integer		::	npad
	integer		::	order, ns
	integer		::	iflag_src
	integer		::	ierr
	

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

	allocate(ur1%u(-4:nnz+5, -4:nnx+5), stat=ierr)
	if(ierr.ne.0)then
		write(*,*)"Can not allocate working memory about ur1, stop!!!"
		stop
	endif

	allocate(ur2%u(-4:nnz+5, -4:nnx+5), stat=ierr)
	if(ierr.ne.0)then
		write(*,*)"Can not allocate working memory about ur2, stop!!!"
		stop
	endif

	allocate(ur3%u(-4:nnz+5, -4:nnx+5), stat=ierr)
	if(ierr.ne.0)then
		write(*,*)"Can not allocate working memory about ur3, stop!!!"
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


	if(iflag_src .eq. 1)then
		allocate(top%u(order, nx_with_apert, ns), stat=ierr)
		if(ierr.ne.0)then
			write(*,*)"Can not allocate working memory about top, stop!!!"
			stop
		endif

		allocate(bot%u(order, nx_with_apert, ns), stat=ierr)
		if(ierr.ne.0)then
			write(*,*)"Can not allocate working memory about bot, stop!!!"
			stop
		endif

		allocate(lef%u(order, nz_with_apert, ns), stat=ierr)
		if(ierr.ne.0)then
			write(*,*)"Can not allocate working memory about lef, stop!!!"
			stop
		endif

		allocate(rig%u(order, nz_with_apert, ns), stat=ierr)
		if(ierr.ne.0)then
			write(*,*)"Can not allocate working memory about rig, stop!!!"
			stop
		endif

	else if(iflag_src .eq. 0)then
		allocate(ust%u(-4:nnz+5, -4:nnx+5, ns), stat=ierr)
		if(ierr.ne.0)then
			write(*,*)"Can not allocate working memory about ust, stop!!!"
			stop
		endif

	endif

    allocate(image_tmp(nz_with_apert, nx_with_apert), stat=ierr)
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
	ur1%u=0.0
	ur2%u=0.0
	ur3%u=0.0

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

	if(iflag_src .eq. 1)then
		top%u=0.0
		bot%u=0.0
		lef%u=0.0
		rig%u=0.0
	else if(iflag_src .eq. 0)then
		ust%u=0.0
	endif

	image_tmp=0.0

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
	pur1 => ur1
	pur2 => ur2
	pur3 => ur3

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


	end subroutine

!***********************************************************************
	subroutine rtm_2d_single_shot_extrapolation_deallocate(iflag_src)

	use module_rtm_extrapolation_2d
	implicit none
	integer		::	iflag_src

	!deallocate wavefiled
	deallocate(us1%u)
	deallocate(us2%u)
	deallocate(us3%u)
	deallocate(ur1%u)
	deallocate(ur2%u)
	deallocate(ur3%u)

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

	if(iflag_src .eq. 1)then
		deallocate(top%u)
		deallocate(bot%u)
		deallocate(lef%u)
		deallocate(rig%u)
	else if(iflag_src .eq. 0)then
		deallocate(ust%u)
	endif

	deallocate(image_tmp)

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


	end subroutine

!***********************************************************************


!***********************************************************************
	subroutine born_2d_single_shot_extrapolation_allocate(nnx, nnz, &
					nx_with_apert, nz_with_apert, &
					nx_bound_l, nx_bound_r, &
					nz_bound_u, nz_bound_d, &
					npad, order, ns, iflag_src)

	use module_born_extrapolation_2d
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
	allocate(us01%u(-4:nnz+5, -4:nnx+5), stat=ierr)
	if(ierr.ne.0)then
		write(*,*)"Can not allocate working memory about us1, stop!!!"
		stop
	endif

	allocate(us02%u(-4:nnz+5, -4:nnx+5), stat=ierr)
	if(ierr.ne.0)then
		write(*,*)"Can not allocate working memory about us2, stop!!!"
		stop
	endif

	allocate(us03%u(-4:nnz+5, -4:nnx+5), stat=ierr)
	if(ierr.ne.0)then
		write(*,*)"Can not allocate working memory about us3, stop!!!"
		stop
	endif

	allocate(ur01%u(-4:nnz+5, -4:nnx+5), stat=ierr)
	if(ierr.ne.0)then
		write(*,*)"Can not allocate working memory about ur1, stop!!!"
		stop
	endif

	allocate(ur02%u(-4:nnz+5, -4:nnx+5), stat=ierr)
	if(ierr.ne.0)then
		write(*,*)"Can not allocate working memory about ur2, stop!!!"
		stop
	endif

	allocate(ur03%u(-4:nnz+5, -4:nnx+5), stat=ierr)
	if(ierr.ne.0)then
		write(*,*)"Can not allocate working memory about ur3, stop!!!"
		stop
	endif


	!allocate pml
    allocate(ux0_a1%u(nnz,nx_bound_l+nx_bound_r), stat=ierr)
	if(ierr.ne.0)then
		write(*,*)"Can not allocate working memory about ux_a1, stop!!!"
		stop
	endif

    allocate(ux0_a2%u(nnz,nx_bound_l+nx_bound_r), stat=ierr)
	if(ierr.ne.0)then
		write(*,*)"Can not allocate working memory about ux_a2, stop!!!"
		stop
	endif

    allocate(ux0_a3%u(nnz,nx_bound_l+nx_bound_r), stat=ierr)
	if(ierr.ne.0)then
		write(*,*)"Can not allocate working memory about ux_a3, stop!!!"
		stop
	endif

    allocate(ux0_b1%u(nnz,nx_bound_l+nx_bound_r), stat=ierr)
	if(ierr.ne.0)then
		write(*,*)"Can not allocate working memory about ux_b1, stop!!!"
		stop
	endif

    allocate(ux0_b2%u(nnz,nx_bound_l+nx_bound_r), stat=ierr)
	if(ierr.ne.0)then
		write(*,*)"Can not allocate working memory about ux_b2, stop!!!"
		stop
	endif

    allocate(ux0_bp1%u(nnz,nx_bound_l+nx_bound_r), stat=ierr)
	if(ierr.ne.0)then
		write(*,*)"Can not allocate working memory about ux_bp1, stop!!!"
		stop
	endif

    allocate(ux0_bp2%u(nnz,nx_bound_l+nx_bound_r), stat=ierr)
	if(ierr.ne.0)then
		write(*,*)"Can not allocate working memory about ux_bp2, stop!!!"
		stop
	endif

    allocate(ux0_bp3%u(nnz,nx_bound_l+nx_bound_r), stat=ierr)
	if(ierr.ne.0)then
		write(*,*)"Can not allocate working memory about ux_bp3, stop!!!"
		stop
	endif

    allocate(ux0_c1%u(nnz,nx_bound_l+nx_bound_r), stat=ierr)
	if(ierr.ne.0)then
		write(*,*)"Can not allocate working memory about ux_c1, stop!!!"
		stop
	endif

    allocate(ux0_c2%u(nnz,nx_bound_l+nx_bound_r), stat=ierr)
	if(ierr.ne.0)then
		write(*,*)"Can not allocate working memory about ux_c2, stop!!!"
		stop
	endif

    allocate(ux0_c3%u(nnz,nx_bound_l+nx_bound_r), stat=ierr)
	if(ierr.ne.0)then
		write(*,*)"Can not allocate working memory about ux_c3, stop!!!"
		stop
	endif


    !z direction
    allocate(uz0_a1%u(nz_bound_u+nz_bound_d,nnx), stat=ierr)
	if(ierr.ne.0)then
		write(*,*)"Can not allocate working memory about uz_a1, stop!!!"
		stop
	endif

    allocate(uz0_a2%u(nz_bound_u+nz_bound_d,nnx), stat=ierr)
	if(ierr.ne.0)then
		write(*,*)"Can not allocate working memory about uz_a2, stop!!!"
		stop
	endif

    allocate(uz0_a3%u(nz_bound_u+nz_bound_d,nnx), stat=ierr)
	if(ierr.ne.0)then
		write(*,*)"Can not allocate working memory about uz_a3, stop!!!"
		stop
	endif

    allocate(uz0_b1%u(nz_bound_u+nz_bound_d,nnx), stat=ierr)
	if(ierr.ne.0)then
		write(*,*)"Can not allocate working memory about uz_b1, stop!!!"
		stop
	endif

    allocate(uz0_b2%u(nz_bound_u+nz_bound_d,nnx), stat=ierr)
	if(ierr.ne.0)then
		write(*,*)"Can not allocate working memory about uz_b2, stop!!!"
		stop
	endif

    allocate(uz0_bp1%u(nz_bound_u+nz_bound_d,nnx), stat=ierr)
	if(ierr.ne.0)then
		write(*,*)"Can not allocate working memory about uz_bp1, stop!!!"
		stop
	endif

    allocate(uz0_bp2%u(nz_bound_u+nz_bound_d,nnx), stat=ierr)
	if(ierr.ne.0)then
		write(*,*)"Can not allocate working memory about uz_bp2, stop!!!"
		stop
	endif

    allocate(uz0_bp3%u(nz_bound_u+nz_bound_d,nnx), stat=ierr)
	if(ierr.ne.0)then
		write(*,*)"Can not allocate working memory about uz_bp3, stop!!!"
		stop
	endif

    allocate(uz0_c1%u(nz_bound_u+nz_bound_d,nnx), stat=ierr)
	if(ierr.ne.0)then
		write(*,*)"Can not allocate working memory about uz_c1, stop!!!"
		stop
	endif

    allocate(uz0_c2%u(nz_bound_u+nz_bound_d,nnx), stat=ierr)
	if(ierr.ne.0)then
		write(*,*)"Can not allocate working memory about uz_c2, stop!!!"
		stop
	endif

    allocate(uz0_c3%u(nz_bound_u+nz_bound_d,nnx), stat=ierr)
	if(ierr.ne.0)then
		write(*,*)"Can not allocate working memory about uz_c3, stop!!!"
		stop
	endif


	allocate(usr1%u(-4:nnz+5, -4:nnx+5), stat=ierr)
	if(ierr.ne.0)then
		write(*,*)"Can not allocate working memory about us1, stop!!!"
		stop
	endif

	allocate(usr2%u(-4:nnz+5, -4:nnx+5), stat=ierr)
	if(ierr.ne.0)then
		write(*,*)"Can not allocate working memory about us2, stop!!!"
		stop
	endif

	allocate(usr3%u(-4:nnz+5, -4:nnx+5), stat=ierr)
	if(ierr.ne.0)then
		write(*,*)"Can not allocate working memory about us3, stop!!!"
		stop
	endif

	allocate(urr1%u(-4:nnz+5, -4:nnx+5), stat=ierr)
	if(ierr.ne.0)then
		write(*,*)"Can not allocate working memory about ur1, stop!!!"
		stop
	endif

	allocate(urr2%u(-4:nnz+5, -4:nnx+5), stat=ierr)
	if(ierr.ne.0)then
		write(*,*)"Can not allocate working memory about ur2, stop!!!"
		stop
	endif

	allocate(urr3%u(-4:nnz+5, -4:nnx+5), stat=ierr)
	if(ierr.ne.0)then
		write(*,*)"Can not allocate working memory about ur3, stop!!!"
		stop
	endif


	!allocate pml
    allocate(uxr_a1%u(nnz,nx_bound_l+nx_bound_r), stat=ierr)
	if(ierr.ne.0)then
		write(*,*)"Can not allocate working memory about ux_a1, stop!!!"
		stop
	endif

    allocate(uxr_a2%u(nnz,nx_bound_l+nx_bound_r), stat=ierr)
	if(ierr.ne.0)then
		write(*,*)"Can not allocate working memory about ux_a2, stop!!!"
		stop
	endif

    allocate(uxr_a3%u(nnz,nx_bound_l+nx_bound_r), stat=ierr)
	if(ierr.ne.0)then
		write(*,*)"Can not allocate working memory about ux_a3, stop!!!"
		stop
	endif

    allocate(uxr_b1%u(nnz,nx_bound_l+nx_bound_r), stat=ierr)
	if(ierr.ne.0)then
		write(*,*)"Can not allocate working memory about ux_b1, stop!!!"
		stop
	endif

    allocate(uxr_b2%u(nnz,nx_bound_l+nx_bound_r), stat=ierr)
	if(ierr.ne.0)then
		write(*,*)"Can not allocate working memory about ux_b2, stop!!!"
		stop
	endif

    allocate(uxr_bp1%u(nnz,nx_bound_l+nx_bound_r), stat=ierr)
	if(ierr.ne.0)then
		write(*,*)"Can not allocate working memory about ux_bp1, stop!!!"
		stop
	endif

    allocate(uxr_bp2%u(nnz,nx_bound_l+nx_bound_r), stat=ierr)
	if(ierr.ne.0)then
		write(*,*)"Can not allocate working memory about ux_bp2, stop!!!"
		stop
	endif

    allocate(uxr_bp3%u(nnz,nx_bound_l+nx_bound_r), stat=ierr)
	if(ierr.ne.0)then
		write(*,*)"Can not allocate working memory about ux_bp3, stop!!!"
		stop
	endif

    allocate(uxr_c1%u(nnz,nx_bound_l+nx_bound_r), stat=ierr)
	if(ierr.ne.0)then
		write(*,*)"Can not allocate working memory about ux_c1, stop!!!"
		stop
	endif

    allocate(uxr_c2%u(nnz,nx_bound_l+nx_bound_r), stat=ierr)
	if(ierr.ne.0)then
		write(*,*)"Can not allocate working memory about ux_c2, stop!!!"
		stop
	endif

    allocate(uxr_c3%u(nnz,nx_bound_l+nx_bound_r), stat=ierr)
	if(ierr.ne.0)then
		write(*,*)"Can not allocate working memory about ux_c3, stop!!!"
		stop
	endif


    !z direction
    allocate(uzr_a1%u(nz_bound_u+nz_bound_d,nnx), stat=ierr)
	if(ierr.ne.0)then
		write(*,*)"Can not allocate working memory about uz_a1, stop!!!"
		stop
	endif

    allocate(uzr_a2%u(nz_bound_u+nz_bound_d,nnx), stat=ierr)
	if(ierr.ne.0)then
		write(*,*)"Can not allocate working memory about uz_a2, stop!!!"
		stop
	endif

    allocate(uzr_a3%u(nz_bound_u+nz_bound_d,nnx), stat=ierr)
	if(ierr.ne.0)then
		write(*,*)"Can not allocate working memory about uz_a3, stop!!!"
		stop
	endif

    allocate(uzr_b1%u(nz_bound_u+nz_bound_d,nnx), stat=ierr)
	if(ierr.ne.0)then
		write(*,*)"Can not allocate working memory about uz_b1, stop!!!"
		stop
	endif

    allocate(uzr_b2%u(nz_bound_u+nz_bound_d,nnx), stat=ierr)
	if(ierr.ne.0)then
		write(*,*)"Can not allocate working memory about uz_b2, stop!!!"
		stop
	endif

    allocate(uzr_bp1%u(nz_bound_u+nz_bound_d,nnx), stat=ierr)
	if(ierr.ne.0)then
		write(*,*)"Can not allocate working memory about uz_bp1, stop!!!"
		stop
	endif

    allocate(uzr_bp2%u(nz_bound_u+nz_bound_d,nnx), stat=ierr)
	if(ierr.ne.0)then
		write(*,*)"Can not allocate working memory about uz_bp2, stop!!!"
		stop
	endif

    allocate(uzr_bp3%u(nz_bound_u+nz_bound_d,nnx), stat=ierr)
	if(ierr.ne.0)then
		write(*,*)"Can not allocate working memory about uz_bp3, stop!!!"
		stop
	endif

    allocate(uzr_c1%u(nz_bound_u+nz_bound_d,nnx), stat=ierr)
	if(ierr.ne.0)then
		write(*,*)"Can not allocate working memory about uz_c1, stop!!!"
		stop
	endif

    allocate(uzr_c2%u(nz_bound_u+nz_bound_d,nnx), stat=ierr)
	if(ierr.ne.0)then
		write(*,*)"Can not allocate working memory about uz_c2, stop!!!"
		stop
	endif

    allocate(uzr_c3%u(nz_bound_u+nz_bound_d,nnx), stat=ierr)
	if(ierr.ne.0)then
		write(*,*)"Can not allocate working memory about uz_c3, stop!!!"
		stop
	endif


	if(iflag_src .eq. 1)then
		allocate(top0%u(order, nx_with_apert, ns), stat=ierr)
		if(ierr.ne.0)then
			write(*,*)"Can not allocate working memory about top, stop!!!"
			stop
		endif

		allocate(bot0%u(order, nx_with_apert, ns), stat=ierr)
		if(ierr.ne.0)then
			write(*,*)"Can not allocate working memory about bot, stop!!!"
			stop
		endif

		allocate(lef0%u(order, nz_with_apert, ns), stat=ierr)
		if(ierr.ne.0)then
			write(*,*)"Can not allocate working memory about lef, stop!!!"
			stop
		endif

		allocate(rig0%u(order, nz_with_apert, ns), stat=ierr)
		if(ierr.ne.0)then
			write(*,*)"Can not allocate working memory about rig, stop!!!"
			stop
		endif

		allocate(topr%u(order, nx_with_apert, ns), stat=ierr)
		if(ierr.ne.0)then
			write(*,*)"Can not allocate working memory about top, stop!!!"
			stop
		endif

		allocate(botr%u(order, nx_with_apert, ns), stat=ierr)
		if(ierr.ne.0)then
			write(*,*)"Can not allocate working memory about bot, stop!!!"
			stop
		endif

		allocate(lefr%u(order, nz_with_apert, ns), stat=ierr)
		if(ierr.ne.0)then
			write(*,*)"Can not allocate working memory about lef, stop!!!"
			stop
		endif

		allocate(rigr%u(order, nz_with_apert, ns), stat=ierr)
		if(ierr.ne.0)then
			write(*,*)"Can not allocate working memory about rig, stop!!!"
			stop
		endif

	else if(iflag_src .eq. 0)then
		allocate(ust0%u(-4:nnz+5, -4:nnx+5, ns), stat=ierr)
		if(ierr.ne.0)then
			write(*,*)"Can not allocate working memory about ust, stop!!!"
			stop
		endif

		allocate(ustr%u(-4:nnz+5, -4:nnx+5, ns), stat=ierr)
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
	us01%u=0.0
	us02%u=0.0
	us03%u=0.0
	ur01%u=0.0
	ur02%u=0.0
	ur03%u=0.0

	!x direction
    ux0_a1%u=0.0
    ux0_a2%u=0.0
    ux0_a3%u=0.0
    ux0_b1%u=0.0
    ux0_b2%u=0.0
    ux0_bp1%u=0.0
    ux0_bp2%u=0.0
    ux0_bp3%u=0.0
    ux0_c1%u=0.0
    ux0_c2%u=0.0
    ux0_c3%u=0.0

    !z direction
    uz0_a1%u=0.0
    uz0_a2%u=0.0
    uz0_a3%u=0.0
    uz0_b1%u=0.0
    uz0_b2%u=0.0
    uz0_bp1%u=0.0
    uz0_bp2%u=0.0
    uz0_bp3%u=0.0
    uz0_c1%u=0.0
    uz0_c2%u=0.0
    uz0_c3%u=0.0

	!wavefields
	usr1%u=0.0
	usr2%u=0.0
	usr3%u=0.0
	urr1%u=0.0
	urr2%u=0.0
	urr3%u=0.0

	!x direction
    uxr_a1%u=0.0
    uxr_a2%u=0.0
    uxr_a3%u=0.0
    uxr_b1%u=0.0
    uxr_b2%u=0.0
    uxr_bp1%u=0.0
    uxr_bp2%u=0.0
    uxr_bp3%u=0.0
    uxr_c1%u=0.0
    uxr_c2%u=0.0
    uxr_c3%u=0.0

    !z direction
    uzr_a1%u=0.0
    uzr_a2%u=0.0
    uzr_a3%u=0.0
    uzr_b1%u=0.0
    uzr_b2%u=0.0
    uzr_bp1%u=0.0
    uzr_bp2%u=0.0
    uzr_bp3%u=0.0
    uzr_c1%u=0.0
    uzr_c2%u=0.0
    uzr_c3%u=0.0

	if(iflag_src .eq. 1)then
		top0%u=0.0
		bot0%u=0.0
		lef0%u=0.0
		rig0%u=0.0
		topr%u=0.0
		botr%u=0.0
		lefr%u=0.0
		rigr%u=0.0
	else if(iflag_src .eq. 0)then
		ust0%u=0.0
		ustr%u=0.0
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
	pus01 => us01
	pus02 => us02
	pus03 => us03
	pur01 => ur01
	pur02 => ur02
	pur03 => ur03

	!x direction
    pux0_a1 => ux0_a1
    pux0_a2 => ux0_a2
    pux0_a3 => ux0_a3
    pux0_b1 => ux0_b1
    pux0_b2 => ux0_b2
    pux0_bp1 => ux0_bp1
    pux0_bp2 => ux0_bp2
    pux0_bp3 => ux0_bp3
    pux0_c1 => ux0_c1
    pux0_c2 => ux0_c2
    pux0_c3 => ux0_c3

    !z direction
    puz0_a1 => uz0_a1
    puz0_a2 => uz0_a2
    puz0_a3 => uz0_a3
    puz0_b1 => uz0_b1
    puz0_b2 => uz0_b2
    puz0_bp1 => uz0_bp1
    puz0_bp2 => uz0_bp2
    puz0_bp3 => uz0_bp3
    puz0_c1 => uz0_c1
    puz0_c2 => uz0_c2
    puz0_c3 => uz0_c3

	!wavefields
	pusr1 => usr1
	pusr2 => usr2
	pusr3 => usr3
	purr1 => urr1
	purr2 => urr2
	purr3 => urr3

	!x direction
    puxr_a1 => uxr_a1
    puxr_a2 => uxr_a2
    puxr_a3 => uxr_a3
    puxr_b1 => uxr_b1
    puxr_b2 => uxr_b2
    puxr_bp1 => uxr_bp1
    puxr_bp2 => uxr_bp2
    puxr_bp3 => uxr_bp3
    puxr_c1 => uxr_c1
    puxr_c2 => uxr_c2
    puxr_c3 => uxr_c3

    !z direction
    puzr_a1 => uzr_a1
    puzr_a2 => uzr_a2
    puzr_a3 => uzr_a3
    puzr_b1 => uzr_b1
    puzr_b2 => uzr_b2
    puzr_bp1 => uzr_bp1
    puzr_bp2 => uzr_bp2
    puzr_bp3 => uzr_bp3
    puzr_c1 => uzr_c1
    puzr_c2 => uzr_c2
    puzr_c3 => uzr_c3


	end subroutine

!***********************************************************************
	subroutine born_2d_single_shot_extrapolation_deallocate(iflag_src)

	use module_born_extrapolation_2d
	implicit none
	integer		::	iflag_src

	!deallocate wavefiled
	deallocate(us01%u)
	deallocate(us02%u)
	deallocate(us03%u)
	deallocate(ur01%u)
	deallocate(ur02%u)
	deallocate(ur03%u)

	!deallocate pml
    !x direction
    deallocate(ux0_a1%u)
    deallocate(ux0_a2%u)
    deallocate(ux0_a3%u)
    deallocate(ux0_b1%u)
    deallocate(ux0_b2%u)
    deallocate(ux0_bp1%u)
    deallocate(ux0_bp2%u)
    deallocate(ux0_bp3%u)
    deallocate(ux0_c1%u)
    deallocate(ux0_c2%u)
    deallocate(ux0_c3%u)

    !z direction
    deallocate(uz0_a1%u)
    deallocate(uz0_a2%u)
    deallocate(uz0_a3%u)
    deallocate(uz0_b1%u)
    deallocate(uz0_b2%u)
    deallocate(uz0_bp1%u)
    deallocate(uz0_bp2%u)
    deallocate(uz0_bp3%u)
    deallocate(uz0_c1%u)
    deallocate(uz0_c2%u)
    deallocate(uz0_c3%u)

	!deallocate wavefiled
	deallocate(usr1%u)
	deallocate(usr2%u)
	deallocate(usr3%u)
	deallocate(urr1%u)
	deallocate(urr2%u)
	deallocate(urr3%u)

	!deallocate pml
    !x direction
    deallocate(uxr_a1%u)
    deallocate(uxr_a2%u)
    deallocate(uxr_a3%u)
    deallocate(uxr_b1%u)
    deallocate(uxr_b2%u)
    deallocate(uxr_bp1%u)
    deallocate(uxr_bp2%u)
    deallocate(uxr_bp3%u)
    deallocate(uxr_c1%u)
    deallocate(uxr_c2%u)
    deallocate(uxr_c3%u)

    !z direction
    deallocate(uzr_a1%u)
    deallocate(uzr_a2%u)
    deallocate(uzr_a3%u)
    deallocate(uzr_b1%u)
    deallocate(uzr_b2%u)
    deallocate(uzr_bp1%u)
    deallocate(uzr_bp2%u)
    deallocate(uzr_bp3%u)
    deallocate(uzr_c1%u)
    deallocate(uzr_c2%u)
    deallocate(uzr_c3%u)

	if(iflag_src .eq. 1)then
		deallocate(top0%u)
		deallocate(bot0%u)
		deallocate(lef0%u)
		deallocate(rig0%u)
		deallocate(topr%u)
		deallocate(botr%u)
		deallocate(lefr%u)
		deallocate(rigr%u)
	else if(iflag_src .eq. 0)then
		deallocate(ust0%u)
		deallocate(ustr%u)
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

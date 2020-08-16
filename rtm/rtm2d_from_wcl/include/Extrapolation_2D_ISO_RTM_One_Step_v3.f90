!*=======================================================================*
!*版本说明
!*2020.04.28
!*2D声波方程波场外推，采用one_step方法实现；
!*采用openMp方式实现；
!*
!*
!*subroutine  Extrapolation_Add_PML_one_step(vv, x1, x2, x3, &
!*			ux_a1, ux_a2, ux_a3, &
!*			ux_b1, ux_b2, ux_bp1, ux_bp2, ux_bp3, &
!*			ux_c1, ux_c2, ux_c3, &
!*			uz_a1, uz_a2, uz_a3, &
!*			uz_b1, uz_b2, uz_bp1, uz_bp2, uz_bp3, &
!*			uz_c1, uz_c2, uz_c3, &
!*			coe, nnx, nnz, nx, nz, &
!*			nx_bound_l, nx_bound_r, &
!*			nz_bound_u, nz_bound_d, &
!*			dx, dz, dt, r, it)
!*函数含义：带吸收边界的波场传播
!*参数含义：
!*		vv：输入的速度场
!*		x1: 输入的t-1时刻的波场
!*		x2: 输入的t时刻的波场
!*		x3: 输出的t+1时刻的波场
!*		X方向的PML方程参数：
!*			ux_a1, ux_a2, ux_a3, &	
!*			ux_b1, ux_b2, ux_bp1, ux_bp2, ux_bp3, &
!*			ux_c1, ux_c2, ux_c3, &
!*		Z方向的PML方程参数：
!*			uz_a1, uz_a2, uz_a3, &
!*			uz_b1, uz_b2, uz_bp1, uz_bp2, uz_bp3, &
!*			uz_c1, uz_c2, uz_c3, &
!*		coe: 差分系数
!*		nnx: X方向的总的网格点数
!*		nnz: Z方向的总的网格点数
!*		nx_bound_l: X方向左侧的PML网格点数
!*		nx_bound_r: X方向右侧的PML网格点数
!*		nz_bound_u: Z方向上侧的PML网格点数
!*		nz_bound_d: Z方向下侧的PML网格点数
!*		dx: X方向的采样间隔
!*		dz: Z方向的采样间隔
!*		dt: 时间采样间隔
!*		r: PML吸收参数，默认为1.0e-3
!***********************************************************************
    subroutine  Extrapolation_2D_Add_PML_One_Step(vv, u1, u2, u3, &
			ux_a1, ux_a2, ux_a3, &
			ux_b1, ux_b2, ux_bp1, ux_bp2, ux_bp3, &
			ux_c1, ux_c2, ux_c3, &
			uz_a1, uz_a2, uz_a3, &
			uz_b1, uz_b2, uz_bp1, uz_bp2, uz_bp3, &
			uz_c1, uz_c2, uz_c3, &
			funa, dfuna, &
			line_xl, line_xr, line_zu, line_zd, &
			coe_2nd_dx2, coe_2nd_dz2, &
			coe_1st_dx, coe_1st_dz, &
			npad, nnx, nnz, nx, nz, &
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
	integer	::	npad
    real    ::	dx, dz, dt
    real    ::	r
	integer ::	it
  
    real	::	vv(nnz, nnx)
    real	::	u1(-npad+1:nnz+npad, -npad+1:nnx+npad)
    real	::	u2(-npad+1:nnz+npad, -npad+1:nnx+npad)
    real	::	u3(-npad+1:nnz+npad, -npad+1:nnx+npad)

    !wavefield for pml field
    !forward
    !x direction
    real	::	ux_a1(nnz,nx_bound_l+nx_bound_r)
    real	::	ux_a2(nnz,nx_bound_l+nx_bound_r)
    real	::	ux_a3(nnz,nx_bound_l+nx_bound_r)
    real	::	ux_b1(nnz,nx_bound_l+nx_bound_r)
    real	::	ux_b2(nnz,nx_bound_l+nx_bound_r)
    real	::	ux_bp1(nnz,nx_bound_l+nx_bound_r)
    real	::	ux_bp2(nnz,nx_bound_l+nx_bound_r)
    real	::	ux_bp3(nnz,nx_bound_l+nx_bound_r)
    real	::	ux_c1(nnz,nx_bound_l+nx_bound_r)
    real	::	ux_c2(nnz,nx_bound_l+nx_bound_r)
    real	::	ux_c3(nnz,nx_bound_l+nx_bound_r)

    !z direction
    real	::	uz_a1(nz_bound_u+nz_bound_d,nnx)
    real	::	uz_a2(nz_bound_u+nz_bound_d,nnx)
    real	::	uz_a3(nz_bound_u+nz_bound_d,nnx)
    real	::	uz_b1(nz_bound_u+nz_bound_d,nnx)
    real	::	uz_b2(nz_bound_u+nz_bound_d,nnx)
    real	::	uz_bp1(nz_bound_u+nz_bound_d,nnx)
    real	::	uz_bp2(nz_bound_u+nz_bound_d,nnx)
    real	::	uz_bp3(nz_bound_u+nz_bound_d,nnx)
    real	::	uz_c1(nz_bound_u+nz_bound_d,nnx)
    real	::	uz_c2(nz_bound_u+nz_bound_d,nnx)
    real	::	uz_c3(nz_bound_u+nz_bound_d,nnx)

    real	::	funa(nnz, nnx)
    real	::	dfuna(nnz, nnx)
	real	::	coe_2nd_dx2(npad), coe_2nd_dz2(npad)
	real	::	coe_1st_dx(npad), coe_1st_dz(npad)
	integer	::	line_xl(nx_bound_l, 2), line_xr(nx_bound_r, 2)
	integer	::	line_zu(nz_bound_u, 2), line_zd(nz_bound_d, 2)


    !Local variables
	integer	::	ii
    integer ::	ix, iz
    integer ::	inx, inz
    integer ::	iix, iiz
    real    ::	v2, dt2
	real	::	x1_deri, z1_deri
	real	::	x2_deri, z2_deri
	real	::	a, da

      
	dt2=dt*dt


!$omp parallel private(ix, iz, iix, iiz, inx, inz, ii, &
!$omp  		v2, z2_deri, x2_deri, z1_deri, x1_deri, a, da)

!$omp do 

    do ix=1, nnx
        do iz=1, nnz
			if(iz>nz_bound_u.and.iz<nnz-nz_bound_d+1.and.&
				ix>nx_bound_l.and.ix<nnx-nx_bound_r+1)then

			inx=ix
			inz=iz
			v2=vv(inz, inx)*vv(inz, inx)

			z2_deri = 0.0
			x2_deri = 0.0
			do ii = 1, npad
				z2_deri = z2_deri + coe_2nd_dz2(ii)*(u2(inz+ii, inx) + u2(inz-ii, inx) - 2.0*u2(inz, inx))
				x2_deri = x2_deri + coe_2nd_dx2(ii)*(u2(inz, inx+ii) + u2(inz, inx-ii) - 2.0*u2(inz, inx))
			enddo

			u3(inz, inx) = 2.0*u2(inz, inx) - u1(inz, inx) + dt2*v2*(z2_deri+x2_deri)
			
			else if(ix<=nx_bound_l)then

			iix = ix
			inx = ix

			iiz = iz
			inz = iz
			v2=vv(inz, inx)*vv(inz, inx)

			z2_deri = 0.0
			x2_deri = 0.0
			x1_deri = 0.0
			do ii = 1, npad
				z2_deri = z2_deri + coe_2nd_dz2(ii)*(u2(inz+ii, inx) + u2(inz-ii, inx) - 2.0*u2(inz, inx))
				x2_deri = x2_deri + coe_2nd_dx2(ii)*(u2(inz, inx+ii) + u2(inz, inx-ii) - 2.0*u2(inz, inx))
				x1_deri = x1_deri + coe_1st_dx(ii) *(u2(inz, inx+ii) - u2(inz, inx-ii)) 
			enddo


			a = funa(inz, inx)
			da = dfuna(inz, inx)

			!equation 1
			ux_a3(iiz, iix) = 2.0*ux_a2(iiz, iix) - ux_a1(iiz, iix) + &
				dt2*(v2*x2_deri - 2.0*a*(ux_a2(iiz, iix) - ux_a1(iiz, iix))/dt - a*a*ux_a2(iiz, iix)) 
			!equation 2
			ux_bp3(iiz, iix) = 2.0*ux_bp2(iiz, iix) - ux_bp1(iiz, iix) + &
				dt2*(-1.0*v2*da*x1_deri - 2.0*a*(ux_bp2(iiz, iix) - ux_bp1(iiz, iix))/dt - a*a*ux_bp2(iiz, iix)) 
			ux_b2(iiz, iix) = ux_b1(iiz, iix) + dt*(ux_bp2(iiz, iix) - a*ux_b1(iiz, iix))
			!equation 3
			ux_c3(iiz, iix) = 2.0*ux_c2(iiz, iix) - ux_c1(iiz, iix) + dt2*v2*z2_deri
              
			u3(inz, inx) = ux_a3(iiz, iix) + ux_b2(iiz, iix) + ux_c3(iiz, iix)


			!x_positive  
			else if(ix>=nnx-nx_bound_r+1)then

			iix = nx_bound_l + ix - (nnx-nx_bound_r)
			inx = ix

			iiz = iz
			inz = iz
			v2=vv(inz, inx)*vv(inz, inx)

			z2_deri = 0.0
			x2_deri = 0.0
			x1_deri = 0.0
			do ii = 1, npad
				z2_deri = z2_deri + coe_2nd_dz2(ii)*(u2(inz+ii, inx) + u2(inz-ii, inx) - 2.0*u2(inz, inx))
				x2_deri = x2_deri + coe_2nd_dx2(ii)*(u2(inz, inx+ii) + u2(inz, inx-ii) - 2.0*u2(inz, inx))
				x1_deri = x1_deri + coe_1st_dx(ii) *(u2(inz, inx+ii) - u2(inz, inx-ii)) 
			enddo

			a = funa(inz, inx)
			da = dfuna(inz, inx)

			!equation 1
			ux_a3(iiz, iix) = 2.0*ux_a2(iiz, iix) - ux_a1(iiz, iix) + &
				dt2*(v2*x2_deri - 2.0*a*(ux_a2(iiz, iix) - ux_a1(iiz, iix))/dt - a*a*ux_a2(iiz, iix)) 
			!equation 2
			ux_bp3(iiz, iix) = 2.0*ux_bp2(iiz, iix) - ux_bp1(iiz, iix) + &
				dt2*(-1.0*v2*da*x1_deri - 2.0*a*(ux_bp2(iiz, iix) - ux_bp1(iiz, iix))/dt - a*a*ux_bp2(iiz, iix)) 
			ux_b2(iiz, iix) = ux_b1(iiz, iix) + dt*(ux_bp2(iiz, iix) - a*ux_b1(iiz, iix))
			!equation 3
			ux_c3(iiz, iix) = 2.0*ux_c2(iiz, iix) - ux_c1(iiz, iix) + dt2*v2*z2_deri
              
			u3(inz, inx) = ux_a3(iiz, iix) + ux_b2(iiz, iix) + ux_c3(iiz, iix)

			!z_negative
			else if(iz<=nz_bound_u)then

			iiz = iz
			inz = iz

			iix = ix
			inx = ix
			v2=vv(inz, inx)*vv(inz, inx)
			z2_deri = 0.0
			x2_deri = 0.0
			z1_deri = 0.0
			do ii = 1, npad
				z2_deri = z2_deri + coe_2nd_dz2(ii)*(u2(inz+ii, inx) + u2(inz-ii, inx) - 2.0*u2(inz, inx))
				x2_deri = x2_deri + coe_2nd_dx2(ii)*(u2(inz, inx+ii) + u2(inz, inx-ii) - 2.0*u2(inz, inx))
				z1_deri = z1_deri + coe_1st_dz(ii) *(u2(inz+ii, inx) - u2(inz-ii, inx)) 
			enddo

			a = funa(inz, inx)
			da = dfuna(inz, inx)

			!equation 1
			uz_a3(iiz, iix) = 2.0*uz_a2(iiz, iix) - uz_a1(iiz, iix) + &
				dt2*(v2*z2_deri - 2.0*a*(uz_a2(iiz, iix) - uz_a1(iiz, iix))/dt - a*a*uz_a2(iiz, iix)) 
			!equation 2
			uz_bp3(iiz, iix) = 2.0*uz_bp2(iiz, iix) - uz_bp1(iiz, iix) + &
				dt2*(-1.0*v2*da*z1_deri - 2.0*a*(uz_bp2(iiz, iix) - uz_bp1(iiz, iix))/dt - a*a*uz_bp2(iiz, iix)) 
			uz_b2(iiz, iix) = uz_b1(iiz, iix) + dt*(uz_bp2(iiz, iix) - a*uz_b1(iiz, iix))
			!equation 3
			uz_c3(iiz, iix) = 2.0*uz_c2(iiz, iix) - uz_c1(iiz, iix) + dt2*v2*x2_deri
              
			u3(inz, inx) = uz_a3(iiz, iix) + uz_b2(iiz, iix) + uz_c3(iiz, iix)


			!z_positive
			else if(iz>=nnz-nz_bound_d+1)then

			iiz = nz_bound_u+iz-(nnz-nz_bound_d)
			inz = iz

			iix = ix
			inx = ix
			v2=vv(inz, inx)*vv(inz, inx)
			z2_deri = 0.0
			x2_deri = 0.0
			z1_deri = 0.0
			do ii = 1, npad
				z2_deri = z2_deri + coe_2nd_dz2(ii)*(u2(inz+ii, inx) + u2(inz-ii, inx) - 2.0*u2(inz, inx))
				x2_deri = x2_deri + coe_2nd_dx2(ii)*(u2(inz, inx+ii) + u2(inz, inx-ii) - 2.0*u2(inz, inx))
				z1_deri = z1_deri + coe_1st_dz(ii) *(u2(inz+ii, inx) - u2(inz-ii, inx)) 
			enddo

			a = funa(inz, inx)
			da = dfuna(inz, inx)

			!equation 1
			uz_a3(iiz, iix) = 2.0*uz_a2(iiz, iix) - uz_a1(iiz, iix) + &
				dt2*(v2*z2_deri - 2.0*a*(uz_a2(iiz, iix) - uz_a1(iiz, iix))/dt - a*a*uz_a2(iiz, iix)) 
			!equation 2
			uz_bp3(iiz, iix) = 2.0*uz_bp2(iiz, iix) - uz_bp1(iiz, iix) + &
				dt2*(-1.0*v2*da*z1_deri - 2.0*a*(uz_bp2(iiz, iix) - uz_bp1(iiz, iix))/dt - a*a*uz_bp2(iiz, iix)) 
			uz_b2(iiz, iix) = uz_b1(iiz, iix) + dt*(uz_bp2(iiz, iix) - a*uz_b1(iiz, iix))
			!equation 3
			uz_c3(iiz, iix) = 2.0*uz_c2(iiz, iix) - uz_c1(iiz, iix) + dt2*v2*x2_deri
			u3(inz, inx) = uz_a3(iiz, iix) + uz_b2(iiz, iix) + uz_c3(iiz, iix)

			endif
		enddo
	enddo
!$omp end do 
!$omp end parallel



    return
    end subroutine

!***********************************************************************

	subroutine absorbing_function_2d(funa, dfuna, vv, nnz, nnx, dx, dz, &
				line_xl, line_xr, line_zu, line_zd, &
				nx_bound_l, nx_bound_r, &
				nz_bound_u, nz_bound_d, &
				r)

    implicit none
	!Data dictionary: declare variable types, definitons, & units
    !dummy variables
    integer	::	nnx, nnz
	integer ::	nx_bound_l, nx_bound_r
    integer ::	nz_bound_u, nz_bound_d
    real	::	dx, dz
    real	::	vv(nnz, nnx)
    real	::	funa(nnz, nnx)
    real	::	dfuna(nnz, nnx)
	integer	::	line_xl(nx_bound_l, 2), line_xr(nx_bound_r, 2)
	integer	::	line_zu(nz_bound_u, 2), line_zd(nz_bound_d, 2)
	real	::	r


    !local variables
    integer	::	ix, iz
    integer	::	iix, iiz
    integer	::	inx, inz
	real	::	tmp
	real	::	dis, dis3
	real	::	x, z
	real	::	kzuxl, kzdxl, kzuxr, kzdxr
	real	::	kxlzu, kxrzu, kxlzd, kxrzd


	!*=================================================================
	!*          now do the wavefield extrapolation
	!*=================================================================
	!left
	tmp=3.0/2.0*log(1.0/r)

	kzuxl=(dz*(nz_bound_u-1))/(dx*(nx_bound_l-1))
	kzdxl=(dz*(nz_bound_d-1))/(dx*(nx_bound_l-1))
	kzuxr=(dz*(nz_bound_u-1))/(dx*(nx_bound_r-1))
	kzdxr=(dz*(nz_bound_d-1))/(dx*(nx_bound_r-1))
	kxlzu=(dx*(nx_bound_l-1))/(dz*(nz_bound_u-1))
	kxrzu=(dx*(nx_bound_r-1))/(dz*(nz_bound_u-1))
	kxlzd=(dx*(nx_bound_l-1))/(dz*(nz_bound_d-1))
	kxrzd=(dx*(nx_bound_r-1))/(dz*(nz_bound_d-1))


	!left
	do ix=1, nx_bound_l
		line_xl(ix,1) = int((ix-1)*kzuxl)
		line_xl(ix,2) = int(nnz-(ix-1)*kzdxl+1)
		if(line_xl(ix, 1) .lt. 1)	line_xl(ix, 1) = 1
		if(line_xl(ix, 2) .gt. nnz)	line_xl(ix, 2) = nnz
	enddo

	!right
	do ix=1, nx_bound_r
		line_xr(ix, 1) = int(nz_bound_u-(ix-1)*kzuxr-1) 
		line_xr(ix, 2) = int(nnz-nz_bound_d+(ix-1)*kzdxr+1)
		if(line_xr(ix, 1) .lt. 1)	line_xr(ix, 1) = 1
		if(line_xr(ix, 2) .gt. nnz)	line_xr(ix, 2) = nnz
	enddo

	!up
	do iz=1, nz_bound_u
		line_zu(iz, 1) = int((iz-1)*kxlzu)
		line_zu(iz, 2) = int(nnx-(iz-1)*kxrzu+1)
		if(line_zu(iz, 1) .lt. 1)	line_zu(iz, 1) = 1
		if(line_zu(iz, 2) .gt. nnx)	line_zu(iz, 2) = nnx
	enddo

	!down
	do iz=1, nz_bound_d
		line_zd(iz, 1) = int(nx_bound_l-(iz-1)*kxlzd-1)
		line_zd(iz, 2) = int(nnx-nx_bound_r +(iz-1)*kxrzd+1)
		if(line_zd(iz, 1) .lt. 1)	line_zd(iz, 1) = 1
		if(line_zd(iz, 2) .gt. nnx)	line_zd(iz ,2) = nnx
	enddo
		



	!pml function d(x)  and d'(x)
	!left
	dis=nx_bound_l*dx
	dis3=dis*dis*dis
	do ix=1, nx_bound_l
		inx=ix
		do inz=line_xl(ix, 1), line_xl(ix, 2)
			x=(nx_bound_l-ix+1)*dx
			funa(inz, inx)=tmp*vv(inz, inx)*x*x/dis3
			dfuna(inz, inx)=tmp*2.0*vv(inz, inx)*x/dis3
		enddo	
	enddo

	!right
	dis=nx_bound_r*dx
	dis3=dis*dis*dis
	do ix=1, nx_bound_r
		inx= nnx-nx_bound_r+ix
		do inz=line_xr(ix, 1), line_xr(ix, 2) 
			x=ix*dx
			funa(inz, inx)=tmp*vv(inz, inx)*x*x/dis3
			dfuna(inz, inx)=tmp*2.0*vv(inz, inx)*x/dis3
		enddo	
	enddo

	!top
	dis=nz_bound_u*dz
	dis3=dis*dis*dis
	do iz=1, nz_bound_u
		inz=iz
		do inx=line_zu(iz, 1), line_zu(iz, 2)
			z=(nz_bound_u-iz+1)*dz
			funa(inz, inx)=tmp*vv(inz, inx)*z*z/dis3
			dfuna(inz, inx)=tmp*2.0*vv(inz, inx)*z/dis3
		enddo	
	enddo

	!down
	dis=nz_bound_d*dz
	dis3=dis*dis*dis
	do iz=1, nz_bound_d
		inz= nnz-nz_bound_d+iz
		do inx=line_zd(iz, 1), line_zd(iz, 2)
			z=iz*dz
			funa(inz, inx)=tmp*vv(inz, inx)*z*z/dis3
			dfuna(inz, inx)=tmp*2.0*vv(inz, inx)*z/dis3
		enddo
	enddo

 	do inx=nx_bound_l+1, nnx-nx_bound_r
        do inz=nz_bound_u+1, nnz-nz_bound_d
			funa(inz, inx)=0.0
			dfuna(inz, inx)=0.0
		enddo
    enddo


	print*, 'here. debug'

    end subroutine

!***********************************************************************








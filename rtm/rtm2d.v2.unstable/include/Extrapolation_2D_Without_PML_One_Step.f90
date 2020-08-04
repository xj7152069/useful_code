!*=======================================================================*
!*版本说明
!*2020.04.28
!*2D声波方程波场外推，采用one_step方法实现；
!*采用openMp方式实现；
!*
!*包含的子函数：
!*subroutine Extrapolation_Without_PML_one_step(vv, x1, x2, x3, coe, nnx, nnz,&
!*				nx_bound_l,nx_bound_r, nz_bound_u,nz_bound_d,&
!*				dx, dz, dt)
!*函数含义：无吸收边界的波场传播
!*参数含义：
!*		vv：输入的速度场
!*		x1: 输入的t-1时刻的波场
!*		x2: 输入的t时刻的波场
!*		x3: 输出的t+1时刻的波场
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
!*

	subroutine Extrapolation_2D_Without_PML_One_Step(vv, x1, x2, x3, coe, nnx, nnz,&
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
    real	::	x3(-4:nnz+5, -4:nnx+5)


    !local variables
    integer	::	ix, iz
	real	::	ddt, vv2

	!*=================================================================
	!*          now do the wavefield extrapolation
	!*=================================================================
    ddt=dt
    
!$omp parallel private(vv2) 
!$omp do

 	do ix=nx_bound_l+1+5, nnx-nx_bound_r-5
        do iz=nz_bound_u+1+5, nnz-nz_bound_d-5

			vv2=vv(iz, ix)*vv(iz, ix)

			x3(iz,ix) = coe(1)*x2(iz, ix)&
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

!*************************************************************************


!*=======================================================================*
!*版本说明
!*2020.04.28
!*2D声波方程波场外推，采用one_step方法实现；
!*采用openMp方式实现；
!*
!*
!*包含的子函数：
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
!*************************************************************************
    subroutine  Extrapolation_2D_Add_PML_One_Step_old(vv, x1, x2, x3, &
			ux_a1, ux_a2, ux_a3, &
			ux_b1, ux_b2, ux_bp1, ux_bp2, ux_bp3, &
			ux_c1, ux_c2, ux_c3, &
			uz_a1, uz_a2, uz_a3, &
			uz_b1, uz_b2, uz_bp1, uz_bp2, uz_bp3, &
			uz_c1, uz_c2, uz_c3, &
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
    real	::	x3(-4:nnz+5, -4:nnx+5)
    real	::	coe(13)

    !wavefield within pml field  


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
             

			!absorbing boundary condition
			!x_negative
			else if(ix<=nx_bound_l)then
				axis_x_negative_pml=real(ix-(nx_bound_l+1))*dx
				axis_x_thickness_pml=real(nx_bound_l)*dx
				grid_x_negative_pml=ix
              
				d_attenuation_x_present=(3.0*vv(iz,ix)/(2.0*&
					axis_x_thickness_pml))*(axis_x_negative_pml/&
					axis_x_thickness_pml)**2*log(1.0/r)

				deri_d_attenuation_x=(3.0*vv(iz,ix)/(2.0*&
					(axis_x_thickness_pml)))*(2.0*axis_x_negative_pml/&
					axis_x_thickness_pml**2)*log(1.0/r)


				!pay attention to the subsripts of x2 etc.
				ux_a3(iz,ix)=(ux_a2(iz,ix)*(2.0*dx**2+2.0*&
					d_attenuation_x_present*ddt*dx**2-&
					d_attenuation_x_present**2*dx**2*ddt**2)-dx**2*&
					ux_a1(iz,ix)+vv(iz,ix)**2*ddt**2*&
					(x2(iz,ix+1)+x2(iz,ix-1)-2.0*x2(iz,ix)))/&
					(dx**2+2.0*d_attenuation_x_present*ddt*dx**2)


				ux_bp3(iz,ix)=(ux_bp2(iz,ix)*(2.0*dx+2.0*&
					d_attenuation_x_present*ddt*dx-ddt**2*dx*&
					d_attenuation_x_present**2)-dx*ux_bp1(iz,ix)&
					-vv(iz,ix)**2*ddt**2*deri_d_attenuation_x*&
					(x2(iz,ix+1)-x2(iz,ix)))/&
					(dx+2.0*d_attenuation_x_present*ddt*dx)

              
				ux_b2(iz,ix)=ux_b1(iz,ix)*&
					(1-d_attenuation_x_present*ddt)&
					+ux_bp3(iz,ix)*ddt
              
				ux_c3(iz,ix)=2.0*ux_c2(iz,ix)-ux_c1(iz,ix)+ &
					vv(iz,ix)**2*ddt**2*(x2(iz+1,ix)+x2(iz-1,ix)-2.0*x2(iz,ix))/(dz**2)
              
				x3(iz,ix)=ux_a3(iz,ix)+ux_b2(iz,ix)+ux_c3(iz,ix)

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
				ux_a3(iz,grid_x_positive_pml)=&
					(ux_a2(iz,grid_x_positive_pml)*(2*dx**2+2*&
					d_attenuation_x_present*ddt*dx**2-&
					d_attenuation_x_present**2*dx**2*ddt**2)-dx**2*&
					ux_a1(iz,grid_x_positive_pml)+vv(iz,ix)**2*&
					ddt**2*(x2(iz,ix+1)+&
					x2(iz,ix-1)-2*x2(iz,ix)))/&
					(dx**2+2*d_attenuation_x_present*ddt*dx**2)
              
				ux_bp3(iz,grid_x_positive_pml)=&
					(ux_bp2(iz,grid_x_positive_pml)*&
					(2*dx+2*d_attenuation_x_present*ddt*dx-ddt**2*dx*&
					d_attenuation_x_present**2)-dx*&
					ux_bp1(iz,grid_x_positive_pml)&
					-vv(iz,ix)**2*ddt**2*deri_d_attenuation_x*&
					(x2(iz,ix+1)-x2(iz,ix)))/&
					(dx+2*d_attenuation_x_present*ddt*dx)
              
				ux_b2(iz,grid_x_positive_pml)=&
					ux_b1(iz,grid_x_positive_pml)*&
					(1-d_attenuation_x_present*ddt)&
					+ux_bp3(iz,grid_x_positive_pml)*ddt
              
				ux_c3(iz,grid_x_positive_pml)=&
					2*ux_c2(iz,grid_x_positive_pml)-&
					ux_c1(iz,grid_x_positive_pml)+vv(iz,ix)**2*&
					ddt**2*(x2(iz+1,ix)+x2(iz-1,ix)-2*&
					x2(iz,ix))/dz**2
              
				x3(iz,ix)=ux_a3(iz,grid_x_positive_pml)+ux_b2(iz,grid_x_positive_pml)+&
					ux_c3(iz,grid_x_positive_pml)

        
			!z_negative
			else if(iz<=nz_bound_u)then
				axis_z_negative_pml=real(iz-(nz_bound_u+1))*dz
				axis_z_thickness_pml=real(nz_bound_u)*dz
				grid_z_negative_pml=iz
              
				d_attenuation_z_present=(3.0*vv(iz,ix)/(2.0*&
					axis_z_thickness_pml))*(axis_z_negative_pml/&
					axis_z_thickness_pml)**2*log(1.0/r)

				deri_d_attenuation_z=(3.0*vv(iz,ix)/(2.0*&
					(axis_z_thickness_pml)))*(2.0*axis_z_negative_pml/&
					axis_z_thickness_pml**2)*log(1.0/r)


				!pay attention to the subsripts of x2 etc.
				uz_a3(iz,ix)=(uz_a2(iz,ix)*(2.0*dz**2+2.0*&
					d_attenuation_z_present*ddt*dz**2-&
					d_attenuation_z_present**2*dz**2*ddt**2)-dz**2*&
					uz_a1(iz,ix)+vv(iz,ix)**2*ddt**2*&
					(x2(iz+1,ix)+&
					x2(iz-1,ix)-2.0*x2(iz,ix)))/&
					(dz**2+2.0*d_attenuation_z_present*ddt*dz**2)

				uz_bp3(iz,ix)=(uz_bp2(iz,ix)*(2.0*dz+2.0*&
					d_attenuation_z_present*ddt*dz-ddt**2*dz*&
					d_attenuation_z_present**2)-dz*uz_bp1(iz,ix)&
					-vv(iz,ix)**2*ddt**2*deri_d_attenuation_z*&
					(x2(iz+1,ix)-x2(iz,ix)))/&
					(dz+2.0*d_attenuation_z_present*ddt*dz)
              
				uz_b2(iz,ix)=uz_b1(iz,ix)*&
					(1-d_attenuation_z_present*ddt)&
					+uz_bp3(iz,ix)*ddt
              
				uz_c3(iz,ix)=2.0*uz_c2(iz,ix)-uz_c1(iz,ix)+&
					vv(iz,ix)**2*ddt**2* &
					(x2(iz,ix+1)+x2(iz,ix-1)-2.0*x2(iz,ix))/(dx**2)
              
				x3(iz,ix)=uz_a3(iz,ix)+uz_b2(iz,ix)+uz_c3(iz,ix)

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
				uz_a3(grid_z_positive_pml,ix)=&
					(uz_a2(grid_z_positive_pml,ix)*(2*dz**2+2*&
					d_attenuation_z_present*ddt*dz**2-&
					d_attenuation_z_present**2*dz**2*ddt**2)-dz**2*&
					uz_a1(grid_z_positive_pml,ix)+vv(iz,ix)**2*&
					ddt**2*(x2(iz+1,ix)+&
					x2(iz-1,ix)-2*x2(iz,ix)))/&
					(dz**2+2*d_attenuation_z_present*ddt*dz**2)

				uz_bp3(grid_z_positive_pml,ix)=&
					(uz_bp2(grid_z_positive_pml,ix)*(2*dz+2*&
					d_attenuation_z_present*ddt*dz-ddt**2*dz*&
					d_attenuation_z_present**2)-dz*&
					uz_bp1(grid_z_positive_pml,ix)&
					-vv(iz,ix)**2*ddt**2*deri_d_attenuation_z*&
					(x2(iz+1,ix)-x2(iz,ix)))/&
					(dz+2*d_attenuation_z_present*ddt*dz)
              
				uz_b2(grid_z_positive_pml,ix)=&
					uz_b1(grid_z_positive_pml,ix)*&
					(1-d_attenuation_z_present*ddt)&
					+uz_bp3(grid_z_positive_pml,ix)*ddt
              
				uz_c3(grid_z_positive_pml,ix)=&
					2*uz_c2(grid_z_positive_pml,ix)-&
					uz_c1(grid_z_positive_pml,ix)+vv(iz,ix)**2*&
					ddt**2*(x2(iz,ix+1)+x2(iz,ix-1)-2*&
					x2(iz,ix))/dx**2
              
				x3(iz,ix)=uz_a3(grid_z_positive_pml,ix)+&
					uz_b2(grid_z_positive_pml,ix)+&
					uz_c3(grid_z_positive_pml,ix)

			endif
		enddo
    enddo

!$omp end do
!$omp end parallel

    return
    end subroutine

!***********************************************************************










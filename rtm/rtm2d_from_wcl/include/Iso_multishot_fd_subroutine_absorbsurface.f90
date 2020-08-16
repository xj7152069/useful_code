!Correct the famous 'negative' to 'positive' mistake
!in pml positive direction of extrapolation
!Add shear wave for better  stablity
!Remove the sigma bug
!Chang sigma's computation method according to Fletcher's way(2009)
!Modified to a VTI version
!Modified for nanjing
!Modified in xibei research institute 
!to a mixed accuracy version for better effciency 2012-9
!*******constant module and externel library*************
        module constant
        implicit none 
        real,parameter::pi=3.1415926
        real,parameter::R=1.0E-3
        end module constant

        module header_module
        type :: segy
         integer(kind = 4) :: tracl
         integer(kind = 4) :: tracr
         integer(kind = 4) :: fldr
         integer(kind = 4) :: tracf
         integer(kind = 4) :: ep
         integer(kind = 4) :: cdp
         integer(kind = 4) :: cdpt
         integer(kind = 2) :: trid
         integer(kind = 2) :: nvs
         integer(kind = 2) :: nhs
         integer(kind = 2) :: duse
         integer(kind = 4) :: offset
         integer(kind = 4) :: gelev
         integer(kind = 4) :: selev
         integer(kind = 4) :: sdepth
         integer(kind = 4) :: gdel
         integer(kind = 4) :: sedl
         integer(kind = 4) :: swdep
         integer(kind = 4) :: gwdep
         integer(kind = 2) :: scalel
         integer(kind = 2) :: scalco
         integer(kind = 4) :: sx
         integer(kind = 4) :: sy
         integer(kind = 4) :: gx
         integer(kind = 4) :: gy
         integer(kind = 2) :: counit
         integer(kind = 2) :: wevel
         integer(kind = 2) :: swevel
         integer(kind = 2) :: sut
         integer(kind = 2) :: gut
         integer(kind = 2) :: sstat
         integer(kind = 2) :: gstat
         integer(kind = 2) :: tstat
         integer(kind = 2) :: laga
         integer(kind = 2) :: lagb
         integer(kind = 2) :: delrt
         integer(kind = 2) :: muts
         integer(kind = 2) :: mute
         integer(kind = 2) :: ns
         integer(kind = 2) :: dt
         integer(kind = 2) :: gain
         integer(kind = 2) :: igc
         integer(kind = 2) :: igi
         integer(kind = 2) :: corr
         integer(kind = 2) :: sfs
         integer(kind = 2) :: sfe
         integer(kind = 2) :: slen
         integer(kind = 2) :: styp
         integer(kind = 2) :: stas
         integer(kind = 2) :: stae
         integer(kind = 2) :: tatyp
         integer(kind = 2) :: afilf
         integer(kind = 2) :: afils
         integer(kind = 2) :: nofilf
         integer(kind = 2) :: nofils
         integer(kind = 2) :: lcf
         integer(kind = 2) :: hcf
         integer(kind = 2) :: lcs
         integer(kind = 2) :: hcs
         integer(kind = 2) :: year
         integer(kind = 2) :: day
         integer(kind = 2) :: hour
         integer(kind = 2) :: minute
         integer(kind = 2) :: sec
         integer(kind = 2) :: timbas
         integer(kind = 2) :: trwf
         integer(kind = 2) :: grnors
         integer(kind = 2) :: grnofr
         integer(kind = 2) :: grnlof
         integer(kind = 2) :: gaps
         integer(kind = 2) :: otrav
         real(kind = 4) :: d1
         real(kind = 4) :: f1
         real(kind = 4) :: d2
         real(kind = 4) :: f2
         real(kind = 4) :: ungpow
         real(kind = 4) :: unscale
         integer(kind = 4) :: ntr
         integer(kind = 2) :: mark
         integer(kind = 2) :: shortpad
         integer(kind = 2) :: unass(14)
         end type
        end module

!********subroutine*****************************************
        subroutine read_data(vz_fn,vz,nx_v,nz_v,myid)
        use constant
        implicit none
        !Dummy variables
        character(len=256)::vz_fn
        real::vz(nx_v,nz_v)
        integer::nx_v,nz_v,nx,nz
        integer::myid
        !Local variables
        integer::i,j
        !Read data from files
        open(10,file=vz_fn,action='read',&
        form='unformatted',access='direct',status='old',recl=nz_v)
        do i=1,nx_v
          read(10,rec=i)vz(i,:)
        enddo
        close (10)
        

!        open(10,file='vz_test',action='write',&
!        form='unformatted',access='direct',status='replace',recl=nz_v)
!        do i=1,nx_v
!          write(10,rec=i)vz(i,:)
!        enddo
!        close (10)

        return
        end subroutine read_data




!********subroutine*****************************************
        subroutine write_currshot_disk(record,shot_fn1,record_acc,shot_fn3,&
								currshot_name,currshot_no,&
                                 nsx,nsz0,currshot_range,currshot_xmin,&
                                 currshot_xmax,nt,dx_v,dz_v,dt)

        use header_module               
        implicit none
        !Dummy variables
        integer::currshot_no,nsx,nsz0,&
                 currshot_range,currshot_xmin,&
                 currshot_xmax,nt
        real::record(nt,currshot_range)
        real::record_acc(nt,currshot_range)
        character(len=256)::shot_fn1
        character(len=256)::shot_fn3
        character(len=256)::currshot_name
        real::dx_v,dz_v,dt
        !Local variables
        type(segy)::su_head
        real::offx_tmp,offy_tmp
        integer::i,j,k,irec

        !Initiallization of the su_head
        su_head=segy(0,0,0,0,0,0,0,0,0,0,0,&
                    0,0,0,0,0,0,0,0,0,0,&
                    0,0,0,0,0,0,0,0,0,0,&
                    0,0,0,0,0,0,0,0,0,0,&
                    0,0,0,0,0,0,0,0,0,0,&
                    0,0,0,0,0,0,0,0,0,0,&
                    0,0,0,0,0,0,0,0,0,0,&
                    0,0,0,0,0,0,0,0,0,0)

        su_head%fldr=currshot_no!shot_no
!        su_head%tracl=currshot_range!number of traces per shot in inline direction
        su_head%trid=1!data_type
        su_head%sx=(nsx-1)*dx_v
        su_head%ns=nt
        su_head%dt=dt*1.0E6!us SEG standard
!        su_head%d1=dx_v!Self_defined
!        su_head%f1=0.0
!        su_head%d2=dz_v


		!Pressure
        open(unit=8,file=trim(adjustl(shot_fn1))//trim(adjustl(currshot_name))//'.su',&
        form='unformatted',access='direct',status='replace',recl=(nt+60))


		!Acceleration
        open(unit=9,file=trim(adjustl(shot_fn3))//trim(adjustl(currshot_name))//'.su',&
        form='unformatted',access='direct',status='replace',recl=(nt+60))


        do i=1,currshot_range

            su_head%gx=(currshot_xmin-1+i-1)*dx_v
			su_head%offset=su_head%gx-su_head%sx

            write(8,rec=i)su_head,(record(k,i),k=1,nt)!**zy
            write(9,rec=i)su_head,(record_acc(k,i),k=1,nt)!**zy


        enddo
        close(unit=8)
		close(unit=9)



        return
        end subroutine write_currshot_disk


        subroutine merge_shot_files(nshot,shot_fn1,shotall_fn1,&
									shot_fn2,shotall_fn2,currshot_xmin,currshot_xmax,nt)

         implicit none
          !Dummy variables
         character(len=256)::shot_fn1,shotall_fn1,shot_fn2,shotall_fn2
         integer::nshot,currshot_xmin,currshot_xmax,nt
         !Local variables
         real,allocatable::record(:,:)
         integer::currshot_ran_x,&
                   is,i,j,k,ierr
         character(len=256)::currshot_noc
         character(len=256)::currshot_name
  
        !Initiallization
        currshot_ran_x=(currshot_xmax-currshot_xmin)+1
        allocate(record(nt+60,currshot_ran_x),STAT=ierr)
        record=0.0



		!================merge No.1 shot fuiles=====================
        open(11,file=trim(adjustl(shotall_fn1))//'.su',access='direct',form='unformatted',&
              action='write',status='replace',recl=currshot_ran_x*(nt+60))  

        do is=1,nshot
          record=0.0
          write(currshot_noc,'(I4)')is
          !Need some test here
          write(currshot_name,*)trim(adjustl(shot_fn1))//trim(adjustl(currshot_noc))//'.su'
          write(*,*)trim(currshot_name),' is bening merged'
          
          open(10,file=currshot_name,access='direct',form='unformatted',&
               action='read',recl=currshot_ran_x*(nt+60))

          !Try to read in this way for better efficiency
           read(10,rec=1)((record(k,i),k=1,nt+60),i=1,currshot_ran_x)
          close (10,status='delete')
!          close (10)

          write(11,rec=is)((record(k,i),k=1,nt+60),i=1,currshot_ran_x)
        enddo

        close(11)


		!================merge No.2 shot fuiles=====================

        open(11,file=trim(adjustl(shotall_fn2))//'.su',access='direct',form='unformatted',&
              action='write',status='replace',recl=currshot_ran_x*(nt+60))  

        do is=1,nshot
          record=0.0
          write(currshot_noc,'(I4)')is
          !Need some test here
          write(currshot_name,*)trim(adjustl(shot_fn2))//trim(adjustl(currshot_noc))//'.su'
          write(*,*)trim(currshot_name),' is bening merged'
          
          open(10,file=currshot_name,access='direct',form='unformatted',&
               action='read',recl=currshot_ran_x*(nt+60))

          !Try to read in this way for better efficiency
           read(10,rec=1)((record(k,i),k=1,nt+60),i=1,currshot_ran_x)
          close (10,status='delete')
!          close (10)

          write(11,rec=is)((record(k,i),k=1,nt+60),i=1,currshot_ran_x)
        enddo

        close(11)


        deallocate(record,STAT=ierr)
        return
        end subroutine merge_shot_files



        subroutine get_currshot_parameters(vz,&
                   currshot_vz,nx_v,nz_v,&
                   currshot_range,currshot_range_all,currshot_xmin,currshot_xmax,&
                   nz,nx_bound_l,nx_bound_r,nz_bound_u,nz_bound_d,nsx,nsz0,&
                   currshot_nsx,currshot_nsz)


        use constant
        implicit none
        !Dummy variables
        real,intent(in)::vz(nx_v,nz_v)
        integer,intent(in)::nx_v,nz_v,currshot_range,&
							currshot_range_all,nz,&
							nx_bound_l,nx_bound_r,&
							nz_bound_u,nz_bound_d,&
							nsx,nsz0,currshot_xmin,&
							currshot_xmax
        real,intent(in out)::currshot_vz(-4:currshot_range_all+5,-4:nz+5)
        integer,intent(in out)::currshot_nsx,currshot_nsz
        !Local  variables
        integer::i,j,ka,i_ori

        !Compute currshot source_position
        currshot_nsx=nsx-currshot_xmin+nx_bound_l!**zy
        currshot_nsz=nsz0+nz_bound_u


        !Assigning&Padding the  parameters
        do i=1,currshot_range
          i_ori=i+currshot_xmin
          if(i_ori<1)i_ori=1
          if(i_ori>nx_v)i_ori=nx_v
          currshot_vz(i+nx_bound_l,nz_bound_u+1:nz_v+nz_bound_u)=vz(i_ori,1:nz_v)
        enddo


            !Using a simpler way to pa velocity
            do    j=-4,nz_bound_u
                do i=-4,currshot_range_all+5
                    currshot_vz(i,j)=currshot_vz(i,nz_bound_u+1)
                enddo
            enddo


            do  j=nz_bound_u+nz_v+1,nz+5
                do i=-4,currshot_range_all+5
                    currshot_vz(i,j)=currshot_vz(i,nz_bound_u+nz_v)
                enddo
            enddo




            do    j=-4,nz+5
                do i=-4,nx_bound_l
                    currshot_vz(i,j)=currshot_vz(nx_bound_l+1,j)
                enddo
            enddo



            do    j=-4,nz+5
                do i=nx_bound_l+currshot_range+1,currshot_range_all+5
                    currshot_vz(i,j)=currshot_vz(nx_bound_l+currshot_range,j)
                enddo
            enddo

        return
        end subroutine get_currshot_parameters

!***********************************************************************************

        subroutine read_par(par_fn,vz_fn,shot_fn1,shot_fn2,shot_fn3,shot_fn4,&
							nshot,offset_min,offset_max,sx0,sz0,rz0,dsx,nx_v,nz_v,&
							dx_v,dz_v,nt,dt,f0,bound_x_l,bound_x_r,bound_z_u,bound_z_d,&
                            order_2nd,order_1st)

        implicit none

        !Dummy variables
        !*Parameter_card_name*
        character(len=256)::par_fn
        !*Read from parameter_card*
        character(len=256)::vz_fn
        character(len=256)::shot_fn1
        character(len=256)::shot_fn2
        character(len=256)::shot_fn3
        character(len=256)::shot_fn4
        integer::nx_v,nz_v
        integer::nt
        integer::order_2nd
        integer::order_1st
        integer::nshot
        real::offset_min,offset_max
        real::dx_v,dz_v
!        real::pml_thick!(m)
        real::bound_x_l,bound_x_r,&
                    bound_z_u,bound_z_d!(m)
        real::sx0,sz0,rz0!(m)
        real::dsx!(m)
        real::f0!(Hz)
        real::dt!(s)
        !*variables of MPI interface*
        integer::myid
      
        !Local variables
        integer::i,j,ierr
        character(len=256)::par_jumpper
        

        open(10,file=par_fn,action='read',status='old',form='formatted',&
        access='sequential')
        
        read(10,'(A)')par_jumpper
        read(10,'(A)')vz_fn
        read(10,'(A)')par_jumpper
        read(10,'(A)')shot_fn1
        read(10,'(A)')par_jumpper
        read(10,'(A)')shot_fn2!**zy
        read(10,'(A)')par_jumpper
        read(10,'(A)')shot_fn3!**zy
        read(10,'(A)')par_jumpper
        read(10,'(A)')shot_fn4!**zy
        read(10,'(A)')par_jumpper
        read(10,*)nshot
        read(10,'(A)')par_jumpper
        read(10,*)offset_min,offset_max
        read(10,'(A)')par_jumpper
        read(10,*)sx0,sz0,rz0
        read(10,'(A)')par_jumpper
        read(10,'(F)')dsx
        read(10,'(A)')par_jumpper
        read(10,*)nx_v,nz_v
        read(10,'(A)')par_jumpper
        read(10,*)dx_v,dz_v
        read(10,'(A)')par_jumpper
        read(10,'(I)')nt
        read(10,'(A)')par_jumpper
        read(10,'(F)')dt
        read(10,'(A)')par_jumpper
        read(10,'(F)')f0
        read(10,'(A)')par_jumpper
        read(10,*)bound_x_l,bound_x_r
        read(10,'(A)')par_jumpper
        read(10,*)bound_z_u,bound_z_d!**zy
        read(10,'(A)')par_jumpper
        read(10,'(I)')order_2nd
        read(10,'(A)')par_jumpper
        read(10,'(I)')order_1st
        close(10)

		return
        end subroutine read_par

!**********************************************************************************
  subroutine  extrapolation_back(vz,&
                                       nsx,nsz0,nrz0,nx_v,nz_v,nx,nz,dx,dz,&
                                       nt,dt,f0,nx_bound_l,nx_bound_r,nz_bound_u,nz_bound_d,&
                                       order_2nd,order_1st,&
                                       u1,u2,&
									   coe_2nd_10,coe_2nd_8,coe_2nd_6,coe_2nd_4,coe_2nd_2,coe_1st&
                                       ) 
  !The basic version of extrapolation function
  !and much modification need to be applied to it
  use constant
  implicit none
  !Dummy variables


  real::vz(-4:nx+5,-4:nz+5)
  integer::nsx,nsz0,nrz0
  integer::nx_v,nz_v,nx,nz,nt,nx_bound_l,nx_bound_r,nz_bound_u,nz_bound_d
  integer::order_2nd,order_1st
  real::dx,dz,dt,f0
  real::u1(-4:nx+5,-4:nz+5),u2(-4:nx+5,-4:nz+5)
  real::coe_2nd_10(order_2nd/2)
  real::coe_2nd_8(4)
  real::coe_2nd_6(3)
  real::coe_2nd_4(2)
  real::coe_2nd_2(1)
  real::coe_1st((order_1st+1)/2)

  real::sum_x_u_d2,sum_z_u_d2


  !Local variables
  integer::i,j,k,ierr


!$OMP PARALLEL PRIVATE(sum_x_u_d2,&
!$OMP	sum_z_u_d2)
!$OMP DO 

    do j=nz_bound_u+1,nz-nz_bound_d
      do i=nx_bound_l+1,nx-nx_bound_r

			sum_x_u_d2=(coe_2nd_10(1)*&
			(u2(i+1,j)+u2(i-1,j)-2*u2(i,j))+&
			coe_2nd_10(2)*&
			(u2(i+2,j)+u2(i-2,j)-2*u2(i,j))+&
			coe_2nd_10(3)*&
			(u2(i+3,j)+u2(i-3,j)-2*u2(i,j))+&
			coe_2nd_10(4)*&
			(u2(i+4,j)+u2(i-4,j)-2*u2(i,j))+&
			coe_2nd_10(5)*&
			(u2(i+5,j)+u2(i-5,j)-2*u2(i,j)))/dx**2


               sum_z_u_d2=(coe_2nd_10(1)*&
             (u2(i,j+1)+u2(i,j-1)-2*u2(i,j))+&
             coe_2nd_10(2)*&
             (u2(i,j+2)+u2(i,j-2)-2*u2(i,j))+&
             coe_2nd_10(3)*&
             (u2(i,j+3)+u2(i,j-3)-2*u2(i,j))+&
             coe_2nd_10(4)*&
             (u2(i,j+4)+u2(i,j-4)-2*u2(i,j))+&
             coe_2nd_10(5)*&
             (u2(i,j+5)+u2(i,j-5)-2*u2(i,j)))/dz**2


			u1(i,j)=2*u2(i,j)-u1(i,j)+dt**2*vz(i,j)**2*(sum_x_u_d2+sum_z_u_d2)


    	enddo
  	enddo

!$OMP END DO
!$OMP END PARALLEL
  end subroutine extrapolation_back





!**********************************************************************************
  subroutine  extrapolation_one_step(vz,&
                                       nsx,nsz0,nrz0,nx_v,nz_v,nx,nz,dx,dz,&
                                       nt,dt,f0,nx_bound_l,nx_bound_r,nz_bound_u,nz_bound_d,&
                                       order_2nd,order_1st,&
                                       u1,u2,&
									   coe_2nd_10,coe_2nd_8,coe_2nd_6,coe_2nd_4,coe_2nd_2,coe_1st,&
                                       u1_1_x,u1_2_x,u2_1_x,u2_2_x,&
                                       u2_tmp_1_x,u2_tmp_2_x,u3_1_x,u3_2_x,&
                                       u1_1_z,u1_2_z,u2_1_z,u2_2_z,&
                                       u2_tmp_1_z,u2_tmp_2_z,u3_1_z,u3_2_z)
                                      
  !The basic version of extrapolation function
  !and much modification need to be applied to it
  use constant
  implicit none
  !Dummy variables


  real::vz(-4:nx+5,-4:nz+5)
  integer::nsx,nsz0,nrz0
  integer::nx_v,nz_v,nx,nz,nt,nx_bound_l,nx_bound_r,nz_bound_u,nz_bound_d
  integer::order_2nd,order_1st
  real::dx,dz,dt,f0
  real::u1(-4:nx+5,-4:nz+5),u2(-4:nx+5,-4:nz+5)
  real::coe_2nd_10(order_2nd/2)
  real::coe_2nd_8(4)
  real::coe_2nd_6(3)
  real::coe_2nd_4(2)
  real::coe_2nd_2(1)
  real::coe_1st((order_1st+1)/2)


  real::u1_1_x(1:nx_bound_l+nx_bound_r,1:nz)
  real::u1_2_x(1:nx_bound_l+nx_bound_r,1:nz)
  real::u2_1_x(1:nx_bound_l+nx_bound_r,1:nz)
  real::u2_2_x(1:nx_bound_l+nx_bound_r,1:nz)
  real::u2_tmp_1_x(1:nx_bound_l+nx_bound_r,1:nz)
  real::u2_tmp_2_x(1:nx_bound_l+nx_bound_r,1:nz)
  real::u3_1_x(1:nx_bound_l+nx_bound_r,1:nz)
  real::u3_2_x(1:nx_bound_l+nx_bound_r,1:nz)

  real::u1_1_z(1:nx,1:nz_bound_u+nz_bound_d)
  real::u1_2_z(1:nx,1:nz_bound_u+nz_bound_d)
  real::u2_1_z(1:nx,1:nz_bound_u+nz_bound_d)
  real::u2_2_z(1:nx,1:nz_bound_u+nz_bound_d)
  real::u2_tmp_1_z(1:nx,1:nz_bound_u+nz_bound_d)
  real::u2_tmp_2_z(1:nx,1:nz_bound_u+nz_bound_d)
  real::u3_1_z(1:nx,1:nz_bound_u+nz_bound_d)
  real::u3_2_z(1:nx,1:nz_bound_u+nz_bound_d)


  !Local variables
  integer::i,j,k,ierr
  !*Local variables for PML wavefield*
  integer::grid_x_positive_pml,grid_x_negative_pml,&
  grid_z_positive_pml,grid_z_negative_pml
  real::axis_x_positive_pml,axis_x_thickness_pml,&
  axis_x_negative_pml,axis_z_positive_pml,&
  axis_z_thickness_pml,axis_z_negative_pml
  real::d_attenuation_x_present,deri_d_attenuation_x,&
  d_attenuation_z_present,deri_d_attenuation_z
  real::sum_x_u_d1,sum_x_u_d2,sum_z_u_d2,sum_z_u_d1


!$OMP PARALLEL PRIVATE(sum_x_u_d2,sum_z_u_d2,sum_x_u_d1,sum_z_u_d1,&
!$OMP d_attenuation_x_present,deri_d_attenuation_x,&
!$OMP d_attenuation_z_present,deri_d_attenuation_z,&
!$OMP axis_x_negative_pml,axis_x_thickness_pml,grid_x_negative_pml,&
!$OMP axis_x_positive_pml,grid_x_positive_pml,&
!$OMP axis_z_negative_pml,axis_z_thickness_pml,grid_z_negative_pml,&
!$OMP axis_z_positive_pml,grid_z_positive_pml) 
!$OMP DO 

    do j=1,nz
      do i=1,nx

			sum_x_u_d2=(coe_2nd_10(1)*&
			(u2(i+1,j)+u2(i-1,j)-2*u2(i,j))+&
			coe_2nd_10(2)*&
			(u2(i+2,j)+u2(i-2,j)-2*u2(i,j))+&
			coe_2nd_10(3)*&
			(u2(i+3,j)+u2(i-3,j)-2*u2(i,j))+&
			coe_2nd_10(4)*&
			(u2(i+4,j)+u2(i-4,j)-2*u2(i,j))+&
			coe_2nd_10(5)*&
			(u2(i+5,j)+u2(i-5,j)-2*u2(i,j)))/dx**2

              

               sum_z_u_d2=(coe_2nd_10(1)*&
             (u2(i,j+1)+u2(i,j-1)-2*u2(i,j))+&
             coe_2nd_10(2)*&
             (u2(i,j+2)+u2(i,j-2)-2*u2(i,j))+&
             coe_2nd_10(3)*&
             (u2(i,j+3)+u2(i,j-3)-2*u2(i,j))+&
             coe_2nd_10(4)*&
             (u2(i,j+4)+u2(i,j-4)-2*u2(i,j))+&
             coe_2nd_10(5)*&
             (u2(i,j+5)+u2(i,j-5)-2*u2(i,j)))/dz**2




			if(i>=nx_bound_l+1.and.i<=nx-nx_bound_r.and.&
           		j>=nz_bound_u+1.and.j<=nz-nz_bound_d)then



			u1(i,j)=2*u2(i,j)-u1(i,j)+dt**2*vz(i,j)**2*(sum_x_u_d2+sum_z_u_d2)


	      !PML wavefield extrapolation
	        else if(i<=nx_bound_l)then
	
	             sum_x_u_d1=(coe_1st(1)*&
	           (u2(i+1,j)-u2(i-1,j))+&
	           coe_1st(2)*&
	           (u2(i+2,j)-u2(i-2,j))+&
	           coe_1st(3)*&
	           (u2(i+3,j)-u2(i-3,j))+&
	           coe_1st(4)*&
	           (u2(i+4,j)-u2(i-4,j))+&
	           coe_1st(5)*&
	           (u2(i+5,j)-u2(i-5,j)))/dx
	
	    
	              
				axis_x_negative_pml=real(i-(nx_bound_l+1))*dx
		        axis_x_thickness_pml=real(nx_bound_l)*dx
		        grid_x_negative_pml=i
		
		
		        !For test
		        d_attenuation_x_present=(3.0*vz(i,j)/(2.0*axis_x_thickness_pml))&
		        *(axis_x_negative_pml/axis_x_thickness_pml)**2*log(1/R)
		
		        deri_d_attenuation_x=(3.0*vz(i,j)/(2.0*(axis_x_thickness_pml)))*&
		        (2.0*axis_x_negative_pml/axis_x_thickness_pml**2)*log(1/R)
		
		        u1_1_x(i,j)=(u1_2_x(i,j)*(2*dx**2+2*d_attenuation_x_present*dt*dx**2-&
		        d_attenuation_x_present**2*dx**2*dt**2)-&
		        dx**2*u1_1_x(i,j)+vz(i,j)**2*dt**2*(sum_x_u_d2*dx**2))/&
		        (dx**2+2*d_attenuation_x_present*dt*dx**2)
		
		
		        u2_tmp_1_x(i,j)=(u2_tmp_2_x(i,j)*(2*dx+2*d_attenuation_x_present*dt*dx-&
		        dt**2*dx*d_attenuation_x_present**2)-&
		        dx*u2_tmp_1_x(i,j)-vz(i,j)**2*dt**2*deri_d_attenuation_x*(sum_x_u_d1*dx))/&
		        (dx+2*d_attenuation_x_present*dt*dx)
		
		        u2_1_x(i,j)=u2_2_x(i,j)*(1-d_attenuation_x_present*dt)+u2_tmp_1_x(i,j)*dt
		
		
		        u3_1_x(i,j)=2*u3_2_x(i,j)-u3_1_x(i,j)+vz(i,j)**2*dt**2*(sum_z_u_d2*dz**2)/dz**2
		
		        u1(i,j)=u1_1_x(i,j)+u2_1_x(i,j)+u3_1_x(i,j)
	        
	
	
	        else if(i>=nx-nx_bound_r+1)then
	
	
	             sum_x_u_d1=(coe_1st(1)*&
	           (u2(i+1,j)-u2(i-1,j))+&
	           coe_1st(2)*&
	           (u2(i+2,j)-u2(i-2,j))+&
	           coe_1st(3)*&
	           (u2(i+3,j)-u2(i-3,j))+&
	           coe_1st(4)*&
	           (u2(i+4,j)-u2(i-4,j))+&
	           coe_1st(5)*&
	           (u2(i+5,j)-u2(i-5,j)))/dx
	
	    
	
	
				axis_x_positive_pml=real(i-(nx-nx_bound_r))*dx
		        axis_x_thickness_pml=real(nx_bound_r)*dx
		        grid_x_positive_pml=i-(nx-nx_bound_r)+nx_bound_l
		
		
		        d_attenuation_x_present=(3.0*vz(i,j)/(2.0*axis_x_thickness_pml))*&
		        (axis_x_positive_pml/axis_x_thickness_pml)**2*log(1/R)
		        
		        deri_d_attenuation_x=(3.0*vz(i,j)/(2.0*axis_x_thickness_pml))*&
		        (2.0*axis_x_positive_pml/axis_x_thickness_pml**2)*log(1/R)
		              
		        u1_1_x(grid_x_positive_pml,j)=(u1_2_x(grid_x_positive_pml,j)*&
		        (2*dx**2+2*d_attenuation_x_present*dt*dx**2-&
		        d_attenuation_x_present**2*dx**2*dt**2)-&
		        dx**2*u1_1_x(grid_x_positive_pml,j)+vz(i,j)**2*dt**2*(sum_x_u_d2*dx**2))/&
		        (dx**2+2*d_attenuation_x_present*dt*dx**2)
		
		        u2_tmp_1_x(grid_x_positive_pml,j)=&
				(u2_tmp_2_x(grid_x_positive_pml,j)*&
		        (2*dx+2*d_attenuation_x_present*dt*dx-&
		        dt**2*dx*d_attenuation_x_present**2)-&
		        dx*u2_tmp_1_x(grid_x_positive_pml,j)-&
				vz(i,j)**2*dt**2*deri_d_attenuation_x*(sum_x_u_d1*dx))/&
		        (dx+2*d_attenuation_x_present*dt*dx)
		
		        u2_1_x(grid_x_positive_pml,j)=u2_2_x(grid_x_positive_pml,j)*&
		        (1-d_attenuation_x_present*dt)+u2_tmp_1_x(grid_x_positive_pml,j)*dt
		
		        !Pay attenuation to the subscripts  of vz
		        u3_1_x(grid_x_positive_pml,j)=2*u3_2_x(grid_x_positive_pml,j)-&
		        u3_1_x(grid_x_positive_pml,j)+vz(i,j)**2*dt**2*(sum_z_u_d2*dz**2)/dz**2
		
		        u1(i,j)=u1_1_x(grid_x_positive_pml,j)+&
		        u2_1_x(grid_x_positive_pml,j)+u3_1_x(grid_x_positive_pml,j)
	        
	
	
	        else if(j<=nz_bound_u)then

	             sum_z_u_d1=(coe_1st(1)*&
	           (u2(i,j+1)-u2(i,j-1))+&
	           coe_1st(2)*&
	           (u2(i,j+2)-u2(i,j-2))+&
	           coe_1st(3)*&
	           (u2(i,j+3)-u2(i,j-3))+&
	           coe_1st(4)*&
	           (u2(i,j+4)-u2(i,j-4))+&
	           coe_1st(5)*&
	           (u2(i,j+5)-u2(i,j-5)))/dx
	
	
				axis_z_negative_pml=real(j-(nz_bound_u+1))*dz
	         	axis_z_thickness_pml=real(nz_bound_u)*dz
				grid_z_negative_pml=j
	
	
		        d_attenuation_z_present=(3.0*vz(i,j)/(2.0*axis_z_thickness_pml))&
		        *(axis_z_negative_pml/axis_z_thickness_pml)**2*log(1/R)
		
		        deri_d_attenuation_z=(3.0*vz(i,j)/(2.0*(axis_z_thickness_pml)))*&
		        (2.0*axis_z_negative_pml/axis_z_thickness_pml**2)*log(1/R)
		
		
		        u1_1_z(i,j)=(u1_2_z(i,j)*(2*dz**2+2*d_attenuation_z_present*dt*dz**2-&
		        d_attenuation_z_present**2*dz**2*dt**2)-&
		        dz**2*u1_1_z(i,j)+vz(i,j)**2*dt**2*(sum_z_u_d2*dz**2))/&
		        (dz**2+2*d_attenuation_z_present*dt*dz**2)
		
		
		        u2_tmp_1_z(i,j)=(u2_tmp_2_z(i,j)*(2*dz+2*d_attenuation_z_present*dt*dz-&
		        dt**2*dz*d_attenuation_z_present**2)-&
		        dz*u2_tmp_1_z(i,j)-vz(i,j)**2*dt**2*deri_d_attenuation_z*(sum_z_u_d1*dz))/&
		        (dz+2*d_attenuation_z_present*dt*dz)
		
		        u2_1_z(i,j)=u2_2_z(i,j)*(1-d_attenuation_z_present*dt)+u2_tmp_1_z(i,j)*dt
		
		        u3_1_z(i,j)=2*u3_2_z(i,j)-u3_1_z(i,j)+vz(i,j)**2*dt**2*sum_x_u_d2
		
		        u1(i,j)=u1_1_z(i,j)+u2_1_z(i,j)+u3_1_z(i,j)






	
	        else if(j>=nz-nz_bound_d+1)then
	
	        
	             sum_z_u_d1=(coe_1st(1)*&
	           (u2(i,j+1)-u2(i,j-1))+&
	           coe_1st(2)*&
	           (u2(i,j+2)-u2(i,j-2))+&
	           coe_1st(3)*&
	           (u2(i,j+3)-u2(i,j-3))+&
	           coe_1st(4)*&
	           (u2(i,j+4)-u2(i,j-4))+&
	           coe_1st(5)*&
	           (u2(i,j+5)-u2(i,j-5)))/dz
	
	         
				axis_z_positive_pml=real(j-(nz-nz_bound_d))*dz
	         	axis_z_thickness_pml=real(nz_bound_d)*dz
				grid_z_positive_pml=j-(nz-nz_bound_d)+nz_bound_u
	
	
	        	d_attenuation_z_present=(3.0*vz(i,j)/(2.0*axis_z_thickness_pml))*&
	        	(axis_z_positive_pml/axis_z_thickness_pml)**2*log(1/R)
	
	        	deri_d_attenuation_z=(3.0*vz(i,j)/(2.0*(axis_z_thickness_pml)))*&
	        	(2.0*axis_z_positive_pml/axis_z_thickness_pml**2)*log(1/R)
	
	        	u1_1_z(i,grid_z_positive_pml)=(u1_2_z(i,grid_z_positive_pml)*&
		        (2*dz**2+2*d_attenuation_z_present*dt*dz**2-&
		        d_attenuation_z_present**2*dz**2*dt**2)-&
		        dz**2*u1_1_z(i,grid_z_positive_pml)+vz(i,j)**2*dt**2*(sum_z_u_d2*dz**2))/&
		        (dz**2+2*d_attenuation_z_present*dt*dz**2)
		
		
		        u2_tmp_1_z(i,grid_z_positive_pml)=(u2_tmp_2_z(i,grid_z_positive_pml)*&
		        (2*dz+2*d_attenuation_z_present*dt*dz-&
		        dt**2*dz*d_attenuation_z_present**2)-dz*u2_tmp_1_z(i,grid_z_positive_pml)&
		        -vz(i,j)**2*dt**2*deri_d_attenuation_z*(sum_z_u_d1*dz))/&
		        (dz+2*d_attenuation_z_present*dt*dz)
		
		        u2_1_z(i,grid_z_positive_pml)=u2_2_z(i,grid_z_positive_pml)*&
		        (1-d_attenuation_z_present*dt)+u2_tmp_1_z(i,grid_z_positive_pml)*dt
		
		        u3_1_z(i,grid_z_positive_pml)=2*u3_2_z(i,grid_z_positive_pml)-&
		        u3_1_z(i,grid_z_positive_pml)+vz(i,j)**2*dt**2*sum_x_u_d2
		
		        u1(i,j)=u1_1_z(i,grid_z_positive_pml)+u2_1_z(i,grid_z_positive_pml)+&
		        u3_1_z(i,grid_z_positive_pml)


      		endif

    	enddo
  	enddo

!$OMP END DO
!$OMP END PARALLEL
  end subroutine extrapolation_one_step





!========================================================================
	subroutine add_source(u1,currshot_nsx,currshot_nsz,nx,nz,&
						f0,it,dt, currshot_no, nshot, shot_ene)

	use constant
	implicit none
	!Dummy variables
	real::u1(-4:nx+5,-4:nz+5)

	real::f0,dt
	integer::currshot_no
	integer::currshot_nsx,currshot_nsz,nx,nz,it
	integer::nshot
	real	::shot_ene(nshot)
	!Local variables
	integer::i,j,k,hw,idx1,idx2,m,n,ierr
	integer,parameter::n1=20,n2=20
	real,parameter::alpha=0.5
	real::wavelet,atten
	real:: coefficent
	real	::	random
	


	hw=ceiling((n1+1)/2.0)
	wavelet=exp(-(pi*f0*(it*dt-1.0/f0))**2)*(1-2*(pi*f0*(it*dt-1.0/f0))**2)
!        wavelet=sin(2*pi*f0*(it*dt-1.0/f0))*EXP(-((2*pi*f0/6)**2)*(it*dt-1.0/f0)**2)

!	if(currshot_no .le. 10)then
!		coefficent=currshot_no*0.5
!	endif
!	if(currshot_no .gt. 10 .and. currshot_no .le.20)then
!		coefficent=(currshot_no-10)*0.5
!	endif
!	if(currshot_no .gt. 20)then
!		coefficent=1.0
!	endif

	coefficent=shot_ene(currshot_no)

!	write(*,*)'currshotnsx=,currshot_nsz, currshot_no ,coefficent',
	write(*,*) 'currshot_no=', currshot_no ,coefficent

!	u1(currshot_nsx,currshot_nsz)=wavelet+u1(currshot_nsx,currshot_nsz)
	u1(currshot_nsx,currshot_nsz)=coefficent*wavelet+u1(currshot_nsx,currshot_nsz)


!	do m=currshot_nsz-hw,currshot_nsz+hw
!		do n=currshot_nsx-hw,currshot_nsx+hw
!              	idx1=m-currshot_nsz
!	      		idx2=n-currshot_nsx
!				atten=exp(-1.0*alpha*(idx1**2+idx2**2))
!
!				if(m>=-4.and.m<=nz+5.and.n>=-4.and.n<=nx+5)then
!					u1(n,m)=wavelet*atten+u1(n,m)
!				endif
!
!		enddo
!    enddo




	return
	end subroutine add_source




!=========================================================================
	subroutine get_acc_record(u,record_acc,order_1st,coe_1st,&
								nx_v,nx,nz_v,nz,nx_bound_l,nz_bound_u,nrz0,nt,it,dz)

	implicit none
	!Dummy variables
	integer::currshot_range,currshot_range_all,&
			nx_v,nx,nz_v,nz,nx_bound_l,nz_bound_u,nrz0,nt,it,&
			order_1st
	real::dz
	real::u(-4:nx+5,-4:nz+5),record_acc(nt,nx_v),&
			coe_1st((order_1st+1)/2)
	!Local variables
	integer::i,j,ixu,izu
	real::sum_z_u_d1

	izu=nrz0+nz_bound_u


	do i=1,nx_v

             sum_z_u_d1=(coe_1st(1)*&
           (u(i+nx_bound_l,izu+1)-u(i+nx_bound_l,izu-1))+&
           coe_1st(2)*&
           (u(i+nx_bound_l,izu+2)-u(i+nx_bound_l,izu-2))+&
           coe_1st(3)*&
           (u(i+nx_bound_l,izu+3)-u(i+nx_bound_l,izu-3))+&
           coe_1st(4)*&
           (u(i+nx_bound_l,izu+4)-u(i+nx_bound_l,izu-4))+&
           coe_1st(5)*&
           (u(i+nx_bound_l,izu+5)-u(i+nx_bound_l,izu-5)))/dz


			record_acc(it,i)=sum_z_u_d1

	enddo



	return

	end subroutine get_acc_record


!=========================================================================

    subroutine coefficient_2nd(order_even,coefficients)
    !use IMSL
    implicit none
     !Dummy variables
     integer::order_even
     real::coefficients(order_even/2)
    !routine variables
    real,allocatable::a(:,:),b(:)
    real::i,j
    real::fact
     integer::k1
     k1=order_even/2
    allocate(a(k1,k1))
    allocate(b(k1))
     !Initiallization
     a=0.0
     b=0.0
    !Computation
    do i=1,k1
      do j=1,k1
        fact=1.0
        a(i,j)=(i**(2*j))/factorical((2*j),fact)
      enddo
    enddo
    !The operation of calling can not be omitted
    call inv()
    deallocate(a)
    deallocate(b)
    return
     contains
       real function factorical(l,fact)
       implicit none
       real::fact
       real::m,l
       do m=1,l
         fact=fact*m
       enddo
      factorical=fact
      !write(*,*)fact
      return
      end function factorical
           
      subroutine inv()
      implicit none
      integer::i1,j1
      real::ftmp1,ftmp5
      real,allocatable::inv_a(:,:),temp(:),ftmp2(:,:),ftmp3(:,:),ftmp4(:),p1(:)
      allocate(p1(k1))
      allocate(temp(k1))
      allocate(ftmp2(k1,k1))
      allocate(ftmp3(k1,k1))
      allocate(inv_a(k1,k1))
      allocate(ftmp4(k1))
      !Initiallization
      p1=0.0
      temp=0.0
      ftmp2=0.0
      ftmp3=0.0
      ftmp4=0.0
      inv_a=0.0
      !Computation
      do i=1,k1
         do j=1,k1
           if (i==j)then
             inv_a(i,j)=1
           endif
         enddo
      enddo
      do i=1,k1
         temp=a(i,:)
        temp(i)=temp(i)-1
        ftmp1=sum(temp*inv_a(:,i))+1
        ftmp4=inv_a(:,i)
        do j1=1,k1
          do i1=1,k1
              ftmp2(j1,i1)=ftmp4(j1)*temp(i1)
            enddo
        enddo
         ftmp3=matmul(ftmp2,inv_a)
         do j=1,k1
           do j1=1,k1
            ftmp3(j,j1)=ftmp3(j,j1)/ftmp1
          enddo
        enddo
        inv_a=inv_a-ftmp3
      enddo
      do i=1,k1
        b=0
        b(i)=0.5
        p1=matmul(inv_a,b)
        coefficients(i)=p1(1)
      enddo
      deallocate(p1)
       deallocate(temp)
       deallocate(ftmp2)
       deallocate(ftmp3)
       deallocate(inv_a)
       deallocate(ftmp4)
       end subroutine inv
    end subroutine  coefficient_2nd


    subroutine coefficient_1st(order_odd,coefficients)
    !use IMSL
    implicit none
    !Dummy variables
    integer::order_odd
    real::coefficients((order_odd+1)/2)
    !routine variables
    real,allocatable::a(:,:),b(:)
    real::i,j
    real::fact
    integer::k1_1
    k1_1=((order_odd+1)/2)!Pay attenuation and need more thicking
    allocate(a(k1_1,k1_1))
    allocate(b(k1_1))
    !Initiallization
    a=0.0
    b=0.0
    !Computation
    do i=1,k1_1
      do j=1,k1_1
        fact=1.0
        a(i,j)=(i**(2*j-1))/factorical((2*j-1),fact)
      enddo
    enddo
    !The operation of calling can not be omitted
     call inv_1st()
     
     deallocate(a)
     deallocate(b)
     return
     contains
       real function factorical(l,fact)
       implicit none
       real::fact
       real::m,l
       do m=1,l
         fact=fact*m
       enddo
      factorical=fact
      return
      end function factorical
           
      subroutine inv_1st()
      implicit none
      integer::i1,j1
      real::ftmp1,ftmp5
      real,allocatable::inv_a(:,:),temp(:),ftmp2(:,:),ftmp3(:,:),ftmp4(:),p1(:)
      allocate(p1(k1_1))
      allocate(temp(k1_1))
      allocate(ftmp2(k1_1,k1_1))
      allocate(ftmp3(k1_1,k1_1))
      allocate(inv_a(k1_1,k1_1))
      allocate(ftmp4(k1_1))
      !Initiallization
      p1=0.0
      temp=0.0
      ftmp2=0.0
      ftmp3=0.0
      inv_a=0.0
      ftmp4=0.0
      do i=1,k1_1
         do j=1,k1_1
           if (i==j)then
             inv_a(i,j)=1
           endif
         enddo
      enddo
     do i=1,k1_1
       temp=a(i,:)
      temp(i)=temp(i)-1
      ftmp1=sum(temp*inv_a(:,i))+1
      ftmp4=inv_a(:,i)
      do j1=1,k1_1
        do i1=1,k1_1
          ftmp2(j1,i1)=ftmp4(j1)*temp(i1)
        enddo
      enddo
       ftmp3=matmul(ftmp2,inv_a)
       do j=1,k1_1
         do j1=1,k1_1
          ftmp3(j,j1)=ftmp3(j,j1)/ftmp1
        enddo
      enddo
      inv_a=inv_a-ftmp3
     enddo
    do i=1,k1_1
      b=0
      b(i)=0.5
      p1=matmul(inv_a,b)
      coefficients(i)=p1(1)!Pay attenuation
    enddo
    !=======Print some basic information on the screen====
!    write(*,*)'=====The order of time derivatives is======'
!    write(*,*)'          4'
!    write(*,*)'=====The order of 1st spatial derivatives is==='
!    write(*,*)k_1
!    write(*,*)'=====The difference coefficients for 1st spatial derivatives are========'
!    write(*,*)p_1

    deallocate(p1)
    deallocate(temp)
    deallocate(ftmp2)
    deallocate(ftmp3)
    deallocate(inv_a)
    deallocate(ftmp4)
     end subroutine inv_1st
    end subroutine  coefficient_1st

        

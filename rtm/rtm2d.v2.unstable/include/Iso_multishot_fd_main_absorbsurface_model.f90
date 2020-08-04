        !Modified to conduct Apture_control during forward modelling
        !Modified to add su head to shot_gather
        !Modified for nanjing
        !Modified to add acc version
        !Modified to a new source implementation 
        !Modified to a acoustic for multiple simlution
        include 'Iso_multishot_fd_subroutine_absorbsurface.f90'


        program Iso_fd_qacoustic_modelling
        implicit none 
        character(len=256)::par_fn=&
!		'./iso_2d_modelling_v1.par'
!		'./iso_2d_modelling_waxian.par'
!		'./iso_2d_modelling_two_layer.par'
!		'./iso_2d_modelling_one_layer_v2.par'
		'./iso_2d_modelling_two_layer_zwy.par'
!        './iso_2d_modelling_sigsbee.par'
!        './iso_2d_modelling_homocline.par'
        call multi_shot_modelling(par_fn)
        end program Iso_fd_qacoustic_modelling

        subroutine multi_shot_modelling(par_fn)
        use mpi
        use constant
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
        integer::order_1st
        integer::order_2nd
        integer::nshot
        real::offset_min,offset_max
        real::dx_v,dz_v
        real::x_bound_l,x_bound_r
        real::z_bound_u,z_bound_d!**zy
        real::sx0,sz0,rz0!(m)
        real::dsx!(m)
        real::f0!(Hz)
        real::dt!(s)
        !*Dummy variables(Out)*
        integer::nx_bound_l,nx_bound_r,nz_bound_u,nz_bound_d
        integer::nsx0,nsz0,nrz0
        integer::currshot_xmin,currshot_xmax
        integer::ndsx,nsx
        integer::currshot_no
        character(len=256)::currshot_name
        real,allocatable::vz(:,:)
        !Local variables
        !*Variables for MPI interface*
        integer::myid,nproc,ierr
        integer::status(MPI_STATUS_SIZE)  
        !*Other local variables*
        integer::i,j,k,err

        !Initiallization of MPI interface
        call MPI_init(ierr)
        call MPI_comm_rank(MPI_COMM_WORLD,myid,ierr)
        call MPI_comm_size(MPI_COMM_WORLD,nproc,ierr)
        call MPI_barrier(MPI_COMM_WORLD,ierr)

        !Read parameter_card
        call read_par(par_fn,vz_fn,shot_fn1,shot_fn2,shot_fn3,shot_fn4,&
                      nshot,offset_min,offset_max,sx0,sz0,rz0,dsx,nx_v,nz_v,&
                      dx_v,dz_v,nt,dt,f0,x_bound_l,x_bound_r,z_bound_u,z_bound_d,&
                      order_2nd,order_1st)


        allocate(vz(nx_v,nz_v),STAT=err)



        !Initiallization
        vz=0.0
        !*Transfrom coordinates to grid_num for the conveniency of computation*!

        nx_bound_l=nint(x_bound_l/dx_v)
        nx_bound_r=nint(x_bound_r/dx_v)

        nz_bound_u=nint(z_bound_u/dz_v)
        nz_bound_d=nint(z_bound_d/dz_v)


!       write(*,*)'nx_bound_l,nx_bound_r,nz_bound_u,nz_bound_d'
!       write(*,*)nx_bound_l,nx_bound_r,nz_bound_u,nz_bound_d

        nsx0=nint(sx0/dx_v)+1!**zy                                    
        nsz0=nint(sz0/dz_v)+1!**zy                                    
        nrz0=nint(rz0/dz_v)+1!**zy
        !*Read parameters from files*!
        call read_data(vz_fn,vz,nx_v,nz_v,myid)


        !Multishot modelling
        if (myid==0)write(*,*)'==========2D_Iso_qusi_acoustic_modelling &
        began=========='



        do i=1+myid,nshot,nproc
          nsx=nsx0+nint((i-1)*dsx/dx_v)
          currshot_xmin=nsx+nint(offset_min/dx_v)
          currshot_xmax=nsx+nint(offset_max/dx_v)
          currshot_no=i


          write(currshot_name,'(I4)')i

          call single_shot_modelling(vz,&
               currshot_xmin,currshot_xmax,nsx,nsz0,nrz0,nx_v,nz_v,dx_v,dz_v,&
               nt,dt,f0,nx_bound_l,nx_bound_r,nz_bound_u,nz_bound_d,&
               order_2nd,order_1st,&
               shot_fn1,shot_fn2,shot_fn3,shot_fn4,currshot_name,currshot_no,myid)
        enddo

        call MPI_barrier(MPI_COMM_WORLD,ierr)
        if (myid==0)then
          write(*,*)'Now merging shot_record files begins'
          call merge_shot_files(nshot,shot_fn1,shot_fn2,&
              shot_fn3,shot_fn4,currshot_xmin,currshot_xmax,nt)
        write(*,*)'==========2D_Iso_acoustic_modelling end&
        =========='
        endif
        
        call MPI_finalize(ierr)

        deallocate(vz,STAT=err)
        
        return
        end subroutine multi_shot_modelling


        subroutine single_shot_modelling(vz,&
                    currshot_xmin,currshot_xmax,nsx,nsz0,nrz0,nx_v,nz_v,dx,dz,&
                    nt,dt,f0,nx_bound_l,nx_bound_r,nz_bound_u,nz_bound_d,&
                    order_2nd,order_1st,&
                    shot_fn1,shot_fn2,shot_fn3,shot_fn4,currshot_name,currshot_no,myid)
        
        use constant
        implicit none

        !Dummy variables
        real::vz(nx_v,nz_v)
        integer::nsx,nsz0,nrz0
        integer::currshot_xmin,currshot_xmax
        integer::nx_v,nz_v,nt,&
                 nx_bound_l,nx_bound_r,&
                 nz_bound_u,nz_bound_d
        integer::order_2nd,order_1st
        integer::myid
        real::dx,dz,dt,f0
        integer::currshot_no
        character(len=256)::shot_fn1,shot_fn2,shot_fn3,shot_fn4,currshot_name
        !Local variables

        !*Buffers for finite_difference*
        real,allocatable::u1(:,:)
        real,allocatable::u2(:,:)

        real,allocatable::record(:,:)
        real,allocatable::record_acc(:,:)


        real,allocatable::coe_2nd_10(:)
        real,allocatable::coe_2nd_8(:)
        real,allocatable::coe_2nd_6(:)
        real,allocatable::coe_2nd_4(:)
        real,allocatable::coe_2nd_2(:)
        real,allocatable::coe_1st(:)



        real,allocatable::currshot_vz(:,:)!**zy

        real,allocatable::top(:,:)!**wcl
        real,allocatable::bot(:,:)!**wcl
        real,allocatable::lef(:,:)!**wcl
        real,allocatable::rig(:,:)!**wcl
        real,allocatable::u3(:,:)!**wcl
        real,allocatable::u4(:,:)!**wcl
		integer::ix,iz!**wcl

        !**Buffers for PML wavefield**
     

        !***Buffers for p wave***
        real,allocatable::u1_1_x(:,:),u1_2_x(:,:),&
                          u2_1_x(:,:),u2_2_x(:,:),&
                          u2_tmp_1_x(:,:),u2_tmp_2_x(:,:),&
                          u3_1_x(:,:),u3_2_x(:,:),&
                          u1_1_z(:,:),u1_2_z(:,:),&
                          u2_1_z(:,:),u2_2_z(:,:),&
                          u2_tmp_1_z(:,:),u2_tmp_2_z(:,:),&
                          u3_1_z(:,:),u3_2_z(:,:)
                          
    
        !*Variables wavefield extrapolation*
        integer::currshot_range,currshot_range_all,nz,&
                  currshot_nsx,currshot_nsz
        !*Other local variables*
        integer::i,i_ori,j,k,it,iit,err
        integer::nwt,decay
		character(len=256)::currt,currtsnap
		character(len=256)::snap_fn1=&
!		'/data4/wcl/result/acoustic_modelling/2d_acoustic_modelling/sigsbee/snap/iso_'
		'/data3/wcl/layer_two/2d/data/snap1/iso_'
		character(len=256)::snap_fn2='./snap2/iso_'!**wcl
		character(len=256)::snap_fn3='./snap3/iso_'!**wcl
		character(len=256)::snap_fn4='./snap4/iso_'!**wcl
		character(len=256)::fnn1='top'!**wcl
		character(len=256)::fnn2='bot'!**wcl
		character(len=256)::fnn3='lef'!**wcl
		character(len=256)::fnn4='rig'!**wcl
        integer::counter=1



        !Computing currshot range
        currshot_range=(currshot_xmax-currshot_xmin)+1
     !   currshot_range_all=currshot_range+(nx_bound_l+nx_bound_r)-2  
      !  nz=nz_v+(nz_bound_u+nz_bound_d)-2!**zy

        currshot_range_all=currshot_range+(nx_bound_l+nx_bound_r)  
        nz=nz_v+(nz_bound_u+nz_bound_d)!**wcl
        !compute the length of wavelet in time dimension
        nwt=2*nint(1.0/(f0*dt))
		write(*,*)'currshot_range',currshot_range
		write(*,*)'currshot_range_all',currshot_range_all
		write(*,*)'nz_v=',nz_v,'nz=',nz

        !Allocate memories for buffers
        allocate(u1(-4:currshot_range_all+5,-4:nz+5),STAT=err)
        allocate(u2(-4:currshot_range_all+5,-4:nz+5),STAT=err)

        allocate(record(nt,currshot_range),STAT=err)!**new
        allocate(record_acc(nt,currshot_range),STAT=err)!**new
        allocate(coe_2nd_10(order_2nd/2),STAT=err)
        allocate(coe_2nd_8(4),STAT=err)
        allocate(coe_2nd_6(3),STAT=err)
        allocate(coe_2nd_4(2),STAT=err)
        allocate(coe_2nd_2(1),STAT=err)

        allocate(coe_1st((order_1st+1)/2),STAT=err)
        allocate(currshot_vz(currshot_range_all+10,nz+10))

        allocate(top(currshot_range,nt))!**wcl
        allocate(bot(currshot_range,nt))!**wcl
        allocate(lef(nz_v,nt))!**wcl
        allocate(rig(nz_v,nt))!**wcl
        allocate(u3(-4:currshot_range_all+5,-4:nz+5),STAT=err)
        allocate(u4(-4:currshot_range_all+5,-4:nz+5),STAT=err)

        !*Buffers for PML wavefield*
        !**Buffers for p wave**
        !***X direction***
      !  allocate(u1_1_x(1:nx_bound_l+nx_bound_r,1:nz))
        allocate(u1_1_x(1:nx_bound_l+nx_bound_r,1:nz),STAT=err)
        allocate(u1_2_x(1:nx_bound_l+nx_bound_r,1:nz))
        allocate(u2_1_x(1:nx_bound_l+nx_bound_r,1:nz))
        allocate(u2_2_x(1:nx_bound_l+nx_bound_r,1:nz))
        allocate(u2_tmp_1_x(1:nx_bound_l+nx_bound_r,1:nz))
        allocate(u2_tmp_2_x(1:nx_bound_l+nx_bound_r,1:nz))
        allocate(u3_1_x(1:nx_bound_l+nx_bound_r,1:nz))
        allocate(u3_2_x(1:nx_bound_l+nx_bound_r,1:nz))
        !***Z direction***
        allocate(u1_1_z(1:currshot_range_all,1:nz_bound_u+nz_bound_d))
        allocate(u1_2_z(1:currshot_range_all,1:nz_bound_u+nz_bound_d))
        allocate(u2_1_z(1:currshot_range_all,1:nz_bound_u+nz_bound_d))
        allocate(u2_2_z(1:currshot_range_all,1:nz_bound_u+nz_bound_d))
        allocate(u2_tmp_1_z(1:currshot_range_all,1:nz_bound_u+nz_bound_d))
        allocate(u2_tmp_2_z(1:currshot_range_all,1:nz_bound_u+nz_bound_d))
        allocate(u3_1_z(1:currshot_range_all,1:nz_bound_u+nz_bound_d))
        allocate(u3_2_z(1:currshot_range_all,1:nz_bound_u+nz_bound_d))



	
        if (err.ne.0)then
          write(*,*)'Allocate memories fails,please check'
          stop
        endif
    
        !Initiallizations
        !*'Zeros'the buffer*
        u1=0.0
        u2=0.0
        record=0.0
        record_acc=0.0
         



        u1_1_x=0.0
        u1_2_x=0.0
        u2_1_x=0.0
        u2_2_x=0.0
        u2_tmp_1_x=0.0
        u2_tmp_2_x=0.0
        u3_1_x=0.0
        u3_2_x=0.0
    
        u1_1_z=0.0
        u1_2_z=0.0
        u2_1_z=0.0
        u2_2_z=0.0
        u2_tmp_1_z=0.0
        u2_tmp_2_z=0.0
        u3_1_z=0.0
        u3_2_z=0.0
    
        coe_2nd_10=0.0
        coe_2nd_8=0.0
        coe_2nd_6=0.0
        coe_2nd_4=0.0
        coe_2nd_2=0.0

        coe_1st=0.0
        currshot_vz=0.0

		top=0.0	!**wcl
		bot=0.0 !**wcl
		lef=0.0 !**wcl
		rig=0.0 !**wcl
		u3=0.0	!**wcl
		u4=0.0	!**wcl

        !Obtain the coefficients for explicit FD extrapolation 
        call coefficient_2nd(order_2nd,coe_2nd_10)

        call coefficient_2nd(8,coe_2nd_8)
        call coefficient_2nd(6,coe_2nd_6)
        call coefficient_2nd(4,coe_2nd_4)
        call coefficient_2nd(2,coe_2nd_2)

        call coefficient_1st(order_1st,coe_1st)


        !Get currshot velocity&anisotropic parameters
        call get_currshot_parameters(vz,&
             currshot_vz,nx_v,nz_v,&
             currshot_range,currshot_range_all,currshot_xmin,currshot_xmax,&
             nz,nx_bound_l,nx_bound_r,nz_bound_u,nz_bound_d,nsx,nsz0,&
             currshot_nsx,currshot_nsz)

        !Test
        if (myid==0)then
          write(*,*)'currshot_range,currshot_range_all,currshot_nsx,nz,nsz0'
          write(*,*)currshot_range,currshot_range_all,currshot_nsx,nz,nsz0
        endif


!		decay=nint((order_2nd*dz/currshot_vz(1,1))/dt)+nwt/2
		decay=nint(nwt/2.0)
		print*,'decay=',decay


        do it=1,nt+decay
          if(modulo(it,2)==1)then


			call add_source(u2,currshot_nsx,currshot_nsz,currshot_range_all,nz,& 
							f0,it,dt)


            call extrapolation_one_step(currshot_vz,&
                                          currshot_nsx,nsz0,nrz0,currshot_range,nz_v,currshot_range_all,&
                                          nz,dx,dz,nt,dt,f0,nx_bound_l,nx_bound_r,nz_bound_u,nz_bound_d,&
                                          order_2nd,order_1st,&
                                          u1,u2,&
										  coe_2nd_10,coe_2nd_8,coe_2nd_6,coe_2nd_4,coe_2nd_2,coe_1st,&
                                          u1_1_x,u1_2_x,u2_1_x,u2_2_x,&
                                          u2_tmp_1_x,u2_tmp_2_x,u3_1_x,u3_2_x,&
                                          u1_1_z,u1_2_z,u2_1_z,u2_2_z,&
                                          u2_tmp_1_z,u2_tmp_2_z,u3_1_z,u3_2_z&
                                          )


            iit=it-decay


			if(iit>=1)then
			

				do ix=1,currshot_range
				      top(ix,iit)=u1(nx_bound_l+ix,nz_bound_u+1)
					  bot(ix,iit)=u1(nx_bound_l+ix,nz_bound_u+nz_v)
				enddo
				do iz=1,nz_v
				      lef(iz,iit)=u1(nx_bound_l+1,nrz0+nz_bound_u+iz)
					  rig(iz,iit)=u1(nx_bound_l+currshot_range,nrz0+nz_bound_u+iz)
				enddo

               record(iit,1:currshot_range)=&
              u1(nx_bound_l+1:currshot_range+nx_bound_l,nrz0+nz_bound_u)   

				call get_acc_record(u1,record_acc,order_1st,coe_1st,&
									currshot_range,currshot_range_all,nz_v,nz,&
									nx_bound_l,nz_bound_u,nrz0,nt,iit,dz)

			endif 



        else




			call add_source(u1,currshot_nsx,currshot_nsz,currshot_range_all,nz,& 
							f0,it,dt)



            call  extrapolation_one_step(currshot_vz,&
                                          currshot_nsx,nsz0,nrz0,currshot_range,nz_v,currshot_range_all,&
                                          nz,dx,dz,nt,dt,f0,nx_bound_l,nx_bound_r,nz_bound_u,nz_bound_d,&
                                          order_2nd,order_1st,&
                                          u2,u1,&
										  coe_2nd_10,coe_2nd_8,coe_2nd_6,coe_2nd_4,coe_2nd_2,coe_1st,&
                                          u1_2_x,u1_1_x,u2_2_x,u2_1_x,&
                                          u2_tmp_2_x,u2_tmp_1_x,u3_2_x,u3_1_x,&
                                          u1_2_z,u1_1_z,u2_2_z,u2_1_z,&
                                          u2_tmp_2_z,u2_tmp_1_z,u3_2_z,u3_1_z&
                                          )
    
    

 			iit=it-decay

			if(iit>=1)then

				do ix=1,currshot_range
				  !    top(ix,iit)=u1(nx_bound_l+ix,nrz0+nz_bound_u+1)
				      top(ix,iit)=u2(nx_bound_l+ix,nz_bound_u+1)
					  bot(ix,iit)=u2(nx_bound_l+ix,nz_bound_u+nz_v)
				enddo
				do iz=1,nz_v
				      lef(iz,iit)=u2(nx_bound_l+1,nrz0+nz_bound_u+iz)
					  rig(iz,iit)=u2(nx_bound_l+currshot_range,nrz0+nz_bound_u+iz)
				enddo

				record(iit,1:currshot_range)=&
                u2(nx_bound_l+1:currshot_range+nx_bound_l,nrz0+nz_bound_u)   

				call get_acc_record(u2,record_acc,order_1st,coe_1st,&
									currshot_range,currshot_range_all,nz_v,nz,&
									nx_bound_l,nz_bound_u,nrz0,nt,iit,dz)
			endif




		endif

 		if (myid==0)then
          	if(modulo(it,10)==0)then
              write(*,*)'shot is',trim(currshot_name),'  myid=',myid,'nt=',nt,'it=',it

				write(currt,'(I6)')it
				write(currtsnap,*)trim(adjustl(snap_fn1)),trim(adjustl(currt))


				open(unit=18,file=currtsnap,form='unformatted',access='direct',&
				status='replace',recl=nz_v)

				do i=1,currshot_range
					write(18,rec=i)u2(i+nx_bound_l,nz_bound_u+1:nz_v+nz_bound_u)
!					write(8,rec=i)u2(i+nx_bound_l,1:nz)
				enddo
				close(unit=18)
          	endif
	
	!		open(unit=33,file=fnn1,form='unformatted',access='direct',&
	!			status='replace',recl=currshot_range)
	!		open(unit=34,file=fnn2,form='unformatted',access='direct',&
	!			status='replace',recl=currshot_range)
	!		open(unit=35,file=fnn3,form='unformatted',access='direct',&
	!			status='replace',recl=nz_v)
	!		open(unit=36,file=fnn4,form='unformatted',access='direct',&
	!			status='replace',recl=nz_v)
	!		do i=1,nt
	!			write(33,rec=i)top(1:201,i)
	!			write(34,rec=i)bot(1:201,i)
	!			write(35,rec=i)lef(1:201,i)
	!			write(36,rec=i)rig(1:201,i)
	!		enddo
	!		close(unit=33)
	!		close(unit=34)
	!		close(unit=35)
	!		close(unit=36)

 		endif


	enddo
    

    !*************************************************************************************  
      
        write(*,*)'shot  ',trim(adjustl(currshot_name)),&
        '  extrapolation finished,and now is writing disk'


        call write_currshot_disk(record,shot_fn1,record_acc,shot_fn3,&
								currshot_name,currshot_no,&!**zy
                                 nsx,nsz0,currshot_range,currshot_xmin,&
                                 currshot_xmax,nt,dx,dz,dt)


        write(*,*)'shot  ',trim(adjustl(currshot_name)),'  is done'
    
        deallocate(u1,STAT=err)
        deallocate(u2,STAT=err)
        deallocate(coe_2nd_10,STAT=err)
        deallocate(coe_2nd_8,STAT=err)
        deallocate(coe_2nd_6,STAT=err)
        deallocate(coe_2nd_4,STAT=err)
        deallocate(coe_2nd_2,STAT=err)
        deallocate(coe_1st,STAT=err)
        deallocate(record,STAT=err)
        deallocate(record_acc,STAT=err)
        deallocate(u1_1_x,STAT=err)
        deallocate(u1_2_x,STAT=err)
        deallocate(u2_1_x,STAT=err)
        deallocate(u2_2_x,STAT=err)
        deallocate(u2_tmp_1_x,STAT=err)
        deallocate(u2_tmp_2_x,STAT=err)
        deallocate(u3_1_x,STAT=err)
        deallocate(u3_2_x,STAT=err)
        deallocate(u1_1_z,STAT=err)
        deallocate(u1_2_z,STAT=err)
        deallocate(u2_1_z,STAT=err)
        deallocate(u2_2_z,STAT=err)
        deallocate(u2_tmp_1_z,STAT=err)
        deallocate(u2_tmp_2_z,STAT=err)
        deallocate(u3_1_z,STAT=err)
        deallocate(u3_2_z,STAT=err)
        deallocate(currshot_vz,STAT=err)
    
		!*************************************************************!


        deallocate(u3,STAT=err)
        deallocate(u4,STAT=err)
        deallocate(top,STAT=err)
        deallocate(bot,STAT=err)
        deallocate(lef,STAT=err)
        deallocate(rig,STAT=err)
    
        if (err.ne.0)then
          write(*,*)'Deallocate memories fails,please check'
          stop
        endif

        return
        end subroutine single_shot_modelling


























         






        






























        

         


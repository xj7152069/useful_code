
!==================================================================
	MODULE INTERP_GLOBAL
		INTEGER,PARAMETER::Ntable=101
		INTEGER,PARAMETER::Ltable=4
		INTEGER,PARAMETER::length=1
		INTEGER ixxx,izzz
		INTEGER lx,lz
		INTEGER mx,mz
		INTEGER nvxm,nvzm
		REAL	ax,az
		REAL	tbl_v(4,Ntable),tbl_vx(4,Ntable)
	END MODULE INTERP_GLOBAL
!==================================================================
!	SUBROUTINE TRAVELTIME_2D(SS1,SS2,TT1,TT2,SLOWNESS,TIME,SLOW_45,&
!					EPSILON_VALUE,SSS,VEL,ELEV,NVX,NX,DVX,DVZ,DX,NVZ,&
!					NZ,DZ,NS_X,NS_Z,SX_COORD,VX_START,NXS_LEFT,&
!					NXS_RIGHT,ntraces,trace_start,trace_locate,&
!					SZ_COORD,ISHOT,DSTEP,NVXS,DXS,DZS,myid)

    PROGRAM test				
		USE INTERP_GLOBAL
         integer::order_even=5
        ! real::coefficients(order_even/2)
         real,allocatable::coefficients(:)
         
         call coefficient_2nd(order_even,coefficients)
         call coefficient_1st(order_even,coefficients)
         
		
		RETURN
    END 

!=========================================================================

    subroutine coefficient_2nd(order_even,coefficients)
    !use IMSL
    implicit none
     !Dummy variables
     integer::order_even
    ! real::coefficients(order_even/2)
     real::coefficients(order_even)
    !routine variables
    real,allocatable::a(:,:),b(:)
    real::i,j
    real::fact
     integer::k1
    ! k1=order_even/2
     k1=order_even
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
    
    call inv(k1)
    
    deallocate(a)
    deallocate(b)
    return
     contains
     
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine inv(k1)
      implicit none
      integer::i1,j1,k1
      real::ftmp1,ftmp5,s
      real::c(k1)
      real::co(k1)
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
      s=0.0
      do i=1,k1
        c=0
        c(i)=0.5
        
        p1=matmul(inv_a,c)
        
        co(i)=p1(1)
        s=s+co(i)
      enddo
      print *,co
      print *,"sum= ", s
       deallocate(p1)
       deallocate(temp)
       deallocate(ftmp2)
       deallocate(ftmp3)
       deallocate(inv_a)
       deallocate(ftmp4)
      end subroutine inv
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    end subroutine  coefficient_2nd

!=========================================================================
  subroutine coefficient_1st(order_even,coefficients)
    !use IMSL
    implicit none
    !Dummy variables
    integer::order_even
    real::coefficients(order_even)
    !routine variables
    real,allocatable::a(:,:),b(:)
    real::i,j
    real::fact
    integer::k1_1
   ! k1_1=((order_odd+1)/2)!Pay attenuation and need more thicking
    k1_1=(order_even)!Pay attenuation and need more thicking
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
     call inv_1st(k1_1)
     
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
           
      subroutine inv_1st(k1_1)
      implicit none
      integer::i1,j1,k1_1
      real::ftmp1,ftmp5,s
      real::c(k1_1)
      real::co(k1_1)
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
     s=0.0
    do i=1,k1_1
      c=0
      c(i)=0.5
      p1=matmul(inv_a,c)
      co(i)=p1(1)!Pay attenuation
      s=s+co(i)
    enddo
    print *,co
    print *,"sum= ", s
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



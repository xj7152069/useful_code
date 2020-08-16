
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

!=========================================================================

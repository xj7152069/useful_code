
!***********************************************************************
    subroutine Laplace_time_2d(image_all, nz, nvx, dx, dz)

	implicit none
	!Data dictionary: declare variable types, definitons, & units
	!Dummy Variables
    integer ::	nz, nvx
    real    ::	dx, dz
    real    ::	image_all(nz, nvx)
	!Local Variables
	integer ::	ierr
	integer ::	iflag
	integer ::	ix, iz
	integer ::	i

    real    ::	b(5)
	real	::	dx2
	real	::	dz2
	real	::	a, a2, a21
	real	::	tmp1, tmp2, tmp3
	real,allocatable :: image_lap(:,:)

	!*======================================================================*
	allocate(image_lap(nz,nvx),stat=ierr)
	if(ierr .ne. 0)then
		write(*,*)"can not allocate working memory, stop!!!"
		stop
	endif
    image_lap=0.0

	iflag=1
    dx2=dx*dx
    dz2=dz*dz
    !second derivative coefficient
    b(1)=1.666667
    b(2)=-0.238095
    b(3)=0.03968254
    b(4)=-0.00496
    b(5)=0.00031746

    a=cos(60.0*3.1415926/180.0)	!pay attention 
    a2=a*a
    a21=1.0/a2

	if(iflag .eq. 1)then
		print*,	'dx2=', dx2, 'dz2=', dz2
		print*, 'a21=', a21
		print*,	'nvx=', nvx, 'nz=', nz
	endif
    do ix=6, nvx-5
		do iz=6, nz-5
			tmp1=0.0
			do i=1,5
				tmp1=tmp1+b(i)*(image_all(iz,ix+i)-&
					2.0*image_all(iz,ix)+image_all(iz,ix-i))/dx2
			enddo 
			tmp2=0.0
			do i=1,5
				tmp2=tmp2+b(i)*(image_all(iz+i,ix)-&
					2.0*image_all(iz,ix)+image_all(iz-i,ix))/dz2
			enddo
			image_lap(iz,ix)=0.25*a21*(tmp1+tmp2)
		enddo
    enddo
    image_all=-image_lap

	deallocate(image_lap)

	return
    end subroutine

!***********************************************************************


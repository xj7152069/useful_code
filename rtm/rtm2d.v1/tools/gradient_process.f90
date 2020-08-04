

!***********************************************************************
	program Gradient_Process

	implicit none
	!Data dictionary: declare variable types, definitons, & units
	!Dummy variables
	integer ::	myid
	integer ::	iter
	integer ::	nz
	integer ::	nvx
	integer	::	nshot
	real	::	dx, dz
	real	::	dt
	real,allocatable	::	gradient(:, :)

	!local variables
	integer	::	ix
	integer ::	iz
	integer ::	ierr
	integer ::	irec
	integer	::	iflag
	integer	::	ishot
	character(len=46)  ::   fn_num

	character(len=256)	::  currt_all
	character(len=128)  ::  fn_in= &
	!	'/data3/wcl/rtm/2d/imaging/waxian_301_test1'
	!	'/data3/wcl/waxian_model/result/rtm/image_test222'
		'/data3/wcl/waxian_model/result/rtm/image_test_wrong_data'

	character(len=128)  ::  fn_out= &
	!	'/data3/wcl/rtm/2d/imaging/waxian_301_test1_lap'
	!	'/data3/wcl/waxian_model/result/rtm/image_test222_lap'
		'/data3/wcl/waxian_model/result/rtm/image_test_wrong_datalap'
	iflag=0
	nshot=201
	nvx=1801
	nz=301
	dx=5.0
	dz=10.0

	allocate(gradient(nz, nvx))
	open(11, file=fn_in, access='direct', recl=nz, status='old')
	open(12, file=fn_out, access='direct', recl=nz, status='replace')

	do ishot=1, nshot
		gradient=0.0
		do ix=1, nvx
			read(11, rec=ix+(ishot-1)*nvx)(gradient(iz, ix), iz=1, nz)
		enddo

    	CALL Laplace(gradient, nz, nvx, dx, dz)

		do ix=1, nvx
			write(12, rec=ix+(ishot-1)*nvx)(gradient(iz, ix), iz=1, nz)
		enddo

	enddo
	close(11)
	close(12)

	!*======================================================================*
!	CALL Gradient_Filter(gradient, nz, nvx)

	!*======================================================================*

	end 

!***********************************************************************
	subroutine	Gradient_Filter(gradient, nz, nvx)

	implicit none
	!Data dictionary: declare variable types, definitons, & units
	!Dummy Variables
	integer ::	nz
	integer ::	nvx

	real	::	gradient(nz, nvx)

	!Local Variables
	integer	::	ix
	integer ::	iz
	integer ::	i
	real	::	trace(nz)

	do ix=1 ,nvx
		do iz=1, nz
			trace(iz)=gradient(iz, ix)
		enddo
		do i=1, 4
			CALL Hamming_Window_All(trace, 1,20, nz-10, nz, nz)
		enddo
		do iz=1, nz
			gradient(iz,ix)=trace(iz)
		enddo
	enddo

	return
	end subroutine

!***********************************************************************
    subroutine Hamming_Window_All(ctrace, nw1, nw2, nw3, nw4, lt)

	implicit none
	!Data dictionary: declare variable types, definitons, & units
	real,parameter	::	pai=3.14159165359
	!Dummy Variables
	integer	::	lt
    integer ::	nw1, nw2, nw3, nw4
    real	::	ctrace(lt)

	!Local Variables
	integer	::	iw
	real	::	hammingw

    do iw=1, lt
		if(iw.ge.nw1.and.iw.le.nw2) then
			hammingw=0.5+0.5*cos(pai*(iw-nw1)/(nw2-nw1)-pai)
			ctrace(iw)=ctrace(iw)*hammingw
        else if(iw.ge.nw3.and.iw.le.nw4) then
			hammingw=0.5+0.5*cos(pai*(nw3-iw)/(nw4-nw3))
			ctrace(iw)=ctrace(iw)*hammingw
        else if(iw.gt.nw4.or.iw.lt.nw1) then
			ctrace(iw)=0.0
        endif
    enddo

    return
    end

!***********************************************************************
    subroutine Laplace(image_all, nz, nvx, dx, dz)

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


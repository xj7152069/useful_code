
!***********************************************************************

	subroutine	Source_and_Receiver_Position_2d(sx, sz, gx, gz, &
					ns_x, ns_z, index_recgx, index_recgz, &
					ntr, cx_min, cz_min, &
					dx, dz, &
					nx_shift, nz_shift, &
					isok)

	use global
	implicit none
	real	::	sx, sz
	real	::	gx(ntr)
	real	::	gz(ntr)
	integer	::	ns_x, ns_z
	integer	::	index_recgx(ntr)
	integer	::	index_recgz(ntr)
	integer	::	ntr
	real	::	cx_min, cz_min
	real	::	dx, dz
	integer	::	nx_shift, nz_shift
	integer	::	isok

	!Local Variables
	integer	::	isxv, iszv, igxv, igzv
	integer	::	isx, isz, igx, igz
	integer	::	nr_x, nr_z
	integer	::	itr

	!*震源位置
    isx = (sx - cx_min)/dx+1.5
	ns_x= isx + nx_shift

    isz = (sz - cz_min)/dz+1.5
	ns_z= isz + nz_shift
	

	!*检波器位置
	index_recgx=-1
	index_recgz=-1

	do itr=1, ntr
		igx = (gx(itr)-cx_min)/dx+1.5
		igz = (gz(itr)-cz_min)/dz+1.5

		nr_x= igx + nx_shift
		nr_z= igz + nz_shift

		index_recgx(itr)=nr_x
		index_recgz(itr)=nr_z
				
	enddo

	isok=0

	end subroutine

!***********************************************************************

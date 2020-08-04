
!***********************************************************************
    subroutine Get_Current_Shot_Velocity_2d(vel, vv, &
					nvz, nvx, nnz, nnx, &
					nz_bound_u, nx_bound_l, &
					nvzz_shift, nvxx_shift, &
					isok)

	use global
	implicit none
	!Data dictionary: declare variable types, definitons, & units
	!Dummy Variables
	real	::	vel(nvz, nvx)
    real	::	vv(nnz, nnx)
	
	integer	::	nvz, nvx
	integer	::	nnz, nnx
	integer	::	nz_bound_u, nx_bound_l
	integer	::	nvzz_shift, nvxx_shift
	integer	::	isok

	!Local Variables
	integer ::	ix, iix, ivx 
	integer ::	iz, iiz, ivz



		do iix=1, nnx
			ivx = iix - nx_bound_l + nvxx_shift
			if(ivx .lt. 1)	ivx=1
			if(ivx .gt. nvx)	ivx=nvx

			do iiz=1, nnz
				ivz = iiz -nz_bound_u + nvzz_shift
				if(ivz .lt. 1)	ivz=1
				if(ivz .gt. nvz)	ivz=nvz
				vv(iiz, iix)=vel(ivz, ivx)
			enddo
		enddo



	isok=0
    return
    end subroutine

!***********************************************************************

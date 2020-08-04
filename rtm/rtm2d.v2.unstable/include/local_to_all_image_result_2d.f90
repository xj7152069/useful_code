



!***********************************************************************
    subroutine local_to_all_image_result_2d(image_vel, image, nz_with_apert, nx_with_apert, nvz, nvx, &
					nx_apert_l, nz_apert_u, &
					nvxx_shift, nvzz_shift, &
					isok)

	use global
	implicit none
	!Data dictionary: declare variable types, definitons, & units
	!Dummy Variables
	real	::	image_vel(nvz, nvx)
	real	::	image(nz_with_apert, nx_with_apert)
	integer	::	nz_with_apert, nx_with_apert
	integer	::	nvz, nvx
	integer	::	nx_apert_l, nz_apert_u
	integer	::	nvxx_shift, nvzz_shift
	integer	::	isok

	!Local Variables
	integer ::	ix, iz
	integer ::	ivx, ivz
    integer	::	irec
	integer	::	shift_x, shift_z


	
!	image=1.0
!	shift_x = nvxx_shift + nx_apert_l
!	shift_z = nvzz_shift + nz_apert_u
	shift_x = nvxx_shift 
	shift_z = nvzz_shift 

    !$omp parallel private(ivx, ivz)
	!$omp do
    		do ix=1, nx_with_apert
        		ivx= ix + shift_x
        		if(ivx .ge. 1 .and. ivx .le. nvx) then
            		do iz=1, nz_with_apert
						ivz= iz + shift_z
						if(ivz .ge. 1 .and. ivz .le. nvz)then
							image_vel(ivz, ivx)=image_vel(ivz, ivx) + image(iz, ix)
						endif
      		    	enddo
        		endif
    		enddo
	!$omp enddo
	!$omp end parallel


	isok=0
    return

    end	subroutine

!***********************************************************************


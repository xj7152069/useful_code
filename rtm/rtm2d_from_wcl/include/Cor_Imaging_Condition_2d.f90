
	subroutine Cor_Imaging_Condition_2d(image, vv, us1, us2, us3, ur1, ur2, ur3,&
							nnz, nnx, nz_with_apert, nx_with_apert, &
							nz_bound_u, nx_bound_l ,dt)

	implicit none
	!Data dictionary: declare variable types, definitons, & units
	!Dummy Variables
	real	::	image(nz_with_apert, nx_with_apert)
	real	::	vv(nnz, nnx)
	real	::	us1(-4:nnz+5, -4:nnx+5)
	real	::	us2(-4:nnz+5, -4:nnx+5)
	real	::	us3(-4:nnz+5, -4:nnx+5)
	real	::	ur1(-4:nnz+5, -4:nnx+5)
	real	::	ur2(-4:nnz+5, -4:nnx+5)
	real	::	ur3(-4:nnz+5, -4:nnx+5)


	integer	::	nnz, nnx
	integer	::	nz_with_apert, nx_with_apert
	integer	::	nz_bound_u, nx_bound_l
	real	::	dt

	!Local Variables
	integer	::	ix, ixx
	integer ::	iz, izz 
	real	::	tmp, tmp1


    !$omp parallel private(ixx, izz, tmp1, tmp)
	!$omp do
		do ix=1, nx_with_apert
		!	ixx = ix + nx_shift
			ixx = ix + nx_bound_l
			do iz=1, nz_with_apert
			!	izz=iz + nz_shift
				izz=iz + nz_bound_u
			
				!**********************************************************!
				tmp1=(us1(izz, ixx)+us3(izz, ixx)-2*us2(izz, ixx))
				tmp1=tmp1/(dt*dt)
				tmp=tmp1*ur3(izz, ixx)
				image(iz, ix)=image(iz, ix) + tmp

				!**********************************************************!
	   	     enddo
		enddo
	!$omp enddo
	!$omp end parallel

    return
    end	subroutine


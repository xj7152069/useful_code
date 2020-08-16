
!***********************************************************************
	subroutine	process_directwave_2d(shotobs, misfit_sum, &
					sx, sz, gx, gz, index_recgx, index_recgz, &
					vsurf, lt, ntr, dt, nwt)

	use global
	implicit none
	!Data dictionary: declare variable types, definitons, & units
	!Dummy Variables
	real	::	shotobs(lt, ntr)
	real	::	misfit_sum
	real	::	gx(ntr)
	real	::	gz(ntr)
	real	::	sx, sz
	integer	::	index_recgx(ntr)
	integer	::	index_recgz(ntr)
	integer	::	lt, ntr, nwt
	real	::	dt
	real	::	vsurf


	!Local Variables
	integer	::	itr, it
	integer	::	flag_mod, flag_offset
	real	::	tmp1, tmp2
	
	!flag_mod 
	real	::	coe1, coe2
	real	::	tmax
	real	::	taper, taper1
	real	::	recgx, recgz
	real	::	offset, offset_min, offset_max
	real	::	tfb, t1, t2
	real	::	win
	integer	::	it1, it2, nwin
	real	::	dis
	real	::	coeff


	flag_mod=1
	if(flag_mod .eq. 1)then

		coe1=0.5
		coe2=0.5
		tmax=lt*dt
		taper=nwt*dt

		offset_min=999999999.9
		offset_max=-999999999.9
		do itr=1, ntr
			recgx=gx(itr)	
			recgz=gz(itr)
            offset=sqrt((recgx-sx)*(recgx-sx) + (recgz-sz)*(recgz-sz))
			if(offset .gt. offset_max)	offset_max=offset
			if(offset .lt. offset_min)	offset_min=offset
		enddo
		coeff = (3.0*taper -taper)/(offset_max - offset_min)


		do itr=1, ntr
			recgx=gx(itr)	
			recgz=gz(itr)
            offset=sqrt((recgx-sx)*(recgx-sx) + (recgz-sz)*(recgz-sz))

			taper1= 3.0*taper - (offset-offset_min)*coeff 
			tfb = offset/vsurf + taper1
			t1  = tfb - taper*0.5
			t2  = tfb + taper*0.5
			if(t1 .gt. tmax)	t1=tmax
			if(t2 .gt. tmax)	t2=tmax

			it1=t1/dt+1
			it2=t2/dt+1
			nwin=it2-it1

			do it=1, lt
				if(it .lt. it1)	then
					win=0.0
				else if(it .lt. it2) then
					dis=1.0*(it-it1)/nwin
					win=coe1-coe2*cos(pi*dis)
				else
					win=1.0
				endif
				shotobs(it, itr)=shotobs(it, itr)*win
			enddo
		enddo
		
	endif


	flag_offset=1
	offset_max = 1000.0
	if(flag_offset .eq. 1)then
		do itr=1, ntr
			recgx=gx(itr)	
			recgz=gz(itr)
            offset=sqrt((recgx-sx)*(recgx-sx) + (recgz-sz)*(recgz-sz))

			if(offset .gt. offset_max)then
				do it=1, lt
					shotobs(it, itr)=0.0
				enddo
			endif
		enddo
	endif


	misfit_sum=0.0
	do itr=1, ntr
		do it=1, lt
			misfit_sum=misfit_sum + shotobs(it, itr)*shotobs(it, itr)
		enddo
	enddo

	end subroutine

!***********************************************************************

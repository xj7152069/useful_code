
!***********************************************************************
    subroutine Trace_Sinc_Interp(trace, lt, dt, trace_real, lt_real, dt_real)

	implicit none
	!Data dictionary: declare variable types, definitons, & units
	!Dummy Variables
	integer ::	lt, lt_real
	real	::	dt, dt_real
	real	::	trace(lt), trace_real(lt_real)

	!Local Variables
	!Data dictionary: declare constants
	integer,parameter	::	lw=4
	real,parameter	::	pi=3.14159265359
	!Data dictionary: declare variable types, definitons, & units
	integer ::	it, iit
	integer ::	iw
	integer	::	n_interped
	integer ::	it_start, it_final
	real	::	rdt
	real	::	tmp_hf
	real	::	tdif
	real	::	x
	real	::	sinc

    n_interped =nint(dt_real/dt)
    rdt=1.0/dt_real

    do it = 1, lt

		!iit=int(it/n_interped)
        !write(*,*)'iit-int'
        !write(*,*)iit
        iit = it/n_interped + 0.5
        !write(*,*)'itt-+0.5'
        !write(*,*)iit

        it_start = iit - lw
        if(it_start.lt.1) it_start = 1

        it_final = iit + lw
        if(it_final.gt.lt_real) it_final = lt_real

        tmp_hf = 0.0
        tdif = it * dt - it_start * dt_real

        do iw = it_start, it_final
			if(tdif.eq.0.0) then
				tmp_hf = tmp_hf + trace_real(iw)
			else
				x = pi * tdif * rdt
				sinc = sin(x)/x
				tmp_hf = tmp_hf + trace_real(iw)*sinc
			endif
			tdif = tdif - dt_real

		enddo
		trace(it)  = tmp_hf
	enddo		

    return
    end subroutine

!***********************************************************************


!***********************************************************************
	subroutine Wavelet_Forming(wavelet, dt, lt, fmain, nwt)
!*==========================================================================
!*  assigning the source function for the shot wavefield extrapolation
!*  fmain: the main frequency of the source function
!*  here the ricker wavelet is used.
!*  2010,2,11, modified by herb.wang
!*  richer wavlet function: f(t) = (1.0-2.0*(pai*fmain*t)**2)exp(-(pai*fmain*t)**2)
!*==========================================================================

	implicit none
	!Dummy Variables
	!Data dictionary: declare variable types, definitons, & units
	integer ::	lt
	integer ::	nwt
	real	::	dt
	real	::	fmain
	real	::	wavelet(lt)

	!Local Variables
	!Data dictionary: declare constants
	real,parameter	::	pi=3.14159265359
	!Data dictionary: declare variable types, definitons, & units
	integer it
	real	tmain
	real	tp1, tp2
	real	time

	!tmain=1.2/fmain
    tmain=1.0/fmain
    nwt=tmain/dt-1
		
    do it=1, lt
		!time = it*sdt-100*sdt  !zhaolei_20100720
		time = (it-1)*dt-tmain
		!time = (it-1)*sdt!zero phase wavelet hephaestus_9_25

		tp1 = pi*fmain*time
		tp2 = tp1*tp1
		wavelet(it) = (1.0-2.0*tp2)*exp(-tp2)!*0.1e7
    enddo

    return
    end subroutine

!***********************************************************************

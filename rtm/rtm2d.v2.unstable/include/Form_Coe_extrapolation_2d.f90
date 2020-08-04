

!***********************************************************************
    subroutine Form_Coe_extrapolation_2d(coe, dx, dz, dt)

	implicit none
	!Data dictionary: declare variable types, definitons, & units
	!Dummy Variables
    real	::	coe(13)
    real	::	dx, dz, dt
	!Local Variables
	real	::	ddt
	real	::	dtx2
	real	::	dtz2
	real	::	omga0, omga1, omga2, omga3, omga4, omga5

	ddt=dt
      
    dtx2=ddt*ddt/dx/dx/2.0
    dtz2=ddt*ddt/dz/dz/2.0

    omga0=-5.8544444444
    omga1=+3.3333333333
    omga2=-0.4761904762
    omga3=+0.0793650794
    omga4=-0.0099206349
    omga5=+0.0006349206

    coe(1) =+2
    coe(2) =-1
    coe(3) =omga0*(dtx2+dtz2)
    coe(4) =omga1*dtx2
    coe(5) =omga2*dtx2
    coe(6) =omga3*dtx2
    coe(7) =omga4*dtx2
    coe(8) =omga5*dtx2
    coe(9) =omga1*dtz2
    coe(10)=omga2*dtz2
    coe(11)=omga3*dtz2
    coe(12)=omga4*dtz2
    coe(13)=omga5*dtz2

    return
    end	subroutine
	
!***********************************************************************

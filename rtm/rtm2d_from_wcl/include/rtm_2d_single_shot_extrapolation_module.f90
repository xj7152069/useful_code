



!***********************************************************************
	  module	module_rtm_extrapolation_2d

	  implicit none

    	!wavefield 
	  	type wavefield
	  		real,allocatable	::	u(:, :)
	  	end	type

		type(wavefield), target	::	us1, us2, us3
		type(wavefield), target	::	ur1, ur2, ur3

	  	type(wavefield), target	::	ux_a1, ux_a2, ux_a3
	  	type(wavefield), target	::	ux_b1, ux_b2
	  	type(wavefield), target	::	ux_bp1, ux_bp2, ux_bp3
	  	type(wavefield), target	::	ux_c1, ux_c2, ux_c3

	  	type(wavefield), target	::	uz_a1, uz_a2, uz_a3
	  	type(wavefield), target	::	uz_b1, uz_b2
	  	type(wavefield), target	::	uz_bp1, uz_bp2, uz_bp3
	  	type(wavefield), target	::	uz_c1, uz_c2, uz_c3


	  	type(wavefield), pointer	::	pus1, pus2, pus3
	  	type(wavefield), pointer	::	pur1, pur2, pur3

	  	type(wavefield), pointer	::	pux_a1, pux_a2, pux_a3
	  	type(wavefield), pointer	::	pux_b1, pux_b2
	  	type(wavefield), pointer	::	pux_bp1, pux_bp2, pux_bp3
	  	type(wavefield), pointer	::	pux_c1, pux_c2, pux_c3


	  	type(wavefield), pointer	::	puz_a1, puz_a2, puz_a3
	  	type(wavefield), pointer	::	puz_b1, puz_b2
	  	type(wavefield), pointer	::	puz_bp1, puz_bp2, puz_bp3
	  	type(wavefield), pointer	::	puz_c1, puz_c2, puz_c3

	  	type(wavefield), pointer	::	pt

	  	!wavefield3
	  	type wavefield3
	  		real,allocatable	::	u(:,:,:)
	 	end type

	  	type(wavefield3), target	::	top, bot, lef, rig
	  	type(wavefield3), target	::	ust

	  	real,allocatable	::	image_tmp(:,:)
	  	
	  	!pml
		real,allocatable	::	coe_1st(:), coe_2nd(:)
		real,allocatable	::	coe_1st_dx(:), coe_1st_dz(:)
		real,allocatable	::	coe_2nd_dx2(:), coe_2nd_dz2(:)
		real,allocatable	::	funa(:,:), dfuna(:,:)
		integer,allocatable	::	line_xl(:,:), line_xr(:,:)
		integer,allocatable	::	line_zu(:,:), line_zd(:,:)

	  	

	  end module 	module_rtm_extrapolation_2d
!***********************************************************************

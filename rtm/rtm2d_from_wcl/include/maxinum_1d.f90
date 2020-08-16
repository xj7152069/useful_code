
!*==========================================================================
	  subroutine mininum_1d(similar, np, ip_min, ss_min, isok)
	  implicit none
	  !Data dictionary: declare variable types, definitons, & units
	  !Dummy variables
	  real		:: similar(np)
	  integer	::	np
	  integer	::	ip_min
	  real		::	ss_min
	  integer	::	isok

	  !Local variables
	  integer	::	ip
	  real		::	tmp

	  tmp=similar(1)
	  ip_min= 1
!debugwcl
!	  	open(212, file='test_adjust_nt', access='direct', recl=ntrace, status='replace')
!	  	write(212, rec=1)(it_adjust(ix), ix=1, ntrace)
!	  	close(212)

	  do ip=2, np
	  	if(abs(similar(ip)) .lt. tmp)then
			tmp=abs(similar(ip))
	  		ip_min=ip
	  	endif
	  enddo

	  ss_min=tmp

	  isok=0
	  return
	  end subroutine
!*==========================================================================
	  subroutine maxinum_1d(similar, np, ip_max, ss_max, isok)
	  implicit none
	  !Data dictionary: declare variable types, definitons, & units
	  !Dummy variables
	  real		:: similar(np)
	  integer	::	np
	  integer	::	ip_max
	  real		::	ss_max
	  integer	::	isok

	  !Local variables
	  integer	::	ip
	  real		::	tmp

	  tmp=similar(1)
	  ip_max= 1

	  do ip=2, np
	  	if(abs(similar(ip)) .gt. tmp)then
			tmp=abs(similar(ip))
	  		ip_max=ip
	  	endif
	  enddo

	  ss_max=tmp

	  isok=0
	  return
	  end subroutine
!*==========================================================================

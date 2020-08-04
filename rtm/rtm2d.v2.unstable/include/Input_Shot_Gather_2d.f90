
!***********************************************************************
	subroutine	Input_Shot_Gather_2d(shot_gather, shotreceiver, &
						sx, sz, gx, gz, &
						ishot, shotfldr, ntrace, ns_total, &
						lt, lt_real, dt, dt_real, isok)

!*==========================================================================
!*
!*     ns_x: the shot point position corresponding to the current shot and 
!*                 receiver range
!*     lt_real, dt_real: the sample point number and the sample rate from the 
!*                       inputed gather
!*
!*==========================================================================
	use global
	use header_module	
	use shotgather_info
	implicit none
	!Data dictionary: declare variable types, definitons, & units
	!Dummy Variables

	type(shotinfo)	::	shotreceiver(ns_total)
    real    ::	shot_gather(lt, ntrace)
	real	::	gx(ntrace), gz(ntrace)
    real	::	sx, sz
	integer	::	ishot, shotfldr
	integer	::	ntrace, ns_total
	integer	::	lt, lt_real
	real	::	dt, dt_real
	integer	::	isok

	!Local Variables
	character(len=para_char_flen)	::	fn_cs
	type(segy)  ::  head
	integer*8	::	krec_save
	integer*8	::	irec
	integer	::	shotfldr_tmp, ntrace_tmp
	real	::	sx_tmp, sz_tmp
	real	::	gx_tmp, gz_tmp
	real	::	sx_cur, sz_cur
	real	::	gx_cur, gz_cur
	integer	::	shot_cur
	integer ::	interp_flag
	integer	::	ierr
	integer	::	itrace
	integer ::	it
	
	real,allocatable	::	trace(:), trace_real(:)


	allocate(trace(lt))
	allocate(trace_real(lt_real))
	trace =0.0
	trace_real =0.0


	fn_cs=trim(shotreceiver(ishot)%filename)
	shotfldr_tmp =shotreceiver(ishot)%fldr
	ntrace_tmp =shotreceiver(ishot)%ntr
	krec_save =shotreceiver(ishot)%krec
	sx_tmp = shotreceiver(ishot)%sx
	sz_tmp = shotreceiver(ishot)%sz

	if(shotfldr .ne. shotfldr_tmp)then
		write(*,*) 'Reading the shotfldr is not right.'
		stop
	endif
	if(ntrace .ne. ntrace_tmp)then
		write(*,*) 'Reading the ntrace is not right.'
		stop
	endif

	!*打开炮集文件
	open(11, file=fn_cs, access='direct', recl=lbyte*(60+lt_real), status='old', iostat=ierr)
	if(ierr /= 0)then
		 write(*,*)  'Shot gather file cannot open right, please check it !!!'
		 write(*,*)  'Shutting down the program'
		 stop
	endif

    irec=krec_save

    if(modulo(dt_real, dt).eq.0.0.and.lt.gt.lt_real) then  !yzhou
		interp_flag = 2
    else
		interp_flag = 1
    end if

	do itrace=1, ntrace
    	read(11, rec=irec, err=5555)head,trace_real
		shot_cur=head%fldr
		sx_cur=head%sx
		gx_cur=head%gx

		if(shot_cur .ne. shotfldr)then
			write(*,*)'Warning : Read the data is error!!!'
			stop
		endif
		if(sx_cur .ne. sx)then
			write(*,*)'Warning : Read the data is error!!!'
			stop
		endif
		if(gx_cur .ne. gx(itrace))then
			write(*,*)'Warning : Read the data is error!!!'
			stop
		endif

		!*********************************************************************!
    	if(interp_flag.eq.2) then
			CALL Trace_Sinc_Interp(trace, lt, dt, trace_real, lt_real, dt_real)
    	else if(interp_flag.eq.1) then
			do it=1, lt
				trace(it) = trace_real(it)
			enddo
    	endif

		!*---------- the shot position of current shot -------------------------*
		do it=1, lt
			shot_gather(it, itrace)=trace(it)
    	enddo

		!*----------------------------------------------------------------------*
    	irec=irec+1
	enddo

5555	continue	
	
	deallocate(trace)
	deallocate(trace_real)

	return
    end	subroutine

!***********************************************************************


!***********************************************************************
    subroutine Read_shotgathers_info(fn_shot_info, &
					ns_start, ns_end, ns_interval, ns_total, nshot, &
					shotreceiver, ntrace_max, lt_real, &
					num_grad, myid, isok)

	use global
	use shotgather_info
	implicit none
	!Data dictionary: declare variable types, definitons, & units
	!Dummy Variables
	character(len=para_char_flen)	::	fn_shot_info

	type(shotinfo)	::	shotreceiver(ns_total)
	integer	::	ns_start, ns_end, ns_interval
	integer	::	ns_total, nshot
	integer	::	ntrace_max
	integer	::	lt_real
	integer	::	num_grad
	integer	::	myid
	integer	::	isok

	!Local Variables
	integer	::	ierr
	integer	::	ishot, itrace
	integer	::	irec
	integer	::	ii
	real	::	xx, yy, zz




	open(21, file=trim(adjustl(fn_shot_info))//'.shotreceiver', status='old', iostat=ierr)
	if(ierr /= 0)then
		write(*,*)	'Shotreceiver file cannot open right, please check it !!!'
		write(*,*)	'Shutting down the program'
		stop
	endif

	ishot=0
	do while (.TRUE.)
		read(21, *, iostat=ierr)ii, xx, yy, zz
        if(ierr /= 0)then
        	exit
		endif
        if(ii .ne. -1)then
            ishot = ishot +1

			if(ishot .gt. nshot)	exit
       		shotreceiver(ishot)%fldr = ii
      		shotreceiver(ishot)%sx = xx
       		shotreceiver(ishot)%sy = yy
      		shotreceiver(ishot)%sz = zz
       		itrace = 0
		!	print*, 'test myid, ishot', myid, ishot

		else if(ii  .eq. -1)then
       		itrace = itrace +1
			if(itrace .gt. ntrace_max)then
				print*, 'Read the trace error!!!!'
				stop
			endif
      		shotreceiver(ishot)%gx(itrace)	=	xx
      		shotreceiver(ishot)%gy(itrace)	=	yy
      		shotreceiver(ishot)%gz(itrace)	=	zz

		endif
	enddo
	close(21)

	print*, 'ishot, nshot', myid, ishot, nshot
	if(ishot .ne. nshot)then
		write(*,*)'Warning: The shot number is error, when reading the shot info!!!'
	endif

	nshot=ishot
	open(21, file=trim(adjustl(fn_shot_info))//'.fileinfo', status='old', iostat=ierr)
	if(ierr /= 0)then
		write(*,*)	'Shotreceiver file cannot open right, please check it !!!'
		write(*,*)	'Shutting down the program'
		stop
	endif
	do ishot=1, nshot
		read(21, *)shotreceiver(ishot)%fldr, shotreceiver(ishot)%krec, shotreceiver(ishot)%ntr
		read(21, '(A)')shotreceiver(ishot)%filename
	enddo
	close(21)

	if(myid  .eq. 0  .and. num_grad .eq. 1)then

		print*, '!=====================================================!'
		print*, 'Read the shot and receiver info :      '
		do ishot=1, nshot
			print*, '      '
			print*, 'fldr', ishot, shotreceiver(ishot)%fldr
			print*, 'filename', ishot, trim(shotreceiver(ishot)%filename)
			print*, 'krec', ishot, shotreceiver(ishot)%krec
			print*, 'sx  ', ishot, shotreceiver(ishot)%sx
			print*, 'sy  ', ishot, shotreceiver(ishot)%sy
			print*, 'sz  ', ishot, shotreceiver(ishot)%sz
			print*, 'ntr ', ishot, shotreceiver(ishot)%ntr
			print*, '      '
		enddo
		print*, '!=====================================================!'

	endif


    end	subroutine
!***********************************************************************


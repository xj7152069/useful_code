
!***********************************************************************
    subroutine index_shotgathers_and_write_info(fn_cs, fn_shot_info, &
					ns_start, ns_end, ns_interval, ns_total, nshot, &
					shotreceiver, ntrace_max, lt_real, &
					myid, isok)

	use global
	use header_module
	use shotgather_info
	implicit none
	!Data dictionary: declare variable types, definitons, & units
	!Dummy Variables
	character(len=para_char_flen)	::	fn_cs
	character(len=para_char_flen)	::	fn_shot_info

	type(shotinfo)	::	shotreceiver(ns_total)
	integer	::	ns_start, ns_end, ns_interval
	integer	::	ns_total, nshot
	integer	::	ntrace_max
	integer	::	lt_real
	integer	::	myid
	integer	::	isok

	!Local Variables
	character(len=para_char_flen)	::	fn_tmp
	type(segy)  ::  head
	real,allocatable	::	trace(:)
	integer*8	::	krec
	integer	::	ierr
	integer	::	ishot, itrace
	integer	::	ifile, nfile
	integer	::	irec
	integer	::	is
	integer	::	nss
	integer	::	icount
	integer	::	it


	allocate(trace(lt_real))

	open(33, file=fn_cs, status='old', iostat=ierr)
	if(ierr /= 0)then
		write(*,*)  'Parameter file cannot open right, please check it !!!'
		stop
	endif
		
	if(myid .eq. 0)then
		print*, 'filename cs: '
		print*, trim(fn_cs)
		print*, '!=====================================================!'
	endif

	ifile=0
	do while (.TRUE.)
		read(33, '(a)', iostat=ierr)fn_tmp
        if(ierr /= 0)then
        	exit
		endif
		if(myid .eq. 0)then
			print*, 'filename tmp: '
			print*, trim(fn_tmp)
		endif
		ifile = ifile +1
	enddo

	nfile = ifile
	if(myid .eq. 0)then
		write(*,*)'The input dada file numbers is ', nfile
		write(*,*) '         '
		print*, '!=====================================================!'
	endif


	!===========================================================================
	rewind(33)
	do ifile=1, nfile

		read(33, '(a)', iostat=ierr)fn_tmp
!		print*, 'filename tmp: ', trim(fn_tmp)
		!*打开炮集文件
		open(331, file=fn_tmp, access='direct', recl=lbyte*(60+lt_real), status='old', iostat=ierr)
		if(ierr /= 0)then
			write(*,*)	'Shot gather file cannot open right, please check it !!!'
			write(*,*)	'Shutting down the program'
			stop
		endif

		if(myid .eq. 0 .and. ifile .eq. 1)then
			print*, 'ns_start, ns_end, ns_interval', ns_start, ns_end, ns_interval
		endif

		irec=0
    	do is=ns_start, ns_end, ns_interval
			ishot = (is - ns_start)/ns_interval + 1

			if(myid .eq. 0 .and. ifile .eq. 1)then
				print*, 'is, ishot', is, ishot
			endif

       		itrace=0				!当前炮的道数
1111    	continue
        	irec=irec+1
        	read(331, rec=irec, err=5555)head, (trace(it), it=1, lt_real)
			nss=head%fldr
        	if(nss .lt. is)  then
				goto 1111

        	elseif(nss .eq. is) then
				itrace=itrace+1

				if(itrace .eq. 1)then
       				shotreceiver(ishot)%fldr = is
					shotreceiver(ishot)%filename = fn_tmp
					shotreceiver(ishot)%krec = irec
      				shotreceiver(ishot)%sx = head%sx
       				shotreceiver(ishot)%sy = head%sy
      				shotreceiver(ishot)%sz = 0.0
				endif
				if(itrace .gt. ntrace_max)then
					print*, 'Read the trace error!!!!'
					stop
				endif
				shotreceiver(ishot)%ntr = itrace
      			shotreceiver(ishot)%gx(itrace)	=	head%gx
      			shotreceiver(ishot)%gy(itrace)	=	head%gy
      			shotreceiver(ishot)%gz(itrace)	=	0.0
				goto 1111

        	elseif(nss .gt. is) then
				irec=irec-1
				goto 2222
        	endif
2222		continue
    	enddo
5555	continue
		close(331)
		!当前文件读取结束
	enddo
	close(33)


	deallocate(trace)
	!===========================================================================

	icount=0
	do ishot=1, ns_total
		if(shotreceiver(ishot)%fldr .ne. 0)	icount=icount+1
	enddo
	nshot=icount

!	if(myid .eq. 0)then

	print*, 'Total shot number : ', ns_total, ' Valid shot number : ', nshot


	open(21, file=trim(adjustl(fn_shot_info))//'.shotreceiver', status='replace')
	do ishot=1, ns_total
		if(shotreceiver(ishot)%fldr .ne. 0)then
			write(21, *)shotreceiver(ishot)%fldr, shotreceiver(ishot)%sx, shotreceiver(ishot)%sy, shotreceiver(ishot)%sz
			do itrace=1, shotreceiver(ishot)%ntr
				write(21, *)-1, shotreceiver(ishot)%gx(itrace), shotreceiver(ishot)%gy(itrace), shotreceiver(ishot)%gz(itrace)
			enddo
		endif
	enddo
	close(21)

	open(21, file=trim(adjustl(fn_shot_info))//'.fileinfo', status='replace')
	do ishot=1, ns_total
		if(shotreceiver(ishot)%fldr .ne. 0)then
			write(21, *)shotreceiver(ishot)%fldr, shotreceiver(ishot)%krec, shotreceiver(ishot)%ntr
			write(21, *)trim(shotreceiver(ishot)%filename)
		endif
	enddo
	close(21)

		!for debug
		print*, 'Write The shot and receiver info :      '
		do ishot=1, ns_total
		if(shotreceiver(ishot)%fldr .ne. 0)then
			print*, '      '
			print*, 'fldr', ishot, shotreceiver(ishot)%fldr
			print*, 'filename', ishot, trim(shotreceiver(ishot)%filename)
			print*, 'krec', ishot, shotreceiver(ishot)%krec
			print*, 'sx  ', ishot, shotreceiver(ishot)%sx
			print*, 'sy  ', ishot, shotreceiver(ishot)%sy
			print*, 'sz  ', ishot, shotreceiver(ishot)%sz
			print*, 'ntr ', ishot, shotreceiver(ishot)%ntr
			print*, '      '
		endif
		enddo

		do ishot=1, ns_total
			if(shotreceiver(ishot)%fldr .eq. 0)then
				print*, 'The absent shot is ', ishot
			endif
		enddo
		print*, '!=====================================================!'

!	endif


    end	subroutine
!***********************************************************************

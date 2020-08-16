

	 module shotgather_info
	 implicit none

	  	type	shotinfo
            character(len=256)	::	filename
			integer*8	::	krec
            integer	::	fldr
            integer	::	ntr
			integer	::	lt
			real	::	sx
			real	::	sy
			real	::	sz
			real,allocatable,dimension(:)	::	gx
			real,allocatable,dimension(:)	::	gy
			real,allocatable,dimension(:)	::	gz
		end type

	 end module shotgather_info




	 subroutine shotreceiver_position(fn_shot_info, shotreceiver, nshot, ntrace)

	 	use shotgather_info
		implicit none
	 	character(len=256)	::	fn_shot_info
	 	type(shotinfo)	shotreceiver(nshot)
		integer	::	nshot
	 	integer	::	ntrace
	 	
	 	integer	::	ierr
		integer	::	ishot
	 	integer	::	ii
	 	integer	::	itrace
	 	real	::	xx, yy, zz

	 	open(21, file=fn_shot_info, status='old')
	 	ishot=0
			print*, 'herere222'
	 	do while (.TRUE.)
			read(21, *, iostat=ierr)ii, xx, yy, zz
          	if(ierr /= 0)then
				exit
			endif
			print*, ii, xx, yy, zz
            if(ii .ge. 0)then
            	ishot = ishot +1
       			shotreceiver(ishot)%fldr = ii
      			shotreceiver(ishot)%sx = xx
       			shotreceiver(ishot)%sy = yy
      			shotreceiver(ishot)%sz = zz
				if(ishot .gt. 1)then
					shotreceiver(ishot-1)%ntr = itrace
				endif
       			itrace = 0
			else if(ii  .eq. -1)then
       			itrace = itrace +1
				print*, 'itrace', itrace
      			shotreceiver(ishot)%gx(itrace)	=	xx
      			shotreceiver(ishot)%gy(itrace)	=	yy
      			shotreceiver(ishot)%gz(itrace)	=	zz
			endif
		enddo
		shotreceiver(ishot).ntr = itrace
		close(21)

		print*, 'dfaafasfasfafasf'
	 end subroutine
	 

	 program test

	 	use shotgather_info
		implicit none
		character(len=256)	::	fn_shot_info= './test.txt'

		type(shotinfo),allocatable	::	shotreceiver(:)

		integer	::	nshot
	 	integer	::	ntrace
	 	integer	::	ishot
	 	integer	::	ii
	 	integer	::	itrace
	 	integer	::	ierr
	 	real	::	xx, yy, zz

	 	character(len=256)	::	currt1
	 	character(len=256)	::	currt2



        nshot=1
        ntrace =10000
	 	allocate(shotreceiver(nshot))
	
     	do ishot =1, nshot
      		allocate(shotreceiver(ishot)%gx(ntrace))
      		allocate(shotreceiver(ishot)%gy(ntrace))
      		allocate(shotreceiver(ishot)%gz(ntrace))

!			print*, sizeof(shotreceiver(ishot)%gx(:))
			
	 	enddo


		write(currt1, '(I8)')nshot
		write(currt2, *)trim(adjustl(fn_shot_info)),'_iter_',trim(adjustl(currt1))
		print*, trim(currt2)




	 	open(21, file=fn_shot_info, status='old')
	 	ishot=0
		!	print*, 'herere222'
	 	do while (.TRUE.)
			read(21, *, iostat=ierr)ii, xx, yy, zz
         	if(ierr /= 0)then
        		exit
			endif
		!	print*, ii, xx, yy, zz
            if(ii .ne. -1)then
            	ishot = ishot +1
				if(ishot .gt. 1)then
					shotreceiver(ishot-1)%ntr = itrace
				endif
				if(ishot .gt. nshot)	exit
				if(ishot .gt. 1)then
					shotreceiver(ishot)%krec = shotreceiver(ishot-1)%krec + itrace
				endif
				if(ishot  .eq. 1)then
					shotreceiver(ishot)%krec = 1
				endif
       			shotreceiver(ishot)%fldr = ii
      			shotreceiver(ishot)%sx = xx
       			shotreceiver(ishot)%sy = yy
      			shotreceiver(ishot)%sz = zz
       			itrace = 0
			else if(ii  .eq. -1)then
       			itrace = itrace +1
				if(itrace .gt. ntrace)then
					print*, 'ntrace errror'
					stop
				endif
		!		print*, 'ishot, itrace', ishot, itrace
      			shotreceiver(ishot)%gx(itrace)	=	xx
      			shotreceiver(ishot)%gy(itrace)	=	yy
      			shotreceiver(ishot)%gz(itrace)	=	zz
			endif
		enddo
		shotreceiver(ishot).ntr = itrace
		close(21)

		do ishot=1, nshot
			print*, '      '
			print*, ishot, shotreceiver(ishot)%fldr
			print*, ishot, shotreceiver(ishot)%krec
			print*, ishot, shotreceiver(ishot)%sx
			print*, ishot, shotreceiver(ishot)%sy
			print*, ishot, shotreceiver(ishot)%sz
			print*, ishot, shotreceiver(ishot)%ntr
			print*, '      '
		enddo



	
     	do ishot =1, nshot
      		deallocate(shotreceiver(ishot)%gx)
      		deallocate(shotreceiver(ishot)%gy)
      		deallocate(shotreceiver(ishot)%gz)
	 	enddo

	 	deallocate(shotreceiver)


			print*, 'herere1111'
!		call shotreceiver_position(fn_shot_info, shotreceiver, nshot, ntrace)

	 end

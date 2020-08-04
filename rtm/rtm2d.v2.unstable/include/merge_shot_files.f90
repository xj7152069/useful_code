

!***********************************************************************
	subroutine  merge_shot_files(fn_cs, fn1, ns_total, lt, shotreceiver)

	use global
	use shotgather_info
	implicit none
	!Data dictionary: declare variable types, definitons, & units
	!Dummy Variables
	character(len=para_char_flen)	::	fn_cs
	character(len=para_char_flen)	::	fn1
	integer	::	ns_total
	integer	::	lt
    type(shotinfo)	::	shotreceiver(ns_total)

	!Local variables
	character(len=para_char_flen)	::	currt
	real,allocatable	::	shotgather(:,:)
	integer	::	ntrace
	integer	::	ishot
	integer	::	itr, it


    open(11, file=fn_cs, access='stream', form='unformatted', status='replace')


	do ishot=1, ns_total

		ntrace=shotreceiver(ishot)%ntr
		print*, 'ishot, ntrace  merging...',  ishot, ntrace

		allocate(shotgather(lt+60, ntrace))

		write(currt, '(I8)')ishot
		open(113, file=trim(adjustl(fn1))//'_ishot'//trim(adjustl(currt))//'.su', access='direct', &
				recl=lbyte*(lt+60)*ntrace, status='old')
		read(113, rec=1)((shotgather(it, itr), it=1, lt+60), itr=1, ntrace)
		close(113)
!		close(113, status='delete',err=333)
333		continue

		write(11)((shotgather(it, itr), it=1, lt+60), itr=1, ntrace)

		deallocate(shotgather)

	enddo


	close(11)
		
	end subroutine

!***********************************************************************




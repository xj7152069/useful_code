
!*=====================================================================*
	subroutine Form_Dynamic_Filename(myid, fn4, fn4_1, ns_start)

	use global
	implicit none
	!Data dictionary: declare variable types, definitons, & units
	!Dummy Variables
	character(len=para_char_flen)	::	fn4
	character(len=para_char_flen)  ::	fn4_1
	integer	::	myid
	integer ::	ns_start
	!Local Variables
	character(len=5)	::	fn4_2
	character(len=3)	::	fn4_3
	character(len=5)	::	fn4_4
	character(len=1)	::	ns1
	character(len=2)	::	ns2
	character(len=3)	::	ns3
	character(len=4)	::	ns4
	character(len=5)	::	ns5
	character(len=1)	::	no1
	character(len=2)	::	no2
	character(len=3)	::	no3
	character(len=4)	::	no4
	character(len=5)	::	no5

	fn4_3='_id'
	if(myid.ge.1.and.myid.le.9)then
		write(no1,'(i1)') myid
		fn4_4='000'//no1
		if(ns_start.ge.1.and.ns_start.le.9)then
			write(ns1,'(i1)')ns_start
			fn4_2='0000'//ns1
			fn4=trim(fn4_1)//fn4_2//fn4_3//fn4_4
		!	write(*,*)'filename=',fn4,'   myid=',myid
		elseif(ns_start.ge.10.and.ns_start.le.99)then
			write(ns2,'(i2)')ns_start
			fn4_2='000'//ns2
			fn4=trim(fn4_1)//fn4_2//fn4_3//fn4_4
		!	write(*,*)'filename=',fn4,'   myid=',myid
        elseif(ns_start.ge.100.and.ns_start.le.999)then
			write(ns3,'(i3)')ns_start
			fn4_2='00'//ns3
			fn4=trim(fn4_1)//fn4_2//fn4_3//fn4_4
		!	write(*,*)'filename=',fn4,'   myid=',myid
        elseif(ns_start.ge.1000.and.ns_start.le.9999)then
			write(ns4,'(i4)')ns_start
			fn4_2='0'//ns4
			fn4=trim(fn4_1)//fn4_2//fn4_3//fn4_4
		!	write(*,*)'filename=',fn4,'   myid=',myid
        elseif(ns_start.ge.10000.and.ns_start.le.99999)then
			write(ns5,'(i5)')ns_start
			fn4_2=ns5
			fn4=trim(fn4_1)//fn4_2//fn4_3//fn4_4
		!	write(*,*)'filename=',fn4,'   myid=',myid
        else
			write(*,*)'form_dynamic_filename error'
			stop
        endif
	elseif(myid.ge.10.and.myid.le.99)then
        write(no2,'(i2)') myid
        fn4_4='00'//no2
        if(ns_start.ge.1.and.ns_start.le.9)then
			write(ns1,'(i1)')ns_start
			fn4_2='0000'//ns1
			fn4=trim(fn4_1)//fn4_2//fn4_3//fn4_4
		!	write(*,*)'filename=',fn4,'   myid=',myid
        elseif(ns_start.ge.10.and.ns_start.le.99)then
			write(ns2,'(i2)')ns_start
			fn4_2='000'//ns2
			fn4=trim(fn4_1)//fn4_2//fn4_3//fn4_4
		!	write(*,*)'filename=',fn4,'   myid=',myid
        elseif(ns_start.ge.100.and.ns_start.le.999)then
			write(ns3,'(i3)')ns_start
			fn4_2='00'//ns3
			fn4=trim(fn4_1)//fn4_2//fn4_3//fn4_4
		!	write(*,*)'filename=',fn4,'   myid=',myid
        elseif(ns_start.ge.1000.and.ns_start.le.9999)then
			write(ns4,'(i4)')ns_start
			fn4_2='0'//ns4
			fn4=trim(fn4_1)//fn4_2//fn4_3//fn4_4
		!	write(*,*)'filename=',fn4,'   myid=',myid
        elseif(ns_start.ge.10000.and.ns_start.le.99999)then
			write(ns5,'(i5)')ns_start
			fn4_2=ns5
			fn4=trim(fn4_1)//fn4_2//fn4_3//fn4_4
		!	write(*,*)'filename=',fn4,'   myid=',myid
        else
			write(*,*)'form_dynamic_filename error'
			stop
        endif
    elseif(myid.ge.100.and.myid.le.999)then
		write(no3,'(i3)') myid
		fn4_4='0'//no3
        if(ns_start.ge.1.and.ns_start.le.9)then
			write(ns1,'(i1)')ns_start
			fn4_2='0000'//ns1
			fn4=trim(fn4_1)//fn4_2//fn4_3//fn4_4
		!	write(*,*)'filename=',fn4,'   myid=',myid
        elseif(ns_start.ge.10.and.ns_start.le.99)then
			write(ns2,'(i2)')ns_start
			fn4_2='000'//ns2
			fn4=trim(fn4_1)//fn4_2//fn4_3//fn4_4
		!	write(*,*)'filename=',fn4,'   myid=',myid
        elseif(ns_start.ge.100.and.ns_start.le.999)then
			write(ns3,'(i3)')ns_start
			fn4_2='00'//ns3
			fn4=trim(fn4_1)//fn4_2//fn4_3//fn4_4
		!	write(*,*)'filename=',fn4,'   myid=',myid
        elseif(ns_start.ge.1000.and.ns_start.le.9999)then
			write(ns4,'(i4)')ns_start
			fn4_2='0'//ns4
			fn4=trim(fn4_1)//fn4_2//fn4_3//fn4_4
		!	write(*,*)'filename=',fn4,'   myid=',myid
        elseif(ns_start.ge.10000.and.ns_start.le.99999)then
			write(ns5,'(i5)')ns_start
			fn4_2=ns5
			fn4=trim(fn4_1)//fn4_2//fn4_3//fn4_4
		!	write(*,*)'filename=',fn4,'   myid=',myid
        else
			write(*,*)'form_dynamic_filename error'
			stop
        endif
	elseif(myid.ge.1000.and.myid.le.9999)then
        write(no4,'(i4)') myid
        fn4_4=no4
        if(ns_start.ge.1.and.ns_start.le.9)then
			write(ns1,'(i1)')ns_start
			fn4_2='0000'//ns1
			fn4=trim(fn4_1)//fn4_2//fn4_3//fn4_4
		!	write(*,*)'filename=',fn4,'   myid=',myid
        elseif(ns_start.ge.10.and.ns_start.le.99)then
			write(ns2,'(i2)')ns_start
			fn4_2='000'//ns2
			fn4=trim(fn4_1)//fn4_2//fn4_3//fn4_4
		!	write(*,*)'filename=',fn4,'   myid=',myid
        elseif(ns_start.ge.100.and.ns_start.le.999)then
			write(ns3,'(i3)')ns_start
			fn4_2='00'//ns3
			fn4=trim(fn4_1)//fn4_2//fn4_3//fn4_4
		!	write(*,*)'filename=',fn4,'   myid=',myid
        elseif(ns_start.ge.1000.and.ns_start.le.9999)then
			write(ns4,'(i4)')ns_start
			fn4_2='0'//ns4
			fn4=trim(fn4_1)//fn4_2//fn4_3//fn4_4
		!	write(*,*)'filename=',fn4,'   myid=',myid
        elseif(ns_start.ge.10000.and.ns_start.le.99999)then
			write(ns5,'(i5)')ns_start
			fn4_2=ns5
			fn4=trim(fn4_1)//fn4_2//fn4_3//fn4_4
		!	write(*,*)'filename=',fn4,'   myid=',myid
        else
			write(*,*)'form_dynamic_filename error'
			stop
        endif
	else
        write(*,*)'form_dynamic_filename error'
        stop
    endif

    return
    end	subroutine


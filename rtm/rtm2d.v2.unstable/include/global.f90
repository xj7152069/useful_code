!*******************************************************************!
!*        Begin The Define of module Need In This Program           ! 
!*==================================================================!
      module global                                                !
       real   ,parameter :: PI = 3.14159265359                     !
       integer,parameter :: lbyte = 1                              !
	   integer,parameter :: para_int_fnum	=	50
	   integer,parameter :: para_long_fnum	=	50
	   integer,parameter :: para_float_fnum	=	50
	   integer,parameter :: para_double_fnum=	50
	   integer,parameter :: para_char_fnum	=	50
	   integer,parameter :: para_char_flen	=	256
      end module global                                            !

!*==================================================================!
      module header_module
       type :: segy
       integer(kind = 4) :: tracl
       integer(kind = 4) :: tracr
       integer(kind = 4) :: fldr
       integer(kind = 4) :: tracf
       integer(kind = 4) :: ep
       integer(kind = 4) :: cdp
       integer(kind = 4) :: cdpt
       integer(kind = 2) :: trid
       integer(kind = 2) :: nvs
       integer(kind = 2) :: nhs
       integer(kind = 2) :: duse
       integer(kind = 4) :: offset
       integer(kind = 4) :: gelev
       integer(kind = 4) :: selev
       integer(kind = 4) :: sdepth
       integer(kind = 4) :: gdel
       integer(kind = 4) :: sedl
       integer(kind = 4) :: swdep
       integer(kind = 4) :: gwdep
       integer(kind = 2) :: scalel
       integer(kind = 2) :: scalco
       integer(kind = 4) :: sx
       integer(kind = 4) :: sy
       integer(kind = 4) :: gx
       integer(kind = 4) :: gy
       integer(kind = 2) :: counit
       integer(kind = 2) :: wevel
       integer(kind = 2) :: swevel
       integer(kind = 2) :: sut
       integer(kind = 2) :: gut
       integer(kind = 2) :: sstat
       integer(kind = 2) :: gstat
       integer(kind = 2) :: tstat
       integer(kind = 2) :: laga
       integer(kind = 2) :: lagb
       integer(kind = 2) :: delrt
       integer(kind = 2) :: muts
       integer(kind = 2) :: mute
       integer(kind = 2) :: ns
       integer(kind = 2) :: dt
       integer(kind = 2) :: gain
       integer(kind = 2) :: igc
       integer(kind = 2) :: igi
       integer(kind = 2) :: corr
       integer(kind = 2) :: sfs
       integer(kind = 2) :: sfe
       integer(kind = 2) :: slen
       integer(kind = 2) :: styp
       integer(kind = 2) :: stas
       integer(kind = 2) :: stae
       integer(kind = 2) :: tatyp
       integer(kind = 2) :: afilf
       integer(kind = 2) :: afils
       integer(kind = 2) :: nofilf
       integer(kind = 2) :: nofils
       integer(kind = 2) :: lcf
       integer(kind = 2) :: hcf
       integer(kind = 2) :: lcs
       integer(kind = 2) :: hcs
       integer(kind = 2) :: year
       integer(kind = 2) :: day
       integer(kind = 2) :: hour
       integer(kind = 2) :: minute
       integer(kind = 2) :: sec
       integer(kind = 2) :: timbas
       integer(kind = 2) :: trwf
       integer(kind = 2) :: grnors
       integer(kind = 2) :: grnofr
       integer(kind = 2) :: grnlof
       integer(kind = 2) :: gaps
       integer(kind = 2) :: otrav
       real(kind = 4) :: d1
       real(kind = 4) :: f1
       real(kind = 4) :: d2
       real(kind = 4) :: f2
       real(kind = 4) :: ungpow
       real(kind = 4) :: unscale
       integer(kind = 4) :: ntr
       integer(kind = 2) :: mark
       integer(kind = 2) :: shortpad
       integer(kind = 2) :: unass(14)
      end type
      end module

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


!*        END Of The Define of module Need In This Program          ! 
!*******************************************************************!


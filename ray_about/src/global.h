!******************************************************************!
!        Begin The Define of module Need In This Program           ! 
!==================================================================!
      module global                                                !
       real   ,parameter :: PI = 3.14159265359                     !
       real   ,parameter :: Boundx = -0.5                           !
       real   ,parameter :: Threshold = 1.0/25
       integer,parameter :: LMAX = 990000000
       integer,parameter :: Mul  = 10                 ! ds = dvz/mul
       integer,parameter :: LTABLE = 4
       integer,parameter :: NTABLE = 101
       integer,parameter :: CELLSIZE = 1                           !
      end module global                                            !

      module types                                                 !
!====== One Cell in which to interpolate complex time and amplitude!
      type :: cell                                                 !
       real    live       !a flag which denote a live cell(1:live) !
       real    time       !traveltime                              !
       real    amp        !a parameter of amplitude                !
       real    road       !a parameter of amplitude                !
      end type cell                                                !
!====== One step along ray                                         !
      type :: raystep                                              !
       real    t           !time                                   !
       real    a           !angle                                  !
       real    x,z         !coordinate x,z                         !
       real    thegam      !a parameter of wide of the beam        !
       real    c,s         !cos(angle) and sin(angle)              !
       real    v,dvdx,dvdz !velocity and its derivatives           !
      end type raystep                                             !
      end module types                                             !
!==================================================================!
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
      type(segy) :: head
      end module

!        END Of The Define of module Need In This Program          ! 
!******************************************************************!


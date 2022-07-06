!==================================================================*
!                                                                  *
!           3D OFFSET PLANE-WAVE PRESTACK TIME MIGRATION           *
!           WITH FINITE DIFFERENCE METHOD PLUS F-K DMOAIN          *
!           ERROR COMPENSATION FOR 3D ACOUSTIC MEDIA.              *
!==================================================================*
!              ALL   RIGHTS   RESERVED                             *
!==================================================================*
!       WPI, School of Ocean and Earth Science, TongJi University  *
!==================================================================*
!       Remarks:                                                   *
!         Author: Feng Bo.                                         *
!         Date  : 2008.12 -- 2022.02                               *
!         Current Version: V20                                     *
!==================================================================*
!       History:                                                   *
!         Stable Version: V16.                                     *
!         Date  : 2010.07                                          *
! Name:     3d_Phx_planewave_acoustic_pstm_for_iCluster_v16.f      *
! Package:  pw3dlib.update.100720.tar                              *
!         Stable Version: V14.                                     *
!         Date  : 2009.08                                          *
! Name:     3d_Phx_planewave_acoustic_pstm_for_iCluster_v14.f      *
! Package:  pw3dlib.update.090820.tar                              *
!==================================================================*
!         Logs.                                                    *
!==================================================================*
!         Updated for V20, 2022-02-1                               *
!           Modified all I-O functions.                            *
!           Simplified algorithm.                                  *
!         Updated for V18, 2010-07-29                              *
!           adding Tau-Axis Break-Point protection feature.        *
!         Updated for V17, 2010-07-24                              *
!           Adding Visco-Acoustic Migration flag.                  *
!           ViscoFlag: =0, acoustic mig; =1, visco-acoustic mig.   *
!         Updated for V16, 2010-07-17                              *
!           Modified Error&Anomaly handling.                       *
!           Adding Open file status checking.                      *
!           Adding Break-Point protection feature.                 *
!           Using Multi-Threads FFT-2D. FFT2D_OMP().               *
!         Updated for V15, 2010-07-15                              *
!           Modified regular-file writting, from data() to file(). *
!           Using Multi-threads FFT-1D. **Profile2d().             *
!==================================================================*
!       Summary:                                                   *
!         This program do a 3D offset plane-wave prestack          *
!       time migration with finite-difference method               *
!       plus W-K domain error compensation for acoustic media.     *
!==================================================================*

PROGRAM PRESTACK_3D_PLANE_WAVE_FD_TIME_MIGRATION
    implicit none
    INCLUDE 'mpif.h'
    integer :: LByte
    REAL    :: PAI2

    !     Parameters reading from parameter file.
    REAL     F1, F2, F3, F4
    INTEGER  Ntau, Iphr1, Iphr2, IBREAK, ViscoFlag
    INTEGER  NtauExtraP, ITauBreak, Ktaur, Ktauw
    INTEGER  IphrR1, IphrR2
    CHARACTER*1024 FN_PAR, FN1, FN2, FN3, FN4, FNTMPW, FNTMPR

    !Parameters for 3D offset-planewave decompostion data.
    INTEGER  IPPAR(25), IFN01
    REAL     FPPAR(20)
    INTEGER  LT, Nmx, Nmy, Nphr
    INTEGER  Ncdp_first, Ncdp_final, Ncdp_step
    INTEGER  Nline_first, Nline_final, Nline_step
    REAL     DT, Dmx, Dmy, Phrmin, Dphr

    !Parameters for 3D interval velocity field(time domain).
    INTEGER  IVPAR(25), IFN02
    REAL     FVPAR(10)
    INTEGER  Nvtau, Nvmy, Nvmx
    INTEGER  Nvtau_first, Nvtau_final, Nvtau_step
    INTEGER  Nvline_first, Nvline_final, Nvline_step
    INTEGER  Nvcdp_first, Nvcdp_final, Nvcdp_step
    REAL     Dvtau, Dvx, Dvy

    !Parameters for common imaging gather of 3D offset plane-wave PSTM.
    INTEGER  icpar(25), IFN03, iflag3
    REAL     cpar(20)
    REAL     Dtau

    INTEGER  Nmx0, Nmx1, Nmy0, Nmy1
    INTEGER  NW1, NW2, NW3, NW4, NW
    REAL     DW, W1, W2, W3, W4

    REAL,ALLOCATABLE :: BUF(:)
    INTEGER*8  L1,  L2,  L3,  L4,  L5,  L6,  L7,  L8,  L9,  L10
    INTEGER*8  L11, L12, L13, L14, L15, L16, L17, L18, L19, L20

    CHARACTER*1024 ProcessorName
    INTEGER  ierr, myid, nproc, NAME_LENGTH

    INTEGER  ierror_flag

    CALL MPI_INIT(ierr)
    CALL MPI_COMM_RANK(MPI_COMM_WORLD, myid,  ierr)
    CALL MPI_COMM_SIZE(MPI_COMM_WORLD, nproc, ierr)

    if( 0 .eq. myid ) then
        CALL GETARG(1, FN_PAR)
    endif
    CALL MPI_BCAST(FN_PAR, 1024, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
    Call READPAR(FN_PAR, FN1, FN2, FN3, FN4, &
        Nvtau, Iphr1, Iphr2, F1, F2, F3, F4, IBREAK, ViscoFlag, &
        NtauExtraP, ITauBreak, Ktaur, Ktauw, FNTMPW, FNTMPR, &
        IphrR1, IphrR2)

    !Init. of parameters.
    Lbyte = 1
    PAI2  = 2*3.14159265359
    Ntau  = NtauExtraP

    !FN1: Offset Plane-Wave data file.
    !CALL GDBFCGFENTERFILENAME(FN1)
    !FN2: Interval velocity field(time domain, slice data).
    !CALL GDBFCGFENTERFILENAME(FN2)
    !FN3: Common-Ph imaging volume, Ntau*Nmx*Nmy*Nphr.
    !CALL GDBFCGFENTERFILENAME(FN3)
    IFN01 = 11
    IFN02 = 12
    IFN03 = 13

    !For Luojia3D data:
    !cdp=[2401,3700], ncdp=1300, dcdp=6.25m.
    !line=[1991,2460], nline=470, dline=12.5m
     LT          = 2901
     DT          = 2.0         !ms
     Nmx         = 1300
     Ncdp_first  = 2401
     Ncdp_final  = 3700
     Ncdp_step   = 1
     Dmx         = 6.25        !m
     Nmy         = 470
     Nline_first = 1991
     Nline_final = 2460
     Nline_step  = 1
     Dmy         = 12.5        !m
     Nphr        = 81
     Phrmin      = 0
     Dphr        = 10.0/1.0E6   !s/m

    !For synthetic Data:
!   LT          = 4000
!   DT          = 1.0         !ms
!   Nmx         = 500
!   Ncdp_first  = 151
!   Ncdp_final  = 650
!   Ncdp_step   = 1
!   Dmx         = 20.0
!   Nmy         = 500
!   Nline_first = 151
!   Nline_final = 650
!   Nline_step  = 1
!   Dmy         = 20.
!   Nphr        = 61
!   Phrmin      = 0
!   Dphr        = 10.0/1.0E6   !s/m

    !New requirement will be included recently.
    IF (Iphr1 .eq. 0 .and. Iphr2 .eq. 0) THEN
        Iphr1 = 1
        Iphr2 = Nphr
    END IF
    IF (Iphr1 .gt. Nphr .or. Iphr2 .gt. Nphr .or. &
        Iphr1 .lt. 1 .or. Iphr2 .lt. 1 .or. Iphr1 .gt. Iphr2) THEN
        write(*,*) 'Wrong index of Ray-parameter.'
        write(*,*) ' Iphr1 = ', Iphr1
        write(*,*) ' Iphr2 = ', Iphr2
        write(*,*) ' Nphr  = ', Nphr
        GOTO 9999
    END IF
    !Update the number of ray-parameters: Nphr
    Nphr = Iphr2 - Iphr1 + 1

    !Nvmy          = 3361
    !Nvline_first  = 35
    !Nvline_final  = 3395
    !Nvline_step   = 1
    !Nvmx          = 4641
    !Nvcdp_first   = 493
    !Nvcdp_final   = 5133
    !Nvcdp_step    = 1

    !For Luojia3D data:
    !cdp=[2401,3700], ncdp=1300, dcdp=6.25m.
    !line=[1991,2460], nline=470, dline=6.25m.
!    Nvtau         = 2800
     Nvtau_first   = 1
     Nvtau_final   = Nvtau
     Nvtau_step    = 1
     Nvmy          = 800
     Nvline_first  = 1801
     Nvline_final  = 2600
     Nvline_step   = 1
     Nvmx          = 2000
     Nvcdp_first   = 2001
     Nvcdp_final   = 4000
     Nvcdp_step    = 1
     Dvtau         = 2.0       !ms
     Dvy           = 12.5      !m
     Dvx           = 6.25      !m

    !For synthetic Data:
!   Nvtau         = 4000
!   Nvtau_first   = 1
!   Nvtau_final   = Nvtau
!   Nvtau_step    = 1
!   Nvmy          = 800
!   Nvline_first  = 1
!   Nvline_final  = 800
!   Nvline_step   = 1
!   Nvmx          = 800
!   Nvcdp_first   = 1
!   Nvcdp_final   = 800
!   Nvcdp_step    = 1
!   Dvtau         = 1.0       !ms
!   Dvy           = 20.0
!   Dvx           = 20.0

    open(IFN02, file=trim(adjustl(FN2)), access='direct', recl=Nvmx*Nvmy)

    !Check the 3d velocity field.
    IF ((Nvline_first .gt. Nline_first) .or. &
        (Nvline_final .lt. Nline_final) .or. &
        (Nvcdp_first  .gt. Ncdp_first) .or. &
        (Nvcdp_final  .lt. Ncdp_final)) THEN
        write(*,*) 'Velocity Field Do Not Match Imaging Area.'
        write(*,*) ' Nline_first = ', Nline_first
        write(*,*) ' Nvline_first = ', Nvline_first
        write(*,*) ' Nline_final = ', Nline_final
        write(*,*) ' Nvline_final = ', Nvline_final
        write(*,*) ' Ncdp_first = ', Ncdp_first
        write(*,*) ' Nvcdp_first = ', Nvcdp_first
        write(*,*) ' Ncdp_final = ', Ncdp_final
        write(*,*) ' Nvcdp_final = ', Nvcdp_final
        write(*,*) ' Program Exit. '
        GOTO 9999
    END IF
    IF (Dmx .le. 1.0E-2 .or. Dvx .le. 1.0E-2 .or. &
        Dmy .le. 1.0E-2 .or. Dvy .le. 1.0E-2) THEN
        write(*,*) ' Velocity field sampling error. '
        write(*,*) ' Or Plane-wave data sampling error. '
        write(*,*) ' Program Exit. '
        GOTO 9999
    END IF
    IF (Dmx/Dvx .lt. 0.95 .or. Dmx/Dvx .gt. 1.05 .or. &
        Dmy/Dvy .lt. 0.95 .or. Dmy/Dvy .gt. 1.05) THEN
        write(*,*) ' Velocity Field Grid Error!'
        write(*,*) ' Dmx = ', Dmx
        write(*,*) ' Dvx = ', Dvx
        write(*,*) ' Dmy = ', Dmy
        write(*,*) ' Dvy = ', Dvy
        write(*,*) ' Program Exit. '
        GOTO 9999
    END IF

    IF ((Nvtau .lt. Ntau)) THEN
        write(*,*) ' Extrapolation Number Error!'
        write(*,*) ' Nvtau = ', Nvtau
        write(*,*) ' Ntau  = ', Ntau
        Ntau = Nvtau
        write(*,*) ' Ntau is set to Nvtau automaticlly'
        write(*,*) ' Ntau  = ', Ntau
    END IF

    DW      = PAI2*1000.0/(LT*DT)
    W1      = PAI2*F1
    NW1     = W1/DW + 0.5
    W2      = PAI2*F2
    NW2     = W2/DW + 0.5
    W3      = PAI2*F3
    NW3     = W3/DW + 0.5
    W4      = PAI2*F4
    NW4     = W4/DW + 0.5
    NW      = NW4-NW1 + 1

    IF (abs(W1) .lt. 1.0E-5) THEN
        IF (myid .eq. 0) THEN
            write(*,*) ' Parameter: F1 Error. '
            write(*,*) ' Make sure that F1 > 1.0Hz  '
            write(*,*) ' Program Exit. '
        END IF
        GOTO 9999
    END IF

    Dtau = Dvtau*0.001     !unit is second.
    Nmx0 = 100
    Nmx1 = Nmx
    Nmx  = Nmx1 + 2*Nmx0
    Nmy0 = 100
    Nmy1 = Nmy
    Nmy  = Nmy1 + 2*Nmy0

    !Print Parameters on the screen.
    IF (myid .eq. 0) THEN
        write(*,*) '*************************************************'
        write(*,*) '                                                 '
        write(*,*) '*****  Parameters of 5d Planewave Data.  ********'
        write(*,*) ' Nline_first  = ', Nline_first
        write(*,*) ' Nline_final  = ', Nline_final
        write(*,*) ' Nline_step   = ', Nline_step
        write(*,*) ' Nmy          = ', Nmy
        write(*,*) ' Ncdp_first   = ', Ncdp_first
        write(*,*) ' Ncdp_final   = ', Ncdp_final
        write(*,*) ' Ncdp_step    = ', Ncdp_step
        write(*,*) ' Nmx          = ', Nmx
        write(*,*) ' LT           = ', LT
        write(*,*) ' DT(ms)       = ', DT
        write(*,*) ' Dmx(m)       = ', Dmx
        write(*,*) ' Dmy(m)       = ', Dmy
        write(*,*) ' Nphr         = ', Nphr
        write(*,*) ' Dphr         = ', Dphr
        write(*,*) ' Phrmin       = ', Phrmin
        write(*,*) '                                                 '

        write(*,*) '*************************************************'
        write(*,*) '                                                 '
        write(*,*) '*****  Parameters of 3d Velocity Field.  ********'
        write(*,*) ' Nvtau_first   = ', Nvtau_first
        write(*,*) ' Nvtau_final   = ', Nvtau_final
        write(*,*) ' Nvtau_step    = ', Nvtau_step
        write(*,*) ' Nvtau         = ', Nvtau
        write(*,*) ' Nvline_first  = ', Nvline_first
        write(*,*) ' Nvline_final  = ', Nvline_final
        write(*,*) ' Nvline_step   = ', Nvline_step
        write(*,*) ' Nvmy          = ', Nvmy
        write(*,*) ' Nvcdp_first   = ', Nvcdp_first
        write(*,*) ' Nvcdp_final   = ', Nvcdp_final
        write(*,*) ' Nvcdp_step    = ', Nvcdp_step
        write(*,*) ' Nvmx          = ', Nvmx
        write(*,*) ' Dvx(m)        = ', Dvx
        write(*,*) ' Dvy(m)        = ', Dvy
        write(*,*) ' Dvtau(ms)     = ', Dvtau
        write(*,*) '                                                 '

        write(*,*) '*************************************************'
        write(*,*) '                                                 '
        write(*,*) '*****  Parameters of 3d CIGs  ********'
        write(*,*) ' Dtau(s) = ', Dtau
        write(*,*) ' Ntau    = ', Ntau
        write(*,*) ' Iphr1   = ', Iphr1
        write(*,*) ' Iphr2   = ', Iphr2
        write(*,*) ' Nmx     = ', Nmx
        write(*,*) ' Nmy     = ', Nmy
        write(*,*) '                                                 '

        write(*,*) '*************************************************'
        write(*,*) '                                                 '
        write(*,*) ' F1, W1, NW1 = ', F1, W1, NW1
        write(*,*) ' F2, W2, NW2 = ', F2, W2, NW2
        write(*,*) ' F3, W3, NW3 = ', F3, W3, NW3
        write(*,*) ' F4, W4, NW4 = ', F4, W4, NW4
        write(*,*) ' NW = ', NW, ' processor number:', nproc
        write(*,*) '                                                 '
        write(*,*) '*************************************************'
    END IF

    L1 = 1
    L2 = L1  + Nmx*Nmy         !VVtau
    L3 = L2  + LT              !Trace
    L4 = L3  + 2*LT            !CTrace
    L5 = L4  + 2*Nmx           !AA_X
    L6 = L5  + 2*Nmx           !BB_X
    L7 = L6  + 2*Nmx           !CC_X
    L8 = L7  + 2*Nmx           !DD_X
    L9 = L8  + 2*Nmx           !XX_X
    L10= L9  + 2*Nmx           !YY_X
    L11= L10 + 2*Nmy           !AA_Y
    L12= L11 + 2*Nmy           !BB_Y
    L13= L12 + 2*Nmy           !CC_Y
    L14= L13 + 2*Nmy           !DD_Y
    L15= L14 + 2*Nmy           !XX_Y
    L16= L15 + 2*Nmy           !YY_Y
    L17= L16 + 2*Nmx*Nmy       !DATW_W
    L18= L17 + 2*Nmx*Nmy       !XX
!   L19= L18 + Ntau*Nmx1*Nmy1  !DATW_R

    ierror_flag = 0
    ALLOCATE(BUF(L18), stat=ierror_flag)
    IF (ierror_flag .ne. 0) THEN
        CALL MPI_GET_PROCESSOR_NAME(ProcessorName, NAME_LENGTH, IERR)
        write(*,*) ' Not Enough Memory at Node: ', TRIM(ProcessorName), '  myid is: ', myid
        write(*,*) ' Program Exit. '
        GOTO 9999
    END IF

    IF (myid .eq. 0) THEN
        write(*,*) '                                                 '
        write(*,*) ' Total Memory Needed:', (L18-1)*4.0/1024/1024,'MB'
        write(*,*) '                                                 '
        write(*,*) '*************************************************'
        write(*,*) '*                                               *'
        write(*,*) '*  3D offset plane-wave FD PSTM begin.          *'
        write(*,*) '*                                               *'
        write(*,*) '*************************************************'
        write(*,*) '                                                 '
    END IF

    IF (nproc .eq. 1) THEN
        write(*,*) '*************************************************'
        write(*,*) '*                                               *'
        write(*,*) '*  Require processor number larger than 1.      *'
        write(*,*) '*  Program Exit Now!                            *'
        write(*,*) '*                                               *'
        write(*,*) '*************************************************'
        GOTO 9999
    END IF

    CALL OFFSET_PLANEWAVE_FD_PSTM_3D(&
        BUF(L1), BUF(L2), BUF(L3),&
        BUF(L4),  BUF(L5),  BUF(L6),  BUF(L7),  BUF(L8),  BUF(L9),&
        BUF(L10), BUF(L11), BUF(L12), BUF(L13), BUF(L14), BUF(L15),&
        BUF(L16), BUF(L17), &
        Nvline_first, Nvline_final, Nvline_step, Nvmy, Dvx,&
        Nvcdp_first,  Nvcdp_final,  Nvcdp_step,  Nvmx, Dvy,&
        Nvtau, &
        Nline_first, Nline_final, Nline_step, Nmy, Nmy0, Nmy1, Dmy,&
        Ncdp_first,  Ncdp_final,  Ncdp_step,  Nmx, Nmx0, Nmx1, Dmx,&
        Nphr, Iphr1, Iphr2, Phrmin, Dphr,&
        NW1, NW2, NW3, NW4, NW, DW, LT, Ntau, Dtau,&
        IFN01, IFN02, IFN03, &
        FN1, FN2, FN3, FN4, FNTMPW, FNTMPR,&
        ViscoFlag, NtauExtraP, IphrR1, IphrR2, &
        ITauBreak, Ktaur, Ktauw, Lbyte, myid, nproc)

    DEALLOCATE(BUF)

    9999  CONTINUE
    IF (myid .eq. 0) THEN
        write(*,*) '                                                 '
        write(*,*) '*************************************************'
        write(*,*) '*                                               *'
        write(*,*) '*  3D offset plane-wave FD PSTM done.           *'
        write(*,*) '*                                               *'
        write(*,*) '*************************************************'
        write(*,*) '                                                 '
    END IF

    CALL MPI_Finalize(ierr)

    STOP
END PROGRAM PRESTACK_3D_PLANE_WAVE_FD_TIME_MIGRATION

!======================================================================*

SUBROUTINE OFFSET_PLANEWAVE_FD_PSTM_3D(&
        VVtau, Trace, CTrace,&
        AA_X, BB_X, CC_X, DD_X, XX_X, YY_X,&
        AA_Y, BB_Y, CC_Y, DD_Y, XX_Y, YY_Y,&
        DATW_W, XX, &
        Nvline_first, Nvline_final, Nvline_step, Nvmy, Dvx,&
        Nvcdp_first,  Nvcdp_final,  Nvcdp_step,  Nvmx, Dvy,&
        Nvtau, &
        Nline_first, Nline_final, Nline_step, Nmy, Nmy0, Nmy1, Dmy,&
        Ncdp_first,  Ncdp_final,  Ncdp_step,  Nmx, Nmx0, Nmx1, Dmx,&
        Nphr, Iphr1, Iphr2, Phrmin, Dphr,&
        NW1, NW2, NW3, NW4, NW, DW, LT, Ntau, Dtau,&
        IFN01, IFN02, IFN03, &
        FN1, FN2, FN3, FN4, FNTMPW, FNTMPR,&
        ViscoFlag, NtauExtraP, IphrR1, IphrR2,&
        ITauBreak, Ktaur, Ktauw, Lbyte, myid, nproc)

    implicit none
    INCLUDE 'mpif.h'
    integer :: LByte
    REAL    :: PAI, PAI2
    INTEGER :: IFN01, IFN02, IFN03

    INTEGER  status(MPI_STATUS_SIZE)
    INTEGER  ierr, myid, nproc
    INTEGER  ViscoFlag, NtauExtraP, ITauBreak, Ktaur, Ktauw

    INTEGER  Nvtau, Nvmx, Nvmy
    INTEGER  Nvline_first, Nvline_final, Nvline_step
    INTEGER  Nvcdp_first,  Nvcdp_final,  Nvcdp_step
    REAL     Dvtau, Dvx, Dvy

    INTEGER  Nline_first, Nline_final, Nline_step, Nmy, Nmy0, Nmy1
    INTEGER  Ncdp_first,  Ncdp_final,  Ncdp_step,  Nmx, Nmx0, Nmx1
    INTEGER  Ntau, LT, Nphr, Iphr1, Iphr2
    INTEGER  NW1, NW2, NW3, NW4, NW, KmyNUM, KmxNUM
    INTEGER  IphrR1, IphrR2
    REAL     Dtau, DIP_max, DIP_min, Dmx, Dmy, Phrmin, Dphr
    REAL     DW, W1, W2, W3, W4, Dky, Dkx, W
    REAL     ALFA, Rkz_max, RPh

    REAL     Trace(LT), VVtau(Nmx, Nmy)
    REAL     ACOE(5), BCOE(5)

    COMPLEX  AA_X(Nmx), BB_X(Nmx), CC_X(Nmx), DD_X(Nmx)
    COMPLEX  XX_X(Nmx), YY_X(Nmx)
    COMPLEX  AA_Y(Nmy), BB_Y(Nmy), CC_Y(Nmy), DD_Y(Nmy)
    COMPLEX  XX_Y(Nmy), YY_Y(Nmy)
    COMPLEX  CTrace(LT), DATW_W(Nmx, Nmy), XX(Nmx, Nmy)

    ! internal parameters.
    INTEGER  ierror1, ierror2, ierror3, ierror4, ierror5, ierror6
    COMPLEX, ALLOCATABLE :: SENDBUF(:)
    COMPLEX, ALLOCATABLE :: TMP2D(:,:)
    COMPLEX, ALLOCATABLE :: DATW(:,:,:)
    REAL,    ALLOCATABLE :: RIMAG(:,:,:)
    REAL,    ALLOCATABLE :: DATW_R(:,:,:)

    INTEGER  LOOPP, Ncount2d
    INTEGER  NumSend, sender, itag, SlaveNodeIndex

    !parameters for writing temp files.
    CHARACTER*1024 FN1, FN2, FN3
    CHARACTER*1024 FN4, FNTMPW, FNTMPR, FNW, FNR, FNCIG
    CHARACTER*4   ciphr
    INTEGER*8     LRIMAG, LDATW

    !2022.02.09, debug only.
    REAL,    ALLOCATABLE :: vel3d(:,:,:)
    REAL,    ALLOCATABLE :: vel3dall(:,:,:)
    INTEGER  imx, ivx, ix, IDX
    INTEGER  imy, ivy, iy
    INTEGER  it, itau, ivtau, Iphr
    INTEGER  Iww, IW
    REAL     v_min, v_max
    REAL*8   stime, etime

    PAI     = 3.14159265359
    PAI2    = 2*3.14159265359
    DIP_max = 79.0
    DIP_min = 25.0
    ALFA    = 0.113
    Dky     = PAI2/(Nmy*Dmy)
    Dkx     = PAI2/(Nmx*Dmx)
    Rkz_max = cos(PAI2*DIP_max/360.0)

    CALL make_order(DIP_max, Acoe, Bcoe, LOOPP)

    ierror1 = 0
    ierror2 = 0
    ierror3 = 0
    ierror4 = 0
    ierror5 = 0
    ALLOCATE(SENDBUF(Nmx*Nmy+2),  stat=ierror1)
    ALLOCATE(TMP2D(Nmy,Nmx),      stat=ierror2)
    ALLOCATE(vel3d(Nmx,Nmy,Ntau), stat=ierror3)
    ALLOCATE(DATW_R(Ntau,Nmx1,Nmy1),    stat=ierror4)
    ALLOCATE(RIMAG(Ntau,Nmx1,Nmy1),     stat=ierror5)
    ALLOCATE(vel3dall(Nvmx,Nvmy,Nvtau), stat=ierror6)
    IF ((ierror1 .ne. 0) .or. (ierror2 .ne. 0) &
        .or. (ierror3 .ne. 0) .or. (ierror4 .ne. 0) &
        .or. (ierror5 .ne. 0) .or. (ierror6 .ne. 0)) THEN
        write(*,*) ' Not Enough Memory for buffers.'
        write(*,*) ' 3D Program STOP Now!'
        STOP
    END IF

    IF (myid .eq. 0) THEN
        write(*,*) '********  Parameters of 5d Planewave Data.  ********'
        write(*,*) ' Nmy          = ', Nmy
        write(*,*) ' Nmx          = ', Nmx
        write(*,*) ' LT           = ', LT
        write(*,*) ' Dmx(m)       = ', Dmx
        write(*,*) ' Dmy(m)       = ', Dmy
        write(*,*) ' Nphr         = ', Nphr
        write(*,*) ' Dphr         = ', Dphr
        write(*,*) ' Phrmin       = ', Phrmin

        write(*,*) '********  Parameters of 3d CIGs  ********'
        write(*,*) ' Dtau(s) = ', Dtau
        write(*,*) ' Ntau    = ', Ntau
        write(*,*) ' Nvtau   = ', Nvtau

        write(*,*) '********  Parameters of Extrapolation  ********'
        write(*,*) ' ALFA   = ', ALFA
        write(*,*) ' Dky    = ', Dky
        write(*,*) ' Dkx    = ', Dkx
        write(*,*) ' DIP_max    = ', DIP_max
        write(*,*) ' DIP_min    = ', DIP_min
        write(*,*) ' Rkz_max    = ', Rkz_max
    END IF

    !Each node read the velocity model into memory.
    DO Itau = 1, Nvtau
!      write(*,*) 'Itau=',Itau, 'Nvtau=', Nvtau
       read(IFN02, rec=Itau) ((vel3dall(ix,iy,Itau),ix=1,Nvmx),iy=1,Nvmy)
    END DO

    !Init. of buffers
    vel3d = 0
    VVtau = 0
    DO Itau = 1, Ntau
       DO Imy = 1, Nmy
          DO Imx = 1, Nmx
             Ivx = Imx + (Ncdp_first - Nvcdp_first - Nmx0)
             Ivy = Imy + (Nline_first - Nvline_first - Nmy0)
             vel3d(Imx, Imy, Itau) = vel3dall(Ivx, Ivy, Itau)
          END DO
       END DO
    END DO
    DEALLOCATE(vel3dall)

!   IF (myid .eq. 0) THEN
!       open(21, file='velo_check.dat', access='direct', recl=Nmx*Nmy)
!       DO Itau = 1, 500
!           DO Imy = 1, Nmy
!           DO Imx = 1, Nmx
!              VVtau(Imx, Imy) = vel3d(Imx, Imy, Itau)
!           END DO
!           END DO
!           write(*,*) 'Itau=',Itau, 'Ntau=', Ntau
!           WRITE(21,REC=Itau)((VVtau(IX,IY),IX=1,Nmx),IY=1,Nmy)
!       END DO
!       close(21)
!   END IF

    IF (myid .eq. 0) THEN
        !Master Node.
        ierror3 = 0
        ALLOCATE(DATW(Nmx,Nmy,NW), stat=ierror3)
        IF (ierror3 .ne. 0) THEN
            write(*,*) ' Not Enough Memory on Master Node!'
            write(*,*) ' Can not Allocate Mem for DATW!'
            write(*,*) ' Program Exit Now!'
            STOP
        END IF
        LRIMAG = 1
        LDATW  = 1
        LRIMAG = LRIMAG*Ntau*Nmx1*Nmy1
        LDATW  = LDATW*NW*Nmx*Nmy
        write(*,*) 'sizeof(RIMAG) = ', LRIMAG*4.0/1024/1024 , 'MB'
        write(*,*) 'sizeof(DATW) = ',  LDATW*8.0/1024/1024 ,  'MB'
    END IF

    DO Iphr = Iphr1, Iphr2

        stime   = MPI_Wtime();
        IF (myid .eq. 0) THEN

            DATW_R  = 0.0
            RIMAG   = 0.0
            RPh     = Phrmin + Dphr*(Iphr-1)
            write(*,*) ' Iphr = ',Iphr,' RPh= ',RPh*1.0E6,'vs/m'

            write(*,*) ' ****READ PlaneWave Data and FFT With OpenMP.'
            CALL READ_PLANEWAVE_DATA_AND_FFT_Profile2d(&
                FN1, Trace, DATW, CTrace,&
                Nmx, Nmy, NW1, NW2, NW3, NW4, NW, LT, DW,&
                Ncdp_first,  Ncdp_final,  Ncdp_step,&
                Nline_first, Nline_final, Nline_step,&
                Nphr, Iphr, IFN01, Nmx0, Nmx1, Nmy0, Nmy1)

            NumSend = 1
            DO Iww = 1, MIN(nproc-1, NW)

                write(*,*) 'Sending: Iww =', Iww, ' Iphr = ', Iphr
                IW  = Iww + NW1 - 1
                W   = (IW-1)*DW

                DO IY=1, Nmy
                    DO IX=1, Nmx
                        IDX = (IY-1)*Nmx + IX
                        SENDBUF(IDX) = DATW(IX, IY, Iww)
                    END DO
                END DO
                SENDBUF(Nmx*Nmy+1) = CMPLX(W,RPh)
                SENDBUF(Nmx*Nmy+2) = CMPLX(Iww,0)

                CALL MPI_Send(SENDBUF, Nmy*Nmx+2, MPI_COMPLEX, &
                    Iww, Iww, MPI_COMM_WORLD, IERR)
                NumSend = NumSend + 1

            END DO

            DO IW = 1, NW

                !Receive the imageing volume from slave node.
                CALL MPI_RECV(SlaveNodeIndex, 1, MPI_INTEGER, &
                    MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)

                sender = status(MPI_SOURCE)
                itag   = status(MPI_TAG)

                IF (NumSend .le. NW) THEN
                    DO IY=1, Nmy
                        DO IX=1, Nmx
                            IDX = (IY-1)*Nmx + IX
                            SENDBUF(IDX) = DATW(IX, IY, NumSend)
                        END DO
                    END DO

                    Iww = NumSend
                    W   = (NumSend-1)*DW
                    SENDBUF(Nmx*Nmy+1) = CMPLX(W,RPh)
                    SENDBUF(Nmx*Nmy+2) = CMPLX(Iww,0)

                    !write(*,*) 'Sending: myid=', myid, 'Iww =', NumSend, ' Iphr = ', Iphr
                    CALL MPI_Send(SENDBUF, Nmy*Nmx+2, MPI_COMPLEX, &
                        sender, NumSend, MPI_COMM_WORLD, IERR)
                    NumSend = NumSend + 1
                ELSE
                    CALL MPI_Send(0, 0, MPI_COMPLEX, sender, 0, MPI_COMM_WORLD, IERR)
                END IF

            END DO

        ELSE !working nodes.
            IF (myid .gt. NW) GOTO 6666  ! Nodes number larger than tasks number.

            DO WHILE(1)

                CALL MPI_RECV(SENDBUF, Nmy*Nmx+2, MPI_COMPLEX, 0, &
                    MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)

                itag   = status(MPI_TAG)
                IF (itag .EQ. 0) GOTO 5555

                DO IY=1, Nmy
                    DO IX=1, Nmx
                        IDX = (IY-1)*Nmx + IX
                        DATW_W(IX, IY) = SENDBUF(IDX)
                    END DO
                END DO
                W   =  REAL(SENDBUF(Nmx*Nmy+1))
                RPh = AIMAG(SENDBUF(Nmx*Nmy+1))
                Iww = int(REAL(SENDBUF(Nmx*Nmy+2)) + 0.5)
                write(*,*) 'Received: myid=', myid, 'Iww =', Iww

                IF ( ViscoFlag .EQ. 0 ) THEN
                    !3D Acoustic Plane-Wave Migration.
                    CALL ACOUSTIC_EXTRAPOLATION_SSF(&
                        myid, Iphr, RPh, W, DW, Iww,&
                        vel3d, NtauExtraP,&
                        DATW_W, DATW_R, VVtau,&
                        Nmx, Dmx, Nmy, Dmy, Ntau, Dtau,&
                        Dkx, Dky, Rkz_max,&
                        Nvline_first, Nvline_final, Nvmy,&
                        Nvcdp_first,  Nvcdp_final,  Nvmx,&
                        Nline_first, Nline_final, Nmy0, Nmy1,&
                        Ncdp_first, Ncdp_final, Nmx0, Nmx1)
                ELSE
                    !3D Visco-Acoustic Plane-Wave Migration.
                    CALL WXFD_EXTRAPOLATION_WXYCOMP_MPI_VISCO(&
                        DATW_W, DATW_R, XX, VVtau, TMP2D,&
                        AA_X, BB_X, CC_X, DD_X, XX_X, YY_X,&
                        AA_Y, BB_Y, CC_Y, DD_Y, XX_Y, YY_Y,&
                        ACOE, BCOE, LOOPP, ALFA,&
                        Nmx, Dmx, Nmy, Dmy, Ntau, Dtau, W, RPh,&
                        Dkx, Dky, Rkz_max, DIP_max, DIP_min,&
                        Nvline_first, Nvline_final, Nvmy,&
                        Nvcdp_first,  Nvcdp_final,  Nvmx,&
                        Nline_first, Nline_final,&
                        Ncdp_first, Ncdp_final, IFN02,&
                        Nmx0, Nmx1, Nmy0, Nmy1, Iphr, Iww, NtauExtraP,&
                        ITauBreak, Ktaur, Ktauw, FNW, FNR,&
                        IphrR1, IphrR2, myid)
                END IF

                CALL MPI_Send(SlaveNodeIndex, 1, MPI_INTEGER, 0, &
                    itag, MPI_COMM_WORLD, IERR)

            END DO

            5555 CONTINUE

        END IF ! end of myid=0

        6666 CONTINUE
        write(*,*) '--------Begin MPI_BARRIER--------, myid is', myid
        CALL MPI_BARRIER(MPI_COMM_WORLD, IERR)

        Ncount2d=Ntau*Nmx1
        DO imy = 1, nmy1
            CALL MPI_REDUCE(DATW_R(1,1,imy), RIMAG(1,1,imy), Ncount2d, &
                MPI_REAL, MPI_SUM, 0, MPI_COMM_WORLD, IERR)
        END DO

        etime   = MPI_Wtime();
        !write the common-ph migrated section.
        IF (myid .eq. 0) THEN
            write(ciphr,'(I4)') Iphr
            FNCIG=trim(adjustl(FN4))//trim(adjustl(ciphr))//'.dat'

            write(*,*) 'Iphr=', Iphr, 'computation time(s) is:', etime-stime
            write(*,*) 'Now, wrting file:', FNCIG
            write(*,*) 'Ntau=', Ntau, 'Nmx1=', Nmx1, 'Nmy1=', Nmy1

            open(IFN03, file=FNCIG, access='direct', recl=Ntau*Nmx1)
            DO imy = 1, Nmy1
                write(IFN03, rec=imy) ((RIMAG(it,imx,imy),it=1,Ntau), imx=1,Nmx1)
            END DO
            close(IFN03)
        END IF

    END DO ! End loop for Iphr

    IF (MYID .EQ. 0) THEN
        DEALLOCATE(DATW)
    END IF

    DEALLOCATE(SENDBUF)
    DEALLOCATE(TMP2D)
    DEALLOCATE(vel3d)
    DEALLOCATE(RIMAG)

    RETURN
END SUBROUTINE OFFSET_PLANEWAVE_FD_PSTM_3D

!=====================================================================*

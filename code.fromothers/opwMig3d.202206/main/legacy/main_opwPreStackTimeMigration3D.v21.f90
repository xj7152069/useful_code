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
    CHARACTER*1024 ProcessorName
    INTEGER  ierr, myid, nproc, NAME_LENGTH

    !     Parameters reading from parameter file.
    REAL     F1, F2, F3, F4
    INTEGER  Iphr1, Iphr2, ViscoFlag, NtauExtraP
    CHARACTER*1024 FN_PAR, FN1, FN2, FN3, FN4

    !Parameters for 3D offset-planewave decompostion data.
    INTEGER  LT, Nmx, Nmy, Nphr
    INTEGER  Ncdp_first, Ncdp_final, Ncdp_step
    INTEGER  Nline_first, Nline_final, Nline_step
    REAL     DT, Dmx, Dmy, Phrmin, Dphr

    !Parameters for 3D interval velocity field(time domain).
    INTEGER  Nvtau, Nvmy, Nvmx
    INTEGER  Nvtau_first, Nvtau_final, Nvtau_step
    INTEGER  Nvline_first, Nvline_final, Nvline_step
    INTEGER  Nvcdp_first, Nvcdp_final, Nvcdp_step
    REAL     Dvtau, Dvx, Dvy

    INTEGER  Nmx0, Nmx1, Nmy0, Nmy1
    INTEGER  Ntau, NW1, NW2, NW3, NW4, NW
    REAL     Dtau, DW, W1, W2, W3, W4

    REAL     PAI2
    INTEGER  LByte, IFN01, IFN02, IFN03, IFN04

    CALL MPI_INIT(ierr)
    CALL MPI_COMM_RANK(MPI_COMM_WORLD, myid,  ierr)
    CALL MPI_COMM_SIZE(MPI_COMM_WORLD, nproc, ierr)
    IF (nproc .eq. 1) THEN
        write(*,*) '*************************************************'
        write(*,*) '*                                               *'
        write(*,*) '*  Number of processor should be larger than 1. *'
        write(*,*) '*  Program Exit Now!                            *'
        write(*,*) '*                                               *'
        write(*,*) '*************************************************'
        GOTO 9999
    END IF

    IF( 0 .eq. myid ) THEN
        CALL GETARG(1, FN_PAR)
    END IF
    CALL MPI_BCAST(FN_PAR, 1024, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
    CALL READPAR(FN_PAR, FN1, FN2, FN3, FN4, &
        Nvtau, Iphr1, Iphr2, F1, F2, F3, F4,&
        ViscoFlag, NtauExtraP)

    !Init. of parameters.
    Lbyte = 1
    PAI2  = 2*3.14159265359
    Ntau  = NtauExtraP
    IFN01 = 11
    IFN02 = 12
    IFN03 = 13
    IFN04 = 14

    !For Luojia3D data:
    !cdp=[2401,3700], ncdp=1300, dcdp=6.25m.
    !line=[1991,2460], nline=470, dline=12.5m.
!   LT          = 2901
!   DT          = 2.0         !ms
!   Nmx         = 1300
!   Ncdp_first  = 2401
!   Ncdp_final  = 3700
!   Ncdp_step   = 1
!   Dmx         = 6.25        !m
!   Nmy         = 470
!   Nline_first = 1991
!   Nline_final = 2460
!   Nline_step  = 1
!   Dmy         = 12.5        !m
!   Nphr        = 81
!   Phrmin      = 0
!   Dphr        = 10.0/1.0E6   !s/m
    !For synthetic Data:
!    LT          = 4000
!    DT          = 1.0         !ms
!    Nmx         = 500
!    Ncdp_first  = 151
!    Ncdp_final  = 650
!    Ncdp_step   = 1
!    Dmx         = 20.0
!    Nmy         = 500
!    Nline_first = 151
!    Nline_final = 650
!    Nline_step  = 1
!    Dmy         = 20.
!    Nphr        = 61
!    Phrmin      = 0
!    Dphr        = 10.0/1.0E6   !s/m

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
    !Nvtau         = 3010
!   Nvtau_first   = 1
!   Nvtau_final   = Nvtau
!   Nvtau_step    = 1
!   Nvmy          = 800
!   Nvline_first  = 1801
!   Nvline_final  = 2600
!   Nvline_step   = 1
!   Nvmx          = 2000
!   Nvcdp_first   = 2001
!   Nvcdp_final   = 4000
!   Nvcdp_step    = 1
!   Dvtau         = 2.0       !ms
!   Dvy           = 12.5      !m
!   Dvx           = 6.25      !m

    !For synthetic Data:
    !Nvtau         = 4000
!    Nvtau_first   = 1
!    Nvtau_final   = Nvtau
!    Nvtau_step    = 1
!    Nvmy          = 800
!    Nvline_first  = 1
!    Nvline_final  = 800
!    Nvline_step   = 1
!    Nvmx          = 800
!    Nvcdp_first   = 1
!    Nvcdp_final   = 800
!    Nvcdp_step    = 1
!    Dvtau         = 1.0       !ms
!    Dvy           = 20.0
!    Dvx           = 20.0

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

        write(*,*) '                                                 '
        write(*,*) '*************************************************'
        write(*,*) '*                                               *'
        write(*,*) '*  3D offset plane-wave FD PSTM begin.          *'
        write(*,*) '*                                               *'
        write(*,*) '*************************************************'
        write(*,*) '                                                 '
    END IF
    GOTO 9999

    CALL OFFSET_PLANEWAVE_FD_PSTM_3D(&
        myid, nproc, Lbyte, ViscoFlag, &
        Nvline_first, Nvline_final, Nvline_step, Nvmy, Dvx,&
        Nvcdp_first,  Nvcdp_final,  Nvcdp_step,  Nvmx, Dvy,&
        Nvtau, NtauExtraP,&
        Nline_first, Nline_final, Nline_step, Nmy, Nmy0, Nmy1, Dmy,&
        Ncdp_first,  Ncdp_final,  Ncdp_step,  Nmx, Nmx0, Nmx1, Dmx,&
        Nphr, Iphr1, Iphr2, Phrmin, Dphr,&
        NW1, NW2, NW3, NW4, NW, DW, LT, Ntau, Dtau,&
        IFN01, IFN02, IFN03, IFN04,&
        FN1, FN2, FN3, FN4)

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
        myid, nproc, Lbyte, ViscoFlag, &
        Nvline_first, Nvline_final, Nvline_step, Nvmy, Dvx,&
        Nvcdp_first,  Nvcdp_final,  Nvcdp_step,  Nvmx, Dvy,&
        Nvtau, NtauExtraP,&
        Nline_first, Nline_final, Nline_step, Nmy, Nmy0, Nmy1, Dmy,&
        Ncdp_first,  Ncdp_final,  Ncdp_step,  Nmx, Nmx0, Nmx1, Dmx,&
        Nphr, Iphr1, Iphr2, Phrmin, Dphr,&
        NW1, NW2, NW3, NW4, NW, DW, LT, Ntau, Dtau,&
        IFN01, IFN02, IFN03, IFN04,&
        FN1, FN2, FN3, FN4)

    implicit none
    INCLUDE 'mpif.h'
    INTEGER  status(MPI_STATUS_SIZE)
    INTEGER  ierr, myid, nproc
    INTEGER  LByte, IFN01, IFN02, IFN03, IFN04

    INTEGER  ViscoFlag, NtauExtraP
    INTEGER  Nvtau, Nvmx, Nvmy
    INTEGER  Nvline_first, Nvline_final, Nvline_step
    INTEGER  Nvcdp_first,  Nvcdp_final,  Nvcdp_step
    REAL     Dvtau, Dvx, Dvy
    INTEGER  Nline_first, Nline_final, Nline_step, Nmy, Nmy0, Nmy1
    INTEGER  Ncdp_first,  Ncdp_final,  Ncdp_step,  Nmx, Nmx0, Nmx1
    INTEGER  Ntau, LT, Nphr, Iphr1, Iphr2
    INTEGER  NW1, NW2, NW3, NW4, NW
    REAL     Dtau, Dmx, Dmy, Phrmin, Dphr
    REAL     DW, W1, W2, W3, W4, Dky, Dkx, W
    REAL     Rkz_max, RPh

    ! internal buffers.
    REAL,    ALLOCATABLE :: Trace(:), VVtau(:,:)
    REAL,    ALLOCATABLE :: vel3d(:,:,:), vel3dall(:,:,:)
    REAL,    ALLOCATABLE :: RIMAG(:,:,:), DATW_R(:,:,:)
    COMPLEX, ALLOCATABLE :: SENDBUF(:), CTrace(:)
    COMPLEX, ALLOCATABLE :: DATW_W(:,:), TMP2D(:,:)
    COMPLEX, ALLOCATABLE :: DATW(:,:,:)

    !parameters for writing temp files.
    CHARACTER*1024 FN1, FN2, FN3, FN4, FNCIG
    CHARACTER*4    ciphr
    REAL*4         sizeGB, LRIMAG, LDATW, LDATW_R, Lvel, LmemAll

    ! internal parameters.
    REAL    :: PAI, PAI2
    INTEGER  Ncount2d, NumSend, sender, itag, SlaveNodeIndex
    INTEGER  ierror1, ierror2, ierror3, ierror4, ierror5, ierror6
    INTEGER  imx, ivx, ix, imy, ivy, iy
    INTEGER  Iphr, it, itau, jtau, IDX, Iww, IW
    REAL*8   stime, etime
    INTEGER  iDebug    !=1, output the common-Phr data file.

    !Init. of parameters.
    iDebug  = 0
    sizeGB  = 1.0*1024*1024*1024
    PAI     = 3.14159265359
    PAI2    = 2*3.14159265359
    Rkz_max = cos(PAI2*85./360.)
    Dky     = PAI2/(Nmy*Dmy)
    Dkx     = PAI2/(Nmx*Dmx)

    !allocate the internal buffers.
    ierror1 = 0
    ierror2 = 0
    ierror3 = 0
    ierror4 = 0
    ierror5 = 0
    ierror6 = 0
    ALLOCATE(Trace(LT),  stat=ierror1)
    ALLOCATE(CTrace(LT), stat=ierror1)
    ALLOCATE(SENDBUF(Nmx*Nmy+2),  stat=ierror1)
    ALLOCATE(VVtau(Nmx,Nmy),      stat=ierror1)
    ALLOCATE(DATW_W(Nmx,Nmy),     stat=ierror1)
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
        !Master Node.
        ierror3 = 0
        ALLOCATE(DATW(Nmx,Nmy,NW), stat=ierror3)
        IF (ierror3 .ne. 0) THEN
            write(*,*) ' Not Enough Memory on Master Node!'
            write(*,*) ' Can not Allocate Mem for DATW!'
            write(*,*) ' Program Exit Now!'
            STOP
        END IF
        LDATW   = 1.*NW*Nmx*Nmy
        LRIMAG  = 1.*Ntau*Nmx1*Nmy1
        LDATW_R = 1.*Ntau*Nmx1*Nmy1
        Lvel    = 1.*Nmx*Nmy*Ntau + 1.*Nvmx*Nvmy*Nvtau
        LmemAll = 8*LDATW + 4*LRIMAG + 4*LDATW_R + 4*Lvel
        write(*,*) 'sizeof(DATW)= ',   8.*LDATW/sizeGB, 'GB'
        write(*,*) 'sizeof(RIMAG)= ',  4.*LRIMAG/sizeGB, 'GB'
        write(*,*) 'sizeof(DATW_R)= ', 4.*LDATW_R/sizeGB, 'GB'
        write(*,*) 'sizeof(vel_loc)= ',4.*Nmx*Nmy*Ntau/sizeGB, 'GB'
        write(*,*) 'sizeof(vel_all)= ',4.*Nvmx*Nvmy*Nvtau/sizeGB, 'GB'
        write(*,*) 'Total Memory Needed: ', LmemAll/sizeGB, 'GB'
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
        write(*,*) ' Dky    = ', Dky
        write(*,*) ' Dkx    = ', Dkx
        write(*,*) ' Rkz_max    = ', Rkz_max
    END IF

    !Each node read the velocity model into memory.
    write(*,*) 'Now, reading the 3-D velocit model. myid is', myid
    OPEN(IFN02, file=trim(adjustl(FN2)), access='direct', recl=Nvmx*Nvmy)
    DO Itau = 1, Nvtau
        !      write(*,*) 'Itau=',Itau, 'Nvtau=', Nvtau
        read(IFN02, rec=Itau) ((vel3dall(ix,iy,Itau),ix=1,Nvmx),iy=1,Nvmy)
    END DO
    CLOSE(IFN02)
    write(*,*) 'Now, reading the 3-D velocit model is done. myid is', myid

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

    IF (myid .eq. 0) THEN
       OPEN(21, file='velo_check_skip50.dat', access='direct', recl=Nmx*Nmy)
       jtau = 1
       DO Itau = 1, Ntau, 50
          DO Imy = 1, Nmy
          DO Imx = 1, Nmx
             VVtau(Imx, Imy) = vel3d(Imx, Imy, Itau)
          END DO
          END DO
          write(*,*) 'Itau=',Itau, 'Jtau=', jtau, 'Ntau=', Ntau
          WRITE(21,REC=jtau)((VVtau(IX,IY),IX=1,Nmx),IY=1,Nmy)
          jtau = jtau + 1
       END DO
       CLOSE(21)
    END IF

    CALL MPI_BARRIER(MPI_COMM_WORLD, IERR)

    DO Iphr = Iphr1, Iphr2

        stime   = MPI_Wtime()
        DATW_R  = 0.0
        RIMAG   = 0.0
        IF (myid .eq. 0) THEN

            !DATW_R  = 0.0
            !RIMAG   = 0.0
            RPh     = Phrmin + Dphr*(Iphr-1)
            write(*,*) ' Iphr = ',Iphr,' RPh= ',RPh*1.0E6,'vs/m'

            write(*,*) ' ****READ PlaneWave Data and FFT With OpenMP.'
            CALL READ_OPW_DATA(&
                FN1, FN3, IFN01, IFN03, iDebug, DATW,&
                Nmx, Nmy, NW1, NW2, NW3, NW4, NW, LT,&
                Ncdp_first,  Ncdp_final,  Ncdp_step,&
                Nline_first, Nline_final, Nline_step,&
                Nphr, Iphr, Nmx0, Nmx1, Nmy0, Nmy1)

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

        etime   = MPI_Wtime()
        !write the common-ph migrated section.
        IF (myid .eq. 0) THEN
            write(ciphr,'(I4)') Iphr
            FNCIG=trim(adjustl(FN4))//trim(adjustl(ciphr))//'.dat'

            write(*,*) 'Iphr=', Iphr, 'computation time(s) is:', etime-stime
            write(*,*) 'Now, wrting file:', FNCIG
            write(*,*) 'Ntau=', Ntau, 'Nmx1=', Nmx1, 'Nmy1=', Nmy1

            OPEN(IFN04, file=FNCIG, access='direct', recl=Ntau*Nmx1)
            DO imy = 1, Nmy1
                write(IFN04, rec=imy) ((RIMAG(it,imx,imy),it=1,Ntau), imx=1,Nmx1)
            END DO
            CLOSE(IFN04)
        END IF
        CALL MPI_BARRIER(MPI_COMM_WORLD, IERR)

    END DO ! End loop for Iphr

    IF (MYID .EQ. 0) THEN
        DEALLOCATE(DATW)
    END IF

    DEALLOCATE(Trace, CTrace, SENDBUF)
    DEALLOCATE(VVtau, DATW_W, TMP2D)
    DEALLOCATE(vel3d, DATW_R, RIMAG)

    RETURN
END SUBROUTINE OFFSET_PLANEWAVE_FD_PSTM_3D

!=====================================================================*

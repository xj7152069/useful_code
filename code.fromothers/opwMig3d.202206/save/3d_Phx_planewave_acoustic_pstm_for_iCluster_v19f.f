*==================================================================*
*                                                                  *
*           3D OFFSET PLANE-WAVE PRESTACK TIME MIGRATION           *
*           WITH FINITE DIFFERENCE METHOD PLUS F-K DMOAIN          *
*           ERROR COMPENSATION FOR 3D ACOUSTIC MEDIA.              *
*==================================================================*
*              ALL   RIGHTS   RESERVED                             *
*==================================================================*
*       School of Ocean and Earth Science, TongJi University       *
*       and IGP, Nanjing.                                          *
*==================================================================*
*       Remarks:                                                   *
*         Author: Feng Bo.                                         *
*         Date  : 2008.12 -- 2010.07                               *
*         Current Version: v18                                     *
*==================================================================*
*       History:                                                   *
*         Stable Version: v16.                                     *
*         Date  : 2010.07                                          *
* Name:     3d_Phx_planewave_acoustic_pstm_for_iCluster_v16.f      *
* Package:  pw3dlib.update.100720.tar                              *
*         Stable Version: v14.                                     *
*         Date  : 2009.08                                          *
* Name:     3d_Phx_planewave_acoustic_pstm_for_iCluster_v14.f      *
* Package:  pw3dlib.update.090820.tar                              *
*==================================================================*
*         Logs.                                                    *
*==================================================================*
*         Updated for V18, 2010-07-29                              *
*           adding Tau-Axis Break-Point protection feature.        *
*         Updated for V17, 2010-07-24                              *
*           Adding Visco-Acoustic Migration flag.                  *
*           ViscoFlag: =0, acoustic mig; =1, visco-acoustic mig.   *
*         Updated for V16, 2010-07-17                              *
*           Modified Error&Anomaly handling.                       *
*           Adding Open file status checking.                      *
*           Adding Break-Point protection feature.                 *
*           Using Multi-Threads FFT-2D. FFT2D_OMP().               *
*         Updated for V15, 2010-07-15                              *
*           Modified regular-file writting, from data() to file(). *
*           Using Multi-threads FFT-1D. **Profile2d().             *
*==================================================================*
*       Summary:                                                   *
*         This program do a 3D offset plane-wave prestack          *
*       time migration with finite-difference method               *
*       plus W-K domain error compensation for acoustic media.     *
*==================================================================*

      PROGRAM PRESTACK_3D_PLANE_WAVE_FD_TIME_MIGRATION 
      INCLUDE 'mpif.h'
      PARAMETER(Lbyte=1)
      PARAMETER(PAI2=2*3.14159265359)

!     Parameters reading from parameter file.
      REAL     F1, F2, F3, F4
      INTEGER  Ntau, Iphr1, Iphr2, IBREAK, ViscoFlag
      INTEGER  NtauExtraP, ITauBreak, Ktaur, Ktauw
      INTEGER  IphrR1, IphrR2
      CHARACTER*256 FN1, FN2, FN3, FN4, FNTMPW, FNTMPR

!     Parameters for 3D offset-planewave decompostion data.
      INTEGER  IPPAR(25), IFN01
      REAL     FPPAR(20)
      INTEGER  LT, Nmx, Nmy, Nphr
      INTEGER  Ncdp_first, Ncdp_final, Ncdp_step
      INTEGER  Nline_first, Nline_final, Nline_step
      REAL     DT, Dmx, Dmy, Phrmin, Dphr

!     Parameters for 3D interval velocity field(time domain).
      INTEGER  IVPAR(25), IFN02
      REAL     FVPAR(10)
      INTEGER  Nvtau, Nvmy, Nvmx
      INTEGER  Nvtau_first, Nvtau_final, Nvtau_step
      INTEGER  Nvline_first, Nvline_final, Nvline_step
      INTEGER  Nvcdp_first, Nvcdp_final, Nvcdp_step
      REAL     Dvtau, Dvx, Dvy

!     Parameters for common imaging gather of 3D offset plane-wave PSTM.
      INTEGER  icpar(25), IFN03, iflag3
      REAL     cpar(20)
      REAL     Dtau

      INTEGER  Nmx0, Nmx1, Nmy0, Nmy1
      INTEGER  NW1, NW2, NW3, NW4, NW
      REAL     DW, W1, W2, W3, W4

      REAL,ALLOCATABLE :: BUF(:)
      INTEGER*8  L1,  L2,  L3,  L4,  L5,  L6,  L7,  L8,  L9,  L10
      INTEGER*8  L11, L12, L13, L14, L15, L16, L17, L18, L19, L20

      CHARACTER*256 ProcessorName
      INTEGER  ierr, myid, np, NAME_LENGTH

*======================================================================*

      CALL MPI_INIT(ierr)
      CALL MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD, np, ierr)

      Call READPAR(FN1, FN2, FN3, FN4,
     +         Ntau, Iphr1, Iphr2, F1, F2, F3, F4, IBREAK, ViscoFlag,
     +         NtauExtraP, ITauBreak, Ktaur, Ktauw, FNTMPW, FNTMPR,
     +         IphrR1, IphrR2)

!     FN1: Offset Plane-Wave data file.
      CALL GDBFCGFENTERFILENAME(FN1)
!     FN2: Interval velocity field(time domain, slice data).
      CALL GDBFCGFENTERFILENAME(FN2)
!     FN3: Common-Ph imaging volume, Ntau*Nmx*Nmy*Nphr.
      CALL GDBFCGFENTERFILENAME(FN3)
*======================================================================*

!     Open the offset plane-wave file(FN1) for reading.  (LT*Nmx*Nmy*Nphr)
      CALL gdbfcinitregulardata(FN1, IPPAR, FPPAR, IFN01)

      IF ( -99999 .eq. IPPAR(6) ) THEN
         CALL MPI_GET_PROCESSOR_NAME(ProcessorName, NAME_LENGTH, IERR)
         write(*,*) ' Open Plane-Wave data file error at Node: ',
     +              TRIM(ProcessorName), '  myid is: ', myid
         write(*,*) ' FLAG = ', IPPAR(6)
         write(*,*) ' Program Exit. '
         GOTO 9999
      END IF

      LT          = IPPAR(5)          ! 1st Dimension grid number.
      DT          = FPPAR(6)          ! sampling rate of 1st Dimension grid.
!     DT          = FPPAR(4)          ! sampling rate of 1st Dimension grid.

      Nmx         = IPPAR(9)          ! 2nd Dimension grid number.
      Ncdp_first  = IPPAR(10)         ! first grid number of 2nd Dimension.
      Ncdp_final  = IPPAR(11)         ! last  grid number of 2nd Dimension.
      Ncdp_step   = IPPAR(12)         ! 1st Dimension grid step.
      Dmx         = FPPAR(8)          ! sampling rate of 2nd Dimension grid.
!     Dmx         = FPPAR(6)          ! sampling rate of 2nd Dimension grid.

      Nmy         = IPPAR(13)         ! 3rd Dimension grid number.
      Nline_first = IPPAR(14)         ! first grid number of 3rd Dimension.
      Nline_final = IPPAR(15)         ! last  grid number of 3rd Dimension.
      Nline_step  = IPPAR(16)         ! 3rd Dimension grid step.
      Dmy         = FPPAR(10)         ! sampling rate of 3rd Dimension grid.
!     Dmy         = FPPAR(8)          ! sampling rate of 3rd Dimension grid.

      Nphr        = IPPAR(17)         ! 4th Dimension grid number.
      Phrmin      = FPPAR(9)          ! origin coordinate of 4th Dimension grid.
      Dphr        = FPPAR(12)         ! sampling rate of 4th Dimension grid.
!     Dphr        = FPPAR(10)         ! sampling rate of 4th Dimension grid.
*======================================================================*

!     New requirement will be included recently.
      IF (Iphr1 .eq. 0 .and. Iphr2 .eq. 0) THEN
          Iphr1 = 1
          Iphr2 = Nphr
      END IF
      IF (Iphr1 .gt. Nphr .or. Iphr2 .gt. Nphr .or.
     +   Iphr1 .lt. 1 .or. Iphr2 .lt. 1 .or. Iphr1 .gt. Iphr2) THEN
          write(*,*) 'Wrong index of Ray-parameter.'
          write(*,*) ' Iphr1 = ', Iphr1
          write(*,*) ' Iphr2 = ', Iphr2
          write(*,*) ' Nphr  = ', Nphr 
          GOTO 6001
      END IF
!     Update the number of ray-parameters: Nphr
      Nphr = Iphr2 - Iphr1 + 1

*======================================================================*
!     Open the velocity file(FN2) for reading.
      CALL GDBFCINITCUTSLICEDATA(FN2, IVPAR, FVPAR, IFN02)

      IF ( -99999 .eq. IVPAR(6) ) THEN
         CALL MPI_GET_PROCESSOR_NAME(ProcessorName, NAME_LENGTH, IERR)
         write(*,*) ' Open velocity file error at Node: ',
     +              TRIM(ProcessorName), '  myid is: ', myid
         write(*,*) ' FLAG = ', IVPAR(6)
         write(*,*) ' Program Exit. '
         GOTO 6001
      END IF

      Nvtau         = IVPAR(2)        ! Sample number of the 3rd dim. 
      Nvtau_first   = IVPAR(3)        ! First sample number of the 3rd dim.
      Nvtau_final   = IVPAR(4)        ! Last  sample number of the 3rd dim.
      Nvtau_step    = IVPAR(5)        ! Sampling   interval of the 3rd dim.
      Nvmy          = IVPAR(6)        ! Sample number of the 2nd dim.
      Nvline_first  = IVPAR(7)        ! First sample number of the 2nd dim.
      Nvline_final  = IVPAR(8)        ! Last  sample number of the 2nd dim.
      Nvline_step   = IVPAR(9)        ! Sampling   interval of the 2nd dim.
      Nvmx          = IVPAR(10)       ! Sample number of the 1st dim.
      Nvcdp_first   = IVPAR(11)       ! First sample number of the 1st dim.
      Nvcdp_final   = IVPAR(12)       ! Last  sample number of the 1st dim.
      Nvcdp_step    = IVPAR(13)       ! Sampling   interval of the 1st dim.

      Dvtau         = FVPAR(3)        ! Real grid interval of the 3rd dim.
      Dvy           = FVPAR(4)        ! Real grid interval of the 2nd dim.
      Dvx           = FVPAR(5)        ! Real grid interval of the 1st dim.
*======================================================================*

!     Check the input parameters.
!     This part should be tested carefully.

!     Check the 3d velocity field.
      IF ((Nvline_first .gt. Nline_first) .or.
     +   (Nvline_final .lt. Nline_final) .or.
     +   (Nvcdp_first  .gt. Ncdp_first) .or.
     +   (Nvcdp_final  .lt. Ncdp_final)) THEN
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
         GOTO 6002
      END IF
      IF (Dmx .le. 1.0E-2 .or. Dvx .le. 1.0E-2 .or.
     +   Dmy .le. 1.0E-2 .or. Dvy .le. 1.0E-2) THEN
         write(*,*) ' Velocity field sampling error. '
         write(*,*) ' Or Plane-wave data sampling error. '
         write(*,*) ' Program Exit. '
         GOTO 6002
      END IF
      IF (Dmx/Dvx .lt. 0.95 .or. Dmx/Dvx .gt. 1.05 .or.
     +   Dmy/Dvy .lt. 0.95 .or. Dmy/Dvy .gt. 1.05) THEN
         write(*,*) ' Velocity Field Grid Error!'
         write(*,*) ' Dmx = ', Dmx
         write(*,*) ' Dvx = ', Dvx
         write(*,*) ' Dmy = ', Dmy
         write(*,*) ' Dvy = ', Dvy
         write(*,*) ' Program Exit. '
         GOTO 6002
      END IF

      IF ((Nvtau .lt. Ntau)) THEN
         write(*,*) ' Extrapolation Number Error!'
         write(*,*) ' Nvtau = ', Nvtau
         write(*,*) ' Ntau  = ', Ntau
         Ntau = Nvtau
         write(*,*) ' Ntau is set to Nvtau automaticlly'
         write(*,*) ' Ntau  = ', Ntau
      END IF
*======================================================================*

      IF (myid .eq. 0) THEN
       IF(IBREAK.eq.0) THEN
!      BREAK POINT SAVING.
         icpar(1)  = 4                   ! Dimension Number
         icpar(2)  = 1                   ! Domain 1:time 2:depth
         icpar(3)  = 1                   ! Data Format 1:float 2:int 3:short...
         icpar(4)  = 1                   ! File Type 1:Seismic 2:Velocity 3:Other
C        cpar(1)   = 0                   ! Value Minimum (Write File Don't use)
C        cpar(2)   = 0                   ! Value Maximum (Write File Don't use)
   
         icpar(5)  = Ntau                ! 1st Dimension grid number.
         icpar(6)  = 1                   ! first grid number of 1st Dimension.
         icpar(7)  = Ntau                ! last  grid number of 1st Dimension.
         icpar(8)  = 1                   ! 1st Dimension grid step.
         cpar(1)   = 0                   ! origin coordinate of 1st Dimension grid.
         cpar(2)   = Dvtau               ! sampling rate of 1st Dimension grid(ms).
C        cpar(3)   = 0                   ! origin coordinate of 1st Dimension grid.
C        cpar(4)   = Dvtau               ! sampling rate of 1st Dimension grid(ms).

         icpar(9)  = Nmx                 ! 2nd Dimension grid number.
         icpar(10) = Ncdp_first          ! first grid number of 2nd Dimension.
         icpar(11) = Ncdp_final          ! last  grid number of 2nd Dimension.
         icpar(12) = Ncdp_step           ! 2nd Dimension grid step.
         cpar(3)   = 0                   ! origin coordinate of 2nd Dimension grid.
         cpar(4)   = Dmx                 ! sampling rate of 2nd Dimension grid.
C        cpar(5)   = 0                   ! origin coordinate of 2nd Dimension grid.
C        cpar(6)   = Dmx                 ! sampling rate of 2nd Dimension grid.

         icpar(13) = Nmy                 ! 3rd Dimension grid number.
         icpar(14) = Nline_first         ! first grid number of 3rd Dimension.
         icpar(15) = Nline_final         ! last  grid number of 3rd Dimension.
         icpar(16) = Nline_step          ! 3rd Dimension grid step.
         cpar(5)   = 0                   ! origin coordinate of 3rd Dimension grid.
         cpar(6)   = Dmy                 ! sampling rate of 3rd Dimension grid.
C        cpar(7)   = 0                   ! origin coordinate of 3rd Dimension grid.
C        cpar(8)   = Dmy                 ! sampling rate of 3rd Dimension grid.

         Ntr = Iphr2 - Iphr1 + 1
         icpar(17) = Ntr                 ! 4th Dimension grid number.
         icpar(18) = Iphr1               ! first grid number of 4th Dimension.
         icpar(19) = Iphr2               ! last  grid number of 4th Dimension.
         icpar(20) = 1                   ! 4th Dimension grid step.
         cpar(7)   = 0                   ! origin coordinate of 4th Dimension grid.
         cpar(8)   = Dphr                ! sampling rate of 4th Dimension grid.
C        cpar(9)   = 0                   ! origin coordinate of 4th Dimension grid.
C        cpar(10)  = Dphr                ! sampling rate of 4th Dimension grid.

!        Open the Common-Ph imaging volume file (FN3), using regular data-4D format.
         iflag3 = 0
!        CALL gdbfcinitwriteregulardata(FN3, icpar, cpar, iflag3, IFN03)
         CALL gdbfcinitwriteregularfile(FN3, icpar, cpar, iflag3, IFN03)
         IF(iflag3 .ne. 0) THEN
            write(*,*) ' gdbfcinitwriteregularfile() failed.'
            write(*,*) ' iflag3 = ', iflag3
            GOTO 6002
         END IF
         write(*,*) '**********Initialize the CIG file done.**********'

!        Link Imaging-volume with velocity.
         iflag = 0
         CALL gsetvelocityregularfile(IFN03,FN2,iflag)
         IF(iflag .ne. 0) THEN
            write(*,*) 'gsetvelocityregularfile() failed.'
            write(*,*)' iflag = ', iflag
            GOTO 6003
         END IF
         write(*,*) '**********Linking   velocity file done.**********'

       ELSE
!      BREAK POINT SAVING.
         iflag3 = 1
         CALL gdbfcinitwriteregularfile(FN3, icpar, cpar, iflag3, IFN03)
         write(*,*) '**********Beak-Point Mode Continue.**********'
         IF(iflag3 .ne. 0) THEN
            write(*,*) ' gdbfcinitwriteregularfile() failed.'
            write(*,*) ' iflag3 = ', iflag3
            GOTO 6002
         END IF
       END IF

      END IF
*======================================================================*

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
         GOTO 6003
      END IF

      Dtau = Dvtau*0.001     !unit is second.

      Nmx0 = 100
      Nmx1 = Nmx
      Nmx  = Nmx1 + 2*Nmx0

      Nmy0 = 40
      Nmy1 = Nmy
      Nmy  = Nmy1 + 2*Nmy0

*======================================================================*
*     Print Parameters on the screen.
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
         write(*,*) ' NW = ', NW, ' processor number:', np
         write(*,*) '                                                 '
         write(*,*) '*************************************************'
      END IF

*======================================================================*
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
      L19= L18 + Ntau*Nmx1*Nmy1  !DATW_R

      ierror_flag = 0
      ALLOCATE(BUF(L19), stat=ierror_flag)
      IF (ierror_flag .ne. 0) THEN
         CALL MPI_GET_PROCESSOR_NAME(ProcessorName, NAME_LENGTH, IERR)
         write(*,*) ' Not Enough Memory at Node: ',
     +              TRIM(ProcessorName), '  myid is: ', myid
         write(*,*) ' Program Exit. '
         GOTO 6003
      END IF

*======================================================================*
      IF (myid .eq. 0) THEN
         write(*,*) '                                                 '
         write(*,*) ' Total Memory Need: ', (L19-1)*4.0/1024/1024, 'MB' 
         write(*,*) '                                                 '
         write(*,*) '*************************************************'
         write(*,*) '*                                               *'
         write(*,*) '*  3D offset plane-wave FD PSTM begin.          *'
         write(*,*) '*                                               *'
         write(*,*) '*************************************************'
         write(*,*) '                                                 '
      END IF
C     pause

      IF (np .eq. 1) THEN
         write(*,*) '*************************************************'
         write(*,*) '*                                               *'
         write(*,*) '*  Require processor number larger than 1.      *'
         write(*,*) '*  Program Exit Now!                            *'
         write(*,*) '*                                               *'
         write(*,*) '*************************************************'
         GOTO 10715
      END IF

*======================================================================*
      CALL OFFSET_PLANEWAVE_FD_PSTM_3D(
     +     BUF(L1), BUF(L2), BUF(L3),
     +     BUF(L4),  BUF(L5),  BUF(L6),  BUF(L7),  BUF(L8),  BUF(L9),
     +     BUF(L10), BUF(L11), BUF(L12), BUF(L13), BUF(L14), BUF(L15),
     +     BUF(L16), BUF(L17), BUF(L18),
     +     Nline_first, Nline_final, Nline_step, Nmy, Nmy0, Nmy1, Dmy,
     +     Ncdp_first,  Ncdp_final,  Ncdp_step,  Nmx, Nmx0, Nmx1, Dmx,
     +     Nphr, Iphr1, Iphr2, Phrmin, Dphr,
     +     NW1, NW2, NW3, NW4, NW, DW, LT, Ntau, Dtau,
     +     IFN01, IFN02, IFN03, FN4, FNTMPW, FNTMPR,
     +     ViscoFlag, NtauExtraP, IphrR1, IphrR2, 
     +     ITauBreak, Ktaur, Ktauw, Lbyte, myid, np)

10715 CONTINUE

*======================================================================*
!     Close the 5d planewave file(FN1).
      IFLAG = 0      ! only close the file.
      CALL gdbfccloseregulardata(IFN01, IFlAG)
      IF (myid .eq. 0) THEN
         write(*,*) '**Close the Plane-Wave file: ', TRIM(FN1), ' done.'
      END IF
!     Close the velocity file(FN2).
      IFLAG = 0      ! only close the file.
      CALL GDBFCCLOSECUTSLICEDATA(IFN02, IFLAG)
      IF (myid .eq. 0) THEN
         write(*,*) '**Close the Velocity   file: ', TRIM(FN2), ' done.'
      END IF
!     Close the CIG file(FN3).
      IF (myid .eq. 0) THEN
         iflag3 = 0
         CALL gdbfccloseregulardata(IFN03, iflag3)
         write(*,*) '**Building Index of file:',    TRIM(FN3), ' done.'
      END IF

*======================================================================*
      DEALLOCATE(BUF)
      GOTO 9999

*======================================================================*
!     Error Handling Area.
*======================================================================*
6003  CONTINUE
!     Error Handling Code.
      IF (myid .eq. 0) THEN
         iflag3 = 0
         CALL gdbfccloseregulardata(IFN03, iflag3)
      END IF

6002  CONTINUE
!     Error Handling Code.
      IFLAG = 0      ! only close the file.
      CALL GDBFCCLOSECUTSLICEDATA(IFN02, IFLAG)

6001  CONTINUE
!     Error Handling Code.
!     Close the 5d planewave file(FN1).
      IFLAG = 0      ! only close the file.
      CALL gdbfccloseregulardata(IFN01, IFlAG)
      GOTO 9999

*======================================================================*
!     Program Exit.
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
      END

*======================================================================*
      SUBROUTINE OFFSET_PLANEWAVE_FD_PSTM_3D(
     +     VVtau, Trace, CTrace,
     +     AA_X, BB_X, CC_X, DD_X, XX_X, YY_X,
     +     AA_Y, BB_Y, CC_Y, DD_Y, XX_Y, YY_Y,
     +     DATW_W, XX, DATW_R,
     +     Nline_first, Nline_final, Nline_step, Nmy, Nmy0, Nmy1, Dmy,
     +     Ncdp_first,  Ncdp_final,  Ncdp_step,  Nmx, Nmx0, Nmx1, Dmx,
     +     Nphr, Iphr1, Iphr2, Phrmin, Dphr,
     +     NW1, NW2, NW3, NW4, NW, DW, LT, Ntau, Dtau,
     +     IFN01, IFN02, IFN03, FN4, FNTMPW, FNTMPR,
     +     ViscoFlag, NtauExtraP, IphrR1, IphrR2,
     +     ITauBreak, Ktaur, Ktauw, Lbyte, myid, np)
          
      INCLUDE 'mpif.h'
      PARAMETER(PAI=3.14159265359, PAI2=2*3.14159265359)
      INTEGER  status(MPI_STATUS_SIZE)
      INTEGER  ierr, myid, np, Lbyte
      INTEGER  ViscoFlag, NtauExtraP, ITauBreak, Ktaur, Ktauw

      INTEGER  Nline_first, Nline_final, Nline_step, Nmy, Nmy0, Nmy1
      INTEGER  Ncdp_first,  Ncdp_final,  Ncdp_step,  Nmx, Nmx0, Nmx1
      INTEGER  Ntau, LT, Nphr, Iphr1, Iphr2
      INTEGER  NW1, NW2, NW3, NW4, NW, KmyNUM, KmxNUM
      INTEGER  IphrR1, IphrR2
      REAL     Dtau, DIP_max, DIP_min, Dmx, Dmy, Phrmin, Dphr
      REAL     DW, W1, W2, W3, W4, Dky, Dkx

      REAL     Trace(LT), VVtau(Nmx, Nmy), DATW_R(Ntau, Nmx1, Nmy1)
      REAL     ACOE(5), BCOE(5)

      COMPLEX  AA_X(Nmx), BB_X(Nmx), CC_X(Nmx), DD_X(Nmx)
      COMPLEX  XX_X(Nmx), YY_X(Nmx)
      COMPLEX  AA_Y(Nmy), BB_Y(Nmy), CC_Y(Nmy), DD_Y(Nmy)
      COMPLEX  XX_Y(Nmy), YY_Y(Nmy)
      COMPLEX  CTrace(LT), DATW_W(Nmx, Nmy), XX(Nmx, Nmy)

      INTEGER  ierror1, ierror2, ierror3, ierror4
      COMPLEX, ALLOCATABLE :: SENDBUF(:)
      COMPLEX, ALLOCATABLE :: TMP2D(:,:)
      COMPLEX, ALLOCATABLE :: DATW(:,:,:)
      REAL,    ALLOCATABLE :: RIMAG(:,:,:)

      INTEGER  LOOPP
      INTEGER  NumSend, sender, itag, SlaveNodeIndex

!     parameters for writing temp files.
      CHARACTER*256 FN4, FNTMPW, FNTMPR, FNW, FNR
      CHARACTER     c1, c2, c3
      INTEGER       fnlen
      INTEGER*8     LRIMAG, LDATW

C     Initialize some variables.
C     fnlen = len_trim(FNTMP)
C     FNTMP(fnlen+1:fnlen+1) = char(0)

      DIP_max = 79.0
      DIP_min = 25.0
      ALFA    = 0.113
      Dky     = PAI2/(Nmy*Dmy)
      Dkx     = PAI2/(Nmx*Dmx)
      Rkz_max = cos(PAI2*DIP_max/360.0)

      CALL make_order(DIP_max, Acoe, Bcoe, LOOPP)

      ierror1 = 0
      ierror4 = 0
      ALLOCATE(SENDBUF(Nmx*Nmy+2), stat=ierror1)
      ALLOCATE(TMP2D(Nmy,Nmx),     stat=ierror4)
      IF ((ierror1 .ne. 0) .or. (ierror4 .ne. 0)) THEN
           write(*,*) ' Not Enough Memory for 2D-buffers.'
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

      write(*,*) '********  Parameters of Extrapolation  ********'
      write(*,*) ' ALFA   = ', ALFA
      write(*,*) ' Dky    = ', Dky
      write(*,*) ' Dkx    = ', Dkx
      write(*,*) ' DIP_max    = ', DIP_max
      write(*,*) ' DIP_min    = ', DIP_min
      write(*,*) ' Rkz_max    = ', Rkz_max
      write(*,*) ' LOOPP, NP  = ', LOOPP, NP
      END IF


      IF (myid .eq. 0) THEN
!     Master Node.
         ierror2 = 0
         ierror3 = 0
         ALLOCATE(RIMAG(Ntau,Nmx1,Nmy1), stat=ierror2)
         ALLOCATE(DATW(Nmx,Nmy,NW),      stat=ierror3)
         IF ((ierror2 .ne. 0) .or. (ierror3 .ne. 0)) THEN
              write(*,*) ' Not Enough Memory on Master Node!'
              write(*,*) ' Can not Allocate Mem for RIMAG&DATW!'
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

      DO 1111 Iphr = Iphr1, Iphr2

         c1=char(Iphr/100+48)
         c2=char((Iphr-Iphr/100*100)/10+48)
         c3=char((Iphr-Iphr/10*10)+48)
         FNW=trim(FNTMPW)//c1//c2//c3//'.dat'
         FNR=trim(FNTMPR)//c1//c2//c3//'.dat'
C        write(*,*) 'Temp Writing File is', trim(FNW)
C        write(*,*) 'Temp Reading File is', trim(FNR)

         DO Itau = 1, Ntau
         DO IY   = 1, Nmy1
         DO IX   = 1, Nmx1
            DATW_R(Itau, IX, IY) = 0.0
         END DO
         END DO
         END DO

         IF (myid .eq. 0) THEN

            DO Itau = 1, Ntau
            DO IY   = 1, Nmy1
            DO IX   = 1, Nmx1
               RIMAG(Itau, IX, IY)  = 0.0
            END DO
            END DO
            END DO

C           DO IW   = 1, NW
C           DO IY   = 1, Nmy
C           DO IX   = 1, Nmx
C              DATW(IX, IY, IW)  = 0.0
C           END DO
C           END DO
C           END DO
C           CALL gdbfcwritecomplex(FNTMP, DATW, LDATW, Iphr,
C    +           1, filestat, 0)
            OPEN(66,FILE=FNW,ACCESS='DIRECT',
     +         STATUS='REPLACE',RECL=4)
            WRITE(66,REC=1) Iphr
            CLOSE(66)

            RPh     = Phrmin + Dphr*(Iphr-1)
            write(*,*) ' Iphr = ',Iphr,' RPh= ',RPh*1.0E6,'vs/m'

            write(*,*) ' ****READ PlaneWave Data and FFT With OpenMP.'
            CALL READ_PLANEWAVE_DATA_AND_FFT_Profile2d(
     +           Trace, DATW, CTrace,
     +           Nmx, Nmy, NW1, NW2, NW3, NW4, NW, LT, DW,
     +           Ncdp_first,  Ncdp_final,  Ncdp_step,
     +           Nline_first, Nline_final, Nline_step,
     +           Iphr, IFN01, Nmx0, Nmx1, Nmy0, Nmy1)

            NumSend = 1
            DO 2222 Iww = 1, MIN(np-1, NW)

               write(*,*) 'Iww =', Iww, ' Iphr = ', Iphr
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

               CALL MPI_Send(SENDBUF, Nmy*Nmx+2, MPI_COMPLEX,
     +                       Iww, Iww, MPI_COMM_WORLD, IERR)
               NumSend = NumSend + 1

2222        CONTINUE

            DO 3333 IW = 1, NW

C              Receive the imageing volume from slave node.
               CALL MPI_RECV(SlaveNodeIndex, 1, MPI_INTEGER,
     +                       MPI_ANY_SOURCE, MPI_ANY_TAG,
     +                       MPI_COMM_WORLD, status, ierr)

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

                  write(*,*) 'Iww =', NumSend, ' Iphr = ', Iphr
                  CALL MPI_Send(SENDBUF, Nmy*Nmx+2, MPI_COMPLEX,
     +                          sender, NumSend, MPI_COMM_WORLD, IERR)
                  NumSend = NumSend + 1
               ELSE
                  CALL MPI_Send(0, 0, MPI_COMPLEX, sender,
     +                          0, MPI_COMM_WORLD, IERR)
               END IF

3333        CONTINUE

         ELSE
!        Slave Nodes.
            IF (myid .gt. NW) GOTO 6666  ! Nodes number larger than tasks number.

            DO WHILE(1)

               CALL MPI_RECV(SENDBUF, Nmy*Nmx+2, MPI_COMPLEX,
     +                       0, MPI_ANY_TAG, 
     +                       MPI_COMM_WORLD, status, ierr)

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

C              goto 5168
               IF ( ViscoFlag .EQ. 0 ) THEN
C              3D Acoustic Plane-Wave Migration.
                  CALL WXFD_EXTRAPOLATION_WXYCOMP_MPI_ACOUSTIC(
     +                 DATW_W, DATW_R, XX, VVtau, TMP2D,
     +                 AA_X, BB_X, CC_X, DD_X, XX_X, YY_X,
     +                 AA_Y, BB_Y, CC_Y, DD_Y, XX_Y, YY_Y,
     +                 ACOE, BCOE, LOOPP, ALFA,
     +                 Nmx, Dmx, Nmy, Dmy, Ntau, Dtau, W, RPh,
     +                 Dkx, Dky, Rkz_max, DIP_max, DIP_min,
     +                 Nline_first, Nline_final,
     +                 Ncdp_first, Ncdp_final, IFN02,
     +                 Nmx0, Nmx1, Nmy0, Nmy1, Iphr, Iww, NtauExtraP,
     +                 ITauBreak, Ktaur, Ktauw, FNW, FNR,
     +                 IphrR1, IphrR2, myid)
               ELSE
C              3D Visco-Acoustic Plane-Wave Migration.
                  CALL WXFD_EXTRAPOLATION_WXYCOMP_MPI_VISCO(
     +                 DATW_W, DATW_R, XX, VVtau, TMP2D,
     +                 AA_X, BB_X, CC_X, DD_X, XX_X, YY_X,
     +                 AA_Y, BB_Y, CC_Y, DD_Y, XX_Y, YY_Y,
     +                 ACOE, BCOE, LOOPP, ALFA,
     +                 Nmx, Dmx, Nmy, Dmy, Ntau, Dtau, W, RPh,
     +                 Dkx, Dky, Rkz_max, DIP_max, DIP_min,
     +                 Nline_first, Nline_final,
     +                 Ncdp_first, Ncdp_final, IFN02,
     +                 Nmx0, Nmx1, Nmy0, Nmy1, Iphr, Iww, NtauExtraP,
     +                 ITauBreak, Ktaur, Ktauw, FNW, FNR,
     +                 IphrR1, IphrR2, myid)
               END IF
C5168          CONTINUE

               CALL MPI_Send(SlaveNodeIndex, 1, MPI_INTEGER,
     +                        0, itag, MPI_COMM_WORLD, IERR)

            END DO

5555        CONTINUE

         END IF

6666     CONTINUE
         CALL MPI_BARRIER(MPI_COMM_WORLD, IERR)

         Ncount2d=Ntau*Nmx1
         DO imy = 1, nmy1
            CALL MPI_REDUCE(DATW_R(1,1,imy), RIMAG(1,1,imy),
     +           Ncount2d, MPI_REAL, MPI_SUM, 
     +           0, MPI_COMM_WORLD, IERR)
         END DO

C6666     CONTINUE
         write(*,*) '--------Begin MPI_BARRIER--------, myid is', myid
         CALL MPI_BARRIER(MPI_COMM_WORLD, IERR)

         IF (myid .eq. 0) THEN
            IF ( ITauBreak .EQ. 0 ) THEN
               fnlen = len_trim(FN4)
               FN4(fnlen+1:fnlen+1) = char(0)
               CALL gdbfcwrite3float(FN4, RIMAG, Nmy1, Nmx1, Ntau, Iphr)
               write(*,*) ' ****begin to write CIG file of IPHR: ', Iphr
               CALL  write_cigs_regular4d(
     +               IFN03, RIMAG, Ntau, Nmx1, Nmy1,
     +               Ncdp_first, Nline_first, Nline_final, 
     +               Iphr1, Iphr2, Iphr)
               write(*,*) ' ****write CIG file done.'
            ELSE
               write(*,*) ' ****begin to write CIG file of IPHR: ', Iphr
               CALL  write_cigs_regular4d_bps(
     +               IFN03, RIMAG, Ntau, Nmx1, Nmy1,
     +               Ncdp_first, Nline_first, Nline_final, 
     +               Iphr1, Iphr2, Iphr, Ktaur)
               write(*,*) ' ****write CIG file done.'
            END IF
         END IF
         
1111  CONTINUE

      IF (MYID .EQ. 0) THEN 
         DEALLOCATE(RIMAG)
         DEALLOCATE(DATW)
      END IF

      DEALLOCATE(SENDBUF)
      DEALLOCATE(TMP2D)

      RETURN
      END

*=====================================================================*
      SUBROUTINE WXFD_EXTRAPOLATION_WXYCOMP_MPI_ACOUSTIC(
     +    DATW_W, DATW_R, XX, VVtau, TMP2D,
     +    AA_X, BB_X, CC_X, DD_X, XX_X, YY_X,
     +    AA_Y, BB_Y, CC_Y, DD_Y, XX_Y, YY_Y,
     +    ACOE, BCOE, LOOPP, ALFA,
     +    NX, Dx, NY, Dy, Ntau, Dtau, W, RPh,
     +    Dkx, Dky, Rkz_max, DIP_max, DIP_min,
     +    Nline_first, Nline_final,
     +    Ncdp_first, Ncdp_final, IFN02,
     +    Nmx0, Nmx1, Nmy0, Nmy1, Iphr, Iww, NtauExtraP,
     +    ITauBreak, Ktaur, Ktauw, FNW, FNR,
     +    IphrR1, IphrR2, myid)

!     Definition of input-parameters
      PARAMETER(PAI=3.14159265359, PAI2=2*3.14159265359)
      INTEGER   Nline_first, Nline_final, NY, Nmy0, Nmy1
      INTEGER   Ncdp_first,  Ncdp_final,  NX, Nmx0, Nmx1
      INTEGER   Ntau, LOOPP, IFN02, myid
      INTEGER   NtauExtraP, ITauBreak, Ktaur, Ktauw
      INTEGER   IphrR1, IphrR2
      REAL      ALFA, Dx, Dy, Dtau, W, RPh, Dkx, Dky, Rkz_max
      REAL      DIP_max, DIP_min, DIP

      REAL      VVtau(NX, NY), DATW_R(Ntau,Nmx1, Nmy1)
      REAL      ACOE(5), BCOE(5)

      COMPLEX   DATW_W(NX, NY), XX(NX, NY), TMP2D(NY, NX)
      COMPLEX   AA_X(NX), BB_X(NX), CC_X(NX), DD_X(NX)
      COMPLEX   AA_Y(NY), BB_Y(NY), CC_Y(NY), DD_Y(NY)
      COMPLEX   XX_X(NX), YY_X(NX), XX_Y(NY), YY_Y(NY)

      INTEGER   Itau, weflag, NtauShallow
      INTEGER   Itaus, Itaue

C     FengBo add temp file saving&reading.
      CHARACTER*256 FNW, FNR
      INTEGER       fnlen, rflag
      LOGICAL       filestat
C     FengBo write velocity, testing method. 2010-07-19.
C     CHARACTER*256 velfile

C     velfile='slice_vel.dat'
C     NtauShallow = 1000

      RPh2     = RPh*RPh
      Sita_max = 85.0
      cos_sita = cos(Sita_max*PAI2/360.0)
      cos_sita2_min = cos_sita*cos_sita

C     OPEN(66,FILE=FNW,ACCESS='DIRECT',STATUS='OLD',RECL=NX*NY*2)
      IF ( ITauBreak .EQ. 0 ) THEN
C     Expolation from the surface.
         Itaus = 1
         Itaue = NtauExtraP
      ELSE
C     Expolation from a certain depth: Ktaur.
         IF(Iphr.ge.IphrR1.and.Iphr.le.IphrR2)THEN
            INQUIRE(FILE=FNR,exist=filestat)
            IF (filestat) THEN
C              Temp wave-field exist! Expolation from depth: Ktaur.
               OPEN(68,FILE=FNR,ACCESS='DIRECT',RECL=NX*NY*2)
               Itaus = Ktaur
               Itaue = NtauExtraP
               READ(68,REC=Iww,IOSTAT=rflag)
     +         ((DATW_W(IX,IY),IX=1,NX),IY=1,NY)
               if(rflag .ne. 0 ) then
                  stop
               end if
            ELSE
C           No temp wave-field exist! Expolation from surface.
               Itaus = 1
               Itaue = NtauExtraP
            END IF
         ELSE
C        No temp wave-field exist! Expolation from surface.
            Itaus = 1
            Itaue = NtauExtraP
         END IF
      END IF

      DO 5555 Itau = Itaus, Itaue

         CALL INPUT_CURRENT_LAYER_VELOCITY(VVtau,
     +        Nline_first, Nline_final,
     +        Ncdp_first, Ncdp_final,
     +        NX, NY, Itau, IFN02,
     +        Nmx0, Nmx1, Nmy0, Nmy1)
         weflag = 0
         DO IY = 1, NY
         DO IX = 1, NX
            VVtau2    = VVtau(IX, IY)*VVtau(IX, IY)
            cos_sita2 = 1.0 - RPh2*VVtau2*0.25
            IF (cos_sita2 .lt. cos_sita2_min) THEN
               weflag = 1
               GOTO 10724
            END IF
         END DO
         END DO

10724    CONTINUE
         IF (weflag .EQ. 0) THEN
            CALL WXFD_EXTRAPOLATION_SINGLE_W(
     +           DATW_W, XX, VVtau, TMP2D,
     +           AA_X, BB_X, CC_X, DD_X, XX_X, YY_X,
     +           AA_Y, BB_Y, CC_Y, DD_Y, XX_Y, YY_Y,
     +           ACOE, BCOE, LOOPP, ALFA,
     +           NX, Dx, NY, Dy, Itau, Dtau, W, RPh, DIP_max,
     +           Dkx, Dky, Rkz_max, Nmx0, Nmx1, Nmy0, Nmy1)
   
            DO IY = 1, Nmy1
            DO IX = 1, Nmx1
               Imy = IY + Nmy0
               Imx = IX + Nmx0
               DATW_R(Itau,IX,IY) = DATW_R(Itau,IX,IY)
     +                            + REAL(DATW_W(Imx, Imy))
            END DO
            END DO
C           Write temp wave-field on depth: Ktauw.
            IF (Itau .eq. Ktauw) THEN
C              write(*,*) 'NX*NY =', NX*NY
C               CALL gdbfcwritecomplex(FN, DATW_W, NX*NY, Iphr,
C     +              Iww, filestat, 1)
               OPEN(66,FILE=FNW,ACCESS='DIRECT',
     +              STATUS='OLD',RECL=NX*NY*2)
               WRITE(66,REC=Iww)((DATW_W(IX,IY),IX=1,NX),IY=1,NY)
               CLOSE(66)
            END IF
         ELSE
            write(*,*) 'Exit! Iphr=',Iphr,' Iww =',Iww,' Itau=',Itau
            RETURN
         END IF

5555  CONTINUE

C     CLOSE(66)
      IF (filestat) THEN
         CLOSE(68)
      END IF

      RETURN
      END

*======================================================================*
      SUBROUTINE WXFD_EXTRAPOLATION_WXYCOMP_MPI_VISCO(
     +    DATW_W, DATW_R, XX, VVtau, TMP2D,
     +    AA_X, BB_X, CC_X, DD_X, XX_X, YY_X,
     +    AA_Y, BB_Y, CC_Y, DD_Y, XX_Y, YY_Y,
     +    ACOE, BCOE, LOOPP, ALFA,
     +    NX, Dx, NY, Dy, Ntau, Dtau, W, RPh,
     +    Dkx, Dky, Rkz_max, DIP_max, DIP_min,
     +    Nline_first, Nline_final,
     +    Ncdp_first, Ncdp_final, IFN02,
     +    Nmx0, Nmx1, Nmy0, Nmy1, Iphr, Iww, NtauExtraP,
     +    ITauBreak, Ktaur, Ktauw, FNW, FNR,
     +    IphrR1, IphrR2, myid)

!     Definition of input-parameters
      PARAMETER(PAI=3.14159265359, PAI2=2*3.14159265359)
      INTEGER   Nline_first, Nline_final, NY, Nmy0, Nmy1
      INTEGER   Ncdp_first,  Ncdp_final,  NX, Nmx0, Nmx1
      INTEGER   Ntau, LOOPP, IFN02, myid
      INTEGER   NtauExtraP, ITauBreak, Ktaur, Ktauw
      INTEGER   IphrR1, IphrR2
      REAL      ALFA, Dx, Dy, Dtau, W, RPh, Dkx, Dky, Rkz_max
      REAL      DIP_max, DIP_min, DIP

      REAL      VVtau(NX, NY), DATW_R(Ntau,Nmx1, Nmy1)
      REAL      ACOE(5), BCOE(5)

      COMPLEX   DATW_W(NX, NY), XX(NX, NY), TMP2D(NY, NX)
      COMPLEX   AA_X(NX), BB_X(NX), CC_X(NX), DD_X(NX)
      COMPLEX   AA_Y(NY), BB_Y(NY), CC_Y(NY), DD_Y(NY)
      COMPLEX   XX_X(NX), YY_X(NX), XX_Y(NY), YY_Y(NY)

      INTEGER   Itau, weflag, NtauShallow
      INTEGER   Itaus, Itaue

C     Visco-Acoustic Features.
      INTEGER   ierror1, ierror2
      REAL      VEL_CONSTANT, Q_CONSTANT
      REAL      Wmain
      COMPLEX, ALLOCATABLE :: Qvalue(:,:)
      COMPLEX, ALLOCATABLE :: RKW(:,:)

C     FengBo add temp file saving&reading.
      CHARACTER*256 FNW, FNR
      INTEGER       fnlen, rflag
      LOGICAL       filestat
C     FengBo write velocity, testing method. 2010-07-19.
C     CHARACTER*256 velfile

      ierror1 = 0
      ierror2 = 0
      ALLOCATE(Qvalue(NX,NY), stat=ierror1)
      ALLOCATE(RKW(NX,NY),    stat=ierror2)
      IF ((ierror1 .ne. 0) .or. (ierror2 .ne. 0)) THEN
           write(*,*) ' Not Enough Memory for 2D-buffers.'
           write(*,*) ' Visco- Program STOP Now!'
           STOP
      END IF

      VEL_CONSTANT = 4000.0
      Q_CONSTANT   = 80.0
      Wmain = 30.0

      RPh2     = RPh*RPh
      Sita_max = 85.0
      cos_sita = cos(Sita_max*PAI2/360.0)
      cos_sita2_min = cos_sita*cos_sita

C     OPEN(66,FILE=FNW,ACCESS='DIRECT',STATUS='OLD',RECL=NX*NY*2)
      IF ( ITauBreak .EQ. 0 ) THEN
C     Expolation from the surface.
         Itaus = 1
         Itaue = NtauExtraP
      ELSE
C     Expolation from a certain depth: Ktaur.
         IF(Iphr.ge.IphrR1.and.Iphr.le.IphrR2)THEN
            INQUIRE(FILE=FNR,exist=filestat)
            IF (filestat) THEN
C              Temp wave-field exist! Expolation from depth: Ktaur.
               OPEN(68,FILE=FNR,ACCESS='DIRECT',RECL=NX*NY*2)
               Itaus = Ktaur
               Itaue = NtauExtraP
               READ(68,REC=Iww,IOSTAT=rflag)
     +         ((DATW_W(IX,IY),IX=1,NX),IY=1,NY)
               if(rflag .ne. 0 ) then
                  stop
               end if
            ELSE
C           No temp wave-field exist! Expolation from surface.
               Itaus = 1
               Itaue = NtauExtraP
            END IF
         ELSE
C        No temp wave-field exist! Expolation from surface.
            Itaus = 1
            Itaue = NtauExtraP
         END IF
      END IF

      DO 5555 Itau = Itaus, Itaue

         write(*,*) 'Itau=',Itau
         CALL INPUT_CURRENT_LAYER_VELOCITY(VVtau,
     +        Nline_first, Nline_final,
     +        Ncdp_first, Ncdp_final,
     +        NX, NY, Itau, IFN02,
     +        Nmx0, Nmx1, Nmy0, Nmy1)
         weflag = 0
         DO IY = 1, NY
         DO IX = 1, NX
            VVtau2    = VVtau(IX, IY)*VVtau(IX, IY)
            cos_sita2 = 1.0 - RPh2*VVtau2*0.25
            IF (cos_sita2 .lt. cos_sita2_min) THEN
               weflag = 1
               GOTO 10724
            END IF
         END DO
         END DO

10724    CONTINUE
         IF (weflag .EQ. 0) THEN
C           CALL COMPUTE_Q_BY_V(VVtau, Qvalue, NX, NY)
            CALL COMPUTE_Q_BY_V_txg(VVtau, Qvalue, NX, NY)

            CALL WXFD_EXTRAPOLATION_SINGLE_W_VISCO(
     +           DATW_W, XX, VVtau, TMP2D,
     +           AA_X, BB_X, CC_X, DD_X, XX_X, YY_X,
     +           AA_Y, BB_Y, CC_Y, DD_Y, XX_Y, YY_Y,
     +           ACOE, BCOE, LOOPP, ALFA,
     +           NX, Dx, NY, Dy, Itau, Dtau, W, RPh, DIP_max,
     +           Dkx, Dky, Rkz_max, Nmx0, Nmx1, Nmy0, Nmy1,
     +           Qvalue, RKW, Wmain)
   
            DO IY = 1, Nmy1
            DO IX = 1, Nmx1
               Imy = IY + Nmy0
               Imx = IX + Nmx0
               DATW_R(Itau,IX,IY) = DATW_R(Itau,IX,IY)
     +                            + REAL(DATW_W(Imx, Imy))
            END DO
            END DO
C           Write temp wave-field on depth: Ktauw.
            IF (Itau .eq. Ktauw) THEN
C              write(*,*) 'NX*NY =', NX*NY
C               CALL gdbfcwritecomplex(FN, DATW_W, NX*NY, Iphr,
C     +              Iww, filestat, 1)
               OPEN(66,FILE=FNW,ACCESS='DIRECT',
     +              STATUS='OLD',RECL=NX*NY*2)
               WRITE(66,REC=Iww)((DATW_W(IX,IY),IX=1,NX),IY=1,NY)
               CLOSE(66)
            END IF
         ELSE
            write(*,*) 'Exit! Iphr=',Iphr,' Iww =',Iww,' Itau=',Itau
            RETURN
         END IF

5555  CONTINUE

C     CLOSE(66)
      IF (filestat) THEN
         CLOSE(68)
      END IF

      DEALLOCATE(Qvalue)
      DEALLOCATE(RKW)

      RETURN
      END

*======================================================================*
      SUBROUTINE WXFD_EXTRAPOLATION_WXYCOMP_MPI_VISCO_OLD(
     +    DATW_W, DATW_R, XX, VVtau, TMP2D,
     +    AA_X, BB_X, CC_X, DD_X, XX_X, YY_X,
     +    AA_Y, BB_Y, CC_Y, DD_Y, XX_Y, YY_Y,
     +    ACOE, BCOE, LOOPP, ALFA,
     +    NX, Dx, NY, Dy, Ntau, Dtau, W, RPh,
     +    Dkx, Dky, Rkz_max, DIP_max, DIP_min,
     +    Nline_first, Nline_final,
     +    Ncdp_first, Ncdp_final, IFN02,
     +    Nmx0, Nmx1, Nmy0, Nmy1, Iphr, Iww, myid)

!     Definition of input-parameters
      PARAMETER(PAI=3.14159265359, PAI2=2*3.14159265359)
      INTEGER   Nline_first, Nline_final, NY, Nmy0, Nmy1
      INTEGER   Ncdp_first,  Ncdp_final,  NX, Nmx0, Nmx1
      INTEGER   Ntau, LOOPP, IFN02, myid
      REAL      ALFA, Dx, Dy, Dtau, W, RPh, Dkx, Dky, Rkz_max
      REAL      DIP_max, DIP_min, DIP

      REAL      VVtau(NX, NY), DATW_R(Ntau,Nmx1, Nmy1)
      REAL      ACOE(5), BCOE(5)

      COMPLEX   DATW_W(NX, NY), XX(NX, NY), TMP2D(NY, NX)
      COMPLEX   AA_X(NX), BB_X(NX), CC_X(NX), DD_X(NX)
      COMPLEX   AA_Y(NY), BB_Y(NY), CC_Y(NY), DD_Y(NY)
      COMPLEX   XX_X(NX), YY_X(NX), XX_Y(NY), YY_Y(NY)

      INTEGER   Itau, weflag, NtauShallow

C     Visco-Acoustic Features.
      INTEGER   ierror1, ierror2
      REAL      VEL_CONSTANT, Q_CONSTANT
      REAL      Wmain
      COMPLEX, ALLOCATABLE :: Qvalue(:,:)
      COMPLEX, ALLOCATABLE :: RKW(:,:)

      ierror1 = 0
      ierror2 = 0
      ALLOCATE(Qvalue(NX,NY), stat=ierror1)
      ALLOCATE(RKW(NX,NY),    stat=ierror2)
      IF ((ierror1 .ne. 0) .or. (ierror2 .ne. 0)) THEN
           write(*,*) ' Not Enough Memory for 2D-buffers.'
           write(*,*) ' Visco- Program STOP Now!'
           STOP
      END IF

      VEL_CONSTANT = 4000.0
      Q_CONSTANT   = 80.0
      Wmain = 30.0

C     NtauShallow = 1000
      RPh2     = RPh*RPh
      Sita_max = 85.0
      cos_sita = cos(Sita_max*PAI2/360.0)
      cos_sita2_min = cos_sita*cos_sita

      DO 5555 Itau = 1, Ntau

         CALL INPUT_CURRENT_LAYER_VELOCITY(VVtau,
     +        Nline_first, Nline_final,
     +        Ncdp_first, Ncdp_final,
     +        NX, NY, Itau, IFN02,
     +        Nmx0, Nmx1, Nmy0, Nmy1)
C        IF ( Itau .LT. NtauShallow ) THEN
C           DIP = DIP_min + (DIP_max-DIP_min)*(Itau-1)/(NtauShallow-1)
C           Rkz_max = cos(PAI2*DIP/360.0)
C        ELSE
C           Rkz_max = cos(PAI2*DIP_max/360.0)
C        END IF
         weflag = 0
         DO IY = 1, NY
         DO IX = 1, NX
            VVtau2    = VVtau(IX, IY)*VVtau(IX, IY)
            cos_sita2 = 1.0 - RPh2*VVtau2*0.25
            IF (cos_sita2 .lt. cos_sita2_min) THEN
               weflag = 1
               GOTO 10724
            END IF
         END DO
         END DO

10724    CONTINUE
         IF (weflag .EQ. 0) THEN
!            CALL COMPUTE_Q_BY_V(VVtau, Qvalue, NX, NY)
            CALL COMPUTE_Q_BY_V_txg(VVtau, Qvalue, NX, NY)

            CALL WXFD_EXTRAPOLATION_SINGLE_W_VISCO(
     +           DATW_W, XX, VVtau, TMP2D,
     +           AA_X, BB_X, CC_X, DD_X, XX_X, YY_X,
     +           AA_Y, BB_Y, CC_Y, DD_Y, XX_Y, YY_Y,
     +           ACOE, BCOE, LOOPP, ALFA,
     +           NX, Dx, NY, Dy, Itau, Dtau, W, RPh, DIP_max,
     +           Dkx, Dky, Rkz_max, Nmx0, Nmx1, Nmy0, Nmy1,
     +           Qvalue, RKW, Wmain)

            DO IY = 1, Nmy1
            DO IX = 1, Nmx1
               Imy = IY + Nmy0
               Imx = IX + Nmx0
               DATW_R(Itau,IX,IY) = DATW_R(Itau,IX,IY)
     +                            + REAL(DATW_W(Imx, Imy))
            END DO
            END DO
         ELSE
            write(*,*) 'Exit! Iphr=',Iphr,' Iww =',Iww,' Itau=',Itau
            RETURN
         END IF

5555  CONTINUE

      DEALLOCATE(Qvalue)
      DEALLOCATE(RKW)

      RETURN
      END

*======================================================================*
      SUBROUTINE write_cigs_regular4d(
     +    IFN03, RIMAG, Ntau, Nmx1, Nmy1,
     +    Ncdp_first, Nline_first, Nline_final,
     +    Iphr1, Iphr2, Iphr)

      INTEGER IFN03, Ntau, Nmx1, Nmy1
      INTEGER Ncdp_first, Nline_first, Nline_final
      INTEGER Iphr1, Iphr2, Iphr

      INTEGER ipar240(20)
      INTEGER icdp, iline, iflag
      REAL    RIMAG(Ntau, Nmx1, Nmy1)

      DO itr = 1, 20
         ipar240(itr) = 0
      END DO

      ipar240(1) = 1                 ! Start 1st Value.
      ipar240(2) = Ncdp_first        ! Start 2nd Value.
      ipar240(3) = Nline_first       ! Start 3rd Value.
      ipar240(4) = Nline_final       ! Final 3rd Value.
      ipar240(5) = Iphr              ! Start 4th Value.

      iflag = 3
!     CALL gdbfcwriteregulardata(IFN03, ipar240, iflag, RIMAG)
      CALL gdbfcwriteregularfile(IFN03, ipar240, iflag, RIMAG)
      IF (iflag .ne. 0) THEN
          write(*,*) ' gdbfcwriteregularfile() failed.'
          write(*,*) ' iflag = ', iflag
          write(*,*) ' Iphr  = ', Iphr
          STOP
      END IF

      RETURN
      END

*=====================================================================*
      SUBROUTINE write_cigs_regular4d_bps(
     +    IFN03, RIMAG, Ntau, Nmx1, Nmy1,
     +    Ncdp_first, Nline_first, Nline_final,
     +    Iphr1, Iphr2, Iphr, Ktaur)
      INTEGER IFN03, Ntau, Nmx1, Nmy1
      INTEGER Ncdp_first, Nline_first, Nline_final
      INTEGER Iphr1, Iphr2, Iphr, Ktaur

      INTEGER ipar240(20)
      INTEGER icdp, iline, iflag
      REAL    RIMAG(Ntau, Nmx1, Nmy1)
      INTEGER imx,imy

      DO itr = 1, 20
         ipar240(itr) = 0
      END DO

      ipar240(1) = Ktaur  ! Start 1st Value.
      ipar240(2) = Ntau   ! End   1st Value.

      DO imy=1,Nmy1
         DO imx=1,Nmx1
          ipar240(3) = imx+Ncdp_first-1    ! Start 2nd Value.
          ipar240(4) = imy+Nline_first-1   ! Start 3rd Value.
          ipar240(5) = Iphr                ! Start 4th Value.
          ipar240(6) = 1                   ! Start 5th Value.

!         write(*,*)'imx',imx,'imy',imy
          iflag = 1
          CALL gdbfcwriteregularfile(IFN03, ipar240, iflag,
     +                               RIMAG(Ktaur,imx,imy))
          IF (iflag .ne. 0) THEN
           write(*,*) ' gdbfcwriteregularfile() failed.'
           write(*,*) ' iflag = ', iflag
           write(*,*) ' Iphr  = ', Iphr
          STOP
          END IF
         ENDDO
      ENDDO

      RETURN
      END

*=====================================================================*
      SUBROUTINE READ_PLANEWAVE_DATA_AND_FFT_Profile2d(
     +      Trace, DATW, CTrace,
     +     Nmx, Nmy, NW1, NW2, NW3, NW4, NW, LT, DW, 
     +     Ncdp_first, Ncdp_final, Ncdp_step,
     +     Nline_first, Nline_final, Nline_step, 
     +     Iphr, IFN01, Nmx0, Nmx1, Nmy0, Nmy1)

      INTEGER NW1, NW, IFN01,LT
      INTEGER Ncdp_first, Ncdp_final, Ncdp_step
      INTEGER Nline_first, Nline_final, Nline_step
      REAL    DW

      REAL    Trace(LT)
      COMPLEX CTrace(LT)
      COMPLEX DATW(Nmx, Nmy, NW)
      INTEGER IPAR240(25)

      REAL,   ALLOCATABLE :: Pdat2d(:,:)
      COMPLEX,ALLOCATABLE :: CTMPTRACE2d(:,:)
      INTEGER ierr, ncdp

      ncdp  = Ncdp_final-Ncdp_first+1
      ierr1 = 0
      ierr2 = 0
      ALLOCATE(CTMPTRACE2d(LT,ncdp),stat=ierr1)
      ALLOCATE(Pdat2d(LT,ncdp),     stat=ierr2)
      IF ((ierr1 .ne. 0).or.(ierr2 .ne. 0)) THEN
           write(*,*) ' Not Enough Memory!'
           STOP
      END IF

      DO iw  = 1, NW
      DO imy = 1, Nmy
      DO imx = 1, Nmx
         DATW(imx,imy,iw) = CMPLX(0.0,0.0)
      END DO
      END DO
      END DO

      IPAR240(1) = 1             ! first grid number of 1st Dimension.
      IPAR240(2) = LT            ! last  grid number of 1st Dimension.
      IPAR240(3) = Ncdp_first    ! location of 2nd Dimension.
      IPAR240(5) = Iphr          ! location of 4th Dimension.
      IPAR240(6) = 1             ! location of 5th Dimension.

      DO 2010 Iline = Nline_first, Nline_final

         IPAR240(4) = Iline      ! location of 3rd Dimension.
         Imy = Iline - Nline_first + 1 + Nmy0
C        write(*,*) 'Iline = ', Iline
         if( mod(Iline,10) .eq. 0 ) write(*,*) 'Iline= ', Iline

         IFLAG = 2
         CALL gdbfcreadregulardata(IFN01, IPAR240, IFLAG, Pdat2d)
         IF (IFLAG.NE.0) THEN
            WRITE(*,*) 'gdbfcreadregulardata() failed!'
            WRITE(*,*) 'IFlAG = ',IFlAG
            WRITE(*,*) 'PROGRAM EXIT.'
            STOP
         END IF
C        write(*,*) 'reading Iline ', Iline, 'done.'
C        goto 2010

*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*
C        FengBo add code, scaning for nan-value and zeroing it.
C        2010-07-17, 21:36.
         do iix=1,ncdp
         do iit=1,lT
            if(isnan(Pdat2d(iit,iix))) then
               Pdat2d(iit,iix) = 0
               WRITE(*,*) 'Iline=',Iline,'iix=',iix,'iit=',iit
            endif
         enddo
         enddo
*<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*

C$OMP PARALLEL PRIVATE(Imx,IX,IT)
C$OMP DO
         DO Icdp = Ncdp_first, Ncdp_final

            Imx = Icdp - Ncdp_first + 1 + Nmx0
            IX  = Icdp - Ncdp_first + 1

            DO 1010 IT = 1, LT
               CTMPTRACE2d(IT,IX)=CMPLX(Pdat2d(IT,IX), 0.0)
1010        CONTINUE

            CALL FFT1D(CTMPTRACE2d(1,IX), LT, 1)

            DO 1011 Iw = 1, NW
               DATW(Imx, Imy, Iw)=CTMPTRACE2d(NW1+Iw-1,IX)
1011        CONTINUE

         END DO
C$OMP END DO
C$OMP END PARALLEL

2010  CONTINUE

      DEALLOCATE(Pdat2d)
      DEALLOCATE(CTMPTRACE2d)

      RETURN
      END

*=====================================================================*
      SUBROUTINE READ_PLANEWAVE_DATA_AND_FFT(Trace, DATW, CTrace,
     +     Nmx, Nmy, NW1, NW2, NW3, NW4, NW, LT, DW, 
     +     Ncdp_first, Ncdp_final, Ncdp_step,
     +     Nline_first, Nline_final, Nline_step, 
     +     Iphr, IFN01, Nmx0, Nmx1, Nmy0, Nmy1)

      INTEGER NW1, NW, IFN01
      INTEGER Ncdp_first, Ncdp_final, Ncdp_step
      INTEGER Nline_first, Nline_final, Nline_step
      REAL    DW

      REAL    Trace(LT)
      COMPLEX CTrace(LT)
      COMPLEX DATW(Nmx, Nmy, NW)
      INTEGER IPAR240(25)

      DO iw  = 1, NW
      DO imy = 1, Nmy
      DO imx = 1, Nmx
         DATW(imx,imy,iw) = CMPLX(0.0,0.0)
      END DO
      END DO
      END DO

      IPAR240(1) = 1             ! first grid number of 1st Dimension.
      IPAR240(2) = LT            ! last  grid number of 1st Dimension.
      IPAR240(5) = Iphr          ! location of 4th Dimension.
      IPAR240(6) = 1             ! location of 5th Dimension.

      IPAR240(4) = Nline_first   ! location of 3rd Dimension.
      DO Iline = Nline_first, Nline_final

         IPAR240(3) = Ncdp_first    ! location of 2nd Dimension.
         Imy = Iline - Nline_first + 1 + Nmy0

         DO Icdp = Ncdp_first, Ncdp_final

            Imx = Icdp - Ncdp_first + 1 + Nmx0
            IFLAG = 1
            CALL gdbfcreadregulardata(IFN01, IPAR240,
     +               IFLAG, Trace)
            IF (IFLAG.NE.0) THEN
               WRITE(*,*) 'gdbfcreadregulardata() failed!'
               WRITE(*,*) 'IFlAG = ',IFlAG
               WRITE(*,*) 'PROGRAM EXIT.'
               STOP
            END IF

C           CALL HAMMING_WINDOW_TRIM(Trace, 1, 1, LT-50, LT, LT)
            DO IT=1, LT
               CTrace(IT)=CMPLX(Trace(IT), 0.0) 
            END DO

            CALL FFT1D(CTrace, LT, 1)
C            CALL HAMMING_WINDOW(CTrace, NW1, NW2, NW3, NW4, LT)

            DO Iw=1, NW
               DATW(Imx, Imy, Iw)=CTrace(NW1+Iw-1)
            END DO

            IPAR240(3) = IPAR240(3) + Ncdp_step

         END DO

         IPAR240(4) = IPAR240(4) + Nline_step

      END DO

      RETURN
      END

*=====================================================================*
      SUBROUTINE make_order(DIP_max, Acoe, Bcoe, LOOPP)
      REAL      ACOE(5), BCOE(5)
      REAL      DIP_max
      INTEGER   LOOPP

      DO I=1,5
        Acoe(I)=0.0
        Bcoe(I)=0.0
      END DO

      IF (DIP_max.LE.45) THEN
        Acoe(1)=0.5
        Bcoe(1)=0.25
        LOOPP=1
      ELSE IF (DIP_max.LE.60) THEN
        Acoe(1)=0.4761
        Bcoe(1)=0.3767
        LOOPP=1
      ELSE IF (DIP_max.LE.70) THEN
        Acoe(1)=0.449565816
        Bcoe(1)=0.426137316
        LOOPP=1
      ELSE IF (DIP_max.LE.79) THEN
        Acoe(1)=0.4575
        Bcoe(1)=0.4575
        LOOPP=1
      ELSE IF (DIP_max.LE.80) THEN
        Acoe(1)=0.457289566
        Bcoe(1)=0.222691983
        Acoe(2)=0.040315157
        Bcoe(2)=0.873981642
        LOOPP=2
      ELSE IF (DIP_max.LE.87) THEN
        Acoe(1)=0.414236605
        Bcoe(1)=0.150843924
        Acoe(2)=0.081312882
        Bcoe(2)=0.744418059
        Acoe(3)=0.00421042
        Bcoe(3)=0.972926132
        LOOPP=3
      ELSE IF (DIP_max.LE.89) THEN
        Acoe(1)=0.367013245
        Bcoe(1)=0.105756624
        Acoe(2)=0.117592008
        Bcoe(2)=0.614520676
        Acoe(3)=0.014853510
        Bcoe(3)=0.919432661
        Acoe(4)=0.000523275
        Bcoe(4)=0.994065088
        LOOPP=4
      ELSE IF (DIP_max.LE.90) THEN
        Acoe(1)=0.318013812
        Bcoe(1)=0.073588213
        Acoe(2)=0.143798076
        Bcoe(2)=0.483340757
        Acoe(3)=0.033860918
        Bcoe(3)=0.824918565
        Acoe(4)=0.004172967
        Bcoe(4)=0.964827992
        Acoe(5)=0.000153427
        Bcoe(5)=0.997370236
        LOOPP=5
      END IF

      RETURN
      END

*======================================================================*
      SUBROUTINE INPUT_CURRENT_LAYER_VELOCITY(VVtau,
     +     Nline_first, Nline_final,
     +     Ncdp_first, Ncdp_final,
     +     Nmx, Nmy, Itau, IFN02,
     +     Nmx0, Nmx1, Nmy0, Nmy1)

      INTEGER Nvmx, Nvmy, Nmx, Nmy
      INTEGER Ncdp_first, Nvcdp_first, Nline_first, Nvline_first
      INTEGER Ncdp_final, Nvcdp_final, Nline_final, Nvline_final

      REAL    VVtau(Nmx, Nmy)

      DO imy = 1, Nmy
      DO imx = 1, Nmx
         VVtau(imx,imy) = 0
      END DO
      END DO

      N4=1
      DO Iline = Nline_first, Nline_final
         imy = Iline - Nline_first + 1 + Nmy0
         imx = Ncdp_first - Ncdp_first + Nmx0 + 1
         IFLAG = 1
         CALL GDBFCREADCUTSLICEDATA(IFN02, N4, Itau, Iline,
     +             Ncdp_first, Ncdp_final, IFLAG, VVtau(imx,imy))
         IF (IFLAG .NE. 0) THEN
            WRITE(*,*) 'GDBFCREADCUTSLICEDATA ERROR.'
            RETURN
         END IF
      END DO

C     Filling the Velocity Field Slice.
      DO imy = Nmy0+1, Nmy0+Nmy1
         DO imx = 1, Nmx0
            VVtau(imx,imy) = VVtau(Nmx0+1,imy)
         END DO
         DO imx = Nmx0+Nmx1+1, Nmx
            VVtau(imx,imy) = VVtau(Nmx0+Nmx1,imy)
         END DO
      END DO

      DO imx = 1, Nmx
         DO imy = 1, Nmy0
            VVtau(imx,imy) = VVtau(imx,Nmy0+1)
         END DO
         DO imy = Nmy0+Nmy1+1, Nmy
            VVtau(imx,imy) = VVtau(imx,Nmy0+Nmy1)
         END DO
      END DO

      RETURN
      END

*======================================================================*
      SUBROUTINE COMPUTE_Q_BY_V(VVtau, Qvalue, NX, NY)

      REAL      VVtau(NX, NY), Qvalue(NX, NY)
      INTEGER   NX, NY, IX, IY

      REAL      scale_qv

      scale_qv = 1/100.0
      DO IY = 1, NY
      DO IX = 1, NX
         Qvalue(IX,IY) = VVtau(IX,IY)*scale_qv
      END DO
      END DO

      RETURN
      END
*======================================================================*
      SUBROUTINE COMPUTE_Q_BY_V_txg(VVtau, Qvalue, NX, NY)

      REAL      VVtau(NX, NY), Qvalue(NX, NY)
      INTEGER   NX, NY, IX, IY

      REAL      scale_qv

      scale_qv = 1/50.0
      DO IY = 1, NY
      DO IX = 1, NX
C         Qvalue(IX,IY) = VVtau(IX,IY)**2.2*1.0E-6*3.516       ! Li's equation.
C         delt_v    = VVtau(IX,IY) - 1750
C         ndeltv    = delt_v/20.
C         Qvalue(IX,IY) = 25 + ndeltv
C         Qvalue(IX,IY) = VVtau(IX,IY)*scale_qv
          delt_v    = ABS(VVtau(IX,IY)-1500.)
          Qvalue(IX,IY) = 20 + delt_v/40
      END DO
      END DO

      RETURN
      END
*======================================================================*
      SUBROUTINE HAMMING_WINDOW(TT_TEMP, NW1, NW2, NW3, NW4, LT)

      PARAMETER(PAI=3.14159165359)
      INTEGER  NW1, NW2, NW3, NW4
      COMPLEX  TT_TEMP(LT)

      DO IW=1, LT/2+1

         IF (IW.GE.NW1.AND.IW.LE.NW2) THEN
           Hammingw=0.54+0.46*cos(PAI*(Iw-NW1)/(NW2-NW1)-PAI)
           TT_TEMP(Iw)=TT_TEMP(Iw)*Hammingw
         ELSE IF (IW.GE.NW3.AND.IW.LE.NW4) THEN
           Hammingw=0.54+0.46*cos(pai*(NW3-Iw)/(NW4-NW3))
           TT_TEMP(Iw)=TT_TEMP(Iw)*Hammingw
         ELSE IF (IW.GT.NW4.OR.IW.LT.NW1) THEN
          TT_TEMP(IW)=0.0
         END IF

       END DO

       RETURN
       END

*=====================================================================*
      SUBROUTINE WXFD_EXTRAPOLATION_SINGLE_W(
     +    DATW_W, XX, VVtau, TMP2D,
     +    AA_X, BB_X, CC_X, DD_X, XX_X, YY_X,
     +    AA_Y, BB_Y, CC_Y, DD_Y, XX_Y, YY_Y,
     +    ACOE, BCOE, LOOPP, ALFA,
     +    NX, Dx, NY, Dy, Itau, Dtau, W, RPh, DIP_max,
     +    Dkx, Dky, Rkz_max, Nmx0, Nmx1, Nmy0, Nmy1)

!     DO W-X-Y Wave Field Extrapolation with Li's compensation.
      PARAMETER(PAI=3.14159265359, PAI2=2*3.14159265359, Wflag=-1.0)

!     Definition of input-parameters
      INTEGER   NY, Nmy0, Nmy1, NX, Nmx0, Nmx1
      INTEGER   LOOPP, Itau
      REAL      ALFA, Dx, Dy, Dtau, W, RPh, Dkx, Dky, Rkz_max, DIP_max

      REAL      VVtau(NX, NY)
      REAL      ACOE(5), BCOE(5)

      COMPLEX   DATW_W(NX, NY), XX(NX, NY), TMP2D(NY, NX)
      COMPLEX   AA_X(NX), BB_X(NX), CC_X(NX), DD_X(NX)
      COMPLEX   AA_Y(NY), BB_Y(NY), CC_Y(NY), DD_Y(NY)
      COMPLEX   XX_X(NX), YY_X(NX), XX_Y(NY), YY_Y(NY)

!     Definition of local-parameters
      REAL      DX2, DY2, W2, VVtau2, VZ
      REAL      cos_sita, cos_sita2_min, cos_sita2, Sita_max
      COMPLEX   PHASE

      RPh2 = RPh*RPh
      DX2  = DX*DX
      DY2  = DY*DY
      W2   = W*W
C     Sita_max = 88.0
C     Sita_max = DIP_max
      Sita_max = 85.0
      cos_sita = cos(Sita_max*PAI2/360.0)
      cos_sita2_min = cos_sita*cos_sita

!     Solve  equation (a) by Phase-Shift.
      DO IY=1, NY
         DO IX=1, NX
            VVtau2         = VVtau(IX, IY)*VVtau(IX, IY)
            cos_sita2      = 1.0 - RPh2*VVtau2*0.25
C           IF (cos_sita2 .lt. cos_sita2_min) THEN
C              weflag = 1
C              RETURN
C           END IF
            cos_sita       = SQRT(cos_sita2)
            Tshift         = Wflag*W*Dtau*cos_sita
            PHASE          = CMPLX(0.0, Tshift)
            DATW_W(IX, IY) = DATW_W(IX, IY)*CEXP(PHASE)
         END DO
      END DO

*========= Solving the equation along the X-direction.=================*

      DO 987 IY = 1, NY

         DO 808 IL=1, LOOPP

         BETAX1 = BCOE(IL)/(4.0*W2*DX2)
         BETAX2 = Wflag*Acoe(IL)*Dtau/(8.0*DX2*W)

         DO IX=1, NX
            VVtau2   = VVtau(IX, IY)*VVtau(IX, IY)
            cos_sita2= 1.0-RPh2*VVtau2*0.25
            cos_sita = SQRT(cos_sita2)
            betax11  = betax1*VVtau2/Cos_Sita2 + ALFA
            betax22  = betax2*VVtau2/cos_sita

            AA_X(IX) = cmplx(betax11,       +betax22)
            BB_X(IX) = cmplx(1-2*betax11, -2*betax22)
            CC_X(IX) = cmplx(betax11,       -betax22)
            DD_X(IX) = cmplx(1-2*betax11, +2*betax22)
         END DO

         YY_X(1) = AA_X(1)*DATW_W(2, IY) + BB_X(1)*DATW_W(1, IY)
         DO IX=2, NX-1
            YY_X(IX) = AA_X(IX)*(DATW_W(IX+1, IY)+DATW_W(IX-1, IY))
     +               + BB_X(IX)*DATW_W(IX, IY)
         END DO
         YY_X(NX) = AA_X(NX)*DATW_W(NX-1, IY) + BB_X(NX)*DATW_W(NX, IY)

         CALL COMPLEX_TRIDAG(CC_X, DD_X, YY_X, XX_X, NX)

         DO IX=1, NX
            DATW_W(IX, IY)=XX_X(IX)
         END DO

C        For the attenuation boundary. Left  side.
         DO IX = 1, Nmx0
            IXX0 = Nmx0 - IX + 1
            DATW_W(IX, IY) = XX_X(IX)*(1.0+cos(PAI*IXX0/Nmx0))*0.5
         END DO

C        For the attenuation boundary. Right side.
         DO IX = Nmx1+Nmx0+1, NX
            IXX0 = IX - Nmx1 - Nmx0
            DATW_W(IX, IY) = XX_X(IX)*(1.0+cos(PAI*IXX0/Nmx0))*0.5
         END DO

808      CONTINUE

987   CONTINUE

*========= Solving the equation along the Y-direction.=================*

      DO 988 IX = 1, NX

         DO 909 IL=1, LOOPP

         BETAY1 = BCOE(IL)/(4.0*W2*DY2)
         BETAY2 = Wflag*ACOE(IL)*Dtau/(8.0*W*DY2)

         DO IY=1, NY
            VVtau2   = VVtau(IX, IY)*VVtau(IX, IY)
            cos_sita2= 1.0-RPh2*VVtau2*0.25
            cos_sita = SQRT(cos_sita2)
            betay11  = betay1*VVtau2/Cos_Sita2 + ALFA
            betay22  = betay2*VVtau2/cos_sita

            AA_Y(IY) = cmplx(betay11,       +betay22)
            BB_Y(IY) = cmplx(1-2*betay11, -2*betay22)
            CC_Y(IY) = cmplx(betay11,       -betay22)
            DD_Y(IY) = cmplx(1-2*betay11, +2*betay22)
         END DO

         YY_Y(1) = AA_Y(1)*DATW_W(IX, 2) + BB_Y(1)*DATW_W(IX, 1)
         DO IY=2, NY-1
            YY_Y(IY) = AA_Y(IY)*(DATW_W(IX, IY+1)+DATW_W(IX, IY-1))
     +               + BB_Y(IY)*DATW_W(IX, IY)
         END DO
         YY_Y(NY) = AA_Y(NY)*DATW_W(IX, NY-1) + BB_Y(NY)*DATW_W(IX, NY)

         CALL COMPLEX_TRIDAG(CC_Y, DD_Y, YY_Y, XX_Y, NY)

         DO Iy=1, NY
            DATW_W(IX, IY)=XX_Y(IY)
         END DO

C        For the attenuation boundary. Left  side.
         DO IY = 1, Nmy0
            IYY0 = Nmy0 - IY + 1
            DATW_W(IX, IY) = XX_Y(IY)*(1.0+cos(PAI*IYY0/Nmx0))*0.5
         END DO

C        For the attenuation boundary. Right side.
         DO IY = Nmy1+Nmy0+1, NY
            IYY0 = IY - Nmy1 - Nmy0
            DATW_W(IX, IY) = XX_Y(IY)*(1.0+cos(PAI*IYY0/Nmx0))*0.5
         END DO

909      CONTINUE

988   CONTINUE

!     Compensate X-Y splitting error.
      VZ = 0
      DO Iy=1, NY
      DO Ix=1, NX
         VZ = VZ + VVtau(Ix,Iy)
      END DO
      END DO
      VZ = VZ/(NX*NY)
      cosalfa = sqrt(1.0 - 0.25*VZ*VZ*RPh2)
      CALL COMPENSATE_PSFD(DATW_W, XX, XX_X, YY_Y, TMP2D,
     +     NX, Dkx, Dx, NY, Dky, Dy,
     +     Dtau, VZ, W, ALFA, Rkz_max, cosalfa)

      RETURN
      END

*======================================================================*
      SUBROUTINE COMPENSATE_PSFD(DATW_W, XX, XX_X, YY_Y, TMP2D,
     +     KxNUM, Dkx, Dx, KYNUM, Dky, Dy, 
     +     Dtau, VZ, W, ALFA, Rkz_max, cosalfa)
      PARAMETER(Wflag=-1.0, ACOE=0.4575, BCOE=0.4575)

      INTEGER  KxNUM, KYNUM
      REAL     Dkx,   Dky, Dx, Dy, Dtau, VZ
      REAL     ALFA
      REAL     AA, BB, Rkz_max

      COMPLEX  XX(KxNUM,KyNUM), DATW_W(KxNUM,KyNUM), TMP2D(KyNUM,KxNUM)
      COMPLEX  XX_X(KxNUM), YY_Y(KyNUM)
      COMPLEX  PHASE1, PHASEX, PHASEY, PHASE1XY

      REAL     a0, a00, tmpa, tmpa2
      REAL     Dx2, Dy2

      DO Iky=1, KyNUM
         DO Ikx=1, KxNUM
            XX(Ikx, Iky)=DATW_W(Ikx, Iky)
         END DO
      END DO
C     CALL FFT2D(XX, XX_X, YY_Y, KxNUM, KyNUM, +1)
      CALL FFT2D_OMP(XX, TMP2D,KxNUM, KyNUM, +1)

      a0  = VZ*VZ*0.25/(W*W*cosalfa*cosalfa)
      a00 = W*cosalfa*Dtau
      Dx2 = Dx*Dx
      Dy2 = Dy*Dy

      DO 6660 Iky  = 1, KyNUM/2 + 1

         Jky    = KyNUM + 2 - Iky
         Rky    = (Iky-1)*Dky
         Rky2   = Rky*Rky                           !Ky^2
         Rkyy2  = 2.0*(1.0-cos(Rky*Dy))/Dy2         !Ky^^2

         CYY    = ACOE*a0*Rkyy2/(1.0-(ALFA*Dy2+BCOE*a0)*Rkyy2)
         PHASEY = CEXP(CMPLX(0.0,Wflag*a00*CYY))

         DO 6661 Ikx  = 1, KxNUM/2 + 1

            Jkx    = KxNUM + 2 - Ikx
            Rkx    = (Ikx-1)*Dkx
            Rkx2   = Rkx*Rkx                        !Kx^2
            Rkxx2  = 2.0*(1.0-cos(Rkx*Dx))/Dx2      !Kx^^2

            tmpa2  = 1.0 - a0*(Rky2+Rkx2)
            IF (tmpa2 .le. 0.0) THEN
               XX(Ikx, Iky) = 0.0
               IF (Ikx.GT.1)              XX(Jkx,Iky)=0.0
               IF (Iky.GT.1)              XX(Ikx,Jky)=0.0
               IF (Iky.GT.1.AND.Ikx.GT.1) XX(Jkx,Jky)=0.0
               goto 6661
            END IF

            tmpa   = sqrt(tmpa2)
            IF (tmpa .le. Rkz_max) THEN
               XX(Ikx, Iky) = 0.0
               IF (Ikx.GT.1)              XX(Jkx,Iky)=0.0
               IF (Iky.GT.1)              XX(Ikx,Jky)=0.0
               IF (Iky.GT.1.AND.Ikx.GT.1) XX(Jkx,Jky)=0.0
               goto 6661
            END IF

            AXY    = tmpa - 1.0
            PHASE1 = CEXP(CMPLX(0.0,Wflag*a00*AXY))

            BXX    = ACOE*a0*Rkxx2/(1.0-(ALFA*Dx2+BCOE*a0)*Rkxx2)
            PHASEX = CEXP(CMPLX(0.0,Wflag*a00*BXX))

            PHASE1XY     = PHASE1*PHASEX*PHASEY
            XX(Ikx, Iky) = XX(Ikx, Iky)*PHASE1XY

            IF (Ikx.GT.1)              XX(Jkx,Iky)=XX(Jkx,Iky)*PHASE1XY
            IF (Iky.GT.1)              XX(Ikx,Jky)=XX(Ikx,Jky)*PHASE1XY
            IF (Iky.GT.1.AND.Ikx.GT.1) XX(Jkx,Jky)=XX(Jkx,Jky)*PHASE1XY

6661     CONTINUE

6660  CONTINUE

C     CALL FFT2D(XX, XX_X, YY_Y, KxNUM, KyNUM, -1)
      CALL FFT2D_OMP(XX, TMP2D,KxNUM, KyNUM, -1)

      DO Iky=1, KyNUM
         DO Ikx=1, KxNUM
            DATW_W(Ikx, Iky)=XX(Ikx, Iky)
         END DO
      END DO

      RETURN
      END

*======================================================================*
      SUBROUTINE WXFD_EXTRAPOLATION_SINGLE_W_VISCO(
     +    DATW_W, XX, VVtau, TMP2D,
     +    AA_X, BB_X, CC_X, DD_X, XX_X, YY_X,
     +    AA_Y, BB_Y, CC_Y, DD_Y, XX_Y, YY_Y,
     +    ACOE, BCOE, LOOPP, ALFA,
     +    NX, Dx, NY, Dy, Itau, Dtau, W, RPh, DIP_max,
     +    Dkx, Dky, Rkz_max, Nmx0, Nmx1, Nmy0, Nmy1,
     +    Qvalue, RKW, Wmain)

!     DO W-X-Y Wave Field Extrapolation with Li's compensation.
      PARAMETER(PAI=3.14159265359, PAI2=2*3.14159265359, Wflag=-1.0)
      PARAMETER(Pai21=1.0/PAI2)

!     Definition of input-parameters
      INTEGER   NY, Nmy0, Nmy1, NX, Nmx0, Nmx1
      INTEGER   LOOPP, Itau
      REAL      ALFA, Dx, Dy, Dtau, W, RPh, Dkx, Dky, Rkz_max, DIP_max

      REAL      VVtau(NX, NY), Qvalue(NX, NY)
      REAL      ACOE(5), BCOE(5)

      COMPLEX   DATW_W(NX, NY), XX(NX, NY), TMP2D(NY, NX)
      COMPLEX   AA_X(NX), BB_X(NX), CC_X(NX), DD_X(NX)
      COMPLEX   AA_Y(NY), BB_Y(NY), CC_Y(NY), DD_Y(NY)
      COMPLEX   XX_X(NX), YY_X(NX), XX_Y(NY), YY_Y(NY)

!     Definition of local-parameters
      REAL      DX2, DY2, W2, VVtau2, VZ, QZ
      REAL      cos_sita, cos_sita2_min, cos_sita2, Sita_max
      COMPLEX   PHASE

!     For Visco-Acoustic Media.
      REAL      Wr, Wmain, Gama, Tmp0, V_dispersion, Coe_Attenuation
      COMPLEX   RKW(NX, NY)
      COMPLEX   PHASE1, PHASE2, CKW

      Wr   = Wmain
      RPh2 = RPh*RPh
      DX2  = DX*DX
      DY2  = DY*DY
      W2   = W*W
C     Sita_max = 88.0
C     Sita_max = DIP_max
      Sita_max = 85.0
      cos_sita = cos(Sita_max*PAI2/360.0)
      cos_sita2_min = cos_sita*cos_sita

C      write(*,*) 'Phase-Shift begin, Wmain=',Wmain,'HZ'
!     Solve  equation (a) by Phase-Shift.
      DO IY=1, NY
         DO IX=1, NX

            Gama           = Pai21*Atan(0.5/Qvalue(IX, IY))
            Tmp0           = (W/Wr)**Gama
            V_dispersion   = VVtau(IX, IY)*Tmp0
            Coe_Attenuation= tan(0.5*pai*Gama)*W*TMP0/VVtau(IX, IY)
            RKW(IX, IY)    = CMPLX(W/V_dispersion, Coe_Attenuation)
            VVtau2         = VVtau(IX, IY)*VVtau(IX, IY)

            cos_sita2      = 1.0 - RPh2*VVtau2*0.25
C           IF (cos_sita2 .lt. cos_sita2_min) THEN
C              weflag = 1
C              RETURN
C           END IF
            cos_sita       = SQRT(cos_sita2)

            Tmp1           = VVtau(IX, IY)*cos_sita*Dtau
            RCK_tau        = Wflag*REAL(RKW(IX, IY))*Tmp1
            PHASE1         = CMPLX(0.0, RCK_tau)

            Coe_Attenua_Tau= AIMAG(RKW(IX, IY))*Tmp1
            PHASE2         = CMPLX(Coe_Attenua_Tau, 0.0)

            DATW_W(IX, IY) = DATW_W(IX, IY)*CEXP(PHASE1)*CEXP(PHASE2)

         END DO
      END DO
C      write(*,*) 'Phase-Shift Pass.'
*========= Solving the equation along the X-direction.=================*

      DO 987 IY = 1, NY

         DO 808 IL=1, LOOPP

         BETAX1 = Bcoe(IL)/(4.0*DX2)
         BETAX2 = Acoe(IL)*Dtau/(8.0*DX2)

         DO IX=1, NX

            VVtau2   = VVtau(IX, IY)*VVtau(IX, IY)
            cos_sita2= 1.0-RPh2*VVtau2*0.25
            cos_sita = SQRT(cos_sita2)

            RWK=REAL(RKW(IX, IY))
            AWK=AIMAG(RKW(IX, IY))
            RWK2=RWK*RWK
            AWK2=AWK*AWK

            RAWK2_minus   = RWK2-AWK2
            RAWK2_plus    = RWK2+AWK2
            RAWK2_plus2   = RAWK2_plus*RAWK2_plus
            RAWK_multiple = RWK*AWK

            Beta11 = BETAX1/cos_sita2
            Atmp1  = (RAWK2_minus/RAWK2_plus2)*Beta11 + ALFA
            Btmp1  = (RAWK_multiple/RAWK2_plus2)*Beta11

            Beta22 = BETAX2*VVtau(IX, IY)/cos_sita
            Atmp2  = (AWK/RAWK2_plus)*Beta22
            Btmp2  = (RWK/RAWK2_plus)*Beta22

         AA_X(IX) = cmplx(Atmp1+Atmp2, 2.0*Btmp1-Btmp2)
         BB_X(IX) = cmplx(1-2.0*(Atmp1+Atmp2), -(4.0*Btmp1-2.0*Btmp2))
         CC_X(IX) = cmplx(Atmp1-Atmp2, 2.0*Btmp1+Btmp2)
         DD_X(IX) = cmplx(1-2.0*(Atmp1-Atmp2), -(4.0*Btmp1+2.0*Btmp2))

         END DO

         YY_X(1) = AA_X(1)*DATW_W(2, IY) + BB_X(1)*DATW_W(1, IY)
         DO IX=2, NX-1
            YY_X(IX) = AA_X(IX)*(DATW_W(IX+1, IY)+DATW_W(IX-1, IY))
     +               + BB_X(IX)*DATW_W(IX, IY)
         END DO
         YY_X(NX) = AA_X(NX)*DATW_W(NX-1, IY) + BB_X(NX)*DATW_W(NX, IY)

         CALL COMPLEX_TRIDAG(CC_X, DD_X, YY_X, XX_X, NX)

         DO IX=1, NX
            DATW_W(IX, IY)=XX_X(IX)
         END DO

C        For the attenuation boundary. Left  side.
         DO IX = 1, Nmx0
            IXX0 = Nmx0 - IX + 1
            DATW_W(IX, IY) = XX_X(IX)*(1.0+cos(PAI*IXX0/Nmx0))*0.5
         END DO

C        For the attenuation boundary. Right side.
         DO IX = Nmx1+Nmx0+1, NX
            IXX0 = IX - Nmx1 - Nmx0
            DATW_W(IX, IY) = XX_X(IX)*(1.0+cos(PAI*IXX0/Nmx0))*0.5
         END DO

808      CONTINUE

987   CONTINUE

*========= Solving the equation along the Y-direction.=================*

      DO 988 IX = 1, NX

         DO 909 IL=1, LOOPP

         BETAY1 = Bcoe(IL)/(4.0*DY2)
         BETAY2 = Acoe(IL)*Dtau/(8.0*DY2)

         DO IY=1, NY

            VVtau2   = VVtau(IX, IY)*VVtau(IX, IY)
            cos_sita2= 1.0-RPh2*VVtau2*0.25
            cos_sita = SQRT(cos_sita2)

            RWK=REAL(RKW(IX, IY))
            AWK=AIMAG(RKW(IX, IY))
            RWK2=RWK*RWK
            AWK2=AWK*AWK

            RAWK2_minus   = RWK2-AWK2
            RAWK2_plus    = RWK2+AWK2
            RAWK2_plus2   = RAWK2_plus*RAWK2_plus
            RAWK_multiple = RWK*AWK

            Beta11 = BETAY1/cos_sita2
            Atmp1  = (RAWK2_minus/RAWK2_plus2)*Beta11 + ALFA
            Btmp1  = (RAWK_multiple/RAWK2_plus2)*Beta11

            Beta22 = BETAY2*VVtau(IX, IY)/cos_sita
            Atmp2  = (AWK/RAWK2_plus)*Beta22
            Btmp2  = (RWK/RAWK2_plus)*Beta22

         AA_Y(IY) = cmplx(Atmp1+Atmp2, 2.0*Btmp1-Btmp2)
         BB_Y(IY) = cmplx(1-2.0*(Atmp1+Atmp2), -(4.0*Btmp1-2.0*Btmp2))
         CC_Y(IY) = cmplx(Atmp1-Atmp2, 2.0*Btmp1+Btmp2)
         DD_Y(IY) = cmplx(1-2.0*(Atmp1-Atmp2), -(4.0*Btmp1+2.0*Btmp2))

         END DO

         YY_Y(1) = AA_Y(1)*DATW_W(IX, 2) + BB_Y(1)*DATW_W(IX, 1)
         DO IY=2, NY-1
            YY_Y(IY) = AA_Y(IY)*(DATW_W(IX, IY+1)+DATW_W(IX, IY-1))
     +               + BB_Y(IY)*DATW_W(IX, IY)
         END DO
         YY_Y(NY) = AA_Y(NY)*DATW_W(IX, NY-1) + BB_Y(NY)*DATW_W(IX, NY)

         CALL COMPLEX_TRIDAG(CC_Y, DD_Y, YY_Y, XX_Y, NY)

         DO Iy=1, NY
            DATW_W(IX, IY)=XX_Y(IY)
         END DO

C        For the attenuation boundary. Left  side.
         DO IY = 1, Nmy0
            IYY0 = Nmy0 - IY + 1
            DATW_W(IX, IY) = XX_Y(IY)*(1.0+cos(PAI*IYY0/Nmx0))*0.5
         END DO

C        For the attenuation boundary. Right side.
         DO IY = Nmy1+Nmy0+1, NY
            IYY0 = IY - Nmy1 - Nmy0
            DATW_W(IX, IY) = XX_Y(IY)*(1.0+cos(PAI*IYY0/Nmx0))*0.5
         END DO

909      CONTINUE

988   CONTINUE

!     Compensate X-Y splitting error.
      VZ = 0
      QZ = 0
      DO IY=1, NY
      DO IX=1, NX
         VZ = VZ + VVtau(IX,IY)
         QZ = QZ + Qvalue(IX, IY)
      END DO
      END DO
      VZ = VZ/(NX*NY)
      QZ = QZ/(NX*NY)

      Gama           = Pai21*Atan(0.5/QZ)
      Tmp0           = (W/Wr)**Gama
      V_dispersion   = VZ*Tmp0
      Coe_Attenuation= tan(0.5*pai*Gama)*W*TMP0/VZ
      CKW            = CMPLX(W/V_dispersion, Coe_Attenuation)

      cosalfa = sqrt(1.0 - 0.25*VZ*VZ*RPh2)
      CALL COMPENSATE_PSFD_VISCO(DATW_W, XX, XX_X, YY_Y, TMP2D,
     +     CKW, NX, Dkx, Dx, NY, Dky, Dy,
     +     Dtau, VZ, W, ALFA, RKz_max, cosalfa)

      RETURN
      END

*======================================================================*
      SUBROUTINE COMPENSATE_PSFD_VISCO(DATW_W, XX, XX_X, YY_Y, TMP2D,
     +     CKW, KxNUM, Dkx, Dx, KYNUM, Dky, Dy,
     +     Dtau, VZ, W, ALFA, RKz_max, cosalfa)
      PARAMETER(Wflag=-1.0, ACOE=0.4575, BCOE=0.4575)

      INTEGER  KxNUM, KYNUM
      REAL     Dkx, Dky, Dx, Dy, Dtau, VZ, W, ALFA, RKz_max, cosalfa
      COMPLEX  CKW

      COMPLEX  XX(KxNUM,KyNUM), DATW_W(KxNUM,KyNUM), TMP2D(KyNUM,KxNUM)
      COMPLEX  XX_X(KxNUM), YY_Y(KyNUM)

!     Definition of local-parameters
      COMPLEX  a0, a00, CKW2cos2, CYY, BXX
      COMPLEX  tmpsqrt, tmpsqrt2, tmpsqrt_vkw
      COMPLEX  PHASE1, PHASEX, PHASEY, PHASE1XY
      COMPLEX  ImagPart

      REAL     Dx2, Dy2, W2, VZ2, RKz_max2
      REAL     tmpa, tmpa2, tmpmin, AR, BI
      REAL     CKW2cos2R2, CKW2cos2I2, CKW2cos2Amp2

      Dx2  = Dx*Dx
      Dy2  = Dy*Dy
      W2   = W*W
      VZ2  = VZ*VZ
      cos2 = cosalfa*cosalfa

      a0   = 0.25/(CKW*CKW*cosalfa*cosalfa)
      a00  = VZ*CKW*cosalfa*Dtau

      CKW2cos2     = CKW*CKW*cos2
      CKW2cos2R2   = REAL(CKW2cos2)*REAL(CKW2cos2)
      CKW2cos2I2   = AIMAG(CKW2cos2)*AIMAG(CKW2cos2)
      CKW2cos2Amp2 = CKW2cos2R2+CKW2cos2I2
      CKW2cos2Amp2 = W2*cos2/VZ2

      RKz_max2 = RKz_max*RKz_max
      tmpmin   = Rkz_max2*W2/VZ2
      ImagPart = CMPLX(0.0,1.0)

      DO Iky=1, KyNUM
      DO Ikx=1, KxNUM
         XX(Ikx, Iky) = DATW_W(Ikx, Iky)
      END DO
      END DO

C     CALL FFT2D(XX, XX_X, YY_Y, KxNUM, KyNUM, +1)
      CALL FFT2D_OMP(XX, TMP2D,KxNUM, KyNUM, +1)

      DO 6660 Iky  = 1, KyNUM/2 + 1

         Jky    = KyNUM + 2 - Iky
         Rky    = (Iky-1)*Dky
         Rky2   = Rky*Rky                           !Ky^2
         Rkyy2  = 2.0*(1.0-cos(Rky*Dy))/Dy2         !Ky^^2

         CYY    = ACOE*a0*Rkyy2/(1.0-(ALFA*Dy2+BCOE*a0)*Rkyy2)
         PHASEY = CEXP(Wflag*a00*CYY*ImagPart)

         DO 6661 Ikx  = 1, KxNUM/2 + 1

            Jkx    = KxNUM + 2 - Ikx
            Rkx    = (Ikx-1)*Dkx
            Rkx2   = Rkx*Rkx                        !Kx^2
            Rkxx2  = 2.0*(1.0-cos(Rkx*Dx))/Dx2      !Kx^^2

            tmpa2  = CKW2cos2Amp2 - 0.25*(Rky2+Rkx2)

            IF (tmpa2 .le. 0.0) THEN
               XX(Ikx, Iky) = 0.0
               IF (Ikx.GT.1)              XX(Jkx,Iky)=0.0
               IF (Iky.GT.1)              XX(Ikx,Jky)=0.0
               IF (Iky.GT.1.AND.Ikx.GT.1) XX(Jkx,Jky)=0.0
               GOTO 6661
            END IF

            IF (tmpa2 .le. tmpmin) THEN
               XX(Ikx, Iky) = 0.0
               IF (Ikx.GT.1)              XX(Jkx,Iky)=0.0
               IF (Iky.GT.1)              XX(Ikx,Jky)=0.0
               IF (Iky.GT.1.AND.Ikx.GT.1) XX(Jkx,Jky)=0.0
               GOTO 6661
            END IF

C           GOTO 6661

            tmpsqrt2    = CKW2cos2 - 0.25*(Rky2+Rkx2)
            tmpsqrt     = CSQRT(tmpsqrt2)
            tmpsqrt_vkw = tmpsqrt - CKW*cosalfa
            AR = REAL(tmpsqrt_vkw)
            BI = (+1.0)*ABS(AIMAG(tmpsqrt_vkw))
            PHASE1 = CEXP(CMPLX(0.0,Wflag*VZ*Dtau*AR))*EXP(VZ*Dtau*BI)

            BXX    = ACOE*a0*Rkxx2/(1.0-(ALFA*Dx2+BCOE*a0)*Rkxx2)
            PHASEX = CEXP(Wflag*a00*BXX*ImagPart)

            PHASE1XY     = PHASE1*PHASEX*PHASEY
            XX(Ikx, Iky) = XX(Ikx, Iky)*PHASE1XY

            IF (Ikx.GT.1)              XX(Jkx,Iky)=XX(Jkx,Iky)*PHASE1XY
            IF (Iky.GT.1)              XX(Ikx,Jky)=XX(Ikx,Jky)*PHASE1XY
            IF (Iky.GT.1.AND.Ikx.GT.1) XX(Jkx,Jky)=XX(Jkx,Jky)*PHASE1XY

6661     CONTINUE

6660  CONTINUE

C     CALL FFT2D(XX, XX_X, YY_Y, KxNUM, KyNUM, -1)
      CALL FFT2D_OMP(XX, TMP2D,KxNUM, KyNUM, -1)

      DO Iky=1, KyNUM
      DO Ikx=1, KxNUM
         DATW_W(Ikx, Iky) = XX(Ikx, Iky)
      END DO
      END DO

      RETURN
      END

*======================================================================*
      SUBROUTINE COMPLEX_TRIDAG(A, B, R, U, N)

      INTEGER   N, NMAX
      COMPLEX   A(N), B(N), R(N), U(N)
      PARAMETER (NMAX=5000)
      INTEGER   j

      complex   bet, gam(NMAX)

      if(NMAX.LT.N)   then
         write(*,*)' NMAX.LT.N,  ERROR in COMPLEX_TRIDAG()'
         stop
      endif

      if(REAL(b(1)).eq.0.0.AND.AIMAG(b(1)).eq.0.0)
     +              pause 'tridag: rewrite equations'
      BET=B(1)
      U(1)=R(1)/BET
      DO 11 J=2, N
         GAM(j)=A(J-1)/BET
         BET=B(J)-A(J)*GAM(J)
         IF (REAL(BET).EQ.0.AND.AIMAG(BET).EQ.0.0)
     +    PAUSE 'TRIDAG FAILED'
         U(J)=(R(J)-A(J)*U(J-1))/BET
11    CONTINUE

      DO 12 J=N-1,1,-1
         U(J)=U(J)-GAM(J+1)*U(J+1)
12    CONTINUE

      RETURN
      END

*======================================================================*
      SUBROUTINE  WRITE_2D_BUFFER_REAL(BUF, LT, NX, FN)
      IMPLICIT NONE
      CHARACTER*256 FN
      INTEGER  LT, NX, IT, IX
      REAL     BUF(LT,NX)

      OPEN(13, FILE=FN, access='direct',recl=LT*NX)
      WRITE(13,rec=1) ((BUF(IT,IX),IT=1,LT),IX=1,NX)
      CLOSE(13)

      RETURN
      END

*======================================================================*
      SUBROUTINE READPAR(FN1, FN2, FN3, FN4,
     +         Ntau, Iphr1, Iphr2, F1, F2, F3, F4, IBREAK, ViscoFlag,
     +         NtauExtraP, ITauBreak, Ktaur, Ktauw, FNTMPW, FNTMPR,
     +         IphrR1, IphrR2 )
      
      INTEGER Ntau, Iphr1, Iphr2, IBREAK, ViscoFlag
      INTEGER NtauExtraP, ITauBreak, Ktaur, Ktauw
      INTEGER IphrR1, IphrR2
      REAL    F1, F2, F3, F4, DIP_max
      CHARACTER*256 FN1, FN2, FN3, FN4, CTMP, FNTMPW, FNTMPR
 
      OPEN(13, FILE='3d_offset_planewave_pstm.par')

C     READ (13,'(A)') FN1              ! Offset-planewave decompostion data.
C     READ (13,'(A)') FN2              ! Interval velocity field(time domain).
C     READ (13,'(A)') FN3              ! Common-ph imaging volume, Ntau*Ncdp*Nline*Nphr
C     READ (13,'(A)') FN4              ! Temp dir and file preffix.
C     READ (13,*)     Ntau             ! Extrapolation number.
C     READ (13,*)     Iphr1            ! The first Ray-paramter number.
C     READ (13,*)     Iphr2            ! The last  Ray-paramter number.
C     READ (13,*)     F1, F2, F3, F4   ! The frequency range for imaging.
C     READ (13,*)     IBREAK           ! Break-Point flag. =0, new; =1, continue.

      READ (13,'(A)') FN1              ! Offset-planewave decompostion data.
      READ (13,'(A)') FN2              ! Interval velocity field(time domain).
      READ (13,'(A)') FN3              ! Common-ph imaging volume, Ntau*Ncdp*Nline*Nphr
      READ (13,'(A)') FN4              ! Temp dir and file preffix.
      READ (13,*) CTMP, Ntau           ! Trace Length of Imaging-volume.
      READ (13,*) CTMP, Iphr1          ! The first Ray-paramter number.
      READ (13,*) CTMP, Iphr2          ! The last  Ray-paramter number.
      READ (13,*) CTMP, F1, F2, F3, F4 ! The frequency range for imaging.
      READ (13,*) CTMP, IBREAK         ! Break-Point flag. =0, new; =1, continue.
      READ (13,*) CTMP, ViscoFlag      ! Visco-acoustic flag. =0, acoustic; =1, visco-.
      READ (13,*) CTMP, NtauExtraP     ! Extrapolation number.
      READ (13,*) CTMP, ITauBreak      ! Tau-Axis Break-Point flag. =0, from surface; =1, from Ktau.
      READ (13,*) CTMP, Ktaur          ! Tau-Axis Extrapolation Start-Reading Point.
      READ (13,*) CTMP, Ktauw          ! Tau-Axis Extrapolation Start-Writing Point.
      READ (13,'(A)') FNTMPW           ! Temp dir and file preffix.
      READ (13,'(A)') FNTMPR           ! Temp dir and file preffix.
      READ (13,*) CTMP, IphrR1         ! The first Ray-paramter number.
      READ (13,*) CTMP, IphrR2         ! The last  Ray-paramter number.

      CLOSE(13)
C     WRITE(*,*) trim(FN1)
C     WRITE(*,*) trim(FN2)
C     WRITE(*,*) trim(FN3)
C     WRITE(*,*) trim(FN4)
C     WRITE(*,*) Ntau
C     WRITE(*,*) Iphr1
C     WRITE(*,*) Iphr2
C     WRITE(*,*) F1, F2, F3, F4
C     WRITE(*,*) IBREAK
C     WRITE(*,*) ViscoFlag

      RETURN
      END
 
*=====================================================================*
       SUBROUTINE  FFT1D(A,N,INV)

       implicit none

       integer  n,i,inv
       COMPLEX A(N)

       IF(INV.EQ.(-1)) THEN
         DO I=1,N
            a(i)=conjg(a(i))/float(n)
         ENDDO
         CALL FFTCC(a,N)
         do i=1,n
            a(i)=conjg(a(i))
         enddo
        else
           CALL FFTCC(a,N)
        ENDIF

        RETURN
        END

*=====================================================================*
      SUBROUTINE FFT2D_OMP(A2D, AT2D, N1, N2, INV)
C     Perform 2D FFT with flag: INV. =0, forward; =1, backward.
C     A2D(N1,N2)  is the input 2d buffer, overwritten by output.
C     AT2D(N1,N2) is a temporal buffer.
C     Author: FengBo, 2010-07-17.
      IMPLICIT NONE
      INTEGER N1, N2, INV
      COMPLEX A2D(N1,N2), AT2D(N2,N1)

      INTEGER I1, I2

C     1D FFT to the fast-dimension.
C$OMP PARALLEL
C$OMP DO
      DO I2 = 1, N2
         CALL FFT1D(A2D(1,I2), N1, INV)
      END DO
C$OMP END DO
C$OMP END PARALLEL

C     Transpose A2D to  AT2D
      DO I1 = 1, N1
      DO I2 = 1, N2
         AT2D(I2,I1) = A2D(I1,I2)
      END DO
      END DO

C     1D FFT to the slow-dimension, using the transposed buffer AT2D
C$OMP PARALLEL
C$OMP DO
      DO I1 = 1, N1
         CALL FFT1D(AT2D(1,I1), N2, INV)
      END DO
C$OMP END DO
C$OMP END PARALLEL

C     Transpose AT2D to  A2D
      DO I1 = 1, N1
      DO I2 = 1, N2
         A2D(I1,I2) = AT2D(I2,I1)
      END DO
      END DO

      RETURN
      END
*======================================================================*
        SUBROUTINE FFT2D(AA, Ahx, Amx, KhxNum, KmxNum, INV)

        INTEGER KmxNum, KhxNum, INV
        COMPLEX AA(KhxNum, KmxNum), Amx(KmxNum), Ahx(KhxNum)

        DO 567 Imx=1, KmxNum

         DO 667 Ihx=1, KhxNum
          Ahx(Ihx)=AA(Ihx, Imx)
667      CONTINUE

         CALL FFT1D(Ahx, KhxNum, INV)

         DO 767 Ihx=1, KhxNum
          AA(Ihx, Imx)=Ahx(Ihx)
767      CONTINUE

567     CONTINUE

        DO 568 Ihx=1, KhxNum

         DO 668 Imx=1, KmxNum
          Amx(Imx)=AA(Ihx, Imx)
668      CONTINUE

         CALL FFT1D(Amx, KmxNum, INV)

         DO 868 Imx=1, KmxNum
          AA(Ihx, Imx)=Amx(Imx)
868      CONTINUE

568     CONTINUE

        RETURN
        END

*===============================================================C 
C                                                               C
C   ROUTINE NAME :  FFTCC                                       C
C                                                               C
C --------------------------------------------------------      C
C                                                               C
C   PURPOSE :  COMPUTE THE FAST FOURIER TRANSFORM OF A          C
C              COMPLEX VALUED SEQUENCE                          C
C                                                               C
C   USAGE :  CALL FFTCC(A,N)                                    C
C                                                               C
C   ARGUMENTS : A    COMPLEX VECTOR OF LENGTH N,                C
C                    ON INPUT A CONTAINS THE                    C
C                    COMPLEX VALUED SEQUENCE TO BE              C
C                    ON OUTPUT A IS REPLACED BY THE             C
C                    FOURIER TRANSFORM.                         C
C               N    INPUT NUMBER OF DATA TO BY                 C
C                    TRANSFORMED. N MAY BE ANY POSITIVE         C
C                    INTEGER.                                   C
C                                                               C
C  NOTATION :   INFORMATION ON SPECIAL NOTATION AND             C
C               CONVENTIONS IS AVAILABLE IN THE MANUAL          C
C               INTRODUCTION OR THROUGH ROUTIME HELP            C
C                                                               C
C  REMARKS : 1. FFTCC COMPUTES THE FOURIER TRANSFORM,X,         C
C              ACCORDING TO THE FOLLOWING FORMULA;              C
C                                                               C
C              X(K+1)= SUM FORM J=0 TO N-1 OF                   C
C                      A(J+1)*CEXP((0.0,(2.0*PI*J*K)/N))        C
C             FOR K=0,1,...,N-1 AND PI=3.1415926...             C
C                                                               C
C             NOTE THAT X OVERWRITES A ON OUTPUT.               C
C                                                               C
C           2. FFTCC CAN BE USED TO COMPUTE                     C
C                                                               C
C             X(K+1)= (1/N)*SUM FORM J = 0 TO N-1 OF            C
C                     A(J+1)*CEXP((0.0,(-2.0*PI*J*K)/N))        C
C             FOR K=0,1,...,N-1 AND PI=3.1415926...             C 
C                                                               C
C          BY PERFORMING THE FOLLOWING STEPS;                   C
C                                                               C
C             DO 10 I=1,N                                       C
C        10   A(I)= CONJG(A(I))                                 C
C             CALL FFTCC(A,N)                                   C
C             DO 20 I=1,N                                       C
C        20   A(I)= CONJG(A(I))/N                               C
C                                                               C
C --------------------------------------------------------------C
      SUBROUTINE FFTCC(A,N)
C     SPECIFICATIONS FOR ARGUMENTS
      INTEGER IWK(4096*2)
      REAL WK(4096*2)
      COMPLEX A(N)
C     SPECIFICATIONS FOR LOCAL VARIABLES
      REAL A0,A1,A2,A3,A4,B0,B1,B2,B3,B4,ZERO,ONE,TWO,
     X     Z0(2),Z1(2),Z2(2),Z3(2),Z4(2),HALF
      COMPLEX ZA0,ZA1,ZA2,ZA3,ZA4,AK2
      EQUIVALENCE (ZA0,Z0(1)),(ZA1,Z1(1)),(ZA2,Z2(1)),
     X            (ZA3,Z3(1)),(A0,Z0(1)),(B0,Z0(2)),
     X            (A1,Z1(1)),(B1,Z1(2)),(A2,Z2(1)),
     X            (B2,Z2(2)),(A3,Z3(1)),(B3,Z3(2)),
     X            (ZA4,Z4(1)),(Z4(1),A4),(Z4(2),B4)
      DATA ZERO/0.0/,HALF/0.5/,ONE/1.0/,TWO/2.0/
C    FIRST EXECUTABLE STATEMENT
C    RAD = MATH. CONS. 2*PI
      RAD=2*3.1415926
      C30=0.5*SQRT(3.0)
      IF(N.EQ.1) GO TO 9005
      K=N
      M=0
      J=2
      JJ=4
      JF=0
C    DETERMINE THE SQUARE FACTORS OF N
      IWK(1)=1
    5 I=K/JJ
      IF(I*JJ.NE.K) GO TO 10
      M=M+1
      IWK(M+1)=J
      K=I
      GO TO 5
   10 J=J+2
      IF(J.EQ.4) J=3
      JJ=J*J
      IF(JJ.LE.K) GO TO 5
      KT=M
C    DETERMINE THE REMAINING FACTORS OF N
      J=2
   15 I=K/J
      IF(I*J.NE.K) GO TO 20
      M=M+1
      IWK(M+1)=J
      K=I
      GO TO 15
   20 J=J+1
      IF(J.EQ.3) GO TO 15
      J=J+1
      IF(J.LE.K) GO TO 15
      K=IWK(M+1)
      IF(IWK(KT+1).GT.IWK(M+1)) K=IWK(KT+1)
      IF(KT.LE.0) GO TO 30
      KTP=KT+2
      DO 25 I=1,KT
      J=KTP-I
      M=M+1
      IWK(M+1)=IWK(J)
   25 CONTINUE
   30 MP=M+1
      IC=MP+1
      ID=IC+MP
      ILL=ID+MP
      IRD=ILL+MP+1
      ICC=IRD+MP
      ISS=ICC+MP
      ICK=ISS+MP
      ISK=ICK+K
      ICF=ISK+K
      ISF=ICF+K
      IAP=ISF+K
      KD2=(K-1)/2+1
      IBP=IAP+KD2
      IAM=IBP+KD2
      IBM=IAM+KD2
      MM1=M-1
      I=1
   35 L=MP-I
      J=IC-I
      IWK(ILL+L)=0
      IF((IWK(J-1)+IWK(J)).EQ.4) IWK(ILL+L)=1
      IF(IWK(ILL+L).EQ.0) GO TO 40
      I=I+1
      L=L-1
      IWK(ILL+L)=0
   40 I=I+1
      IF(I.LE.MM1) GO TO 35
      IWK(ILL+1)=0
      IWK(ILL+MP)=0
      IWK(IC)=1
      IWK(ID)=N
      DO 45 J = 1 ,M
      K=IWK(J+1)
      IWK(IC+J)=IWK(IC+J-1)*K
      IWK(ID+J)=IWK(ID+J-1)/K
      WK(IRD+J)=RAD/IWK(IC+J)
      C1=RAD/K
      IF(K.LE.2) GO TO 45
      WK(ICC+J)=COS(C1)
      WK(ISS+J)=SIN(C1)
   45 CONTINUE
      MM=M
      IF(IWK(ILL+M).EQ.1) MM=M-1
      IF(MM.LE.1) GO TO 50
      SM=IWK(IC+MM-2)*WK(IRD+M)
      CM=COS(SM)
      SM=SIN(SM)
   50 KB=0
      KN=N
      JJ=0
      I=1
      C1=1.0
      S1=0.0
      L1=1
   55 IF(IWK(ILL+I+1).EQ.1) GO TO 60
      KF=IWK(I+1)
      GO TO 65
   60 KF=4
      I=I+1
   65 ISP=IWK(ID+I)
      IF(L1.EQ.1) GO TO 70
      S1=JJ*WK(IRD+I)
      C1=COS(S1)
      S1=SIN(S1)
C    FACTORS OF 2, 3, AND 4 ARE
C    HANDLESEPARATELY.
   70 IF(KF.GT.4) GO TO 140
      GO TO (75,75,90,115), KF
   75 K0=KB+ISP
      K2=K0+ISP
      IF(L1.EQ.1) GO TO 85
   80 K0=K0-1
      IF(K0.LT.KB) GO TO 190
      K2=K2-1
      ZA4=A(K2+1)
      A0=A4*C1-B4*S1
      B0=A4*S1+B4*C1
      A(K2+1)=A(K0+1)-ZA0
      A(K0+1)=A(K0+1)+ZA0
      GO TO 80
   85 K0=K0-1
      IF(K0.LT.KB) GO TO 190
      K2=K2-1
      AK2=A(K2+1)
      A(K2+1)=A(K0+1)-AK2
      A(K0+1)=A(K0+1)+AK2
      GO TO 85
   90 IF(L1.EQ.1) GO TO 95
      C2=C1*C1-S1*S1
      S2=2.0*C1*S1
   95 JA=KB+ISP-1
      KA=JA+KB
      IKB=KB+1
      IJA=JA+1
      DO 110 II=IKB,IJA
      K0=KA-II+1
      K1=K0+ISP
      K2=K1+ISP
      ZA0=A(K0+1)
      IF(L1.EQ.1) GO TO 100
      ZA4=A(K1+1)
      A1=A4*C1-B4*S1
      B1=A4*S1+B4*C1
      ZA4=A(K2+1)
      A2=A4*C2-B4*S2
      B2=A4*S2+B4*C2
      GO TO 105
  100 ZA1=A(K1+1)
      ZA2=A(K2+1)
  105 A(K0+1)=CMPLX(A0+A1+A2,B0+B1+B2)
      A0=-HALF*(A1+A2)+A0
      A1=(A1-A2)*C30
      B0=-HALF*(B1+B2)+B0
      B1=(B1-B2)*C30
      A(K1+1)=CMPLX(A0-B1,B0+A1)
      A(K2+1)=CMPLX(A0+B1,B0-A1)
  110 CONTINUE
      GO TO 190
  115 IF(L1.EQ.1) GO TO 120
      C2=C1*C1-S1*S1
      S2=2.0*C1*S1
      C3=C1*C2-S1*S2
      S3=S1*C2+C1*S2
  120 JA=KB+ISP-1
      KA=JA+KB
      IKB=KB+1
      IJA=JA+1
      DO 135 II=IKB,IJA
      K0=KA-II+1
      K1=K0+ISP
      K2=K1+ISP
      K3=K2+ISP
      ZA0=A(K0+1)
      IF(L1.EQ.1) GO TO 125
      ZA4=A(K1+1)
      A1=A4*C1-B4*S1
      B1=A4*S1+B4*C1
      ZA4=A(K2+1)
      A2=A4*C2-B4*S2
      B2=A4*S2+B4*C2
      ZA4=A(K3+1)
      A3=A4*C3-B4*S3
      B3=A4*S3+B4*C3
      GO TO 130
  125 ZA1=A(K1+1)
      ZA2=A(K2+1)
      ZA3=A(K3+1)
  130 A(K0+1)=CMPLX(A0+A1+A2+A3,B0+B1+B2+B3)
      A(K1+1)=CMPLX(A0+A2-A1-A3,B0+B2-B1-B3)
      A(K2+1)=CMPLX(A0-A2-B1+B3,B0-B2+A1-A3)
      A(K3+1)=CMPLX(A0-A2+B1-B3,B0-B2-A1+A3)
  135 CONTINUE
      GO TO 190
  140 JK=KF-1
      KH=JK/2
      K3=IWK(ID+I-1)
      K0=KB+ISP
      IF(L1.EQ.1) GO TO 150
      K=JK-1
      WK(ICF+1)=C1
      WK(ISF+1)=S1
      DO 145 J=1,K
      WK(ICF+J+1)=WK(ICF+J)*C1-WK(ISF+J)*S1
      WK(ISF+J+1)=WK(ICF+J)*S1+WK(ISF+J)*C1
  145 CONTINUE
  150 IF(KF.EQ.JF) GO TO 160
      C2=WK(ICC+I)
      WK(ICK+1)=C2
      WK(ICK+JK)=C2
      S2=WK(ISS+I)
      WK(ISK+1)=S2
      WK(ISK+JK)=-S2
      DO 155 J=1,KH
      K=JK-J
      WK(ICK+K)=WK(ICK+J)*C2-WK(ISK+J)*S2
      WK(ICK+J+1)=WK(ICK+K)
      WK(ISK+J+1)=WK(ICK+J)*S2+WK(ISK+J)*C2
      WK(ISK+K)=-WK(ISK+J+1)
  155 CONTINUE
  160 K0=K0-1
      K1=K0
      K2=K0+K3
      ZA0=A(K0+1)
      A3=A0
      B3=B0
      DO 175 J=1,KH
      K1=K1+ISP
      K2=K2-ISP
      IF(L1.EQ.1) GO TO 165
      K=KF-J
      ZA4=A(K1+1)
      A1=A4*WK(ICF+J)-B4*WK(ISF+J)
      B1=A4*WK(ISF+J)+B4*WK(ICF+J)
      ZA4=A(K2+1)
      A2=A4*WK(ICF+K)-B4*WK(ISF+K)
      B2=A4*WK(ISF+K)+B4*WK(ICF+K)
      GO TO 170
  165 ZA1=A(K1+1)
      ZA2=A(K2+1)
  170 WK(IAP+J)=A1+A2
      WK(IAM+J)=A1-A2
      WK(IBP+J)=B1+B2
      WK(IBM+J)=B1-B2
      A3=A1+A2+A3
      B3=B1+B2+B3
  175 CONTINUE
      A(K0+1)=CMPLX(A3,B3)
      K1=K0
      K2=K0+K3
      DO 185 J=1,KH
      K1=K1+ISP
      K2=K2-ISP
      JK=J
      A1=A0
      B1=B0
      A2=0.0
      B2=0.0
      DO 180 K=1,KH
      A1=A1+WK(IAP+K)*WK(ICK+JK)
      A2=A2+WK(IAM+K)*WK(ISK+JK)
      B1=B1+WK(IBP+K)*WK(ICK+JK)
      B2=B2+WK(IBM+K)*WK(ISK+JK)
      JK=JK+J
      IF(JK.GE.KF) JK=JK-KF
  180 CONTINUE
      A(K1+1)=CMPLX(A1-B2,B1+A2)
      A(K2+1)=CMPLX(A1+B2,B1-A2)
  185 CONTINUE
      IF(K0.GT.KB) GO TO 160
      JF=KF
  190 IF(I.GE.MM) GO TO 195
      I=I+1
      GO TO 55
  195 I=MM
      L1=0
      KB=IWK(ID+I-1)+KB
      IF(KB.GE.KN) GO TO 215
  200 JJ=IWK(IC+I-2)+JJ
      IF(JJ.LT.IWK(IC+I-1)) GO TO 205
      I=I-1
      JJ=JJ-IWK(IC+I)
      GO TO 200
  205 IF(I.NE.MM) GO TO 210
      C2=C1
      C1=CM*C1-SM*S1
      S1=SM*C2+CM*S1
      GO TO 70
  210 IF(IWK(ILL+I).EQ.1)I=I+1
      GO TO 55
  215 I=1
      JA=KT-1
      KA=JA+1
      IF(JA.LT.1) GO TO 225
      DO 220 II=1,JA
      J=KA-II
      IWK(J+1)=IWK(J+1)-1
      I=IWK(J+1)+I
  220 CONTINUE
C   THE RESULT IS NOW PERMUTED TO NORMAL ORDER.
  225 IF(KT.LE.0) GO TO 270
      J=1
      I=0
      KB=0
  230 K2=IWK(ID+J)+KB
      K3=K2
      JJ=IWK(IC+J-1)
      JK=JJ
      K0=KB+JJ
      ISP=IWK(IC+J)-JJ
  235 K=K0+JJ
  240 ZA4=A(K0+1)
      A(K0+1)=A(K2+1)
      A(K2+1)=ZA4
      K0=K0+1
      K2=K2+1
      IF(K0.LT.K) GO TO 240
      K0=K0+ISP
      K2=K2+ISP
      IF(K0.LT.K3) GO TO 235
      IF(K0.GE.K3+ISP) GO TO 245
      K0=K0-IWK(ID+J)+JJ
      GO TO 235
  245 K3=IWK(ID+J)+K3
      IF(K3-KB.GE.IWK(ID+J-1)) GO TO 250
      K2=K3+JK
      JK=JK+JJ
      K0=K3-IWK(ID+J)+JK
      GO TO 235
  250 IF(J.GE.KT) GO TO 260
      K=IWK(J+1)+I
      J=J+1
  255 I=I+1
      IWK(ILL+I)=J
      IF(I.LT.K) GO TO 255
      GO TO 230
  260 KB=K3
      IF(I.LE.0) GO TO 265
      J=IWK(ILL+I)
      I=I-1
      GO TO 230
  265 IF(KB.GE.N) GO TO 270
      J=1
      GO TO 230
  270 JK=IWK(IC+KT)
      ISP=IWK(ID+KT)
      M=M-KT
      KB=ISP/JK-2
      IF(KT.GE.M-1) GO TO 9005
      ITA=ILL+KB+1
      ITB=ITA+JK
      IDM1=ID-1
      IKT=KT+1
      IM=M+1
      DO 275 J=IKT,IM
      IWK(IDM1+J)=IWK(IDM1+J)/JK
  275 CONTINUE
      JJ=0
      DO 290 J=1,KB
      K=KT
  280 JJ=IWK(ID+K+1)+JJ
      IF(JJ.LT.IWK(ID+K)) GO TO 285
      JJ=JJ-IWK(ID+K)
      K=K+1
      GO TO 280
  285 IWK(ILL+J)=JJ
      IF(JJ.EQ.J) IWK(ILL+J)=-J
  290 CONTINUE
C    DETERMINE THE PERMUTATION CYCLES
C    OF LENGTH GREATER THAN OR EQUAL TO TWO.
      DO 300 J=1,KB
      IF(IWK(ILL+J).LE.0) GO TO 300
      K2=J
  295 K2=IABS(IWK(ILL+K2))
      IF(K2.EQ.J) GO TO 300
      IWK(ILL+K2)=-IWK(ILL+K2)
      GO TO 295
  300 CONTINUE
C    REORDER A FOLLOWING THE PERMUATION CYCLES
      I=0
      J=0
      KB=0
      KN=N
      ITI=1
      ITT=1
  305 J=J+1
      IF(IWK(ILL+J).LT.0) GO TO 305
      K=IWK(ILL+J)
      K0=JK*K+KB
  310 ZA4=A(K0+I+1)
      WK(ITA+I)=A4
      WK(ITB+I)=B4
      I=I+1
      IF(I.LT.JK) GO TO 310
      I=0
  315 K=-IWK(ILL+K)
      JJ=K0
      K0=JK*K+KB
  320 A(JJ+I+1)=A(K0+I+1)
      I=I+1
      IF(I.LT.JK) GO TO 320
      I=0
      IF(K.NE.J) GO TO 315
  325 A(K0+I+1)=CMPLX(WK(ITA+I),WK(ITB+I))
      I=I+1
      IF(I.LT.JK) GO TO 325
      I=0
      IF(J.LT.K2) GO TO 305
      J=0
      KB=KB+ISP
      IF(KB.LT.KN) GO TO 305
 9005 RETURN
      END

*=====================================================================*

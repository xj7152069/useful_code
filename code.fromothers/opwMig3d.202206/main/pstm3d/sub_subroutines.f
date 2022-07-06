      SUBROUTINE INPUT_CURRENT_LAYER_VELOCITY(VVtau,
     +     Nvline_first, Nvline_final,
     +     Nvcdp_first, Nvcdp_final,
     +     Nvmx, Nvmy,
     +     Nline_first, Nline_final,
     +     Ncdp_first, Ncdp_final,
     +     Nmx, Nmy, Itau, IFN02,
     +     Nmx0, Nmx1, Nmy0, Nmy1)

      INTEGER Nvmx, Nvmy, Nmx, Nmy
      INTEGER Ncdp_first, Nvcdp_first, Nline_first, Nvline_first
      INTEGER Ncdp_final, Nvcdp_final, Nline_final, Nvline_final

      REAL    VVtau(Nmx, Nmy)
      REAL, ALLOCATABLE :: Vslice2d(:,:)

      ALLOCATE(Vslice2d(Nvmx, Nvmy), stat=ierror)

      DO imy = 1, Nmy
      DO imx = 1, Nmx
         VVtau(imx,imy) = 0
      END DO
      END DO

C     read a whole 2d depth-slice at once.
      read(IFN02, rec=Itau) ((Vslice2d(ix,iy),ix=1,Nvmx),
     +    iy=1,Nvmy)

      DO Iline = Nline_first, Nline_final
         DO Icdp = Ncdp_first, Ncdp_final
            imy = Iline - Nline_first + 1 + Nmy0
            imx = Icdp - Ncdp_first + Nmx0 + 1
            ivy = Iline - Nvline_first + 1
            ivx = Icdp - Nvcdp_first + 1
            VVtau(imx,imy) = Vslice2d(ivx, ivy)
         END DO
      END DO

C     N4=1
C     DO Iline = Nline_first, Nline_final
C        imy = Iline - Nline_first + 1 + Nmy0
C        imx = Ncdp_first - Ncdp_first + Nmx0 + 1
C        IFLAG = 1
C        CALL GDBFCREADCUTSLICEDATA(IFN02, N4, Itau, Iline,
C    +             Ncdp_first, Ncdp_final, IFLAG, VVtau(imx,imy))
C        IF (IFLAG .NE. 0) THEN
C           WRITE(*,*) 'GDBFCREADCUTSLICEDATA ERROR.'
C           RETURN
C        END IF
C     END DO

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

      DEALLOCATE(Vslice2d)

      RETURN
      END SUBROUTINE INPUT_CURRENT_LAYER_VELOCITY

*======================================================================*
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
      END SUBROUTINE make_order

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
      CHARACTER*1024 FN
      INTEGER  LT, NX, IT, IX
      REAL     BUF(LT,NX)

      OPEN(13, FILE=FN, access='direct',recl=LT*NX)
      WRITE(13,rec=1) ((BUF(IT,IX),IT=1,LT),IX=1,NX)
      CLOSE(13)

      RETURN
      END SUBROUTINE  WRITE_2D_BUFFER_REAL

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
      SUBROUTINE READPAR(FN_PAR, FN1, FN2, FN3, FN4,
     +         Nvtau, Iphr1, Iphr2,
     +         F1, F2, F3, F4,
     +         ViscoFlag, NtauExtraP)

      INTEGER Nvtau, Iphr1, Iphr2, ViscoFlag, NtauExtraP
      REAL    F1, F2, F3, F4
      CHARACTER*1024 FN_PAR, FN1, FN2, FN3, FN4, CTMP

      OPEN(13, FILE=FN_PAR)
      READ (13,'(A)') FN1              ! Offset-planewave decompostion data.
      READ (13,'(A)') FN2              ! Interval velocity field(time domain).
      READ (13,'(A)') FN3              ! the stacked image.
      READ (13,'(A)') FN4              ! Common-ph imaging volume, Ntau*Ncdp*Nline*Nphr
      READ (13,*) CTMP, Nvtau          ! for interval velocity model.
      READ (13,*) CTMP, Iphr1          ! The first Ray-paramter number.
      READ (13,*) CTMP, Iphr2          ! The last  Ray-paramter number.
      READ (13,*) CTMP, F1, F2, F3, F4 ! The frequency range for imaging.
      READ (13,*) CTMP, ViscoFlag      ! Visco-acoustic flag. =0, acoustic; =1, visco-.
      READ (13,*) CTMP, NtauExtraP     ! Extrapolation number.
      CLOSE(13)

      RETURN
      END SUBROUTINE READPAR

*=====================================================================*

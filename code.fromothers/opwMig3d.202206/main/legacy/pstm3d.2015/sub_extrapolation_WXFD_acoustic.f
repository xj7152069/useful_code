
      SUBROUTINE WXFD_EXTRAPOLATION_WXYCOMP_MPI_ACOUSTIC(
     +    vel3d,
     +    DATW_W, DATW_R, XX, VVtau, TMP2D,
     +    AA_X, BB_X, CC_X, DD_X, XX_X, YY_X,
     +    AA_Y, BB_Y, CC_Y, DD_Y, XX_Y, YY_Y,
     +    ACOE, BCOE, LOOPP, ALFA,
     +    NX, Dx, NY, Dy, Ntau, Dtau, W, RPh,
     +    Dkx, Dky, Rkz_max, DIP_max, DIP_min,
     +    Nvline_first, Nvline_final, Nvmy,
     +    Nvcdp_first,  Nvcdp_final,  Nvmx,
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

      INTEGER   Itau, weflag
      INTEGER   Itaus, Itaue
      REAL      vel3d(NX, NY, Ntau)

C     FengBo add temp file saving&reading.
      CHARACTER*1024 FNW, FNR

      RPh2     = RPh*RPh
      Sita_max = 85.0
      cos_sita = cos(Sita_max*PAI2/360.0)
      cos_sita2_min = cos_sita*cos_sita
      Itaus = 1
      Itaue = NtauExtraP

C     Check the input time-domain layer velocity.
C     open(21, file='velo_check.dat', access='direct',
C    +     recl=NX*NY)

      DO 5555 Itau = Itaus, Itaue

C        CALL INPUT_CURRENT_LAYER_VELOCITY(VVtau,
C    +        Nvline_first, Nvline_final,
C    +        Nvcdp_first, Nvcdp_final,
C    +        Nvmx, Nvmy,
C    +        Nline_first, Nline_final,
C    +        Ncdp_first, Ncdp_final,
C    +        NX, NY, Itau, IFN02,
C    +        Nmx0, Nmx1, Nmy0, Nmy1)
         DO IY = 1, NY
         DO IX = 1, NX
              VVtau(IX, IY) = vel3d(IX, IY, Itau)
         END DO
         END DO

C        WRITE(21,REC=Itau)((VVtau(IX,IY),IX=1,NX),IY=1,NY)

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
         ELSE
            write(*,*) 'Exit! Iphr=',Iphr,' Iww =',Iww,' Itau=',Itau
            RETURN
         END IF

5555  CONTINUE

C     close(21)
      RETURN
      END

*======================================================================*
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
C$OMP PARALLEL PRIVATE(VVtau2,cos_sita2,cos_sita,Tshift,PHASE)
C$OMP DO
      DO IY=1, NY
         DO IX=1, NX
            VVtau2         = VVtau(IX, IY)*VVtau(IX, IY)
            cos_sita2      = 1.0 - RPh2*VVtau2*0.25
            cos_sita       = SQRT(cos_sita2)
            Tshift         = Wflag*W*Dtau*cos_sita
            PHASE          = CMPLX(0.0, Tshift)
            DATW_W(IX, IY) = DATW_W(IX, IY)*CEXP(PHASE)
         END DO
      END DO
C$OMP END DO
C$OMP END PARALLEL
      RETURN

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

      RETURN
!     Compensate X-Y splitting error.
      VZ = 0
      DO Iy=1, NY
      DO Ix=1, NX
         VZ = VZ + VVtau(Ix,Iy)
      END DO
      END DO
      VZ = VZ/(NX*NY)
      cosalfa = sqrt(1.0 - 0.25*VZ*VZ*RPh2)
C     2022.02.08, No Compensation needed.
C     CALL COMPENSATE_PSFD(DATW_W, XX, XX_X, YY_Y, TMP2D,
C    +     NX, Dkx, Dx, NY, Dky, Dy,
C    +     Dtau, VZ, W, ALFA, Rkz_max, cosalfa)

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

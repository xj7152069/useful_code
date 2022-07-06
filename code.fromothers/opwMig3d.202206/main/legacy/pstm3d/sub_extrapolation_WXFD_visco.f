
      SUBROUTINE WXFD_EXTRAPOLATION_WXYCOMP_MPI_VISCO(
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

      INTEGER   Itau, weflag, NtauShallow
      INTEGER   Itaus, Itaue

C     Visco-Acoustic Features.
      INTEGER   ierror1, ierror2
      REAL      VEL_CONSTANT, Q_CONSTANT
      REAL      Wmain
      COMPLEX, ALLOCATABLE :: Qvalue(:,:)
      COMPLEX, ALLOCATABLE :: RKW(:,:)

C     FengBo add temp file saving&reading.
      CHARACTER*1024 FNW, FNR

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
      Itaus = 1
      Itaue = NtauExtraP

      DO 5555 Itau = Itaus, Itaue

         write(*,*) 'Itau=',Itau
         CALL INPUT_CURRENT_LAYER_VELOCITY(VVtau,
     +        Nvline_first, Nvline_final,
     +        Nvcdp_first, Nvcdp_final,
     +        Nvmx, Nvmy,
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
         ELSE
            write(*,*) 'Exit! Iphr=',Iphr,' Iww =',Iww,' Itau=',Itau
            RETURN
         END IF

5555  CONTINUE

      DEALLOCATE(Qvalue)
      DEALLOCATE(RKW)

      RETURN
      END

*=====================================================================*
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
      END SUBROUTINE COMPUTE_Q_BY_V

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
      END SUBROUTINE COMPUTE_Q_BY_V_txg
*======================================================================*

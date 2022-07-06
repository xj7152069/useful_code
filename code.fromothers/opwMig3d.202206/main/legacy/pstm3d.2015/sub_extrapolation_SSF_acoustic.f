
      SUBROUTINE ACOUSTIC_EXTRAPOLATION_SSF(
     +    myid, Iphr, RPh, W, DW, Iww,
     +    vel3d, NtauExtraP,
     +    DATW_W, DATW_R, VVtau,
     +    NX, Dx, NY, Dy, Ntau, Dtau,
     +    Dkx, Dky, Rkz_max,
     +    Nvline_first, Nvline_final, Nvmy,
     +    Nvcdp_first,  Nvcdp_final,  Nvmx,
     +    Nline_first, Nline_final, Nmy0, Nmy1,
     +    Ncdp_first, Ncdp_final, Nmx0, Nmx1)

      implicit none
!     Definition of input-parameters
      INTEGER   Nline_first, Nline_final, NY, Nmy0, Nmy1
      INTEGER   Ncdp_first,  Ncdp_final,  NX, Nmx0, Nmx1
      INTEGER   Nvline_first, Nvline_final, Nvmy
      INTEGER   Nvcdp_first,  Nvcdp_final,  Nvmx
      INTEGER   Ntau, NtauExtraP, myid, Iphr, Iww
      REAL      Dx, Dy, Dtau, W, DW, RPh, Dkx, Dky, Rkz_max
      REAL      VVtau(NX, NY), DATW_R(Ntau,Nmx1, Nmy1)
      COMPLEX   DATW_W(NX, NY)
      REAL      vel3d(NX, NY, Ntau)

C     Definitions of internal parameters.
      REAL      PAI, PAI2
      REAL      RPh2, Sita_max, cos_sita2_min
      REAL      VZ, VVtau2, cos_sita, cos_sita2, cosalfa
      INTEGER   weflag, Itaus, Itaue, Itau
      INTEGER   Ix, Iy, Imx, Imy
      COMPLEX   ImagUnit, Cimage

C     FengBo add temp file saving&reading.
      INTEGER  ierror1, ierror2, ierror3
      COMPLEX, ALLOCATABLE :: DATW_WK_SAVE(:,:)
      COMPLEX, ALLOCATABLE :: DATW_WK(:,:), TMP2D(:,:)

C     Init. of parameters.
      ALLOCATE(DATW_WK_SAVE(NX,NY), stat=ierror1)
      ALLOCATE(DATW_WK(NX,NY), stat=ierror2)
      ALLOCATE(TMP2D(NY,NX), stat=ierror3)

      ImagUnit = CMPLX(0.0,1.)
      PAI      = 3.14159265359
      PAI2     = 2.*PAI
      RPh2     = RPh*RPh
      Sita_max = 85.0
      cos_sita = cos(Sita_max*PAI2/360.0)
      cos_sita2_min = cos_sita*cos_sita
      Itaus = 1
      Itaue = NtauExtraP

      DO 5555 Itau = Itaus, Itaue

         DO Iy = 1, NY
         DO Ix = 1, NX
              VVtau(Ix, Iy) = vel3d(Ix, Iy, Itau)
         END DO
         END DO

C        V(z) is the average velocity.
         VZ = 0
         DO Iy=1, NY
         DO Ix=1, NX
             VZ = VZ + VVtau(Ix,Iy)
         END DO
         END DO
         VZ = VZ/(1.*NX*NY)
         cosalfa = sqrt(1.0 - 0.25*VZ*VZ*RPh2)

         weflag = 0
         DO Iy = 1, NY
         DO Ix = 1, NX
            VVtau2    = VVtau(Ix, Iy)*VVtau(Ix, Iy)
            cos_sita2 = 1.0 - RPh2*VVtau2*0.25
            IF (cos_sita2 .lt. cos_sita2_min) THEN
               weflag = 1
               GOTO 10724
            END IF
         END DO
         END DO

10724    CONTINUE
         IF (weflag .EQ. 0) THEN
            CALL SSF_EXTRAPOLATION_SINGLE_W(
     +           Itau, DATW_W, DATW_WK, TMP2D,
     +           NX, Dkx, Dx, NY, Dky, Dy,
     +           Dtau, VZ, W, DW, Rkz_max, cosalfa)

            DO Iy = 1, Nmy1
            DO Ix = 1, Nmx1
               Imy = Iy + Nmy0
               Imx = Ix + Nmx0
C              Cimage = DATW_W(Imx, Imy)/(-W*ImagUnit)
               Cimage = DATW_W(Imx, Imy)
               DATW_R(Itau,Ix,Iy) = DATW_R(Itau,Ix,Iy)
     +                            + REAL(Cimage)
C    +                            + REAL(DATW_W(Imx, Imy))
            END DO
            END DO
         ELSE
            write(*,*) 'Exit! Iphr=',Iphr,' Iww =',Iww,' Itau=',Itau
            RETURN
         END IF

5555  CONTINUE

      DEALLOCATE(DATW_WK_SAVE, DATW_WK, TMP2D)
C     DEALLOCATE(DATW_WK)
C     DEALLOCATE(TMP2D)

      RETURN
      END

*======================================================================*
      SUBROUTINE SSF_EXTRAPOLATION_SINGLE_W(
     +     Itau, DATW_WXY, DATW_WK, TMP2D,
     +     NX, Dkx, Dx, NY, Dky, Dy,
     +     Dtau, VZ, W, DW, Rkz_max, cosalfa)
      implicit none
C     Definitions of subroutine parameters.
      COMPLEX  DATW_WXY(NX,NY), DATW_WK_SAVE(NX,NY)
      COMPLEX  DATW_WK(NX,NY), TMP2D(NY,NX)
      INTEGER  NX, NY, Itau
      REAL     Dkx, Dky, Dx, Dy, Dtau, VZ, W, DW
      REAL     Rkz_max, cosalfa

C     Definitions of internal parameters.
      REAL     Wflag
      INTEGER  Iky, Ikx, Jkx, Jky
      REAL     a0, tmpa, tmpa2
      REAL     kz_tau, kr, kr2
      REAL     Rkx, Rkx2, Rky, Rky2
      COMPLEX  Phase

C     Init. of parameters.
      Wflag   = -1.0
      a0      = VZ*VZ*0.25/(W*W*cosalfa*cosalfa)

C     (1) Phase-Shift extrapolation using constant reference velocity V(z).
C     2D forward FFT, form u(x,y;w) to U(kx,Ky;w)
      IF ( 1 .eq. Itau ) THEN
          DO Iky=1, NY
              DO Ikx=1, NX
                  DATW_WK(Ikx, Iky)=DATW_WXY(Ikx, Iky)
              END DO
          END DO
          CALL FFT2D_OMP(DATW_WK, TMP2D,NX, NY, +1)
      ELSE
          DO Iky=1, NY
              DO Ikx=1, NX
                  DATW_WK(Ikx, Iky)=DATW_WK_SAVE(Ikx, Iky)
              END DO
          END DO
      END IF

      DO 6660 Iky  = 1, NY/2 + 1

         Jky    = NY + 2 - Iky
         Rky    = (Iky-1)*Dky
         Rky2   = Rky*Rky                           !Ky^2

         DO 6661 Ikx  = 1, NX/2 + 1

            Jkx    = NX + 2 - Ikx
            Rkx    = (Ikx-1)*Dkx
            Rkx2   = Rkx*Rkx                        !Kx^2

            kr2    = Rky2 + Rkx2                    !Kr^2=Kx^2+Ky^2
            tmpa2  = 1.0 - a0*kr2
            IF (tmpa2 .le. 0.0) THEN
               DATW_WK(Ikx, Iky) = 0.0
               IF (Ikx.GT.1) DATW_WK(Jkx,Iky)=0.0
               IF (Iky.GT.1) DATW_WK(Ikx,Jky)=0.0
               IF (Iky.GT.1.AND.Ikx.GT.1) DATW_WK(Jkx,Jky)=0.0
               goto 6661
            END IF

            tmpa   = sqrt(tmpa2)
C           IF (tmpa .le. Rkz_max) THEN
C              XX(Ikx, Iky) = 0.0
C              IF (Ikx.GT.1)              XX(Jkx,Iky)=0.0
C              IF (Iky.GT.1)              XX(Ikx,Jky)=0.0
C              IF (Iky.GT.1.AND.Ikx.GT.1) XX(Jkx,Jky)=0.0
C              goto 6661
C           END IF

            kz_tau       = W*cosalfa*tmpa
            Phase     = CEXP(CMPLX(0.0,Wflag*kz_tau*Dtau))
            DATW_WK(Ikx, Iky) = DATW_WK(Ikx, Iky)*Phase

            IF (Ikx.GT.1) DATW_WK(Jkx,Iky)=DATW_WK(Jkx,Iky)*Phase
            IF (Iky.GT.1) DATW_WK(Ikx,Jky)=DATW_WK(Ikx,Jky)*Phase
            IF (Iky.GT.1.AND.Ikx.GT.1) THEN
                DATW_WK(Jkx,Jky)=DATW_WK(Jkx,Jky)*Phase
            END IF

6661     CONTINUE

6660  CONTINUE

C     save the (kx-ky) domain wavefield.
      DO Iky=1, NY
          DO Ikx=1, NX
              DATW_WK_SAVE(Ikx, Iky)=DATW_WK(Ikx, Iky)
          END DO
      END DO

C     2D inverse FFT, form U(kx,Ky;w) to u(x,y;w)
      CALL FFT2D_OMP(DATW_WK, TMP2D,NX, NY, -1)
      DO Iky=1, NY
         DO Ikx=1, NX
            DATW_WXY(Ikx, Iky)=DATW_WK(Ikx, Iky)
         END DO
      END DO

C     (2) correct the horizontal velocity variation using SSF operator.
      RETURN
      END SUBROUTINE SSF_EXTRAPOLATION_SINGLE_W

*======================================================================*
      SUBROUTINE SSF_EXTRAPOLATION_COMPLEX_W(
     +     Itau, DATW_WXY, DATW_WK, TMP2D,
     +     NX, Dkx, Dx, NY, Dky, Dy,
     +     Dtau, VZ, W, DW, Rkz_max, cosalfa)
      implicit none
C     Definitions of subroutine parameters.
      COMPLEX  DATW_WXY(NX,NY), DATW_WK_SAVE(NX,NY)
      COMPLEX  DATW_WK(NX,NY), TMP2D(NY,NX)
      INTEGER  NX, NY, Itau
      REAL     Dkx, Dky, Dx, Dy, Dtau, VZ, W, DW
      REAL     Rkz_max, cosalfa

C     Definitions of internal parameters.
      REAL     Wflag
      INTEGER  Iky, Ikx, Jkx, Jky
      REAL     a0, tmpa, tmpa2
      REAL     kz_tau, kr, kr2
      REAL     Rkx, Rkx2, Rky, Rky2
      COMPLEX  Phase, CW
      COMPLEX  Ca0, Ctmpa, Ctmpa2, Ckz_tau
      COMPLEX  ImagUnit, Cimage

C     Init. of parameters.
      ImagUnit = CMPLX(0.0,1.)
      Wflag   = -1.0
      CW      = CMPLX(W,DW)
      a0      = VZ*VZ*0.25/(W*W*cosalfa*cosalfa)
      Ca0     = VZ*VZ*0.25/(CW*CW*cosalfa*cosalfa)

C     (1) Phase-Shift extrapolation using constant reference velocity V(z).
C     2D forward FFT, form u(x,y;w) to U(kx,Ky;w)
      IF ( 1 .eq. Itau ) THEN
          DO Iky=1, NY
              DO Ikx=1, NX
                  DATW_WK(Ikx, Iky)=DATW_WXY(Ikx, Iky)
              END DO
          END DO
          CALL FFT2D_OMP(DATW_WK, TMP2D,NX, NY, +1)
      ELSE
          DO Iky=1, NY
              DO Ikx=1, NX
                  DATW_WK(Ikx, Iky)=DATW_WK_SAVE(Ikx, Iky)
              END DO
          END DO
      END IF

      DO 6660 Iky  = 1, NY/2 + 1

         Jky    = NY + 2 - Iky
         Rky    = (Iky-1)*Dky
         Rky2   = Rky*Rky                           !Ky^2

         DO 6661 Ikx  = 1, NX/2 + 1

            Jkx    = NX + 2 - Ikx
            Rkx    = (Ikx-1)*Dkx
            Rkx2   = Rkx*Rkx                        !Kx^2

            kr2    = Rky2 + Rkx2                    !Kr^2=Kx^2+Ky^2
            tmpa2  = 1.0 - a0*kr2
            Ctmpa2  = 1.0 - Ca0*kr2
C           IF (tmpa2 .le. 0.0) THEN
C              DATW_WK(Ikx, Iky) = 0.0
C              IF (Ikx.GT.1) DATW_WK(Jkx,Iky)=0.0
C              IF (Iky.GT.1) DATW_WK(Ikx,Jky)=0.0
C              IF (Iky.GT.1.AND.Ikx.GT.1) DATW_WK(Jkx,Jky)=0.0
C              goto 6661
C           END IF

            tmpa   = sqrt(tmpa2)
            Ctmpa  = csqrt(Ctmpa2)
C           IF (tmpa .le. Rkz_max) THEN
C              XX(Ikx, Iky) = 0.0
C              IF (Ikx.GT.1)              XX(Jkx,Iky)=0.0
C              IF (Iky.GT.1)              XX(Ikx,Jky)=0.0
C              IF (Iky.GT.1.AND.Ikx.GT.1) XX(Jkx,Jky)=0.0
C              goto 6661
C           END IF

            kz_tau       = W*cosalfa*tmpa
            Phase     = CEXP(CMPLX(0.0,Wflag*kz_tau*Dtau))
C           Ckz_tau      = CW*cosalfa*Ctmpa
C           Phase     = CEXP(Wflag*ImagUnit*Ckz_tau*Dtau)
            DATW_WK(Ikx, Iky) = DATW_WK(Ikx, Iky)*Phase

            IF (Ikx.GT.1) DATW_WK(Jkx,Iky)=DATW_WK(Jkx,Iky)*Phase
            IF (Iky.GT.1) DATW_WK(Ikx,Jky)=DATW_WK(Ikx,Jky)*Phase
            IF (Iky.GT.1.AND.Ikx.GT.1) THEN
                DATW_WK(Jkx,Jky)=DATW_WK(Jkx,Jky)*Phase
            END IF

6661     CONTINUE

6660  CONTINUE

C     save the (kx-ky) domain wavefield.
      DO Iky=1, NY
          DO Ikx=1, NX
              DATW_WK_SAVE(Ikx, Iky)=DATW_WK(Ikx, Iky)
          END DO
      END DO

C     2D inverse FFT, form U(kx,Ky;w) to u(x,y;w)
      CALL FFT2D_OMP(DATW_WK, TMP2D,NX, NY, -1)
      DO Iky=1, NY
         DO Ikx=1, NX
            DATW_WXY(Ikx, Iky)=DATW_WK(Ikx, Iky)
         END DO
      END DO

C     (2) correct the horizontal velocity variation using SSF operator.
      RETURN
      END SUBROUTINE SSF_EXTRAPOLATION_COMPLEX_W

*======================================================================*

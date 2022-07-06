SUBROUTINE READ_OPW_DATA(&
        FN1, FN4, IFN01, IFN04, iDebug, DATW,&
        Nmx, Nmy, NW1, NW2, NW3, NW4, NW, LT,&
        Ncdp_first, Ncdp_final, Ncdp_step,&
        Nline_first, Nline_final, Nline_step,&
        Iphr, Nmx0, Nmx1, Nmy0, Nmy1)

    implicit none
    CHARACTER*1024 FN1, FN4
    COMPLEX DATW(Nmx, Nmy, NW)
    INTEGER Nmx, Nmy, NW1, NW2, NW3, NW4, NW, LT
    INTEGER Ncdp_first, Ncdp_final, Ncdp_step
    INTEGER Nline_first, Nline_final, Nline_step
    INTEGER Iphr, Nmx0, Nmx1, Nmy0, Nmy1
    INTEGER IFN01, IFN04, iDebug

    ! internal buffers.
    REAL,   ALLOCATABLE :: Pdat2d(:,:)
    COMPLEX,ALLOCATABLE :: CTMPTRACE2d(:,:)

    ! internal parameters.
    INTEGER ierr1, ierr2, ncdp
    INTEGER irec
    INTEGER iw, imy, imx, Iline, Icdp
    INTEGER it, ix
    CHARACTER*1024 FN_opw
    CHARACTER*4    ciphr
    REAL    ampAbsMax

    ! for Luojia data only, 2022.02.17.
    ampAbsMax = 1E7

    ncdp  = Ncdp_final-Ncdp_first+1
    ALLOCATE(CTMPTRACE2d(LT,ncdp),stat=ierr1)
    ALLOCATE(Pdat2d(LT,ncdp),     stat=ierr2)

!    WRITE(*,*) 'test of FFT1D, LT=',LT
!    DO Icdp = Ncdp_first, Ncdp_final
!
!        Imx = Icdp - Ncdp_first + 1 + Nmx0
!        IX  = Icdp - Ncdp_first + 1
!        Imy = 1
!        WRITE(*,*) 'Icdp=', Icdp, 'IX=', IX
!
!        DO IT = 1, LT
!            CTMPTRACE2d(IT,IX)=0
!        END DO
!
!        !LT=2001
!        CALL FFT1D(CTMPTRACE2d(1,IX), LT, 1)
!
!        DO Iw = 1, NW
!            DATW(Imx, Imy, Iw)=CTMPTRACE2d(NW1+Iw,IX)
!        END DO
!    END DO

    DO iw  = 1, NW
        DO imy = 1, Nmy
            DO imx = 1, Nmx
                DATW(imx,imy,iw) = CMPLX(0.0,0.0)
            END DO
        END DO
    END DO

    ! 2022.02, FengBo use fortran IO function to access 4D-OPW data.
    ! open the plane-wave gather file: FN1.
    OPEN(IFN01, file=trim(adjustl(FN1)), access='direct', recl=LT*ncdp)

    WRITE(*,*) 'Iphr=', Iphr, 'LT=', LT, 'ncdp=', ncdp
    WRITE(ciphr,'(I4)') Iphr
    FN_opw=trim(adjustl(FN4))//'opwData_Iphr_'//trim(adjustl(ciphr))//'.dat'
    IF (iDebug .eq. 1) THEN
        WRITE(*,*) 'open file: ',trim(adjustl(FN_opw)),' to write.'
        OPEN(IFN04, file=FN_opw, access='direct', recl=LT*ncdp)
    END IF

    DO Iline = Nline_first, Nline_final

        Imy = Iline - Nline_first + 1 + Nmy0
        if( mod(Iline,50) .eq. 0 ) WRITE(*,*) 'Iline= ', Iline

        irec=(Iphr-1)*(Nline_final-Nline_first+1)+(Iline-Nline_first+1)
        !WRITE(*,*) 'Iline= ', Iline, 'irec=', irec
        READ(IFN01, rec=irec) ((Pdat2d(it,ix),it=1,LT),ix=1,ncdp)

        do ix=1,ncdp
            do it=1,lT
                if(isnan(Pdat2d(it,ix))) then
                    Pdat2d(it,ix) = 0
                    WRITE(*,*) 'Iline=',Iline,'ix=',ix,'it=',it
                endif
                !if(ABS(Pdat2d(it,ix))>ampAbsMax) then
                !    WRITE(*,*) 'Amp=', ABS(Pdat2d(it,ix)), 'Iline=',Iline,'ix=',ix,'it=',it
                !    Pdat2d(it,ix) = 0
                !endif
            enddo
        enddo
        IF (iDebug .eq. 1) THEN
            WRITE(IFN04, rec=Iline-Nline_first+1) ((Pdat2d(it,ix),it=1,LT),ix=1,ncdp)
        END IF

        !$OMP PARALLEL PRIVATE(Icdp,Imx,IX,IT,Iw,Imy)
        !$OMP DO
        DO Icdp = Ncdp_first, Ncdp_final

            Imx = Icdp - Ncdp_first + 1 + Nmx0
            IX  = Icdp - Ncdp_first + 1
            !WRITE(*,*) 'Iline= ', Iline, 'Icdp=', Icdp

            DO IT = 1, LT
                CTMPTRACE2d(IT,IX)=CMPLX(Pdat2d(IT,IX), 0.0)
            END DO

            CALL FFT1D(CTMPTRACE2d(1,IX), LT, 1)

            DO Iw = 1, NW
                DATW(Imx, Imy, Iw)=CTMPTRACE2d(NW1+Iw,IX)
!               DATW(Imx, Imy, Iw)=CTMPTRACE2d(NW1+Iw-1,IX)
            END DO
        END DO
        !$OMP END DO
        !$OMP END PARALLEL
!        goto 3456
!3456    continue

    END DO !end do iline

    IF (iDebug .eq. 1) THEN
        CLOSE(IFN04)
    END IF
    CLOSE(IFN01)

    DEALLOCATE(Pdat2d)
    DEALLOCATE(CTMPTRACE2d)

    RETURN
END SUBROUTINE READ_OPW_DATA

!======================================================================*

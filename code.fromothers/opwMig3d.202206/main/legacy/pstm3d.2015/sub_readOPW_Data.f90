SUBROUTINE READ_PLANEWAVE_DATA_AND_FFT_Profile2d(&
        FN1, Trace, DATW, CTrace,&
        Nmx, Nmy, NW1, NW2, NW3, NW4, NW, LT, DW, &
        Ncdp_first, Ncdp_final, Ncdp_step,&
        Nline_first, Nline_final, Nline_step,&
        Nphr, Iphr, IFN01, Nmx0, Nmx1, Nmy0, Nmy1)

    implicit none
    CHARACTER*1024 FN1
    REAL    Trace(LT)
    COMPLEX CTrace(LT)
    COMPLEX DATW(Nmx, Nmy, NW)
    INTEGER Nmx, Nmy, NW1, NW2, NW3, NW4, NW, LT
    REAL    DW
    INTEGER Ncdp_first, Ncdp_final, Ncdp_step
    INTEGER Nline_first, Nline_final, Nline_step
    INTEGER Nphr, Iphr, IFN01, Nmx0, Nmx1, Nmy0, Nmy1

    REAL,   ALLOCATABLE :: Pdat2d(:,:)
    COMPLEX,ALLOCATABLE :: CTMPTRACE2d(:,:)
    INTEGER ierr1, ierr2, ncdp
    INTEGER irec
    INTEGER iw, imy, imx, Iline, iix, iit, Icdp
    INTEGER it, ix

    ncdp  = Ncdp_final-Ncdp_first+1
    ALLOCATE(CTMPTRACE2d(LT,ncdp),stat=ierr1)
    ALLOCATE(Pdat2d(LT,ncdp),     stat=ierr2)

    DO iw  = 1, NW
        DO imy = 1, Nmy
            DO imx = 1, Nmx
                DATW(imx,imy,iw) = CMPLX(0.0,0.0)
            END DO
        END DO
    END DO

    ! 2022.02, FengBo use fortran IO function to access 4D-OPW data.
    ! open the plane-wave gather file: FN1.
    open(IFN01, file=trim(adjustl(FN1)), access='direct', recl=LT*ncdp)

    WRITE(*,*) 'LT=',LT,'ncdp=',ncdp
!   open(78, file='./comPhr.dat',access='direct', recl=LT*ncdp)

    DO Iline = Nline_first, Nline_final

        Imy = Iline - Nline_first + 1 + Nmy0
        if( mod(Iline,10) .eq. 0 ) write(*,*) 'Iline= ', Iline

        irec=(Iphr-1)*(Nline_final-Nline_first+1)+(Iline-Nline_first+1)
        read(IFN01, rec=irec) ((Pdat2d(it,ix),it=1,LT),ix=1,ncdp)

        do iix=1,ncdp
            do iit=1,lT
                if(isnan(Pdat2d(iit,iix))) then
                    Pdat2d(iit,iix) = 0
                    WRITE(*,*) 'Iline=',Iline,'iix=',iix,'iit=',iit
                endif
            enddo
        enddo
!       write(78, rec=Iline-Nline_first+1) ((Pdat2d(it,ix),it=1,LT),ix=1,ncdp)

        !$OMP PARALLEL PRIVATE(Imx,IX,IT,Iw)
        !$OMP DO
        DO Icdp = Ncdp_first, Ncdp_final

            Imx = Icdp - Ncdp_first + 1 + Nmx0
            IX  = Icdp - Ncdp_first + 1

            DO IT = 1, LT
                CTMPTRACE2d(IT,IX)=CMPLX(Pdat2d(IT,IX), 0.0)
            END DO

            CALL FFT1D(CTMPTRACE2d(1,IX), LT, 1)

            DO Iw = 1, NW
                DATW(Imx, Imy, Iw)=CTMPTRACE2d(NW1+Iw-1,IX)
            END DO
        END DO
        !$OMP END DO
        !$OMP END PARALLEL

    END DO !end do iline

!   close(78)
    close(IFN01)

    DEALLOCATE(Pdat2d)
    DEALLOCATE(CTMPTRACE2d)

    RETURN
END SUBROUTINE READ_PLANEWAVE_DATA_AND_FFT_Profile2d

!======================================================================*

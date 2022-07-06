SUBROUTINE readParMig3d(FN_PAR, FN1, FN2, FN3, FN4,&
      LT, DT, Ncdp_first, Ncdp_final, Ncdp_step, Nmx, Dmx,&
      Nline_first, Nline_final, Nline_step, Nmy, Dmy,&
      Nphr, Phrmin, Dphr,&
      Nvtau, Dvtau,&
      Nvcdp_first, Nvcdp_final, Nvcdp_step, Nvmx, Dvx,&
      Nvline_first, Nvline_final, Nvline_step, Nvmy, Dvy,&
      ViscoFlag, NtauExtraP, F1, F2, F3, F4,&
      Iphr1, Iphr2)

      implicit none
      CHARACTER*1024 FN_PAR, FN1, FN2, FN3, FN4, CTMP
      REAL    F1, F2, F3, F4
      INTEGER LT, Ncdp_first, Ncdp_final, Ncdp_step, Nmx
      INTEGER Nline_first, Nline_final, Nline_step, Nmy
      INTEGER Nvtau
      INTEGER Nvcdp_first, Nvcdp_final, Nvcdp_step, Nvmx
      INTEGER Nvline_first, Nvline_final, Nvline_step, Nvmy
      REAL    DT, Dmx, Dmy, Dvtau, Dvx, Dvy, Phrmin, Dphr
      INTEGER Nphr, Iphr1, Iphr2, ViscoFlag, NtauExtraP

      OPEN(13, FILE=FN_PAR)
      READ (13,'(A)') FN1              ! Offset-planewave decompostion data.
      READ (13,'(A)') FN2              ! Interval velocity field(time domain).
      READ (13,'(A)') FN3              ! the stacked image.
      READ (13,'(A)') FN4              ! Common-ph imaging volume, Ntau*Ncdp*Nline*Nphr
      ! parameters of the offset plane-wave data.
      READ (13,*) CTMP, LT             ! for interval velocity model.
      READ (13,*) CTMP, DT             ! unit of DT is ms.
      READ (13,*) CTMP, Ncdp_first     ! minimum cdp number of the project.
      READ (13,*) CTMP, Ncdp_final     ! maximum cdp number of the project.
      READ (13,*) CTMP, Dmx            ! the cdp spacing (m)
      Nmx = Ncdp_final - Ncdp_first + 1
      READ (13,*) CTMP, Nline_first    ! minimum line number of the project.
      READ (13,*) CTMP, Nline_final    ! maximum line number of the project.
      READ (13,*) CTMP, Dmy            ! the line spacing (m)
      Nmy = Nline_final - Nline_first + 1
      READ (13,*) CTMP, Nphr           ! number of offset plane-wave decomposition
      READ (13,*) CTMP, Phrmin         ! minimum Phr value (us/m).
      READ (13,*) CTMP, Dphr           ! Phr interval (us/m)
      Phrmin = Phrmin/1.0E6 
      Dphr   = Dphr/1.0E6
      ! parameters of the time-domain interval velocity model.
      READ (13,*) CTMP, Nvtau          ! for interval velocity model.
      READ (13,*) CTMP, Dvtau          ! unit of Dvtau is ms.
      READ (13,*) CTMP, Nvcdp_first    ! minimum cdp number of the velocity model.
      READ (13,*) CTMP, Nvcdp_final    ! maximum cdp number of the velocity model.
      READ (13,*) CTMP, Dvx            ! the cdp spacing (m)
      Nvmx = Nvcdp_final - Nvcdp_first + 1
      READ (13,*) CTMP, Nvline_first   ! minimum line number of the velocity model.
      READ (13,*) CTMP, Nvline_final   ! maximum line number of the velocity model.
      READ (13,*) CTMP, Dvy            ! the line spacing (m)
      Nvmy = Nvline_final - Nvline_first + 1
      ! parameters for the 3-D migration.
      READ (13,*) CTMP, NtauExtraP     ! Extrapolation number.
      READ (13,*) CTMP, F1, F2, F3, F4 ! The frequency range for imaging.
      READ (13,*) CTMP, Iphr1          ! The first Ray-parameter number.
      READ (13,*) CTMP, Iphr2          ! The last  Ray-parameter number.
      CLOSE(13)

      ! Init. of parameters.
      ViscoFlag   = 0
      Nvcdp_step  = 1
      Nvline_step = 1
      Ncdp_step   = 1
      Nline_step  = 1

      RETURN

END SUBROUTINE readParMig3d

!======================================================================*

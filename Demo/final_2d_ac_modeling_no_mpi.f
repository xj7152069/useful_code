*********************************************************************
*                         AC2DMOD UNDER TOPO                        *
*********************************************************************
*              MARINE GEOLOGY AND GEOPHISICS DEPARTMENT, 	    *
*                        TONG JI UNIVERSITY                         *
*                             2005.7.19       			    *
*                       All right reserved                          *
*********************************************************************
       PROGRAM MODEL2D
       !include "mpif.h"

       PARAMETER(Lmax=530000000)
       DIMENSION buf(Lmax)
       CHARACTER*256 FNAME, VEL_FNAME, FN3

        INTEGER NTRACE,NMAX,NMAZ,NSHOT,NZ,LT,rdt
        REAL    DX,DZ,XS_START,DX_SHOT,DX_TRACE,dx_trace_new
        REAL    XG_START
        REAL    dt_new 

*==================================================================
       integer myid,np,ierr

       !call mpi_init(ierr)

        CALL READPAR(VEL_FNAME, FN3, 
     +        NSHOT, NTRACE, NMAX, NMAZ, 
     +        DX, DZ, XS_START, DS_FR_FS, 
     +        DX_SHOT, DX_TRACE, DT, LT, 
     +        PEAK, ipc, dt_new)

c       pause'check paras'

        dx_trace_new=dx_trace
        IF(dx_trace_new.lt.dx_trace) THEN
        write(*,*)'dx_trace_new can not be less than dt_new'
        write(*,*)'stop and modify dx_trace_new !'
        stop
        ENDIF

        IF(dt_new.lt.dt) THEN
        write(*,*)'dt_new can not be less than dt'
        write(*,*)'stop and modify dt_new !'
        stop
        ENDIF

        NZ=NMAZ
        na=20
        ntrace=ntrace+2*na
        xg_start=DS_FR_FS-na*dx_trace

        rdt=dt_new/DT
        lt_new=LT/rdt
        write(*,*) 'dt_new,DT,LT,rdt,lt_new'
        write(*,*)  dt_new,DT,LT,rdt,lt_new

        kp=1
        kboud=2
        f0=25
        t0=40

        t0=t0/1000.
        dt=dt/1000.

        x=(ntrace-1)*dx_trace
        nx=INT(x/dx+0.5)+1

        x2=(na-1)*dx_trace
        nx2=INT(x2/dx+0.5)+1
      
***********************************************************************
*       20 GRID LINES ARE ADDED ABOVE THE TOP OF THE VELOCITY MODEL
*       IN ORDER TO USING CHARASTERISTIC ABSORBING BOUNDARIES
*       BOTH THE SHOTS AND RECEIVERS LOCATE AT THE 21 GIED_LINE
*********************************************************************** 
	nz=nmaz+20
****************************************************************
*       THE DELAYED TIME OF THE WAVELET MUST BE INCLUDED 
*       IN THE CALCULATION TIME
**************************************************************** 

       lt=lt+int(t0/dt+0.5)
       IJ=INT(T0/DT+0.5)

       L1=1
       L2=L1+5
       L3=L2+nx*nz
       L4=L3+(5+NX+5)*(5+NZ+5)
       L5=L4+(5+NX+5)*(5+NZ+5)
       L6=L5+(5+NX+5)*(5+NZ+5)
       L7=L6+ntrace*lt
       L8=L7+nmax
       L9=L8+nz

       write(*,*)"*****************************"
       write(*,*)"Memomery left is", Lmax-L9
       write(*,*)"Memomery for last buf is",ntrace*lt 
       write(*,*)"*****************************"

       if((ntrace*lt).gt.(Lmax-L9)) THEN
       write(*,*)"Memomery left is not enough for last buf!"
       write(*,*)"Please applying for more memomery !"
       stop
       endif

       CALL multi_mod(buf(L1),buf(L2),buf(L3),buf(L4),
     +         buf(L5),buf(L6),buf(L7),buf(L8),buf(L9),
     +          kboud,kp,f0,t0,nmax,nmaz,nx,nx2,nz,dx,dz,
     +          xs_start,dx_shot,nshot,xg_start,dx_trace,
     +          ntrace,dt,lt,vel_fname,ij,na,fn3,
     +          peak,ipc,lt_new,rdt,dx_trace_new)

       !call mpi_finalize(ierr) 
       END
*=====================================================================
        subroutine multi_mod(C,v,u,u1,u2,seism,whole_stat_elev,v1,
     +             seism_out,kboud,kp,f0,t0,nmax,nmaz,nx,nx2,
     +             nz,dx,dz,xs_start,
     +             dx_shot,nshot,xg_start,dx_trace,ntrace,dt,lt,
     +             vel_fname,ij,na,fn3,peak,ipc,lt_new,rdt,dx_trace_new)
        
        integer kboud,kp,nmax,nmaz,nz,nshot,ntrace,lt,lt_new,rdt,ng_left
        integer ng_right,iw,ncx_shot,ncz_shot,deltax
        integer ierr,myid,np
	INTEGER NSPX_SHOT,NSPZ_SHOT,ij,na
	real f0,t0,dx,dz,dx_shot,xs_start,xg_start,dx_trace,dx_trace_new,dt
        !include "mpif.h"
        DIMENSION  C(5)
	DIMENSION  V(NX, NZ),V1(NZ)
        DIMENSION  u(-4:nx+5, -4:nz+5)
        DIMENSION  u1(-4:nx+5, -4:nz+5)
        DIMENSION  u2(-4:nx+5, -4:nz+5)
        DIMENSION  seism(ntrace, lt)
        DIMENSION  seism_out(ntrace, lt)
        DIMENSION  whole_stat_eleV(nmax)

        DIMENSION  TMPVV(17913, 7465)

        CHARACTER*256 VEL_FNAME, FN3
        !call mpi_comm_rank(MPI_COMM_WORLD,myid,ierr)
        !call mpi_comm_size(MPI_COMM_WORLD,np,ierr)
        write(*,*)'np= ', np, 'myid=', myid

        if(myid.eq.0) then
        write(*,*)'vel_fname,fn3'
        write(*,*) vel_fname,fn3
        write(*,*)'kboud,kp,f0,t0,nmax,nmaz,nx'
        write(*,*) kboud,kp,f0,t0,nmax,nmaz,nx 
        write(*,*)'nmax,nmaz,nx'
        write(*,*) nmax,nmaz,nx 
        write(*,*)'nz,dx,dz,xs_start'
        write(*,*) nz,dx,dz,xs_start
        write(*,*)'dx_shot,nshot,xg_start'
        write(*,*) dx_shot,nshot,xg_start
        write(*,*)'dx_trace,ntrace,dt,lt'
        write(*,*) dx_trace,ntrace,dt,lt
        write(*,*)'ij,na,peak'
        write(*,*) ij,na,peak
        write(*,*)'ipc,lt_new'
        write(*,*) ipc,lt_new
        endif

c       pause
 
        OPEN(11, FILE=vel_fname, ACCESS='DIRECT', RECL=ipc*nmaz) 
        OPEN(12, FILE=fn3)

         if(myid.eq.0) write(*,*)'PEAK',peak
          do i = 1,nmax
            read(12,*)a,elev
            whole_stat_elev(i) = int((peak-elev)/dz+0.5)+1
c           whole_stat_elev(i) = 1
          enddo
        CLOSE(12)

*====== VIP:CHECK IF THE SHOT OR RECEIVER ARE LOCATED IN THE AIR
c       OPEN(91, FILE='tmp', ACCESS='DIRECT', RECL=ipc*nmaz)
c         do ix = 1,nmax
c           read(11,rec=ix)(TMPVV(ix,iz),iz=1,nmaz)
c         enddo
c         do ix = 1,nmax
c           kkz=whole_stat_elev(ix)
c           write(*,*)ix,'kkz=',kkz
c           TMPVV(ix,kkz)=8000.
c         enddo
c         do ix = 1,nmax
c           write(91,rec=ix)(TMPVV(ix,iz),iz=1,nmaz)
c         enddo
c       CLOSE(11)
c       CLOSE(91)
c       STOP'tmp'

*====== NOW THE SHOT LOOP BEGINING......
*       write(*,*)" NOW THE SHOT LOOP BEGINING......"

        !DO 2222 IS=1+myid, nshot, np
        DO 2222 IS=1, nshot
 
        !if(myid.eq.0)then
        write(*,*) 'Myid=', myid, 'Is=', Is
        write(*,*) 'shot No.=== ', is
        write(*,*) 'xs_start ', xs_start
        write(*,*) 'xg_start ', xg_start
        !endif

*===================================================================
*       ZERO THE WORKING BUFFERS
*===================================================================
 
        CALL ZERO_BUF(NX, NZ, NTRACE, LT, U, U1, U2, SEISM, V) 

*-------------------------------------------------------------------
*       LATERAL GRID NO. OF THE CURRENT SHOT
*-------------------------------------------------------------------

	NSPX_SHOT=INT((xs_start+(is-1)*dx_shot)/dx+1)

	IW=int(dx_trace/dx+0.5)
*-------------------------------------------------------------------
*       LEFT GRID NO. OF A SHOT
*-------------------------------------------------------------------
	ng_left=INT(xg_start/dx-0.5)+NSPX_SHOT
*-------------------------------------------------------------------
*       RIGHT GRID NO. OF A SHOT
*-------------------------------------------------------------------
	ng_right=ng_left+nx-1

        if(myid.eq.0)write(*,*)'nspx_shot =',nspx_shot
c       if(myid.eq.0)write(*,*)'nspz_shot =',whole_stat_elev(nspx_shot)
        if(myid.eq.0)write(*,*)'ng_left   =',ng_left
        if(myid.eq.0)write(*,*)'ng_right  =',ng_right
c       pause

*=======get the velocity of the current shot

        CALL VEL_SHOT(V,V1,NMAX,NMAZ,NX,nx2,NZ,NG_LEFT,NG_RIGHT,myid,IS)

*...... SINGAL SHOT MODELING BEGINING......

*-------SHOT POINT POSITION IN A SINGLE SHOT VELOCITY FIELD

        NCX_SHOT=int(ABS(xg_start)/dx+0.5)+1
        NCZ_SHOT=int(whole_stat_elev(NCX_SHOT+ng_left-1))+20
       
c       if(myid.eq.0)write(*,*)'NCX_SHOT=',NCX_SHOT
c       if(myid.eq.0)write(*,*)'NCZ_SHOT=',NCZ_SHOT
c       pause

        CALL SINGLE_SHOT_MOD(C,V, NX, NZ, U, U1, U2, SEISM, IRZ, 
     +           NTRACE, DX, DZ, NCX_SHOT, NCZ_SHOT, LT, DT, 
     +           kboud,kp,f0,t0,iw,is,whole_stat_elev,nmax,ng_left)

*====== output the calculated result of the current shot
      
c        OPEN(21, FILE='output', ACCESS='DIRECT', RECL=ipc*(LT-ij))
c           DO IX=1, NTRACE-2*NA
c           KREC=(IS-1)*(NTRACE-2*NA)+IX
c           WRITE(21, REC=KREC) (SEISM(IX+NA, IT), IT=IJ+1, LT)
c           END DO
c        CLOSE(21)

      CALL WRITE_DISK(IS,NTRACE,LT,lt_new,SEISM,T0,DT,IJ,
     +   NA,rdt,dx_trace_new,SEISM_OUT,dx_trace,ipc)
 
2222    CONTINUE
 
*=======PROGRAM FINISHED!
    
c       CLOSE(10)
        CLOSE(11)

        RETURN
        END

*===================================================================
*       THE WORKING BUFFERS AER SET TO ZEROS
*===================================================================

       SUBROUTINE ZERO_BUF(NX, NZ, NTRACE, LT, U, U1, U2, SEISM, V)

       INTEGER NX, NZ, NTRACE, LT

       DIMENSION V(NX, NZ)
       DIMENSION SEISM(NTRACE, LT)
       DIMENSION U(-4:NX+5, -4:NZ+5)
       DIMENSION U1(-4:NX+5, -4:NZ+5)
       DIMENSION U2(-4:NX+5, -4:NZ+5)

       DO IX=1, NTRACE
          DO IT=1, LT
             SEISM(IX, IT)=0.0
          END DO
       END DO

       DO IX=1, NX
          DO IZ=1, NZ
             V(IX, IZ)=0.0
          END DO
       END DO

       DO IX=-4, NX+5
          DO IZ=-4, NZ+5
             U(IX, IZ)=0.0
             U1(IX, IZ)=0.0
             U2(IX, IZ)=0.0
          END DO
       END DO

       RETURN
       END

*===================================================================
*
*      GET THE VELOCITY FIELD OF THE CURRENT SHOT
*
*===================================================================

       SUBROUTINE VEL_SHOT(V,V1,NMAX,NMAZ,NX,nx2,NZ,
     +       NG_LEFT,NG_RIGHT,myid,IS)
       INTEGER   NMAX, NMAZ, NX, nx2, NZ
       INTEGER   NG_LEFT,NG_RIGHT 
 
       DIMENSION V1(NZ),V(NX, NZ)

c      character*256 vel_fname
 
       IF(NG_LEFT.lt.1) THEN

       if(myid.eq.0)then
         write(*,*) 'The left boundary of the current shot is'
         write(*,*) 'beyound the left boundary'
         write(*,*) 'of the velocity field'
         write(*,*) 'the first trace in velocity field will be'
         write(*,*) 'used to fill the negative part'
c        write(*,*) 'the calculation stop here !!'
c        stop
       endif
       END IF

       IF(NG_RIGHT.GT.NMAX) THEN
         write(*,*) 'The right boundary of the current shot is'
         write(*,*) 'beyound the rightest boundary'
         write(*,*) 'of the velocity field'
         write(*,*) 'the last trace in velocity field will be'
         write(*,*) 'used to fill the excess part'
c        write(*,*) 'the calculation stop here !!'
c        stop
       END IF

       IF((NG_LEFT.lt.1).and.(NG_RIGHT.GT.NMAX)) THEN
         write(*,*) 'The left and right boundary of the current shot is'
         write(*,*) 'beyound the left and right boundary'
         write(*,*) 'of the velocity field simultaneously'
         write(*,*) 'Attention !!!'
         !stop
       END IF

       IF(NG_LEFT.gt.NMAX) THEN
         write(*,*) 'The left boundary of the current shot is'
         write(*,*) 'beyound the right boundary'
         write(*,*) 'of the velocity field'
         write(*,*) 'You are totally wrong !!!'
         stop
       END IF

        DO IX=1,NX
        DO IZ=1,NZ
           V(IX,IZ)=0.
        END DO
        END DO

       NNN=NX-ng_right
       write(*,*) 'ng_left=',ng_left
       write(*,*) 'ng_right=',ng_right
       write(*,*) 'NNN=',NNN

******** WHEN ng_left < 1, PROCESSING AS FOLLOWS:
        IF(ng_left.lt.1) THEN
        write(*,*) 'ng_left=',ng_left
          DO IX=1,ng_right
           read(11,rec=IX)(v1(j),j=1,nmaz)
           DO IZ=1, NZ
             IF(IZ.LE.20) THEN
               V(IX+NNN,IZ)=V1(1)
             ELSE
               V(IX+NNN,IZ)=V1(IZ-20)
             ENDIF
           END DO
          END DO
******** CP V(NNN+1) TO V(1:NNN) 
        DO IX=1,NNN
        DO IZ=1,NZ
           V(IX,IZ)=V(NNN+1,IZ)
        END DO
        END DO

        ENDIF
******** END OF PROCESSING WHEN ng_left < 1

******** WHEN ng_left > 1 and ng_right < NMAX, PROCESSING AS FOLLOWS:
        IF(ng_left.ge.1.and.ng_right.le.NMAX) THEN
          II=0
          DO IX=ng_left,ng_right
           II=II+1
           read(11,rec=IX)(v1(j),j=1,nmaz)
           DO IZ=1, NZ
             IF(IZ.LE.20) THEN
               V(II,IZ)=V1(1)
             ELSE
               V(II,IZ)=V1(IZ-20)
             ENDIF
           END DO
          END DO
        ENDIF
******** END OF PROCESSING WHEN ng_left > 1 and ng_right < NX 

******** WHEN ng_right < NX, PROCESSING AS FOLLOWS:
        IF(ng_right.gt.NMAX) THEN
c         write(*,*) 'ng_right=',ng_right
          II=0
          DO IX=ng_left,NMAX
           II=II+1
           read(11,rec=IX)(v1(j),j=1,nmaz)
           DO IZ=1, NZ
             IF(IZ.LE.20) THEN
               V(II,IZ)=V1(1)
             ELSE
               V(II,IZ)=V1(IZ-20)
             ENDIF
           END DO
          END DO

c         write(*,*) 'II=',II
******** FILL ZERO TRACE WITH VTMP(NMAX,NZ) 
         DO IX=II+1,NX
         DO IZ=1,NZ
         V(IX,IZ)=V(II,IZ)
         END DO
         END DO
        ENDIF
******** END OF PROCESSING WHEN ng_left < 1
c       OPEN(33,FILE='velocity',access='direct',RECL=1*NZ)
c       write(*,*) '*** start and end record when write disk='
c       write(*,*) (is-1)*nx+1,is*nx
c       ii=0
c       do ix=(is-1)*nx+1,is*nx
c       ii=ii+1
c       write(33,rec=ix) (v(ii,iz),iz=1,NZ)
c       end do
c       close(33)
c       write(*,*)'velocity',NZ,NX

        RETURN
        END

*========================================================
*      SINGAL SHOT MODELING SUBROUTINE
*========================================================

       SUBROUTINE SINGLE_SHOT_MOD(C,V, NX, NZ, U, U1, U2, SEISM, IRZ,
     +            NTRACE, DX, DZ, NCX_SHOT, NCZ_SHOT, LT, DT, 
     +            kboud,kp,f0,t0,iw,is, whole_stat_elev,nmax,ng_left)

       REAL      DT,DX,DZ,f0,t0 
       INTEGER   NTRACE,LT,NX,NZ,NCX_SHOT, NCZ_SHOT,kboud,kp,iw 
       INTEGER   IIJJ,L,is, ng_left, kkjj
       DIMENSION C(5)
       DIMENSION U( -4:NX+5 , -4:NZ+5 )
       DIMENSION U1( -4:NX+5 , -4:NZ+5 )
       DIMENSION U2( -4:NX+5 , -4:NZ+5 )
       DIMENSION V(NX, NZ)
       DIMENSION SEISM(NTRACE, LT)
       DIMENSION IRZ(NX)
       DIMENSION whole_stat_elev(nmax)

       mm=10
       mm=mm/2
       if(mm.eq.2) then
	   sg=0.
       else
	   sg=1.
       endif
c===================================================================c
c                  Differencial Parameter Specified                 c
c===================================================================c
       if(mm.eq.1) then
	    c(1)=2
	    c(2)=0.
            c(3)=0.
	    c(4)=0.
	    c(5)=0.
       elseif(mm.eq.2) then
	    c(1)=8./3.
            c(2)=-1./6.
	    c(3)=0.
	    c(4)=0.
	    c(5)=0.
       elseif(mm.eq.3) then
	    c(1)=3.
            c(2)=-3./10
	    c(3)=1./45.
	    c(4)=0.
	    c(5)=0.
       elseif(mm.eq.4) then
	    c(1)=16./5.
            c(2)=-2./5.
	    c(3)=16./315.
	    c(4)=-1./280.
	    c(5)=0.
       elseif(mm.eq.5) then
            c(1)=3.3333333
            c(2)=-0.4761905
            c(3)=0.07936508
            c(4)=-0.009920635
            c(5)=0.0006349206
c	    c0=-2.927222164
c	    c(1)=1.66666665
c            c(2)=-0.23809525
c	    c(3)=0.03968254
c	    c(4)=-0.004960318
c	    c(5)=0.0003174603
       endif
            c0=0.
            do i=1,mm
               c0=c0-c(i)
            enddo
            do i=1,mm
               c(i)=0.5*c(i)
            enddo

        dtx=dt/dx
	dtz=dt/dz
	dtxz=dtx*dtz

        dr1=dtx*dtx/2.
	dr2=dtz*dtz/2.

        dtx4=dtx*dtx*dtx*dtx
	dtz4=dtz*dtz*dtz*dtz
	dtxz4=dtx*dtx*dtz*dtz
	
c===================Time Extrapolation Now Begining===================c

        do 30 l=1,LT
	  
         if(mod(l,100).eq.0) then
         write(*,*)'shot=',is,'     TIME=',l
         endif

c------------Grid Points Circling Begining in a Time Slice----------c	  

	  do 40 i=1,nx
	  do 55 k=1,nz
		  
	    vv2=v(i,k)*v(i,k)
            drd1=dr1*vv2
            drd2=dr2*vv2
c-------------------------------------------------------------------c
c                    Calculating Wavelet                            c
c-------------------------------------------------------------------c
	    tt=(l-1)*dt
	    tt=tt-t0
	    call gauss(f0,tt,fx)
c-------------------------------------------------------------------c
c                     Source Specified                              c
c-------------------------------------------------------------------c
	    if(kp.eq.1) then
	        if(i.eq.NCX_SHOT.and.k.eq.NCZ_SHOT) then
	          qx=1.
	        else
	          qx=0.
	        endif
	    else
	       if(i.eq.NCX_SHOT-1.and.(k.eq.NCZ_SHOT.or.
     +          k.eq.NCZ_SHOT+1)) then
	         qx=1.
	       else if(i.eq.NCX_SHOT.and.(k.eq.NCZ_SHOT-1.or.
     +	           k.eq.NCZ_SHOT+2)) then
	         qx=1.
	       else if(i.eq.NCX_SHOT+1.and.(k.eq.NCZ_SHOT-1.or.
     +               k.eq.NCZ_SHOT+2)) then
	         qx=1.
	       else if(i.eq.NCX_SHOT+2.and.(k.eq.NCZ_SHOT.or.
     +               k.eq.NCZ_SHOT+1)) then
	         qx=1.
	       else
	         qx=0.
	       endif
            endif
      
      if(kboud.eq.1) then

c-------------------------------------------------------------------c
c         Absorbing Boundaries Starting                             c
c-------------------------------------------------------------------c
      
        if(k.gt.nz-5.and.i.gt.5.and.i.lt.nx-5) then
	  p1=drd1*(u1(i+1,k)-2*u1(i,k)+u1(i-1,k))
	  p2=-dtz*v(i,k)*(u1(i,k)-u2(i,k)-u1(i,k-1)+u2(i,k-1))
	  u(i,k)=2*u1(i,k)-u2(i,k)+p1+p2
	else if(k.lt.6.and.i.gt.5.and.i.lt.nx-5) then
	  p1=drd1*(u1(i+1,k)-2*u1(i,k)+u1(i-1,k))
	  p2=dtz*v(i,k)*(u1(i,k+1)-u2(i,k+1)-u1(i,k)+u2(i,k))
	  u(i,k)=2*u1(i,k)-u2(i,k)+p1+p2
	else if(i.lt.6.and.k.gt.5.and.k.lt.nz-5) then
	  p1=drd2*(u1(i,k+1)-2*u1(i,k)+u1(i,k-1))
	  p2=dtx*v(i,k)*(u1(i+1,k)-u2(i+1,k)-u1(i,k)+u2(i,k))
	  u(i,k)=2*u1(i,k)-u2(i,k)+p1+p2
	else if(i.gt.nx-5.and.k.gt.5.and.k.lt.nz-5) then
	  p1=drd2*(u1(i,k+1)-2*u1(i,k)+u1(i,k-1))
	  p2=-dtx*v(i,k)*(u1(i,k)-u2(i,k)-u1(i-1,k)+u2(i-1,k))
	  u(i,k)=2*u1(i,k)-u2(i,k)+p1+p2
	else if(i.lt.6.and.k.lt.6) then
	  p1=dtx*v(i,k)*(u1(i+1,k)-u2(i+1,k)-u1(i,k)+u2(i,k))
	  p2=dtz*v(i,k)*(u1(i,k+1)-u2(i,k+1)-u1(i,k)+u2(i,k))
	  p3=-0.75*vv2*dtxz
	  p3=p3*(u1(i+1,k+1)-u1(i+1,k)-u1(i,k+1)+u1(i,k))
	  u(i,k)=2*u1(i,k)-u2(i,k)+p1+p2+p3
	else if(i.gt.nx-5.and.k.lt.6) then
	  p1=-dtx*v(i,k)*(u1(i,k)-u2(i,k)-u1(i-1,k)+u2(i-1,k))
	  p2=dtz*v(i,k)*(u1(i,k+1)-u2(i,k+1)-u1(i,k)+u2(i,k))
	  p3=0.75*vv2*dtxz
	  p3=p3*(u1(i,k+1)-u1(i,k)-u1(i-1,k+1)+u1(i-1,k))
	  u(i,k)=2*u1(i,k)-u2(i,k)+p1+p2+p3
	else if(i.gt.nx-5.and.k.gt.nz-5) then
	  p1=-dtx*v(i,k)*(u1(i,k)-u2(i,k)-u1(i-1,k)+u2(i-1,k))
	  p2=-dtz*v(i,k)*(u1(i,k)-u2(i,k)-u1(i,k-1)+u2(i,k-1))
	  p3=-0.75*vv2*dtxz
	  p3=p3*(u1(i,k)-u1(i-1,k)-u1(i,k-1)+u1(i-1,k-1))
	  u(i,k)=2*u1(i,k)-u2(i,k)+p1+p2+p3
        else if(i.lt.6.and.k.gt.nz-5) then
	  p1=dtx*v(i,k)*(u1(i+1,k)-u2(i+1,k)-u1(i,k)+u2(i,k))
	  p2=-dtz*v(i,k)*(u1(i,k)-u2(i,k)-u1(i,k-1)+u2(i,k-1))
	  p3=0.75*vv2*dtxz
	  p3=p3*(u1(i+1,k)-u1(i,k)-u1(i+1,k-1)+u1(i,k-1))
	  u(i,k)=2*u1(i,k)-u2(i,k)+p1+p2+p3
	else
c-------------------------------------------------------------------c
            temp1=0.
	    temp2=0.
		  do 100 kk=1,mm
		    temp1=temp1+c(kk)*(u1(i+kk,k)+u1(i-kk,k))
		    temp2=temp2+c(kk)*(u1(i,k+kk)+u1(i,k-kk))
100		  continue
		  
		  temp1=(temp1+c0*u1(i,k))*vv2*dtx*dtx
		  temp2=(temp2+c0*u1(i,k))*vv2*dtz*dtz
		  
	      p1=dtx4*(u1(i+2,k)-4.*u1(i+1,k)+6.*u1(i,k)
     +         -4.*u1(i-1,k)+u1(i-2,k))
	      
	      p2=dtz4*(u1(i,k+2)-4.*u1(i,k+1)+6.*u1(i,k)
     +	         -4.*u1(i,k-1)+u1(i,k-2))
	      
	      p3=u1(i+1,k+1)-2.*u1(i+1,k)+u1(i+1,k-1)-2.*u1(i,k+1)
     +	         +4.*u1(i,k)-2.*u1(i,k-1)+u1(i-1,k+1)-2.*u1(i-1,k)
     +         +u1(i-1,k-1)
	      p3=p3*2*dtxz4

            p=p1+p2+p3
	      
            q1=vv2*vv2*p/12.
	      q2=temp1+temp2
              
            u(i,k)=2.*u1(i,k)-u2(i,k)+sg*q1+q2+qx*fx
         endif

c-------------------------------------------------------------------c
c                Absorbing Boundaries End!                          c
c-------------------------------------------------------------------c
      
        else
		 
c==================================================================c 
c           Characteristic analysis Absorbing Boundaries Begin!    c
c==================================================================c

                  temp1=0.
		  temp2=0.
		  do 101 kk=1,mm
		    temp1=temp1+c(kk)*(u1(i+kk,k)+u1(i-kk,k))
		    temp2=temp2+c(kk)*(u1(i,k+kk)+u1(i,k-kk))
101		  continue
		  
		  temp1=(temp1+c0*u1(i,k))*vv2*dtx*dtx
		  temp2=(temp2+c0*u1(i,k))*vv2*dtz*dtz
		  
	      p1=dtx4*(u1(i+2,k)-4.*u1(i+1,k)+6.*u1(i,k)
     +         -4.*u1(i-1,k)+u1(i-2,k))
	      
	      p2=dtz4*(u1(i,k+2)-4.*u1(i,k+1)+6.*u1(i,k)
     +	         -4.*u1(i,k-1)+u1(i,k-2))
	      
	      p3=u1(i+1,k+1)-2.*u1(i+1,k)+u1(i+1,k-1)-2.*u1(i,k+1)
     +	         +4.*u1(i,k)-2.*u1(i,k-1)+u1(i-1,k+1)-2.*u1(i-1,k)
     +         +u1(i-1,k-1)
	      p3=p3*2*dtxz4

              p=p1+p2+p3
	      
              q1=vv2*vv2*p/12.
	      q2=temp1+temp2
              
              u(i,k)=2.*u1(i,k)-u2(i,k)+sg*q1+q2+qx*fx

	if(k.ge.nz-20.and.i.gt.21.and.i.lt.nx-20) then
	  pq=0.01*(k-nz+20)
	  p3=u(i,k)-u1(i,k)+0.5*dtz*v(i,k)*(u1(i,k+1)-u1(i,k-1))
	  p3=pq*p3
	  u(i,k)=u(i,k)-p3
	elseif(k.le.21.and.i.gt.21.and.i.lt.nx-20) then
	  pq=0.01*(21-k)
	  p4=u(i,k)-u1(i,k)-0.5*dtz*v(i,k)*(u1(i,k+1)-u1(i,k-1))
	  p4=pq*p4
	  u(i,k)=u(i,k)-p4
	elseif(i.le.21.and.k.gt.21.and.k.lt.nz-20) then
	  pq=0.01*(21-i)
	  p4=u(i,k)-u1(i,k)-0.5*dtx*v(i,k)*(u1(i+1,k)-u1(i-1,k))
	  p4=pq*p4
	  u(i,k)=u(i,k)-p4
	elseif(i.ge.nx-20.and.k.gt.21.and.k.lt.nz-20) then
	  pq=0.01*(i-nx+20)
	  p3=u(i,k)-u1(i,k)+0.5*dtx*v(i,k)*(u1(i+1,k)-u1(i-1,k))
	  p3=pq*p3
	  u(i,k)=u(i,k)-p3
	elseif(i.le.21.and.k.le.21) then
	  pq=0.01*(21-i)
	  p4=u(i,k)-u1(i,k)-0.5*dtx*v(i,k)*(u1(i+1,k)-u1(i-1,k))
	  p4=pq*p4
	  pq=0.01*(21-k)
	  p3=u(i,k)-u1(i,k)-0.5*dtz*v(i,k)*(u1(i,k+1)-u1(i,k-1))
	  p3=pq*p3
	  u(i,k)=u(i,k)-p3-p4
	elseif(i.ge.nx-20.and.k.le.21) then
	  pq=0.01*(21-k)
	  p4=u(i,k)-u1(i,k)-0.5*dtz*v(i,k)*(u1(i,k+1)-u1(i,k-1))
	  p4=pq*p4
	  pq=0.01*(i-nx+20)
	  p3=u(i,k)-u1(i,k)+0.5*dtx*v(i,k)*(u1(i+1,k)-u1(i-1,k))
	  p3=pq*p3
	  u(i,k)=u(i,k)-p3-p4
	elseif(i.ge.nx-20.and.k.ge.nz-20) then
	  pq=0.01*(k-nz+20)
	  p4=u(i,k)-u1(i,k)+0.5*dtz*v(i,k)*(u1(i,k+1)-u1(i,k-1))
	  p4=pq*p4
	  pq=0.01*(i-nx+20)
	  p3=u(i,k)-u1(i,k)+0.5*dtx*v(i,k)*(u1(i+1,k)-u1(i-1,k))
	  p3=pq*p3
	  u(i,k)=u(i,k)-p3-p4
	elseif(i.le.21.and.k.ge.nz-20) then
	  pq=0.01*(k-nz+20)
	  p4=u(i,k)-u1(i,k)+0.5*dtz*v(i,k)*(u1(i,k+1)-u1(i,k-1))
	  p4=pq*p4
	  pq=0.01*(21-i)
	  p3=u(i,k)-u1(i,k)-0.5*dtx*v(i,k)*(u1(i+1,k)-u1(i-1,k))
	  p3=pq*p3
	  u(i,k)=u(i,k)-p3-p4
	endif
	
c===================================================================c
c        Characteristic analysis Absorbing Boundaries Finished!     c
c===================================================================c
       endif
     
55	  continue
40	  continue

c-------------------------------------------------------------------c
 
        do 110 i=1,nx
	     do 110 k=1,nz
	        u2(i,k)=u1(i,k)
              u1(i,k)=u(i,k)
110     continue
c-------------------------------------------------------------------c

	 do i=1,ntrace 
	    iijj=1+(i-1)*iw

            if((iijj+ng_left-1).lt.1) then
            whole_stat_elev(iijj+ng_left-1)=whole_stat_elev(1)
            elseif((iijj+ng_left-1).gt.nmax) then
            whole_stat_elev(iijj+ng_left-1)=whole_stat_elev(nmax)
            endif

c           write(*,*)i,iijj+ng_left-1,kkjj-20
            kkjj=int(whole_stat_elev(iijj+ng_left-1))+20
            seism(i,l)=u(iijj,kkjj)
         enddo

c       if(l.eq.1) pause

30	continue

       return   
	end

*=============================================================*

	SUBROUTINE RICK(F,TT,FX)

	real f,tt,fx,sp,PIE
	PIE=3.1415926
	sp=PIE*f*tt
	fx=1000.*exp(-sp*sp)*(1.-2.*sp*sp)
	end

	subroutine gauss(f,tt,fx)
	real f,tt,fx,sp,PIE
	PIE=3.1415926
	sp=PIE*f*tt
	fx=-2000.*pie*pie*f*f*tt*exp(-sp*sp)*(3.-2.*sp*sp)

	ENd


*=============================================================

      SUBROUTINE WRITE_DISK(IS,NTRACE,LT,lt_new,SEISM,T0,DT,IJ,
     +   NA,rdt,dx_trace_new,SEISM_OUT,dx_trace,ipc)

        INTEGER IS, NTRACE, LT, lt_new, KREC, IJ, NA
        INTEGER IT1, rdt, ipc
        REAL T0, DT,dx_trace_new,dx_trace
        INTEGER RE
        DIMENSION SEISM(NTRACE,LT)
        DIMENSION SEISM_OUT(NTRACE,LT) 
        CHARACTER a,b,c,d

c       write(*,*) ntrace,lt
c       open(30,file='ttt',access='direct',recl=1*lt)
c       do ix=1,ntrace
c       write(30,rec=ix) (seism(ix,it),it=1,lt)
c       end do
c       close(30)
c       write(*,*) 'ttt'

        RE=int(dx_trace_new/dx_trace)

        write(*,*) 'new resample time interval'
        write(*,*) rdt

        NTRACE_NEW=0
        DO IX=1+NA, NTRACE-NA,RE
          NTRACE_NEW=NTRACE_NEW+1
          IT1=0
          DO IT=IJ+1, LT, rdt
            IT1=IT1+1
            SEISM_OUT(NTRACE_NEW, IT1)=SEISM(IX, IT)
          ENDDO
        ENDDO
c       write(*,*) 'traces after resample=',NTRACE_NEW

        a=char((is-is/10000*10000)/1000+48)
        b=char((is-is/1000*1000)/100+48)
        c=char((is-is/100*100)/10+48)
        d=char((is-is/10*10)+48)

        OPEN(10,FILE=
     +'shot_cnooc_fault_'//a//b//c//d//'.dat', ACCESS='DIRECT',
     +   status='replace',RECL=ipc*lt_new)        
 
        DO IX=1,NTRACE_NEW 
           WRITE(10, REC=IX) (SEISM_OUT(IX, IT), IT=1, lt_new)
        END DO

        CLOSE(10)
        RETURN
        END
        
c---------------------------------------------------------
        SUBROUTINE READPAR(VEL_FNAME, FN3,
     +        NSHOT, NTRACE, NMAX, NMAZ,
     +        DX, DZ, XS_START, DS_FR_FS,
     +        DX_SHOT, DX_TRACE, DT, LT,
     +        PEAK, ipc,dt_new)

        INTEGER NTRACE,NMAX,NMAZ,NSHOT,NZ,LT
        REAL    DX,DZ,XS_START,DX_SHOT,DX_TRACE,dx_trace_new,SP_Z
        REAL    XG_START
        REAL    dt_new 

        character*256 VEL_FNAME,FN3

        OPEN(15,FILE='final_2d_ac_modeling_no_mpi.par')

        read(15,'(a)') vel_fname
        write(*,*) ' velocity file name'
        write(*,'(a)') vel_fname

        read(15,'(a)') FN3 
        write(*,*) ' Station vs Elevation file name'
        write(*,'(a)')  FN3

        READ (15,*) nshot,ntrace
        write(*,*)  'total shot number in this line'
        write(*,*)  'trace number in each shot'
        write(*,*)  nshot,ntrace

        READ(15, *) nmax,nmaz
        write(*,*) 'maximum number sample number'
        write(*,*) 'in X-,Z- direction'
        write(*,*) nmax,nmaz

        READ(15,*) DX, DZ
        WRITE(*,*) 'X-,Z- DIRECTION SAMPLE RATE OF VELOCITY MODEL(m)'
        WRITE(*,*) DX, DZ

        READ(15, *) XS_START, DS_FR_FS
        write(*,*) 'FIRST SHOT&RECEIVER LOCATION(m)'
        write(*,*)  XS_START, DS_FR_FS

        READ(15,*) DX_SHOT, DX_TRACE
        WRITE(*,*)'SHOT&RECEIVER interval(m)'
        WRITE(*,*) DX_SHOT, DX_TRACE

        READ(15,*) DT, LT
        WRITE(*,*)'expectful sample interval and sample points'
        WRITE(*,*) DT, LT

        READ(15,*) PEAK
        WRITE(*,*)'PEAK Elevation in this line '
        WRITE(*,*) PEAK

        READ(15,*) ipc 
        WRITE(*,*)'1 or 4 '
        WRITE(*,*) ipc

        READ(15,*) dt_new
        WRITE(*,*)'NEW TIME SAMPLE INTERVAL (ms)'
        WRITE(*,*) dt_new

        CLOSE(15)

        RETURN
        END
c---------------------------------------------------------

*=============================================================*
*                                                             *
*            End of the file                                  *
*                                                             *
*=============================================================*

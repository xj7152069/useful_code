
!==================================================================
	MODULE INTERP_GLOBAL
		INTEGER,PARAMETER::Ntable=101
		INTEGER,PARAMETER::Ltable=4
		INTEGER,PARAMETER::length=1
		INTEGER ixxx,izzz
		INTEGER lx,lz
		INTEGER mx,mz
		INTEGER nvxm,nvzm
		REAL	ax,az
		REAL	tbl_v(4,Ntable),tbl_vx(4,Ntable)
	END MODULE INTERP_GLOBAL
!==================================================================

    PROGRAM RAY_TRACE_2D				
		USE INTERP_GLOBAL

		INTEGER NX, NZ, NS_X, NS_Z, NVX, NVZ, NVXS
		INTEGER NXS_LEFT, NXS_RIGHT
		INTEGER ISHOT,myid
        REAL    DX, DZ, DVX, DVZ
		REAL	SX_COORD, VX_START
		REAL	DSTEP,DXS,DZS
		REAL	DS_X, DS_Z
		REAL	SZ_COORD!,maxpath
        real    SX_LEFT,SX_RIGHT,DEPTH,nvzs
		integer ntraces,trace_start,trace_locate

        REAL,ALLOCATABLE::TIME(:,:)
        REAL,ALLOCATABLE::TIME2(:,:)
        REAL,ALLOCATABLE::pathmat(:,:)

!*	START COORDINATE IN X/Z DIMENSION OF THE VELOCITY MODEL	*!
!*	0.0,0.0
!=====================================Par=====================================
!model sample begin with (1,1); source should begin with (2,2)
        NS_X=40
        NS_Z=2
        NR_X=40
        NR_Z=60

!model size (depth = (nvz-1)*dvz)
        nvx=80
        nvz=60
        DEPTH=600.0

!sample gep
        DVX=10.0
		DVZ=10.0
        dxs=10.0
        dzs=10.0

!source distance to the left and right boundary
		SX_LEFT=(NS_X-1)*DVX
		SX_RIGHT=(NVX-1)*DVX-SX_LEFT

        VX_START=0.0
        VZ_START=0.0

        dstep=10.0

!=============================================================================

        sx_coord=(NS_X-1)*DVX
		sz_coord=(NS_Z-1)*DVZ
        rx_coord=(NR_X-1)*DVX
		rz_coord=(NR_Z-1)*DVZ

		DX=DVX
		DZ=DVZ

		nvxs=ceiling(nvx/(dxs/dvx))
		nvzs=ceiling(nvz/(dzs/dvz))

		NXS_LEFT=SX_LEFT/DVX+0.5
		NXS_RIGHT=SX_RIGHT/DVX+0.5

		NX=NXS_LEFT+NXS_RIGHT+1
		NZ=DEPTH/DZ
		
		NV=(NVX-1)*(NVZ-1)
		NVS=(NVXS-1)*(NVZS-1)

		XMIN=0.0
		ZMIN=0.0

!=============================================================================
		ALLOCATE(TIME(NX,NZ))
        ALLOCATE(TIME2(NZ,NX))
        ALLOCATE(pathmat(NZ,NX))


        Open(12,File = "time.bin" , access="stream" , form = "unformatted" )
        Read( 12 ) time2
        close(12)

        DO IX=1, NX
			DO IZ=1, NZ
			    TIME(ix,iz)=TIME2(iz,ix)
			END DO
		END DO

!		CUR_RAY     =trim(adjustl(FN4))//'ray'//trim(adjustl(CUR_IS))//'.txt'
!		CUR_PATH    =trim(adjustl(FN5))//'path'//trim(adjustl(CUR_IS))//'.dat'
!		CUR_TIME    =trim(adjustl(FN6))//'traveltime'//trim(adjustl(CUR_IS))//'.dat'
        CALL PreInterpolators_Bspline()

!=============================================================================
!		DO i=1,ntraces
!			READ(13,REC=i+trace_start-1) geo
		
            Open(12,File = "time.ray.bin" , access="stream" , form = "unformatted",status='replace')
		    OPEN(14,FILE='test.ray',STATUS='REPLACE')
		    OPEN(18,FILE='test.pathmat',access="stream", STATUS='REPLACE')
            OPEN(19,FILE='test.path',STATUS='REPLACE')

			NR_X=NINT((RX_COORD-VX_START)/DVX)+1

!			IF(NR_X.LT.1) CYCLE
!			IF(NR_X.GT.NVX) CYCLE

!			ITR=trace_locate+i-1
            ITR=0
			pathmat(:,:)=0.0
!			if(ITR.eq.1) then

            IF(NR_X>=1 .and. NR_X<=NVX) then
			    CALL TRACING(NX, NZ, DX, DZ, TIME, XMIN, ZMIN, DSTEP,&
						SX_COORD,SZ_COORD,RX_COORD,RZ_COORD,VX_START,&
						NVXS, DXS, DZS, ITR, time2, pathmat)
            end if
            write( 12 ) time2
            write( 18 ) pathmat
        
!		END DO
!=============================================================================
        CLOSE(14)
		CLOSE(18)
		CLOSE(19)
        close(12)

	END

!===========================================================================
SUBROUTINE TRACING(NX, NZ, DX, DZ, TIME, XMIN, ZMIN, DSTEP,&
				SX_COORD,SZ_COORD,RX_COORD,RZ_COORD,VX_START,&
				NVXS, DXS, DZS, ITR, time2, pathmat)
USE INTERP_GLOBAL
IMPLICIT NONE

INTEGER NX, NZ, NVXS
INTEGER ITR
REAL	DX, DZ
REAL	XMIN, ZMIN
REAL	SX_COORD,SZ_COORD,RX_COORD,RZ_COORD
REAL	VX_START
REAL	TMP
REAL	DSTEP
REAL	DXS,DZS
REAL,DIMENSION(NX,NZ)::TIME
REAL,DIMENSION(NZ,NX)::TIME2
REAL,DIMENSION(NZ,NX)::pathmat

REAL	X, Z, T, DTDX, DTDZ	!CURRENT POINT PARAMETER
REAL	X1,Z1,X2,Z2,X3,Z3,X4,Z4
INTEGER IVX1,IVZ1,IVX2,IVZ2,IVX,IVZ,ISX,ISZ,IV,I,ik,k
REAL	PATH,PATH1,PATH2,TMP1,TMP2,TMP3,TMP4
REAL 	m,l,a,b

X=RX_COORD
Z=RZ_COORD

PATH=0.0
PATH1=0.0
PATH2=0.0
ik=1

CALL VELINTERPED(X, Z, T, XMIN, ZMIN, DX, DZ,&
				DTDX, DTDZ, NX, NZ, TIME)

TMP=T

m=DTDX
l=DTDZ
                
IF(.NOT.(abs(rx_coord-sx_coord)+abs(rz_coord-sz_coord))>0.0)THEN !!!!!!!!!!!!!!!!!!!
	TMP=-1
	GOTO 1201
ELSE
	WRITE(14,*) 1, X, Z!, 'geophone'                                              !raycoordinate

	X1=X
	Z1=Z
	IVX1=INT((X1-VX_START)/DXS)+1
	IVZ1=INT(Z1/DZS)+1
!	PRINT*,'1716',X1,Z1,IVX1,IVZ1

	DO WHILE(.TRUE.)
		X=X-DSTEP*DTDX/SQRT(DTDX*DTDX+DTDZ*DTDZ)
		Z=Z-DSTEP*DTDZ/SQRT(DTDX*DTDX+DTDZ*DTDZ)
		TMP3=SQRT((X-SX_COORD)*(X-SX_COORD)+&
        		(Z-SZ_COORD)*(Z-SZ_COORD))
									
		IF(Z.LT.0.0) THEN
			TMP=-1
			GOTO 1201
		ELSE
			WRITE(14,*) 0, X, Z

            mx=INT(X/DX)+1
            mz=INT(Z/DZ)+1
            IF(mx>=1 .and. mx<=(nx) .and. mz>=1 .and. mz<=(nz) .and. time2(mz,mx)<50)THEN
                time2(mz,mx)=time2(mz,mx)*10                                           
            end if

!raycoordinate

			X2=X
			Z2=Z
			IVX2=INT((X2-VX_START)/DXS)+1
			IVZ2=INT(Z2/DZS)+1
!			PRINT*,'1870',X2,Z2,IVX2,IVZ2
			IF((IVX1.EQ.IVX2).AND.(IVZ1.EQ.IVZ2)) THEN
				PATH=PATH+DSTEP
			ELSE IF(IVX1.EQ.IVX2) THEN
				PATH1=ABS((Z1-(MAX(IVZ1,IVZ2)-1)*DZS)*&
					SQRT(DTDX*DTDX+DTDZ*DTDZ)/DTDZ)
				PATH=PATH+PATH1
				IVX=IVX1
				IVZ=IVZ1
				IV=(IVZ-1)*(NVXS-1)+IVX

                pathmat(IVZ,IVX)=pathmat(IVZ,IVX)+PATH
                write(19,*) IVZ,IVX,PATH
				ik=ik+1
!				WRITE(17,*) ITR, IV, PATH
!				PRINT*,'1881',ITR,IV,PATH
				PATH2=DSTEP-PATH1
				PATH=PATH2
			ELSE IF(IVZ1.EQ.IVZ2) THEN
				PATH1=ABS((X1-(MAX(IVX1,IVX2)-1)*DXS-VX_START)*&
					SQRT(DTDX*DTDX+DTDZ*DTDZ)/DTDX)
				PATH=PATH+PATH1
				IVX=IVX1
				IVZ=IVZ1
				IV=(IVZ-1)*(NVXS-1)+IVX	

                pathmat(IVZ,IVX)=pathmat(IVZ,IVX)+PATH
!                write(*,*) IVZ,IVX,PATH
                write(19,*) IVZ,IVX,PATH
				ik=ik+1
!				WRITE(17,*) ITR, IV, PATH
!				PRINT*,'1892',ITR,IV,PATH
				PATH2=DSTEP-PATH1
				PATH=PATH2
			ELSE
				TMP1=ABS((Z1-(MAX(IVZ1,IVZ2)-1)*DZS)*&
						SQRT(DTDX*DTDX+DTDZ*DTDZ)/DTDZ)
				TMP2=ABS((X1-(MAX(IVX1,IVX2)-1)*DXS-VX_START)*&
						SQRT(DTDX*DTDX+DTDZ*DTDZ)/DTDX)
				PATH1=MIN(TMP1,TMP2)
				PATH=PATH+PATH1

				IVX=IVX1
				IVZ=IVZ1
				IV=(IVZ-1)*(NVXS-1)+IVX	

                pathmat(IVZ,IVX)=pathmat(IVZ,IVX)+PATH
!                write(*,*) IVZ,IVX,PATH
                write(19,*) IVZ,IVX,PATH
				ik=ik+1
!				WRITE(17,*) ITR, IV, PATH
!				PRINT*,'1907',ITR,IV,PATH

				PATH=MAX(TMP1,TMP2)-MIN(TMP1,TMP2)

				X3=X-TMP1*DTDX/SQRT(DTDX*DTDX+DTDZ*DTDZ)
				Z3=Z-TMP1*DTDZ/SQRT(DTDX*DTDX+DTDZ*DTDZ)
				X4=X-TMP2*DTDX/SQRT(DTDX*DTDX+DTDZ*DTDZ)
				Z4=Z-TMP2*DTDZ/SQRT(DTDX*DTDX+DTDZ*DTDZ)
        
				IVX=INT((0.5*(X3+X4)-VX_START)/DXS)+1
				IVZ=INT(0.5*(Z3+Z4)/DZS)+1
				IV=(IVZ-1)*(NVXS-1)+IVX	

                pathmat(IVZ,IVX)=pathmat(IVZ,IVX)+PATH
!                write(*,*) IVZ,IVX,PATH
                write(19,*) IVZ,IVX,PATH
				ik=ik+1
!				WRITE(17,*) ITR, IV, PATH
!				PRINT*,'1920',ITR,IV,PATH
          
				PATH2=DSTEP-MAX(TMP1,TMP2)
				PATH=PATH2

			END IF
		END IF                       

		CALL VELINTERPED(X, Z, T, XMIN, ZMIN, DX, DZ,&
							DTDX, DTDZ, NX, NZ, TIME)
   
!		PRINT*,'X,Z',X,Z
!		PRINT*,'DTDX',DTDX,'DTDZ',DTDZ
!       print*,'1799',X,Z,T

		X1=X
		Z1=Z
		IVX1=INT((X1-VX_START)/DXS)+1
		IVZ1=INT(Z1/DZS)+1

		IF(TMP3.LT.(1*DSTEP)) THEN
			PATH=PATH+TMP3
			WRITE(14,*) 2,sx_coord,sz_coord!,'source'
            mx=INT(sx_coord/DX)+1
            mz=INT(sz_coord/DZ)+1
            IF(mx>=1 .and. mx<=(nx) .and. mz>=1 .and. mz<=(nz) .and. time2(mz,mx)<50)THEN
                time2(mz,mx)=time2(mz,mx)*10                                           
            end if
			EXIT
		END IF
							   

		IF(.NOT.(abs(x1-sx_coord)+abs(z1-sz_coord))>0.0)THEN
			TMP=-1
			GOTO 1201
		END IF
		
		b=(m*dtdx+l*dtdz)/((sqrt(m*m+l*l))*(sqrt(dtdx*dtdx+dtdz*dtdz)))-0.000001
		a=ACOS(b) !Calucate Angle(radian)

		IF(.NOT.(0.0<=a.AND.a<=1.256)) THEN
			TMP=-1
			GOTO 1201
		END IF

		m=dtdx
		l=dtdz

	END DO
	
	IVX=IVX1
	IVZ=IVZ1
	IV=(IVZ-1)*(NVXS-1)+IVX	

    pathmat(IVZ,IVX)=pathmat(IVZ,IVX)+PATH
!    write(*,*) IVZ,IVX,PATH
    write(19,*) IVZ,IVX,PATH
	ik=ik+1
	
!	WRITE(17,*) ITR, IV, PATH
END IF
	
1201	CONTINUE
	
WRITE(14,*) 3, TMP, TMP*1000!, 'ray_time(s|ms)'
!print* , 'itr,tmp', itr,tmp	

END

!===========================================================================

	SUBROUTINE VELINTERPED(x,z,v,Xmin,Zmin,dvx,dvz,&
					dvdx,dvdz,nvx,nvz,Vel)
		USE INTERP_GLOBAL
		INTEGER nvx,nvz
		REAL	a,b,c,d
		REAL    x,z,v,Xmin,Zmin,dvx,dvz
		REAL    dvdx,dvdz
		REAL    Vel(nvx, nvz)

!********* determine offsets into v and interpolation coefficients *
		ax=(x-Xmin)/dvx+1
		ixxx=ax
		lx=(ax-ixxx)*(Ntable-1)+1.5
		mx=ixxx-Ltable/2+1
		az=(z-Zmin)/dvz+1
		izzz=az
		lz=(az-izzz)*(Ntable-1)+1.5
		mz=izzz-Ltable/2+1		
		nvxm=nvx-Ltable+1
		nvzm=nvz-Ltable+1
		
!********* if totally within input array, use fast method **********      
		IF(mx>=1 .and. mx<=(nvxm) .and. mz>=1 .and. mz<=(nvzm))THEN
			v = Vel(mx  ,mz  )*tbl_v(1,lz)*tbl_v(1,lx)+&
				Vel(mx+1,mz  )*tbl_v(2,lz)*tbl_v(1,lx)+&
				Vel(mx+2,mz  )*tbl_v(3,lz)*tbl_v(1,lx)+&
				Vel(mx+3,mz  )*tbl_v(4,lz)*tbl_v(1,lx)+&
				Vel(mx  ,mz+1)*tbl_v(1,lz)*tbl_v(2,lx)+&
				Vel(mx+1,mz+1)*tbl_v(2,lz)*tbl_v(2,lx)+&
				Vel(mx+2,mz+1)*tbl_v(3,lz)*tbl_v(2,lx)+&
				Vel(mx+3,mz+1)*tbl_v(4,lz)*tbl_v(2,lx)+&
				Vel(mx  ,mz+2)*tbl_v(1,lz)*tbl_v(3,lx)+&
				Vel(mx+1,mz+2)*tbl_v(2,lz)*tbl_v(3,lx)+&
				Vel(mx+2,mz+2)*tbl_v(3,lz)*tbl_v(3,lx)+&
				Vel(mx+3,mz+2)*tbl_v(4,lz)*tbl_v(3,lx)+&
				Vel(mx  ,mz+3)*tbl_v(1,lz)*tbl_v(4,lx)+&
				Vel(mx+1,mz+3)*tbl_v(2,lz)*tbl_v(4,lx)+&
				Vel(mx+2,mz+3)*tbl_v(3,lz)*tbl_v(4,lx)+&
				Vel(mx+3,mz+3)*tbl_v(4,lz)*tbl_v(4,lx)
				
			dvdx=(Vel(mx ,mz )*tbl_vx(1,lx)*tbl_v(1,lz)+&
				Vel(mx+1,mz  )*tbl_vx(2,lx)*tbl_v(1,lz)+&
				Vel(mx+2,mz  )*tbl_vx(3,lx)*tbl_v(1,lz)+&
				Vel(mx+3,mz  )*tbl_vx(4,lx)*tbl_v(1,lz)+&
				Vel(mx  ,mz+1)*tbl_vx(1,lx)*tbl_v(2,lz)+&
				Vel(mx+1,mz+1)*tbl_vx(2,lx)*tbl_v(2,lz)+&
				Vel(mx+2,mz+1)*tbl_vx(3,lx)*tbl_v(2,lz)+&
				Vel(mx+3,mz+1)*tbl_vx(4,lx)*tbl_v(2,lz)+&
				Vel(mx  ,mz+2)*tbl_vx(1,lx)*tbl_v(3,lz)+&
				Vel(mx+1,mz+2)*tbl_vx(2,lx)*tbl_v(3,lz)+&
				Vel(mx+2,mz+2)*tbl_vx(3,lx)*tbl_v(3,lz)+&
				Vel(mx+3,mz+2)*tbl_vx(4,lx)*tbl_v(3,lz)+&
				Vel(mx  ,mz+3)*tbl_vx(1,lx)*tbl_v(4,lz)+&
				Vel(mx+1,mz+3)*tbl_vx(2,lx)*tbl_v(4,lz)+&
				Vel(mx+2,mz+3)*tbl_vx(3,lx)*tbl_v(4,lz)+&
				Vel(mx+3,mz+3)*tbl_vx(4,lx)*tbl_v(4,lz))/dvx

			dvdz=(Vel(mx ,mz )*tbl_v(1,lx)*tbl_vx(1,lz)+&
				Vel(mx+1,mz  )*tbl_v(2,lx)*tbl_vx(1,lz)+&
				Vel(mx+2,mz  )*tbl_v(3,lx)*tbl_vx(1,lz)+&
				Vel(mx+3,mz  )*tbl_v(4,lx)*tbl_vx(1,lz)+&
				Vel(mx  ,mz+1)*tbl_v(1,lx)*tbl_vx(2,lz)+&
				Vel(mx+1,mz+1)*tbl_v(2,lx)*tbl_vx(2,lz)+&
				Vel(mx+2,mz+1)*tbl_v(3,lx)*tbl_vx(2,lz)+&
				Vel(mx+3,mz+1)*tbl_v(4,lx)*tbl_vx(2,lz)+&
				Vel(mx  ,mz+2)*tbl_v(1,lx)*tbl_vx(3,lz)+&
				Vel(mx+1,mz+2)*tbl_v(2,lx)*tbl_vx(3,lz)+&
				Vel(mx+2,mz+2)*tbl_v(3,lx)*tbl_vx(3,lz)+&
				Vel(mx+3,mz+2)*tbl_v(4,lx)*tbl_vx(3,lz)+&
				Vel(mx  ,mz+3)*tbl_v(1,lx)*tbl_vx(4,lz)+&
				Vel(mx+1,mz+3)*tbl_v(2,lx)*tbl_vx(4,lz)+&
				Vel(mx+2,mz+3)*tbl_v(3,lx)*tbl_vx(4,lz)+&
				Vel(mx+3,mz+3)*tbl_v(4,lx)*tbl_vx(4,lz))/dvz
			
			
!********** else handle end effects with constant extrapolation ********
		ELSE

			v=0.0
			dvdx=0.0
			dvdz=0.0

			DO jx=1,4
				jmx=mx+jx-1
				IF(jmx.lt.1)jmx=1
				IF(jmx.gt.nvx)jmx=nvx

				DO jz=1,4
					jmz=mz+jz-1
					IF(jmz.lt.1)jmz=1
					IF(jmz.gt.nvz)jmz=nvz

					v=v+Vel(jmx,jmz)*tbl_v(jx,lx)*tbl_v(jz,lz)
					dvdx=dvdx+Vel(jmx,jmz)*tbl_vx(jx,lx)*tbl_v(jz,lz)
					dvdz=dvdz+Vel(jmx,jmz)*tbl_v(jx,lx)*tbl_vx(jz,lz)
	
				ENDDO
			ENDDO 
			dvdx=dvdx/dvx
			dvdz=dvdz/dvz

		ENDIF

		RETURN
	END

!***********************************************************************
!                Interpolation coefficient(Bspline)                    !  
!======================================================================*
Subroutine PreInterpolators_Bspline()
      use INTERP_GLOBAL
      INTEGER:: i
      REAL::x
       
DO i=1, NTABLE
	x=(i-1)*1.0/(NTABLE-1)
    tbl_v  (1,i)= (- x*x*x + 3.0 * x*x - 3.0 * x + 1.0) / 6.0
	tbl_v  (2,i)= 0.5 * x*x*x - x*x + 2.0 / 3.0
    tbl_v  (3,i)=-0.5 * x*x*x + 0.5 * x*x + 0.5 * x + 1.0 / 6.0
    tbl_v  (4,i)=x*x*x / 6.0
       
	tbl_vx (1,i)= -0.5 * x*x + x - 0.5
    tbl_vx (2,i)=  1.5 * x*x - 2.0 * x
    tbl_vx (3,i)=-1.5 * x*x + x + 0.5
    tbl_vx (4,i)=0.5 * x*x
ENDDO
END SUBROUTINE PreInterpolators_Bspline
	



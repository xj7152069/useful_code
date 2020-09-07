
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
    PROGRAM main				
		USE INTERP_GLOBAL
        INTEGER NX, NZ, NS_X, NS_Z, NR_X, NR_Z, NVX, NVZ
        REAL    DXS, DZS, DVX, DVZ
		REAL	VX_START,VZ_START
		REAL	DSTEP,DEPTH
        REAL    SX_LEFT,SX_RIGHT

        REAL,ALLOCATABLE::TIME2(:,:)
        REAL,ALLOCATABLE::PATHMAT(:,:)
        REAL,ALLOCATABLE::TIME(:,:)
        
!=====================================PAR=====================================
!MODEL SAMPLE BEGIN WITH (1,1); SOURCE SHOULD BEGIN WITH (2,2)
        NS_X=180
        NS_Z=2
!        NR_X=501
        NR_Z=2

!MODEL SIZE (DEPTH = (NVZ-1)*DVZ)
        NVX=501
        NVZ=250
        DEPTH=2500.0

!SAMPLE GEP
        DVX=10.0
		DVZ=10.0
        DXS=10.0
        DZS=10.0

!SOURCE DISTANCE TO THE LEFT AND RIGHT BOUNDARY
        VX_START=0.0
        VZ_START=0.0

		SX_LEFT=(NS_X-1)*DVX
		SX_RIGHT=(NVX-1)*DVX-SX_LEFT

        DSTEP=10.0

!=============================================================================
        ALLOCATE(TIME2(NVZ,NVX))
        ALLOCATE(PATHMAT(NVZ,NVX))
        
        OPEN(12,FILE = "time.bin" , ACCESS="STREAM" , FORM = "UNFORMATTED" )
        READ( 12 ) TIME2
        CLOSE(12)

		ALLOCATE(TIME(NVX,NVZ))
        DO IX=1, NVX
			DO IZ=1, NVZ
			    TIME(IX,IZ)=TIME2(IZ,IX)
			END DO
		END DO

        open(12,file = "time.ray.bin" , access="stream" , form = "unformatted" )
	    open(14,file='ray.txt',status='replace')
	    open(18,file='pathmat.bin',access="stream", status='replace')
        open(19,file='pathpoint.txt',status='replace')
        open(20,file='ray.bin',access="stream", status='replace')

        NR_Z=1
        DO NR_X=1, 501, 1
            PATHMAT(:,:)=0.0
			call CAL_TRACE_2D(NS_X, NS_Z, NR_X, NR_Z, NVX, NVZ, &
				SX_LEFT, SX_RIGHT, VX_START, VZ_START, &
				DVX, DVZ, DXS, DZS, DEPTH, DSTEP, TIME, TIME2, PATHMAT)
            WRITE( 18 ) PATHMAT
		END DO     
        NR_Z=250    
        DO NR_X=1, 501, 1
            PATHMAT(:,:)=0.0
			call CAL_TRACE_2D(NS_X, NS_Z, NR_X, NR_Z, NVX, NVZ, &
				SX_LEFT, SX_RIGHT, VX_START, VZ_START, &
				DVX, DVZ, DXS, DZS, DEPTH, DSTEP, TIME, TIME2, PATHMAT)
            WRITE( 18 ) PATHMAT
		END DO   
        NR_X=1    
        DO NR_Z=1, 250, 1
            PATHMAT(:,:)=0.0
			call CAL_TRACE_2D(NS_X, NS_Z, NR_X, NR_Z, NVX, NVZ, &
				SX_LEFT, SX_RIGHT, VX_START, VZ_START, &
				DVX, DVZ, DXS, DZS, DEPTH, DSTEP, TIME, TIME2, PATHMAT)
            WRITE( 18 ) PATHMAT
		END DO  
        NR_X=501    
        DO NR_Z=1, 250, 1
            PATHMAT(:,:)=0.0
			call CAL_TRACE_2D(NS_X, NS_Z, NR_X, NR_Z, NVX, NVZ, &
				SX_LEFT, SX_RIGHT, VX_START, VZ_START, &
				DVX, DVZ, DXS, DZS, DEPTH, DSTEP, TIME, TIME2, PATHMAT)
            WRITE( 18 ) PATHMAT
		END DO  

            WRITE( 12 ) TIME2
            

        CLOSE(20)
        CLOSE(14)
		CLOSE(18)
		CLOSE(19)
        close(12)

    END
!==================================================================
    SUBROUTINE CAL_TRACE_2D(NS_X, NS_Z, NR_X, NR_Z, NVX, NVZ, &
				SX_LEFT, SX_RIGHT, VX_START, VZ_START, &
				DVX, DVZ, DXS, DZS, DEPTH, DSTEP, TIME, TIME2, PATHMAT)

		INTEGER NX, NZ, NS_X, NS_Z, NVX, NVZ, NVXS
		INTEGER NXS_LEFT, NXS_RIGHT
		INTEGER ISHOT,MYID
        REAL    DX, DZ, DVX, DVZ
		REAL	VX_START,VZ_START
		REAL	DSTEP,DXS,DZS
		REAL	DS_X, DS_Z
		REAL	SZ_COORD,SX_COORD!,MAXPATH
        REAL    SX_LEFT,SX_RIGHT,DEPTH,NVZS
		INTEGER NTRACES,TRACE_START,TRACE_LOCATE

        DIMENSION::TIME(1:NVX,1:NVZ)
        DIMENSION::TIME2(1:NVZ,1:NVX)
        DIMENSION::PATHMAT(1:NVZ,1:NVX)
!*	START COORDINATE IN X/Z DIMENSION OF THE VELOCITY MODEL	*!
!*	0.0,0.0
!=====================================PAR=====================================
!MODEL SAMPLE BEGIN WITH (1,1); SOURCE SHOULD BEGIN WITH (2,2)
!        NS_X=1 
!        NS_Z=1
!        NR_X=501
!        NR_Z=1

!MODEL SIZE (DEPTH = (NVZ-1)*DVZ)
!        NVX=501
!        NVZ=250
!        DEPTH=2500.0

!SAMPLE GEP
!        DVX=10.0
!		DVZ=10.0
!        DXS=10.0
!        DZS=10.0

!SOURCE DISTANCE TO THE LEFT AND RIGHT BOUNDARY
!		SX_LEFT=10.0
!		SX_RIGHT=4990.0

!        VX_START=0.0
!        VZ_START=0.0

!        DSTEP=10.0

!=============================================================================

        SX_COORD=(NS_X-1)*DVX
		SZ_COORD=(NS_Z-1)*DVZ
        RX_COORD=(NR_X-1)*DVX
		RZ_COORD=(NR_Z-1)*DVZ

		DX=DVX
		DZ=DVZ

		NVXS=CEILING(NVX/(DXS/DVX))
		NVZS=CEILING(NVZ/(DZS/DVZ))

		NXS_LEFT=SX_LEFT/DVX+0.5
		NXS_RIGHT=SX_RIGHT/DVX+0.5

		NX=NXS_LEFT+NXS_RIGHT+1
		NZ=DEPTH/DZ
		
		NV=(NVX-1)*(NVZ-1)
		NVS=(NVXS-1)*(NVZS-1)

		XMIN=0.0
		ZMIN=0.0

!=============================================================================
!		CUR_RAY     =TRIM(ADJUSTL(FN4))//'RAY'//TRIM(ADJUSTL(CUR_IS))//'.TXT'
!		CUR_PATH    =TRIM(ADJUSTL(FN5))//'PATH'//TRIM(ADJUSTL(CUR_IS))//'.DAT'
!		CUR_TIME    =TRIM(ADJUSTL(FN6))//'TRAVELTIME'//TRIM(ADJUSTL(CUR_IS))//'.DAT'
        CALL PREINTERPOLATORS_BSPLINE()

!=============================================================================
!		DO I=1,NTRACES
!			READ(13,REC=I+TRACE_START-1) GEO

			NR_X=NINT((RX_COORD-VX_START)/DVX)+1

!			IF(NR_X.LT.1) CYCLE
!			IF(NR_X.GT.NVX) CYCLE

!			ITR=TRACE_LOCATE+I-1
!            ITR=0
			
!			IF(ITR.EQ.1) THEN

            IF(NR_X>=1 .AND. NR_X<=NVX) THEN
			    CALL TRACING(NX, NZ, DX, DZ, TIME, XMIN, ZMIN, DSTEP,&
						SX_COORD,SZ_COORD,RX_COORD,RZ_COORD,VX_START,&
						NVXS, DXS, DZS, TIME2, PATHMAT)
            END IF
        
!		END DO
!=============================================================================

	END

!===========================================================================
SUBROUTINE TRACING(NX, NZ, DX, DZ, TIME, XMIN, ZMIN, DSTEP,&
				SX_COORD,SZ_COORD,RX_COORD,RZ_COORD,VX_START,&
				NVXS, DXS, DZS, time2, pathmat)
USE INTERP_GLOBAL
IMPLICIT NONE

INTEGER NX, NZ, NVXS
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
                
IF(.NOT.(((rx_coord-sx_coord)*dtdx)>=0.0))THEN
	TMP=-1
	GOTO 1201
ELSE
	WRITE(14,*) 1, X, Z!, 'geophone'            
    WRITE( 20 ) 1.0, X, Z!, 'geophone'                                       !raycoordinate

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
            WRITE( 20 ) 0.0, X, Z

            mx=INT(X/DX)+1
            mz=INT(Z/DZ)+1
            IF(mx>=1 .and. mx<=(nx) .and. mz>=1 .and. mz<=(nz))THEN
                time2(mz,mx)=100                                           
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
			WRITE(14,*) 2 ,sx_coord,sz_coord!,'source'
            WRITE( 20 ) 2.0,sx_coord,sz_coord!,'source'
            mx=INT(sx_coord/DX)+1
            mz=INT(sz_coord/DZ)+1
            IF(mx>=1 .and. mx<=(nx) .and. mz>=1 .and. mz<=(nz) .and. time2(mz,mx)<50)THEN
                time2(mz,mx)=time2(mz,mx)*10                                           
            end if
			EXIT
		END IF
							   

		IF(.NOT.(((x1-sx_coord)*dtdx)>=0.0))THEN
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
WRITE( 20 ) 3.0, TMP, TMP*1000!, 'ray_time(s|ms)'
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
	



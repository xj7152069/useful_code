
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
!	SUBROUTINE TRAVELTIME_2D(SS1,SS2,TT1,TT2,SLOWNESS,TIME,SLOW_45,&
!					EPSILON_VALUE,SSS,VEL,ELEV,NVX,NX,DVX,DVZ,DX,NVZ,&
!					NZ,DZ,NS_X,NS_Z,SX_COORD,VX_START,NXS_LEFT,&
!					NXS_RIGHT,ntraces,trace_start,trace_locate,&
!					SZ_COORD,ISHOT,DSTEP,NVXS,DXS,DZS,myid)

    PROGRAM TRAVELTIME_2D				
		USE INTERP_GLOBAL

		INTEGER NX, NZ, NS_X, NS_Z, NVX, NVZ, NVXS
		INTEGER NXS_LEFT, NXS_RIGHT
		INTEGER ISHOT,myid
        REAL    DX, DZ, DVX, DVZ
		REAL	SX_COORD, VX_START
		REAL	DSTEP,DXS,DZS
		REAL	DS_X, DS_Z
		REAL	SZ_COORD
        real    SX_LEFT,SX_RIGHT,DEPTH,nvzs
		integer ntraces,trace_start,trace_locate

		REAL,ALLOCATABLE::ELEV(:)
		REAL,ALLOCATABLE::SS1(:)
		REAL,ALLOCATABLE::SS2(:)
		REAL,ALLOCATABLE::SSS(:,:)
		REAL,ALLOCATABLE::VEL(:,:)
		REAL,ALLOCATABLE::TT1(:)
		REAL,ALLOCATABLE::TT2(:)
		REAL,ALLOCATABLE::SLOWNESS(:,:)
		REAL,ALLOCATABLE::TIME(:,:)
        REAL,ALLOCATABLE::TIME2(:,:)
		REAL,ALLOCATABLE::SLOW_45(:,:)
		REAL,ALLOCATABLE::EPSILON_VALUE(:,:)

!=====================================Par=====================================
!model sample begin with (1,1); source should begin with (2,2)
        NS_X=35
        NS_Z=2

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

!=============================================================================

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

		ALLOCATE(ELEV(NVX))
		ALLOCATE(SS1(NX+1))
		ALLOCATE(SS2(NX+1))
		ALLOCATE(SSS(NVZ,NVX))
		ALLOCATE(VEL(NVZ,NVX))
		ALLOCATE(TT1(NX))
		ALLOCATE(TT2(NX))
		ALLOCATE(SLOWNESS(NX+1,NZ+1))
		ALLOCATE(TIME(NX,NZ))
        ALLOCATE(TIME2(NZ,NX))
		ALLOCATE(SLOW_45(NX+1,NZ+1))
		ALLOCATE(EPSILON_VALUE(NX+1,NZ+1))

        Open(12,File = "model.dat" , access="stream" , form = "unformatted")
        Read( 12 ) vel
        close(12)

        DO IX=1, NX
			DO IZ=1, NZ
			    SSS(iz,ix)=1.0/(vel(iz,ix)+0*iz) 
!                SSS(iz,ix)=1.0/vel(iz,ix)  
			END DO
		END DO
        
        open(12,File='slowness.bin',Access='stream',Form='Unformatted',status='replace')
        Write( 12 ) SSS
        close(12)
    !===== ZEROING THE WORKING BUFFER
		
        CALL ZERO_BUF(NX, NZ, SS1, SS2, TT1, TT2, TIME,&
				SLOWNESS, SLOW_45) !!1

    !===== READ THE CURRENT-SHOT VELOCITY INTO THE WORKING BUFFER SLOWNESS

        CALL READ_CURRENT_SHOT_VELO(NS_X, NXS_LEFT, NXS_RIGHT, SLOW_45,&
                           NVX, NX, NVZ, NZ, SLOWNESS, SSS, EPSILON_VALUE) !!1

		NS_X = NXS_LEFT+1
		DS_X = NS_X*DX
		DS_Z = NS_Z*DZ 
            
        CALL TRAVELTIME_CAL(SS1, SS2, TT1, TT2, SLOWNESS, TIME, SLOW_45,&
					EPSILON_VALUE, NX, DX, NZ, DZ, NS_X, NS_Z, DS_X, DS_Z)

        DO IX=1, NX
			DO IZ=1, NZ
			    TIME2(iz,ix)=TIME(ix,iz)   
			END DO
		END DO
        open(12,File='time.bin',Access='stream',Form='Unformatted',status='replace')
        Write( 12 ) TIME2
        close( 12 )
!		R1=4.0
!		R2=4.0
!		CALL smooth2f(NZ,NX,R1,R2,TIME)

!        CALL WRITE_DISK(NX, NZ, ISHOT, TIME)

!		CALL RAYTRACING(NX, NZ, NVX, NVZ, DVX, DVZ, DX, DZ, TIME, &
!					ELEV,VEL,SX_COORD,ntraces,trace_start,trace_locate,&
!					VX_START, SZ_COORD, ISHOT, DSTEP, NVXS, DXS,&
!					DZS,myid)
		
		RETURN
    END 

!=============================================================================

    SUBROUTINE TRAVELTIME_CAL(SS1, SS2, TT1, TT2, SLOWNESS, TIME,&
				SLOW_45, EPSILON_VALUE, NX, DX, NY, DY, &
                NS_X, NS_Y, DS_X, DS_Y)
 
		INTEGER NX, NY, NS_X, NS_Y
		REAL    DX, DY, DS_X, DS_Y
       
		DIMENSION SS1(0:NX), SS2(0:NX), TT1(NX), TT2(NX)
		DIMENSION SLOWNESS(0:NX, 0:NY), TIME(NX, NY)
        DIMENSION EPSILON_VALUE(0:NX, 0:NY)
        DIMENSION SLOW_45(0:NX, 0:NY)
        
!===== CALCULATING THE NS_Y LAYER TRAVEL TIME

		CALL START_LAYER(TIME, NX, DX, NY, DY,&
            NS_X, NS_Y, SS1, TT1, SLOWNESS, EPSILON_VALUE)

!===== CALCULATING THE START TIME VALUE

		CALL UPWARD_CALCUL(TIME, NX, DX, NY, DY, &
                 NS_X, NS_Y, DS_X, DS_Y, SS1, SLOW_45,&
                 SS2, TT1, TT2, SLOWNESS, EPSILON_VALUE)       

!===== FORWARD CALCULATING THE MINIMUM TRAVELTIME  

		CALL FORWARD_CALCUL(TIME, NX, DX, NY, DY,&
                  NS_X, NS_Y, DS_X, DS_Y, SS1, SLOW_45,&
                  SS2, TT1, TT2, SLOWNESS, EPSILON_VALUE) 
     
!===== BACKWARD CALCULATING FOR REPLACING THE FORWARD CALCULATED 
!===== MINIMUM TRAVELTIME

		CALL BACKWARD_CALCUL(TIME, NX, DX, NY, DY,&
                   DS_X, DS_Y, SS1, SLOW_45, SS2, TT1,&
                   TT2, SLOWNESS, EPSILON_VALUE)

!===== WRITE THE CALCULATED MINIMUM TRAVELTIME ONTO THE DISKFILE

		RETURN
    END

!=================================================================
	SUBROUTINE START_LAYER(TIME, NX, DX, NY, DY,&
				NS_X, NS_Y, SS1, TT1, SLOWNESS, EPSILON_VALUE)
      
		REAL      DX, DY
		INTEGER   NX, NY, NS_X, NS_Y
		DIMENSION SS1(0:NX), TT1(NX)
		DIMENSION TIME(NX, NY), SLOWNESS(0:NX, 0:NY)
		DIMENSION EPSILON_VALUE(0:NX,0:NY) 
      
!===== CALCULATING THE START TIME VALUE

		DO IX=0, NX
			IF(IX.EQ.0)THEN
				SS1(IX)=SLOWNESS(IX, NS_Y-1)/&
					SQRT(1+2*EPSILON_VALUE(1, NS_Y-1))
			ELSE
				SS1(IX)=SLOWNESS(IX, NS_Y-1)/&
					SQRT(1+2*EPSILON_VALUE(IX, NS_Y-1))
			END IF
		END DO
		TEMP_T=0.0
		DO IX=NS_X-1, 1, -1
			TEMP_T=TEMP_T+DX*SS1(IX)
			TT1(IX)=TEMP_T
		END DO

		TT1(NS_X)=0.0

		TEMP_T=0.0
		DO IX=NS_X+1, NX
			TEMP_T=TEMP_T+DX*SS1(IX-1)
			TT1(IX)=TEMP_T
		END DO

		DO IX=1, NX
			TIME(IX, NS_Y)=TT1(IX)
		END DO

		RETURN
	END

!=================================================================

	SUBROUTINE UPWARD_CALCUL(TIME, NX, DX, NY, DY,&
				NS_X, NS_Y, DS_X, DS_Y, SS1,SLOW_45, SS2, &
				TT1, TT2, SLOWNESS, EPSILON_VALUE )

		REAL    T1, T2, X1, X2, Y1, Y2
		REAL    DX, DY, DS_X, DS_Y ,EPSILON_V, S_PHI
		INTEGER NX, NY, NS_X, NS_Y
		DIMENSION  SS1(0:NX), SS2(0:NX), TT1(NX), TT2(NX)
		DIMENSION TIME(NX, NY), SLOWNESS(0:NX, 0:NY),EPSILON_VALUE(0:NX,0:NY)
		DIMENSION SLOW_45(0:NX, 0:NY)

      
!=================================================
		DO IX=1, NX
			TT1(IX)=TIME(IX, NS_Y)
		END DO

		DO 7777 IY=NS_Y-1 ,1, -1
         
!======  GET THE SLOWNESS OF EACH LAYER
      
			DO IX=0, NX
				SS1(IX)=SLOWNESS(IX, IY)
				SS2(IX)=SLOWNESS(IX, IY-1)
			END DO

			DO IX=1, NX
				TT2(IX)=100000.0
			END DO

			IF(IY.EQ.(NS_Y-1))  THEN
				TEMP_T=AMIN1(SLOWNESS(NS_X-1, NS_Y-1),& 
						SLOWNESS(NS_X, NS_Y-1)) * DY
				TT2(NS_X)=TEMP_T
				TEMP_T=SQRT(DX*DX + DY*DY)*SLOWNESS(NS_X-1, NS_Y-1)
				TT2(NS_X-1)=TEMP_T
				TEMP_T=SQRT(DX*DX + DY*DY)*SLOWNESS(NS_X, NS_Y-1)
				TT2(NS_X+1)=TEMP_T
			END IF
   
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        
			DO IX=2, NX
         
				T1=TT1(IX-1)
				T2=TT1(IX)
				X1=ABS((IX-1)*DX - DS_X)
				X2=ABS(IX*DX     - DS_X)

				IF(IX.EQ.NS_X) THEN
					TT2(NS_X)=T2+AMIN1(SS1(IX-1), SS1(IX))*DY
				END IF

             
				IF(SS1(IX).LT.SS1(IX-1)) THEN     

					TS=T2+SS1(IX)*DY
					IF(TS.LT.TT2(IX)) TT2(IX)=TS
				END IF
        
				SNESS=SS1(IX-1)
           
				W=(T2*T2-T1*T1)/(X2*X2-X1*X1)
				EPSILON_V = EPSILON_VALUE(IX,IY)
				SLOW = SLOW_45(IX,IY)
!           	WRITE(*,*)'EPSILON_V',EPSILON_V,IX,IY
				JJ = - 1
!========================== PHI = PI/2 - PUXI + SITA ==================
				CALL ROOT1(X0, W, SNESS, X1, X2,T1,T2, DY,EPSILON_V,&
							S_PHI, JJ, SLOW)
						
				T0=SQRT(W*(X0*X0 - X1*X1) + T1*T1)
!           	TS=T0+SNESS*SQRT((X2-X0)*(X2-X0) + DY*DY)
				TS=T0+ S_PHI*SQRT((X2-X0)*(X2-X0) + DY*DY)
!           	WRITE(*,*)SNESS,S_PHI
			
				IF(TS.LT.TT2(IX)) TT2(IX)=TS

			END DO

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

			DO IX=NX-1, 1, -1
 
				T1=TT1(IX+1)
				T2=TT1(IX)
				X1=ABS((IX+1)*DX-DS_X)
				X2=ABS(IX*DX-DS_X)

				IF(SS1(IX-1).LT.SS1(IX)) THEN     

					TS=T2+SS1(IX-1)*DY
					IF(TS.LT.TT2(IX)) TT2(IX)=TS

!         		ELSE
				END IF
        
				SNESS=SS1(IX)
           
				W=(T2*T2-T1*T1)/(X2*X2-X1*X1)
				EPSILON_V = EPSILON_VALUE(IX,IY)
				SLOW = SLOW_45(IX,IY)
!========================== PHI = PI/2 - PUXI - SITA ==================
				CALL ROOT1(X0, W, SNESS, X1, X2,T1,T2,DY,EPSILON_V,&
							S_PHI, JJ, SLOW)
						
				T0=SQRT(W*(X0*X0-X1*X1)+T1*T1)
!           	TS=T0+SNESS*SQRT((X2-X0)*(X2-X0)+DY*DY)
				TS=T0+S_PHI*SQRT((X2-X0)*(X2-X0)+DY*DY)
			
				IF(TS.LT.TT2(IX)) TT2(IX)=TS

			END DO

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

			DO IX=2, NX
           
				T1=TT1(IX-1)
				T2=TT2(IX-1)
				Y1=ABS(IY*DY-DS_Y)
				Y2=ABS((IY-1)*DY-DS_Y)

				IF(SS2(IX-1).LT.SS1(IX-1)) THEN

					TS=T2+SS2(IX-1)*DX
					IF(TS.LT.TT2(IX)) TT2(IX)=TS

!        		ELSE
				END IF

				SNESS=SS1(IX-1)
            
				W=(T2*T2-T1*T1)/(Y2*Y2-Y1*Y1)
				EPSILON_V = EPSILON_VALUE(IX,IY) 
				SLOW = SLOW_45(IX,IY)
				JJ = 1  
!========================== PHI =  PUXI + SITA ==================
				CALL ROOT1(Y0, W, SNESS, Y1, Y2,T1,T2,DX,EPSILON_V,&
							S_PHI, JJ, SLOW)
            
				T0=SQRT(W*(Y0*Y0-Y1*Y1)+T1*T1)
!            	TS=T0+SNESS*SQRT((Y2-Y0)*(Y2-Y0)+DX*DX)
				TS=T0+S_PHI*SQRT((Y2-Y0)*(Y2-Y0)+DX*DX)

				IF(TS.LT.TT2(IX)) TT2(IX)=TS
         
			END DO

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

			DO IX=NX-1, 1, -1
          
				T1=TT1(IX+1)
				T2=TT2(IX+1)
				Y1=ABS(IY*DY-DS_Y)
				Y2=ABS((IY-1)*DY-DS_Y)

				IF(SS2(IX).LT.SS1(IX)) THEN

					TS=T2+SS2(IX)*DX
					IF(TS.LT.TT2(IX)) TT2(IX)=TS

!        		ELSE
				END IF
        
				SNESS=SS1(IX)
          
				W=(T2*T2-T1*T1)/(Y2*Y2-Y1*Y1)
				EPSILON_V = EPSILON_VALUE(IX,IY)
				SLOW = SLOW_45(IX,IY)
!========================== PHI =  PUXI + SITA ==================
				CALL ROOT1(Y0, W, SNESS, Y1,Y2,T1, T2, DX,EPSILON_V,&
							S_PHI, JJ, SLOW)
						
				T0=SQRT(W*(Y0*Y0-Y1*Y1)+T1*T1)
!          		TS=T0+SNESS*SQRT((Y2-Y0)*(Y2-Y0)+DX*DX)
				TS=T0+S_PHI*SQRT((Y2-Y0)*(Y2-Y0)+DX*DX)

				IF(TS.LT.TT2(IX)) TT2(IX)=TS

			END DO

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
          
			DO IX=1, NX
				TT1(IX)=TT2(IX)
			END DO

			DO IX=1, NX
				TIME(IX, IY)=TT2(IX)
			END DO

7777  	CONTINUE

		RETURN
	END


!====================================================================

	SUBROUTINE FORWARD_CALCUL(TIME, NX, DX, NY, DY,&
				NS_X, NS_Y, DS_X, DS_Y, SS1,SLOW_45, SS2,&
				TT1, TT2, SLOWNESS, EPSILON_VALUE)

		REAL      T1, T2, X1, X2, Y1, Y2
		REAL      DX, DY, DS_X, DS_Y,EPSILON_V
		INTEGER   NX, NY, NS_X, NS_Y
		DIMENSION SS1(0:NX), SS2(0:NX), TT1(NX), TT2(NX)
		DIMENSION TIME(NX, NY), SLOWNESS(0:NX, 0:NY),EPSILON_VALUE(0:NX, 0:NY)
		DIMENSION SLOW_45(0:NX, 0:NY)

!============================================================= 

		DO IX=1,NX
			TT1(IX)=TIME(IX,1)
		END DO

		DO 8888 IY=2, NY
        
!====== GET THE SLOWNESS OF EACH LAYER
      
			DO IX=0, NX
				SS1(IX)=SLOWNESS(IX, IY-1)
				SS2(IX)=SLOWNESS(IX, IY)
			END DO

!===== ASSIGNING THE (BEYOND) MAXIMUM TRAVELTIME FOR THE
!===== CALCULATING TRAVELTIME

			IF(IY.LE.NS_Y) THEN

				DO IX=1, NX
					TT2(IX)=TIME(IX, IY)
				END DO     

			ELSE

				DO IX=1, NX
					TT2(IX)=1000000.0
				END DO

			ENDIF			 

!====== CALCULATING THE TRAVEL TIME AROUND THE SOURCE
      
			IF(IY.EQ.(NS_Y+1)) THEN

				TEMP_T=DY*AMIN1(SLOWNESS(NS_X, NS_Y),&
						SLOWNESS(NS_X-1, NS_Y))
				TT2(NS_X)=TEMP_T

				TEMP_T=SQRT(DX*DX+DY*DY)*SLOWNESS(NS_X-1, NS_Y)
				TT2(NS_X+1)=TEMP_T

				TEMP_T=SQRT(DX*DX+DY*DY)*SLOWNESS(NS_X, NS_Y) 
				TT2(NS_X-1)=TEMP_T

			END IF

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

			DO IX=2, NX

				T1=TT1(IX-1)
				T2=TT1(IX)
				X1=ABS((IX-1)*DX-DS_X)
				X2=ABS((IX)*DX-DS_X)

				IF(IY.GT.NS_Y+1.AND.IX.EQ.NS_X) THEN
					TT2(NS_X)=T2+AMIN1(SS1(IX-1), SS1(IX))*DY
				ENDIF

				IF(SS1(IX).LT.SS1(IX-1)) THEN     

					TS=T2+SS1(IX)*DY
					IF(TS.LT.TT2(IX)) TT2(IX)=TS

!         ELSE
				END IF

				SNESS=SS1(IX-1)
          
				W=(T2*T2-T1*T1)/(X2*X2-X1*X1)
				EPSILON_V = EPSILON_VALUE(IX,IY)
				SLOW = SLOW_45(IX,IY)
				JJ = - 1
!========================== PHI = PI/2 - PUXI - SITA ==================
				CALL ROOT1(X0, W, SNESS, X1, X2,T1,T2,DY,EPSILON_V,&
							S_PHI, JJ, SLOW)
				T0=SQRT(W*(X0*X0-X1*X1)+T1*T1)
!            	TS=T0+SNESS*SQRT((X2-X0)*(X2-X0)+DY*DY)
				TS=T0+S_PHI*SQRT((X2-X0)*(X2-X0)+DY*DY)


				IF(TS.LT.TT2(IX)) TT2(IX)=TS

			END DO

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

			DO IX=NX-1, 1, -1

				T1=TT1(IX+1)
				T2=TT1(IX)
				X1=ABS((IX+1)*DX-DS_X)
				X2=ABS(IX*DX-DS_X)

				IF(SS1(IX-1).LT.SS1(IX)) THEN     

					TS=T2+SS1(IX-1)*DY
					IF(TS.LT.TT2(IX)) TT2(IX)=TS

!         		ELSE
				END IF
       
				SNESS=SS1(IX)
          
				W=(T2*T2-T1*T1)/(X2*X2-X1*X1)
				EPSILON_V = EPSILON_VALUE(IX,IY)
				SLOW = SLOW_45(IX,IY)
!========================== PHI = PI/2 - PUXI - SITA ==================
				CALL ROOT1(X0, W, SNESS, X1, X2,T1,T2,DY,EPSILON_V,&
							S_PHI,JJ, SLOW)
				T0=SQRT(W*(X0*X0-X1*X1)+T1*T1)
!          		TS=T0+SNESS*SQRT((X2-X0)*(X2-X0)+DY*DY)
				TS=T0+S_PHI*SQRT((X2-X0)*(X2-X0)+DY*DY)


				IF(TS.LT.TT2(IX)) TT2(IX)=TS

			END DO
         
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

			DO IX=2, NX

				T1=TT1(IX-1)
				T2=TT2(IX-1)
				Y1=ABS((IY-1)*DY-DS_Y)
				Y2=ABS((IY)*DY-DS_Y)

				IF(SS2(IX-1).LT.SS1(IX-1)) THEN

					TS=T2+SS2(IX-1)*DX
					IF(TS.LT.TT2(IX)) TT2(IX)=TS

!        		ELSE
				END IF

				SNESS=SS1(IX-1)
           
				W=(T2*T2-T1*T1)/(Y2*Y2-Y1*Y1)
				EPSILON_V = EPSILON_VALUE(IX,IY)
				SLOW = SLOW_45(IX,IY)
				JJ = 1 
!========================== PHI = PUXI - SITA ==================
				CALL ROOT1(Y0, W, SNESS, Y1, Y2,T1,T2,DX,EPSILON_V,&
							S_PHI,JJ, SLOW)
				T0=SQRT(W*(Y0*Y0-Y1*Y1)+T1*T1)
!           	TS=T0+SNESS*SQRT((Y2-Y0)*(Y2-Y0)+DX*DX)
				TS=T0+S_PHI*SQRT((Y2-Y0)*(Y2-Y0)+DX*DX)

				IF(TS.LT.TT2(IX)) TT2(IX)=TS

			END DO
          
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

			DO IX=NX-1, 1, -1

				T1=TT1(IX+1)
				T2=TT2(IX+1)
				Y1=ABS((IY-1)*DY-DS_Y)
				Y2=ABS((IY)*DY-DS_Y)

				IF(SS2(IX).LT.SS1(IX)) THEN

					TS=T2+SS2(IX)*DX
					IF(TS.LT.TT2(IX)) TT2(IX)=TS

!          		ELSE
				END IF

				SNESS=SS1(IX)
           
				W=(T2*T2-T1*T1)/(Y2*Y2-Y1*Y1)
				EPSILON_V = EPSILON_VALUE(IX,IY)
				SLOW = SLOW_45(IX,IY)
!========================== PHI = PUXI - SITA ==================
				CALL ROOT1(Y0, W, SNESS, Y1, Y2,T1,T2,DX,EPSILON_V,&
							S_PHI,JJ, SLOW)
				T0=SQRT(W*(Y0*Y0-Y1*Y1)+T1*T1)
!            	TS=T0+SNESS*SQRT((Y2-Y0)*(Y2-Y0)+DX*DX)
				TS=T0+S_PHI*SQRT((Y2-Y0)*(Y2-Y0)+DX*DX)

				IF(TS.LT.TT2(IX)) TT2(IX)=TS

			END DO
          
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
           
			DO IX=1, NX
				TT1(IX)=TT2(IX)
			END DO

			DO IX=1, NX
				TIME(IX, IY)=TT2(IX)
			END DO

8888   	CONTINUE
 
		RETURN
	END

!===========================================================================
    
	SUBROUTINE BACKWARD_CALCUL(TIME, NX, DX, NY, DY,& 
				DS_X, DS_Y, SS1,SLOW_45, SS2, TT1, TT2, &
				SLOWNESS, EPSILON_VALUE)

		REAL    T1, T2, X1, X2, Y1, Y2
		REAL    DX, DY, DS_X, DS_Y, EPSILON_V
		INTEGER NX, NY
		DIMENSION SS1(0:NX), SS2(0:NX), TT1(NX), TT2(NX)
		DIMENSION TIME(NX, NY), SLOWNESS(0:NX, 0:NY),EPSILON_VALUE(0:NX,0:NY)
		DIMENSION SLOW_45(0:NX, 0:NY)

!==================================================

		DO IX=1, NX
			TT1(IX)=TIME(IX, NY)
		END DO

		DO 9999 IY=NY-1, 1, -1

!======  GET THE SLOWNESS OF EACH LAYER
      
			DO IX=0, NX
				SS1(IX)=SLOWNESS(IX, IY)
				SS2(IX)=SLOWNESS(IX, IY-1)
			END DO

			DO IX=1, NX
				TT2(IX)=TIME(IX, IY)
			END DO

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        
			DO IX=2, NX
        
				T1=TT1(IX-1)
				T2=TT1(IX)
				X1=ABS((IX-1)*DX-DS_X)
				X2=ABS((IX)*DX-DS_X)

				IF(SS1(IX).LT.SS1(IX-1)) THEN     

					TS=T2+SS1(IX)*DY
					IF(TS.LT.TT2(IX)) TT2(IX)=TS

!         		ELSE
				END IF
        
				SNESS=SS1(IX-1)
           
				W=(T2*T2-T1*T1)/(X2*X2-X1*X1)
				EPSILON_V = EPSILON_VALUE(IX,IY)
				SLOW = SLOW_45(IX,IY)
				JJ = - 1
!========================== PHI = PI/2 - PUXI + SITA =================
				CALL ROOT1(X0, W, SNESS, X1, X2,T1,T2,DY,EPSILON_V,&
							S_PHI,JJ ,SLOW)
				T0=SQRT(W*(X0*X0-X1*X1)+T1*T1)
!           	TS=T0+SNESS*SQRT((X2-X0)*(X2-X0)+DY*DY)
				TS=T0+S_PHI*SQRT((X2-X0)*(X2-X0)+DY*DY)

   
				IF(TS.LT.TT2(IX)) TT2(IX)=TS

			END DO

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

			DO IX=NX-1, 1, -1
 
				T1=TT1(IX+1)
				T2=TT1(IX)
				X1=ABS((IX+1)*DX-DS_X)
				X2=ABS((IX)*DX-DS_X)

				IF(SS1(IX-1).LT.SS1(IX)) THEN     

					TS=T2+SS1(IX-1)*DY
					IF(TS.LT.TT2(IX)) TT2(IX)=TS

!         		ELSE
				END IF
        
				SNESS=SS1(IX)
           
				W=(T2*T2-T1*T1)/(X2*X2-X1*X1)
				EPSILON_V = EPSILON_VALUE(IX,IY)
				SLOW = SLOW_45(IX,IY)
!========================== PHI = PI/2 - PUXI + SITA =================
				CALL ROOT1(X0, W, SNESS, X1, X2,T1,T2,DY,EPSILON_V,&
							S_PHI,JJ, SLOW)
				T0=SQRT(W*(X0*X0-X1*X1)+T1*T1)
!           	TS=T0+SNESS*SQRT((X2-X0)*(X2-X0)+DY*DY)
				TS=T0+S_PHI*SQRT((X2-X0)*(X2-X0)+DY*DY)


				IF(TS.LT.TT2(IX)) TT2(IX)=TS

			END DO

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

			DO IX=2, NX
           
				T1=TT1(IX-1)
				T2=TT2(IX-1)
				Y1=ABS(IY*DY-DS_Y)
				Y2=ABS((IY-1)*DY-DS_Y)

				IF(SS2(IX-1).LT.SS1(IX-1)) THEN

					TS=T2+SS2(IX-1)*DX
					IF(TS.LT.TT2(IX)) TT2(IX)=TS

!        		ELSE
				END IF

				SNESS=SS1(IX-1)
            
				W=(T2*T2-T1*T1)/(Y2*Y2-Y1*Y1)
				EPSILON_V = EPSILON_VALUE(IX,IY)
				SLOW = SLOW_45(IX,IY)
				JJ = 1 
!========================== PHI = PUXI + SITA ==================
				CALL ROOT1(Y0, W, SNESS, Y1,Y2,T1,T2, DX,EPSILON_V,&
							S_PHI,JJ, SLOW)
				T0=SQRT(W*(Y0*Y0-Y1*Y1)+T1*T1)
!            	TS=T0+SNESS*SQRT((Y2-Y0)*(Y2-Y0)+DX*DX)
				TS=T0+S_PHI*SQRT((Y2-Y0)*(Y2-Y0)+DX*DX)


				IF(TS.LT.TT2(IX)) TT2(IX)=TS
	 
			END DO

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

			DO IX=NX-1, 1, -1
          
				T1=TT1(IX+1)
				T2=TT2(IX+1)
				Y1=ABS(IY*DY-DS_Y)
				Y2=ABS((IY-1)*DY-DS_Y)

				IF(SS2(IX).LT.SS1(IX)) THEN

					TS=T2+SS2(IX)*DX
					IF(TS.LT.TT2(IX)) TT2(IX)=TS

!        		ELSE
				END IF
        
				SNESS=SS1(IX)
          
				W=(T2*T2-T1*T1)/(Y2*Y2-Y1*Y1)
				EPSILON_V = EPSILON_VALUE(IX,IY)
				SLOW = SLOW_45(IX,IY)
!========================== PHI = PUXI + SITA ==================
				CALL ROOT1(Y0,W,SNESS, Y1, Y2, T1, T2,DX,EPSILON_V,&
							S_PHI,JJ, SLOW)
				T0=SQRT(W*(Y0*Y0-Y1*Y1)+T1*T1)
!          		TS=T0+SNESS*SQRT((Y2-Y0)*(Y2-Y0)+DX*DX)
				TS=T0+S_PHI*SQRT((Y2-Y0)*(Y2-Y0)+DX*DX)

      
				IF(TS.LT.TT2(IX)) TT2(IX)=TS
         
			END DO

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
          
			DO IX=1, NX
				TT1(IX)=TT2(IX)
			END DO

			DO IX=1, NX
				TIME(IX, IY)=TT2(IX)
			END DO

9999  	CONTINUE

		RETURN
	END

!=============================================================================
       
	SUBROUTINE ROOT1(RTBIS,W,S,X1,X2,T1,T2,DY,EPSILON_V,&
				S_PHI,JJ, SLOW_45)

		REAL W, S, T1, T2, DY

		INTEGER JMAX
		REAL RTBIS,X1,X2,XACC, FUNC, EPSILON_V, TS
		EXTERNAL FUNC_NEW
		PARAMETER(XACC=1.0E-8)
		PARAMETER (JMAX=40)
		INTEGER J
		REAL DX,F,FMID,XMID,S_PHI,SLOW_45
		REAL A1, A2 , A3 
      
		A1 = S*S/(1+2*EPSILON_V)
		A3 =  4.0*SLOW_45*SLOW_45 - 2.0*S*S - 2*A1  
		A2 = S*S + A3 - A1

		FMID=FUNC_NEW(X2, W, S, X1, X2, T1, DY,A1,A2,A3,S_PHI,JJ)
		F=FUNC_NEW(X1, W,  S, X1, X2, T1, DY,A1,A2,A3,S_PHI,JJ)

		IF(F*FMID.GE.0.) THEN

!        	RTBIS=X1
			RTBIS=AMIN1(X1, X2)
			IF(JJ.EQ.-1)S_PHI=S
			IF(JJ.EQ.1) THEN
				PHI = ATAN(DY/ABS(X1-X2))
				COS_PHI = COS(PHI)
				COS_PHI2 = COS_PHI*COS_PHI
				SIN_PHI = SIN(PHI)
				S_PHI = SQRT(A1 + A2*COS_PHI2 - A3*COS_PHI2*COS_PHI2)
			ENDIF

			RETURN

		END IF

		IF(F.LT.0.)THEN
			RTBIS=X1
			DX=X2-X1
		ELSE
			RTBIS=X2
			DX=X1-X2
		END IF

		DO 11 J=1,JMAX

			DX=DX*.5
			XMID=RTBIS+DX

			FMID=FUNC_NEW(XMID, W, S, X1, X2, T1, DY,A1,A2,A3,S_PHI,JJ)

			IF(FMID.LE.0.) RTBIS=XMID

			IF(ABS(DX).LT.XACC .OR. FMID.EQ.0.) THEN
				RETURN
			END IF

11		CONTINUE

		STOP 'TOO MANY BISECTIONS IN RTBIS'
	END
!==========================================================================


	REAL FUNCTION FUNC_NEW(Z0, W, S, Z1, Z2, T1, DX,A1,A2,A3,S_PHI,JJ)

		REAL Z0, W, S, Z1, Z2, T1, DX,A1,A2,A3 
		REAL PHI, T0,TEMP, T02, TEMP2
		REAL S_PHI, COS_PHI, COS_PHI2, SIN_PHI, SIN_PHI2
		PHI = ATAN(DX/ABS(Z2-Z0))
		IF(JJ.EQ.-1)PHI = 3.1415926/2 - PHI
		COS_PHI = COS(PHI)
		COS_PHI2 = COS_PHI*COS_PHI
		COS_PHI3 = COS_PHI2*COS_PHI
		SIN_PHI = SIN(PHI)
		SIN_PHI2 = SIN_PHI*SIN_PHI
		SIN_PHI3 = SIN_PHI2*SIN_PHI
		S_PHI = SQRT(A1 + A2*COS_PHI2 - A3*COS_PHI2*COS_PHI2) 
!      	S_PHI = S*SQRT(1-0.408*SIN_PHI2*COS_PHI2-0.378*SIN_PHI**4)
		T0 = SQRT(T1*T1 + W*(Z0*Z0-Z1*Z1))
		TEMP = SQRT(DX*DX+(Z2-Z0)*(Z2-Z0))
		T02 = T1*T1 + W*(Z0*Z0-Z1*Z1)
		TEMP2 = DX*DX+(Z2-Z0)*(Z2-Z0)
      
!      	FUNC_NEW=W*Z0/T0 - S*(Z2-Z0)/TEMP  
		FUNC_NEW=W*Z0/T0 - S_PHI*(Z2-Z0)/TEMP + &
			JJ*(2*A3*COS_PHI3*SIN_PHI-A2*COS_PHI*SIN_PHI)*DX/S_PHI/TEMP
!     + JJ*(S*S*(0.408*COS_PHI*SIN_PHI*(SIN_PHI2-COS_PHI2) - 2*0.378*
!     + SIN_PHI3*COS_PHI)/S_PHI)*DX/TEMP
!      S_PHI = S
	END


!============================================================================
             
	SUBROUTINE ZERO_BUF(NX, NZ, SS1, SS2, TT1, TT2, TIME,&
					SLOWNESS, SLOW_45)

		INTEGER NX, NZ
		DIMENSION SS1(0:NX), SS2(0:NX), TT1(NX), TT2(NX)
		DIMENSION TIME(NX, NZ), SLOWNESS(0:NX, 0:NZ)
		DIMENSION SLOW_45(0:NX, 0:NZ)

!=======================================================================
      
!		WRITE(*,*)  'NX, NZ=', NX, NZ
      
!=======================================================================      
		DO IX=0, NX
			SS1(IX)=0.0
			SS2(IX)=0.0
		END DO

		DO IX=1, NX
			TT1(IX)=0.0
			TT2(IX)=0.0
		END DO

		DO IX=1, NX
			DO IZ=1, NZ
				TIME(IX, IZ)=100000.0
			END DO
		END DO	

		DO IX=0, NX
			DO IZ=0, NZ
				SLOWNESS(IX, IZ)=0.0
				SLOW_45(IX, IZ)=0.0
			END DO
		END DO

		RETURN
	END

!=================================================================

    SUBROUTINE READ_CURRENT_SHOT_VELO(NS_X, NXS_LEFT, NXS_RIGHT, &
				SLOW_45, NVX, NX, NVZ, NZ, SLOWNESS, SSS, EPSILON_VALUE)
		INTEGER NX, NZ, NS_X, NXS_LEFT, NXS_RIGHT
		INTEGER IX, IZ, IXX, IIX
		DIMENSION SLOWNESS(0:NX, 0:NZ)
		DIMENSION SSS(NVZ,NVX)
		DIMENSION EPSILON_VALUE(0:NX,0:NZ)
		DIMENSION SLOW_45(0:NX, 0:NZ)
		
		DO IX=0, NX
			DO IZ=0, NZ
				SLOWNESS(IX, IZ)=0.1   ! ASSIGNNING A SMALL VALUE AVIODING
					                   ! ERROR VELOCITY INPUTTED
			END DO
		END DO
		
		DO IX=0, NX
			DO IZ=0, NZ
				EPSILON_VALUE(IX, IZ)=0.0
			ENDDO
		ENDDO		
		
		IIX=1
		DO IX=NS_X-NXS_LEFT, NS_X+NXS_RIGHT
			IXX=IX
			IF(IXX.LT.1) IXX=1
			IF(IXX.GT.NVX) IXX=NVX
			DO IZ=1,NZ
				SLOWNESS(IIX, IZ)=SSS(IZ,IXX)
				SLOW_45(IIX, IZ)=SSS(IZ,IXX)
			END DO
			IIX=IIX+1
		END DO
	
		DO IX=1, NX
			DO IZ=1, NZ
				IF(IZ.EQ.1)SLOWNESS(IX, IZ-1) = SLOWNESS(IX, IZ)
				IF(IZ.EQ.1) SLOW_45(IX, IZ-1) = SLOW_45(IX, IZ)
				IF(IX.EQ.1)SLOWNESS(IX-1, IZ) = SLOWNESS(IX, IZ)
				IF(IX.EQ.1) SLOW_45(IX-1, IZ) = SLOW_45(IX, IZ)
			END DO
		END DO		
	
	END
!==========================================================================

!===========================================================================

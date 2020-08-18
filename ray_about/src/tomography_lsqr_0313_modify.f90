!===============================================================!
!                                                               !
!             TRAVLE TIME CALCULATION                           !
!                     2 - D                                     ! 
!        AUTHOR  :  HERB. WANG                                  !
!        DATE    :  97,9                                        !
!        MODIFIED BY: HERB. WANG FOR ICLUSTER                   !  
!        DATE    :  08,3                                        !   
!        REFERENCE :                                            !
!               SCHNEIDER, W. A., ET AL , 1992, A DYNAIC        !
!        PROGRAMMING APPROACH TO FIRST ARRIVAL TRAVELTIME       !
!        COMPUTATION IN MEDIA WITH ARBITRARY DISTRIBUTED        !
!        VELOCITIES. GEOPHYSICS, VOL.57, NO.1, P39-50.          ! 
!                                                               !
!     NOTE: THIS PROGRAM CAN CALCULATE THE TRAVELTIME WHEREVER  !
!     THE SHOT POINT STAY.                                      ! 
!===============================================================!


	MODULE INTERP_GLOBAL
		INTEGER,PARAMETER::Ntable=101
		INTEGER,PARAMETER::Ltable=4
		INTEGER ixxx,izzz
		INTEGER lx,lz
		INTEGER mx,mz
		INTEGER nvxm,nvzm
		REAL	ax,az
		REAL	tbl_v(4,Ntable),tbl_vx(4,Ntable)
	END MODULE INTERP_GLOBAL

!========================== MAIN ===============================!

	PROGRAM TOMOGRAPHY_MAIN
		USE lsqrDataModule
		USE lsqrTestModule
		USE INTERP_GLOBAL

		INCLUDE 'mpif.h'

		INTEGER NX, NZ, NVX, NVZ, NVXS, NVXZ
		INTEGER NS_X, NS_X2, NS_Z
		INTEGER SX_LEFT,SX_RIGHT,NXS_LEFT, NXS_RIGHT
		REAL    DX, DZ, DVX, DVZ,DOFFSET
		REAL	DXS,DZS
		REAL    VX_START
		REAL	SX_START
		REAL	SX_COORD
		REAL	XAPERTURE_MAX,DEPTH
		REAL	OFFSETMAX,OFFSETMIN
		INTEGER ISHOT,NSHOT,NOFFSET,IX
		INTEGER ITR, NIT, IIT
		INTEGER(ip) NTR, NV, NVS, NVSS
		REAL	DSHOT
		REAL	DSTEP

		INTEGER ierr, myid, np

		REAL,ALLOCATABLE::ELEV(:)
		REAL,ALLOCATABLE::SS1(:)
		REAL,ALLOCATABLE::SS2(:)
		REAL,ALLOCATABLE::SSS(:,:)
		REAL,ALLOCATABLE::TT1(:)
		REAL,ALLOCATABLE::TT2(:)
		REAL,ALLOCATABLE::SLOWNESS(:,:)
		REAL,ALLOCATABLE::TIME(:,:)
		REAL,ALLOCATABLE::SLOW_45(:,:)
		REAL,ALLOCATABLE::EPSILON_VALUE(:,:)
		
		REAL,ALLOCATABLE::TOBS(:)
		REAL,ALLOCATABLE::DV(:,:)
		REAL,ALLOCATABLE::LIGHT(:,:)
		REAL,ALLOCATABLE::TCAL(:)
		REAL,ALLOCATABLE::J(:)
		REAL,ALLOCATABLE::DSSS(:,:)
		REAL,ALLOCATABLE::DVSS(:,:)
		REAL(dp),ALLOCATABLE::TRES(:)	!b
		REAL(dp),ALLOCATABLE::DSS(:)	!x
		
		CHARACTER(LEN=256) FN1,FN2,FN3,FN4,FN5,FN6,FN7,FN8,FN9,FN10
		CHARACTER(LEN=256) CUR_IS,CUR_IIT,ORI_VEL,CUR_VEL
		CHARACTER(LEN=256) CUR_LIGHT, CUR_LIGHT_IB
		CHARACTER(LEN=256) DSS_FN,DV_FN,DVS_FN,DVSC_FN,DVSM_FN

		!*  MPI START
		CALL MPI_INIT(ierr)
		CALL MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)
		CALL MPI_COMM_SIZE(MPI_COMM_WORLD, np, ierr)


		CALL READPAR(FN1,FN2,FN3,FN4,FN5,FN6,FN7,FN8,FN9,FN10,NVX,NVZ,DVX,DVZ,& 
						NVXS, NVZS, VX_START, XAPERTURE_MAX,DEPTH,&
						OFFSETMAX,OFFSETMIN,SX_START,DOFFSET,NSHOT,DSHOT,NIT,myid)

!		SX_START=3750.0
			CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
		
		CALL PRE_COMPUTED_INTERPOLATORS()
		
		DX=DVX
		DZ=DVZ
		DSTEP=0.5*MIN(DX,DZ)

		IF(NVXS.GT.NVX) THEN
			IF(myid.eq.0) THEN
				PRINT*,'NVXS IS GREATER THAN NVX'
				PRINT*,'SET NVXS = NVX'
			END IF
			NVXS=NVX
		END IF
		IF(NVZS.GT.NVZ) THEN
			IF(myid.eq.0) THEN
				PRINT*,'NVZS IS GREATER THAN NVZ'
				PRINT*,'SET NVZS = NVZ'
			END IF
			NVZS=NVZ
		END IF

		DXS=ceiling(1.0*(NVX-1)/(NVXS-1))*DVX
		DZS=ceiling(1.0*(NVZ-1)/(NVZS-1))*DVZ

		SX_LEFT=XAPERTURE_MAX
		SX_RIGHT=XAPERTURE_MAX

		NXS_LEFT=SX_LEFT/DVX+0.5
		NXS_RIGHT=SX_RIGHT/DVX+0.5

		NX=NXS_LEFT+NXS_RIGHT+1
		NZ=DEPTH/DZ
		
		NV=(NVX-1)*(NVZ-1)
		NVS=(NVXS-1)*(NVZS-1)

		IF(OFFSETMAX.GT.XAPERTURE_MAX) THEN
			IF(myid.eq.0) THEN
				PRINT*,'OFFSETMAX IS GREAT THAN XAPERTURE_MAX'
				PRINT*,'SET OFFSETMAX = XAPERTURE_MAX'
			END IF
			OFFSETMAX = XAPERTURE_MAX
		END IF

		NOFFSET=((OFFSETMAX-OFFSETMIN)/DOFFSET)+1		!Double side receiver
!		NOFFSET=(OFFSETMAX-OFFSETMIN)/DOFFSET+1			!Single side receiver		
		NTR=NOFFSET*NSHOT

		IF(myid.eq.0) THEN
			PRINT*,'DXS=',DXS,'DZS=',DZS
			PRINT*,'NX=',NX,'NZ',NZ
			PRINT*,'NV=',NV,'NVS=',NVS
			PRINT*,'NTR=',NTR
		END IF
		
		ALLOCATE(ELEV(NVX))
		ALLOCATE(SS1(NX+1))
		ALLOCATE(SS2(NX+1))
		ALLOCATE(SSS(NVZ,NVX))
		ALLOCATE(TT1(NX))
		ALLOCATE(TT2(NX))
		ALLOCATE(SLOWNESS(NX+1,NZ+1))
		ALLOCATE(TIME(NX,NZ))
		ALLOCATE(SLOW_45(NX+1,NZ+1))
		ALLOCATE(EPSILON_VALUE(NX+1,NZ+1))
		
		ELEV(1:NVX)=50.0                       !altitude

!		NIT=25
		IIT=1

		IF(myid.eq.0) THEN

			ALLOCATE(TOBS(NTR))	
			ALLOCATE(TCAL(NTR))
			ALLOCATE(TRES(NTR))		!b
			ALLOCATE(J(NIT))
			ALLOCATE(DV(NVZ,NVX))
			ALLOCATE(DSS(NVS))		!x
			ALLOCATE(DSSS(NVZS,NVXS))
			ALLOCATE(DVSS(NVZS,NVXS))
			ALLOCATE(LIGHT(NVZS,NVXS))

			J(1:NIT)=0.0
			
			CALL READ_TIME_OBS(TOBS,FN5,NTR,NSHOT)
                        

			CALL allocatememory()

		END IF
		
		DO WHILE(IIT.LE.NIT)
				
			IF(myid.eq.0) PRINT*,'IIT=',IIT
			
			OPEN(12,FILE=FN2,ACCESS='DIRECT',RECL=NZ,STATUS='REPLACE')!Traveltime Field

			!*	READ Velocity Field For Raytracing
			IF(IIT.EQ.1) THEN
				ORI_VEL=FN1
			ELSE
				ORI_VEL=CUR_VEL
			END IF
			IF(myid.eq.0) PRINT*,'ORI_VEL=',ORI_VEL

			OPEN(11,FILE=ORI_VEL,ACCESS='DIRECT',RECL=NVZ,STATUS='OLD')
			
			CALL READ_INITIAL_SLOWNESS(NVX, NVZ, SSS)
			
			CLOSE(11)
			
			CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)

!			GOTO 1137

			CALL TRAVELTIME_BASED_RAYTRACING(SS1,SS2,TT1,TT2,&
					SLOWNESS,TIME,SLOW_45,EPSILON_VALUE,&
					SSS,ELEV,NVX,NX,DVX,DX,NVZ,NZ,DZ,NS_X,NS_Z,&
					VX_START,SX_START,NSHOT,DSHOT,NXS_LEFT,&
					NXS_RIGHT,OFFSETMAX,OFFSETMIN,NOFFSET,DOFFSET,DSTEP,&
					NVXS,DXS,DZS,FN4,FN7,FN8,myid,np)
			CLOSE(12)

1137		CONTINUE

!			GOTO 1133

			WRITE(CUR_IIT,'(I4)') IIT
			CUR_VEL=trim(adjustl(FN6))//trim(adjustl(CUR_IIT))//'.dat'
			CUR_LIGHT=trim(adjustl(FN3))//trim(adjustl(CUR_IIT))//'.dat'

			CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)

			IF(myid.eq.0) THEN			!Master Process solve equation

				DSS_FN=trim(adjustl(FN10))//'dss'//trim(adjustl(CUR_IIT))//'.dat'
				DVS_FN=trim(adjustl(FN10))//'dvs'//trim(adjustl(CUR_IIT))//'.dat'
				DV_FN=trim(adjustl(FN10))//'dv'//trim(adjustl(CUR_IIT))//'.dat'
				DVSM_FN=trim(adjustl(FN10))//'dvsmoothed'//trim(adjustl(CUR_IIT))//'.dat'
				DVSC_FN=trim(adjustl(FN10))//'dvconsted'//trim(adjustl(CUR_IIT))//'.dat'

				DSS(1:NVS)=0.0

				PRINT*,'Start to read ray path'

				OPEN(21,FILE=CUR_LIGHT,ACCESS='DIRECT',RECL=NVZS,STATUS='REPLACE')
				CALL READ_RAY_DATA(TRES,TCAL,TOBS,NTR,NSHOT,NVSS,&
									NVXS,NVZS,FN7,FN8)

				PRINT*,'Start to illumination'

				CALL ILLUMINATION(NVXS,NVZS,NTR,NVSS,TRES,LIGHT)
				
				CLOSE(21)

				PRINT*,'Start to solve equation'

				CALL solve(NTR,NVSS,DSS,TRES)		!Solve equation by lsqr
                              			
				DO ITR=1,NTR
					J(IIT)=J(IIT)+TCAL(ITR)*TCAL(ITR)
					TRES(ITR)=0.0
					TCAL(ITR)=0.0
				END DO

				PRINT*,'Start to update velocity'

				OPEN(41,FILE=DSS_FN,ACCESS='DIRECT',RECL=NVZS)
				OPEN(42,FILE=DVS_FN,ACCESS='DIRECT',RECL=NVZS)
				OPEN(43,FILE=DV_FN,ACCESS='DIRECT',RECL=NVZ)
				OPEN(44,FILE=DVSM_FN,ACCESS='DIRECT',RECL=NVZ)
				OPEN(45,FILE=DVSC_FN,ACCESS='DIRECT',RECL=NVZ)
			
				CALL UPDATE_SLOWNESS(NV,NVX,NVZ,NVS,NVXS,NVZS,DSS,DSSS,DVSS,SSS,DV)

				CLOSE(41)
				CLOSE(42)
				CLOSE(43)
				CLOSE(44)
				CLOSE(45)

				print*,'Start to write Updated velocity'

				PRINT*,'CUR_VEL=',CUR_VEL	
				OPEN(16,FILE=CUR_VEL,ACCESS='DIRECT',RECL=NVZ,STATUS='REPLACE')

				CALL WRITE_VELOCITY(NVX,NVZ,SSS)

				CLOSE(16)
			
			END IF

			IIT=IIT+1

			CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
			
		END DO

1133	continue

		DEALLOCATE(ELEV)
		DEALLOCATE(SS1)
		DEALLOCATE(SS2)		
		DEALLOCATE(SSS)
		DEALLOCATE(TT1)
		DEALLOCATE(TT2)
		DEALLOCATE(SLOWNESS)
		DEALLOCATE(TIME)
		DEALLOCATE(SLOW_45)
		DEALLOCATE(EPSILON_VALUE)
		
		IF(myid.eq.0) THEN

			OPEN(51,FILE=FN9,ACCESS='DIRECT',RECL=NIT)
			WRITE(51,REC=1)(J(IIT),IIT=1,NIT)
			CLOSE(51)

			DEALLOCATE(TOBS)
			DEALLOCATE(TCAL)	
			DEALLOCATE(TRES)
			DEALLOCATE(J)
			DEALLOCATE(DV)
			DEALLOCATE(DSS)
			DEALLOCATE(DSSS)
			DEALLOCATE(DVSS)
			DEALLOCATE(LIGHT)

			CALL deallocatememory()
		END IF

		CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
			
		CALL MPI_Finalize(ierr)
			
	END

!==================================================================
	
	SUBROUTINE TRAVELTIME_BASED_RAYTRACING(SS1,SS2,TT1,TT2,&
					SLOWNESS,TIME,SLOW_45,EPSILON_VALUE,&
					SSS,ELEV,NVX,NX,DVX,DX,NVZ,NZ,DZ,NS_X,NS_Z,&
					VX_START,SX_START,NSHOT,DSHOT,NXS_LEFT,&
					NXS_RIGHT,OFFSETMAX,OFFSETMIN,NOFFSET,DOFFSET,DSTEP,&
					NVXS,DXS,DZS,FN4,FN7,FN8,myid,np)
	

		USE INTERP_GLOBAL
		INCLUDE 'mpif.h'

		INTEGER NX, NZ, NS_X, NS_X2, NS_Z
		INTEGER NXS_LEFT, NXS_RIGHT
		INTEGER NOFFSET
		INTEGER NSHOT, ISHOT
		INTEGER IV
		INTEGER NVX,NVXS
		REAL	DXS,DZS

		INTEGER myid,np,ntask
		CHARACTER*256 ProcessorName
		INTEGER  status(MPI_STATUS_SIZE)
		INTEGER  isend, itask, inode, ierr
		INTEGER  R_table_Max, S_table_Max
		INTEGER  R_table(2), S_table(2)

        REAL    DX, DZ, DVX, DOFFSET, DSHOT
		REAL	SX_COORD, OFFSETMAX,OFFSETMIN
		REAL	SX_START
		REAL	SZ_COORD
		REAL	DSTEP
		CHARACTER(LEN=256) FN4,FN7,FN8
		CHARACTER(LEN=256) CUR_RAY,CUR_PATH,CUR_TIME,CUR_IS

		DIMENSION ELEV(NVX)
        DIMENSION SS1(0:NX), SS2(0:NX), TT1(NX), TT2(NX)
		DIMENSION SSS(NVZ,NVX)
        DIMENSION SLOWNESS(0:NX,0:NZ)
        DIMENSION SLOW_45(0:NX,0:NZ)
        DIMENSION EPSILON_VALUE(0:NX,0:NZ)

		R_table_Max = 2
		S_table_Max = 2
		ntask=NSHOT

		IF (myid .eq. 0) THEN
!Master Node.
			write(*,*) 'Task number: ', ntask, 'Processor number: ', np

			DO 1111 itask = 1, ntask

				CALL MPI_RECV(R_table, R_table_Max, MPI_INTEGER,&
								MPI_ANY_SOURCE, MPI_ANY_TAG,&
								MPI_COMM_WORLD, status, ierr)

				isend      = status(MPI_SOURCE)
				S_table(1) = 1
				S_table(2) = itask

				CALL MPI_Send(S_table, S_table_Max, MPI_INTEGER,&
								isend, itask,&
								MPI_COMM_WORLD, ierr)

1111     	CONTINUE

			DO 2222 inode = 1, np-1

				CALL MPI_RECV(R_table, R_table_Max, MPI_INTEGER,&
								MPI_ANY_SOURCE, MPI_ANY_TAG,&
								MPI_COMM_WORLD, status, ierr)

				isend      = status(MPI_SOURCE)
				S_table(1) = 0
				S_table(2) = 0

!				write(*,*) 'inode = ', inode, 'isend = ', isend
				CALL MPI_Send(S_table, S_table_Max, MPI_INTEGER,&
								isend, itask,&
								MPI_COMM_WORLD, ierr)

2222     	CONTINUE

			write(*,*) 'Done. myid is ', myid

		ELSE
!Working Node.

			S_table(1) = 0
			S_table(2) = 0
			CALL MPI_Send(S_table, S_table_Max, MPI_INTEGER,&
							0, 0, MPI_COMM_WORLD, ierr)

			DO WHILE(1)

				CALL MPI_RECV(R_table, R_table_Max, MPI_INTEGER,&
								0, MPI_ANY_TAG,&
								MPI_COMM_WORLD, status, ierr)

				IF (R_table(1) .ne. 0) THEN
					S_table(1) = R_table(1)
					S_table(2) = R_table(2)
					itask      = R_table(2)

					write(*,*) 'Working node: ', myid,&
							'received itask: ', itask

					SX_COORD=SX_START+(itask-1)*DSHOT

					NS_X=NINT((SX_COORD-VX_START)/DVX)+1
					NS_X2=NS_X

					IF(NS_X2.LT.1) NS_X2=1
					IF(NS_X2.GT.NVX) NS_X2=NVX

					SZ_COORD = ELEV(NS_X2)
					NS_Z =NINT(ABS(SZ_COORD)/DZ) + 1

					WRITE(*,*)'ISHOT=',itask

					WRITE(CUR_IS,'(I4)') itask
					CUR_RAY =trim(adjustl(FN4))//trim(adjustl(CUR_IS))//'.txt'
					CUR_PATH=trim(adjustl(FN7))//trim(adjustl(CUR_IS))//'.dat'
					CUR_TIME=trim(adjustl(FN8))//trim(adjustl(CUR_IS))//'.dat'

					OPEN(14,FILE=CUR_RAY,STATUS='REPLACE')
					OPEN(17,FILE=CUR_PATH,STATUS='REPLACE')
					OPEN(18,FILE=CUR_TIME,STATUS='REPLACE')

					CALL TRAVELTIME_2D(SS1,SS2,TT1,TT2,SLOWNESS,TIME,SLOW_45,&
							EPSILON_VALUE,SSS,ELEV,NVX,NX,DVX,DX,NVZ,&
							NZ,DZ,NS_X,NS_Z,SX_COORD,VX_START,NXS_LEFT,&
							NXS_RIGHT,OFFSETMAX,OFFSETMIN,NOFFSET,DOFFSET,&
							SZ_COORD,itask,DSTEP,NVXS,DXS,DZS)

					CLOSE(14)
					CLOSE(17)
					CLOSE(18)

					CALL MPI_Send(S_table, S_table_Max, MPI_INTEGER,&
									0, 0, MPI_COMM_WORLD, ierr)
				ELSE
					GOTO 3333
				END IF

			END DO
              
3333     	CONTINUE

			write(*,*) 'Done. myid is ', myid
		END IF

	END

!==================================================================
	SUBROUTINE TRAVELTIME_2D(SS1,SS2,TT1,TT2,SLOWNESS,TIME,SLOW_45,&
					EPSILON_VALUE,SSS,ELEV,NVX,NX,DVX,DX,NVZ,&
					NZ,DZ,NS_X,NS_Z,SX_COORD,VX_START,NXS_LEFT,&
					NXS_RIGHT,OFFSETMAX,OFFSETMIN,NOFFSET,DOFFSET,&
					SZ_COORD,ISHOT,DSTEP,NVXS,DXS,DZS)
					
		USE INTERP_GLOBAL
		INTEGER NX, NZ, NS_X, NS_Z, NVX, NVZ, NVXS
		INTEGER NXS_LEFT, NXS_RIGHT
		INTEGER NOFFSET
		INTEGER ISHOT
        REAL    DX, DZ, DVX, DOFFSET
		REAL	SX_COORD, VX_START, OFFSETMAX,OFFSETMIN
		REAL	DSTEP,DXS,DZS
		REAL	DS_X, DS_Z
		REAL	SZ_COORD
		REAL	R1,R2

		DIMENSION ELEV(NVX)
        DIMENSION SS1(0:NX), SS2(0:NX), TT1(NX), TT2(NX)
        DIMENSION SLOWNESS(0:NX, 0:NZ), TIME(NX, NZ) 
        DIMENSION EPSILON_VALUE(0:NX, 0:NZ)
        DIMENSION SLOW_45(0:NX, 0:NZ)    
        DIMENSION SSS(NVZ,NVX)
    !===== ZEROING THE WORKING BUFFER
		
        CALL ZERO_BUF(NX, NZ, SS1, SS2, TT1, TT2, TIME,&
				SLOWNESS, SLOW_45)

    !===== READ THE CURRENT-SHOT VELOCITY INTO THE WORKING BUFFER SLOWNESS
                
        CALL READ_CURRENT_SHOT_VELO(NS_X, NXS_LEFT, NXS_RIGHT, SLOW_45,&
                           NVX, NX, NVZ, NZ, SLOWNESS, SSS, EPSILON_VALUE)
                
		NS_X = NXS_LEFT+1
		DS_X = NS_X*DX
		DS_Z = NS_Z*DZ
              
                       
        CALL TRAVELTIME_CAL(SS1, SS2, TT1, TT2, SLOWNESS, TIME, SLOW_45,&
					EPSILON_VALUE, NX, DX, NZ, DZ, NS_X, NS_Z, DS_X, DS_Z)

		R1=4.0
		R2=4.0
		CALL smooth2f(NZ,NX,R1,R2,TIME)

        CALL WRITE_DISK(NX, NZ, ISHOT, TIME)
			
		CALL RAYTRACING(NX, NZ, NVX, DVX, DX, DZ, TIME, &
					ELEV, NOFFSET, DOFFSET, OFFSETMAX,OFFSETMIN,SX_COORD,&
					VX_START, SZ_COORD, ISHOT, DSTEP, NVXS, DXS, DZS)
		
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

		PAUSE 'TOO MANY BISECTIONS IN RTBIS'
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
		
		DO IX=1, NX
			DO IZ=1, NZ
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
!======  WRITE THE CALCULATED TRAVEL TIME ONTO THE DISK FILE FN2

	SUBROUTINE WRITE_DISK(NX, NZ, ISHOT, TIME)

		INTEGER NX, NZ, KREC_SAVE, ISHOT
		DIMENSION TIME(NX, NZ)

		KREC_SAVE=(ISHOT-1)*NX

		DO IX=1, NX
			WRITE(12, REC=KREC_SAVE+IX) (TIME(IX, IZ), IZ=1, NZ)
		END DO

       RETURN
    END

!===========================================================================

	SUBROUTINE READPAR(FN1,FN2,FN3,FN4,FN5,FN6,FN7,FN8,FN9,FN10,NVX,NVZ,DVX,DVZ,& 
						NVXS, NVZS, VX_START, XAPERTURE_MAX,DEPTH,&
						OFFSETMAX,OFFSETMIN,SX_START,DOFFSET,NSHOT,DSHOT,NIT,myid)

		INTEGER NVX,NVZ,NSHOT
		INTEGER NVXS, NVZS
		INTEGER NIT,myid
		REAL	DSHOT
		REAL	DVX,DVZ,DOFFSET,SX_START
		REAL    XAPERTURE_MAX,DEPTH,OFFSETMAX,OFFSETMIN
		REAL    VX_START
		CHARACTER(256) FN1,FN2,FN3,FN4,FN5,FN6,FN7,FN8,FN9,FN10,PAR_JUMPPER

		OPEN(11, FILE='modify_TOMOGRAPHY.PAR')

		READ(11,'(A)') PAR_JUMPPER
		READ(11,'(A)') FN1
		
		READ(11,'(A)') PAR_JUMPPER
		READ(11,'(A)') FN2
		
		READ(11,'(A)') PAR_JUMPPER
		READ(11,'(A)') FN3		

		READ(11,'(A)') PAR_JUMPPER
		READ(11,'(A)') FN4

		READ(11,'(A)') PAR_JUMPPER
		READ(11,'(A)') FN5

		READ(11,'(A)') PAR_JUMPPER
		READ(11,'(A)') FN6

		READ(11,'(A)') PAR_JUMPPER
		READ(11,'(A)') FN7

		READ(11,'(A)') PAR_JUMPPER
		READ(11,'(A)') FN8

		READ(11,'(A)') PAR_JUMPPER
		READ(11,'(A)') FN9
           
                READ(11,'(A)') PAR_JUMPPER
                READ(11,'(A)') FN10

		READ(11,'(A)') PAR_JUMPPER
		READ(11,*) NVX,NVZ

		READ(11,'(A)') PAR_JUMPPER
		READ(11,*) DVX,DVZ

		READ(11,'(A)') PAR_JUMPPER
		READ(11,*) NVXS,NVZS
		
		READ(11,'(A)') PAR_JUMPPER
		READ(11,*) VX_START 

                READ(11,'(A)') PAR_JUMPPER
		READ(11,*) XAPERTURE_MAX, DEPTH

		READ(11,'(A)') PAR_JUMPPER
		READ(11,*) OFFSETMIN,OFFSETMAX,DOFFSET

		READ(11,'(A)') PAR_JUMPPER
		READ(11,*) NSHOT,DSHOT

		READ(11,'(A)') PAR_JUMPPER
		READ(11,*) SX_START
           
                READ(11,'(A)') PAR_JUMPPER
                READ(11,*) NIT

		IF(myid.eq.0) THEN

			WRITE(*,*) 'THE INITIAL VELOCITY'
			WRITE(*,'(A)') FN1

			WRITE(*,*) 'THE TRAVELTIME FILENAME'
			WRITE(*,'(A)') FN2

			WRITE(*,*) 'THE LIGHT DISTRIBUTION FILENAME'
			WRITE(*,'(A)') FN3

			WRITE(*,*) 'THE RAYTRACING RESULT'
			WRITE(*,'(A)') FN4

			WRITE(*,*) 'THE PICKED TRAVELTIME FILENAME'
			WRITE(*,'(A)') FN5

			WRITE(*,*) 'THE UPDATED VELOCITY'
			WRITE(*,'(A)') FN6
		        
			WRITE(*,*) 'THE RAY PATH FILENAME'
			WRITE(*,'(A)') FN7

			WRITE(*,*) 'THE CALCULATE TRACELTIME'
			WRITE(*,'(A)') FN8

			WRITE(*,*) 'THE CONVERGENCE CURVE FILENAME'
			WRITE(*,'(A)') FN9
                        
                        WRITE(*,*) 'THE LSQR PARAMETER(DSS,DVS,DV,DVSM,DVSC)'
		        WRITE(*,'(A)') FN10
                           
			WRITE(*,*) 'NX,NZ OF THE VELOCITY MODEL'
			WRITE(*,*) NVX,NVZ

			WRITE(*,*) 'DX,DZ OF THE VELOCITY MODEL'
			WRITE(*,*) DVX,DVZ

			WRITE(*,*) 'SOLVING EQUATION NX AND NZ'
			WRITE(*,*) NVXS, NVZS

			WRITE(*,*) 'VELOCITY START COORDINATE IN X DIMENSION'
			WRITE(*,*) VX_START

			WRITE(*,*) 'XAPERTURE_MAX AND DEPTH'
			WRITE(*,*) XAPERTURE_MAX, DEPTH

			WRITE(*,*) 'MIN OFFSET MAX OFFSET AND DOFFSET'
			WRITE(*,*) OFFSETMIN,OFFSETMAX,DOFFSET

			WRITE(*,*) 'NSHOT DSHOT'
			WRITE(*,*) NSHOT,DSHOT

			WRITE(*,*) 'SX_START'
			WRITE(*,*)  SX_START

                        WRITE(*,*) 'INTERATIONS'
                        WRITE(*,*) NIT

		END IF
		
		CLOSE(11)

	END

!===========================================================================

	SUBROUTINE RAYTRACING(NX, NZ, NVX, DVX, DX, DZ, TIME, &
					ELEV, NOFFSET, DOFFSET, OFFSETMAX,OFFSETMIN,SX_COORD,&
					VX_START, SZ_COORD, ISHOT, DSTEP, NVXS, DXS, DZS)
		USE INTERP_GLOBAL
		INTEGER NX, NZ, NVX, NVXS
		INTEGER NOFFSET, NOFFSET_H_R,NOFFSET_H_L, IOFFSET
		INTEGER NR_X
		INTEGER ISHOT
		REAL	DX, DZ, DVX, DOFFSET, SX_COORD
		REAL	VX_START, XMIN, ZMIN
		REAL	RX_COORD, RZ_COORD
		REAL	DSTEP,DXS,DZS
		REAL	SZ_COORD,OFFSETMIN,OFFSETMAX
		DIMENSION TIME(NX, NZ)
		DIMENSION ELEV(NVX)
		INTEGER	ITR

!		XMIN ZMIN:MINIMUM COORDINATE OF THE TIME PROFILE
		XMIN=SX_COORD-(NX-1)*DX/2
		ZMIN=0.0		

		DO IOFFSET=1,NOFFSET
		
!			PRINT*,'IOFFSET=',IOFFSET
      RX_COORD=SX_COORD+OFFSETMIN+DOFFSET*(IOFFSET-1)	!Double side
!			RX_COORD=SX_COORD+DOFFSET*(IOFFSET-1)+OFFSETMIN		!Single side receiver
			NR_X=NINT((RX_COORD-VX_START)/DVX)+1

			IF(NR_X.LT.1) CYCLE                             
			IF(NR_X.GT.NVX) CYCLE
			RZ_COORD=ELEV(NR_X)

			ITR=(ISHOT-1)*NOFFSET+IOFFSET

			CALL TRACING(NX, NZ, DX, DZ, TIME, XMIN, ZMIN, DSTEP,&
						SX_COORD,SZ_COORD,RX_COORD,RZ_COORD,VX_START,&
						NVXS, DXS, DZS, ITR)

		END DO

	END

!===========================================================================

	SUBROUTINE TRACING(NX, NZ, DX, DZ, TIME, XMIN, ZMIN, DSTEP,&
						SX_COORD,SZ_COORD,RX_COORD,RZ_COORD,VX_START,&
						NVXS, DXS, DZS, ITR)
		USE INTERP_GLOBAL
		INTEGER NX, NZ, NVXS

		INTEGER ITR
		REAL	DX, DZ
		REAL	XMIN, ZMIN
		REAL	SX_COORD,SZ_COORD,RX_COORD,RZ_COORD
		REAL	VX_START
		REAL	TMP
		REAL	DSTEP
		REAL	DXS,DZS
		DIMENSION TIME(NX, NZ)

		REAL	X, Z, T, DTDX, DTDZ	!CURRENT POINT PARAMETER
		REAL	X1,Z1,X2,Z2,X3,Z3,X4,Z4
		INTEGER IVX1,IVZ1,IVX2,IVZ2,IVX,IVZ,IV
		REAL	PATH1,PATH2,PATH,TMP1,TMP2,TMP3

		X=RX_COORD
		Z=RZ_COORD

		PATH=0.0
		PATH1=0.0
		PATH2=0.0

		CALL VELINTERPED(X, Z, T, XMIN, ZMIN, DX, DZ,&
						DTDX, DTDZ, NX, NZ, TIME)
!		PRINT*,'DTDX',DTDX,'DTDZ',DTDZ,'DTDZ/DTDX',DTDZ/DTDX
		
		TMP=T
!		WRITE(18,*) ITR,TMP
!		GOTO 1201

!		IF((DTDZ.GT.0).OR.(ABS(DTDZ/DTDX).LT.0.05)) THEN 
		IF(DTDZ.GT.0) THEN 
			TMP=-1.0
			GOTO 1201
		ELSE

			WRITE(14,*) X, Z                                                    !raycoordinate

			X1=X
			Z1=Z
			IVX1=INT((X1-VX_START)/DXS)+1
			IVZ1=INT(Z1/DZS)+1

	
!			DO WHILE (T.GE.0.001)	!IF CURRENT IS SHOT POINT,T=0
			DO WHILE(1)

				X=X-DSTEP*DTDX/SQRT(DTDX*DTDX+DTDZ*DTDZ)
				Z=Z-DSTEP*DTDZ/SQRT(DTDX*DTDX+DTDZ*DTDZ)
				
				TMP3=SQRT((X-SX_COORD)*(X-SX_COORD)+&
						(Z-SZ_COORD)*(Z-SZ_COORD))

				IF((Z.LE.0.0).AND.(TMP3.GT.DSTEP)) THEN
					TMP=-1
					GOTO 1201
				END IF
				
				WRITE(14,*) X, Z                                           !raycoordinate

				X2=X
				Z2=Z
				IVX2=INT((X2-VX_START)/DXS)+1
				IVZ2=INT(Z2/DZS)+1

				IF((IVX1.EQ.IVX2).AND.(IVZ1.EQ.IVZ2)) THEN
					PATH=PATH+DSTEP
				ELSE IF(IVX1.EQ.IVX2) THEN
					PATH1=ABS((Z1-(MAX(IVZ1,IVZ2)-1)*DZS)*&
						SQRT(DTDX*DTDX+DTDZ*DTDZ)/DTDZ)
					PATH=PATH+PATH1
					IVX=IVX1
					IVZ=IVZ1
					IV=(IVZ-1)*(NVXS-1)+IVX	
					WRITE(17,*) ITR, IV, PATH
					PATH2=DSTEP-PATH1
					PATH=PATH2
				ELSE IF(IVZ1.EQ.IVZ2) THEN
					PATH1=ABS((X1-(MAX(IVX1,IVX2)-1)*DXS-VX_START)*&
						SQRT(DTDX*DTDX+DTDZ*DTDZ)/DTDX)
					PATH=PATH+PATH1
					IVX=IVX1
					IVZ=IVZ1
					IV=(IVZ-1)*(NVXS-1)+IVX	
					WRITE(17,*) ITR, IV, PATH
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
					WRITE(17,*) ITR, IV, PATH

					PATH=MAX(TMP1,TMP2)-MIN(TMP1,TMP2)

					X3=X-TMP1*DTDX/SQRT(DTDX*DTDX+DTDZ*DTDZ)
					Z3=Z-TMP1*DTDZ/SQRT(DTDX*DTDX+DTDZ*DTDZ)
					X4=X-TMP2*DTDX/SQRT(DTDX*DTDX+DTDZ*DTDZ)
					Z4=Z-TMP2*DTDZ/SQRT(DTDX*DTDX+DTDZ*DTDZ)

					IVX=INT((0.5*(X3+X4)-VX_START)/DXS)+1
					IVZ=INT(0.5*(Z3+Z4)/DZS)+1
					IV=(IVZ-1)*(NVXS-1)+IVX	
					WRITE(17,*) ITR, IV, PATH

					PATH2=DSTEP-MAX(TMP1,TMP2)
					PATH=PATH2

				END IF

				CALL VELINTERPED(X, Z, T, XMIN, ZMIN, DX, DZ,&
							DTDX, DTDZ, NX, NZ, TIME)
!				PRINT*,'X,Z',X,Z
!				PRINT*,'DTDX',DTDX,'DTDZ',DTDZ

				X1=X
				Z1=Z
				IVX1=INT((X1-VX_START)/DXS)+1
				IVZ1=INT(Z1/DZS)+1

				IF(TMP3.LT.DSTEP) THEN
				!	PATH=PATH+TMP3
					EXIT
				END IF

			END DO

			IVX=IVX1
			IVZ=IVZ1
			IV=(IVZ-1)*(NVXS-1)+IVX	
			WRITE(17,*) ITR, IV, PATH

		END IF

1201	CONTINUE

		WRITE(18,*) ITR, TMP

		RETURN
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
		IF(mx.gt.0.and.mx.le.nvxm.and.mz.gt.0.and.mz.le.nvzm)THEN
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

	
	
!********** pre-computed interpolators, for 0th, 1st, and 2nd derivatives *****

	SUBROUTINE PRE_COMPUTED_INTERPOLATORS()
		USE INTERP_GLOBAL

		DO i=1, Ntable
			x=(i-1)*1.0/(Ntable-1)
			tbl_v(1,i)=-0.5*x*x*x+x*x-0.5*x
			tbl_v(2,i)=1.5*x*x*x-2.5*x*x+1
			tbl_v(3,i)=-1.5*x*x*x+2*x*x+0.5*x
			tbl_v(4,i)=0.5*x*x*x-0.5*x*x
       
			tbl_vx(1,i)=-1.5*x*x+2*x-0.5
			tbl_vx(2,i)=4.5*x*x-5*x
			tbl_vx(3,i)=-4.5*x*x+4*x+0.5
			tbl_vx(4,i)=1.5*x*x-x

		ENDDO

		RETURN
	END
	
!***********************************************************************!	
	SUBROUTINE READ_INITIAL_SLOWNESS(NVX, NVZ, SSS)
		INTEGER NVX, NVZ
		INTEGER IVX, IVZ
		REAL	SSS(NVZ,NVX)
		
		DO IVX=1,NVX
			READ(11,REC=IVX)(SSS(IVZ,IVX),IVZ=1,NVZ)
		END DO

		DO IVX=1,NVX
			DO IVZ=1,NVZ
				SSS(IVZ,IVX)=1.0/SSS(IVZ,IVX)
			END DO
		END DO
			
	END

!***********************************************************************!
	SUBROUTINE UPDATE_SLOWNESS(NV,NVX,NVZ,NVS,NVXS,NVZS,DSS,DSSS,DVSS,SSS,DV)
		USE lsqrDataModule
		USE lsqrTestModule
		USE INTERP_GLOBAL

		INTEGER NVX, NVZ
		INTEGER NVXS,NVZS
		INTEGER IVX,IVZ,IVXS,IVZS,IV
		INTEGER IIVX,IIVZ
		INTEGER(ip) NV,NVS
		INTEGER DXS,DZS
		REAL	X,Z,DX,DZ,V,DTDX,DTDZ
		REAL	DSSS(NVZS,NVXS)
		REAL	DVSS(NVZS,NVXS)
		REAL	R1,R2
		REAL	SSS(NVZ,NVX)
		REAL	DV(NVZ,NVX)
		REAL(dp) DSS(NVS)
		REAL,ALLOCATABLE::TMP(:,:)

		ALLOCATE(TMP(NVXS,NVZS))

		DXS=ceiling(1.0*(NVX-1)/(NVXS-1))
		DZS=ceiling(1.0*(NVZ-1)/(NVZS-1))

		!*	1D slowness increment to 2D
		DO IVZ=1,NVZS-1
			DO IVX=1,NVXS-1
				IV=(IVZ-1)*(NVXS-1)+IVX	
				DSSS(IVZ,IVX)=DSS(IV)	
			END DO
		END DO
		DO IVZ=1,NVZS-1
			DSSS(IVZ,NVXS)=DSSS(IVZ,NVXS-1)
		END DO
		DO IVX=1,NVXS
			DSSS(NVZS,IVX)=DSSS(NVZS-1,IVX)
		END DO

		DO IVX=1,NVXS
			WRITE(41,REC=IVX)(DSSS(IVZ,IVX),IVZ=1,NVZS)
		END DO

		!*	sampling original slowness to NVZS*NVXS
		DO IVX=1,NVXS-1
			IIVX=(IVX-1)*DXS+1
			DO IVZ=1,NVZS-1
				IIVZ=(IVZ-1)*DZS+1
				DVSS(IVZ,IVX)=SSS(IIVZ,IIVX)
			END DO
		END DO
		DO IVZ=1,NVZS-1
			IIVZ=(IVZ-1)*DZS+1
			DVSS(IVZ,NVXS)=SSS(IIVZ,NVX)
		END DO
		DO IVX=1,NVXS-1
			IIVX=(IVX-1)*DXS+1
			DVSS(NVZS,IVX)=SSS(NVZ,IIVX)
		END DO
		DVSS(NVZS,NVXS)=SSS(NVZ,NVX)

		!*	compute velocity increment
		DO IVX=1,NVXS
			DO IVZ=1,NVZS
				DSSS(IVZ,IVX)=1.0/(DSSS(IVZ,IVX)+DVSS(IVZ,IVX))
				DVSS(IVZ,IVX)=1.0/(DVSS(IVZ,IVX))
				DVSS(IVZ,IVX)=DSSS(IVZ,IVX)-DVSS(IVZ,IVX)
			END DO
		END DO

		DO IVX=1,NVXS
			WRITE(42,REC=IVX)(DVSS(IVZ,IVX),IVZ=1,NVZS)
		END DO

		!*	smooth velocity increment

		R1=20.0
		R2=20.0
		CALL smooth2f(NVXS,NVZS,R1,R2,DVSS)

		DO IVX=1,NVXS
			DO IVZ=1,NVZS
				TMP(IVX,IVZ)=DVSS(IVZ,IVX)
			END DO
		END DO

		!*	Interplotion velocity increment to NVZ*NVX
		XMIN=1.0
		ZMIN=1.0
		DX=DXS*1.0
		DZ=DZS*1.0
		DO IVX=1,NVX
			DO IVZ=1,NVZ
				X=IVX*1.0
				Z=IVZ*1.0
				CALL VELINTERPED(X,Z,V,XMIN,ZMIN,DX,DZ,&
								DTDX,DTDZ,NVXS,NVZS,TMP)
				DV(IVZ,IVX)=V	
			END DO
		END DO

		DO IVX=1,NVX
			WRITE(43,REC=IVX)(DV(IVZ,IVX),IVZ=1,NVZ)
		END DO

		!*	smooth velocity increment

		R1=5.0
		R2=5.0
		CALL smooth2f(NVX,NVZ,R1,R2,DV)

		DO IVX=1,NVX
			WRITE(44,REC=IVX)(DV(IVZ,IVX),IVZ=1,NVZ)
		END DO

		!*	Constraint velocity increment 
		DO IVX=1,NVX
			DO IVZ=1,NVZ
!				IF(ABS(DV(IVZ,IVX)).GT.300.0) THEN
!					DV(IVZ,IVX)=300.0*DV(IVZ,IVX)/ABS(DV(IVZ,IVX))
!				END IF
				DV(IVZ,IVX)=0.3*DV(IVZ,IVX)                        !0.3=update step-size(trying times)
			END DO
		END DO

		DO IVX=1,NVX
			WRITE(45,REC=IVX)(DV(IVZ,IVX),IVZ=1,NVZ)
		END DO

		!*	update velocity
		DO IVX=1,NVX
			DO IVZ=1,NVZ
				SSS(IVZ,IVX)=DV(IVZ,IVX)+1.0/SSS(IVZ,IVX)
			END DO
		END DO

		DEALLOCATE(TMP)

		RETURN
	END
!***********************************************************************!		
	SUBROUTINE WRITE_VELOCITY(NVX,NVZ,SSS)		
		INTEGER NVX, NVZ
		INTEGER IVX, IVZ
		DIMENSION SSS(NVZ,NVX)
					
		DO IVX=1,NVX
			WRITE(16,REC=IVX)(SSS(IVZ,IVX),IVZ=1,NVZ)
		END DO
		RETURN
	END	

!***********************************************************************!		
	SUBROUTINE READ_RAY_DATA(TRES,TCAL,TOBS,NTR,NSHOT,NVSS,&
								NVXS,NVZS,FN7,FN8)

		USE lsqrDataModule
		USE lsqrTestModule

		REAL(dp)	TRES(NTR)
		REAL		TOBS(NTR)
		REAL		TCAL(NTR)
		INTEGER(ip)	NTR, NV, IV, NVSS
		INTEGER		NVXS,NVZS,IVX,IVZ
		INTEGER		ITR,ISHOT,NSHOT
		INTEGER tmpi,tmpj
		REAL	tmpp,tmpt
		
		CHARACTER(LEN=256) FN7,FN8,CUR_IS,CUR_RAY,CUR_PATH,CUR_TIME

!		PRINT*,'NSHOT=',NSHOT
!		PRINT*,'FN7=',FN7
!		PRINT*,'FN8=',FN8
!		PRINT*,'NVXS=',NVXS
!		PRINT*,'NVZS=',NVZS

		!*	initialize array
		CALL initialization()

		NVSS=0
		ITR=0

		DO ISHOT=1,NSHOT
			WRITE(CUR_IS,'(I4)') ISHOT
			CUR_TIME=trim(adjustl(FN8))//trim(adjustl(CUR_IS))//'.dat'
			OPEN(88,FILE=CUR_TIME,STATUS='OLD')
			DO WHILE(.NOT.EOF(88))
!				ITR=ITR+1
!				READ(88,*) TCAL(ITR)
				READ(88,*) ITR,tmpt
				TCAL(ITR)=tmpt
			END DO
			CLOSE(88)
		END DO

		DO ISHOT=1,NSHOT
			WRITE(CUR_IS,'(I4)') ISHOT
			CUR_PATH=trim(adjustl(FN7))//trim(adjustl(CUR_IS))//'.dat'
			OPEN(77,FILE=CUR_PATH,STATUS='OLD')
			DO WHILE(.NOT.EOF(77))
				READ(77,*) tmpi,tmpj,tmpp
				IF(TCAL(tmpi).NE.-1) THEN
					nzz=nzz+1
					IF(nzz.GT.maxnz) THEN
						PRINT*,'nnz is bigger than maxnz'
						STOP
					END IF
					indx(nzz)=tmpi
					jndx(nzz)=tmpj
					path(nzz)=tmpp
					IF(NVSS.LT.tmpj) NVSS=tmpj
				END IF
			END DO
			CLOSE(77)
		END DO
		
		DO ITR=1,NTR
			IF(TCAL(ITR).EQ.-1) THEN
				TCAL(ITR)=0.0
			ELSE
				TCAL(ITR)=TOBS(ITR)-TCAL(ITR)
			END IF
			TRES(ITR)=TCAL(ITR)
		END DO
		
		RETURN
	END

!***********************************************************************!		
	SUBROUTINE READ_TIME_OBS(TOBS,FN5,NTR,NSHOT)
		USE lsqrDataModule
		USE lsqrTestModule
		INTEGER(ip) NTR
		INTEGER		NSHOT,ISHOT
		REAL		TMP
		REAL		TOBS(NTR)
		CHARACTER(LEN=256) FN5,CUR_IS,CUR_TIME
		ITR=0		
		DO ISHOT=1,NSHOT
			WRITE(CUR_IS,'(I4)') ISHOT
			CUR_TIME=trim(adjustl(FN5))//trim(adjustl(CUR_IS))//'.dat'
!			PRINT*,'ISHOT=',ISHOT
!			PRINT*,'CUR_TIME=',CUR_TIME

			OPEN(55,FILE=CUR_TIME,STATUS='OLD')
			DO WHILE(.not. eof(55))
!				ITR=ITR+1
!				READ(55,*) TOBS(ITR)
				READ(55,*) ITR,TMP
				TOBS(ITR)=TMP
			END DO
			CLOSE(55)
		END DO
	END
!***********************************************************************!		
	SUBROUTINE ILLUMINATION(NVXS,NVZS,NTR,NVSS,TRES,LIGHT)
		USE lsqrDataModule
		USE lsqrTestModule

		INTEGER(ip)	NTR, NVSS, izz
		INTEGER		NVXS,NVZS,IVX,IVZ
		REAL(dp)	TRES(NTR)
		REAL		LIGHT(NVZS,NVXS)
		INTEGER		tmpi,tmpj
		REAL		tmpp

		LIGHT(1:NVZS,1:NVXS)=0.0

		DO izz=1,nzz
			tmpi=indx(izz)
			tmpj=jndx(izz)
			tmpp=path(izz)

			IVZ=tmpj/(NVXS-1)+1
			IVX=MOD(tmpj,NVXS-1)
			IF(IVX.EQ.0) THEN
				IVX=NVXS-1
				IVZ=tmpj/(NVXS-1)
			END IF
			LIGHT(IVZ,IVX)=LIGHT(IVZ,IVX)+tmpp
		END DO

		DO IVX=1,NVXS
			WRITE(21,REC=IVX)(LIGHT(IVZ,IVX),IVZ=1,NVZS)
		END DO

	END SUBROUTINE ILLUMINATION

!**************=======================****************
!**************     END OF PROGRAM    ****************
!**************=======================****************

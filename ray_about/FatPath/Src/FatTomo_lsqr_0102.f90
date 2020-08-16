!/***************************************
!*     Author:luofei
!*     Date:2016-12-06 20:26
!*     Filename:RayTomo_lsqr.f90
!*     Description:
!*     胖层析反演,第一版,框架搭建
!*
!*     Last modified:2019-08-12 11:35
!****************************************/

!========================== MAIN ===============================!

PROGRAM TOMOGRAPHY_MAIN
	USE lsqrDataModule
	USE lsqrTestModule
	USE gli
	USE AnalysisGLI

IMPLICIT NONE
INCLUDE 'mpif.h'

INTEGER NX, NZ, NVX, NVZ, NVXS, NVZS,ivx,ivz,ielev,iline,ishot
INTEGER NS_X, NS_X2, NS_Z
REAL    DVX, DVZ
REAL	DXS,DZS
REAL    VX_START,VZ_START
REAL	XAPERTURE_LEFT, XAPERTURE_RIGHT,DEPTH
INTEGER SHOT_START,SHOTS,SHOT_INTERVAL,LINE_START,LINES,LINE_INTERVAL,IX
INTEGER ITR,NIT, IIT,tflag
INTEGER(ip) NTR, NV, NVS, NVSS
REAL	DSHOT,Total
REAL	DSTEP,t_period,i_sigma_rule

INTEGER ierr, myid, np,n

REAL,ALLOCATABLE::ELEV(:)
REAL,ALLOCATABLE::sss(:,:)
	
REAL,ALLOCATABLE::TOBS(:)
REAL,ALLOCATABLE::DV(:,:)
REAL,ALLOCATABLE::LIGHT(:,:)
REAL,ALLOCATABLE::TCAL(:)
REAL,ALLOCATABLE::J(:)
REAL,ALLOCATABLE::DSSS(:,:)
REAL,ALLOCATABLE::DVSS(:,:)
REAL(dp),ALLOCATABLE::TRES(:)	!b
REAL(dp),ALLOCATABLE::DSS(:)	!x
INTEGER,ALLOCATABLE::table(:,:,:)  !table(:,:,1) Every shot position in the geo file
								   !table(:,:,2) Each shot traces
								   !table(:,:,3) The start location of every shot in FirstBreak files

CHARACTER(LEN=256) FN1,FN2,FN3,FN4,FN5,FN6,FN7,FN8,FN9,FN10,fn_tmp
CHARACTER(LEN=256) CUR_IS,CUR_IIT,ORI_VEL,CUR_VEL,CUR_RAY_DAT,CUR_FAT_DAT
CHARACTER(LEN=256) CUR_LIGHT, CUR_LIGHT_IB
CHARACTER(LEN=256) DSS_FN,DV_FN,DVS_FN,DVSC_FN,DVSM_FN

LOGICAL::lexist             !True if fn_tmp exists

INTEGER::state              !=1,RayTracing+LSQR
							!=2,RayTracing
							!=3,LSQR

!*  MPI START
CALL MPI_INIT(ierr)
CALL MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)
CALL MPI_COMM_SIZE(MPI_COMM_WORLD, np, ierr)

CALL READPAR(FN1,FN2,FN3,FN4,FN5,FN6,FN7,FN8,FN9,FN10,NVX,NVZ,DVX,DVZ,& 
			DXS, DZS, VX_START, VZ_START, XAPERTURE_LEFT, XAPERTURE_RIGHT, DEPTH,&
			SHOT_START,SHOTS,SHOT_INTERVAL,LINE_START,LINES,LINE_INTERVAL,tflag,&
			NIT,DSTEP,t_period,i_sigma_rule,myid,state)

ALLOCATE(table(lines,shots,3))

fn_tmp=TRIM(ADJUSTL(fn2))//'_tmp.dat'

INQUIRE(FILE=fn_tmp,EXIST=lexist)

IF((myid==0).AND.(.NOT.lexist))THEN
	CALL InvGLI(fn2,fn_tmp)             !Convert GLI into Internal File
END IF

IF(myid==0)THEN
	CALL ScanningGLI(line_start,lines,line_interval,&
					shot_start,shots,shot_interval,&
					fn_tmp,table,ntr)
END IF

CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)

nvxs=ceiling(nvx/(dxs/dvx))
nvzs=ceiling(nvz/(dzs/dvz))

nv  = nvx*nvz
nvs = nvxs*nvzs

IF(myid.eq.0) THEN
	PRINT*,'DXS=',DXS,'DZS=',DZS
	PRINT*,'NX=',NX,'NZ',NZ
	PRINT*,'NV=',NV,'NVS=',NVS
	PRINT*,'NTR=',NTR
END IF
		
ALLOCATE(ELEV(NVX))
ALLOCATE(sss(NVZ,NVX))

open(250,file=FN10,access="direct",status="old",recl=nvx)
	read(250,rec=1) elev(1:nvx)
close(250)

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
	TOBS	=0.0
	TCAL	=0.0
	TRES	=0.0
	dv  	=0.0
	dsss	=0.0
	dvss	=0.0
	light	=0.0

	CALL Read_Time_obs(fn_tmp,ntr,tobs,lines,shots,table)

	CALL allocatememory()

END IF
		
CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)

	
OPEN(51,FILE=FN9,ACCESS='stream',status="replace")
DO WHILE(IIT.LE.NIT)
				
	IF(myid.eq.0) PRINT*,'IIT=',IIT

	!*	READ Velocity Field For Raytracing
	IF(IIT.EQ.1) THEN
		ORI_VEL=FN1
	ELSE
		ORI_VEL=CUR_VEL
	END IF

	IF(myid.eq.0)	PRINT*,'ORI_VEL=',ORI_VEL

	OPEN(11,FILE=ORI_VEL,ACCESS='DIRECT',STATUS='OLD',RECL=NVZ)
			
		CALL READ_INITIAL_SLOWNESS(NVX, NVZ, SSS)
			
	CLOSE(11)
			
	CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
	

	IF(state.ne.3)then
		OPEN(13,FILE=fn_tmp,STATUS='OLD',ACCESS='DIRECT',RECL=9)

			CALL TRAVELTIME_BASED_RAYTRACING(lines,shots,table,DSTEP,t_period,i_sigma_rule,tflag,&
						sss,ELEV,NVX,DVX,DVZ,NVZ,VX_START,VZ_START,XAPERTURE_LEFT, XAPERTURE_RIGHT, DEPTH,&
						NVXS,DXS,DZS,FN4,FN5,FN6,myid,np)

			CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
		CLOSE(13)
	END IF
	CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)

	IF(state.ne.2) THEN
		WRITE(CUR_IIT,'(I4)') IIT
		CUR_VEL=trim(adjustl(FN7))//'updatevel'//trim(adjustl(CUR_IIT))//'.dat'
		CUR_LIGHT=trim(adjustl(FN8))//'light'//trim(adjustl(CUR_IIT))//'.dat'
		CUR_LIGHT_IB=trim(adjustl(FN8))//'light_ib'//trim(adjustl(CUR_IIT))//'.dat'


		IF(myid.eq.0) THEN			!Master Process solve equation

			DSS_FN=trim(adjustl(FN8))//'dss'//trim(adjustl(CUR_IIT))//'.dat'
			DVS_FN=trim(adjustl(FN8))//'dvs'//trim(adjustl(CUR_IIT))//'.dat'
			DV_FN=trim(adjustl(FN8))//'dv'//trim(adjustl(CUR_IIT))//'.dat'
			DVSM_FN=trim(adjustl(FN8))//'dvsm'//trim(adjustl(CUR_IIT))//'.dat'
			DVSC_FN=trim(adjustl(FN8))//'dvsc'//trim(adjustl(CUR_IIT))//'.dat'

			DSS(1:NVS)=0.0
			Total     =0.0

			PRINT*,'Start to read ray path'

			CALL READ_RAY_DATA(TRES,TCAL,TOBS,NTR,lines,shots,NVSS,&
								NVXS,NVZS,FN5,FN6,Total)

			! OPEN(21,FILE=CUR_LIGHT,ACCESS='DIRECT',RECL=NVZS,STATUS='REPLACE')
			! OPEN(31,FILE=CUR_LIGHT_IB,ACCESS='DIRECT',RECL=NVZS,STATUS='REPLACE')

			! 	PRINT*,'Start to blanace illumination'

			! 	CALL BALANCE_ILLUMINATION(NVXS,NVZS,NTR,NVSS,TRES,LIGHT)
			
			! CLOSE(21)
			! CLOSE(31)

			PRINT*,'Start to solve equation'

			CALL solve(NTR,NVSS,DSS,TRES)		!Solve equation by lsqr
			
			J(IIT)=Total
			J(IIT)=SQRT(J(IIT)/J(1))
			WRITE(51) J(IIT)
			TRES     =0.0
			TCAL     =0.0
			Total    =0.0

			PRINT*,'Start to update velocity'

			OPEN(41,FILE=DSS_FN,ACCESS='DIRECT',RECL=NVZS)
			OPEN(42,FILE=DVS_FN,ACCESS='DIRECT',RECL=NVZS)
			OPEN(43,FILE=DV_FN,ACCESS='DIRECT',RECL=NVZ)
			OPEN(44,FILE=DVSM_FN,ACCESS='DIRECT',RECL=NVZ)
			OPEN(45,FILE=DVSC_FN,ACCESS='DIRECT',RECL=NVZS)
			
				CALL UPDATE_SLOWNESS(NV,NVX,NVZ,NVS,NVXS,NVZS,DSS,DSSS,DVSS,SSS,DV)

			CLOSE(41)
			CLOSE(42)
			CLOSE(43)
			CLOSE(44)
			CLOSE(45)

			print*,'Start to write Updated velocity'

			do ivx=1,nvx
				ielev = nint(elev(ivx)/dvz)
				do ivz=1,ielev
					sss(ivz,ivx)=2000.0
				end do
			end do

			PRINT*,'CUR_VEL=',CUR_VEL	
			OPEN(16,FILE=CUR_VEL,ACCESS='DIRECT',RECL=NVZ,STATUS='REPLACE')

				CALL WRITE_VELOCITY(NVX,NVZ,SSS)

			CLOSE(16)
			
			
		END IF
	END IF

	IIT=IIT+1

	CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
	
END DO
CLOSE(51)

DEALLOCATE(ELEV)	
DEALLOCATE(SSS)
deallocate(table)
		
IF(myid.eq.0) THEN

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
			
END PROGRAM TOMOGRAPHY_MAIN

!==================================================================
	SUBROUTINE TRAVELTIME_BASED_RAYTRACING(lines,shots,table,DSTEP,t_period,i_sigma_rule,tflag,&
					sss,ELEV,NVX,DVX,DVZ,NVZ,VX_START, VZ_START,XAPERTURE_LEFT, XAPERTURE_RIGHT, DEPTH,&
					NVXS,DXS,DZS,FN4,FN5,FN6,myid,np)

		USE global
		USE gli
		use FatPath

		implicit none

		INCLUDE 'mpif.h'

		INTEGER::shots,ntraces,trace_start,trace_locate,lines,iline,ishot
		INTEGER IV
		INTEGER NVX,NVZ,NVXS
		REAL	DXS,DZS,VX_START,VZ_START
		REAL,INTENT(IN)::XAPERTURE_LEFT, XAPERTURE_RIGHT, DEPTH
		INTEGER,DIMENSION(lines,shots,3)::table

		INTEGER myid,np,ntask,tflag
		CHARACTER*256 ProcessorName
		INTEGER  status(MPI_STATUS_SIZE)
		INTEGER  isend, itask, inode, ierr
		INTEGER  R_table_Max, S_table_Max
		integer  R_table(9), S_table(9)

        REAL   DVX, DVZ

		REAL::sx_coord,sz_coord
		REAL	DSTEP,t_period,i_sigma_rule
		CHARACTER(LEN=256) FN4,FN5,FN6
		CHARACTER(LEN=256) CUR_RAY,CUR_PATH,CUR_TIME,cur_iline,cur_ishot

		real,allocatable::time(:,:)
		REAL,DIMENSION(NVX):: ELEV
		REAL,DIMENSION(NVZ,NVX)::sss

		integer::ns_x,ns_z,nxs_left,nxs_right,i,nx,nz
		real::sx_min,sx_max,sx_left,sx_right
		real::dx,dz


		R_table_Max = 9
		S_table_Max = 9
		ntask=shots*lines

		IF (myid .eq. 0) THEN
!Master Node.
			write(*,*) 'Task number: ', ntask, 'Processor number: ', np
			do iline = 1,lines
				do ishot = 1,shots
					itask = ishot + (iline-1)*shots

					READ(13,REC=table(iline,ishot,1)) geo

					sx_coord=geo%x
					sz_coord=geo%z
					ntraces=table(iline,ishot,2)
					trace_start=table(iline,ishot,1)+1
					trace_locate=table(iline,ishot,3)

					CALL MPI_RECV(R_table, R_table_Max, MPI_INTEGER,&
									MPI_ANY_SOURCE, MPI_ANY_TAG,&
									MPI_COMM_WORLD, status, ierr)

					isend      = status(MPI_SOURCE)
					s_table(1) = 1
					s_table(2) = 0
					s_table(3) = iline
					s_table(4) = sx_coord
					s_table(5) = ishot
					s_table(6) = sz_coord
					s_table(7) = ntraces
					s_table(8) = trace_start
					s_table(9) = trace_locate

!				print*,itask,sx_coord,sz_coord,ntraces,trace_start,trace_locate

					CALL MPI_Send(S_table, S_table_Max, MPI_integer,&
									isend, itask,&
									MPI_COMM_WORLD, ierr)

				end do
			end do

			DO 2222 inode = 1, np-1

				CALL MPI_RECV(R_table, R_table_Max, MPI_integer,&
								MPI_ANY_SOURCE, MPI_ANY_TAG,&
								MPI_COMM_WORLD, status, ierr)

				isend      = status(MPI_SOURCE)
    			s_table(1) = 0
    			s_table(2) = 0
    			s_table(3) = 0
    			s_table(4) = 0
				s_table(5) = 0
				s_table(6) = 0
				s_table(7) = 0
				s_table(8) = 0
				s_table(9) = 0

	!			write(*,*) 'inode = ', inode, 'isend = ', isend
				CALL MPI_Send(S_table, S_table_Max, MPI_INTEGER,&
								isend, itask,&
								MPI_COMM_WORLD, ierr)

2222     	CONTINUE

			write(*,*) 'Done. myid is ', myid

		ELSE
!Working Node.

    		s_table(1) = 0
    		s_table(2) = 0
    		s_table(3) = 0
    		s_table(4) = 0
			s_table(5) = 0
			s_table(6) = 0
			s_table(7) = 0
			s_table(8) = 0
			s_table(9) = 0
			CALL MPI_Send(S_table, S_table_Max, MPI_INTEGER,&
							0, 0, MPI_COMM_WORLD, ierr)

			DO WHILE(.True.)

				CALL MPI_RECV(R_table, R_table_Max, MPI_INTEGER,&
								0, MPI_ANY_TAG,&
								MPI_COMM_WORLD, status, ierr)

				IF (R_table(3) .ne. 0) THEN

					iline        = r_table(3)
					sx_coord     = r_table(4)
					ishot		 = r_table(5)
					sz_coord     = r_table(6)
					ntraces      = r_table(7)
					trace_start  = r_table(8)
					trace_locate = r_table(9)

					write(*,*) 'Working node: ', myid, 'received itask: ', ishot+(iline-1)*shots

					!确定震源的位置
					NS_X=INT((SX_COORD-VX_START)/DVX)+1

					IF(NS_X.LT.1) NS_X=1
					IF(NS_X.GT.NVX) NS_X=NVX

					NS_Z =INT(ABS(SZ_COORD)/DVZ) + 1
					!*********************************

					!计算旅行时场孔径
					sx_min = 999999
					sx_max = 0
					do i=1,ntraces
						READ(13,REC=i+trace_start-1) geo
						if(geo%x<sx_min) sx_min = geo%x
						if(geo%x>sx_max) sx_max = geo%x
					end do

					sx_left = sx_coord - sx_min * 1.0

					IF(sx_left<=0.0) THEN
						sx_left = 100.0
					else if(sx_left>=xaperture_left) then
						sx_left = xaperture_left+100.0
					else
						sx_left = sx_left+100.0
					END IF

					sx_right = sx_max * 1.0 - sx_coord

					IF(sx_right<=0.0) THEN
						sx_right = 100.0
					else if(sx_right>=xaperture_right) then
						sx_right = xaperture_right+100.0
					else
						sx_right = sx_right+100.0
					END IF

					nx = int((sx_left+sx_right)/dvx) + 1
					nz = int(depth/dvz) + 1

					nxs_left = int(sx_left/dvx)
					nxs_right= int(sx_right/dvx)

					dx = dvx
					dz = dvz

					ALLOCATE(TIME(NX,NZ))

					!*************************************************

					WRITE(CUR_iline,'(I4)') iline
					WRITE(CUR_ishot,'(I4)') ISHOT
					CUR_RAY     =trim(adjustl(FN4))//'ray'//trim(adjustl(CUR_iline))//"_"//trim(adjustl(CUR_ishot))//'.txt'
					CUR_PATH    =trim(adjustl(FN5))//'path'//trim(adjustl(CUR_iline))//"_"//trim(adjustl(CUR_ishot))//'.dat'
					CUR_TIME    =trim(adjustl(FN6))//'traveltime'//trim(adjustl(cur_iline))//"_"//trim(adjustl(CUR_ishot))//'.dat'



					OPEN(14,FILE=CUR_RAY,STATUS='REPLACE')
					OPEN(17,FILE=CUR_PATH,ACCESS='STREAM',STATUS='REPLACE')
					OPEN(18,FILE=CUR_TIME,ACCESS='STREAM',STATUS='REPLACE')


                   	CALL TRAVELTIME_2D(TIME,sss,ELEV,NVX,NX,DVX,DVZ,DX,NVZ,&
							NZ,DZ,NS_X,NS_Z,SX_COORD,VX_START,NXS_LEFT,&
							NXS_RIGHT,ntraces,trace_start,trace_locate,&
							SZ_COORD,itask,DSTEP,NVXS,DXS,DZS,myid)


					call raytracing_fresnell(NX, NZ, NVX, NVZ, DVX, DVZ, DX, DZ, TIME, &
									ELEV,sss,SX_COORD,ntraces,trace_start,trace_locate,tflag,&
									VX_START, SZ_COORD, ISHOT, DSTEP, NVXS, DXS,NXS_LEFT, NXS_RIGHT,&
									DZS,t_period,i_sigma_rule,myid)

					CLOSE(14)
					CLOSE(17)
					CLOSE(18)

					deallocate(TIME)
					
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
	SUBROUTINE TRAVELTIME_2D(TIME,sss,ELEV,NVX,NX,DVX,DVZ,DX,NVZ,&
					NZ,DZ,NS_X,NS_Z,SX_COORD,VX_START,NXS_LEFT,&
					NXS_RIGHT,ntraces,trace_start,trace_locate,&
					SZ_COORD,ISHOT,DSTEP,NVXS,DXS,DZS,myid)
					
		USE global

		implicit none

		INTEGER NX, NZ, NS_X, NS_Z, NVX, NVZ, NVXS
		INTEGER NXS_LEFT, NXS_RIGHT
		INTEGER ISHOT,myid
        REAL    DX, DZ, DVX, DVZ
		REAL	SX_COORD, VX_START
		REAL	DSTEP,DXS,DZS
		REAL	DS_X, DS_Z
		REAL	SZ_COORD
		REAL	R1,R2
		integer ntraces,trace_start,trace_locate

		REAL,ALLOCATABLE::SS1(:)
		REAL,ALLOCATABLE::SS2(:)
		REAL,ALLOCATABLE::TT1(:)
		REAL,ALLOCATABLE::TT2(:)
		REAL,ALLOCATABLE::SLOWNESS(:,:)
		REAL,ALLOCATABLE::SLOW_45(:,:)
		REAL,ALLOCATABLE::EPSILON_VALUE(:,:)

		real,dimension(NX,NZ),intent(inout)::time
		REAL,DIMENSION(NVX),intent(in):: ELEV
		REAL,DIMENSION(NVZ,NVX),intent(in)::sss
		
		ALLOCATE(SS1(NX+1))
		ALLOCATE(SS2(NX+1))
		ALLOCATE(TT1(0:(NX+1)))
		ALLOCATE(TT2(0:(NX+1)))
		ALLOCATE(SLOWNESS(NX+1,NZ+1))
		ALLOCATE(SLOW_45(NX+1,NZ+1))
		ALLOCATE(EPSILON_VALUE(NX+1,NZ+1))


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


		DEALLOCATE(SS1)
		DEALLOCATE(SS2)
		DEALLOCATE(TT1)
		DEALLOCATE(TT2)
		DEALLOCATE(SLOWNESS)
		DEALLOCATE(SLOW_45)
		DEALLOCATE(EPSILON_VALUE)
	
		RETURN
    END 

!=============================================================================

    SUBROUTINE TRAVELTIME_CAL(SS1, SS2, TT1, TT2, SLOWNESS, TIME,&
				SLOW_45, EPSILON_VALUE, NX, DX, NY, DY, &
                NS_X, NS_Y, DS_X, DS_Y)
 
		INTEGER NX, NY, NS_X, NS_Y
		REAL    DX, DY, DS_X, DS_Y
       
		DIMENSION SS1(0:NX), SS2(0:NX), TT1(0:(NX+1)), TT2(0:(NX+1))
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
		DIMENSION SS1(0:NX), TT1(0:(NX+1))
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
		DIMENSION  SS1(0:NX), SS2(0:NX), TT1(0:(NX+1)), TT2(0:(NX+1))
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
		DIMENSION SS1(0:NX), SS2(0:NX), TT1(0:(NX+1)), TT2(0:(NX+1))
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
		DIMENSION SS1(0:NX), SS2(0:NX), TT1(0:(NX+1)), TT2(0:(NX+1))
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
		DIMENSION SS1(0:NX), SS2(0:NX), TT1(0:(NX+1)), TT2(0:(NX+1))
		DIMENSION TIME(NX, NZ), SLOWNESS(0:NX, 0:NZ)
		DIMENSION SLOW_45(0:NX, 0:NZ)

!=======================================================================
      
!		WRITE(*,*)  'NX, NZ=', NX, NZ
      
!=======================================================================      
		DO IX=0, NX
			SS1(IX)=0.0
			SS2(IX)=0.0
		END DO

		DO IX=0, NX+1
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
						DVXS, DVZS, VX_START, VZ_START,	XAPERTURE_LEFT, XAPERTURE_RIGHT, DEPTH,&
						SHOT_START,SHOTS,SHOT_INTERVAL,LINE_START,LINES,LINE_INTERVAL,tflag,&
						NIT,DSTEP,t_period,i_sigma_rule,myid,state)

		INTEGER NVX,NVZ
		INTEGER SHOT_START,SHOTS,SHOT_INTERVAL,LINE_START,LINES,LINE_INTERVAL
		INTEGER myid,NIT,state,tflag
		REAL	t_period,i_sigma_rule,DSTEP
		REAL	DVX,DVZ,DVXS,DVYS
		REAL    XAPERTURE_LEFT,XAPERTURE_RIGHT,DEPTH
		REAL    VX_START,VZ_START
		CHARACTER(256) FN1,FN2,FN3,FN4,FN5,FN6,FN7,FN8,FN9,FN10,PAR_JUMPPER

		OPEN(11, FILE='InvTomo.par')

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
		READ(11,*) DVXS,DVZS
		
		READ(11,'(A)') PAR_JUMPPER
		READ(11,*) VX_START,VZ_START

        READ(11,'(A)') PAR_JUMPPER
		READ(11,*) XAPERTURE_LEFT,XAPERTURE_RIGHT, DEPTH

		READ(11,'(A)') PAR_JUMPPER
		READ(11,*) LINE_START,LINES,LINE_INTERVAL

		READ(11,'(A)') PAR_JUMPPER
		READ(11,*) SHOT_START,SHOTS,SHOT_INTERVAL

		READ(11,'(A)') PAR_JUMPPER
		READ(11,*) NIT

		READ(11,'(A)') PAR_JUMPPER
		READ(11,*) state

		READ(11,'(A)') PAR_JUMPPER
		READ(11,*) tflag

		READ(11,'(A)') PAR_JUMPPER
		READ(11,*) DSTEP

		READ(11,'(A)') PAR_JUMPPER
		READ(11,*) t_period

		READ(11,'(A)') PAR_JUMPPER
		READ(11,*) i_sigma_rule

		IF(myid.eq.0) THEN

			WRITE(*,*) 'THE INITIAL VELOCITY'
			WRITE(*,'(A)') FN1
			
			WRITE(*,*) 'GLI'
			WRITE(*,'(A)') FN2

			WRITE(*,*) 'THE TRAVELTIME FILENAME'
			WRITE(*,'(A)') FN3

			WRITE(*,*) 'THE RAYTRACING RESULT(.TXT)'
			WRITE(*,'(A)') FN4
		        
			WRITE(*,*) 'THE Ray PATH FILENAME(.Dat)'
			WRITE(*,'(A)') FN5

			WRITE(*,*) 'THE CALCULATE TRACELTIME'
			WRITE(*,'(A)') FN6

			WRITE(*,*) 'THE UPDATED VELOCITY'
			WRITE(*,'(A)') FN7
			
			WRITE(*,*) 'THE LSQR PARAMETER(DSS,DVS,DV,DVSM,DVSC)'
			WRITE(*,'(A)') FN8

			WRITE(*,*) 'THE CONVERGENCE CURVE FILENAME'
			WRITE(*,'(A)') FN9
		        
			WRITE(*,*) 'THE Elevation file'
			WRITE(*,'(A)') FN10

			WRITE(*,*) 'NX,NZ OF THE VELOCITY MODEL'
			WRITE(*,*) NVX,NVZ

			WRITE(*,*) 'DX,DZ OF THE VELOCITY MODEL'
			WRITE(*,*) DVX,DVZ

			WRITE(*,*) 'DVXS,DVZS WITH SOLVING EQUATION'
			WRITE(*,*) DVXS,DVZS

			WRITE(*,*) 'VELOCITY START COORDINATE IN X/Z DIMENSION'
			WRITE(*,*) VX_START,VZ_START

			WRITE(*,*) 'XAPERTURE_LEFT,XAPERTURE_RIGHT AND DEPTH'
			WRITE(*,*) XAPERTURE_LEFT,XAPERTURE_RIGHT,DEPTH

			WRITE(*,*) 'line_start,lines,line_interval'
			WRITE(*,*) LINE_START,LINES,LINE_INTERVAL

			WRITE(*,*) 'shot_start,shots,shot_interval'
			WRITE(*,*) SHOT_START,SHOTS,SHOT_INTERVAL

			WRITE(*,*) 'INTERATIONS'
			WRITE(*,*) NIT

			WRITE(*,*) '=1,RayTracing+LSQR;=2,RayTracing;=3,LSQR'
			WRITE(*,*) state

			WRITE(*,*) 'Calculate fresnell zone(1=Dynamic;2=Approximate)'
			WRITE(*,*) tflag

			WRITE(*,*) 'DSTEP'
			WRITE(*,*) DSTEP

			WRITE(*,*) 'Period (1/HZ)'
			WRITE(*,*) t_period

			WRITE(*,*) 'Sigma'
			WRITE(*,*) i_sigma_rule

		END IF
		
		CLOSE(11)

	END

!***********************************************************************!	
	SUBROUTINE READ_INITIAL_SLOWNESS(NVX, NVZ, SSS)
		INTEGER NVX, NVZ
		INTEGER IVX, IVZ
		REAL	SSS(NVZ,NVX)
		
		DO IVX=1,NVX
			READ(11,REC=IVX)(sss(IVZ,IVX),IVZ=1,NVZ)
		END DO

		DO IVX=1,NVX
			DO IVZ=1,NVZ
				SSS(IVZ,IVX)=1.0/sss(IVZ,IVX)
			END DO
		END DO
			
	END

!***********************************************************************!
	SUBROUTINE UPDATE_SLOWNESS(NV,NVX,NVZ,NVS,NVXS,NVZS,DSS,DSSS,DVSS,SSS,DV)
		USE lsqrDataModule
		USE lsqrTestModule
		USE global
		use LinearInterTravelTime

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

		DXS=ceiling(nvx*1.0/nvxs)
		DZS=ceiling(nvz*1.0/nvzs)

		!*	1D slowness increment to 2D
		DO IVZ=1,NVZS
			DO IVX=1,NVXS
				IV=(IVZ-1)*NVXS+IVX	
				DSSS(IVZ,IVX)=DSS(IV)	
			END DO
		END DO

		DO IVX=1,NVXS
			WRITE(41,REC=IVX)(DSSS(IVZ,IVX),IVZ=1,NVZS)
		END DO

		!*	sampling original slowness to NVZS*NVXS
		DO IVX=1,NVXS
			IIVX=(IVX-1)*DXS+1
			DO IVZ=1,NVZS
				IIVZ=(IVZ-1)*DZS+1
				DVSS(IVZ,IVX)=SSS(IIVZ,IIVX)
			END DO
		END DO

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
		
		!*	Constraint velocity increment 
		DO IVX=1,NVXS
			DO IVZ=1,NVZS
				IF(ABS(DVSS(IVZ,IVX)).GT.3000.0) THEN
					DVSS(IVZ,IVX)=0.0
				END IF
			END DO
		END DO
		DO IVX=1,NVXS
			WRITE(45,REC=IVX)(DVSS(IVZ,IVX),IVZ=1,NVZS)
		END DO

		!*	smooth velocity increment

		R1=10.0
		R2=10.0
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
	SUBROUTINE READ_RAY_DATA(TRES,TCAL,TOBS,NTR,lines,shots,NVSS,&
								NVXS,NVZS,FN5,FN6,Total)

		USE lsqrDataModule
		USE lsqrTestModule

		REAL(dp)	TRES(NTR)
		REAL		TOBS(NTR)
		REAL		TCAL(NTR)
		INTEGER(ip)	NTR, NV, IV, NVSS
		INTEGER		NVXS,NVZS,IVX,IVZ,iline,ishot
		INTEGER		ITR,NSHOT,lines,shots
		INTEGER tmpi,tmpj
		REAL	tmpp,tmpt,tmpt1,Total
		INTEGER j
		
		CHARACTER(LEN=256) FN5,FN6,CUR_RAY,CUR_PATH,CUR_TIME,cur_line,cur_shot


		!*	initialize array
		CALL initialization()

		NVSS=0
		ITR=0

		do iline = 1,lines
			do ishot = 1,shots
				WRITE(cur_line,'(I4)') iline
				WRITE(cur_shot,'(I4)') ishot
				CUR_TIME=trim(adjustl(FN6))//'traveltime'//trim(adjustl(cur_line))//"_"//trim(adjustl(cur_shot))//'.dat'
				OPEN(88,FILE=CUR_TIME,ACCESS='STREAM',STATUS='OLD')
				DO while(.true.)
					READ(88,iostat=j) ITR,tmpt1,tmpt
					IF(J/=0) EXIT
					TCAL(ITR)=tmpt
					Total = Total + (TOBS(ITR)-tmpt1)*(TOBS(ITR)-tmpt1)
				END DO
				CLOSE(88)
			END DO
		end do

		do iline = 1,lines
			do ishot = 1,shots
				WRITE(cur_line,'(I4)') iline
				WRITE(cur_shot,'(I4)') ishot
				CUR_PATH=trim(adjustl(FN5))//'path'//trim(adjustl(cur_line))//"_"//trim(adjustl(cur_shot))//'.dat'
				OPEN(77,FILE=CUR_PATH,ACCESS='STREAM',STATUS='OLD')
				DO
					READ(77,iostat=j) tmpi,tmpj,tmpp
					if(j/=0) exit
					nzz=nzz+1
					IF(nzz.GT.maxnz) THEN
						PRINT*,'nnz is bigger than maxnz'
						STOP	
					END IF
					indx(nzz)=tmpi
					jndx(nzz)=tmpj
					path(nzz)=tmpp
					IF(NVSS.LT.tmpj) NVSS=tmpj
				END DO
				CLOSE(77)
			END DO
		end do
		
		DO ITR=1,NTR
			IF(TCAL(ITR).LE.0.0)THEN
				TCAL(ITR)=0.0
			ELSE
				TCAL(ITR)=TOBS(ITR)-TCAL(ITR)
			END IF
			TRES(ITR)=TCAL(ITR)
		END DO
		
		RETURN
	END

!***********************************************************************!		
!*************************************************************************!
!                     Read  FirstBreak by Picking                         !
!*************************************************************************!
SUBROUTINE Read_Time_obs(fn_tmp,ntr,tobs,lines,shots,table)
	USE global
	USE gli

IMPLICIT NONE

INTEGER::nu,l,s,itr,i
INTEGER,INTENT(IN)::ntr,lines,shots,table(lines,shots,3)
REAL,INTENT(INOUT)::tobs(ntr)
CHARACTER(LEN=256),INTENT(IN)::fn_tmp

itr=1
INQUIRE(IOLENGTH=nu) geo
OPEN(55,FILE=fn_tmp,ACCESS='DIRECT',STATUS='OLD',RECL=nu*length)
DO l=1,lines
	DO s=1,shots
		DO i=1,table(l,s,2)
			READ(55,REC=i+table(l,s,1)) geo
			tobs(itr)=geo%fb*0.001
			itr=itr+1
		END DO
	END DO
END DO
CLOSE(55)
END SUBROUTINE Read_Time_obs
!***********************************************************************!		
	SUBROUTINE BALANCE_ILLUMINATION(NVXS,NVZS,NTR,NVSS,TRES,LIGHT)
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

			IVZ=tmpj/(NVXS)+1
			IVX=MOD(tmpj,NVXS)
			IF(IVX.EQ.0) THEN
				IVX=NVXS
				IVZ=tmpj/(NVXS)
			END IF
			LIGHT(IVZ,IVX)=LIGHT(IVZ,IVX)+tmpp
		END DO

		DO IVX=1,NVXS
			WRITE(21,REC=IVX)(LIGHT(IVZ,IVX),IVZ=1,NVZS)
		END DO

!		CALL illumination(NVXS,NVZS,nzz,NTR,NVSS,indx,jndx,path,TRES)

		LIGHT(1:NVZS,1:NVXS)=0.0

		DO izz=1,nzz
			tmpi=indx(izz)
			tmpj=jndx(izz)
			tmpp=path(izz)

			IVZ=tmpj/(NVXS)+1
			IVX=MOD(tmpj,NVXS)
			IF(IVX.EQ.0) THEN
				IVX=NVXS
				IVZ=tmpj/(NVXS)
			END IF
			LIGHT(IVZ,IVX)=LIGHT(IVZ,IVX)+tmpp
		END DO

		DO IVX=1,NVXS
			WRITE(31,REC=IVX)(LIGHT(IVZ,IVX),IVZ=1,NVZS)
		END DO

	END SUBROUTINE BALANCE_ILLUMINATION
!===========================================================================
!===========================================================================
!**************=======================****************
!**************     END OF PROGRAM    ****************
!**************=======================****************


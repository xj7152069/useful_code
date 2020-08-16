module FatPath
implicit none
public::raytracing_fresnell

contains

SUBROUTINE raytracing_fresnell(NX, NZ, NVX, NVZ, DVX, DVZ, DX, DZ, TIME, &
								ELEV,sss,SX_COORD,ntraces,trace_start,trace_locate,tflag,&
								VX_START, SZ_COORD, ISHOT, DSTEP, NVXS, DXS,NXS_LEFT, NXS_RIGHT,&
								DZS,t_period,i_sigma_rule,myid)

USE global
USE gli

INTEGER NX, NZ, NVX, NVXS,NVZ,ivx,ivz
INTEGER::i,myid,tflag
INTEGER NR_X,NXS_LEFT, NXS_RIGHT
INTEGER ISHOT
REAL	DX, DZ, DVX, DVZ,SX_COORD
REAL	VX_START, XMIN, ZMIN
REAL	RX_COORD, RZ_COORD
REAL	DSTEP,DXS,DZS,t_period,i_sigma_rule
REAL	SZ_COORD
real,dimension(nx,nz)::time
real,dimension(nvz,nvx)::sss
real,dimension(nvx)::ELEV
INTEGER	ITR
INTEGER,INTENT(IN)::ntraces,trace_start,trace_locate
real,allocatable::vel(:,:)


XMIN=SX_COORD-NXS_LEFT*DX
ZMIN=0.0

allocate(vel(nvz,nvx))
do ivx = 1,nvx
	do ivz = 1,nvz
		vel(ivz,ivx) = 1.0/sss(ivz,ivx)
	end do
end do

DO i=1,ntraces
	READ(13,REC=i+trace_start-1) geo
	
	RX_COORD=geo%x
	NR_X=NINT((RX_COORD-VX_START)/DVX)+1

	RZ_COORD=geo%z

    ITR=trace_locate+i-1

	IF(NR_X.LT.1) THEN
		write(18) itr,0.0,-1.0
		CYCLE
	END IF
	IF(NR_X.GT.NVX) THEN
		write(18) itr,0.0,-1.0
		CYCLE
	END IF


	CALL TRACING_rolling_fresnel_Modify(NX, NZ, NVX, NVZ, DVX, DVZ, DX,	DZ, TIME, ELEV,XMIN, ZMIN, DSTEP,&
							SX_COORD,SZ_COORD,RX_COORD,RZ_COORD,VX_START,VEL,tflag,&
							NVXS, DXS, DZS, t_period,i_sigma_rule,ITR,myid)

END DO

deallocate(vel)

END SUBROUTINE raytracing_fresnell
!***********************************************************************!
!===============================================================================!
SUBROUTINE TRACING_rolling_fresnel_Modify(NX, NZ, NVX, NVZ, DVX, DVZ, DX, DZ, TIME, ELEV,XMIN, ZMIN, DSTEP,&
					SX_COORD,SZ_COORD,RX_COORD,RZ_COORD,VX_START,VEL,tflag,&
					NVXS, DXS, DZS, t_period,i_sigma_rule,ITR,myid)
USE global
use LinearInterTravelTime

		
implicit none
INTEGER NX, NZ, NVXS, NVX,NVZ,nv

INTEGER ITR,NN,i,n_point,myid,j,tflag
REAL	DX, DZ, DVX, DVZ
REAL	XMIN, ZMIN
REAL	SX_COORD,SZ_COORD,RX_COORD,RZ_COORD
REAL	VX_START
REAL	TMP
REAL	DSTEP
REAL	DXS,DZS,t_period,i_sigma_rule
real,dimension(nx,nz)::time
real,dimension(nvz,nvx)::vel
real,dimension(nvx)::elev
real,dimension(4,Ntable)::tbl_v,tbl_vx

REAL	X, Z, T, DTDX, DTDZ,t1!CURRENT POINT PARAMETER
REAL	X1,Z1,X2,Z2,X3,Z3,X4,Z4
INTEGER IVX1,IVZ1,IVX2,IVZ2,IVX,IVZ,IV
REAL	PATH1,PATH2,PATH,TMP1,TMP2,TMP3
real 	m,n,a,b

REAL,ALLOCATABLE::x_loc_ray(:),z_loc_ray(:),xp_ray(:),zp_ray(:),tau_ray(:)
real,allocatable::kernel(:),ray(:)

NN=2*(nvx*NINT(dvx/MIN(dvx,dvz))+nvz*NINT(dvz/MIN(dvx,dvz)))
nv=nvx*nvz

ALLOCATE(x_loc_ray(NN))
ALLOCATE(z_loc_ray(NN))
ALLOCATE(xp_ray(NN))
ALLOCATE(zp_ray(NN))
ALLOCATE(tau_ray(NN))
ALLOCATE(kernel(nv))
ALLOCATE(ray(nv))

do j=1,nn
	x_loc_ray(j)=0.0
	z_loc_ray(j)=0.0
	xp_ray(j)   =0.0
	zp_ray(j)   =0.0
	tau_ray(j)  =0.0
end do

do j=1,nv
	kernel(j)=0.0
	ray(j)=0.0
end do

X=RX_COORD
Z=RZ_COORD

PATH=0.0
PATH1=0.0
PATH2=0.0
i=1
n_point=0


CALL VELINTERPED(X, Z, T, XMIN, ZMIN, DX, DZ,&
				DTDX, DTDZ, NX, NZ, TIME)

TMP=T
t1 =T

m=dtdx
n=dtdz
		
x_loc_ray(i)=X
z_loc_ray(i)=Z
xp_ray(i)   =-DTDX/SQRT(DTDX*DTDX+DTDZ*DTDZ)
zp_ray(i)   =-DTDZ/SQRT(DTDX*DTDX+DTDZ*DTDZ)
tau_ray(i)  =t1-T
i=i+1

DO WHILE(1)	
	X=X-DSTEP*DTDX/SQRT(DTDX*DTDX+DTDZ*DTDZ)
	Z=Z-DSTEP*DTDZ/SQRT(DTDX*DTDX+DTDZ*DTDZ)
				   
                                   				
	TMP3=SQRT((X-SX_COORD)*(X-SX_COORD)+&
			(Z-SZ_COORD)*(Z-SZ_COORD))
					

	IF(Z.LT.0.0) THEN
		tmp=-1
		EXIT
	else
		CALL VELINTERPED(X, Z, T, XMIN, ZMIN, DX, DZ,&
						DTDX, DTDZ, NX, NZ, TIME)

		IF(TMP3.LT.DSTEP) THEN
			x_loc_ray(i)=X
			z_loc_ray(i)=Z
			xp_ray(i)   =(x-sx_coord)/dstep;
			zp_ray(i)   =(z-sz_coord)/dstep;
			tau_ray(i)  =t1-T
			i=i+1
			! WRITE(14,*) sx_coord,sz_coord 
			EXIT
		END IF
		
		! WRITE(14,*) X, Z                            !raycoordinate


		b=(m*dtdx+n*dtdz)/((sqrt(m*m+n*n))*(sqrt(dtdx*dtdx+dtdz*dtdz)))-0.000001
		a=acos(b)

		if(.NOT.(0.0<=a.AND.a<=1.256))then
			tmp=-1
			exit
		end if

		m=dtdx
		n=dtdz

		x_loc_ray(i)=X
		z_loc_ray(i)=Z
		xp_ray(i)   =-DTDX/SQRT(DTDX*DTDX+DTDZ*DTDZ)
		zp_ray(i)   =-DTDZ/SQRT(DTDX*DTDX+DTDZ*DTDZ)
		tau_ray(i)  =t1-T
		i=i+1
	end if            
END DO

n_point=i-1
	
WRITE(18) ITR,t1,TMP


if(tmp.gt.0.0) then

	DO i=1,n_point
		WRITE(14,*) x_loc_ray(i), z_loc_ray(i)
	END DO
	WRITE(14,*) sx_coord,sz_coord

	call calculate_kernel_fresnel(x_loc_ray,z_loc_ray,xp_ray,zp_ray,tau_ray,tflag,&
			 n_point,nvx,nvz,dvx,dvz,t_period,i_sigma_rule,sx_coord,sz_coord,dstep,itr,&
			 elev,vel,ray,kernel,myid)

	do ivx=1,nv
		if(kernel(ivx).ne.0.0)then
			write(17) itr,ivx,kernel(ivx)
		end if
	end do
end if

DEALLOCATE(x_loc_ray)
DEALLOCATE(z_loc_ray)
DEALLOCATE(xp_ray)
DEALLOCATE(zp_ray)
DEALLOCATE(tau_ray)
DEALLOCATE(kernel)
DEALLOCATE(ray)

RETURN
END SUBROUTINE TRACING_rolling_fresnel_Modify

end module FatPath

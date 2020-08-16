module LinearInterTravelTime
implicit none
public::VELINTERPED
public::PRE_COMPUTED_INTERPOLATORS
public::PreInterpolators_Bspline

contains


SUBROUTINE VELINTERPED(x,z,v,Xmin,Zmin,dvx,dvz,&
						dvdx,dvdz,nvx,nvz,Vel)
		USE global
		implicit none

		INTEGER nvx,nvz
	    INTEGER ixxx,izzz
	    INTEGER lx,lz
	    INTEGER mx,mz
	    INTEGER nvxm,nvzm
	    REAL	ax,az
		REAL	a,b,c,d
		REAL    x,z,v,Xmin,Zmin,dvx,dvz
		REAL    dvdx,dvdz
		REAL,dimension(nvx,nvz)::Vel

		integer i,j,jx,jmx,jz,jmz

		real,dimension(4,Ntable)::tbl_v,tbl_vx

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

		call PRE_COMPUTED_INTERPOLATORS(tbl_v,tbl_vx)
		
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
SUBROUTINE	 PRE_COMPUTED_INTERPOLATORS(tbl_v,tbl_vx)
	USE global
	implicit none

	real,dimension(4,Ntable),intent(inout)::tbl_v,tbl_vx
	integer::i
	real::x

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
END SUBROUTINE PRE_COMPUTED_INTERPOLATORS
!***********************************************************************
!                Interpolation coefficient(Bspline)                    !  
!======================================================================*
Subroutine PreInterpolators_Bspline(tbl_v,tbl_vx)
      use global
      INTEGER:: i
	  REAL::x
	  real,dimension(4,Ntable),intent(inout)::tbl_v,tbl_vx
       
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

end module LinearInterTravelTime

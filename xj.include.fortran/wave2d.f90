!/*********(version 1.0)***********/
!*
!   fortran module file: 
!/
!/********************************/
module  wave2d_type

    type wave2d
        real(4)  :: xs2(0:4)=[1.666667,-0.238095,0.039683,-0.004960,0.000317]
        real(4)  :: xs1(0:9)=[-0.0007926,0.00991800,-0.0595200,0.238080,-0.833333,&
        0.833333,-0.238080,0.0595200,-0.00991800,0.0007926]
        real(4) :: dx,dz,dt,PML_wide,R
        integer :: nx,nz,suface
        real(4),allocatable :: sw(:,:)
        real(4),allocatable :: sx11(:,:)
        real(4),allocatable :: sx12(:,:)
        real(4),allocatable :: sx13(:,:)
        real(4),allocatable :: sxp21(:,:)
        real(4),allocatable :: sxp22(:,:)
        real(4),allocatable :: sxp23(:,:)
        real(4),allocatable :: sx21(:,:)
        real(4),allocatable :: sx22(:,:) 
        real(4),allocatable :: sx31(:,:)
        real(4),allocatable :: sx32(:,:)
        real(4),allocatable :: sx33(:,:)
        real(4),allocatable :: vel(:,:)
        real(4),allocatable :: s1(:,:)
        real(4),allocatable :: s2(:,:)
        real(4),allocatable :: s3(:,:)
    end	type

end module wave2d_type

module  wave2d_basic
    use wave2d_type
implicit none   

contains 
    subroutine wave2d_creater(w,nz,nx)
        integer :: nz,nx
        type(wave2d),intent(inout) :: w

        w%nz=nz
        w%nx=nx
        w%dx=5.0
        w%dz=5.0
        w%dt=0.0005
        w%PML_wide=30
        w%suface=1
        w%R=1000

        allocate(w%sw(0:nz, 0:nx))

        allocate(w%s1(0:nz, 0:nx))

        allocate(w%s2(0:nz, 0:nx))

        allocate(w%s3(0:nz, 0:nx))

        allocate(w%sx11(0:nz, 0:nx))

        allocate(w%sx12(0:nz, 0:nx))

        allocate(w%sx13(0:nz, 0:nx))
        allocate(w%sx21(0:nz, 0:nx))
        allocate(w%sx22(0:nz, 0:nx))

        allocate(w%sx31(0:nz, 0:nx))

        allocate(w%sx32(0:nz, 0:nx))

        allocate(w%sx33(0:nz, 0:nx))

        allocate(w%sxp21(0:nz, 0:nx))

        allocate(w%sxp22(0:nz, 0:nx))

        allocate(w%sxp23(0:nz, 0:nx))
        allocate(w%vel(0:nz, 0:nx))

        w%sw(:,:)=0.0
        w%s1(:,:)=0.0
        w%s2(:,:)=0.0
        w%s3(:,:)=0.0
        w%sx11(:,:)=0.0
        w%sx12(:,:)=0.0
        w%sx13(:,:)=0.0
        w%sxp21(:,:)=0.0
        w%sxp22(:,:)=0.0
        w%sxp23(:,:)=0.0
        w%sx21(:,:)=0.0
        w%sx22(:,:)=0.0
        w%sx31(:,:)=0.0
        w%sx32(:,:)=0.0
        w%sx33(:,:)=0.0
        w%vel(:,:)=3000.0
    end subroutine wave2d_creater

    subroutine wave2d_del(w)
        type(wave2d),intent (inout) :: w
        deallocate(w%sw)
        deallocate(w%s1)
        deallocate(w%s2)
        deallocate(w%s3)
        deallocate(w%sx11)
        deallocate(w%sx12)
        deallocate(w%sx13)
        deallocate(w%sx21)
        deallocate(w%sx22)
        deallocate(w%sxp21)
        deallocate(w%sxp22)
        deallocate(w%sxp23)
        deallocate(w%sx31)
        deallocate(w%sx32)
        deallocate(w%sx33)
        deallocate(w%vel)
    end subroutine wave2d_del

    subroutine wave2d_timeslicecal(w)
        type(wave2d),intent(inout) :: w
        real(4) :: DX,DY,DT,xshd
        integer :: X,Y,suface_PML
        real(4) :: fadx,fady,faddx,faddy,snx1,sny1,snx2,sny2,t2,t5
        integer :: i,j,n,t3,t4
        real(4) :: u1,u2,u,ux,uy
        real(4) :: C_X
        real(4) :: C_Y
        real(4) :: DT2,DT3,DX2,DY2,mo2
    
        DX=w%dx
        DY=w%dz
        DT=w%dt
        xshd=w%PML_wide
        X=w%nx
        Y=w%nz
        suface_PML=w%suface

        DT2=DT*DT
        DT3=DT2*DT
        DX2=DX*DX
        DY2=DY*DY
        t5=real(Y)/real(X)
        u1=0
        u2=0
        u=0
        ux=0
        uy=0
        C_X=log(w%R)*3.0/2.0/(xshd)/(xshd)/(xshd)/DX2/DX
        C_Y=log(w%R)*3.0/2.0/(xshd)/(xshd)/(xshd)/DY2/DY

        iloop: do i=5, X-5

            if (i>=0.5*(X)) then
                t3=1
            else
                t3=-1
            end if
    
            jloop: do j=5, Y-5
                
            !根据系数求得二阶偏微分的离散算子
                n0loop: do n=0, 4
                    u=u+2*w%xs2(n);
                    u1=u1+w%xs2(n)*(w%s2(j-n-1,i)+w%s2(j+n+1,i))
                    u2=u2+w%xs2(n)*(w%s2(j,i-n-1)+w%s2(j,i+n+1))
                end do n0loop
    
                snx1=0.0
                snx2=0.0
                sny1=0.0
                sny2=0.0
                fadx=0
                faddx=0
                fady=0
                faddy=0
    
                if (suface_PML==1) then
                    if (i>=X-xshd-5 .and. j<t5*i .and. j>-t5*i+Y) then
                        snx1=i-(X-xshd-5)
                    end if
                    if (i<=xshd+5 .and. j>t5*i .and. j<-t5*i+Y) then	
                        snx2=xshd+5-i
                    end if
                    if(j>=Y-xshd-5 .and. j>=t5*i .and. j>=-t5*i+Y) then
                        sny1=j-(Y-xshd-5)
                    end if
                    if(j<=xshd+5 .and. j<=t5*i .and. j<=-t5*i+Y) then
                        sny2=xshd+5-j
                    end if
                else
                    if (i>=X-xshd-5 .and. j<=t5*i) then
                        snx1=i-(X-xshd-5)
                    end if
                    if (i<=xshd+5 .and. j<=-t5*i+Y) then		
                        snx2=xshd+5-i
                    end if
                    if(j>=Y-xshd-5 .and. j>t5*i .and. j>-t5*i+Y) then		
                        sny1=j-(Y-xshd-5)
                    end if
                    if (j<=xshd+5) then		
                        sny2=0
                    end if
                end if
    
                if (sny1 /= 0) then
                    fady=w%vel(j,i)*C_Y*sny1*sny1*DY2
                    faddy=w%vel(j,i)*C_Y*2*sny1*DY
                end if
                if (sny2 /= 0) then
                    fady=w%vel(j,i)*C_Y*sny2*sny2*DY2
                    faddy=w%vel(j,i)*C_Y*2*sny2*DY
                end if
                if (snx1 /= 0) then
                    fadx=w%vel(j,i)*C_X*snx1*snx1*DX2
                    faddx=w%vel(j,i)*C_X*2*snx1*DX
                end if
                if (snx2 /= 0) then 
                    fadx=w%vel(j,i)*C_X*snx2*snx2*DX2
                    faddx=w%vel(j,i)*C_X*2*snx2*DX
                end if
    
                if(j>=0.5*(Y)) then
                    t4=1
                else
                    t4=-1
                end if
    
                mo2=w%vel(j,i)*w%vel(j,i)
    
                if(snx1/=0 .or. snx2/=0) then
                !根据系数求得一阶偏微分的离散算子
                    n1loop: do n=0, 9
                        if(n<5) then
                            ux=ux+w%s2(j,i+n-5)*w%xs1(n)
                            uy=uy+w%s2(j+n-5,i)*w%xs1(n)
                        else
                            ux=ux+w%s2(j,i+n-4)*w%xs1(n)
                            uy=uy+w%s2(j+n-4,i)*w%xs1(n)
                        end if
                    end do n1loop
    
                    !equation 1
                    w%sx11(j,i)=mo2*DT2*(u2-u*w%s2(j,i))*(1.0/(DX2)) &
                    -fadx*fadx*DT2*w%sx12(j,i)+(2*w%sx12(j,i) &
                    -w%sx13(j,i))+DT*(2*fadx*(w%sx13(j,i)-w%sx12(j,i)))
    
                    !equation 2
                    w%sxp21(j,i) = 2.0*w%sxp22(j,i) - w%sxp23(j,i) &
                    +DT2*(-1.0*mo2*faddx*(1.0/(DX))*(ux*t3) - &
                    2.0*fadx*(w%sxp22(j,i) - w%sxp23(j,i))/DT - fadx*fadx*w%sxp22(j,i))
                    w%sx21(j,i) = w%sx22(j,i) + DT*(w%sxp22(j,i) - fadx*w%sx22(j,i))
    
                    !equation 3
                    w%sx31(j,i)=DT2*mo2*(1.0/(DY2))*(u1-u*w%s2(j,i)) &
                    +2*w%sx32(j,i)-w%sx33(j,i)
                    
                else if (sny1/=0 .or. sny2/=0) then
                !根据系数求得一阶偏微分的离散算子
                    n2loop: do n=0, 9
                        if(n<5) then
                            ux=ux+w%s2(j,i+n-5)*w%xs1(n)
                            uy=uy+w%s2(j+n-5,i)*w%xs1(n)
                        else
                            ux=ux+w%s2(j,i+n-4)*w%xs1(n)
                            uy=uy+w%s2(j+n-4,i)*w%xs1(n)
                        end if
                    end do n2loop
    
                    !equation 1
                    w%sx11(j,i)=mo2*DT2*(u1-u*w%s2(j,i))*(1.0/(DY2)) &
                    -fady*fady*DT2*w%sx12(j,i)+(2*w%sx12(j,i) &
                    -w%sx13(j,i))+DT*(2*fady*(w%sx13(j,i)-w%sx12(j,i)))
    
                    !equation 2 : 包含三阶偏微分,需将其拆解为一阶偏微分(p)的二阶导数离散求解
                    w%sxp21(j,i) = 2.0*w%sxp22(j,i) - w%sxp23(j,i) &
                    +DT2*(-1.0*mo2*faddy*(1.0/(DY))*(uy*t4) - 2.0*fady &
                    *(w%sxp22(j,i) - w%sxp23(j,i))/DT - fady*fady*w%sxp22(j,i))
                    w%sx21(j,i) = w%sx22(j,i) + DT*(w%sxp22(j,i) - fady*w%sx22(j,i))
    
                    !equation 3
                    w%sx31(j,i)=DT2*mo2*(1.0/(DX2))*(u2-u*w%s2(j,i)) &
                    +2*w%sx32(j,i)-w%sx33(j,i)
                end if
    
                if(snx1==0 .and. snx2==0 .and. sny1==0 .and. sny2==0) then
                    w%s3(j,i)=(mo2*DT2*(u2-u*w%s2(j,i))*(1.0/(DX2)) &
                    +DT2*mo2*(1.0/(DY2))*(u1-u*w%s2(j,i)))+2*w%s2(j,i)-w%s1(j,i)
                else
                    w%s3(j,i)=w%sx11(j,i)+w%sx21(j,i)+w%sx31(j,i)
                end if
                u1=0
                u2=0
                u=0
                ux=0
                uy=0
            end do jloop
        end do iloop
    
        w%sw(:,:)=w%s1(:,:)
        w%s1(:,:)=w%s2(:,:)
        w%s2(:,:)=w%s3(:,:)

        w%sw(:,:)=w%sx13(:,:)
        w%sx13(:,:)=w%sx12(:,:)
        w%sx12(:,:)=w%sx11(:,:)

        w%sw(:,:)=w%sxp23(:,:)
        w%sxp23(:,:)=w%sxp22(:,:)
        w%sxp22(:,:)=w%sxp21(:,:)

        w%sw(:,:)=w%sx22(:,:)
        w%sx22(:,:)=w%sx21(:,:)

        w%sw(:,:)=w%sx33(:,:)
        w%sx33(:,:)=w%sx32(:,:)
        w%sx32(:,:)=w%sx31(:,:)
    end subroutine wave2d_timeslicecal

    subroutine wavelet02(s, N, DT, hz)
        integer :: N
        real(4),intent(inout) :: s(0:N)
        real(4) :: DT, hz
        real(4) :: pi=3.1415926
        real(4) :: f,det
        integer :: k
        det=0.05*(30.0/hz)

        do k=0, N
            f=(pi)*(pi)*hz*hz*(k*DT-det)* &
            exp((-pi*pi*hz*hz*(k*DT-det)*(DT*k-det))) &
            *(3.0-2.0*pi*pi*hz*hz*(k*DT-det)*(DT*k-det))
            s(k)=f
        end do
    end subroutine wavelet02

    subroutine wavelet01(s, N, DT, hz)
        integer :: N
        real(4),intent(inout) :: s(0:N)
        real(4) :: DT, hz
        real(4) :: pi=3.1415926
        real(4) :: f,det
        integer :: k
        det=0.05*(30.0/hz)

        do k=0, N
            f=exp((-pi*pi*hz*hz*(k*DT-det)*(DT*k-det))) &
            *(1.0-2.0*pi*pi*hz*hz*(k*DT-det)*(DT*k-det))
            s(k)=f
        end do
    end subroutine wavelet01

end module wave2d_basic

program wave2d_test
    use wave2d_type
    use wave2d_basic
    
    implicit none
    
        integer :: k,nz,nx,nt
        real(4) :: dt
        real(4) :: s(0:3000)
        type(wave2d) :: w
        nt=3000
        nz=500
        nx=500
        dt=0.0005
    
        call wave2d_creater(w,nz,nx)
        call wavelet01(s, nt, dt, 30.0)
    
        open( 12 , File = 'testmovie.bin' , Access = 'stream' , Form = 'Unformatted'  )
    
        do k=0, nt
            w%s2(nz/2,nx/2)=w%s2(nz/2,nx/2)+s(k)
            call wave2d_timeslicecal(w)
            if(mod(k,10)==0) then
                Write( 12  ) w%s2
            end if 
            if(mod(k,100)==0) then
                print *, k
            end if 
        end do
    
        close( 12 )
        call wave2d_del(w)
    
end program wave2d_test

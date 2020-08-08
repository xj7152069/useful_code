!/*********(version 1.0)***********/
!*
!   fortran module file: 
!/
!/********************************/

module  wave2d_basic
implicit none   

    type wave2d
        real,target  :: xs2(5)=[1.666667,-0.238095,0.039683,-0.004960,0.000317]
        real,target  :: xs1(10)=[-0.0007926,0.00991800,-0.0595200,0.238080,-0.833333,&
        0.833333,-0.238080,0.0595200,-0.00991800,0.0007926]
        real :: dx,dz,dt,PML_wide,R
        integer :: nx,nz,suface
        real,allocatable,target :: sx11(:,:)
        real,allocatable,target :: sx12(:,:)
        real,allocatable,target :: sx13(:,:)
        real,allocatable,target :: sxp21(:,:)
        real,allocatable,target :: sxp22(:,:)
        real,allocatable,target :: sxp23(:,:)
        real,allocatable,target :: sx21(:,:)
        real,allocatable,target :: sx22(:,:) 
        real,allocatable,target :: sx31(:,:)
        real,allocatable,target :: sx32(:,:)
        real,allocatable,target :: sx33(:,:)
        real,allocatable,target :: vel(:,:)
        real,allocatable,target :: s1(:,:)
        real,allocatable,target :: s2(:,:)
        real,allocatable,target :: s3(:,:)
    end	type wave2d

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

        allocate(w%s1(0:nz, 0:nx), stat=ierr)
        if(ierr.ne.0)then
            write(*,*)"Can not allocate working memory about us1, stop!!!"
            stop
        endif

        allocate(w%s2(0:nz, 0:nx), stat=ierr)
        if(ierr.ne.0)then
            write(*,*)"Can not allocate working memory about us1, stop!!!"
            stop
        endif

        allocate(w%s3(0:nz, 0:nx), stat=ierr)
        if(ierr.ne.0)then
            write(*,*)"Can not allocate working memory about us1, stop!!!"
            stop
        endif

        allocate(w%sx11(0:nz, 0:nx), stat=ierr)
        if(ierr.ne.0)then
            write(*,*)"Can not allocate working memory about us1, stop!!!"
            stop
        endif

        allocate(w%sx12(0:nz, 0:nx), stat=ierr)
        if(ierr.ne.0)then
            write(*,*)"Can not allocate working memory about us1, stop!!!"
            stop
        endif

        allocate(w%sx13(0:nz, 0:nx), stat=ierr)
        if(ierr.ne.0)then
            write(*,*)"Can not allocate working memory about us1, stop!!!"
            stop
        endif

        allocate(w%sx21(0:nz, 0:nx), stat=ierr)
        if(ierr.ne.0)then
            write(*,*)"Can not allocate working memory about us1, stop!!!"
            stop
        endif

        allocate(w%sx22(0:nz, 0:nx), stat=ierr)
        if(ierr.ne.0)then
            write(*,*)"Can not allocate working memory about us1, stop!!!"
            stop
        endif

        allocate(w%sx31(0:nz, 0:nx), stat=ierr)
        if(ierr.ne.0)then
            write(*,*)"Can not allocate working memory about us1, stop!!!"
            stop
        endif

        allocate(w%sx32(0:nz, 0:nx), stat=ierr)
        if(ierr.ne.0)then
            write(*,*)"Can not allocate working memory about us1, stop!!!"
            stop
        endif

        allocate(w%sx33(0:nz, 0:nx), stat=ierr)
        if(ierr.ne.0)then
            write(*,*)"Can not allocate working memory about us1, stop!!!"
            stop
        endif

        allocate(w%sxp21(0:nz, 0:nx), stat=ierr)
        if(ierr.ne.0)then
            write(*,*)"Can not allocate working memory about us1, stop!!!"
            stop
        endif

        allocate(w%sxp22(0:nz, 0:nx), stat=ierr)
        if(ierr.ne.0)then
            write(*,*)"Can not allocate working memory about us1, stop!!!"
            stop
        endif

        allocate(w%sxp23(0:nz, 0:nx), stat=ierr)
        if(ierr.ne.0)then
            write(*,*)"Can not allocate working memory about us1, stop!!!"
            stop
        endif

        allocate(w%vel(0:nz, 0:nx), stat=ierr)
        if(ierr.ne.0)then
            write(*,*)"Can not allocate working memory about us1, stop!!!"
            stop
        endif

        s1(:,:)=0.0
        s2(:,:)=0.0
        s3(:,:)=0.0
        sx11(:,:)=0.0
        sx12(:,:)=0.0
        sx13(:,:)=0.0
        sxp21(:,:)=0.0
        sxp22(:,:)=0.0
        sxp23(:,:)=0.0
        sx21(:,:)=0.0
        sx22(:,:)=0.0
        sx31(:,:)=0.0
        sx32(:,:)=0.0
        sx33(:,:)=0.0
        vel(:,:)=3000.0
    end subroutine wave2d_creater

    subroutine wave2d_del(w)
        type(wave2d),intent (out) :: w
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

    subroutine timeslicecal(w)
        type(wave2d),intent (inout) :: w
        real :: DX,DY,DT,xshd
        integer :: X,Y,suface_PML
        real :: fadx,fady,faddx,faddy,snx1,sny1,snx2,sny2,t2,t5
        integer :: i,j,n,t3,t4
        real :: u1,u2,u,ux,uy
        real :: C_X
        real :: C_Y
        real :: DT2,DT3,DX2,DY2,mo2
    
        real,pointer :: xs1_in(:)=>w%xs1, xs2_in(:)=>w%xs2    
        real,pointer :: sx11_in(:,:)=>w%sx11, sx12_in(:,:)=>w%sx12, sx13_in(:,:)=>w%sx13 
        real,pointer :: sxp21i(:,:)=>w%sxp21, sxp22i(:,:)=>w%sxp22, sxp23i(:,:)=>w%sxp23
        real,pointer :: sx21_in(:,:)=>w%sx21, sx22_in(:,:)=>w%sx22 
        real,pointer :: sx31_in(:,:)=>w%sx31, sx32_in(:,:)=>w%sx32, sx33_in(:,:)=>w%sx33 
        real,pointer :: p2_in(:,:)=>w%vel, swap=>null()
        real,pointer :: s1_in(:,:)=>w%s1, s2_in(:,:)=>w%s2, s3_in(:,:)=>w%s3 
    
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
        t5=float(Y)/X
        u1=0
        u2=0
        u=0
        ux=0
        uy=0
        C_X=log(w%R)*3/2/(xshd)/(xshd)/(xshd)/DX2/DX
        C_Y=log(w%R)*3/2/(xshd)/(xshd)/(xshd)/DY2/DY

        iloop: do i=5, X-5

            if (i>=0.5*(X)) then
                t3=1
            else
                t3=-1
            end if
    
            jloop: do j=5, Y-5
                
            !根据系数求得二阶偏微分的离散算子
                nloop: do n=0, 5
                    u=u+2*xs2_in(n);
                    u1=u1+xs2_in(n)*(s2_in(j-n-1,i)+s2_in(j+n+1,i))
                    u2=u2+xs2_in(n)*(s2_in(j,i-n-1)+s2_in(j,i+n+1))
                end do nloop
    
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
                    fady=p2_in(j,i)*C_Y*sny1*sny1*DY2
                    faddy=p2_in(j,i)*C_Y*2*sny1*DY
                end if
                if (sny2 /= 0) then
                    fady=p2_in(j,i)*C_Y*sny2*sny2*DY2
                    faddy=p2_in(j,i)*C_Y*2*sny2*DY
                end if
                if (snx1 /= 0) then
                    fadx=p2_in(j,i)*C_X*snx1*snx1*DX2
                    faddx=p2_in(j,i)*C_X*2*snx1*DX
                end if
                if (snx2 /= 0) then 
                    fadx=p2_in(j,i)*C_X*snx2*snx2*DX2
                    faddx=p2_in(j,i)*C_X*2*snx2*DX
                end if
    
                if(j>=0.5*(Y)) then
                    t4=1
                else
                    t4=-1
                end if
    
                mo2=p2_in(j,i)*p2_in(j,i)
    
                if(snx1/=0 .or. snx2/=0) then
                !根据系数求得一阶偏微分的离散算子
                    nloop: do n=0, 10
                        if(n<5) then
                            ux=ux+s2_in(j,i+n-5)*xs1_in(n)
                            uy=uy+s2_in(j+n-5,i)*xs1_in(n)
                        else
                            ux=ux+s2_in(j,i+n-4)*xs1_in(n)
                            uy=uy+s2_in(j+n-4,i)*xs1_in(n)
                        end if
                    end do nloop
    
                    !equation 1
                    sx11_in(j,i)=mo2*DT2*(u2-u*s2_in(j,i))*(1.0/(DX2)) &
                    -fadx*fadx*DT2*sx12_in(j,i)+(2*sx12_in(j,i) &
                    -sx13_in(j,i))+DT*(2*fadx*(sx13_in(j,i)-sx12_in(j,i)))
    
                    !equation 2
                    sxp21i(j,i) = 2.0*sxp22i(j,i) - sxp23i(j,i) &
                    +DT2*(-1.0*mo2*faddx*(1.0/(DX))*(ux*t3) - &
                    2.0*fadx*(sxp22i(j,i) - sxp23i(j,i))/DT - fadx*fadx*sxp22i(j,i))
                    sx21_in(j,i) = sx22_in(j,i) + DT*(sxp22i(j,i) - fadx*sx22_in(j,i))
    
                    !equation 3
                    sx31_in(j,i)=DT2*mo2*(1.0/(DY2))*(u1-u*s2_in(j,i)) &
                    +2*sx32_in(j,i)-sx33_in(j,i)
                    
                else if (sny1/=0 .or. sny2/=0) then
                !根据系数求得一阶偏微分的离散算子
                    nloop: do n=0, 10
                        if(n<5) then
                            ux=ux+s2_in(j,i+n-5)*xs1_in(n)
                            uy=uy+s2_in(j+n-5,i)*xs1_in(n)
                        else
                            ux=ux+s2_in(j,i+n-4)*xs1_in(n)
                            uy=uy+s2_in(j+n-4,i)*xs1_in(n)
                        end if
                    end do nloop
    
                    !equation 1
                    sx11_in(j,i)=mo2*DT2*(u1-u*s2_in(j,i))*(1.0/(DY2)) &
                    -fady*fady*DT2*sx12_in(j,i)+(2*sx12_in(j,i) &
                    -sx13_in(j,i))+DT*(2*dy*(sx13_in(j,i)-sx12_in(j,i)))
    
                    !equation 2 : 包含三阶偏微分,需将其拆解为一阶偏微分(p)的二阶导数离散求解
                    sxp21i(j,i) = 2.0*sxp22i(j,i) - sxp23i(j,i) &
                    +DT2*(-1.0*mo2*faddy*(1.0/(DY))*(uy*t4) - 2.0*fady &
                    *(sxp22i(j,i) - sxp23i(j,i))/DT - fady*fady*sxp22i(j,i))
                    sx21_in(j,i) = sx22_in(j,i) + DT*(sxp22i(j,i) - fady*sx22_in(j,i))
    
                    !equation 3
                    sx31_in(j,i)=DT2*mo2*(1.0/(DX2))*(u2-u*s2_in(j,i)) &
                    +2*sx32_in(j,i)-sx33_in(j,i)
                end if
    
                if(snx1==0 .and. snx2==0 .and. sny1==0 .and. sny2==0) then
                    s3_in(j,i)=(mo2*DT2*(u2-u*s2_in(j,i))*(1.0/(DX2)) &
                    +DT2*mo2*(1.0/(DY2))*(u1-u*s2_in(j,i)))+2*s2_in(j,i)-s1_in(j,i)
                else
                    s3_in(j,i)=sx11_in(j,i)+sx21_in(j,i)+sx31_in(j,i)
                end if
                u1=0
                u2=0
                u=0
                ux=0
                uy=0
            end do jloop
        end do iloop
    
        swap=>s1_in
        s1_in=>s2_in
        s2_in=>s3_in
        s3_in=>swap

        swap=>sx13_in
        sx13_in=>sx12_in
        sx12_in=>sx11_in
        sx11_in=>swap

        swap=>sxp23i
        sxp23i=>sxp22i
        sxp22i=>sxp21i
        sxp21i=>swap

        swap=>sx22_in
        sx22_in=>sx21_in
        sx21_in=>swap

        swap=>sx33_in
        sx33_in=>sx32_in
        sx32_in=>sx31_in
        sx31_in=>swap
    end subroutine timeslicecal

    subroutine wavelet02(s, N, DT, hz)
        integer :: N
        real,intent(inout) :: s(N)
        real :: DT, hz
        real :: pi=3.1415926
        real :: f,det
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
        real,intent(inout) :: s(N)
        real :: DT, hz
        real :: pi=3.1415926
        real :: f,det
        integer :: k
        det=0.05*(30.0/hz)

        do k=0, N
            f=exp((-pi*pi*hz*hz*(k*DT-det)*(DT*k-det))) &
            *(1.0-2.0*pi*pi*hz*hz*(k*DT-det)*(DT*k-det))
            s(k)=f
        end do
    end subroutine wavelet01

end module wave2d_basic


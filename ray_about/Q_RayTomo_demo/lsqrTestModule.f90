	module lsqrTestModule
		use lsqrDataModule
		use lsqrblasInterface
		use lsqrModule
		use lsqrCheckModule

		implicit none

		public	::	solve	!slove equations by lsqr; parameter definition
		public	::  allocatememory, deallocatememory
		public	::	initialization
		private	::	Aprod1, Aprod2

		integer(ip),parameter :: maxnz=100000000	!max number of none-zero element
		integer(ip)	::nzz							!true number of none-zero element
		integer(ip),allocatable	::indx(:)			!Row(Hang)
		integer(ip),allocatable ::jndx(:)			!Column(Lie)
		real(dp),allocatable	::path(:)

	contains

		subroutine allocatememory()

			allocate(indx(maxnz))
			allocate(jndx(maxnz))
			allocate(path(maxnz))

			return

		end subroutine allocatememory

		subroutine deallocatememory()

			deallocate(indx)
			deallocate(jndx)
			deallocate(path)

			return

		end subroutine deallocatememory

		subroutine initialization()
			
			integer(ip)	::	i

			nzz=0
			do i=1,maxnz
				indx(i)=0
				jndx(i)=0
				path(i)=0
			end do

			return

		end subroutine initialization

		subroutine Aprod1(m,n,x,y)	!y=y+A*x
			
			integer(ip), intent(in)  :: m, n
			real(dp), intent(in)     :: x(n)
			real(dp), intent(inout)  :: y(m)

			integer(ip)	::	i,j,k

			do k=1,nzz
				i=indx(k)
				j=jndx(k)

				y(i)=y(i)+path(k)*x(j)

			end do

			return

		end subroutine Aprod1

		subroutine Aprod2(m,n,x,y) !x=x+A'*y

			integer(ip), intent(in)  :: m, n
			real(dp), intent(inout)  :: x(n)
			real(dp), intent(in)     :: y(m)

			integer(ip) ::  i,j,k

			do k=1,nzz
				i=indx(k)
				j=jndx(k)

				x(j)=x(j)+path(k)*y(i)

			end do
			
			return

		end subroutine Aprod2

		subroutine solve(m, n, x, b)

			integer(ip),intent(in)	::	m, n	!m hang; n lie
			real(dp), intent(out)   ::  x(n)
			real(dp), intent(in)	::	b(m)
			
			integer(ip)	::	i, nout
			integer(ip)	::	itn, itnlim, istop, inform
			real(dp)	::	atol, btol, conlim, damp
			real(dp)	::	Anorm, Acond, Arnorm, rnorm, xnorm
			logical     ::	wantse
			real(dp),allocatable	::	se(:)
!real(dp)::t
!real(dp),parameter::one=1.0_dp
			allocate(se(n))

			print*,'m=',m,'n=',n

			print*,'nzz=',nzz

			open(11,file='lsqr.log')

			nout=11

!do i=1,n
!	t=t+one
!	x(i)=sqrt(t)
!	print*,x(i)
!end do
!			do i=1,m
!				if(b(i).ne.0.0) then
!					write(11,*) b(i)
!					PRINT*,i, b(i)
!				end if
!			end do
!			stop

!			do i=1,nzz
!				write(11,*) i,indx(i),jndx(i),path(i)
!				print*, i,indx(i),jndx(i),path(i)
!			end do
!			stop

			wantse = .false.
			atol   = 5e-3_dp                            ! Fixed accuracy.
			btol   = atol
!			conlim = 1e4_dp
			conlim = 1/(10*sqrt(eps))
!			itnlim = 4*(m + n + 50)
			itnlim = 30000
!			damp   = 1e-6_dp
			damp   = 0_dp

!print*,'lsqr-143'
!			call Acheck(m, n, Aprod1, Aprod2, nout, inform)	! Check Aprod1, Aprod2.
!print*,'lsqr-143'

			call LSQR  (m, n, Aprod1, Aprod2, b, damp, wantse,   &
						x, se,                                   &
						atol, btol, conlim, itnlim, nout,        &
						istop, itn, Anorm, Acond, rnorm, Arnorm, xnorm )

!			call xcheck(m, n, Aprod1, Aprod2, b, damp, x, &    ! Check x
!						Anorm, atol, nout,                &
!						inform)

!			do i=1,n
!				if(x(i).ne.0.0) write(nout,*) x(i)
!			end do

			close(11)
			deallocate(se)

			return

		end subroutine solve

	end module lsqrTestModule

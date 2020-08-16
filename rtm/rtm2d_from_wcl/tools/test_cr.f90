

include '../include/global.f90'

	  program   test

		use global
	  	implicit none

		character(len=para_char_flen)	:: file_cs = &
			'/ssd/ssd_3.2T/wcl/data/layer/2d/cs/layer_shot31_nt8000.dat'
		integer	::	lt=8000
		integer	::	nx=601

		real	::	dx=10.0
		real	::	dt=0.0005

		real,allocatable	::	shot(:,:)
		real,allocatable	::	tt(:,:)

		integer	::	it, ix
		real	::	a, b, c
		real	::	v1, v2, v3, v4
		real	::	x1, x2, x3, x4, x
		real	::	d1, d2, d3, d4, dd
		real	::	t0, t1, t2, t3, t4

		allocate(shot(lt, nx))
		allocate(tt(nx, 4))
		shot=0.0
		tt=0.0

		open(11, file=file_cs, form='unformatted', access='direct', recl=lbyte*lt*nx, status='old')
		read(11, rec=1)((shot(it, ix),it=1, lt), ix=1, nx)
		close(11)





		v1=2000.0
		d1=1000.0

		v2=2600.00
		d2=1100.00

		v3=3100.00
		d3=1100.00

		v4=4000.0
		d4=800.0

		do ix=1, nx
			x=(ix-1)*dx
			x=x/2.0
			dd=sqrt(x*x+d1*d1)
			tt(ix, 1)=dd*2/v1/dt

			x1=d1*x/(d1+d2)
			x2=x-x1
			t1=sqrt(d1*d1+x1*x1)/v1
			t2=sqrt(d2*d2+x2*x2)/v2
			tt(ix, 2)=(t1+t2)*2.0/dt


			x1=d1*x/(d1+d2+d3)
		!	x2=(d1+d2)*x/(d1+d2+d3)
		!	x2=x2-x1
			x2=d2*x/(d1+d2+d3)
			x3=x-x1-x2
			t1=sqrt(d1*d1+x1*x1)/v1
			t2=sqrt(d2*d2+x2*x2)/v2
			t3=sqrt(d3*d3+x3*x3)/v3
			tt(ix, 3)=(t1+t2+t3)*2.0/dt


			x1=d1*x/(d1+d2+d3+d4)
			x2=d2*x/(d1+d2+d3+d4)
			x3=d3*x/(d1+d2+d3+d4)
			x4=x-x1-x2-x3
			t1=sqrt(d1*d1+x1*x1)/v1
			t2=sqrt(d2*d2+x2*x2)/v2
			t3=sqrt(d3*d3+x3*x3)/v3
			t4=sqrt(d4*d4+x4*x4)/v4
			tt(ix, 4)=(t1+t2+t3+t4)*2.0/dt

		enddo


		do ix=1, nx
			do it=1, lt
				
				t0=tt(ix,1)-50
				t2=tt(ix,1)+50
				if( it .gt. t0 .and. it .lt. t2)then
				
				else

					t0=tt(ix,2)-50
					t2=tt(ix,2)+50
					if( it .gt. t0 .and. it .lt. t2)then

					else

						t0=tt(ix,3)-50
						t2=tt(ix,3)+50
						if( it .gt. t0 .and. it .lt. t2)then

						else
							t0=tt(ix,4)-50
							t2=tt(ix,4)+50
							if( it .gt. t0 .and. it .lt. t2)then

							else
								shot(it, ix)=0.0
							endif

						endif

					endif
				endif
			enddo
		enddo

		open(11, file='test1111', form='unformatted', access='direct', recl=lbyte*lt*nx, status='replace')
		write(11, rec=1)((shot(it, ix),it=1, lt), ix=1, nx)
		close(11)

		deallocate(shot)


	  end

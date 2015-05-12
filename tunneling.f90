program tunneling
	use, intrinsic :: iso_c_binding
	use modules
	implicit none
	include 'fftw3.f03'

	type(C_PTR) :: plan
	integer :: N, i 
	complex(8), dimension( :), allocatable :: PSI
	real(8), dimension(:), allocatable :: P, V
	real(8) :: k, alpha, P_total, t!, sigma
	real(8), parameter :: pi = 4._8*datan(1._8), delta_t = 0.1, sigma = 1, p_0 = 1, x_0 = 700, L = 100
	complex(8), parameter :: nim = (0,-1)


	N = 1024
	t = 0	
	
	allocate ( PSI(N) )
	allocate ( P(N) )
	allocate ( V(N) )


	k = 3
	alpha = 20d0
	P_total = 0d0
	V = 0d0

	write(*,*) "set xrange [0:",N,"]"
	write(*,*) "set yrange [-0.05:0.1]"

	



	call InitializeWavePacket

	do i=1,N
		PSI(i) = sqrt(1/(sqrt(pi)*alpha)) * exp(-(i/1d0-100)**2/(2*alpha**2)) * cmplx(cos(k*(i/1d0)*1d0) , sin(k*(i/1d0)*1d0))
		P(i) = real(PSI(i))**2 + aimag(PSI(i))**2
		P_total = P_total + P(i)
		if ( i > 0.4 * N .and. i< 0.6*N) then
			V(i) = 3
		end if
	end do


	write(*,*) "plot 'probabilities.dat' u 1:2 w l notitle, 'probabilities.dat' u 1:3 w l notitle"

	do while ( t <  1000)
		t = t + delta_t
		call NextTimeStep(PSI,N, delta_t, V, L)
	end do

contains
	
	subroutine InitializeWavePacket


end program



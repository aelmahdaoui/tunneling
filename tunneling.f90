program tunneling
	use, intrinsic :: iso_c_binding
	implicit none
	include 'fftw3.f03'

	type(C_PTR) :: plan
	integer :: N, i 
	complex(8), dimension( :), allocatable :: PSI
	real(8), dimension(:), allocatable :: P, V
	real(8) :: k, alpha, P_total, t, delta_t
	real(8), parameter :: pi = 4._8*datan(1._8)

	call ReadInputFile
	allocate ( PSI(N) )
	allocate ( P(N) )
	allocate ( V(N) )

	call InitializeWave
	call InitializePotential

	t = 0
	do while ( t <  10000*delta_t)
		t = t + delta_t
		call NextTimeStep
		if (int(t/delta_t) == 5000) then
			call WriteWaveToFile		!Snapshot taken
		end if
	end do

contains

	subroutine ReadInputFile	
		open(1, file='TunnelingInput.txt', status='old', action='read')
	  	read(1,*), N
		read(1,*), alpha
		read(1,*), k
		read(1,*), delta_t
		close(1)
	end subroutine
	
	subroutine InitializeWave
		do i=1,N
			PSI(i) = sqrt(1/(sqrt(pi)*alpha)) * exp(-(i/1d0-100)**2/(2*alpha**2)) * cmplx(cos(k*(i/1d0)*1d0) , sin(k*(i/1d0)*1d0))
		end do
	end subroutine

	subroutine InitializePotential
		V = 0d0
		do i=1,N
			if ( i > 0.4 * N .and. i< 0.6*N) then
				V(i) = 3
			end if
		end do
	end subroutine

	subroutine NextTimeStep
		call ApplyPotentialOperator
		call ApplyFourierTransform
		call ApplyKineticOperator
		call ApplyInverseFourierTransform
		call NormalizeWave
	end subroutine

	subroutine ApplyPotentialOperator
		do i=1,N
			PSI(i) = exp(cmplx(0.0, -delta_t * V(i))) * PSI(i)
		end do
	end subroutine

	subroutine ApplyFourierTransform
		plan = fftw_plan_dft_1d(N, PSI, PSI, FFTW_FORWARD,FFTW_ESTIMATE)
		call fftw_execute_dft(plan, PSI, PSI)
		call fftw_destroy_plan(plan)
	end subroutine

	subroutine ApplyKineticOperator
		do i=1,N
			if ( i < N/2) then
				PSI(i) = exp(cmplx(0.0, -delta_t * 0.5 * (2*pi*i/N)**2)) * PSI(i)		
			else 
				PSI(i) = exp(cmplx(0.0, -delta_t * 0.5 * (2*pi*i/N - 2*pi)**2)) * PSI(i)
			end if
		end do
	end subroutine

	subroutine ApplyInverseFourierTransform
		plan = fftw_plan_dft_1d(N, PSI, PSI, FFTW_BACKWARD,FFTW_ESTIMATE)
		call fftw_execute_dft(plan, PSI, PSI)
		call fftw_destroy_plan(plan)
	end subroutine	

	subroutine NormalizeWave
		P_total = 0
		do i = 1,N
			P(i) = real(PSI(i))**2 + aimag(PSI(i))**2
			P_total = P_total + P(i)
		end do
		PSI(:) = PSI(:)/sqrt(P_total)		
	end subroutine
	
	subroutine WriteWaveToFile
		open (unit=3,file="probabilities.dat",action="write", status="replace")
		do i = 1,N
			P(i) = real(PSI(i))**2 + aimag(PSI(i))**2
			write(3,*) i, P(i), V(i)
		end do
		close(3)
	end subroutine
end program



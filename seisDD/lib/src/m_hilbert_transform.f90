!! module for hilbert transform
!! created by Yanhua O. Yuan ( yanhuay@princeton.edu)

module m_hilbert_transform
    use constants, only : CUSTOM_REAL
    implicit none 

    complex :: CI=(0.0,1.0)

contains

    subroutine hilbert(trace,nt)
        !! return hilbert transform of real signal trace(nt)
        !! i.e. imaginary part of analytic signal 
        !! a = cmplx(trace, hilbert(trace))

        integer,                  intent(in)    :: nt
        real(kind=CUSTOM_REAL),         intent(inout) :: trace(nt)

        complex, allocatable, dimension(:) :: C
        integer :: NPT,IMID   

        ! extend nt to a power of 2
        IF ( nt <= 0 ) STOP 'FATAL ERROR in HILBERT: nt must be positive'

        NPT = 2**( INT( LOG10( REAL( nt ) ) / 0.30104 ) + 1 )
        ! IF ( NPT /= nt) print*,'pad trace from length ', nt, ' to ',NPT
        IF (NPT > 16784) STOP 'FATAL ERROR in HILBERT: nt(NPT) exceeds 16784 '

        allocate(C(NPT))
        C=cmplx(0.,0.)
        C(1:nt)=cmplx(trace(1:nt),0.0)

        ! Fourier transform 
        call CFFT(C,NPT,1)
        ! scaling 
        C=C/NPT

        !  Multiply by i * sgn( f )
        IMID = NPT / 2
        C( 1:IMID-1 ) = -CI * C( 1:IMID-1 )   ! pos. spectrum (-i)
        C( IMID     ) = 0.0                   ! d.c. component
        C(IMID+1:NPT) = CI * C( IMID+1:NPT )   ! neg. spectrum (i)

        ! inverse Fourier transform
        call  CFFT(C,NPT,-1) 

        ! output
        trace(1:nt)=real(C(1:nt))

        deallocate(C)

    end subroutine hilbert

    SUBROUTINE CFFT(trace,N,iforw)
        !!! complex FFT
        IMPLICIT NONE  
        integer,         intent(in)    :: N, iforw
        complex,         intent(inout) :: trace(N)
        INTEGER :: I1, I2A, I2B, I3, I3Rev, IP1, IP2, ISign
        REAL    ::  theta, sinth
        COMPLEX :: TEMP, W, WSTP

        ISIGN = -IFORW
        I3REV = 1

        DO I3 = 1, N
        IF ( I3 < I3REV ) THEN   ! switch values
            TEMP          = trace( I3 )
            trace( I3 )    = trace( I3REV )
            trace( I3REV ) = TEMP
        ENDIF
        ! following loop is just to compute I3REV
        IP1 = N / 2
        DO WHILE (  I3REV > IP1 )
        IF ( IP1 <= 1 ) EXIT
        I3REV = I3REV - IP1
        IP1   = IP1 / 2
        END DO
        I3REV = I3REV + IP1
        END DO

        IP1 = 1

        DO WHILE  ( IP1 < N )
        IP2   = IP1 * 2
        THETA = 6.283185307 / FLOAT( ISIGN * IP2 )
        SINTH = SIN( THETA / 2.)
        WSTP  = CMPLX( -2. * SINTH * SINTH, SIN( THETA ) )
        W     = 1.

        DO I1 = 1, IP1
        DO I3 = I1, N, IP2
        I2A = I3
        I2B = I2A + IP1
        TEMP      = W * trace( I2B )
        trace( I2B ) = trace( I2A ) - TEMP
        trace( I2A ) = trace( I2A ) + TEMP
        END DO

        W = W * WSTP + W
        END DO

        IP1 = IP2
        END DO

        RETURN
    END SUBROUTINE CFFT

end module m_hilbert_transform

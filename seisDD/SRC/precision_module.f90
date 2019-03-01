module precision
  ! Uncomment the appropriate line to select the desired precision
  
  integer, parameter  :: RK = KIND(1.0) ! Single prec. (RK: for Real Kind, for my workstation's single prec. fftw)
  !integer, parameter  :: RK = KIND(1.0D0) ! Double prec. (for Shaheen fftw2_dp, or for my workstation's double prec. fftw)
  integer, parameter  :: RK1= KIND(1.0)  !! adjustable, set to single or double, used in marine_module
  integer, parameter  :: RKw= KIND(1.0D0)!! adjustable, used in FWI
  
  !integer, parameter  :: PK = 4  ! pointer on 32-bit machine has 4 bytes
  integer, parameter  :: PK = 8  ! pointer on 64-bit machine has 8 bytes
end module precision

! Then in every program unit,
! USE precision
!
! Remember that the use statement must immediately follow the 
! unit invocation (program, subroutine, function). Then declare 
! all real variables using the form:
! INTEGER(PK):: plan
! REAL (RK)  :: r, d
! COMPLEX(RK):: c, z
! and so forth. It will be necessary to recompile every time 
! you change from real to double or vice versa, but this module 
! makes it much simpler. 

! --------- This type is adaptive, unlike Num. Recp.'s SP, DP, SPC, DPC,
! which are intercorrelated with other functions, so can't be easily changed.
! Here, make the type that's most likely to be changed adptive, e.g.,
! FFTW's precision.  On Shaheen Blue Gene it's double while on Linux
! cluster its single.

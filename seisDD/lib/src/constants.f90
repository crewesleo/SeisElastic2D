module constants

implicit none

integer,parameter :: IIN=40
integer, parameter :: IOUT=41
integer,parameter :: MAX_STRING_LEN= 65535
integer,parameter :: MAX_KERNEL_NUM=20
integer,parameter :: MAX_FILENAME_LEN=65535
integer,parameter :: MAX_LINES=100000000
integer,parameter :: MAX_MISFIT_TYPE=10
integer, parameter :: SIZE_REAL=4
integer, parameter :: SIZE_DOUBLE=8
integer, parameter :: CUSTOM_REAL=SIZE_REAL
integer, parameter :: CUSTOM_COMPLEX=2*SIZE_REAL 

!===================================================================
! filter parameters for xapiir subroutine (filter type is BP)
real(kind=CUSTOM_REAL), parameter :: TRBDNDW=0.3
real(kind=CUSTOM_REAL), parameter :: APARM=30.0
integer, parameter :: IORD=5
integer, parameter :: PASSES=2

! -------------------------------------------------------------
! array dimensions
! note that some integer arrays (iM,iL,iR) are NWINDOWS * NWINDOWS
! THESE SHOULD PROBABLY BE USER PARAMETERS, SINCE THEY WILL AFFECT
! THE SPEED OF THE PROGRAM (ALTHOUGH NOT THE OUTPUT).
!integer, parameter :: NDIM=10000
integer, parameter :: NWINDOWS=2500

! -------------------------------------------------------------
! miscellaneous - do not modify!
! -------------------------------------------------------------

! mathematical constants
real(kind=CUSTOM_REAL), parameter :: PI=3.1415926535897
real(kind=CUSTOM_REAL), parameter :: E= 2.7182818284590

! filter types
integer, parameter :: HANNING=1
integer, parameter :: HAMMING=2
integer, parameter :: COSINE= 3

! modified constants 

! add by YY
! constants
real(kind=CUSTOM_REAL), parameter :: TWOPI=2.0 * PI
complex (CUSTOM_REAL), parameter :: CCI=cmplx(0.,1.)
real(kind=CUSTOM_REAL), parameter :: LARGE_VAL=huge(0.0)
real(kind=CUSTOM_REAL), parameter :: SMALL_VAL=tiny(0.0)
! phase correction control parameters, set this between (PI, 2PI),
! use a higher value for conservative phase wrapping
real(kind=CUSTOM_REAL), parameter :: PHASE_STEP=1.5 * PI
! FFT parameters
integer, parameter :: LNPT=13, NPT=2**LNPT
real(kind=CUSTOM_REAL), parameter :: FORWARD_FFT=1.0
real(kind=CUSTOM_REAL), parameter :: REVERSE_FFT=-1.0
! CUTOFF for phase unwrapping 
real(kind=CUSTOM_REAL), parameter :: CUTOFF=PI
! water level for effective spectrum 
real(kind=CUSTOM_REAL), parameter :: WTR=0.05
real(kind=CUSTOM_REAL), parameter ::wtr_env=0.05
! water level for mtm 
real(kind=CUSTOM_REAL), parameter ::wtr_mtm=1.e-10
! multitaper
real(kind=CUSTOM_REAL), parameter :: mt_threshold=0.9! eigenvalue threshold
!integer, parameter :: MW=10 ! number of segments of uncorrelated frequency points 
real(kind=CUSTOM_REAL), parameter :: NW=3
integer, parameter :: NTAPER =int(2*NW-1)
! error estimation 
logical :: USE_ERROR_CC=.false.
! minimum error for dt and dlnA
real(kind=CUSTOM_REAL), parameter :: DT_SIGMA_MIN=1.0
real(kind=CUSTOM_REAL), parameter :: DLNA_SIGMA_MIN=0.5
logical :: USE_ERROR_MT=.false.
! taper power 
integer :: ipwr_w=10
real(kind=CUSTOM_REAL), parameter :: ipwr_t=10 ! for time-domain cosine taper 
! CG orhtogonality threshold for conscutive gradients
real(kind=CUSTOM_REAL), parameter :: CG_threshold=0.1

! WT parameters 
integer, parameter :: nvm=6
character(len=MAX_STRING_LEN) :: WT_directory='./WT_basis'
 
! gradient rho, Vp, Vs dim 
integer, parameter :: gdim_rho=3
integer, parameter :: gdim_vp=4
integer, parameter :: gdim_vs=5
! model rho, Vp, Vs dim 
integer, parameter :: mdim_rho=4
integer, parameter :: mdim_vp=5
integer, parameter :: mdim_vs=6

! Display 
logical :: DISPLAY_DETAILS=.false.
character(len=MAX_STRING_LEN) :: output_dir='OUTPUT_FILES'
! -------------------------------------------------------------

end module constants


module seismo_parameters
! yanhuay@princeton.edu

use constants, only: IIN, IOUT, MAX_STRING_LEN,MAX_FILENAME_LEN,MAX_KERNEL_NUM, &
    MAX_LINES,MAX_MISFIT_TYPE, SIZE_REAL, SIZE_DOUBLE, CUSTOM_REAL,CUSTOM_COMPLEX,&
    LARGE_VAL, SMALL_VAL,PI

implicit none

!----------------------------------------------------------------------
! constants
! number of Gauss-Lobatto-Legendre (GLL) points (i.e., polynomial degree + 1)
INTEGER, PARAMETER :: NGLLX=5
! number of Gauss-Lobatto-Jacobi (GLJ) points in the axial elements (i.e.,
! polynomial degree + 1)
! the code does NOT work if NGLLZ /= NGLLX because it then cannot handle a
! non-structured mesh
! due to non matching polynomial degrees along common edges
INTEGER, PARAMETER :: NGLLZ=5
INTEGER, PARAMETER :: NGLLY=1

!! solver
CHARACTER (LEN=20) :: solver='specfem2D'
CHARACTER (LEN=MAX_STRING_LEN) :: LOCAL_PATH='OUTPUT_FILES'
CHARACTER (LEN=MAX_STRING_LEN) :: IBOOL_NAME='ibool.bin'

!! FORWARD MODELNG INFO
INTEGER, PARAMETER :: NSTEP=4800 
REAL(KIND=CUSTOM_REAL), PARAMETER :: deltat=0.06 
REAL(KIND=CUSTOM_REAL), PARAMETER :: t0=0.0 
REAL(KIND=CUSTOM_REAL), PARAMETER :: f0=0.084 
INTEGER, PARAMETER :: NREC=2 
INTEGER, PARAMETER :: NSRC=1 

!! PRE-PROCESSING
! wavelet
INTEGER, PARAMETER :: Wscale=0 
!window
INTEGER, PARAMETER :: is_window=0 
INTEGER, PARAMETER :: window_type=3
REAL(KIND=CUSTOM_REAL), PARAMETER :: Vmax=3900 
REAL(KIND=CUSTOM_REAL), PARAMETER :: Vmin=3100
REAL(kind=CUSTOM_REAL), PARAMETER :: wavelet_len=3.0/f0
REAL(kind=CUSTOM_REAL), PARAMETER :: taper_len=1.2/f0

! damping
INTEGER, PARAMETER :: is_laplace=0
REAL(KIND=CUSTOM_REAL), PARAMETER :: X_decay=1.0
REAL(KIND=CUSTOM_REAL), PARAMETER :: T_decay=1.0
REAL(KIND=CUSTOM_REAL) :: lambda_x=1.0/X_decay 
REAL(KIND=CUSTOM_REAL) :: lambda_t=1.0/T_decay
! mute
INTEGER, PARAMETER :: mute_near=0 
REAL(KIND=CUSTOM_REAL), PARAMETER :: offset_near=0 
INTEGER, PARAMETER :: mute_far=0 
REAL(KIND=CUSTOM_REAL), PARAMETER :: offset_far=0 

!! event scale
REAL(KIND=CUSTOM_REAL), PARAMETER :: lambda_min=Vmin/f0
REAL(KIND=CUSTOM_REAL), PARAMETER :: lambda=Vmax/f0
REAL(KIND=CUSTOM_REAL), PARAMETER :: wavenumber=2.0*pi/lambda
REAL(KIND=CUSTOM_REAL), PARAMETER :: omega=2*pi*f0

!! sensitivity
LOGICAL :: sensitivity=.false.

!! DD
REAL(KIND=CUSTOM_REAL), PARAMETER :: cc_threshold=0.9
REAL(KIND=CUSTOM_REAL), PARAMETER :: DD_min=SMALL_VAL
REAL(KIND=CUSTOM_REAL), PARAMETER :: DD_max=LARGE_VAL
REAL(KIND=CUSTOM_REAL) :: ratio_data_syn=0.01

!! OPTIMIZATION
CHARACTER (LEN=2) :: opt_scheme='QN'
INTEGER, PARAMETER :: CGSTEPMAX=10 
CHARACTER (LEN=2) :: CG_scheme='PR'
INTEGER, PARAMETER :: BFGS_STEPMAX=4 
REAL(KIND=CUSTOM_REAL), PARAMETER :: initial_step_length=0.04 
INTEGER, PARAMETER :: max_step=5
REAL(KIND=CUSTOM_REAL), PARAMETER :: min_step_length=0.01 
LOGICAL :: backtracking=.false.

!! Scaling factor
!REAL(KIND=CUSTOM_REAL), PARAMETER :: scaling_factor=1.0

!! CONVERGENCE?
INTEGER, PARAMETER :: iter_start=1
INTEGER, PARAMETER :: iter_end=20
REAL(KIND=CUSTOM_REAL) :: misfit_ratio_initial=0.001
REAL(KIND=CUSTOM_REAL) :: misfit_ratio_previous=0.0001

!! POST-PROCESSING
LOGICAL :: smooth=.false.
LOGICAL :: MASK_SOURCE=.false.
LOGICAL :: MASK_STATION=.false.
LOGICAL :: MASK_DAMP=.true.
LOGICAL :: MASK_MODEL=.false.
REAL(KIND=CUSTOM_REAL), PARAMETER :: source_radius=8.0
REAL(KIND=CUSTOM_REAL), PARAMETER :: station_radius=4.0
INTEGER, PARAMETER :: mask_z=0
INTEGER, PARAMETER :: mask_zend=0

!! DISPLAY 
LOGICAL :: DISPLAY_DETAILS=.false.
!! COMPUTE ATTENUATION KERNELS by PWY
LOGICAL :: VISCOELASTIC=.false.
!! Using joint misfit or not for viscoelastic FWI
LOGICAL :: JOINT_MISFIT=.false.
!! Joint misfit tradeoff parameter
!! misfit_lambda*WD + (1-misfit_lambda)*RA
REAL(KIND=CUSTOM_REAL), PARAMETER :: misfit_lambda=0.0

!!!!!!!!!!!!!!!!! gloabl variables !!!!!!!!!!!!!!!!!!!!!!
INTEGER :: myrank,nproc,iproc

!! ADJOINT?
LOGICAL :: compute_adjoint

!! VISCOELASTIC?
!LOGICAL :: VISCOELASTIC
!! data
INTEGER :: ndata
INTEGER,PARAMETER :: MAX_DATA_NUM = 4
INTEGER, DIMENSION(:), ALLOCATABLE :: which_proc_receiver
INTEGER :: nrec_proc  ! trace from a single proc
REAL(KIND=CUSTOM_REAL), DIMENSION(:,:), ALLOCATABLE :: seism_obs
REAL(KIND=CUSTOM_REAL), DIMENSION(:,:), ALLOCATABLE :: seism_syn
REAL(KIND=CUSTOM_REAL), DIMENSION(:,:), ALLOCATABLE :: seism_adj
REAL(KIND=CUSTOM_REAL), DIMENSION(:,:), ALLOCATABLE :: seism_adj_AD
REAL(KIND=CUSTOM_REAL), DIMENSION(:,:), ALLOCATABLE :: seism_adj_DD
REAL(KIND=CUSTOM_REAL), DIMENSION(:), ALLOCATABLE :: st_xval,st_yval,st_zval
REAL(KIND=CUSTOM_REAL) :: x_source, y_source, z_source
INTEGER(KIND=4) :: r4head(60)
INTEGER(KIND=2) :: header2(2)

!! measurement
CHARACTER(LEN=MAX_STRING_LEN) :: measurement_list
CHARACTER(LEN=MAX_STRING_LEN) :: misfit_type_list
CHARACTER(LEN=MAX_STRING_LEN) :: measurement_attenuation

!! window 
INTEGER,DIMENSION(:), ALLOCATABLE  :: win_start, win_end

!! source-timefunction
LOGICAL :: conv_stf=.false.
CHARACTER (LEN=MAX_FILENAME_LEN) :: stf_file='source.txt'
REAL(KIND=CUSTOM_REAL) :: tshift_stf, integral_stf 
REAL(KIND=CUSTOM_REAL), DIMENSION(:), ALLOCATABLE :: stf
INTEGER :: stf_len

!! misfit
INTEGER :: num_AD, num_DD
INTEGER, DIMENSION(:,:), ALLOCATABLE :: is_pair
REAL(KIND=CUSTOM_REAL), DIMENSION(:), ALLOCATABLE :: misfit_proc
REAL(KIND=CUSTOM_REAL) :: misfit_AD,misfit_DD, misfit

!! kernels
INTEGER :: nspec
INTEGER, DIMENSION(:), ALLOCATABLE :: nspec_proc
INTEGER :: nker, npre
REAL(KIND=CUSTOM_REAL), DIMENSION(:), ALLOCATABLE :: g_new
REAL(KIND=CUSTOM_REAL), DIMENSION(:), ALLOCATABLE :: p_new
LOGICAL :: precond=.false.
REAL(KIND=CUSTOM_REAL) :: wtr_precond=1.0d-6

!! models
INTEGER :: nmod
REAL(KIND=CUSTOM_REAL), DIMENSION(:), ALLOCATABLE :: m_new
REAL(KIND=CUSTOM_REAL), DIMENSION(:), ALLOCATABLE :: m_try

!! linesearch
INTEGER :: is_done, is_cont, is_brak
REAL(KIND=CUSTOM_REAL) :: step_length, next_step_length,optimal_step_length

!! mask source
REAL(KIND=CUSTOM_REAL), DIMENSION(:,:,:,:), ALLOCATABLE :: xstore
REAL(KIND=CUSTOM_REAL), DIMENSION(:,:,:,:), ALLOCATABLE :: ystore
REAL(KIND=CUSTOM_REAL), DIMENSION(:,:,:,:), ALLOCATABLE :: zstore
REAL(KIND=CUSTOM_REAL), DIMENSION(:,:,:,:), ALLOCATABLE :: mask
!----------------------------------------------------------------------

end module seismo_parameters

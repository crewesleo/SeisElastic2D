!=====================================================================
!
!                   S P E C F E M 2 D  Version 7 . 0
!
!=====================================================================

! setup/precision.h.  Generated from precision.h.in by configure.

!
! solver in single or double precision depending on the machine
!
! set to MPI_REAL to run in single precision
! set to MPI_DOUBLE_PRECISION to run in double precision
!
! ALSO CHANGE FILE constants.h ACCORDINGLY
!
  integer, parameter :: CUSTOM_MPI_TYPE = MPI_REAL

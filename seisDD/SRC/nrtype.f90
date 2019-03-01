MODULE nrtype

  IMPLICIT NONE
  INTEGER, PARAMETER :: I4  = 4  ! a 4 byte number = how many record length units?
  ! On GNU,   recl unit = 1 byte, therefore I4 = 4
  ! On Intel, recl unit = 4 bytes, so I4 would = 1
  !  If -assume byterecl, then recl unit = 1 byte, and I4 = 4
  
  INTEGER, PARAMETER :: I8B = SELECTED_INT_KIND(18)
  INTEGER, PARAMETER :: I4B = SELECTED_INT_KIND(9)
  INTEGER, PARAMETER :: I2B = SELECTED_INT_KIND(4)
  INTEGER, PARAMETER :: I1B = SELECTED_INT_KIND(2)
  INTEGER, PARAMETER :: SP = KIND(1.0)
  INTEGER, PARAMETER :: DP = KIND(1.0D0)
  INTEGER, PARAMETER :: SPC = KIND((1.0,1.0))
  INTEGER, PARAMETER :: DPC = KIND((1.0D0,1.0D0)) ! print out those kind values for Windows
  ! I8B = 8,  I4B = 4,  I2B = 2,   I1B = 1,
  ! SP  = 4,   DP = 8,  SPC = 4,   DPC = 8.  So actually SP=SPC, DP=DPC.
  ! I.e., complex kind is the same as its underlying real kind. 
  ! (Of course total size is twice.) 
  
  INTEGER, PARAMETER  :: LGT = KIND(.true.)
  REAL(SP), PARAMETER :: PI=3.141592653589793238462643383279502884197_sp
  REAL(SP), PARAMETER :: PIO2=1.57079632679489661923132169163975144209858_sp
  REAL(SP), PARAMETER :: TWOPI=6.283185307179586476925286766559005768394_sp
  REAL(SP), PARAMETER :: SQRT2=1.41421356237309504880168872420969807856967_sp
  REAL(SP), PARAMETER :: EULER=0.5772156649015328606065120900824024310422_sp
  REAL(DP), PARAMETER :: PI_D=3.141592653589793238462643383279502884197_dp
  REAL(DP), PARAMETER :: PIO2_D=1.57079632679489661923132169163975144209858_dp
  REAL(DP), PARAMETER :: TWOPI_D=6.283185307179586476925286766559005768394_dp
  TYPE sprs2_sp
    INTEGER(I4B) :: n,len
    REAL(SP), DIMENSION(:), POINTER :: val
    INTEGER(I4B), DIMENSION(:), POINTER :: irow
    INTEGER(I4B), DIMENSION(:), POINTER :: jcol
  END TYPE sprs2_sp
  TYPE sprs2_dp
    INTEGER(I4B) :: n,len
    REAL(DP), DIMENSION(:), POINTER :: val
    INTEGER(I4B), DIMENSION(:), POINTER :: irow
    INTEGER(I4B), DIMENSION(:), POINTER :: jcol
  END TYPE sprs2_dp
  
  INTERFACE logic2byte  ! logical 1 or 4 bytes, to 1 byte integer
    module procedure  &
      logic_byte0,  logic_byte1,  logic_byte2, &
      logic1_byte0, logic1_byte1, logic1_byte2
  END INTERFACE  ! _{0,1,2}: dimensionality
  
CONTAINS !!!!!!!!!!!!!!

subroutine logic1_byte0(b, l)
integer(1), intent(out) :: b
logical(1), intent(in)  :: l
if (l) then
  b = 1
else
  b = 0
end if
end subroutine
subroutine logic_byte0(b, l)
integer(1), intent(out) :: b
logical,    intent(in)  :: l
if (l) then
  b = 1
else
  b = 0
end if
end subroutine  

subroutine logic1_byte1(b, l)
integer(1), intent(out) :: b(:)
logical(1), intent(in)  :: l(:)
integer k
do k = 1, size(l)
  if (l(k)) then
    b(k) = 1
  else
    b(k) = 0
  end if
end do
end subroutine
subroutine logic_byte1(b, l)
integer(1), intent(out) :: b(:)
logical,    intent(in)  :: l(:)
integer k
do k = 1, size(l)
  if (l(k)) then
    b(k) = 1
  else
    b(k) = 0
  end if
end do
end subroutine  
 
subroutine logic1_byte2(b, l)
integer(1), intent(out) :: b(:,:)
logical(1), intent(in)  :: l(:,:)
integer j, k
do j = 1, size(l,2)
  do k = 1, size(l,1)
    if (l(k,j)) then
      b(k,j) = 1
    else
      b(k,j) = 0
    end if
  end do
end do
end subroutine
subroutine logic_byte2(b, l)
integer(1), intent(out) :: b(:,:)
logical,    intent(in)  :: l(:,:)
integer j, k
do j = 1, size(l,2)
  do k = 1, size(l,1)
    if (l(k,j)) then
      b(k,j) = 1
    else
      b(k,j) = 0
    end if
  end do
end do
end subroutine
 
END MODULE nrtype

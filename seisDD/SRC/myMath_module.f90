MODULE myMath_module
! basic matrix and vector functions, following MATLAB names
  USE nrtype;   USE nrutil, ONLY : assert 
  IMPLICIT NONE
  
  INTERFACE Bspline_conv
    module procedure  Bspline1D_mirr_conv_S, Bspline1D_mirr_conv_D, &
            Bspline2D_mirr_conv_S,  Bspline3D_mirr_conv_S
  END INTERFACE ! Bspline1D_mirr_conv_S(A, W, tmp1, tmp2), ..._D(A, W, tmp), 
           !Bspline2D_mirr_conv_S(A, W)
           
  INTERFACE percentile
    module procedure  percentile1D
    module procedure  percentile2D
  END INTERFACE

  INTERFACE diff     ! like the MATLAB diff.  diff(A,dim) has 1 fewer elements along dim,
  ! e.g., diff([1 2 5]) = [1 3].  only 1d implemented.
    module procedure diff_rv
  END INTERFACE diff
  
  INTERFACE diff0   ! diff0(A,dim) has same size as A, 'coz  diff wrt A(0)==0, thus
  ! retaining A(1), as A(1)-A(0) = A(1), e.g.,
  ! diff([3 5 9]) = [3 2 4]
    module procedure diff0_rv, diff0_rm 
  END INTERFACE diff0
  
  INTERFACE ind2sub  ! index --> subscript:  sub = ind2sub(idx,siz)
    module procedure ind2sub_v, ind2sub_s 
  END INTERFACE ind2sub
  
  ! All accumulative operations involving single precision should be
  ! accumulated on double precision.  Otherwise the error could be too big!
  ! 1d vectors may be fine, but when the size is large, need to use double!
  
  INTERFACE norm_1    ! norm_1:  \Sum_i |x_i|
    module procedure norm_1_rv, norm_1_dv,  &
                     norm_1_rm, norm_1_dm,  &
                     norm_1_r3d,norm_1_d3d
  END INTERFACE
  
  INTERFACE norm_inf  !   max_i |x_i|
    module procedure norm_inf_rv, norm_inf_dv,  &
                     norm_inf_rm, norm_inf_dm,  &
                     norm_inf_r3d,norm_inf_d3d
  END INTERFACE
  
  INTERFACE norm2    ! square of Euclidean norm
    module procedure norm2_rv, norm2_dv, norm2_cv, norm2_zv, &
                     norm2_rm, norm2_dm, norm2_cm, norm2_zm, &
                     norm2_r3d,norm2_d3d,norm2_c3d,norm2_z3d
  END INTERFACE
  
  INTERFACE norm    ! Euclidean norm
    module procedure norm_rv, norm_dv, norm_cv, norm_zv, &
                     norm_rm, norm_dm, norm_cm, norm_zm, &
                     norm_r3d,norm_d3d,norm_c3d,norm_z3d
  END INTERFACE
  
  INTERFACE dot_prod    ! dot product.  0,1,2,3 dimensions.
    module procedure dot_prod_rv, dot_prod_dv, dot_prod_cv, dot_prod_zv, &
                     dot_prod_rm, dot_prod_dm, dot_prod_cm, dot_prod_zm, &
                     dot_prod_r3d,dot_prod_d3d,dot_prod_c3d,dot_prod_z3d
  END INTERFACE
  
  INTERFACE swap_pointers   ! swap_pointers(pa, pb)  ! pntr always pass by reference
    module procedure swap_pointers_1d, swap_pointers_1r, swap_pointers_1i, &
                     swap_pointers_2r, swap_pointers_2i, &
                     swap_pointers_3r, swap_pointers_3i
  END INTERFACE
  
  INTERFACE swapDimInplace  ! swapDimInplace(A,dim1,dim2), swapDimInplace(A)
    module procedure swapDimInplace2r, swapDimInplace3r, &
                     swapDimInplace2c, swapDimInplace3c, &
                     swapDimInplace2i, &
                     swapDimInplace2d, swapDimInplace3d, &
                     swapDimInplace2z, swapDimInplace3z
  END INTERFACE swapDimInplace ! A's length along dim1 & dim2 must be equal
  ! otherwise inplace swapping can't work because of different shape
  
  INTERFACE permute  ! rearrange dimensions of N-D array
  ! e.g., 
  ! call permute(Y, X, (/2 3 1/) ) 
  ! Y <-- X(2nd dim, 3rd dim, 1st dim), where Y's type is array pointer:
  !    ..., POINTER :: Y(:,:,:)
    module procedure permute_r, permute_c, permute_d, permute_z
  ! only 3d array implemented, 'coz 2d matrix can transpose
  END INTERFACE permute
  
  INTERFACE cat1
    module procedure cat1_iv, cat1_rv, cat1_cv, cat1_dv, cat1_zv, &
                     cat1_im, cat1_rm, cat1_cm, cat1_dm, cat1_zm
  END INTERFACE cat1 
  INTERFACE cat2
    module procedure cat2_iv, cat2_rv, cat2_cv, cat2_dv, cat2_zv, &
                     cat2_im, cat2_rm, cat2_cm, cat2_dm, cat2_zm
  END INTERFACE cat2
  INTERFACE catfun1
    module procedure &
      catfun1_iv, catfun1_rv, catfun1_cv, catfun1_dv, catfun1_zv, &
      catfun1_im, catfun1_rm, catfun1_cm, catfun1_dm, catfun1_zm
  END INTERFACE catfun1 
  INTERFACE catfun2
    module procedure &
      catfun2_iv, catfun2_rv, catfun2_cv, catfun2_dv, catfun2_zv, &
      catfun2_im, catfun2_rm, catfun2_cm, catfun2_dm, catfun2_zm
  END INTERFACE catfun2
  
  INTEGER, PARAMETER :: Nprimes   = 1000
  INTEGER, DIMENSION(Nprimes) :: primes ! 1st 1000 prime numbers, in 100 lines X 10 numbers/line
  LOGICAL  :: primes_read_already = .false.
  PRIVATE  :: Nprimes, primes, primes_read_already
CONTAINS
!==================================

real function percentile2D(p, a) result(val)
! p-percentile of a(:,:)   0<=p<=100
real,     intent(in)   :: p, a(:,:)
integer   :: k, j, M, km,kn
real      :: ak
k   = max(1, floor(p/100*size(a) + 0.5))
if (p>=100) then
  val = maxval(a)
else if (k==1) then
  val = minval(a)
else
  M = size(a,1)
  kn= (k-1)/M    ! 0-based index
  km=  k - M *kn
  kn = kn + 1    ! 1-based
  ak=a(km,kn)
  
  j  = count(a < ak)  ! # how many a(:) < ak
  ! write(*,*) 'Q2D:   j, k, ak', j, k, ak
  if (j>=k) then
    val = selec_rank1D(k, pack(a, a < ak))
  else
    j   = k - (size(a) - count(a > ak))
    if(j>0) then
      val = selec_rank1D(j, pack(a, a > ak))
    else
      val = ak
    end if
  end if
end if
end function

real function percentile1D(p, a) result(val)
! p-percentile of a(:)   0<=p<=100
real,     intent(in)   :: p, a(:)
integer   :: k
k   = max(1, floor(p/100*size(a) + 0.5))
if (p>=100) then
  val = maxval(a)
else if (k==1) then
  val = minval(a)
else
  val = selec_rank1D(k, a)    
end if
end function

recursive real function selec_rank1D(k, a) result(val)
! if in a sorted array, val = a(k)
integer,  intent(in)   :: k ! # position in array
real,     intent(in)   :: a(:)
integer   :: j
real      :: ak
ak = a(k)
j  = count(a < ak)  ! # how many a(:) < ak
! write(*,*) 'Q1D:   j, k, ak', j, k, ak
if (j>=k) then
  val = selec_rank1D(k, pack(a, a < ak))
else
  j   = k - (size(a) - count(a > ak))
  if(j>0) then
    val = selec_rank1D(j, pack(a, a > ak))
  else
    val = ak
  end if
end if
end function
!!!!!!!!!!!!

subroutine laplace2_ext(B, A)  ! perform laplacian, then the boundary values are extended out from the inner
real,      intent(out) :: B(:,:)
real,      intent(in)  :: A(:,:)  ! same size
integer    :: nx, nz, ix, iz
nx = size(A,1)
nz = size(A,2)
!$OMP parallel do PRIVATE(ix)   
do iz = 2, nz-1
  do ix = 2, nx-1
    B(ix,iz) = A(ix-1,iz) + A(ix+1,iz) + A(ix,iz-1) + A(ix,iz+1) - 4.*A(ix,iz)
  end do
  B(1, iz) = B(2,   iz)
  B(nx,iz) = B(nx-1,iz)
end do
B(:,1) = B(:,   2)
B(:,nz)= B(:,nz-1)
end subroutine
! . . . . . . . 

function diff_rv(x) result(d)
real, intent(in) :: x(:)
real  :: d(size(x)-1)
d = x(2:) - x(:size(x)-1)
end function

function diff0_rv(x) result(d)
real, intent(in) :: x(:)
real  :: d(size(x))
d(1)  = x(1)
d(2:) = x(2:) - x(:size(x)-1)
end function

function diff0_rm(x, wei) result(d)
real,    intent(in) :: x(:,:)
integer, intent(in) :: wei   ! dimension, weidu
real  :: d(size(x,1), size(x,2))
integer  :: k, M, N
M = size(x,1);    N = size(x,2)
if (wei==1) then
  do k = 1, N
    d(1,k)   = x(1,k)
    d(2:, k) = x(2:, k) - x(:M-1, k)
  end do
else
  d(:,1) = x(:,1)
  do k = 2, N
    d(:,k) = x(:,k) - x(:,k-1)
  end do
end if

end function

subroutine cumsum2d(B,A,kdim)
! B <--cumulative sum A, along dimension kdim
real, intent(out) :: B(:,:)  ! same size as A, space pre-assigned
real, intent(in)  :: A(:,:)
integer, optional, intent(in) :: kdim ! default 1
integer  :: d, j,k, shp(2)

d  = 1;  if(present(kdim))  d = kdim
shp= shape(A)

if (d==1) then
  do k = 1, shp(2)
    B(1,k) = A(1,k)
    do j=2, shp(1)
      B(j,k) = B(j-1,k) + A(j,k)
    end do
  end do
else
 B(:,1) = A(:,1)
 do k = 2, shp(2)
   B(:,k) = B(:,k-1) + A(:,k)
 end do
end if
end subroutine
! ******* ******* ******* SMOOTHING ******* ******* ******* *****
subroutine rect1D_mirror_conv_D(A, il, ir, B)  ! rect(x \in [-il,ir]) (*) A --> B. 
!    Note: A is modified as cumsum(A)
real(8), intent(inout) :: A(:)   ! A is Double
real(8), intent(out)   :: B(:)   ! same size as A
integer, intent(in)    :: il, ir ! index displacement.  I.e, a rect of width 4, center is regarded as index 2.  Then il = -1, ir = 2
! if -il == ir, then rect width is odd; otherwise, even.
integer  :: i, na, sml, smr      ! sm = sum of algebraic index and mirror-reflected index
na   = size(A)

! handle the boundaries first:
! See the  MRI, MRH, ... in Bspline1Dconv()
! left  boundary MRH, if (rect is MRI) | (rect is MRH & -il > ir)
sml = 2    ! MRI,  by default
if (-il >= ir)  sml = 1  ! MRH
! right boundary MRH, if (rect is MRI) | (rect is MRH & -il < ir)
smr = 2*na ! MRI,  at the right boundary, by default
if (ir >= -il)  smr = smr + 1    ! MRH

do i = 1, -il
  B(i) = sum(A(sml:sml-(i+il))) + sum(A(1:i+ir))
end do
do i = na-ir+1, na
  B(i) = sum(A(smr-(i+ir):smr-(na+1))) + sum(A(i+il:na))
end do

! do the valid (inner) part
do i = 2, na      ! cumsum, in-place
  A(i) = A(i) + A(i-1)
end do
B(1-il)= A(1-il+ir)
do i   = 2-il, na-ir
  B(i) = A(i+ir) - A(i+il-1)
end do
end subroutine

subroutine Bspline1D_mirr_conv_S(A, W, tmp1, tmp2)
! A <-- A (*) rect's, of width W(1), W(2),...,W(n).  A  single precision.
! If conv w/ a single rect, then write W as (/ w /), in an array form.
! tmp1,2:   same length as A,  double precision, to avoid error accumulation
! bundary:  mirror reflection
! Lemma:    if A is mirror-symmetric (MR), then A (*) rect is also MR.
! Two cases of MR:  
! 1)  MRH:  MR w.r.t. half integer
! 2)  MRI:  MR w.r.t. integer
! If W(i) is even, then the rect is MRH; otherwise, MRI
! Lemma:    MRH (*) MRI = MRH;   MRI (*) MRI = MRI;  MRH (*) MRH = MRI
! If  rect is MRI, then let A be MRH, so the output length stays the same
! If  rect is MRH, then let A be MRI, so the output length = length(A) - 1  [we can extend 1 pixel easily]
!             or set A's left- and right-boundary to be MRI & MRH.  See rect1D_mirror_conv_D().
real,    intent(inout) :: A(:)
integer, intent(in)    :: W(:)
real(8), intent(inout),target :: tmp1(:), tmp2(:)
real(8), pointer       :: p1(:), p2(:)
integer  :: wi, i, c, il, ir, cenR  ! when rect of even width, take center+0.5 as location 0
p1 => tmp1; p2 => tmp2
p1 = A   ! a copying
cenR = 0 ! init.  when wi is even, then flip-flop
do i = 1, size(W)
  wi = W(i)
  if (wi==1)   cycle !!!!!!!!! skip triviality
  c  = ceiling(wi/2.)
  if (c+c==wi) then  ! even
    if (cenR==1) c = c + 1
    cenR = 1 - cenR  ! flip-flop
  end if
  il = 1 - c;  ir = wi - c
  call rect1D_mirror_conv_D(p1, il, ir, p2)
  call swap_pointers(p1, p2)
end do
A = p1  ! copying back
end subroutine

subroutine Bspline1D_mirr_conv_D(A, W, tmp)  ! like Bspline1D_mirr_conv_S, except A is double
real(8), intent(inout),target :: A(:), tmp(:)
integer, intent(in)    :: W(:)
real(8), pointer       :: p1(:), p2(:)
integer  :: nw, wi, i, c, il, ir, cenR  ! when rect of even width, take center+0.5 as location 0
p1 => A; p2 => tmp
cenR = 0  ! init.  when wi is even, then flip-flop
nw   = size(W)
do i = 1, nw
  wi = W(i)
  c  = ceiling(wi/2.)
  if (c+c==wi) then  ! even
    if (cenR==1) c = c + 1
    cenR = 1 - cenR  ! flip-flop
  end if
  il = 1 - c;  ir = wi - c
  call rect1D_mirror_conv_D(p1, il, ir, p2)
  call swap_pointers(p1, p2)
end do
if ((nw/2)*2 /= nw) A = p1   ! need to copy back to A if nw is odd.
end subroutine
 
subroutine Bspline2D_mirr_conv_S(A, W) 
! see the comment in Bspline1D_mirr_conv_S()
real,    intent(inout) :: A(:,:)
integer, intent(in)    :: W(:)
integer  :: shp(2), j
real(8), allocatable   :: tmp1(:), tmp2(:)
shp = shape(A)

!$omp parallel private(j, TMP1, TMP2)
allocate(tmp1(shp(1)), tmp2(shp(1)) )
!$omp do
do j = 1, shp(2)
  tmp1  = A(:,j)
  call Bspline1D_mirr_conv_D(tmp1, W, tmp2)
  A(:,j)= tmp1
end do
!$omp end do
deallocate(tmp1, tmp2)

allocate(tmp1(shp(2)), tmp2(shp(2)) )
!$omp do
do j= 1, shp(1)
  tmp1  = A(j,:)
  call Bspline1D_mirr_conv_D(tmp1, W, tmp2)
  A(j,:)= tmp1
end do
!$omp end do
deallocate(tmp1, tmp2)
!$omp end parallel
end subroutine

subroutine Bspline3D_mirr_conv_S(A, W)
real,    intent(inout) :: A(:,:,:)
integer, intent(in)    :: W(:)
integer  :: shp(3), j, k, M
real(8), allocatable   :: tmp1(:), tmp2(:)
shp = shape(A)
M   = maxval(shp)

!$omp parallel private(j, k, TMP1, TMP2)  !===
allocate(tmp1(M), tmp2(M))

!$omp do
do k= 1, shp(3)
  do j= 1, shp(2)
    tmp1(1:shp(1))  = A(:,j,k)
    call Bspline1D_mirr_conv_D( tmp1(1:shp(1)), W, tmp2(1:shp(1)) )
    A(:,j,k)  = tmp1(1:shp(1))
  enddo
  do j= 1, shp(1)
    tmp1(1:shp(2))  = A(j,:,k)
    call Bspline1D_mirr_conv_D(tmp1(1:shp(2)), W, tmp2(1:shp(2)) )
    A(j,:,k)  = tmp1(1:shp(2))
  end do  
end do
!$omp end do

!$omp do
do k = 1, shp(2)
  do j= 1, shp(1)
    tmp1(1:shp(3))  = A(j,k,:)
    call Bspline1D_mirr_conv_D(tmp1(1:shp(3)), W, tmp2(1:shp(3)) )
    A(j,k,:)  = tmp1(1:shp(3))
  end do
end do    
!$omp end do
  
deallocate(tmp1, tmp2)
!$omp end parallel      !===
end subroutine
!******************************** PRIME # ************************
subroutine read_prime_numbers
 integer :: i, beg
 open(10, file='prime_numbers_1000.txt', err=20, form='formatted', status='old')
 do i = 1, 100
   beg = (i-1)*10 + 1
   read(10, *) primes(beg:beg+9)
 enddo
 close(10)
 return
20 stop '  Read error: prime_numbers_1000.txt!'
end subroutine
!-------------
subroutine PRIMEFACTORS(num, factors, nf)
  INTEGER, INTENT(IN) ::num  ! input number
  INTEGER, INTENT(OUT)::nf,  factors(:) ! # of factors, array to store factors
  ! e.g., factors = {2,2,7,7}, nf = 4
  integer  :: ip, n, tmp     ! ip: index to prime #,  n: store 'num'
  if (.not. primes_read_already) then
    call read_prime_numbers
    primes_read_already = .true.
  endif
  ip = 1
  nf = 0
  n  = num
  do while(n > 1)
    if (ip > Nprimes)  stop 'needed prime factor > 7919, too large!'
    tmp = primes(ip)
    if ( (n/tmp)*tmp == n ) then ! a new factor
      nf = nf + 1;
      factors(nf) = tmp
      n  = n/tmp
    else
      ip = ip+1
    endif
  enddo
end subroutine
!-------------
function If_Has_Small_Prime_Factors(num, nf) result(whe) 
  ! num:  input number in question
  ! nf:   the smallest nf prime factors
  ! whe:  whether num consists of only the nf smallest prime factors?
  integer, intent(in)  :: num, nf
  logical  :: whe
  integer  :: ip, n, tmp
  if (.not. primes_read_already) then
    call read_prime_numbers
    primes_read_already = .true.
  endif
  n = num
  do ip = 1, nf
    tmp = primes(ip)
    do while ((n/tmp)*tmp == n)
      n = n/tmp
    end do
  end do
  whe = n==1
end function
!---------------------------------
function findBoundingIntv(s, y) result(md)  !!!!!!!!!!!!!!!!!!!!!!!!!!
! s is a non-decreasing series
! find 'md' such that s(md) < y <= s(md+1), using 2-beta search
integer, intent(in)  :: s(:), y
integer  :: md, bg, en
bg  = 1;    en = size(s)
if (y <= s(1)) then
  md = 0     ! as out-of-bound signal
  return
else if (y > s(en)) then
  md = en    ! as out-of-bound signal
  return
endif
DO
  if (en == bg+1)  then
    md = bg
    return
  endif
  md = (bg + en)/2
  if (s(md+1) < y) then
    bg = md + 1   ! take right half
  else if (y <= s(md)) then
    en = md       ! take left  half
  else             
    return
  endif
END DO
end function
!-------------
subroutine binnedIncrSeriesIdx(s, y1, y2, idx, num)
! s is a non-decreasing series
! find all the idx satisfying y1 <= s(idx) < y2,  num = length(idx)
!!! i.e., y2 is non-inclusive 
integer, intent(in)  :: s(:), y1, y2
integer, intent(OUT) :: idx(:), num   ! make sure idx has enough space
integer  :: b, a, k
! 1.  find a, such that s(a)<y1 and s(a+1)>=y1
a   = findBoundingIntv(s, y1)
! 2.  find b, such that s(b)<y2 and s(b+1)>=y2
b   = findBoundingIntv(s, y2)
idx = (/ (k, k=a+1, b) /)
num = b-a   ! # elements inside [y1, y2)
end subroutine
!------
subroutine binnedSeriesNums(s, wintv, nss, cmsm0)
! s is a series
! intervals are [1, wintv] [1+wintv, 2*wintv] [1+2*wintv, 3*wintv], ...
! find how many series elements fall in each interval --> nss
! cumulative sum of 'nss' starting @ 0                --> cmsm0.  
!   is the displacement of each interval's stored s's relative to the beginning
integer, intent(in)  :: s(:), wintv
integer, intent(OUT) :: nss(:)       ! make sure enough length
integer, intent(OUT),   optional ::  cmsm0(:)  
integer  :: k, i, L
L    = size(s)
nss  = 0  ! init.
do k = 1, L
  i  = (s(k)-1)/wintv + 1  ! find interval index
  nss(i) = nss(i) + 1
enddo
if (present(cmsm0)) then
  L = size(nss)
  cmsm0(1) = 0
  if (L < 2) return
  cmsm0(2) = nss(1)
  do k = 3, L
    cmsm0(k) = cmsm0(k-1) + nss(k-1)
  enddo
endif
end subroutine
!---------------------
subroutine stripSmallPrimeFactors(A, C, B)
! A = C B, where C = 2^m 3^n 5^o 7^p, B is a prime number > 7
integer, intent(in)  :: A
integer, intent(out) :: C, B
integer  :: prm(4) = (/2,3,5,7/), i, p
C = 1;  B = A
outer: do i = 1, 4  ! go through the 4 smallest prime numbers
  p = prm(i)
  inner: do
    if ( (B/p)*p == B ) then
      B = B/p;  C = C*p
    else
      exit inner
    endif
  enddo inner 
enddo outer
end subroutine

function roundUp2SmallPrimeProduct(A) result (C)
! round A up to a smallest integer C that is 2^m 3^n 5^o 7^p
integer, intent(in) :: A
integer  :: C
C = A
do while(.not.If_Has_Small_Prime_Factors(C, 4))
  C = C + 1
enddo
end function

!*****************************************
subroutine find_indices(ind, nfound, crit)
  integer, intent(out) :: ind(:),  nfound   ! nfound limited by size(ind)
  logical, intent(in)  :: crit(:)
  
  integer  :: nind, k
  nind = size(ind)
  nfound = 0
  do k = 1, size(crit)  ! prefer head
    if(.not. crit(k)) cycle
    nfound = nfound + 1
    ind(nfound) = k
    if (nfound==nind) exit
  end do
end subroutine

FUNCTION SUB2IND(siz, sub) result (ind)  !!!!!!!!!!!!!!!!!!!!!!!!!!
! siz    D x 1   array of the shape of multi-dim tensor
! sub    D x 1   subscript
! return scalar index
! ***using FORTRAN array convention***
integer, intent(in) :: siz(:), sub(:)
integer  :: ind, D, k
D   = size(siz)
ind = sub(D)-1
do k = D-1, 1 ,-1
  ind = ind*siz(k) + sub(k)-1
enddo
ind = ind + 1
END FUNCTION

!===========
function ind2sub_s(ndx, siz)
! scalar ndx
! siz    D x 1   array of the shape of multi-dim tensor
! ndx    1 x 1   index into the tensor, treated as 1d vector (1st dim is the fastest)
! return D x 1   subscript vector
! ***using FORTRAN array convention***
!    following MATLAB's naming.
integer, intent(in) :: siz(:), ndx
integer             :: ind2sub_s(size(siz))
integer             :: D, k, i, n
D    = size(siz)
ind2sub_s(1) = 1
ind2sub_s(2) = siz(1)
do k = 3, D
  ind2sub_s(k) = siz(k-1) * ind2sub_s(k-1)
enddo  ! [1, cumprod]
i    = ndx
do k = D, 2, -1
  n  = (i-1)/ind2sub_s(k)  ! C convention
  i  = i - n*ind2sub_s(k)  
  ind2sub_s(k) = n + 1     ! to Fortran convention
enddo
ind2sub_s(1)   = i   ! @ k = 1
end function
!-----

function ind2sub_v(ndx, siz)
! vector ndx
! siz    D x 1   array of the shape of multi-dim tensor
! ndx    N x 1   N indices into the same sized tensor, treated as 1d vector (1st dim is the fastest)
! return D x N   subscript matrix.  Each column is a D-dimensional subscript vector, and 
!        there're N of them.
! ***using FORTRAN array convention***
!    following MATLAB's naming.
integer, intent(in) :: siz(:), ndx(:)
integer             :: ind2sub_v(size(siz), size(ndx))
integer             :: D, N, kd, kn, i, x

D   = size(siz);    N = size(ndx)
ind2sub_v(1,1) = 1
ind2sub_v(2,1) = siz(1)
do kd = 3, D
  ind2sub_v(kd,1) = siz(kd-1) * ind2sub_v(kd-1,1)
enddo  ! [1; cumprod] saved @ the 1st column of ind2sub_v
do kn = N, 1, -1
  i   = ndx(kn)
  do kd = D, 2, -1
    x   = (i-1)/ind2sub_v(kd,1)  ! C convention
    i   = i - x*ind2sub_v(kd,1) 
    ind2sub_v(kd,kn) = x + 1     ! to Fortran convention
  enddo
  ind2sub_v(1, kn)   = i   ! @ kd = 1
enddo
end function




! * * * * * * * * * * * * * * * * * *
!           1-Norm
! * * * * * * * * * * * * * * * * * *
function norm_1_rv(v)
REAL(SP) :: v(:), norm_1_rv
integer  :: i
norm_1_rv = 0.
do i = 1, size(v)
  if (v(i)>0.) then
    norm_1_rv = norm_1_rv + v(i)
  else
    norm_1_rv = norm_1_rv - v(i)
  endif
end do
end function

function norm_1_dv(v)  result(nrm)
REAL(DP) :: v(:), nrm
integer  :: i
nrm = 0.
do i = 1, size(v)
  if (v(i)>0.) then
    nrm = nrm + v(i)
  else
    nrm = nrm - v(i)
  endif
end do
end function

function norm_1_rm(v)  result(nrm)  ! matrix v
REAL(SP) :: v(:,:), nrm
REAL(DP) :: nrm_c   ! local accumulator, double prec.
integer  :: i, j
nrm_c = 0.
do i = 1, size(v,2)
  do j = 1, size(v,1)
  if (v(j,i)>0.) then
    nrm_c = nrm_c + v(j,i)
  else
    nrm_c = nrm_c - v(j,i)
  end if
  end do
end do
nrm = nrm_c
end function

function norm_1_dm(v)  result(nrm)  ! matrix v
REAL(DP) :: v(:,:), nrm
integer  :: i, j
nrm = 0.
do i = 1, size(v,2)
  do j = 1, size(v,1)
  if (v(j,i)>0.) then
    nrm = nrm + v(j,i)
  else
    nrm = nrm - v(j,i)
  endif
  end do
end do
end function

function norm_1_r3d(v)  result(nrm_o)  ! matrix v
REAL(SP) :: v(:,:,:), nrm_o
REAL(DP) :: nrm
integer  :: i, j, k
nrm = 0.
do k = 1, size(v,3)
  do i = 1, size(v,2)
    do j = 1, size(v,1)
    if (v(j,i,k)>0.) then
      nrm = nrm + v(j,i,k)
    else
      nrm = nrm - v(j,i,k)
    end if
    end do
  end do
end do
nrm_o = nrm
end function

function norm_1_d3d(v)  result(nrm)  ! matrix v
REAL(DP) :: v(:,:,:), nrm
integer  :: i, j, k
nrm = 0.
do k = 1, size(v,3)
  do i = 1, size(v,2)
    do j = 1, size(v,1)
      if (v(j,i,k)>0.) then
        nrm = nrm + v(j,i,k)
      else
        nrm = nrm - v(j,i,k)
      end if
    end do
  end do
end do  
end function

! * * * * * * * * * * * * *
!       inf-Norm
! * * * * * * * * * * * * *
function norm_inf_rv(v) result(nrm)
REAL(SP) :: v(:), nrm
integer  :: i
nrm = 0.
do i = 1, size(v)
  if (v(i)> nrm) then
    nrm  = v(i)
  else if (-v(i)> nrm) then
    nrm  = -v(i)  
  end if  
end do
end function

function norm_inf_dv(v) result(nrm)
REAL(DP) :: v(:), nrm
integer  :: i
nrm = 0.
do i = 1, size(v)
  if (v(i)> nrm) then
    nrm  = v(i)
  else if (-v(i)> nrm) then
    nrm  = -v(i)  
  end if  
end do
end function

function norm_inf_rm(v) result(nrm)
REAL(SP) :: v(:,:), nrm
integer  :: i, j
nrm = 0.
do j = 1, size(v,2)
  do i = 1, size(v,1)
    if (v(i,j)> nrm) then
      nrm  = v(i,j)
    else if (-v(i,j)> nrm) then
      nrm  = -v(i,j)  
    end if  
  end do
end do
end function

function norm_inf_dm(v) result(nrm)
REAL(DP) :: v(:,:), nrm
integer  :: i, j
nrm = 0.
do j = 1, size(v,2)
  do i = 1, size(v,1)
    if (v(i,j)> nrm) then
      nrm  = v(i,j)
    else if (-v(i,j)> nrm) then
      nrm  = -v(i,j)  
    end if  
  end do
end do
end function

function norm_inf_r3d(v) result(nrm)
REAL(SP) :: v(:,:,:), nrm
integer  :: i, j, k
nrm = 0.
do k = 1, size(v,3)
  do j = 1, size(v,2)
    do i = 1, size(v,1)
     if (v(i,j,k)> nrm) then
      nrm  = v(i,j,k)
     else if (-v(i,j,k)> nrm) then
      nrm  = -v(i,j,k)  
     end if  
    end do
  end do
end do
end function

function norm_inf_d3d(v) result(nrm)
REAL(DP) :: v(:,:,:), nrm
integer  :: i, j, k
nrm = 0.
do k = 1, size(v,3)
  do j = 1, size(v,2)
    do i = 1, size(v,1)
     if (v(i,j,k)> nrm) then
      nrm  = v(i,j,k)
     else if (-v(i,j,k)> nrm) then
      nrm  = -v(i,j,k)  
     end if  
    end do
  end do
end do
end function

! * * * * * * * * * * * * *
! square of Euclidean norm
FUNCTION norm2_rv(v)  ! real
REAL(SP)     :: v(:)
REAL(SP)     :: norm2_rv
norm2_rv = dot_product(v,v)
END FUNCTION

FUNCTION norm2_dv(v)  ! double
REAL(DP)     :: v(:)
REAL(DP)     :: norm2_dv
norm2_dv = dot_product(v,v)
END FUNCTION

FUNCTION norm2_cv(v)  ! complex
COMPLEX(SPC) :: v(:)
REAL(SP)     :: norm2_cv
norm2_cv = dot_product(v,v)
END FUNCTION

FUNCTION norm2_zv(v)  ! complex double
COMPLEX(DPC) :: v(:)
REAL(DP)     :: norm2_zv
norm2_zv = dot_product(v,v)
END FUNCTION
!-------------------
FUNCTION norm2_rm(m)  result(norm2_o)  !! matrix
REAL(SP)     :: m(:,:), norm2_o
REAL(DP)     :: norm2 !!! double
integer      :: i,j, shp(2)
shp   = shape(m)
norm2 = 0.
do j  = 1, shp(2)
  do i = 1, shp(1)
    norm2 = norm2 + DBLE(m(i,j))**2
  end do
end do
norm2_o = norm2
END FUNCTION

FUNCTION norm2_dm(m)  result(norm2)  
REAL(DP)     :: m(:,:)
REAL(DP)     :: norm2
integer      :: i, j, shp(2)
shp   = shape(m)
norm2 = 0.
do j  = 1, shp(2)
  do i = 1, shp(1)
    norm2 = norm2 + m(i,j)**2
  enddo
enddo
END FUNCTION

FUNCTION norm2_cm(m)  result(norm2_o) 
COMPLEX(SPC) :: m(:,:)
REAL(SP)     :: norm2_o
REAL(DP)     :: norm2  !!! double
integer      :: i, j, shp(2)
shp   = shape(m)
norm2 = 0.
do j  = 1, shp(2)
  do i = 1, shp(1)
    norm2 = norm2 + DBLE(m(i,j))*DBLE(conjg(m(i,j)))
  enddo
enddo
norm2_o = norm2
END FUNCTION

FUNCTION norm2_zm(m)  result(norm2) 
COMPLEX(DPC) :: m(:,:)
REAL(DP)     :: norm2
integer      :: i, j, shp(2)
shp   = shape(m)
norm2 = 0
do j  = 1, shp(2)
  do i = 1, shp(1)
    norm2 = norm2 + m(i,j)*conjg(m(i,j))
  enddo
enddo 
END FUNCTION
!-------------------
FUNCTION norm2_r3d(t)    result(norm2_o)  !!! 3d tensor
REAL(SP)     :: t(:,:,:), norm2_o
REAL(DP)     :: norm2    !!! double
integer      :: i, j, k, shp(3)
shp   = shape(t)
norm2 = 0.
do k = 1, shp(3)
  do j = 1, shp(2)
    do i = 1, shp(1)
      norm2 = norm2 + DBLE(t(i,j,k))**2
    enddo
  enddo
enddo
norm2_o = norm2
END FUNCTION

FUNCTION norm2_d3d(t)   result(norm2) 
REAL(DP)     :: t(:,:,:)
REAL(DP)     :: norm2
integer      :: i, j, k, shp(3)
shp   = shape(t)
norm2 = 0
do k = 1, shp(3)
  do j = 1, shp(2)
    do i = 1, shp(1)
      norm2 = norm2 + t(i,j,k)**2
    enddo
  enddo
enddo
END FUNCTION

FUNCTION norm2_c3d(t)   result(norm2_o) 
COMPLEX(SPC) :: t(:,:,:)
REAL(SP)     :: norm2_o
REAL(DP)     :: norm2   !!! double
integer      :: i, j, k, shp(3)
shp   = shape(t)
norm2 = 0.
do k = 1, shp(3)
  do j = 1, shp(2)
    do i = 1, shp(1)
      norm2 = norm2 + DBLE(t(i,j,k))*DBLE(conjg(t(i,j,k)))
    enddo
  enddo
enddo
norm2_o = norm2
END FUNCTION

FUNCTION norm2_z3d(t)  result(norm2) 
COMPLEX(DPC) :: t(:,:,:)
REAL(DP)     :: norm2
integer      :: i, j, k, shp(3)
shp   = shape(t)
norm2 = 0
do k = 1, shp(3)
  do j = 1, shp(2)
    do i = 1, shp(1)
      norm2 = norm2 + t(i,j,k)*conjg(t(i,j,k))
    enddo
  enddo
enddo
END FUNCTION

! * * * * * *  Euclidean norm * * * * *
FUNCTION norm_rv(v)  result(norm) ! real
REAL(SP)     :: v(:)
REAL(SP)     :: norm
norm = sqrt(norm2(v))
END FUNCTION

FUNCTION norm_dv(v)  result(norm) ! double
REAL(DP)     :: v(:)
REAL(DP)     :: norm
norm = sqrt(norm2(v))
END FUNCTION

FUNCTION norm_cv(v)  result(norm)  ! complex
COMPLEX(SPC) :: v(:)
REAL(SP)     :: norm
norm = sqrt(norm2(v))
END FUNCTION

FUNCTION norm_zv(v)  result(norm)  ! complex double
COMPLEX(DPC) :: v(:)
REAL(DP)     :: norm
norm = sqrt(norm2(v))
END FUNCTION
!-------------------
FUNCTION norm_rm(m)  result(norm)  !! matrix
REAL(SP)     :: m(:,:)
REAL(SP)     :: norm
norm = sqrt(norm2(m))
END FUNCTION

FUNCTION norm_dm(m)  result(norm)
REAL(DP)     :: m(:,:)
REAL(DP)     :: norm
norm = sqrt(norm2(m))
END FUNCTION

FUNCTION norm_cm(m)  result(norm)
COMPLEX(SPC) :: m(:,:)
REAL(SP)     :: norm
norm = sqrt(norm2(m))
END FUNCTION

FUNCTION norm_zm(m)  result(norm)
COMPLEX(DPC) :: m(:,:)
REAL(DP)     :: norm
norm = sqrt(norm2(m))
END FUNCTION
!-------------------
FUNCTION norm_r3d(t)    result(norm)  !!! 3d tensor
REAL(SP)     :: t(:,:,:)
REAL(SP)     :: norm
norm = sqrt(norm2(t))
END FUNCTION

FUNCTION norm_d3d(t)   result(norm) 
REAL(DP)     :: t(:,:,:)
REAL(DP)     :: norm
norm = sqrt(norm2(t))
END FUNCTION

FUNCTION norm_c3d(t)   result(norm) 
COMPLEX(SPC) :: t(:,:,:)
REAL(SP)     :: norm
norm = sqrt(norm2(t))
END FUNCTION

FUNCTION norm_z3d(t)  result(norm) 
COMPLEX(DPC) :: t(:,:,:)
REAL(DP)     :: norm
norm = sqrt(norm2(t))
END FUNCTION

! * * * * * * * * * * * * *
! * * *  dot product  * * *

FUNCTION dot_prod_rv(u,v) result(dtp) ! real
REAL(SP)     :: u(:), v(:)
REAL(SP)     :: dtp
dtp = dot_product(u,v)
END FUNCTION

FUNCTION dot_prod_dv(u,v) result(dtp)  ! double
REAL(DP)     :: u(:), v(:)
REAL(DP)     :: dtp
dtp = dot_product(u,v)
END FUNCTION

FUNCTION dot_prod_cv(u,v) result(dtp)  ! complex
COMPLEX(SPC) :: u(:), v(:)
REAL(SP)     :: dtp
dtp = dot_product(u,v)
END FUNCTION

FUNCTION dot_prod_zv(u,v) result(dtp)  ! complex double
COMPLEX(DPC) :: u(:), v(:)
REAL(DP)     :: dtp
dtp = dot_product(u,v)
END FUNCTION
!-------------------
FUNCTION dot_prod_rm(m,n)  result(dtp_o)  !! matrix
REAL(SP)     :: m(:,:), n(:,:), dtp_o
REAL(DP)     :: dtp   !!! double
integer      :: i, j, shp(2)
shp   = shape(m)
dtp   = 0.
do j  = 1, shp(2)
  do i = 1, shp(1)
    dtp = dtp + DBLE(m(i,j)) * DBLE(n(i,j))
  enddo
enddo
dtp_o = dtp
END FUNCTION

FUNCTION dot_prod_dm(m,n)  result(dtp)  
REAL(DP)     :: m(:,:), n(:,:), dtp
integer      :: i, j, shp(2)
shp   = shape(m)
dtp   = 0.
do j  = 1, shp(2)
  do i = 1, shp(1)
    dtp = dtp + m(i,j) * n(i,j) 
  enddo
enddo
END FUNCTION

FUNCTION dot_prod_cm(m,n)  result(dtp_o) 
COMPLEX(SPC) :: m(:,:), n(:,:), dtp_o
COMPLEX(DPC) :: dtp  !!! double
integer      :: i, j, shp(2)
shp   = shape(m)
dtp   = 0.
do j  = 1, shp(2)
  do i = 1, shp(1)
    dtp = dtp + n(i,j)*conjg(m(i,j))
  enddo
enddo
dtp_o = dtp
END FUNCTION

FUNCTION dot_prod_zm(m,n)  result(dtp) 
COMPLEX(DPC) :: m(:,:),n(:,:), dtp
integer      :: i, j, shp(2)
shp   = shape(m)
dtp   = 0.
do j  = 1, shp(2)
  do i = 1, shp(1)
    dtp = dtp + n(i,j)*conjg(m(i,j))
  enddo
enddo 
END FUNCTION
!-------------------
FUNCTION dot_prod_r3d(s,t)    result(dtp_o)  !!! 3d tensor
REAL(SP)     :: s(:,:,:), t(:,:,:), dtp_o
REAL(DP)     :: dtp      !!! double
integer      :: i, j, k, shp(3)
shp  = shape(t)
dtp  = 0.
do k = 1, shp(3)
  do j = 1, shp(2)
    do i = 1, shp(1)
      dtp = dtp + DBLE(t(i,j,k)) * DBLE(s(i,j,k))
    enddo
  enddo
enddo
dtp_o = dtp
END FUNCTION

FUNCTION dot_prod_d3d(s,t)   result(dtp) 
REAL(DP)     :: s(:,:,:), t(:,:,:), dtp
integer      :: i, j, k, shp(3)
shp  = shape(t)
dtp   = 0
do k = 1, shp(3)
  do j = 1, shp(2)
    do i = 1, shp(1)
      dtp = dtp + t(i,j,k) * s(i,j,k) 
    enddo
  enddo
enddo
END FUNCTION

FUNCTION dot_prod_c3d(s,t)   result(dtp_o) 
COMPLEX(SPC) :: s(:,:,:), t(:,:,:), dtp_o
COMPLEX(DPC) :: dtp      !!! double
integer      :: i, j, k, shp(3)
shp  = shape(t)
dtp  = 0.
do k = 1, shp(3)
  do j = 1, shp(2)
    do i = 1, shp(1)
      dtp = dtp + t(i,j,k)*conjg(s(i,j,k))
    enddo
  enddo
enddo
dtp_o = dtp
END FUNCTION

FUNCTION dot_prod_z3d(s,t)  result(dtp) 
COMPLEX(DPC) :: s(:,:,:), t(:,:,:), dtp
integer      :: i, j, k, shp(3)
shp  = shape(t)
dtp   = 0
do k = 1, shp(3)
  do j = 1, shp(2)
    do i = 1, shp(1)
      dtp = dtp + t(i,j,k)*conjg(s(i,j,k))
    enddo
  enddo
enddo
END FUNCTION

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SUBROUTINE swap_pointers_1d(pa, pb)
REAL(DP), pointer  :: pa(:), pb(:), tmp(:)
tmp => pa;  pa => pb;  pb => tmp
END SUBROUTINE

SUBROUTINE swap_pointers_1r(pa, pb)
REAL(SP), pointer  :: pa(:), pb(:), tmp(:)
tmp => pa;  pa => pb;  pb => tmp
END SUBROUTINE

SUBROUTINE swap_pointers_1i(pa, pb)
integer,  pointer  :: pa(:), pb(:), tmp(:)
tmp => pa;  pa => pb;  pb => tmp
END SUBROUTINE

SUBROUTINE swap_pointers_2r(pa, pb)
REAL(SP), pointer  :: pa(:,:), pb(:,:), tmp(:,:)
tmp => pa;  pa => pb;  pb => tmp
END SUBROUTINE

SUBROUTINE swap_pointers_2i(pa, pb)
integer,  pointer  :: pa(:,:), pb(:,:), tmp(:,:)
tmp => pa;  pa => pb;  pb => tmp
END SUBROUTINE

SUBROUTINE swap_pointers_3r(pa, pb)
REAL(SP), pointer  :: pa(:,:,:), pb(:,:,:), tmp(:,:,:)
tmp => pa;  pa => pb;  pb => tmp
END SUBROUTINE

SUBROUTINE swap_pointers_3i(pa, pb)
integer,  pointer  :: pa(:,:,:), pb(:,:,:), tmp(:,:,:)
tmp => pa;  pa => pb;  pb => tmp
END SUBROUTINE
!~~~~~~~~~~~~~~~~~~ swapDimInplace ~~~~~~~~~~~~~~~~~~
SUBROUTINE swapDimInplace2r(A)
REAL(SP), intent(inout)  :: A(:,:)
REAL(SP)                 :: T
INTEGER(I4B) :: I,J, S
S = SIZE(A,1)
CALL ASSERT(S==SIZE(A,2), 'INPUT MUST BE SQUARE MATRIX TO ALLOW INPLACE TRANSPOSE')
DO I = 1,S-1
  DO J = I+1,S
    T = A(J,I);  A(J,I) = A(I,J);  A(I,J) = T;
  ENDDO
ENDDO
END SUBROUTINE
SUBROUTINE swapDimInplace2c(A)
COMPLEX(SPC), intent(inout)  :: A(:,:)
COMPLEX(SPC)                 :: T
INTEGER(I4B) :: I,J, S
S = SIZE(A,1)
CALL ASSERT(S==SIZE(A,2), 'INPUT MUST BE SQUARE MATRIX TO ALLOW INPLACE TRANSPOSE')
DO I = 1,S-1
  DO J = I+1,S
    T = A(J,I);  A(J,I) = A(I,J);  A(I,J) = T;
  ENDDO
ENDDO
END SUBROUTINE
SUBROUTINE swapDimInplace2i(A)
INTEGER(I4B), intent(inout)   :: A(:,:)
INTEGER(I4B)  :: T, I,J, S
S = SIZE(A,1)
CALL ASSERT(S==SIZE(A,2), 'INPUT MUST BE SQUARE MATRIX TO ALLOW INPLACE TRANSPOSE')
DO I = 1,S-1
  DO J = I+1,S
    T = A(J,I);  A(J,I) = A(I,J);  A(I,J) = T;
  ENDDO
ENDDO
END SUBROUTINE
SUBROUTINE swapDimInplace2d(A)
REAL(DP), intent(inout)  :: A(:,:)
REAL(DP)  :: T
INTEGER(I4B) :: I,J, S
S = SIZE(A,1)
CALL ASSERT(S==SIZE(A,2), 'INPUT MUST BE SQUARE MATRIX TO ALLOW INPLACE TRANSPOSE')
DO I = 1,S-1
  DO J = I+1,S
    T = A(J,I);  A(J,I) = A(I,J);  A(I,J) = T;
  ENDDO
ENDDO
END SUBROUTINE
SUBROUTINE swapDimInplace2z(A)
COMPLEX(DPC), intent(inout) :: A(:,:)
COMPLEX(DPC) :: T
INTEGER(I4B) :: I,J, S
S = SIZE(A,1)
CALL ASSERT(S==SIZE(A,2), 'INPUT MUST BE SQUARE MATRIX TO ALLOW INPLACE TRANSPOSE')
DO I = 1,S-1
  DO J = I+1,S
    T = A(J,I);  A(J,I) = A(I,J);  A(I,J) = T;
  ENDDO
ENDDO
END SUBROUTINE

SUBROUTINE swapDimInplace3r(A, d1, d2)  ! A: 3D array, swap it's d1 & d2 diemsions in-place
REAL(SP), intent(inout)   :: A(:,:,:)
REAL(SP)  :: T
INTEGER(I4B) :: d1,d2, I,J,K, S(3)
CALL ASSERT(d1/=d2, d1<4, d2<4, 'Dimension indices to swap must be valid (unequal, 1-->3)')
S = SHAPE(A)
IF (d1/=3 .AND. d2/=3) THEN  ! swap 1 <--> 2
  CALL ASSERT(S(1)==S(2), 'dimension to swap in-place must be of equal sizes')
  DO K = 1,S(3)
   DO I = 1,S(2)-1
    DO J = I+1,S(1)
    T = A(J,I,K);  A(J,I,K) = A(I,J,K);  A(I,J,K) = T;
    ENDDO
   ENDDO
  ENDDO
ELSE IF(d1/=2 .AND. d2/=2) THEN ! swap 1 <--> 3
  CALL ASSERT(S(1)==S(3), 'dimension to swap in-place must be of equal sizes')
  DO K = 1,S(3)-1
   DO I = 1,S(2)
    DO J = K+1,S(1)
    T = A(J,I,K);  A(J,I,K) = A(K,I,J);  A(K,I,J) = T;
    ENDDO
   ENDDO
  ENDDO
ELSE  ! swap 3 <--> 2
  CALL ASSERT(S(2)==S(3), 'dimension to swap in-place must be of equal sizes')
  DO K = 1,S(3)-1
   DO I = K+1,S(2)
    DO J = 1,S(1)
    T = A(J,I,K);  A(J,I,K) = A(J,K,I);  A(J,K,I) = T;
    ENDDO
   ENDDO
  ENDDO
END IF
END SUBROUTINE

SUBROUTINE swapDimInplace3c(A, d1, d2)  ! A: 3D array, swap it's d1 & d2 diemsions in-place
COMPLEX(SPC), intent(inout)   :: A(:,:,:)
COMPLEX(SPC)  :: T
INTEGER(I4B)  :: d1,d2, I,J,K, S(3)
CALL ASSERT(d1/=d2, d1<4, d2<4, 'Dimension indices to swap must be valid (unequal, 1-->3)')
S = SHAPE(A)
IF (d1/=3 .AND. d2/=3) THEN  ! swap 1 <--> 2
  CALL ASSERT(S(1)==S(2), 'dimension to swap in-place must be of equal sizes')
  DO K = 1,S(3)
   DO I = 1,S(2)-1
    DO J = I+1,S(1)
    T = A(J,I,K);  A(J,I,K) = A(I,J,K);  A(I,J,K) = T;
    ENDDO
   ENDDO
  ENDDO
ELSE IF(d1/=2 .AND. d2/=2) THEN ! swap 1 <--> 3
  CALL ASSERT(S(1)==S(3), 'dimension to swap in-place must be of equal sizes')
  DO K = 1,S(3)-1
   DO I = 1,S(2)
    DO J = K+1,S(1)
    T = A(J,I,K);  A(J,I,K) = A(K,I,J);  A(K,I,J) = T;
    ENDDO
   ENDDO
  ENDDO
ELSE  ! swap 3 <--> 2
  CALL ASSERT(S(2)==S(3), 'dimension to swap in-place must be of equal sizes')
  DO K = 1,S(3)-1
   DO I = K+1,S(2)
    DO J = 1,S(1)
    T = A(J,I,K);  A(J,I,K) = A(J,K,I);  A(J,K,I) = T;
    ENDDO
   ENDDO
  ENDDO
END IF
END SUBROUTINE
SUBROUTINE swapDimInplace3d(A, d1, d2)  ! A: 3D array, swap it's d1 & d2 diemsions in-place
REAL(DP), intent(inout)   :: A(:,:,:)
REAL(DP)  :: T
INTEGER(I4B) :: d1,d2, I,J,K, S(3)
CALL ASSERT(d1/=d2, d1<4, d2<4, 'Dimension indices to swap must be valid (unequal, 1-->3)')
S = SHAPE(A)
IF (d1/=3 .AND. d2/=3) THEN  ! swap 1 <--> 2
  CALL ASSERT(S(1)==S(2), 'dimension to swap in-place must be of equal sizes')
  DO K = 1,S(3)
   DO I = 1,S(2)-1
    DO J = I+1,S(1)
    T = A(J,I,K);  A(J,I,K) = A(I,J,K);  A(I,J,K) = T;
    ENDDO
   ENDDO
  ENDDO
ELSE IF(d1/=2 .AND. d2/=2) THEN ! swap 1 <--> 3
  CALL ASSERT(S(1)==S(3), 'dimension to swap in-place must be of equal sizes')
  DO K = 1,S(3)-1
   DO I = 1,S(2)
    DO J = K+1,S(1)
    T = A(J,I,K);  A(J,I,K) = A(K,I,J);  A(K,I,J) = T;
    ENDDO
   ENDDO
  ENDDO
ELSE  ! swap 3 <--> 2
  CALL ASSERT(S(2)==S(3), 'dimension to swap in-place must be of equal sizes')
  DO K = 1,S(3)-1
   DO I = K+1,S(2)
    DO J = 1,S(1)
    T = A(J,I,K);  A(J,I,K) = A(J,K,I);  A(J,K,I) = T;
    ENDDO
   ENDDO
  ENDDO
END IF
END SUBROUTINE

SUBROUTINE swapDimInplace3z(A, d1, d2)  ! A: 3D array, swap it's d1 & d2 diemsions in-place
COMPLEX(DPC), intent(inout)  :: A(:,:,:)
COMPLEX(DPC)  :: T
INTEGER(I4B)   :: d1,d2, I,J,K, S(3)
CALL ASSERT(d1/=d2, d1<4, d2<4, 'Dimension indices to swap must be valid (unequal, 1-->3)')
S = SHAPE(A)
IF (d1/=3 .AND. d2/=3) THEN  ! swap 1 <--> 2
  CALL ASSERT(S(1)==S(2), 'dimension to swap in-place must be of equal sizes')
  DO K = 1,S(3)
   DO I = 1,S(2)-1
    DO J = I+1,S(1)
    T = A(J,I,K);  A(J,I,K) = A(I,J,K);  A(I,J,K) = T;
    ENDDO
   ENDDO
  ENDDO
ELSE IF(d1/=2 .AND. d2/=2) THEN ! swap 1 <--> 3
  CALL ASSERT(S(1)==S(3), 'dimension to swap in-place must be of equal sizes')
  DO K = 1,S(3)-1
   DO I = 1,S(2)
    DO J = K+1,S(1)
    T = A(J,I,K);  A(J,I,K) = A(K,I,J);  A(K,I,J) = T;
    ENDDO
   ENDDO
  ENDDO
ELSE  ! swap 3 <--> 2
  CALL ASSERT(S(2)==S(3), 'dimension to swap in-place must be of equal sizes')
  DO K = 1,S(3)-1
   DO I = K+1,S(2)
    DO J = 1,S(1)
    T = A(J,I,K);  A(J,I,K) = A(J,K,I);  A(J,K,I) = T;
    ENDDO
   ENDDO
  ENDDO
END IF
END SUBROUTINE





!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! pppppppppppppppp  PERMUTE pppppppppppppppp
SUBROUTINE permute_r(Y, X, dims)
REAL(SP), POINTER    :: Y(:,:,:)
REAL(SP), intent(in) :: X(:,:,:)
INTEGER(I4B)      :: dims(3), shpx(3), shpy(3), indy(3), i,j,k
shpx = shape(X)
shpy = shpx(dims)
IF(ASSOCIATED(Y)) DEALLOCATE(Y)
ALLOCATE (Y(shpy(1),shpy(2),shpy(3)))
do k = 1, shpx(3)
  do j = 1, shpx(2)
    do i = 1, shpx(1)
      shpy = (/ i,j,k /) ! act as tmp
      indy = shpy(dims)
      Y(indy(1),indy(2),indy(3)) = X(i,j,k)
    enddo
  enddo
enddo
END SUBROUTINE
SUBROUTINE permute_c(Y, X, dims)
COMPLEX(SPC), POINTER    :: Y(:,:,:)
COMPLEX(SPC), intent(in) :: X(:,:,:)
INTEGER(I4B)          :: dims(3), shpx(3), shpy(3), indy(3), i,j,k
shpx = shape(X)
shpy = shpx(dims)
IF(ASSOCIATED(Y)) DEALLOCATE(Y)
ALLOCATE (Y(shpy(1),shpy(2),shpy(3)))
do k = 1, shpx(3)
  do j = 1, shpx(2)
    do i = 1, shpx(1)
      shpy = (/ i,j,k /) ! act as tmp
      indy = shpy(dims)
      Y(indy(1),indy(2),indy(3)) = X(i,j,k)
    enddo
  enddo
enddo
END SUBROUTINE
SUBROUTINE permute_d(Y, X, dims)
REAL(DP), POINTER    :: Y(:,:,:)
REAL(DP), intent(in) :: X(:,:,:)
INTEGER(I4B)      :: dims(3), shpx(3), shpy(3), indy(3), i,j,k
shpx = shape(X)
shpy = shpx(dims)
IF(ASSOCIATED(Y)) DEALLOCATE(Y)
ALLOCATE (Y(shpy(1),shpy(2),shpy(3)))
do k = 1, shpx(3)
  do j = 1, shpx(2)
    do i = 1, shpx(1)
      shpy = (/ i,j,k /) ! act as tmp
      indy = shpy(dims)
      Y(indy(1),indy(2),indy(3)) = X(i,j,k)
    enddo
  enddo
enddo
END SUBROUTINE
SUBROUTINE permute_z(Y, X, dims)
COMPLEX(DPC), POINTER    :: Y(:,:,:)
COMPLEX(DPC), intent(in) :: X(:,:,:)
INTEGER(I4B)          :: dims(3), shpx(3), shpy(3), indy(3), i,j,k
shpx = shape(X)
shpy = shpx(dims)
IF(ASSOCIATED(Y)) DEALLOCATE(Y)
ALLOCATE (Y(shpy(1),shpy(2),shpy(3)))
do k = 1, shpx(3)
  do j = 1, shpx(2)
    do i = 1, shpx(1)
      shpy = (/ i,j,k /) ! act as tmp
      indy = shpy(dims)
      Y(indy(1),indy(2),indy(3)) = X(i,j,k)
    enddo
  enddo
enddo
END SUBROUTINE





! pppppppppppppppppppppppppppppppppppppppppp

! * * * * * * * * * * * CAT * * * * * * * * * * *

! Need to implement all cases of 1d vector, 2d matrix, etc.
! 
! Difficulty: has to know the output shape beforehand, because the output argument's
! shape has to be specified in the beginning of the subroutine.  For example, concate.
! two col. vectors along the column or along the row would produce different rank (2 or 1)
! of the result.
! 
! Solution: split 'cat' into three categories:
! call cat1(D, A, B)    concat. along dim1   (subroutine: D <== A [+] B)
! call cat2(D, A, B)    concat. along dim2
! call cat3(D, A, B)    concat. along dim3
! Within each category, ranks and types are overloaded
! 
! Also implement functions:
! D = catfun1(A,B)
! D = catfun2(A,B)
! D = catfun3(A,B)

! * * * * * * * * * * * CAT1 * * * * * * * * * * *

SUBROUTINE cat1_iv(D, A, B)  ! integer (column) vector
INTEGER(I4B), INTENT(IN)  :: A(:), B(:)
INTEGER(I4B), INTENT(out) :: D(:)
integer(I4B) :: na, nb, nd
na = size(A);   nb = size(B);  nd = size(D);
call assert(nd == na+nb, ' cat dim1: Ldest /= La + Lb')
D(1:na)     = A
D(na+1 :nd) = B
end SUBROUTINE cat1_iv

SUBROUTINE cat1_rv(D, A, B)  ! real (column) vector
REAL(SP), INTENT(IN)  :: A(:), B(:)
REAL(SP), INTENT(out) :: D(:)
integer(I4B) :: na, nb, nd
na = size(A);   nb = size(B);  nd = size(D);
call assert(nd == na+nb, ' cat dim1: Ldest /= La + Lb')
D(1:na)     = A
D(na+1 :nd) = B
end SUBROUTINE cat1_rv

SUBROUTINE cat1_cv(D, A, B) ! complex (column) vector
COMPLEX(SPC), INTENT(IN)  :: A(:), B(:)
COMPLEX(SPC), INTENT(out) :: D(:)
integer(I4B) :: na, nb, nd
na = size(A);   nb = size(B);  nd = size(D);
call assert(nd == na+nb, ' cat dim1: Ldest /= La + Lb')
D(1:na)     = A
D(na+1 :nd) = B
end SUBROUTINE cat1_cv

SUBROUTINE cat1_dv(D, A, B)  ! double (column) vector
REAL(DP), INTENT(IN)  :: A(:), B(:)
REAL(DP), INTENT(out) :: D(:)
integer(I4B) :: na, nb, nd
na = size(A);   nb = size(B);  nd = size(D);
call assert(nd == na+nb, ' cat dim1: Ldest /= La + Lb')
D(1:na)     = A
D(na+1 :nd) = B
end SUBROUTINE cat1_dv

SUBROUTINE cat1_zv(D, A, B) ! complex double (column) vector
COMPLEX(DPC), INTENT(IN)  :: A(:), B(:)
COMPLEX(DPC), INTENT(out) :: D(:)
integer(I4B) :: na, nb, nd
na = size(A);   nb = size(B);  nd = size(D);
call assert(nd == na+nb, ' cat dim1: Ldest /= La + Lb')
D(1:na)     = A
D(na+1 :nd) = B
end SUBROUTINE cat1_zv


SUBROUTINE cat1_im(D, A, B)  ! integer matrix
INTEGER(I4B), INTENT(IN)  :: A(:,:), B(:,:)
INTEGER(I4B), INTENT(out) :: D(:,:)
integer(I4B) :: sa(2), sb(2), sd(2)
sa = shape(A);   sb = shape(B);  sd = shape(D);
call assert(sd(1)== sa(1)+sb(1), sd(2)==sa(2), sd(2)==sb(2), &
               ' cat dim 1, size inconsist.') 
D(1:sa(1),:)    = A
D(sa(1)+1 :sd(1),:)= B
end SUBROUTINE cat1_im

SUBROUTINE cat1_rm(D, A, B)  ! real matrix
REAL(SP), INTENT(IN)  :: A(:,:), B(:,:)
REAL(SP), INTENT(out) :: D(:,:)
integer(I4B) :: sa(2), sb(2), sd(2)
sa = shape(A);   sb = shape(B);  sd = shape(D);
call assert(sd(1)== sa(1)+sb(1), sd(2)==sa(2), sd(2)==sb(2), &
               ' cat dim 1, size inconsist.') 
D(1:sa(1),:)    = A
D(sa(1)+1 :sd(1),:)= B
end SUBROUTINE cat1_rm

SUBROUTINE cat1_cm(D, A, B) ! complex matrix
COMPLEX(SPC), INTENT(IN)  :: A(:,:), B(:,:)
COMPLEX(SPC), INTENT(out) :: D(:,:)
integer(I4B) :: sa(2), sb(2), sd(2)
sa = shape(A);   sb = shape(B);  sd = shape(D);
call assert(sd(1)== sa(1)+sb(1), sd(2)==sa(2), sd(2)==sb(2), &
               ' cat dim 1, size inconsist.')
D(1:sa(1),:)    = A
D(sa(1)+1 :sd(1),:)= B
end SUBROUTINE cat1_cm

SUBROUTINE cat1_dm(D, A, B)  ! double matrix
REAL(DP), INTENT(IN)  :: A(:,:), B(:,:)
REAL(DP), INTENT(out) :: D(:,:)
integer(I4B) :: sa(2), sb(2), sd(2)
sa = shape(A);   sb = shape(B);  sd = shape(D);
call assert(sd(1)== sa(1)+sb(1), sd(2)==sa(2), sd(2)==sb(2), &
               ' cat dim 1, size inconsist.') 
D(1:sa(1),:)    = A
D(sa(1)+1 :sd(1),:)= B
end SUBROUTINE cat1_dm

SUBROUTINE cat1_zm(D, A, B) ! complex double matrix
COMPLEX(DPC), INTENT(IN)  :: A(:,:), B(:,:)
COMPLEX(DPC), INTENT(out) :: D(:,:)
integer(I4B) :: sa(2), sb(2), sd(2)
sa = shape(A);   sb = shape(B);  sd = shape(D);
call assert(sd(1)== sa(1)+sb(1), sd(2)==sa(2), sd(2)==sb(2), &
               ' cat dim 1, size inconsist.')
D(1:sa(1),:)    = A
D(sa(1)+1 :sd(1),:)= B
end SUBROUTINE cat1_zm

!!!!!!!!!!!!!!  CAT2  !!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE cat2_iv(D, A, B)  ! integer (column) vector
INTEGER(I4B), INTENT(IN)  :: A(:), B(:)
INTEGER(I4B), INTENT(out) :: D(:,:)
integer(I4B) :: sd(2), na, nb
sd = shape(D)
na = size(B);   nb = size(B);
call assert(sd(1) == na, sd(1) == nb, sd(2)==2, ' cat dim 2, shape must match.')
D(:,1) = A
D(:,2) = B
end SUBROUTINE cat2_iv

SUBROUTINE cat2_rv(D, A, B)  ! real (column) vector
REAL(SP), INTENT(IN)  :: A(:), B(:)
REAL(SP), INTENT(out) :: D(:,:)
integer(I4B) :: sd(2), na, nb
sd = shape(D)
na = size(B);   nb = size(B);
call assert(sd(1) == na, sd(1) == nb, sd(2)==2, ' cat dim 2, shape must match.')

D(:,1) = A
D(:,2) = B
end SUBROUTINE cat2_rv

SUBROUTINE cat2_cv(D, A, B)  ! complex (column) vector
COMPLEX(SPC), INTENT(IN)  :: A(:), B(:)
COMPLEX(SPC), INTENT(out) :: D(:,:)
integer(I4B) :: sd(2), na, nb
sd = shape(D)
na = size(B);   nb = size(B);
call assert(sd(1) == na, sd(1) == nb, sd(2)==2, ' cat dim 2, shape must match.')

D(:,1) = A
D(:,2) = B
end SUBROUTINE cat2_cv

SUBROUTINE cat2_dv(D, A, B)  ! double (column) vector
REAL(DP), INTENT(IN)  :: A(:), B(:)
REAL(DP), INTENT(out) :: D(:,:)
integer(I4B) :: sd(2), na, nb
sd = shape(D)
na = size(B);   nb = size(B);
call assert(sd(1) == na, sd(1) == nb, sd(2)==2, ' cat dim 2, shape must match.')

D(:,1) = A
D(:,2) = B
end SUBROUTINE cat2_dv

SUBROUTINE cat2_zv(D, A, B)  ! double complex (column) vector
COMPLEX(DPC), INTENT(IN)  :: A(:), B(:)
COMPLEX(DPC), INTENT(out) :: D(:,:)
integer(I4B) :: sd(2), na, nb
sd = shape(D)
na = size(B);   nb = size(B);
call assert(sd(1) == na, sd(1) == nb, sd(2)==2, ' cat dim 2, shape must match.')

D(:,1) = A
D(:,2) = B
end SUBROUTINE cat2_zv

SUBROUTINE cat2_im(D, A, B)  ! integer matrix
INTEGER(I4B), INTENT(IN)  :: A(:,:), B(:,:)
INTEGER(I4B), INTENT(out) :: D(:,:)
integer(I4B) :: sa(2), sb(2), sd(2)
sa = shape(A);   sb = shape(B);  sd = shape(D);
call assert(sd(2)== sa(2)+sb(2), sd(1)==sa(1), sd(1)==sb(1), &
               ' cat dim 2, size inconsist.') 
D(:, 1:sa(2))    = A
D(:, sa(2)+1 :)  = B
end SUBROUTINE cat2_im

SUBROUTINE cat2_rm(D, A, B)  ! real matrix
REAL(SP), INTENT(IN)  :: A(:,:), B(:,:)
REAL(SP), INTENT(out) :: D(:,:)
integer(I4B) :: sa(2), sb(2), sd(2)
sa = shape(A);   sb = shape(B);  sd = shape(D);
call assert(sd(2)== sa(2)+sb(2), sd(1)==sa(1), sd(1)==sb(1), &
               ' cat dim 2, size inconsist.') 
D(:, 1:sa(2))    = A
D(:, sa(2)+1 :)  = B
end SUBROUTINE cat2_rm

SUBROUTINE cat2_cm(D, A, B)  ! complex matrix
COMPLEX(SPC), INTENT(IN)  :: A(:,:), B(:,:)
COMPLEX(SPC), INTENT(out) :: D(:,:)
integer(I4B) :: sa(2), sb(2), sd(2)
sa = shape(A);   sb = shape(B);  sd = shape(D);
call assert(sd(2)== sa(2)+sb(2), sd(1)==sa(1), sd(1)==sb(1), &
               ' cat dim 2, size inconsist.') 
D(:, 1:sa(2))    = A
D(:, sa(2)+1 :)  = B
end SUBROUTINE cat2_cm

SUBROUTINE cat2_dm(D, A, B)  ! double matrix
REAL(DP), INTENT(IN)  :: A(:,:), B(:,:)
REAL(DP), INTENT(out) :: D(:,:)
integer(I4B) :: sa(2), sb(2), sd(2)
sa = shape(A);   sb = shape(B);  sd = shape(D);
call assert(sd(2)== sa(2)+sb(2), sd(1)==sa(1), sd(1)==sb(1), &
               ' cat dim 2, size inconsist.') 
D(:, 1:sa(2))    = A
D(:, sa(2)+1 :)  = B
end SUBROUTINE cat2_dm

SUBROUTINE cat2_zm(D, A, B)  ! complex matrix
COMPLEX(DPC), INTENT(IN)  :: A(:,:), B(:,:)
COMPLEX(DPC), INTENT(out) :: D(:,:)
integer(I4B) :: sa(2), sb(2), sd(2)
sa = shape(A);   sb = shape(B);  sd = shape(D);
call assert(sd(2)== sa(2)+sb(2), sd(1)==sa(1), sd(1)==sb(1), &
               ' cat dim 2, size inconsist.') 
D(:, 1:sa(2))    = A
D(:, sa(2)+1 :)  = B
end SUBROUTINE cat2_zm
! * * * * * * * * * * * * * * * * * * * * * *

! * * * * * * * CATFUN1 * * * * * * *

FUNCTION catfun1_iv(A, B)  ! integer (column) vector
INTEGER(I4B)  :: A(:), B(:)
INTEGER(I4B)  :: catfun1_iv(size(A)+size(B))
catfun1_iv(1:size(A))   = A
catfun1_iv(size(A)+1 :) = B
end FUNCTION

FUNCTION catfun1_rv(A, B)  ! real (column) vector
REAL(SP)  :: A(:), B(:)
REAL(SP)  :: catfun1_rv(size(A)+size(B))
catfun1_rv(1:size(A))   = A
catfun1_rv(size(A)+1 :) = B
end FUNCTION

FUNCTION catfun1_cv(A, B) ! complex (column) vector
COMPLEX(SPC)  :: A(:), B(:)
COMPLEX(SPC)  :: catfun1_cv(size(A)+size(B))
catfun1_cv(1:size(A))   = A
catfun1_cv(size(A)+1 :) = B
end FUNCTION

FUNCTION catfun1_dv(A, B)  ! double (column) vector
REAL(DP)  :: A(:), B(:)
REAL(DP)  :: catfun1_dv(size(A)+size(B))
catfun1_dv(1:size(A))   = A
catfun1_dv(size(A)+1 :) = B
end FUNCTION

FUNCTION catfun1_zv(A, B) ! double complex (column) vector
COMPLEX(DPC)  :: A(:), B(:)
COMPLEX(DPC)  :: catfun1_zv(size(A)+size(B))
catfun1_zv(1:size(A))   = A
catfun1_zv(size(A)+1 :) = B
end FUNCTION

FUNCTION catfun1_im(A, B)  ! integer matrix
INTEGER(I4B)  :: A(:,:), B(:,:)
INTEGER(I4B)  :: catfun1_im(size(A,1)+size(B,1), size(A,2))
catfun1_im(1:size(A,1),:)   = A
catfun1_im(size(A,1)+1 :,:) = B
end FUNCTION

FUNCTION catfun1_rm(A, B)  ! real matrix
REAL(SP)  :: A(:,:), B(:,:)
REAL(SP)  :: catfun1_rm(size(A,1)+size(B,1), size(A,2))
catfun1_rm(1:size(A,1),:)   = A
catfun1_rm(size(A,1)+1 :,:) = B
end FUNCTION

FUNCTION catfun1_cm(A, B) ! complex matrix
COMPLEX(SPC)  :: A(:,:), B(:,:)
COMPLEX(SPC)  :: catfun1_cm(size(A,1)+size(B,1), size(A,2))
catfun1_cm(1:size(A,1),:)   = A
catfun1_cm(size(A,1)+1 :,:) = B
end FUNCTION

FUNCTION catfun1_dm(A, B)  ! double matrix
REAL(DP)  :: A(:,:), B(:,:)
REAL(DP)  :: catfun1_dm(size(A,1)+size(B,1), size(A,2))
catfun1_dm(1:size(A,1),:)   = A
catfun1_dm(size(A,1)+1 :,:) = B
end FUNCTION

FUNCTION catfun1_zm(A, B) ! double complex matrix
COMPLEX(DPC)  :: A(:,:), B(:,:)
COMPLEX(DPC)  :: catfun1_zm(size(A,1)+size(B,1), size(A,2))
catfun1_zm(1:size(A,1),:)   = A
catfun1_zm(size(A,1)+1 :,:) = B
end FUNCTION

!!!!!!!!!!!!!!!!! CATFUN2 !!!!!!!!!!!!!!!!!!!!

FUNCTION catfun2_iv(A, B)  ! integer (column) vector
INTEGER(I4B)  :: A(:), B(:)
INTEGER(I4B)  :: catfun2_iv(size(A),2)
catfun2_iv(:,1) = A
catfun2_iv(:,2) = B
end FUNCTION
 
FUNCTION catfun2_rv(A, B)  ! real (column) vector
REAL(SP)  :: A(:), B(:)
REAL(SP)  :: catfun2_rv(size(A),2)
catfun2_rv(:,1) = A
catfun2_rv(:,2) = B
end FUNCTION

FUNCTION catfun2_cv(A, B)  ! complex (column) vector
COMPLEX(SPC)  :: A(:), B(:)
COMPLEX(SPC)  :: catfun2_cv(size(A),2)
catfun2_cv(:,1) = A
catfun2_cv(:,2) = B
end FUNCTION

FUNCTION catfun2_dv(A, B)  ! double (column) vector
REAL(DP)  :: A(:), B(:)
REAL(DP)  :: catfun2_dv(size(A),2)
catfun2_dv(:,1) = A
catfun2_dv(:,2) = B
end FUNCTION

FUNCTION catfun2_zv(A, B)  ! double complex (column) vector
COMPLEX(DPC)  :: A(:), B(:)
COMPLEX(DPC)  :: catfun2_zv(size(A),2)
catfun2_zv(:,1) = A
catfun2_zv(:,2) = B
end FUNCTION

FUNCTION catfun2_im(A, B)  ! integer matrix
INTEGER(I4B)  :: A(:,:), B(:,:)
INTEGER(I4B)  :: catfun2_im(size(A,1), size(A,2)+size(B,2))
catfun2_im(:, 1:size(A,2))    = A
catfun2_im(:, size(A,2)+1 :)  = B
end FUNCTION

FUNCTION catfun2_rm(A, B)  ! real matrix
REAL(SP)  :: A(:,:), B(:,:)
REAL(SP)  :: catfun2_rm(size(A,1), size(A,2)+size(B,2))
catfun2_rm(:, 1:size(A,2))    = A
catfun2_rm(:, size(A,2)+1 :)  = B
end FUNCTION

FUNCTION catfun2_cm(A, B)  ! complex matrix
COMPLEX(SPC)  :: A(:,:), B(:,:)
COMPLEX(SPC)  :: catfun2_cm(size(A,1), size(A,2)+size(B,2))
catfun2_cm(:, 1:size(A,2))    = A
catfun2_cm(:, size(A,2)+1 :)  = B
end FUNCTION 

FUNCTION catfun2_dm(A, B)  ! double matrix
REAL(DP)  :: A(:,:), B(:,:)
REAL(DP)  :: catfun2_dm(size(A,1), size(A,2)+size(B,2))
catfun2_dm(:, 1:size(A,2))    = A
catfun2_dm(:, size(A,2)+1 :)  = B
end FUNCTION

FUNCTION catfun2_zm(A, B)  ! double complex matrix
COMPLEX(DPC)  :: A(:,:), B(:,:)
COMPLEX(DPC)  :: catfun2_zm(size(A,1), size(A,2)+size(B,2))
catfun2_zm(:, 1:size(A,2))    = A
catfun2_zm(:, size(A,2)+1 :)  = B
end FUNCTION 
! ****----------- CAT END ---------- ****

END MODULE


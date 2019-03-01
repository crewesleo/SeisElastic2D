module seism_uti_module
use nrtype; use nrutil; use myMath_module; use IEEE_ARITHMETIC

IMPLICIT NONE
REAL(SP), PARAMETER :: RECL_MAX = 1.E8, TRUNC_F_RATIO = 2.5 !2.7638
! ricker spectrum truncation ratio: f_hi/f_peak
integer ::STATbuff(13), status  !used to get file size
private ::STATbuff, status

INTERFACE read_binfile  ! read_binfile(filename, data, <recn>); <recn> <= 1  --> read all
! otherwise, data is a subarray in filename, w/ record length being data's length
  module procedure &
    read_binfile_b1, read_binfile_b2, &
    read_binfile_i1, read_binfile_i2, read_binfile_i3, &
    read_binfile_r1, read_binfile_r2, read_binfile_r3, &
    read_binfile_d1, read_binfile_d2, read_binfile_d3
END INTERFACE

INTERFACE write_binfile ! write_binfile(filename, data, <rec_num>), rec_num<0 --> append                  
   module procedure &  
    write_binfile_b1, write_binfile_b2, &
    write_binfile_i1, write_binfile_i2, write_binfile_i3, &
    write_binfile_r0, write_binfile_r1, write_binfile_r2, write_binfile_r3, &
    write_binfile_d1, write_binfile_d2, write_binfile_d3, &
    write_binfile_c1, write_binfile_c2, write_binfile_c3
END INTERFACE ! b: byte

INTERFACE read_tex
    module procedure read_tex_i1, read_tex_r1  ! 1 columned text file
END INTERFACE

INTERFACE write_tex    ! write_tex(filename, data, flg)  ! flg:  optional
!  (In MPI debugging setting, typically individual thread has its own filename.)
!  flg:  fastdim_or_clearup flag
!  flg <= 0  --> clear up the file [while the default is append]
!      == +-1 or +-2 --> which dim of the matrix is printed out the fastest, 1 or 2? 
!               (i.e., on the same line)  [default: 2, visually intuitive.]
! b:  integer(1)
! l1: logical(1)
    module procedure  write_tex_char, &
      write_tex_i0, write_tex_r0, write_tex_l10,write_tex_l0, write_tex_b0,&
      write_tex_i1, write_tex_r1, write_tex_c1, write_tex_z1, write_tex_l11,write_tex_l1,write_tex_b1,&
      write_tex_i2, write_tex_r2, write_tex_c2, write_tex_z2
END INTERFACE

INTERFACE ricker   ! RICKER(RF,F0,T)  or  RICKER(RF,Fmx,N),  RF is pointer
  module procedure ricker_pk, ricker_hi ! use real() or integer() of the
END INTERFACE                           ! 3rd argument to disambiguate


CONTAINS  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!subroutine nega_diff2( y, x )  ! y = -x"  ! no longer needed
!real,      intent(out) :: y(:)
!real,      intent(in)  :: x(:)
!integer    :: i, sz
!
!sz    = size(x)
!do i  = 2, sz-1
!  y(i)= 2*x(i) - x(i-1) - x(i+1)
!end do
!y(1)  = y(2)
!y(sz) = y(sz-1)
!end subroutine
!------------------------------
subroutine freeze_water_velocity(A, z_fix_v, waterBtm, v0)
!          freeze water velocity
! if z_fix_v <=0, do nothing
! else,  if size(waterBtm) > 0,  water region: z_fix_v <= z <= waterBtm
! else,  water region:  1<= z <= z_fix_v
!  if v0 absent, it's 0, i.e., to 0 out the GRADIENT associated w/ water
!  x is fast dim.
real,         intent(inout) :: A(:,:)
integer,      intent(in)    :: z_fix_v
real,         intent(in)    :: waterBtm(:)
real,optional,intent(in)    :: v0
real          :: valu
integer       :: ix,  iz,  nx,  nz

if( z_fix_v < 1)           return   ! no water layer

valu = 0.;    if(present(v0))  valu = v0
if(size(waterBtm)==0)      then
  A(:, 1 : z_fix_v) = valu
  return
end if

nx   = size(A,1);       nz = size(A,2)
call assert(nx == size(waterBtm), 'waterBtm length must == nx')
do ix= 1, nx
   A( ix, z_fix_v : floor(waterBtm(ix)) ) = valu
end do
end subroutine

subroutine misfit_adaptGain(e, e_nrm2, c, o, ng)
real,   intent(out) :: e(:,:), e_nrm2 ! e:  alpha c  -  o
! e_nrm2: ||e||^2.  For each column, find alpha to minimize e_nrm2.
real,   intent(in)  :: c(:,:), o(:,:) ! calculated, observed
integer,intent(in)  :: ng      ! only deals with 1st ng columns

integer :: i
real    :: alpha
e_nrm2 = 0.
do i     = 1, ng
    alpha  = dot_prod(c(:,i), o(:,i)) / norm2(c(:,i))
    e(:,i) = alpha* c(:,i) - o(:,i)
    e_nrm2 = e_nrm2 + norm2(e(:,i))
end do
end subroutine
!-------------------------------
subroutine showArryInfo(arr, nam)
! show array statistics, for diagnosis
real, intent(in) :: arr(:,:)
character(len=*) :: nam
write(*,*) trim(nam) // ':, min, max, mean, norm_1/N', &
  minval(arr), maxval(arr), sum(arr)/size(arr), norm_1(arr)/size(arr)
call flush(6)
end subroutine

function anyNan(A) result(has)
logical  :: has
real,intent(in) :: A(:,:)
integer  :: i,j, nr  ! # rows
has  = .false.
nr   = size(A,1)
do j = 1, size(A,2)
  do i = 1, nr
    has = ieee_is_Nan(A(i,j));  if(has) return 
  end do
end do
end function

subroutine findPaddingWidth(pdl,pdr,ny)
integer, intent(out) :: pdl,pdr
integer, intent(in)  :: ny
! padding width to the left & to the right: pdl, pdr >= 25
! total width pdl + ny + pdr  consists of small prime factors (for ease of FFT)
! for phase Shift migration, padding & tapering to kill space wrap-around.
integer  :: L, w
L   = roundUp2SmallPrimeProduct(ny+50) - ny
w   = L/2
pdl = w
pdr = w + (L-w-w)
end subroutine
!-----------------
subroutine weights4Padding(wtl,wtr,pdl,pdr)
real, pointer, intent(inout):: wtl(:), wtr(:)
integer,       intent(in)   :: pdl,pdr
! weights for left & right padding to give a Tukey window
! for phase Shift migration, padding & tapering to kill space wrap-around.
integer         :: L
real, parameter :: pi = 3.14159265
IF(ASSOCIATED(wtl))   DEALLOCATE(wtl)
IF(ASSOCIATED(wtr))   DEALLOCATE(wtr)
L   = pdl
allocate(wtl(L))
wtl = (1+ cos( arth(-L,1,L)*(pi/(L+1)) ))/2
L   = pdr
allocate(wtr(L))
wtr = (1+ cos( arth(1,1,L)*(pi/(L+1)) ))/2
end subroutine
!-----------------
SUBROUTINE pad_slow_ref_y_3d(slw, ref, szx, szy,  hdx, hdy, slwin, refin)
! slw, ref      nx  x ny  x nz  slowness and reflectance
! slwin, refin  nx1 x ny1 x nz  input slowness and reflectance
! get a chunk from slwin/refin, starting @ (hdx, hdy, 1), of size (szx, szy, nz) which is smaller
! put this chunk in slw/ref starting @ (1,1,1), and then
! pad along y, from szy+1 to ny, for all x \in [1,szx]
! if  hdx, hdy, slwin, refin not present. then assume the chunk has already been copied to slw/ref
! The way the padding works is smoothed slowness: succeedingly a grid point (x,y,z)'s slowness
! becomes the average of the two neighbors' average (x, y-1, z) (x, y-1, z+1)
real, intent(inout) :: slw(:,:,:), ref(:,:,:)
integer, intent(in) :: szx, szy
integer, intent(in), optional :: hdx, hdy
real,    intent(in), optional :: slwin(:,:,:), refin(:,:,:)
integer  :: ny, nz, iy
if (present(refin)) then
  slw(1:szx, 1:szy, :) = slwin(hdx:hdx-1+szx, hdy:hdy-1+szy, :)
  ref(1:szx, 1:szy, :) = refin(hdx:hdx-1+szx, hdy:hdy-1+szy, :)
endif
nz = size(slw, 3);  ny = size(slw, 2)
do iy = szy+1, ny
  slw(1:szx, iy, 1:nz-1) = (slw(1:szx, iy-1, 1:nz-1) + slw(1:szx, iy-1, 2:nz))/2.
  slw(1:szx, iy, nz)     =  slw(1:szx, iy, nz-1)
enddo
ref(1:szx, szy+1:ny, nz) = 0.
ref(1:szx, szy+1:ny, 1:nz-1) = (slw(1:szx, szy+1:ny, 1:nz-1) - slw(1:szx, szy+1:ny, 2:nz))&
   /(slw(1:szx, szy+1:ny, 1:nz-1) + slw(1:szx, szy+1:ny, 2:nz))
end SUBROUTINE

SUBROUTINE pad_slow_ref_y_2d(slw, ref, szy,  hdy, slwin, refin)
! similar to pad_slow_ref_y_3d, except that nx==1
real, intent(inout) :: slw(:,:,:), ref(:,:,:)
integer, intent(in) :: szy
integer, intent(in), optional :: hdy
real,    intent(in), optional :: slwin(:,:,:), refin(:,:,:)

integer  :: ny, nz, iy
if (present(refin)) then
  slw(1, 1:szy, :) = slwin(1, hdy:hdy-1+szy, :)
  ref(1, 1:szy, :) = refin(1, hdy:hdy-1+szy, :)
endif
nz = size(slw, 3);  ny = size(slw, 2)
do iy = szy+1, ny
  slw(1, iy, 1:nz-1) = (slw(1, iy-1, 1:nz-1) + slw(1, iy-1, 2:nz))/2.
  slw(1, iy, nz)     =  slw(1, iy, nz-1)
enddo
ref(1, szy+1:ny, nz) = 0.
ref(1, szy+1:ny, 1:nz-1) = (slw(1, szy+1:ny, 1:nz-1) - slw(1, szy+1:ny, 2:nz))&
   /(slw(1, szy+1:ny, 1:nz-1) + slw(1, szy+1:ny, 2:nz))
end subroutine

SUBROUTINE pad_slow_ref_x(slw, ref, szx)
! slw, ref      nx  x ny  x nz  slowness and reflectance
! A chunk of size (szx, :, :) is in slw/ref
! Pad along x, from szx+1 to nx, for all y, z
! (assume having padded along y first.  Thus the chunk is already in memory.)
real, intent(inout) :: slw(:,:,:), ref(:,:,:)
integer, intent(in) :: szx
integer  :: nx, nz, ix

nx = size(slw, 1);  nz = size(slw, 3)
do ix = szx+1, nx
  slw(ix, :, 1:nz-1) = (slw(ix-1, :, 1:nz-1) + slw(ix-1, :, 2:nz))/2.
  slw(ix, :, nz)     =  slw(ix, :, nz-1)
enddo
ref(szx+1:nx, :, nz) = 0.
ref(szx+1:nx, :, 1:nz-1) = (slw(szx+1:nx, :, 1:nz-1) - slw(szx+1:nx, :, 2:nz))/(slw(szx+1:nx, :, 1:nz-1) + slw(szx+1:nx, :, 2:nz))
end SUBROUTINE

SUBROUTINE pad_slow_ref(slw, ref, szx, szy,  hdx, hdy, slwin, refin)
! arguments similar to those in  pad_slow_ref_y_3d
! pad along y and then along x (i.e., for real 3d data)
real, intent(inout) :: slw(:,:,:), ref(:,:,:)
integer, intent(in) :: szx, szy
integer, intent(in) :: hdx, hdy
real,    intent(in) :: slwin(:,:,:), refin(:,:,:)

call pad_slow_ref_y_3d(slw, ref, szx, szy,  hdx, hdy, slwin, refin)
call pad_slow_ref_x(slw, ref, szx)
end SUBROUTINE



!-----------------------------------------------------
function file_size(filename) !! file size (in bytes)
character(len=*), intent(in) :: filename
integer  :: file_size, fstat_buff(13), status
!! through compiler flag -D_GNUf90_, if using gnu compiler:
#ifdef _GNUf90_
  CALL STAT(filename, fstat_buff, status)
  file_size = fstat_buff(8)
  if (STATUS /= 0)  file_size = 0  ! file doesn't exist
#else
  inquire(FILE=filename, SIZE=file_size)  
  if (file_size < 0)  file_size = 0
#endif
end function


!-------------------------------
SUBROUTINE ricker_t(s, f0, dt)  ! Ricker in time domain, peak freq f0 Hz
! fill the time trace s w/ Ricker wavelet from time 0+
real, intent(out) :: s(:)
real, intent(in)  :: f0, dt
integer :: np, i, L
real, allocatable :: t(:), env(:)
real  ::  sm_ric,  sm_env
np = min(ceiling(1.23456789/(f0*dt)),  (size(s)-1)/2)
L  = 2*np+1 
allocate(t(L), env(L))
t  = (/ (i, i=-np, np) /)
t  = (t * (PI *f0 *dt) )**2
env     = exp(-t)   ! envelope
sm_env  = sum(env)
s(1:L)  = (1. - 2.*t)*env
sm_ric  = sum(s(1:L))
s(1:L)  = s(1:L) - sm_ric/sm_env * env ! 0-mean
s(L+1:) = 0.
#ifdef DEBUG
  write(*,*) 'lengths of source & ricker:', size(s), L
#endif
deallocate(t, env)
end SUBROUTINE
!------------------
SUBROUTINE ricker_pk(rf, f0, T, idx1,idx2)
! 0-phase ricker wavelet in freq. domain, truncated
! at high end where the value is 1% of the peak.
! determined by its peak frequency, and time duration
! f0:  peak freq.
! T:   duration in time domain, which determines df = 1/T
! idx1,idx2   the index range of the ricker to be computed and sent to out
! Usage:
!   ricker_pk(ptr, f0, T, idx1,idx2) 
!   where ptr is a pointer to real array, e.g., REAL, POINTER :: ptr(:)
! **Make sure in the end of main program to deallocate(ptr).

USE nrutil,  ONLY: arth
IMPLICIT  NONE
REAL(SP),    intent(in):: f0, T
REAL(SP),    POINTER, CONTIGUOUS   :: rf(:)  ! <====
INTEGER(I4B), intent(in), optional :: idx1,idx2 ! only return rf(idx1:idx2) chunk
REAL(SP)     :: df, val, alph
INTEGER(I4B) :: NumF, totNf, i, i1, i2
    
df   = 1./T
IF(ASSOCIATED(rf))   DEALLOCATE(rf)
totNf= FLOOR(f0*TRUNC_F_RATIO / df)
if (present(idx1))   then
  i2 = idx2;   i1 = idx1
else
  i2 = totNf;  i1 = 1
endif
NumF = i2 - i1 + 1;  ALLOCATE (rf(NumF))  
rf = (arth(df*i1, df, NumF) / f0)**2
rf = rf * exp(-rf + 1)   ! set max to 1 (in Freq. dom.):
! now, to make rf(totNf+1)==0, to remove truncation step induced ripples:
! @ totNf+1, value 
val= (df*(totNf+1) / f0)**2
val= val * exp(-val + 1)
alph = val / (totNf+1)**2  ! alph * i*i is the parabola to be subtracted
do i = i1, i2
  rf(i-i1+1) = rf(i-i1+1) - alph * (i*i)
enddo
END SUBROUTINE
!-------------------------------
SUBROUTINE ricker_hi(rf, fhi, totNf, idx1,idx2)
! 0-phase ricker wavelet in freq. domain, truncated
! at high end where the value is 1% of the peak.
! determined by its high frequency limit, and # of points
! fhi:    upper freq. limit
! totNf:  # freq. points
! idx1,idx2   the index range of the ricker to be computed and sent to out
! Usage:
!   ricker_hi(ptr, fhi, totNf, idx1,idx2) 
!   where ptr is a pointer to real array, e.g., REAL, POINTER :: ptr(:)
! **Make sure in the end of main program to deallocate(ptr).
USE nrutil,  ONLY: arth
IMPLICIT  NONE
REAL(SP),    intent(in):: fhi
REAL(SP),    POINTER, CONTIGUOUS :: rf(:)  ! <====
INTEGER(I4B),intent(in):: totNf
INTEGER(I4B),intent(in),  optional ::  idx1,idx2 ! only return rf(idx1:idx2) chunk
REAL(SP)     :: df, f0, val, alph
integer(I4B) :: NumF, i, i1, i2
    
f0   = fhi / TRUNC_F_RATIO
df   = fhi / totNf
IF(ASSOCIATED(rf))   DEALLOCATE(rf)
if (present(idx1))   then
  i2 = idx2;   i1 = idx1
else
  i2 = totNf;  i1 = 1
endif
NumF = i2 - i1 + 1;  ALLOCATE (rf(NumF))  
rf = (arth(df*i1, df, NumF) / f0)**2
rf = rf * exp(-rf + 1)   ! set to norm 1 (in Freq. domain):
! now, to make rf(totNf+1)==0, to remove truncation step induced ripples:
! @ totNf+1, value 
val= (df*(totNf+1) / f0)**2
val= val * exp(-val + 1)
alph = val / (totNf+1)**2  ! alph * i*i is the parabola to be subtracted
do i = i1, i2
  rf(i-i1+1) = rf(i-i1+1) - alph * (i*i)
enddo
END SUBROUTINE

!-------------------------------
SUBROUTINE read_binfile_b1(filename, data, recn)
character(len=*), intent(in)              :: filename
integer(1),       intent(out)             :: data(:)
integer,          intent(inout), optional :: recn    
! input: if < 0 --> automatically figure out # of data elements read, and return
!        if > 0 --> read recn elements.
! if not present--> if file size smaller than data size, issue error; else, like the <0 case
integer  :: nf, nd, n
nf = file_size(filename);   nd = size(data) ! bytes
if ((.not. present(recn)) .and. nf < nd) then
  write(*,*) 'read_binfile_b1: ', trim(filename), ' size ', nf, &
             ' is smaller than data size ', nd;     stop
end if
n = -1;  if(present(recn))  n  = recn
if(n<0)  n = min(nf, nd)
if(present(recn))   recn = n

open(10, file=filename, access='direct', err=20, recl = n)
read(10, rec=1) data(1:n)
close(10)
return
20 stop '  read error read_binfile_b1'
END SUBROUTINE

SUBROUTINE read_binfile_i1(filename, data, recn)
character(len=*), intent(in)              :: filename
integer,          intent(out)             :: data(:)
integer,          intent(inout), optional :: recn    
! input: if < 0 --> automatically figure out # of data elements read, and return
!        if > 0 --> read recn elements.
! if not present--> if file size smaller than data size, issue error; else, like the <0 case
integer  :: nf, nd, n
nf = file_size(filename)/4;   nd = size(data) ! in 4 bytes
if ((.not. present(recn)) .and. nf < nd) then
  write(*,*) 'read_binfile_b1: ', trim(filename), ' size/4 ', nf, &
             ' is smaller than data size ', nd;     stop
end if
n = -1;  if(present(recn))  n  = recn
if(n<0)  n = min(nf, nd)
if(present(recn))   recn = n

open(10, file=filename, access='direct', err=20, recl=I4*n)
read(10, rec=1) data(1:n)
close(10)
return
20 stop '  read error read_binfile_i1'
END SUBROUTINE
!-------------------------------
SUBROUTINE read_binfile_b2(filename, data, recn)
character(len=*), intent(in) :: filename
integer(1), intent(out) :: data(:,:)
integer, intent(in), optional :: recn
integer  :: n
n = 1
if(present(recn))  n = recn
if (n<=0)    n = 1
open(10, file=filename, access='direct', err=20, recl=size(data))
read(10, rec=1) data
close(10)
return
20 stop '  read error read_binfile_b2'
END SUBROUTINE

SUBROUTINE read_binfile_i2(filename, data, recn)
character(len=*), intent(in) :: filename
integer, intent(out) :: data(:,:)
integer, intent(in), optional :: recn
integer  :: n
n = 1
if(present(recn))  n = recn
if (n<=0)    n = 1
open(10, file=filename, access='direct', err=20, recl=I4*size(data))
read(10, rec=1) data
close(10)
return
20 stop '  read error read_binfile_i2'
END SUBROUTINE
!-------------------------------
SUBROUTINE read_binfile_i3(filename, data)
character(len=*), intent(in) :: filename
integer, intent(out) :: data(:,:,:)
open(10, file=filename, access='direct', err=20, recl=I4*size(data))
read(10, rec=1) data
close(10)
return
20 stop '  read error read_binfile_i3'
END SUBROUTINE

!-------------------------------
SUBROUTINE read_binfile_r1(filename, data, recn)
character(len=*), intent(in)              :: filename
real(SP),         intent(out)             :: data(:)
integer,          intent(inout), optional :: recn    
! input: if < 0 --> automatically figure out # of data elements read, and return
!        if > 0 --> read recn elements.
! if not present--> if file size smaller than data size, issue error; else, like the <0 case
integer  :: nf, nd, n
nf = file_size(filename)/4;   nd = size(data) ! in 4 bytes
if ((.not. present(recn)) .and. nf < nd) then
  write(*,*) 'read_binfile_b1: ', trim(filename), ' size/4 ', nf, &
             ' is smaller than data size ', nd;     stop
end if
n = -1;  if(present(recn))  n  = recn
if(n<0)  n = min(nf, nd)
if(present(recn))   recn = n

open(10, file=filename, access='direct', err=20, recl=I4*n)
read(10, rec=1) data(1:n)
close(10)
return
20 stop '  read error read_binfile_r1'
END SUBROUTINE
!-------------------------------
SUBROUTINE read_binfile_r2(filename, data, recn)
character(len=*), intent(in) :: filename
real(SP), intent(out) :: data(:,:)
integer,  intent(in), optional :: recn
integer  :: n
n = 1
if(present(recn))  n = recn
if (n<=0)    n = 1
open(10, file=filename, access='direct', err=20, recl=I4*size(data))
read(10, rec=n) data
close(10)
return
20 stop '  read error read_binfile_r2'
END SUBROUTINE
!-------------------------------
SUBROUTINE read_binfile_r3(filename, data)
character(len=*), intent(in) :: filename
real(SP), intent(out) :: data(:,:,:)
open(10, file=filename, access='direct', err=20, recl=I4*size(data))
read(10, rec=1) data
close(10)
return
20 stop '  read error read_binfile_r3'
END SUBROUTINE
!-------------------------------
SUBROUTINE read_binfile_d1(filename, data, recn)
character(len=*), intent(in)              :: filename
real(DP),         intent(out)             :: data(:)
integer,          intent(inout), optional :: recn    
! input: if < 0 --> automatically figure out # of data elements read, and return
!        if > 0 --> read recn elements.
! if not present--> if file size smaller than data size, issue error; else, like the <0 case
integer  :: nf, nd, n
nf = file_size(filename)/8;   nd = size(data) ! in 8 bytes
if ((.not. present(recn)) .and. nf < nd) then
  write(*,*) 'read_binfile_b1: ', trim(filename), ' size/8 ', nf, &
             ' is smaller than data size ', nd;     stop
end if
n = -1;  if(present(recn))  n  = recn
if(n<0)  n = min(nf, nd)
if(present(recn))   recn = n

open(10, file=filename, access='direct', err=20, recl=2*I4*size(data))
read(10, rec=1) data(1:n)
close(10)
return
20 stop '  read error read_binfile_d1'
END SUBROUTINE
!-------------------------------
SUBROUTINE read_binfile_d2(filename, data, recn)
character(len=*), intent(in) :: filename
real(DP), intent(out) :: data(:,:)
integer,  intent(in), optional :: recn
integer  :: n
n = 1
if(present(recn))  n = recn
if (n<=0)    n = 1
open(10, file=filename, access='direct', err=20, recl=2*I4*size(data))
read(10, rec=n) data
close(10)
return
20 stop '  read error read_binfile_d2'
END SUBROUTINE
!-------------------------------
SUBROUTINE read_binfile_d3(filename, data)
character(len=*), intent(in) :: filename
real(DP), intent(out) :: data(:,:,:)
open(10, file=filename, access='direct', err=20, recl=2*I4*size(data))
read(10, rec=1) data
close(10)
return
20 stop '  read error read_binfile_d3'
END SUBROUTINE

!================================
SUBROUTINE write_binfile_r0(filename, data, rcn) ! real, scalar
character(len=*), intent(in) :: filename
real(SP),         intent(in) :: data
integer,optional, intent(in) :: rcn     ! the record # to be written at
integer  :: recn, leng

leng = I4    ! record length
if (present(rcn)) then
  if (rcn > 0) then
    recn = rcn
  else
    recn = file_size(filename)  !!!! file size (in bytes)
    recn = recn / leng + 1
  endif 
else
  recn = 1  ! default
endif    
#ifdef DEBUG1
 write(*,*) trim(filename), ' write rec # ', recn
#endif 
if (recn==1) then
  open(10, file=filename, access='direct', recl=leng,status='replace')
else
  open(10, file=filename, access='direct', recl=leng)  
endif
write(10,rec=recn)   data  ! write to # recn
call flush(10)
close(10)
END SUBROUTINE

SUBROUTINE write_binfile_b1(filename, data, rcn) ! byte, 1d array
character(len=*), intent(in) :: filename
integer(1),       intent(in) :: data(:)
integer,optional, intent(in) :: rcn     ! the record # to be written at, assuming we know
  ! the # of previously written records.  If we don't, then set   RCN = -1, and it's gonna be
  ! found out automatically what's the current record #.  We then APPEND data to the file.
integer  :: recn, leng

leng = size(data)    ! record length, make sure the compiler assumes a byte has length 1
if (present(rcn)) then
  if (rcn > 0) then
    recn = rcn
  else
    recn = file_size(filename)  !!!! file size (in bytes)
    recn = recn / leng + 1
  endif 
else
  recn = 1  ! default
endif
if (recn==1) then
  open(10, file=filename, access='direct', recl=leng,status='replace')
else
  open(10, file=filename, access='direct', recl=leng)  
endif
write(10,rec=recn)   data  ! write to # recn
call flush(10)
close(10)
END SUBROUTINE

SUBROUTINE write_binfile_i1(filename, data, rcn) ! integer, 1d array
character(len=*), intent(in) :: filename
integer(I4B),     intent(in) :: data(:)
integer,optional, intent(in) :: rcn     ! the record # to be written at, assuming we know
  ! the # of previously written records.  If we don't, then set   RCN = -1, and it's gonna be
  ! found out automatically what's the current record #.  We then APPEND data to the file.
integer  :: recn, leng

leng = size(data)*I4    ! record length
if (present(rcn)) then
  if (rcn > 0) then
    recn = rcn
  else
    recn = file_size(filename)  !!!! file size (in bytes)
    recn = recn / leng + 1
  endif 
else
  recn = 1  ! default
endif
if (recn==1) then
  open(10, file=filename, access='direct', recl=leng,status='replace')
else
  open(10, file=filename, access='direct', recl=leng)  
endif
write(10,rec=recn)   data  ! write to # recn
call flush(10)
close(10)
END SUBROUTINE

SUBROUTINE write_binfile_r1(filename, data, rcn) ! real, 1d array
character(len=*), intent(in) :: filename
real(SP),         intent(in) :: data(:)
integer,optional, intent(in) :: rcn     ! the record # to be written at
integer  :: recn, leng

leng = size(data)*I4    ! record length
if (present(rcn)) then
  if (rcn > 0) then
    recn = rcn
  else
    recn = file_size(filename)  !!!! file size (in bytes)
    recn = recn / leng + 1
  endif 
else
  recn = 1  ! default
endif    
if (recn==1) then
  open(10, file=filename, access='direct', recl=leng,status='replace')
else
  open(10, file=filename, access='direct', recl=leng)  
endif
write(10,rec=recn)   data  ! write to # recn
call flush(10)
close(10)
END SUBROUTINE

SUBROUTINE write_binfile_d1(filename, data, rcn) ! double, 1d array
character(len=*), intent(in) :: filename
real(DP),         intent(in) :: data(:)
integer,optional, intent(in) :: rcn     ! the record # to be written at
integer  :: recn, leng

leng = 2*size(data)*I4    ! record length
if (present(rcn)) then
  if (rcn > 0) then
    recn = rcn
  else
    recn = file_size(filename)  !!!! file size (in bytes)
    recn = recn / leng + 1
  endif 
else
  recn = 1  ! default
endif    
if (recn==1) then
  open(10, file=filename, access='direct', recl=leng,status='replace')
else
  open(10, file=filename, access='direct', recl=leng)  
endif
write(10,rec=recn)   data  ! write to # recn
call flush(10)
close(10)
END SUBROUTINE

SUBROUTINE write_binfile_c1(filename, data, rcn) ! complex, 1d array
character(len=*), intent(in) :: filename
complex(SPC),     intent(in) :: data(:)
integer,optional, intent(in) :: rcn     ! the record # to be written at
integer  :: recn, leng

leng = size(data)*2*I4    ! record length
if (present(rcn)) then
  if (rcn > 0) then
    recn = rcn
  else
    recn = file_size(filename)  !!!! file size (in bytes)
    recn = recn / leng + 1
  endif 
else
  recn = 1  ! default
endif    
if (recn==1) then
  open(10, file=filename, access='direct', recl=leng,status='replace')
else
  open(10, file=filename, access='direct', recl=leng)  
endif
write(10,rec=recn)   data  ! write to # recn
call flush(10)
close(10)
END SUBROUTINE
!=======
SUBROUTINE write_binfile_b2(filename, data, rcn) ! byte, 2d array
character(len=*), intent(in) :: filename
integer(1),       intent(in) :: data(:,:)
integer,optional, intent(in) :: rcn     ! the record # to be written at
integer  :: recn, leng

leng = size(data)    ! record length
if (present(rcn)) then
  if (rcn > 0) then
    recn = rcn
  else
    recn = file_size(filename)  !!!! file size (in bytes)
    recn = recn / leng + 1 
  endif 
else
  recn = 1  ! default
endif    
if (recn==1) then
  open(10, file=filename, access='direct', recl=leng,status='replace')
else
  open(10, file=filename, access='direct', recl=leng)  
endif
write(10,rec=recn)   data  ! write to # recn
call flush(10)
close(10)
END SUBROUTINE

SUBROUTINE write_binfile_i2(filename, data, rcn) ! integer, 2d array
character(len=*), intent(in) :: filename
integer(I4B),     intent(in) :: data(:,:)
integer,optional, intent(in) :: rcn     ! the record # to be written at
integer  :: recn, leng

leng = size(data)*I4    ! record length
if (present(rcn)) then
  if (rcn > 0) then
    recn = rcn
  else
    recn = file_size(filename)  !!!! file size (in bytes)
    recn = recn / leng + 1 
  endif 
else
  recn = 1  ! default
endif    
if (recn==1) then
  open(10, file=filename, access='direct', recl=leng,status='replace')
else
  open(10, file=filename, access='direct', recl=leng)  
endif
write(10,rec=recn)   data  ! write to # recn
call flush(10)
close(10)
END SUBROUTINE

SUBROUTINE write_binfile_r2(filename, data, rcn) ! real, 2d array
character(len=*), intent(in) :: filename
real(SP),         intent(in) :: data(:,:)
integer,optional, intent(in) :: rcn     ! the record # to be written at
integer  :: recn, leng

leng = size(data)*I4    ! record length
if (present(rcn)) then
  if (rcn > 0) then
    recn = rcn
  else
    recn = file_size(filename)  !!!! file size (in bytes)
    recn = recn / leng + 1 
  endif 
else
  recn = 1  ! default
endif    
if (recn==1) then
  open(10, file=filename, access='direct', recl=leng,status='replace')
else
  open(10, file=filename, access='direct', recl=leng)  
endif
write(10,rec=recn)   data  ! write to # recn
call flush(10)
close(10)
END SUBROUTINE

SUBROUTINE write_binfile_d2(filename, data, rcn) ! double, 2d array
character(len=*), intent(in) :: filename
real(DP),         intent(in) :: data(:,:)
integer,optional, intent(in) :: rcn     ! the record # to be written at
integer  :: recn, leng

leng = size(data)*2*I4    ! record length
if (present(rcn)) then
  if (rcn > 0) then
    recn = rcn
  else
    recn = file_size(filename)  !!!! file size (in bytes)
    recn = recn / leng + 1 
  endif 
else
  recn = 1  ! default
endif    
if (recn==1) then
  open(10, file=filename, access='direct', recl=leng,status='replace')
else
  open(10, file=filename, access='direct', recl=leng)  
endif
write(10,rec=recn)   data  ! write to # recn
call flush(10)
close(10)
END SUBROUTINE

SUBROUTINE write_binfile_c2(filename, data, rcn) ! complex, 2d array
character(len=*), intent(in) :: filename
complex(SPC),     intent(in) :: data(:,:)
integer,optional, intent(in) :: rcn     ! the record # to be written at
integer  :: recn, leng

leng = size(data)*2*I4    ! recod length
if (present(rcn)) then
  if (rcn > 0) then
    recn = rcn
  else
    recn = file_size(filename)  !!!! file size (in bytes)
    recn = recn / leng + 1
  endif 
else
  recn = 1  ! default
endif    
if (recn==1) then
  open(10, file=filename, access='direct', recl=leng,status='replace')
else
  open(10, file=filename, access='direct', recl=leng)  
endif
write(10,rec=recn)   data  ! write to # recn
call flush(10)
close(10)
END SUBROUTINE
!==
SUBROUTINE write_binfile_i3(filename, data, rcn) ! integer, 3d array
character(len=*), intent(in) :: filename
integer(I4B),     intent(in) :: data(:,:,:)
integer(I4B)   :: shp(3), k
integer,optional, intent(in) :: rcn ! the record # to be written at, assuming record length as datasize
integer  :: recn, leng
logical  :: findSize, doSlice

shp = shape(data)
leng = product(shp) *I4
if (leng < RECL_MAX) THEN
  doSlice = .false.  ! can do 3d in 1 step
else
  leng    = shp(1)*shp(2) *I4
  doSlice = .true.   ! must do 1 slice @ a time to avoid write buffer overflow
endif  

findSize = .false.
if (present(rcn)) then
  if (rcn > 0) then
    recn = rcn
  else          ! rcn < 0
    recn = file_size(filename)  !!!! file size (in bytes)
    recn = recn / leng + 1
    if (recn > 1)  findSize = .true.  ! we found the filesize to get the (correct) recn value
  endif 
else
  recn = 1  ! default
endif 

if (doSlice .and. .not. findSize)   recn = (recn - 1) * shp(3) + 1  
!^ because the input rcn was counted as the # of 3d matrices, but here 1 record is a slice only,
!  so recn needs to be recomputed.

if (recn==1) then
  open(10, file=filename, access='direct', recl=leng,status='replace')
else
  open(10, file=filename, access='direct', recl=leng)  
endif

if (doSlice) then
  do k = 1, shp(3)                
    write(10,rec=recn+k-1) data(:,:,k)
  enddo
else  
  write(10,rec=recn)  data 
endif
call flush(10)
close(10)
END SUBROUTINE

SUBROUTINE write_binfile_r3(filename, data, rcn) ! real, 3d array
character(len=*), intent(in) :: filename
real(SP),         intent(in) :: data(:,:,:)
integer(I4B)   :: shp(3), k
integer,optional, intent(in) :: rcn ! the record # to be written at, assuming record length as datasize
integer  :: recn, leng
logical  :: findSize, doSlice

shp = shape(data)
leng = product(shp) *I4
if (leng < RECL_MAX) THEN
  doSlice = .false.  ! can do 3d in 1 step
else
  leng    = shp(1)*shp(2) *I4
  doSlice = .true.   ! must do 1 slice @ a time to avoid write buffer overflow
endif  

findSize = .false.
if (present(rcn)) then
  if (rcn > 0) then
    recn = rcn
  else          ! rcn < 0
    recn = file_size(filename)  !!!! file size (in bytes)
    recn = recn / leng + 1
    if (recn > 1)  findSize = .true.  ! we found the filesize to get the (correct) recn value
  endif 
else
  recn = 1  ! default
endif 

if (doSlice .and. .not. findSize)   recn = (recn - 1) * shp(3) + 1  
!^ because the input rcn was counted as the # of 3d matrices, but here 1 record is a slice only,
!  so recn needs to be recomputed.

if (recn==1) then
  open(10, file=filename, access='direct', recl=leng,status='replace')
else
  open(10, file=filename, access='direct', recl=leng)  
endif

if (doSlice) then
  do k = 1, shp(3)                
    write(10,rec=recn+k-1) data(:,:,k)
  enddo
else  
  write(10,rec=recn)  data 
endif
call flush(10)
close(10)
END SUBROUTINE

SUBROUTINE write_binfile_d3(filename, data, rcn) ! double, 3d array
character(len=*), intent(in) :: filename
real(DP),         intent(in) :: data(:,:,:)
integer(I4B)   :: shp(3), k
integer,optional, intent(in) :: rcn ! the record # to be written at, assuming record length as datasize
integer  :: recn, leng
logical  :: findSize, doSlice

shp = shape(data)
leng = product(shp) *2*I4
if (leng < RECL_MAX) THEN
  doSlice = .false.  ! can do 3d in 1 step
else
  leng    = shp(1)*shp(2) *2*I4
  doSlice = .true.   ! must do 1 slice @ a time to avoid write buffer overflow
endif  

findSize = .false.
if (present(rcn)) then
  if (rcn > 0) then
    recn = rcn
  else          ! rcn < 0
    recn = file_size(filename)  !!!! file size (in bytes)
    recn = recn / leng + 1
    if (recn > 1)  findSize = .true.  ! we found the filesize to get the (correct) recn value
  endif 
else
  recn = 1  ! default
endif 

if (doSlice .and. .not. findSize)   recn = (recn - 1) * shp(3) + 1  
!^ because the input rcn was counted as the # of 3d matrices, but here 1 record is a slice only,
!  so recn needs to be recomputed.

if (recn==1) then
  open(10, file=filename, access='direct', recl=leng,status='replace')
else
  open(10, file=filename, access='direct', recl=leng)  
endif

if (doSlice) then
  do k = 1, shp(3)                
    write(10,rec=recn+k-1) data(:,:,k)
  enddo
else  
  write(10,rec=recn)  data 
endif
call flush(10)
close(10)
END SUBROUTINE

SUBROUTINE write_binfile_c3(filename, data, rcn) ! complex, 3d array
character(len=*), intent(in) :: filename
complex(SPC),     intent(in) :: data(:,:,:)
integer(I4B)   :: shp(3), k
integer,optional, intent(in) :: rcn ! the record # to be written at, assuming record length as datasize
integer  :: recn, leng
logical  :: findSize, doSlice

shp = shape(data)
leng = product(shp) *2 *I4
if (leng < RECL_MAX) THEN
  doSlice = .false.  ! can do 3d in 1 step
else
  leng    = shp(1)*shp(2) *2 *I4
  doSlice = .true.   ! must do 1 slice @ a time to avoid write buffer overflow
endif  

findSize = .false.
if (present(rcn)) then
  if (rcn > 0) then
    recn = rcn
  else          ! rcn < 0
    recn = file_size(filename)  !!!! file size (in bytes)
    recn = recn / leng + 1
    if (recn > 1)  findSize = .true.  ! we found the filesize to get the (correct) recn value
  endif 
else
  recn = 1  ! default
endif 

if (doSlice .and. .not. findSize)   recn = (recn - 1) * shp(3) + 1  
!^ because the input rcn was counted as the # of 3d matrices, but here 1 record is a slice only,
!  so recn needs to be recomputed.

if (recn==1) then
  open(10, file=filename, access='direct', recl=leng,status='replace')
else
  open(10, file=filename, access='direct', recl=leng)  
endif

if (doSlice) then
  do k = 1, shp(3)                
    write(10,rec=recn+k-1) data(:,:,k)
  enddo
else  
  write(10,rec=recn)  data 
endif
call flush(10)
close(10)
END SUBROUTINE

!=============================
SUBROUTINE read_tex_i1(A, fname, n)
 integer,  allocatable,intent(out) :: A(:)
 character(len=*),     intent(in)  :: fname
 integer,  optional,   intent(in)  :: n
 integer               :: i, nr
 
 open(10, file=fname, form='formatted')
 ! determine # rows
 if (present(n)) then
   nr = n
 else
	 nr = 0
	 do 
	  read(10,*,end=99) i
	  nr = nr + 1
	 end do
	 99  continue
	 rewind(10)
 end if
 
 if (allocated(A)) then
   i = size(A)
   if (i /= nr) then
     deallocate(A)
	 allocate(A(nr))
   end if
 else
   allocate(A(nr))
 end if
 
 ! read in to A
 do i = 1, nr
   read(10,*) A(i)
 end do
 close(10)
END SUBROUTINE
!-------------
SUBROUTINE read_tex_r1(A, fname, n)
 real,     allocatable,intent(out) :: A(:)
 character(len=*),     intent(in)  :: fname
 integer,  optional,   intent(in)  :: n
 integer               :: nr, i
 real                  :: rtmp
 
 open(10, file=fname, form='formatted')
 ! determine # rows
 if (present(n)) then
   nr = n
 else
	 nr = 0
	 do 
	  read(10,*,end=99) rtmp
	  nr = nr + 1
	 end do
99   continue
	 rewind(10)
 end if
 
 if (allocated(A)) then
   i = size(A)
   if (i /= nr) then
	 deallocate(A)
	 allocate(A(nr))
   end if
 else
   allocate(A(nr))
 end if

 ! read in to A
 do i = 1, nr
   read(10,*) A(i)
 end do
 close(10)
END SUBROUTINE
!=============================
SUBROUTINE write_tex_char(fname, data, flg)
character(len=*), intent(in) :: fname
character(len=*), intent(in) :: data
integer,optional, intent(in) :: flg
integer    :: uni
logical    :: REP
uni = 10;  if (fname==''.or. fname=='n/a' .or. fname=='N/A') uni=6 ! stdout

REP = .false.
if (present(flg)) then
  if (flg <= 0)   REP = .true. 
endif
if (REP) then
  open(uni, file=fname, status='replace') 
else
  open(uni, file=fname, position='append')
endif
write(uni, *) trim(data)
call flush(uni)
close(uni)
END SUBROUTINE

SUBROUTINE write_tex_r0(fname, data, flg)
character(len=*), intent(in) :: fname
real, intent(in) :: data
integer,optional, intent(in) :: flg
logical    :: REP
integer    :: uni
uni = 10;  if (fname==''.or. fname=='n/a' .or. fname=='N/A') uni=6 ! stdout

REP = .false.
if (present(flg)) then
  if (flg <= 0)   REP = .true. 
endif
if (REP) then
  open(uni, file=fname, status='replace') 
else
  open(uni, file=fname, position='append')
endif
write(uni, *) data
call flush(uni)
close(uni)
END SUBROUTINE

SUBROUTINE write_tex_i0(fname, data, flg)
character(len=*), intent(in) :: fname
integer, intent(in) :: data
integer,optional, intent(in) :: flg
logical    :: REP
integer    :: uni
uni = 10;  if (fname==''.or. fname=='n/a' .or. fname=='N/A') uni=6 ! stdout

REP = .false.
if (present(flg)) then
  if (flg <= 0)   REP = .true. 
endif
if (REP) then
  open(uni, file=fname, status='replace') 
else
  open(uni, file=fname, position='append')
endif
write(uni, *) data
call flush(uni)
close(uni)
END SUBROUTINE

SUBROUTINE write_tex_l0(fname, data, flg)
character(len=*), intent(in) :: fname
logical, intent(in) :: data
integer,optional, intent(in) :: flg
logical    :: REP
integer    :: uni
uni = 10;  if (fname==''.or. fname=='n/a' .or. fname=='N/A') uni=6 ! stdout

REP = .false.
if (present(flg)) then
  if (flg <= 0)   REP = .true. 
endif
if (REP) then
  open(uni, file=fname, status='replace') 
else
  open(uni, file=fname, position='append')
endif
write(uni, *) data
call flush(uni)
close(uni)
END SUBROUTINE

SUBROUTINE write_tex_l10(fname, data, flg)
character(len=*), intent(in) :: fname
logical(1), intent(in) :: data
integer,optional, intent(in) :: flg
logical    :: REP
integer    :: uni
uni = 10;  if (fname==''.or. fname=='n/a' .or. fname=='N/A') uni=6 ! stdout

REP = .false.
if (present(flg)) then
  if (flg <= 0)   REP = .true. 
endif
if (REP) then
  open(uni, file=fname, status='replace') 
else
  open(uni, file=fname, position='append')
endif
write(uni, *) data
call flush(uni)
close(uni)
END SUBROUTINE

SUBROUTINE write_tex_b0(fname, data, flg)
character(len=*), intent(in) :: fname
integer(1), intent(in) :: data
integer,optional, intent(in) :: flg
logical    :: REP
integer    :: uni
uni = 10;  if (fname==''.or. fname=='n/a' .or. fname=='N/A') uni=6 ! stdout

REP = .false.
if (present(flg)) then
  if (flg <= 0)   REP = .true. 
endif
if (REP) then
  open(uni, file=fname, status='replace') 
else
  open(uni, file=fname, position='append')
endif
write(uni, *) data
call flush(uni)
close(uni)
END SUBROUTINE

SUBROUTINE write_tex_i1(fname, data, flg)
character(len=*), intent(in) :: fname
integer, intent(in) :: data(:)
integer,optional, intent(in) :: flg
logical    :: REP
integer    :: uni
uni = 10;  if (fname==''.or. fname=='n/a' .or. fname=='N/A') uni=6 ! stdout

REP = .false.
if (present(flg)) then
  if (flg <= 0)   REP = .true. 
endif
if (REP) then
  open(uni, file=fname, status='replace') 
else
  open(uni, file=fname, position='append')
endif
write(uni, '(100I5)') data
call flush(uni)
close(uni)
END SUBROUTINE

SUBROUTINE write_tex_r1(fname, data, flg)
character(len=*), intent(in) :: fname
real,   intent(in) :: data(:)
integer,optional, intent(in) :: flg
logical    :: REP
integer    :: uni
uni = 10;  if (fname==''.or. fname=='n/a' .or. fname=='N/A') uni=6 ! stdout

REP = .false.
if (present(flg)) then
  if (flg <= 0)   REP = .true. 
endif
if (REP) then
  open(uni, file=fname, status='replace') 
else
  open(uni, file=fname, position='append')
endif  
write(uni, '(100EN13.3)') data
call flush(uni)
close(uni)
END SUBROUTINE

SUBROUTINE write_tex_c1(fname, data, flg)
character(len=*), intent(in) :: fname
complex, intent(in) :: data(:)
integer, optional, intent(in) :: flg
logical    :: REP
integer    :: uni
uni = 10;  if (fname==''.or. fname=='n/a' .or. fname=='N/A') uni=6 ! stdout

REP = .false.
if (present(flg)) then
  if (flg <= 0)   REP = .true. 
endif
if (REP) then
  open(uni, file=fname, status='replace') 
else
  open(uni, file=fname, position='append')
endif
write(uni, '( 100("(",EN13.3,",",EN13.3,") ") )') data
call flush(uni)
close(uni)
END SUBROUTINE

SUBROUTINE write_tex_z1(fname, data, flg)
character(len=*), intent(in) :: fname
complex(DPC), intent(in) :: data(:)
integer, optional, intent(in) :: flg
logical    :: REP
integer    :: uni
uni = 10;  if (fname==''.or. fname=='n/a' .or. fname=='N/A') uni=6 ! stdout

REP = .false.
if (present(flg)) then
  if (flg <= 0)   REP = .true. 
endif
if (REP) then
  open(uni, file=fname, status='replace') 
else
  open(uni, file=fname, position='append')
endif
write(uni, '( 100("(",EN13.3,",",EN13.3,") ") )') data
call flush(uni)
close(uni)
END SUBROUTINE

SUBROUTINE write_tex_l11(fname, data, flg)
character(len=*), intent(in) :: fname
logical(1), intent(in) :: data(:)
integer,optional, intent(in) :: flg
logical    :: REP
integer    :: uni
uni = 10;  if (fname==''.or. fname=='n/a' .or. fname=='N/A') uni=6 ! stdout

REP = .false.
if (present(flg)) then
  if (flg <= 0)   REP = .true. 
endif
if (REP) then
  open(uni, file=fname, status='replace') 
else
  open(uni, file=fname, position='append')
endif
write(uni, '(100L4)') data
call flush(uni)
close(uni)
END SUBROUTINE

SUBROUTINE write_tex_l1(fname, data, flg)
character(len=*), intent(in) :: fname
logical, intent(in) :: data(:)
integer,optional, intent(in) :: flg
logical    :: REP
integer    :: uni
uni = 10;  if (fname==''.or. fname=='n/a' .or. fname=='N/A') uni=6 ! stdout

REP = .false.
if (present(flg)) then
  if (flg <= 0)   REP = .true. 
endif
if (REP) then
  open(uni, file=fname, status='replace') 
else
  open(uni, file=fname, position='append')
endif
write(uni, '(100L4)') data
call flush(uni)
close(uni)
END SUBROUTINE

SUBROUTINE write_tex_b1(fname, data, flg)
character(len=*), intent(in) :: fname
integer(1), intent(in) :: data(:)
integer,optional, intent(in) :: flg
logical    :: REP
integer    :: uni
uni = 10;  if (fname==''.or. fname=='n/a' .or. fname=='N/A') uni=6 ! stdout

REP = .false.
if (present(flg)) then
  if (flg <= 0)   REP = .true. 
endif
if (REP) then
  open(uni, file=fname, status='replace') 
else
  open(uni, file=fname, position='append')
endif
write(uni, '(100I4)') data
call flush(uni)
close(uni)
END SUBROUTINE

SUBROUTINE write_tex_i2(fname, data, flg)
character(len=*), intent(in) :: fname
integer, intent(in) :: data(:,:)
integer, optional, intent(in) :: flg
integer  :: fstdim, m, n
logical  :: REP
integer    :: uni
uni = 10;  if (fname==''.or. fname=='n/a' .or. fname=='N/A') uni=6 ! stdout

REP = .false.;   fstdim = 2   ! default
if (present(flg)) then
  fstdim = flg
  if (flg <= 0)   then
    REP = .true.; fstdim = -flg
  endif
endif

if (REP) then
  open(uni, file=fname, status='replace') 
else
  open(uni, file=fname, position='append')
end if
if ( fstdim == 2 ) then
  do m = 1, size(data,1)
    write(uni, '(100I5)') (data(m,n),n=1,size(data,2))
  enddo  
else
  do n = 1, size(data,2)
    write(uni, '(100I5)') (data(m,n),m=1,size(data,1))
  enddo  
endif
call flush(uni);  close(uni)
END SUBROUTINE

SUBROUTINE write_tex_r2(fname, data, flg)
character(len=*), intent(in) :: fname
real,   intent(in) :: data(:,:)
integer,optional, intent(in) :: flg
integer  :: fstdim, m, n
logical  :: REP
integer    :: uni
uni = 10;  if (fname==''.or. fname=='n/a' .or. fname=='N/A') uni=6 ! stdout

REP = .false.;   fstdim = 2   ! default
if (present(flg)) then
  fstdim = flg
  if (flg <= 0)   then
    REP = .true.;  fstdim = -flg
  endif
endif

if (REP) then
  open(uni, file=fname, status='replace') 
else
  open(uni, file=fname, position='append')
endif
if ( fstdim == 2 ) then
  do m = 1, size(data,1)
    write(uni, '(100EN13.3)') (data(m,n),n=1,size(data,2))
  enddo  
else
  do n = 1, size(data,2)
    write(uni, '(100EN13.3)') (data(m,n),m=1,size(data,1))
  enddo  
endif
 call flush(uni); close(uni)
END SUBROUTINE

SUBROUTINE write_tex_c2(fname, data, flg)
character(len=*), intent(in) :: fname
complex,intent(in) :: data(:,:)
integer,optional, intent(in) :: flg
integer  :: fstdim, m, n
logical  :: REP
integer    :: uni
uni = 10;  if (fname==''.or. fname=='n/a' .or. fname=='N/A') uni=6 ! stdout

REP = .false.;   fstdim = 2   ! default
if (present(flg)) then
  fstdim = flg
  if (flg <= 0)   then
    REP = .true.;  fstdim = -flg
  endif
endif
if (REP) then
  open(uni, file=fname, status='replace') 
else
  open(uni, file=fname, position='append')
endif

if ( fstdim == 2) then
  do m = 1, size(data,1)
    write(uni, '( 100("(",EN13.3,",",EN13.3,") ") )') (data(m,n),n=1,size(data,2))
  enddo  
else
  do n = 1, size(data,2)
    write(uni, '( 100("(",EN13.3,",",EN13.3,") ") )') (data(m,n),m=1,size(data,1))
  enddo  
endif
 call flush(uni); close(uni)
END SUBROUTINE

SUBROUTINE write_tex_z2(fname, data, flg)
character(len=*), intent(in) :: fname
complex(DPC),intent(in) :: data(:,:)
integer,optional, intent(in) :: flg
integer  :: fstdim, m, n
logical  :: REP
REP = .false.;   fstdim = 2   ! default
if (present(flg)) then
  fstdim = flg
  if (flg <= 0)   then
    REP = .true.;  fstdim = -flg
  endif
endif
if (REP) then
  open(10, file=fname, status='replace') 
else
  open(10, file=fname, position='append')
endif

if ( fstdim == 2) then
  do m = 1, size(data,1)
    write(10, '( 100("(",EN13.3,",",EN13.3,") ") )') (data(m,n),n=1,size(data,2))
  enddo  
else
  do n = 1, size(data,2)
    write(10, '(100EN13.3)') (data(m,n),m=1,size(data,1))
  enddo  
endif
call flush(10)
close(10)
END SUBROUTINE
!-----------------------------

END MODULE seism_uti_module


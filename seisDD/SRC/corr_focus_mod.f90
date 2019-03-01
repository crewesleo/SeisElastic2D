MODULE corr_focus_mod

  USE myMath_module;    USE seism_uti_module
 ! USE myString_mod

  IMPLICIT NONE

  ! PRIVATE  :: 

CONTAINS
!==================================
function SQRT_to_odd1(n2) result(od)
integer  :: n2, od
integer  :: j, k
! Gieven input n2, find closest od >= 3, such that od^2 == n2

j  = floor( (sqrt(real(n2))-1.)/2. )
k  = j+1

j  = 2*j+1  ! odd
k  = 2*k+1  ! odd

if (n2 - j**2 <= k**2 - n2) then
  od = j
else
  od = k
end if
od = max(1, od)

end function
!-----------
subroutine SQRT_to_odd2(od, er, n2)
integer, intent(out) :: od(2), er
integer, intent(in)  :: n2
integer  :: j(4), k(4), e(4), i, loc(1), resi

! Given an input n2, find the best od(2)>=3, such that \sum od(:)^2 == n2

j = (/ -1, 2, 1, 0 /) +  floor( (sqrt(real(n2)/2.)-1.)/2. )

j = 2*j + 1 ! odd
e = 10000   ! init big
od= 0

do i = 1,4
  resi = n2 - j(i)**2 
  if (j(i) <= 1  .OR. resi <= 0)  cycle
  k(i) = SQRT_to_odd1(resi) 
  e(i) = abs(resi - k(i)**2)
end do

loc = minloc(e)
od(1) = j( loc(1) )
od(2) = k( loc(1) )
er    = e( loc(1) )
end subroutine
!------------

subroutine SQRT_to_odd3(od, er, n2)
integer, intent(out) :: od(3), er
integer, intent(in)  :: n2
integer  :: j(4), e(4), i, loc(1), resi

! Given an input n2, find the best od(3), such that \sum od(:)^2 == n2

j = (/ -1, 2, 1, 0 /) +  floor( (sqrt(real(n2)/3.)-1.)/2. )

j = 2*j + 1 ! odd
e = 10000   ! init big
od= 0

do i = 1,4
  resi = n2 - j(i)**2
  if (j(i) <=1 .OR. resi <=0)   cycle
  call SQRT_to_odd2(od(2:3), e(i), resi)
end do

loc = minloc(e)
er  = e ( loc(1) )
od(1) = j(loc(1) )

if (er==10000)  return  ! all are invalid

if (loc(1) /= 4) then
  call SQRT_to_odd2(od(2:3), er, n2-od(1)**2)
end if

end subroutine
!------------
subroutine SQRT_to_odd4(od, er, n2)
integer, intent(out) :: od(4), er
integer, intent(in)  :: n2
integer  :: j(4), e(4), i, loc(1), resi

! Given an input n2, find the best od(4), such that \sum od(:)^2 == n2

j = (/ -1, 2, 1, 0 /) +  floor( (sqrt(real(n2)/4.)-1.)/2. )

j = 2*j + 1 ! odd
e = 10000   ! init big
od= 0

do i = 1,4
  resi = n2 - j(i)**2
  if (j(i) <=1 .OR. resi <=0)   cycle
  call SQRT_to_odd3(od(2:4), e(i), resi)
end do

loc = minloc(e)
er  = e ( loc(1) )
od(1) = j(loc(1) )

if (er==10000)  return  ! all are invalid

if (loc(1) /= 4) then
  call SQRT_to_odd3(od(2:4), er, n2-od(1)**2)
end if

end subroutine
!-------------
subroutine Bspline_wids(wids, M)
! rect width M, to be equated w/ a B-spline with 4 narrower widths --> their spatial var. are eq
integer, intent(out):: wids(4)
integer, intent(in) :: M   ! var = (M^2 - 1)/12

integer  :: er
if (M <=  6) then
   wids = 3
   return
end if
call SQRT_to_odd4(wids, er, M*M + 3)

end subroutine
!-------------

subroutine make_imp_resp_Bspline(V, wids)
! B-spline with 4 narrower widths.  Find the impulse response in V
real, allocatable, intent(out):: V(:)   
integer, intent(in)  :: wids(4)

integer  :: M2, M
real(8), allocatable :: tmp1(:), tmp2(:)

M  = sum(wids) - (size(wids) - 1)  ! because of keeping conv'ing
M2 = M/2

if (allocated(V))    deallocate(V)
allocate(V(M), tmp1(M), tmp2(M))
V       =  0.
V(M2+1) =  1.
call Bspline1D_mirr_conv_S(V, wids, tmp1, tmp2)

deallocate(tmp1, tmp2)

end subroutine

!=============

subroutine index_intersec(hsi,tsi, Hbi,Tbi, hs,ts, Hb,Tb)
! a short vector whose index range on the domain of a big vector is
! supposedly [hs,ts].  The big vector's valid index range is [Hb,Tb].
! Intersect them, so now the [hsi,tsi] of the small vec resides on 
! [Hbi, Tbi] of the big one.

integer, intent(out):: hsi,tsi, Hbi,Tbi
integer, intent(in) :: hs, ts,  Hb, Tb
integer  :: L

L = ts-hs+1   ! small vec's length

if (hs>=Hb) then
   hsi = 1
   Hbi = hs
else
   hsi = 1+Hb-hs
   Hbi = Hb
end if

if (ts<=Tb) then
   tsi = L
   Tbi = ts
else
   tsi = L - (ts-Tb)
   Tbi = Tb
end if

end subroutine
!-------------

subroutine corr_cent(y, u, W, lag12, HEAD)
! y = Xcorr(u, W)
! u: N x 1
! W: M x 1,  M << N, M odd, M's reference point is its center, if HEAD == <0>
!                     else, the reference pnt is W(1)
!
! W sliding across u, produce the corr output y
! lag12  2 x 1,  the lag range of W's reference pnt wrt u
! y:  lag12(2)-lag12(1)+1 x 1
!===  All vectors' memory assigned ====

real,   intent(out):: y(:)
real,   intent(in) :: u(:), W(:)
integer,intent(in) :: lag12(2)
integer,intent(in),optional :: HEAD

integer :: i, hsi,tsi, Hbi,Tbi, N, M, M2, headx ! local var HEAD
! integer, save :: callcounter = 0

headx = 0;    if(present(HEAD))  headx = HEAD

N = size(u)
M = size(W)
if (headx==0) then
  M2 = M/2
  do i = lag12(1), lag12(2)
      call index_intersec(hsi,tsi, Hbi,Tbi, i-M2, i+M2, 1, N)
      y(i-lag12(1)+1) = sum(W(hsi:tsi)*u(Hbi:Tbi))
  end do
else
  do i = lag12(1), lag12(2)
      call index_intersec(hsi,tsi, Hbi,Tbi, i, i+M-1, 1, N)
      y(i-lag12(1)+1) = sum(W(hsi:tsi)*u(Hbi:Tbi))
  end do
end if
!if (callcounter == 0) then
!   call write_binfile('y.bin', y)  !!!!!!!!!
!   call write_binfile('u.bin', u)  !!!!!!!!!
!   call write_binfile('W.bin', W)  !!!!!!!!!
!end if
!callcounter = callcounter + 1

end subroutine
!-------------
subroutine squashNega(f,g, x, alph, MODE)
! idea is: reduce the influence of the negative part of x
! x:  vector
! alph:  1x1
! f:  same size as x
! g:  gradient elem-wise, \partial f_k / \partial x_k, so same size as x
!
! MODE == 0:
! map x \in [-1, 1] to f \in [0, 1]  ~ exp(x)
!     x is normalized, so \in [-1, 1]
!
! MODE == 1, 2:
! f(x) = exp(alph.*x), when x<=0.  When x>0:
! 1:
!     1 + alph.*x      f, f' conti
! 2:
!     alph^2/2 .*(x+ 1/alph).^2 + 1/2,   f, f', f" conti

real, intent(out) :: f(:), g(:)
real, intent(in)  :: x(:), alph
integer,intent(in):: MODE
real  :: bet, exa

select case (MODE)
case (0)
   exa = exp(alph)
   bet = 1./(exa - 1./exa)
   g   = exp(alph*x)
   f   = bet * (g - 1./exa)
   g   = (alph*bet) * g
case (1)
   where (x<0)
      f = exp(alph*x)
      g = alph*f
   elsewhere
      f = 1 + alph*x
      g = alph
   end where
case (2)
   bet = alph*alph
   where (x<0)
      f = exp(alph*x)
      g = alph*f
   elsewhere
      f = 0.5*(1. + bet* (x + 1/alph)**2)
      g =       bet* (x + 1/alph)
   end where
end select

end subroutine
!-------------

subroutine localMx(pk, y, L2, Boff, thresh, pktmp, itmp)
! y    N x 1  waveform
! pk   N x 1  logical, marking the local peaks of y within a
!             [-L2, L2] neighborhood
! Boff 1 x 1  the distance within Boff near the boundary is turned off
! thresh 1x1  peaks must >= thresh
! pktmp  Nx1  same size as pk

logical, intent(out)   :: pk(:)
logical, intent(inout) :: pktmp(:)
real,    intent(in)    :: y(:), thresh
integer, intent(in)    :: L2,   Boff
integer, intent(inout) :: itmp(:)

integer  :: N, i, M, M2, lf, rt, id(1)

N     =  size(y)
pk    = .false.
pktmp = .false.
pktmp(2:N-1) = y(2:N-1) > y(1:N-2) .and. y(2:N-1) > y(3:N) &
               .and. y(2:N-1) > thresh

pktmp(1:Boff)      = .false.
pktmp(N-Boff+1:N)  = .false.  ! preventing peaks near the boundaries

M     = count(pktmp)
itmp  = pack( (/(i, i=1,N) /), pktmp )  ! idx to peaks

do i = 1, M
  lf = max(1, itmp(i)-L2)
  rt = min(N, itmp(i)+L2)
  id = maxloc(y(lf:rt), pktmp(lf:rt))
  if (id(1) == itmp(i)-lf+1)  pk(itmp(i)) = .true.
end do

end subroutine
!-------------

subroutine envelop_mx(y, pk, itmp, GauSig2)
! y   N x 1,  waveform, in-place form an envelop over the peaks of y
! pk  N x 1,  logical, marking the peaks of y
! GauSig2 1 x 1   real,  Gaussian's sigma^2.  The Gaussian is the yenv outside the itmp range
!     if GauSig2 <= 0, then outside is unmodified (0 most likely)
!     if ~present, then outside is constant
real,    intent(inout) :: y(:)
logical, intent(in)    :: pk(:)
integer, intent(inout) :: itmp(:)
real,    intent(in), optional :: GauSig2

integer  :: N, i, M, dx, x1, x2
real     :: y1, y2

integer, save :: callcounter = 0

N     =  size(y)
itmp  =  pack( (/ (i, i=1,N) /), pk )  ! idx to peaks
M     =  count(pk)

if (callcounter==0) then
   call write_binfile('y--0.bin', y)  !!!!!!!!!
else if (callcounter==1) then
   call write_binfile('y--1.bin', y)  !!!!!!!!!
end if

if (present(GauSig2))  then
   if (GauSig2 > 0) then
      y(1:itmp(1)) = y(itmp(1))*exp( (/ (i, i=1-itmp(1),0) /)**2 /(-2.*GauSig2) )
      y(itmp(M):N) = y(itmp(M))*exp( (/ (i, i=0,N-itmp(M)) /)**2 /(-2.*GauSig2) )
   end if  ! otherwise don't modify
else
   y(1:itmp(1)) = y(itmp(1))
   y(itmp(M):N) = y(itmp(M))
end if

do i = 1, M-1
   x1 = itmp(i)
   x2 = itmp(i+1)
   dx = x2 - x1
   y1 = y(x1)
   y2 = y(x2)
   y(x1:x2) = y2 + (y1-y2)/2. * (1.+ &
     cos( (/ (N, N=0, x2-x1) /)/(real(dx)/PI) ) )
end do

if (callcounter==0) then
   call write_binfile('yenv0.bin', y)  !!!!!!!!!
else if (callcounter==1) then
   call write_binfile('yenv1.bin', y)  !!!!!!!!!
end if
callcounter = callcounter + 1

end subroutine
!-------------

subroutine corr_focusing1D(J, g, u,upks,V,Vwids, d,W0,W1,Wsig2, alph, lag12,&
                unrm,gchnk,   dtprd,fdtp, gdtp, simi,simF, dnrm, dtmp, &
                itmp,pk,pktmp,Vtmp, tmpD1,tmpD2)
! J = \Sum_k { \int W f(Xcorr^k, alph) / \int f(Xcorr^k, alph) }
!          k: peak index in u
!          Xcorr^k  corr(d, u^k), where d: data trace; u^k: kth u's chunk
!
! g:  grad J / u,   N x 1
!                     f(., alpha, MODE)   squashing-negativ nonlinearity
! 
! u    N x 1    predicted trace
! upks Npk x 1  peak locations in u
!               around each peak [-M/2, M/2] is a u's chunk
!      NOTE: u's peaks should be M/2 away from the trace beg & end!      
!
! V    M x 1    the narrow weight profile for a /u chunk/
! Vwids Nvwids x 1  rect wids, keep convolving => the B-spline V
!
! d    N x 1    data trace
! W0   lag(2)-lag(1)+1 x 1   Gaussian weigtht profile over xcorr lags
! W1   ditto                 specific weights
! 
! Wsig2 1x 1   for overall weight profile (e.g., Gaussian, B-spline) over Xcorr shifts,
!              encouraging xcorr mass @ 0.
! lag12  2 x 1  the Xcorr is limied to this range of *relative* lags (relative means,
!     assuming u chunk's center @ t=0, where does it shift around on d's trace?)

! MODE: norm. dot-prod followed by a nega-squashing nonlinearity of mode
! 0:   [e(a x)
! 1:   [e(a x), x-b]
! 2:   [e(a x), c(x-b)^2+d], then \int w/ weight

real,   intent(out)  :: J,    g(:)
real,   intent(in)   :: u(:), d(:), V(:),    W0(:),   Wsig2, alph
integer,intent(in)   :: upks(:),  Vwids(:),  lag12(2)
real   ,intent(inout):: unrm(:),  gchnk(:),  W1(:),   dtprd(:,:), &
                        fdtp(:,:),gdtp(:,:),             & ! all pre-allocated memory
                        simi(:,:),simF(:,:), dnrm(:), dtmp(:), Vtmp(:)     
! unrm:    u chunks' norm (Npk #)
! gchnk:   g_chunk  M x 1
! dtprd:   normlz. dot-prod (d, u_chunks)  diff(lag12)+1 x Npk
! fdtp, gdtp: funct & grad of passing dtprd through a nega-squashing nonlinearity
!       same size as dtprd
! simi: init. similarity weights, between d's elements and u's chunks, N x Npk
! simF: Final ... (after double-normalizing)
!       From simF/simi, we know the gain weights.
! dnrm: d's local sliding-window's norm, N x 1
! dtmp: N x 1
! Vtmp: length(M), temp for a u chunk

real(8),intent(inout) :: tmpD1(:), tmpD2(:)  ! for Bspline V (*) ..., length(M)
integer,intent(inout) :: itmp(:) ! peak indices 
logical,intent(inout) :: pk(:), pktmp(:) 

integer ::  N, Npk, i, k, M, M2, lf,rt, da,db,Va,Vb, lg1,lg2,  dh,dt,Wh,Wt  
            ! da,db,Va,Vb represent beg & end indices on d and V  (after intersection)
            ! dh,dt,Wh,Wt  ...      beg & end indices on d and W  (after intersection)

real    ::  Jk, E, eps, thresh

! dtprd & simi: d-trace & u chunk

N     =  size(d)
M     =  size(V);   M2 = M/2
Npk   =  size(upks)

! compute d_nrm:
dnrm  = d*d
call Bspline_conv(dnrm, Vwids, tmpD1, tmpD2)
dnrm  = sqrt(dnrm)  ! local d_norm
! call write_binfile('dnrm.bin',dnrm)  !!!!!!!!
! compute u_nrm:

eps  = 1.e-4
simi = 0.
simF = 0.
dtprd= 0.

do k = 1, Npk
  lf = upks(k)-M2           ! may outside [1,N]
  rt = upks(k)+M2
  call index_intersec(Va,Vb, da,db, lf, rt,  1,N)

  lg1= lag12(1) + upks(k)   ! may outside [1,N]
  lg2= lag12(2) + upks(k)
  call index_intersec(Wh,Wt, dh,dt, lg1,lg2, 1,N)

  Vtmp    = 0.
  Vtmp(Va:Vb) = V(Va:Vb)  * u(da : db)               ! weighted u chunk
  unrm(k)     = sqrt(sum(Vtmp(Va:Vb) * u(da : db) )) ! chunk norm, weighted
  call corr_cent(dtprd(Wh:Wt,k), d, Vtmp, (/dh, dt /) )
!  dtprd(Wh:Wt,k) = dtprd(Wh:Wt,k)/(unrm(k)+eps)/(dnrm(dh:dt) + eps)

  call squashNega(fdtp(Wh:Wt,k),gdtp(Wh:Wt,k), dtprd(Wh:Wt,k), alph, 2) ! MODE 2
  simi(dh:dt,k) = fdtp(Wh:Wt,k)
end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#if 0 

  call write_binfile('dtprd.bin',dtprd)  !!!!!!!!!
  call write_binfile('fdtp.bin', fdtp)  !!!!!!!!!
  call write_binfile('simi0.bin', simi)  !!!!!!!!!

! threshold of peaks
thresh = 0.2*minval(maxval(simi, 1))

lf = 0.5*M2  ! a temp int
do k = 1, Npk
  call localMx(pk, simi(:,k), lf, max(M2,maxval(abs(lag12))), thresh, pktmp, itmp)
  call envelop_mx( simi(:,k), pk, itmp, real(lag12(1)**2) )  ! should prod Gaussian 1st, then mark peaks

  simi(:,k) = simi(:,k) * exp( (/ (i-upks(k), i=1,N) /)**2 /(-2.*Wsig2) )
end do  ! init similarity envelope

call write_binfile('simi.bin',simi)  !!!!!!!!!

simF = simi
! compute simi as attribution weight 
dtmp = sum(simF, 2) + eps
do k = 1, Npk
   simF(:,k) = simF(:,k) / dtmp  ! 1) normalize over u_chunks
end do
call write_binfile('simF1.bin',simF)  !!!!!!!!!

do k = 1, Npk
   simF(:,k) = simF(:,k) / (sum(simF(:,k))+eps)  ! 2) normalize over t
end do
call write_binfile('simF2.bin',simF)  !!!!!!!!!


dtmp = sum(simF, 2) + eps
do k = 1, Npk
   simF(:,k) = simF(:,k) / dtmp  ! 3) normalize over u_chunks
end do
  call write_binfile('simF3.bin',simF)  !!!!!!!!!

#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

J    = 0.
g    = 0.
do k = 1, Npk
   lf= upks(k)-M2
   rt= upks(k)+M2
   call index_intersec(Va,Vb, da,db, lf, rt, 1,N)

  lg1= lag12(1) + upks(k)
  lg2= lag12(2) + upks(k)
  call index_intersec(Wh,Wt, dh,dt, lg1,lg2, 1,N)


   W1 = W0  !* simF(lg1:lg2,k)/(simi(lg1:lg2,k) + eps)

   E = sum(fdtp(Wh:Wt,k)) + eps
   Jk= sum(W1(Wh:Wt) * fdtp(Wh:Wt,k)) / E

   gchnk = 0.

!   call corr_cent(gchnk,d,(W1(Wh:Wt)-Jk)*gdtp(Wh:Wt,k)/(dnrm(dh:dt)+eps),dh+(/-M2,M2/),1)
!   gchnk(Va:Vb) = gchnk(Va:Vb) -sum((W1(Wh:Wt)-Jk)*gdtp(Wh:Wt,k)*dtprd(Wh:Wt,k))/(unrm(k)+eps) *u(da:db)
!   gchnk(Va:Vb) = gchnk(Va:Vb) * V(Va:Vb) / ((unrm(k)+eps) * E)

   call corr_cent(gchnk,d,(W1(Wh:Wt)-Jk)*gdtp(Wh:Wt,k), dh+(/-M2,M2/), 1)
   gchnk(Va:Vb) = gchnk(Va:Vb) * V(Va:Vb) / E

   g(da:db) = g(da:db) + gchnk(Va:Vb)

   J     = J + Jk
end do
! call write_binfile('grd.bin', g)  !!!!!!!!!

end subroutine
!-------------


subroutine corr_focusing_Adj_Src(AdjS1, misF,  datP1, datO1, f0, dt, Terr_max)
! Output:
! AdjS   nt x ng   adjoint source, or the negative-gradient direction to minimize the misfit
! misF   scalar    criterion function (of negative value)
!
! Input:
! datP   nt x ng   predicted data
! datO   nt x ng   observed data
! f0               peak freq. (Hz) 
! dt               delta t (s)
! Terr_max         max time error (s)

real,  intent(out) :: AdjS1(:,:),  misF
real,  intent(in)  :: datP1(:,:), datO1(:,:), f0, dt, Terr_max
!--------------------

real               :: J,  alph, Wsig2, thresh, stdV
integer            :: Vwids(4), lag12(2)
real,   allocatable:: V(:),  unrm(:), gchnk(:), dtprd(:,:),&
                      fdtp(:,:),gdtp(:,:),   W0(:), W1(:),        &
                      simi(:,:),simF(:,:),   dnrm(:),  dtmp(:), Vtmp(:),  &
                      datP(:,:),datO(:,:), AdjS(:,:)
integer,allocatable:: upks(:),  itmp(:) 
real(8),allocatable:: tmpD1(:), tmpD2(:) 
logical,allocatable:: pk(:),    pktmp(:) 

integer :: nt1, nt, ng, ig, i, M2, M, L2, L, Npk, Npkmx

nt1   = size(datP1,1)
ng    = size(datP1,2)
nt   = 3*nt1

allocate( datP(nt,ng), datO(nt,ng), AdjS(nt,ng) )
do i = 1, nt1
   datP(i+nt1,:)   = datP1(i,:);    datO(i+nt1,:) = datO1(i,:)
   datP(i,:) = datP1(nt1+1-i,:);    datO(i,:) = datO1(nt1+1-i,:)  ! mirror reflec
   datP(i+2*nt1,:) = datP(i,:);     datO(i+2*nt1,:) = datO(i,:)
end do

L2    = NINT(Terr_max/dt);  L = 2*L2+1  ! Xcorr lag range:
lag12 = (/ -L2, L2 /)

M     = NINT(2/f0/dt);      M2 = M/2;   M = 2*M2+1  ! wavelet chunk length ~ 2 T, odd
call Bspline_wids(Vwids, M)
call make_imp_resp_Bspline(V, Vwids)
M     = size(V);              M2 = M/2
!stdV  = sqrt( ( sum(Vwids**2)-size(Vwids) )/12. )  
stdV  = 0.5/f0/dt  ! spacing, half period
Npk   = NINT((Nt1-2*M2)/(stdV)) + 1

allocate(upks(Npk), W0(L), W1(L))

do i = 1,Npk
  upks(i) = Nt1+ 1+M2 + NINT( (i-1)*(1.5 *stdV) )
end do

Wsig2 = (L2*2)**2        ! var for W weight profile over Xcorr shifts
W0    = exp( (/ (i, i= -L2, L2) /)**2 /(-2.*Wsig2) )

! write(*,*) 'M, Vwids', M, Vwids; call flush(6)
! call write_binfile('V.bin', V)   !!!!!!!!!!!!!

alph  = 8.
Npkmx = nt1/2 ! # local maxima

allocate(pk(Nt),    pktmp(Nt), itmp(Npkmx),  &
         gchnk(M),  dnrm(Nt),  &
         dtmp(Nt),  Vtmp(M),   tmpD1(Nt), tmpD2(Nt), &
         unrm(Npk),   dtprd(L,Npk), &
         fdtp(L,Npk), gdtp(L,Npk),  &
         simi(Nt,Npk),simF(Nt,Npk)  )

misF = 0.
AdjS = 0.
do ig = 1, ng
    call corr_focusing1D(J, AdjS(:,ig), datP(:,ig),upks,V,Vwids, datO(:,ig), &
                        W0, W1, Wsig2,  alph, lag12, &
                        unrm,gchnk,   dtprd,fdtp, gdtp, simi,simF,dnrm,dtmp, &
                        itmp,pk,pktmp,Vtmp, tmpD1,tmpD2)
    misF = misF -J
end do
Adjs1 = AdjS(nt1+1: 2*nt1,:)

deallocate(pk,     pktmp, itmp, &
          gchnk,   W0,W1, dnrm, &
          dtmp,    V,Vtmp,tmpD1, tmpD2, &
          upks,    &
          unrm, dtprd,  &
          fdtp, gdtp,   &
          simi, simF    )

deallocate(datP, datO, Adjs)

end subroutine
!-------------

END MODULE


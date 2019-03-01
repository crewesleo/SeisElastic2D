!! main subroutines for misfit_adjoint_lib and preprocess_lib
!! created by Yanhua O.Yuan (yanhuay@princeton.edu)

!----------------------------------------------------------------------
subroutine detrend(x,n)
    use constants
    implicit none

    integer, intent(in) :: n
    real(kind=CUSTOM_REAL), dimension(*) :: x

    real(kind=CUSTOM_REAL) :: ds1,ds2,dan,davei,davex,dslope,dai
    integer :: i, an

    an = n
    dan=n
    ds1=0
    ds2=0

    do i=1,n
    ds1 = ds1+ x(i)
    ds2 = ds2 + ds1
    enddo
    davei = 0.5 * (dan+1.0)
    davex = ds1/dan
    dslope = -12.0*(ds2-davei*ds1)/(dan*(dan*dan-1.0))
    do i=1,n
    dai = i-1
    x(i) = x(i)- davex - dslope*(dai-davei)
    enddo

end subroutine detrend
!------------------------------------------------------------------------
subroutine xcorr_calc(dat1,dat2,npts,i1,i2,ishift,dlnA,cc_max)
    ! SAC libarary -- to get optimal shift (not scaled by time interval) 
    ! corresponding to maximal cross-correlation coefficients
    ! cross correlation time T(dat1-dat2) within window i1 and i2
    ! finxed the window for dat2, shift the window for dat1

    use constants
    implicit none

    real(kind=CUSTOM_REAL), dimension(*), intent(in) :: dat1,dat2
    integer, intent(in) :: npts, i1, i2

    ! outputs:
    ! ishift = index lag (d-s) for max cross correlation
    ! cc_max = maximum of cross correlation (normalised by sqrt(synthetic*data))
    integer, intent(out) :: ishift
    real(kind=CUSTOM_REAL), intent(out) :: cc_max,dlnA

    ! local variables
    integer :: nlen
    integer :: i_left, i_right, i, j, id_left, id_right
    real(kind=CUSTOM_REAL) :: cc, norm, norm_s

    ! initialise shift and cross correlation to zero
    ishift = 0
    dlnA = 0.0
    cc_max = 0.0

    ! print*,'ishift,dlnA,cc_max :',ishift,dlnA,cc_max

    if (i1.lt.1 .or. i1.gt.i2 .or. i2.gt.npts) then
        write(*,*) 'Error with window limits: i1, i2, npts ', i1, i2, npts
        return
    endif

    ! length of window (number of points, including ends)
    nlen = i2 - i1 + 1

    ! power of synthetic signal in window
    norm_s = sqrt(sum(dat2(i1:i2)*dat2(i1:i2)))

    ! left and right limits of index (time) shift search
    ! NOTE: This looks OUTSIDE the time window of interest to compute TSHIFT and CC.
    !       How far to look outside, in theory, should be another parameter.
    !       However, it does not matter as much if the data and synthetics are
    !          zeroed outside the windows, as currently done in calc_criteria.
    i_left = -1*int(nlen/2.0)
    i_right = int(nlen/2.0)
    !!  i_left = -nlen
    !!  i_right = nlen

    ! i is the index to shift to be applied to DATA (d)
    do i = i_left, i_right

    ! normalization factor varies as you take different windows of d
    id_left = max(1,i1+i)      ! left index for data window
    id_right = min(npts,i2+i)  ! right index for data window
    norm = norm_s * sqrt(sum(dat1(id_left:id_right)*(dat1(id_left:id_right))))

    ! cc as a function of i
    cc = 0.0
    do j = i1, i2   ! loop over full window length
    if((j+i).ge.1 .and. (j+i).le.npts) cc = cc + dat2(j)*dat1(j+i)  ! d is shifted by i
    enddo
    !print*,'norm,norm_s :',norm,norm_s

    ! normalized by norm of data
    if(norm > SMALL_VAL ) cc = cc/norm

    ! keeping cc-max only
    if (cc .gt. cc_max) then
        cc_max = cc
        ishift = i
        !print*,'i,cc,cc_max :',i,cc,cc_max
    endif
    enddo

    dlnA = 0.5 * log( sum(dat1(i1:i2) * dat1(i1:i2)) / sum(dat2(i1:i2) * dat2(i1:i2)) )

    !print*,'ishift,dlnA,cc_max :',ishift,dlnA,cc_max
    ! EXAMPLE: consider the following indexing:
    ! Two records are from 1 to 100, window is i1=20 to i2=41.
    !    --> nlen = 22, i_left = -11, i_right = 11
    !    i   i1+i   i2+i  id_left  id_right
    !  -11     9     30      9        30
    !   -5    15     36     15        36
    !    0    20     41     20        41
    !    5    25     46     25        46
    !   10    31     52     31        52

end subroutine xcorr_calc
!------------------------------------------------------------------------

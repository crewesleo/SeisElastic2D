!! main subroutines for preprocessing
!! created by Yanhua O. Yuan ( yanhuay@princeton.edu)

subroutine preprocess(d,NSTEP,deltat,t0,f0,&
        is_laplace,dis_sr,S_x,S_t, &
        is_window,window_type,V_slow,V_fast,&
        mute_near,offset_near,mute_far,offset_far,&
        Wscale,&
        istart,iend)

    !! preprocessing library
    use constants
    implicit none

    real(kind=CUSTOM_REAL), dimension(*), intent(out) :: d
    real(kind=CUSTOM_REAL), intent(in) :: deltat,t0,f0
    integer, intent(in) :: NSTEP, Wscale
    integer, intent(in) :: is_laplace, is_window, window_type,mute_near,mute_far
    real(kind=CUSTOM_REAL), intent(in) :: dis_sr, S_x, S_t,V_slow,V_fast,offset_near,offset_far
    integer, intent(out) :: istart,iend

    integer :: itime,NA
    real(kind=CUSTOM_REAL) :: win(NSTEP)

    istart=1
    iend=NSTEP

    !! mute 
    if(mute_near .eq. 1 .and. dis_sr<=offset_near ) then
        d(1:NSTEP) =0.0
        istart=NSTEP
        iend=0
    endif 
    if(mute_far .eq. 1 .and. dis_sr>=offset_far) then
        d(1:NSTEP) =0.0
        istart=NSTEP
        iend=0
    endif

    !! laplace damping spatially and temporally 
    if (is_laplace .eq. 1 .and. istart<iend) then
        ! spatial
        d(1:NSTEP)=d(1:NSTEP)*exp(-dis_sr*S_x)
        ! temporal
        do itime=1,NSTEP
        d(itime)=d(itime)*exp(-((itime-1)*deltat*S_t))
        enddo
    endif

    ! dip-window using slopes 
    if(is_window .eq. 1.and. istart<iend) then
        ! estimate surface-wave arrivals
        istart=max(int((dis_sr/V_fast+t0-1.2/f0)/deltat)+1,istart)
        iend=min(int((dis_sr/V_slow+3.0/f0+t0+1.2/f0)/deltat)+1,iend)
        if(istart<iend) then 
            ! window
            win(1:NSTEP)=0.d0
            call window(NSTEP,istart,iend,window_type,win)
            d(1:NSTEP)=d(1:NSTEP)*win(1:NSTEP)
        endif 
    endif ! window

    ! WT filtering
    if( Wscale .gt. 0.and. istart<iend) then
        call WT(d,NSTEP,Wscale,NA)
        istart=max(istart-NA,1)
        iend=min(iend+NA,NSTEP) 
    endif

end subroutine preprocess

!----------------------------------------------------------------------

subroutine window(npts,istart,iend,window_type,win)
    use constants
    implicit none

    integer, intent(in) :: npts, istart,iend, window_type
    real(kind=CUSTOM_REAL), dimension(*), intent(out) :: win

    integer ::  i, nlen
    real(kind=CUSTOM_REAL) :: sfac1
    real(kind=CUSTOM_REAL), dimension(npts) :: fac

    nlen = iend - istart+1

    ! some constants
    sfac1 = (2./real(nlen))**2   ! for Welch taper

    ! initialization
    win(1:npts) = 0.d0
    fac(1:npts) = 0.d0

    do i = 1, nlen
    if(window_type ==2) then
        fac(i) = 1 - sfac1*((i-1) - real(nlen)/2.)**2
    elseif(window_type ==3) then
        fac(i) = 1. - cos(PI*(i-1)/(nlen-1))**ipwr_t
    elseif(window_type ==4) then
        fac(i) = 0.5 - 0.5*cos(TWOPI*(i-1)/(nlen-1))
    else
        fac(i) = 1. ! boxcar window
    endif
    enddo

    !! return 
    win(istart:iend)=fac(1:nlen)

end subroutine window
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine taper(NSTEP,taper_len,tas)
    use constants
    implicit none

    ! input parameters
    integer, intent(in) :: NSTEP,taper_len
    real, dimension(*), intent(out) :: tas

    integer i, n_width_left,n_width_right
    real ::  omega_left,omega_right, f0, f1

    ! initialization
    tas(1:NSTEP)=1.0

    ! set the number of points corresponding to the taper width
    n_width_left=taper_len
    n_width_right=taper_len

    ! set the taper properties according to type
    omega_left = PI/real(n_width_left)
    omega_right = PI/real(n_width_right)
    f0 = 0.5
    f1 = 0.5

    ! apply the taper symmetrically to left side of the data
    do i = 1, n_width_left
    tas(i) = f0-f1*cos(omega_left*(i-1))
    enddo
    do i = 1, n_width_right
    !  print*,i,omega_right,n_width_right,(f0-f1*cos(omega_right*(i-1)))
    tas(NSTEP+1-i) = f0-f1*cos(omega_right*(i-1))
    end do

end subroutine taper
!------------------------------------------------------------------------
subroutine xcorr(d,s,npts1,npts2,cc_array,ishift,cc_max,sigma)
    ! recommended by Yanhua -- to calculate cross-correlation time series cc_array 
    ! (length of ntps1+npts2-1) for d and s with
    ! lengths of npts1 and npts2 respectively, as well as their maiximal
    ! coefficients cc_max and corresponding shift ishift 
    ! sigma -- A gaussian weighting centered at zero applied to cc_array (if sigma=0
    ! no weighting)
    ! this function can replace xcorr_calc if command the normalization operator
    use constants
    implicit none
    ! inputs
    real(kind=CUSTOM_REAL), dimension(*), intent(in) :: s,d
    integer, intent(in) :: npts1, npts2
    real(kind=CUSTOM_REAL), intent(in) :: sigma

    ! outputs:
    integer, intent(out) :: ishift
    real(kind=CUSTOM_REAL), intent(out) :: cc_max
    real(kind=CUSTOM_REAL), dimension(*), intent(out) :: cc_array

    ! local variables
    integer :: i_left, i_right, i, j,k
    real(kind=CUSTOM_REAL) :: cc

    ! initialise shift and cross correlation to zero
    ishift = 0
    cc_max = 0.0d0

    i_left = -(npts2-1) ! toward right shift (negative)
    i_right =npts1-1    ! toward left shift  (postive)
    ! i is the index to shift to be applied to DATA (d)
    k=0
    do i = i_left, i_right
    k=k+1
    ! cc as a function of i or k
    cc = 0.
    do j = 1,npts2    ! loop over full window length
    if((j+i).ge.1 .and. (j+i).le.npts1) then
        cc = cc + s(j)*d(j+i)  ! d is shifted by i
    endif
    enddo
    cc_array(k)=cc
    ! keeping cc-max only
    if (cc .gt. cc_max) then
        cc_max = cc
        ishift = i
    endif
    !! weighting cc_array with a Gaussian function
    if (sigma>0) then
        cc_array(k)=cc_array(k)*exp(-4*(i**2)/(sigma**2))
    endif

    enddo
end subroutine xcorr
!-------------------------------------------------------------------------
subroutine xconv(d,s,npts1,npts2,c_array,npts3)
    ! Yanhua -- convolution of d with npts1 points and s with npts2 points 
    ! output is c_array with length npts3
    use constants 
    implicit none
    ! inputs
    real(kind=CUSTOM_REAL), dimension(*), intent(in) :: s,d
    integer, intent(in) :: npts1,npts2,npts3

    ! outputs:
    real(kind=CUSTOM_REAL), dimension(*), intent(out) :: c_array

    ! local variables
    integer :: i_left, i_right, i,j,k
    real(kind=CUSTOM_REAL),dimension(npts1+npts2-1) :: cc_array
    real(kind=CUSTOM_REAL) :: cc

    if (npts3 .lt. 1 .or. npts3 .gt. npts1+npts2-1) then
        write(*,*) 'output length is improper: ', npts3
        return
    endif

    i_left = 2
    i_right = npts1+npts2
    ! i is the index to shift to be applied to DATA (d)
    k=0
    do i = i_left, i_right
    k=k+1
    ! cc as a function of i or k
    cc = 0.
    do j = 1,npts1    ! loop over full window length
    if((i-j) .ge. 1 .and. (i-j) .le. npts2) then
        cc = cc + d(j)*s(i-j)  ! d is shifted by i
    endif
    enddo
    cc_array(k)=cc
    enddo
    i_left=1+ceiling((npts1+npts2-1-npts3)*0.5)
    i_right=npts3+ceiling((npts1+npts2-1-npts3)*0.5)
    c_array(1:npts3)=cc_array(i_left:i_right)

end subroutine xconv
!-----------------------------------------------------------------------
subroutine gaussmf(x,sig,c,npts,f)
    ! Gaussian curve membership function
    ! by Yanhua -- see matlab y = gaussmf(x,[sig c])
    use constants
    implicit none
    ! inputs
    real(kind=CUSTOM_REAL), dimension(*), intent(in) :: x
    real(kind=CUSTOM_REAL), intent(in) :: sig,c
    integer, intent(in) :: npts

    ! outputs:
    real(kind=CUSTOM_REAL), dimension(*), intent(out) :: f

    integer :: i

    do i=1,npts
    f(i)=exp(-(x(i)-c)**2/(2*sig**2))
    enddo
end subroutine gaussmf
!---------------------------------------------------------------------------
subroutine gauspuls(x,npts,fc,sig,c,fe,f )
    ! Gaussian curve membership function
    use constants
    implicit none
    ! inputs
    real(kind=CUSTOM_REAL), dimension(*), intent(in) :: x
    real(kind=CUSTOM_REAL), intent(in) :: fc,sig,c
    integer, intent(in) :: npts

    ! outputs:
    real(kind=CUSTOM_REAL), dimension(*), intent(out) :: fe, f

    ! compute gaussian envelope
    call gaussmf(x,sig,c,npts,fe)
    ! gaussian modulated sine/cosine function
    f(1:npts) = fe(1:npts) * cos(TWOPI*fc*x(1:npts));    

end subroutine gauspuls
! ---------------------------------------------------------------------------
subroutine WT(seism,NSTEP,level,NA)
    use constants
    implicit none

    real(kind=CUSTOM_REAL), dimension(*), intent(out) :: seism
    integer, intent(in) :: NSTEP,level
    integer :: iend,istart,i,j,st
    integer, intent(out) :: NA
    real(kind=CUSTOM_REAL),dimension(:,:),allocatable :: basis, seism_WT
    real(kind=CUSTOM_REAL) :: wc 
    real(kind=CUSTOM_REAL),dimension(:), allocatable :: nonzero_basis
    character(len=500) :: fname,filename

    NA=0
    if(level .gt. 0) then
        if(level==12) NA=45046
        if(level==11) NA=22518
        if(level==10) NA=11254   
        if(level==9) NA=5622
        if(level==8) NA=2806
        if(level==7) NA=1398
        if(level==6) NA=694  
        if(level==5) NA=342 
        if(level==4) NA=166
        if(level==3) NA=78
        if(level==2) NA=34
        if(level==1) NA=12
    endif !! level

    if(level .gt. 0 .and. NA .gt. 0) then
        allocate(nonzero_basis(NA))
        allocate(basis(NSTEP,1))
        allocate(seism_WT(NSTEP,1))
        seism_WT(1:NSTEP,1)=0.0

        ! print*,' load basis functions for wavelet Daubechies',nvm,' scale',level
        write(fname,"(a,i0,a,i0)") 'scale_basis_Daubechies',nvm,'_scale',level
        ! filename=''//trim(WT_directory)//'/'//trim(fname)//'.dat'
        filename='WT_basis/'//trim(fname)//'.dat'
        ! print*,filename
        OPEN (UNIT=20,FILE=filename,STATUS='OLD',action='read',iostat=st)
        if(st>0) then
            print*,'Error opening file: ',filename
            stop
        else
            read(20,*) nonzero_basis
        end if
        close(20)

        ! initialization
        iend=1
        istart=2-NA
        ! find shifted basis 
        do while (istart<=NSTEP)
        basis(1:NSTEP,1)=0.0
        j=0
        do i=istart,iend
        j=j+1
        if(i>=1 .and. i<=NSTEP) then
            basis(i,1)=nonzero_basis(j)
        endif
        enddo
        !! WT
        wc=dot_product(seism(1:NSTEP),basis(1:NSTEP,1))
        !   print*, 'istart=',istart,' wc=',wc
        ! inverse wavelet transform to construct data
        seism_WT(1:NSTEP,1)=seism_WT(1:NSTEP,1)+wc*basis(1:NSTEP,1)
        !    print*, ' prepare for next basis '
        iend=iend+2**level
        istart=istart+2**level
        !    print*, 'istart=',istart,' iend=',iend
        enddo
        seism(1:NSTEP)=seism_WT(1:NSTEP,1)

        deallocate(nonzero_basis)
        deallocate(basis)
        deallocate(seism_WT)
    endif  !! level

end subroutine WT

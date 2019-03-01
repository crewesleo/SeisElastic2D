!! main subroutines for mt measurement 
!

! =================================================================================================
! subroutine mt_measure()
! Boxcar/Cosine/Multitaper estimates of the phase function between dat1 and dat2 (dat1-dat2)
!
!  Input:
!        dat1_dt(:), dat2_dt(:), t0, dt, npts -- original dat1 and dat2 array
!        tstart, tend -- start and end of the measurement window (can be from Alessia's code)
!  Output:
!        istart -- starting index of the windowed portion of  original trace
!        dat1_dtw(:), dat2_dtw(:), nlen -- windowed and shifted dat1, windowed dat2
!        tshift, dlnA, cc_max -- time shift and amplitude cross-correlation measurements
!        i_right -- the maximum reliable frequency estimate index
!        dtau_w(:), dlnA_w(:) -- estimates of travel-time and amplitude anomaly
!        err_dt(:), err_dlnA(:) -- error bar of the travel-time and amplitude estimates (MT only)
!
! =================================================================================================

subroutine mt_measure(dat1,dat2,npts,deltat,nlen,tshift_cc,dlnA_cc, i_fstart,i_fend,&
        wvec,&
        trans_func,dtau_w,dlnA_w,err_dtau,err_dlnA)
    use constants
    implicit none

    integer, intent(in) :: nlen, npts
    real(kind=CUSTOM_REAL), dimension(*), intent(in) :: dat1,dat2
    real(kind=CUSTOM_REAL), intent(in) ::  deltat,tshift_cc,dlnA_cc 
    integer, intent(in) :: i_fstart,i_fend
    real(kind=CUSTOM_REAL), dimension(*), intent(in) :: wvec
    complex (CUSTOM_REAL), dimension(NPT), intent(out) :: trans_func
    real(kind=CUSTOM_REAL), dimension(*), intent(out) :: dtau_w,dlnA_w,err_dtau,err_dlnA


    !! fft -- double size
    ! FFT parameters 
    real(kind=SIZE_DOUBLE) :: ampmax,wtr_mtm_use,wtr_use
    integer :: i,ictaper,iom

    !! multitaper
    real(kind=SIZE_DOUBLE), dimension(NPT) :: ey1,ey2
    real(kind=SIZE_DOUBLE), dimension(:,:),allocatable :: tas

    !! fft to estimate transfer function
    ! tapered data in time domain
    real(kind=CUSTOM_REAL), dimension(NPT) :: dat1_dtw_h,dat2_dtw_h
    ! fft domain
    complex (SIZE_DOUBLE), dimension(NPT) :: dat2_dtwo, dat1_dtwo, &
        dat2_dtw_ho, dat1_dtw_ho,  &
        top_mtm,bot_mtm
    !! transfer function to phase/amplitude measurement
    real(kind=CUSTOM_REAL), dimension(NPT) :: phi_mtm,dtau_mtm,abs_mtm,dlnA_mtm
    real(kind=CUSTOM_REAL), dimension(:,:),allocatable :: phi_mul,dtau_mul,abs_mul,dlnA_mul

    !! jackknife error
    real(kind=CUSTOM_REAL) :: edt_ave,edt_iom,edA_ave,edA_iom

    character (len=MAX_FILENAME_LEN) :: filename

    !-------------------------------------------------------------

    if ( nlen > npts .or. nlen <1) then
        print*, nlen, npts
        stop 'Check nlen'
    endif

    !-------------------------------------------------------------------------------
    ! multitaper estimation of transfer function 
    !-------------------------------------------------------------------------------
    !print*,'nlen,NW,ntaper -- ',nlen, NW, ntaper
    ! calculate the tapers
    allocate(tas(NPT,ntaper))
    call staper(nlen, dble(NW), ntaper, tas, NPT, ey1, ey2)

    if( DISPLAY_DETAILS) then
        do ictaper=1,ntaper
        write(filename,'(i)') ictaper
        open(1,file=trim(output_dir)//'/'//'taper_t_'//trim(adjustl(filename)),status='unknown')
        open(2,file=trim(output_dir)//'/'//'taper_w_'//trim(adjustl(filename)),status='unknown')
        do i = 1,nlen
        write(1,'(f15.2,E15.5)') (i-1)*deltat,tas(i,ictaper)
        ! write(1,'(f15.2,<ntaper>e15.5)') (i-1)*deltat,tas(i,1:ntaper)
        enddo
        close(1)
        dat1_dtw_ho(:) = cmplx(0.0_SIZE_DOUBLE,0.0_SIZE_DOUBLE)
        dat1_dtw_h(1:nlen) = real(tas(1:nlen,ictaper),CUSTOM_REAL)
        dat1_dtw_ho(1:nlen) = dcmplx(dat1_dtw_h(1:nlen),0.0_SIZE_DOUBLE)
        call fft(LNPT,dat1_dtw_ho,dble(FORWARD_FFT),dble(deltat))
        do  i =  1,50
        write(2,'(2f15.2,E15.5)') wvec(i),abs(dat1_dtw_ho(i))
        enddo
        close(1)
        close(2)
        enddo
    endif

    !! estimate transfer function from mtm 
    ! initialize transfer function terms
    top_mtm(:)   = cmplx(0.0_SIZE_DOUBLE,0.0_SIZE_DOUBLE)
    bot_mtm(:)   = cmplx(0.0_SIZE_DOUBLE,0.0_SIZE_DOUBLE)
    trans_func(:) = cmplx(0.0_SIZE_DOUBLE,0.0_SIZE_DOUBLE)

    do ictaper = 1, ntaper

    dat2_dtw_ho(:) = cmplx(0.0_SIZE_DOUBLE,0.0_SIZE_DOUBLE) ! note: this has to be initialized inside the loop
    dat1_dtw_ho(:) = cmplx(0.0_SIZE_DOUBLE,0.0_SIZE_DOUBLE)

    ! apply time-domain taper
    do i = 1, nlen
    dat2_dtw_h(i) = dat2(i) * real(tas(i,ictaper),CUSTOM_REAL)     ! single-tapered, windowed data2
    dat1_dtw_h(i) = dat1(i) * real(tas(i,ictaper),CUSTOM_REAL)  ! single-tapered, windowed, shifted data1
    enddo

    dat2_dtw_ho(1:nlen) = dcmplx(dat2_dtw_h(1:nlen),0.0_SIZE_DOUBLE)
    dat1_dtw_ho(1:nlen) = dcmplx(dat1_dtw_h(1:nlen),0.0_SIZE_DOUBLE)

    ! apply FFT to get complex spectra
    call fft(LNPT,dat2_dtw_ho,dble(FORWARD_FFT),dble(deltat))
    call fft(LNPT,dat1_dtw_ho,dble(FORWARD_FFT),dble(deltat))

    if(DISPLAY_DETAILS .and. ictaper==1) then
        open(1,file=trim(output_dir)//'/fft_in',status='unknown')
        open(2,file=trim(output_dir)//'/fft_out',status='unknown')
        do i = 1, nlen
        write(1,'(3E15.5)') i,dat1_dtw_h(i),dat2_dtw_h(i)
        enddo
        do  i =  i_fstart,i_fend
        write(2,'(3E15.5)') wvec(i)/TWOPI, abs(dat1_dtw_ho(i)),abs(dat2_dtw_ho(i))
        enddo
        close(1)
        close(2)
    endif

    ! calculate top and bottom of transfer function 
    do i = i_fstart,i_fend
    top_mtm(i) = top_mtm(i) + dat1_dtw_ho(i) * conjg(dat2_dtw_ho(i))   
    bot_mtm(i) = bot_mtm(i) + dat2_dtw_ho(i) * conjg(dat2_dtw_ho(i))  
    enddo

    enddo  ! ictapers

    ! water level for transfer function
    ampmax = maxval(abs(bot_mtm))

    ! LQY -- is this too small ???
    wtr_mtm_use = ampmax * wtr_mtm**2
    wtr_mtm_use = ampmax * 0.01

    ! calculate MT transfer function using water level
    do i =  i_fstart,i_fend
    !      if(abs(bot_mtm(i)) > abs(wtr_mtm_use)) trans_func(i) = cmplx(top_mtm(i) / bot_mtm(i))
    !      if(abs(bot_mtm(i)) < abs(wtr_mtm_use)) trans_func(i) = cmplx(top_mtm(i) / (bot_mtm(i)+wtr_mtm_use))
    trans_func(i) = cmplx(top_mtm(i) / (bot_mtm(i)+wtr_mtm_use))
    enddo

    if(DISPLAY_DETAILS) then
        ! print*
        ! print*,'water level for multi taper transfer function is : ',wtr_use
        ! fft of untapered inputs 
        dat1_dtwo = cmplx(0.0_SIZE_DOUBLE,0.0_SIZE_DOUBLE)
        dat2_dtwo = cmplx(0.0_SIZE_DOUBLE,0.0_SIZE_DOUBLE)
        dat1_dtwo(1:nlen) = dcmplx(dat1(1:nlen),0.0)
        dat2_dtwo(1:nlen) = dcmplx(dat2(1:nlen),0.0)
        call fft(LNPT,dat1_dtwo,dble(FORWARD_FFT),dble(deltat))
        call fft(LNPT,dat2_dtwo,dble(FORWARD_FFT),dble(deltat))
        open(1,file=trim(output_dir)//'/top_bot_mtm',status='unknown')
        open(2,file=trim(output_dir)//'/trans_func',status='unknown')
        open(3,file=trim(output_dir)//'/trans_error',status='unknown')
        do  i =  i_fstart,i_fend
        write(1,'(3E15.5)') wvec(i)/TWOPI, abs(top_mtm(i)),abs(bot_mtm(i))
        write(2,'(2E15.5)') wvec(i)/TWOPI, abs(trans_func(i))
        write(3,'(4E15.5)') wvec(i)/TWOPI, abs(dat1_dtwo(i)),abs(dat2_dtwo(i)),abs(dcmplx(trans_func(i))*dat2_dtwo(i))
        enddo
        close(1)
        close(2)
        close(3)
    endif

    !! phase and amplitude measurement derived from transfer functions
    phi_mtm(:) = 0.d0
    abs_mtm(:) = 0.d0
    dtau_mtm(:) = 0.d0
    dlnA_mtm(:) = 0.d0
    call write_trans(trans_func, wvec, i_fstart,i_fend,tshift_cc,dlnA_cc,phi_mtm, abs_mtm,dtau_mtm, dlnA_mtm)

    ! pass results to main routine
    dtau_w(i_fstart:i_fend) = dtau_mtm(i_fstart:i_fend)
    dlnA_w(i_fstart:i_fend) = dlnA_mtm(i_fstart:i_fend)

    if(DISPLAY_DETAILS) then
        open(1,file=trim(output_dir)//'/dtau_dlnA',status='unknown')
        do  i =  i_fstart,i_fend
        write(1,'(3E15.5)') wvec(i)/TWOPI,dtau_w(i),dlnA_w(i)
        enddo
        close(1)
    endif

    !-------------------------------------------------------------------------------
    ! multitaper error estimation
    !-------------------------------------------------------------------------------

    if (ntaper > 1 .and. USE_ERROR_MT) then

        ! allocate Jacknife MT estimates
        allocate(phi_mul(NPT,ntaper))
        allocate(dtau_mul(NPT,ntaper))
        allocate(abs_mul(NPT,ntaper))
        allocate(dlnA_mul(NPT,ntaper))

        do iom = 1, ntaper

        top_mtm(:) = cmplx(0.,0.)
        bot_mtm(:) = cmplx(0.,0.)

        do ictaper = 1, ntaper
        if(ictaper.eq.iom) cycle

        ! apply ictaper-th taper
        dat2_dtw_h(1:nlen) = dat2(1:nlen) * real(tas(1:nlen,ictaper),CUSTOM_REAL)
        dat1_dtw_h(1:nlen) = dat1(1:nlen) * real(tas(1:nlen,ictaper),CUSTOM_REAL)

        ! complex tapered series
        dat2_dtw_ho(:) = cmplx(0.,0.)
        dat1_dtw_ho(:) = cmplx(0.,0.)
        dat2_dtw_ho(1:nlen) = dcmplx(dat2_dtw_h(1:nlen),0.)
        dat1_dtw_ho(1:nlen) = dcmplx(dat1_dtw_h(1:nlen),0.)

        ! apply f.t. to get complex spectra
        call fft(LNPT,dat2_dtw_ho,dble(FORWARD_FFT),dble(deltat))
        call fft(LNPT,dat1_dtw_ho,dble(FORWARD_FFT),dble(deltat))

        ! calculate top and bottom of Jacknife transfer function
        do i = i_fstart,i_fend
        top_mtm(i) = top_mtm(i) + dat1_dtw_ho(i) * conjg(dat2_dtw_ho(i))
        bot_mtm(i) = bot_mtm(i) + dat2_dtw_ho(i) * conjg(dat2_dtw_ho(i))
        enddo
        enddo ! ictaper

        ! water level
        ampmax = maxval(abs(bot_mtm))
        wtr_use = cmplx(ampmax * wtr_mtm ** 2, 0.)
        wtr_mtm_use = ampmax * 0.01

        !  calculate transfer function using water level
        do i = i_fstart,i_fend
        !if(abs(bot_mtm(i)).gt.abs(wtr_use)) trans_func(i) = top_mtm(i) / bot_mtm(i)
        !if(abs(bot_mtm(i)).le.abs(wtr_use)) trans_func(i) = top_mtm(i) /(bot_mtm(i)+wtr_use)
        trans_func(i) = cmplx(top_mtm(i) /(bot_mtm(i)+wtr_use))
        enddo

        call write_trans(trans_func, wvec,i_fstart,i_fend,tshift_cc,dlnA_cc,&
            phi_mul(:,iom),abs_mul(:,iom),dtau_mul(:,iom), dlnA_mul(:,iom))

        enddo ! iom

        !----------------------
        err_dtau(1:NPT)   = 0.
        err_dlnA(1:NPT)   = 0.

        do i = i_fstart,i_fend

        edt_ave   = 0.
        edA_ave   = 0. 

        do iom = 1, ntaper
        edt_iom = ntaper*dtau_mtm(i) - (ntaper-1)*dtau_mul(i,iom)
        edt_ave = edt_ave + edt_iom

        edA_iom = ntaper*dlnA_mtm(i) - (ntaper-1)*dlnA_mul(i,iom)
        edA_ave = edA_ave + edA_iom

        enddo

        edt_ave   = edt_ave   / (ntaper)
        edA_ave   = edA_ave   / (ntaper)

        do iom = 1, ntaper
        err_dlnA(i)  = err_dlnA(i) + ( dlnA_mul(i,iom) - edA_ave)**2
        err_dtau(i)   = err_dtau(i)  + (dtau_mul(i,iom) - edt_ave)**2
        enddo

        err_dlnA(i)  =  sqrt( err_dlnA(i) / (ntaper * (ntaper-1) ) )
        err_dtau(i)   =  sqrt( err_dtau(i) / (ntaper * (ntaper-1) ) )
        ! set the error bar for the first point corresponding to
        ! static offset to be large, which makes no contribution to
        ! the adjoint source
        !  if (i == 1) err_dt(i) = LARGE_VAL

        enddo ! i_fstart

        deallocate(phi_mul)
        deallocate(dtau_mul)
        deallocate(abs_mul)
        deallocate(dlnA_mul)

        if(DISPLAY_DETAILS) then
            open(1,file=trim(output_dir)//'/err_dtau_dlnA',status='unknown')
            do  i =  i_fstart,i_fend
            write(1,'(3E15.5)') wvec(i)/TWOPI,err_dtau(i),err_dlnA(i)
            enddo
            close(1)
        endif

    endif  ! if use_error_MT

    !     ------------------------------------------------------------------
    !     End error calculation loop
    !     ------------------------------------------------------------------

    deallocate(tas)

end subroutine mt_measure
! ------------------------------------------------------------------------------------
subroutine mtm_adj(syn,npts,deltat,nlen,df,i_fstart,i_fend,dtau_w,dlnA_w,&
        err_dt_cc,err_dlnA_cc, &
        err_dtau_mt,err_dlnA_mt, &
        fp,fq)
    use constants
    implicit none

    integer, intent(in) :: npts,nlen
    real(kind=CUSTOM_REAL), dimension(*), intent(in) :: syn
    real(kind=CUSTOM_REAL), intent(in) ::  deltat,df,err_dt_cc,err_dlnA_cc
    integer, intent(in) :: i_fstart,i_fend
    real(kind=CUSTOM_REAL), dimension(*), intent(in) :: dtau_w,dlnA_w,err_dtau_mt,err_dlnA_mt
    real(kind=CUSTOM_REAL), dimension(*), intent(out) :: fp, fq 

    ! frequency-domain taper 
    integer :: nfrange, window_type =3
    real(kind=CUSTOM_REAL), dimension(NPT) :: w_taper, wp_taper, wq_taper
    real(kind=CUSTOM_REAL) :: ffac, dtau_wtr, dlnA_wtr, err_t, err_A

    !! fft -- double size
    ! FFT parameters 
    real(kind=SIZE_DOUBLE) :: ampmax,wtr_mtm_use,wtr_use
    integer :: i,ictaper,iom

    !! multitaper
    real(kind=SIZE_DOUBLE), dimension(NPT) :: ey1,ey2
    real(kind=SIZE_DOUBLE), dimension(:,:),allocatable :: tas

    !! fft to evaluate adj
    ! tapered data in time domain
    real(kind=CUSTOM_REAL), dimension(NPT) :: syn_dtw_h, syn_vtw_h
    ! fft and ifft (forward and backward)
    complex (SIZE_DOUBLE), dimension(:,:),allocatable :: syn_dtw_ho_all,syn_vtw_ho_all
    complex (SIZE_DOUBLE), dimension(NPT) :: p_bot_mtm, q_bot_mtm, &
        pwc_adj,qwc_adj
    real(kind=SIZE_DOUBLE), dimension(NPT) ::  dtau_pj_t,dlnA_qj_t

    ! frequency-domain tapers
    ! THIS CHOICE WILL HAVE AN EFFECT ON THE ADJOINT SOURCES
    nfrange = i_fend -i_fstart 
    w_taper(:) = 0.
    do i = i_fstart, i_fend    
    ! if(window_type ==1) w_taper(i) = 1. ! boxcar window
    ! if(window_type ==2) w_taper(i) = 1. - (2.0/nw)**2 * ((i-1) - nw/2.0)**2       !welch window
    if(window_type ==3) w_taper(i) = &
        1. - cos(PI*(i-i_fstart)/(nfrange))**ipwr_w !cosine window
    enddo

    ! compute normalization factor for w_taper
    ! note: 2 is needed for the integration from -inf to inf
    ffac = 2.0 * df * sum(w_taper(i_fstart:i_fend) )   
    if (DISPLAY_DETAILS) then
        print* 
        print*, ' normalize frequency-domain taper' 
        print*, ' df = ', df 
        print *, 'Taper normalization factor, ffac = ', ffac
    endif
    ! no error estimate
    ! only adds normalization factor
    wp_taper = w_taper / ffac
    wq_taper = w_taper / ffac

    ! add error estimated 
    if (USE_ERROR_CC) then
        !          ! CC error estimate
        do i = i_fstart,i_fend
        wp_taper(i) = wp_taper(i) / (err_dt_cc ** 2)
        wq_taper(i) = wq_taper(i) / (err_dlnA_cc ** 2)
        enddo 
    endif
    if (USE_ERROR_MT) then
        ! MT jack-knife error estimate
        dtau_wtr = WTR * sum(abs(dtau_w(i_fstart:i_fend)))/ nfrange
        dlnA_wtr = WTR * sum(abs(dlnA_w(i_fstart:i_fend)))/ nfrange
        print*
        print*,'to add mt error :'
        print*,'water level for err_dtau is ',dtau_wtr
        print*,'water level for err_dlnA is ',dlnA_wtr

        do i = i_fstart,i_fend
        err_t = err_dtau_mt(i)
        err_A = err_dlnA_mt(i)
        if (err_t < dtau_wtr)  err_t = err_t + dtau_wtr
        if (err_A < dlnA_wtr)  err_A = err_A + dlnA_wtr
        wp_taper(i) = wp_taper(i) / (err_t ** 2)
        wq_taper(i) = wq_taper(i) / (err_A ** 2)
        enddo
    endif

    if( DISPLAY_DETAILS) then
        open(1,file=trim(output_dir)//'/frequency_taper',status='unknown')
        do  i =  i_fstart,i_fend
        write(1,'(I5,3e15.5)') i,w_taper(i),wp_taper(i),wq_taper(i)
        enddo
        close(1)
    endif

    ! allocate MT variables
    allocate(syn_dtw_ho_all(NPT,ntaper))
    allocate(syn_vtw_ho_all(NPT,ntaper))

    ! calculate the tapers
    allocate(tas(NPT,ntaper))
    call staper(nlen, NW, ntaper, tas, NPT, ey1, ey2)

    p_bot_mtm = 0.
    q_bot_mtm = 0.

    do ictaper = 1,ntaper

    ! tapered synthetic displacement
    syn_dtw_h(1:nlen) = syn(1:nlen) * real(tas(1:nlen,ictaper),CUSTOM_REAL)

    ! compute velocity of tapered syn
    call compute_vel(syn_dtw_h,npts,deltat,nlen,syn_vtw_h)

    ! single-tapered complex synthetic displacement and velocity
    syn_dtw_ho_all(:,ictaper) = 0.
    syn_vtw_ho_all(:,ictaper) = 0.
    syn_dtw_ho_all(1:nlen,ictaper) = dcmplx(syn_dtw_h(1:nlen),0.)
    syn_vtw_ho_all(1:nlen,ictaper) = dcmplx(syn_vtw_h(1:nlen),0.)

    ! apply FFT get complex spectra
    call fft(LNPT,syn_dtw_ho_all(:,ictaper),dble(FORWARD_FFT),dble(deltat))
    call fft(LNPT,syn_vtw_ho_all(:,ictaper),dble(FORWARD_FFT),dble(deltat))

    p_bot_mtm(:) = p_bot_mtm(:) + syn_vtw_ho_all(:,ictaper) &
        * conjg(syn_vtw_ho_all(:,ictaper))
    q_bot_mtm(:) = q_bot_mtm(:) + syn_dtw_ho_all(:,ictaper) &
        * conjg(syn_dtw_ho_all(:,ictaper))

    enddo ! ictaper

    if( DISPLAY_DETAILS) then
        open(2,file=trim(output_dir)//'/adj_bot_pq',status='unknown')
        do  i = i_fstart,i_fend 
        write(2,*) i,abs(p_bot_mtm(i)), abs(q_bot_mtm(i))
        enddo
        close(2)
    endif

    ampmax = maxval(abs(p_bot_mtm))
    !    wtr_use = ampmax * wtr_mtm**2
    wtr_use = ampmax * 0.01

    ! compute p_j, q_j, P_j, Q_j and adjoint source fp, fq
    fp(1:npts) = 0.
    fq(1:npts) = 0.
    do ictaper = 1,ntaper

    ! compute p_j(w) and q_j(w)
    pwc_adj(:) = cmplx(0.,0.)
    qwc_adj(:) = cmplx(0.,0.)

    do i = i_fstart,i_fend
    ! if(abs(p_bot_mtm(i)) > abs(wtr_use)) pwc_adj(i) = &
    !        syn_vtw_ho_all(i,ictaper) / p_bot_mtm(i)
    ! if(abs(p_bot_mtm(i)) > abs(wtr_use)) pwc_adj(i) = &
    !        syn_vtw_ho_all(i,ictaper) / (p_bot_mtm(i) + wtr_use)
    pwc_adj(i) =  syn_vtw_ho_all(i,ictaper) / p_bot_mtm(i)
    qwc_adj(i) = -syn_dtw_ho_all(i,ictaper) / q_bot_mtm(i)
    enddo

    ! compute P_j(w) and Q_j(w)
    ! NOTE: the MT measurement is incorporated here
    !             also note that wp_taper and wq_taper can contain
    !             uncertainty estimations
    ! adds misfit measurement dtau, dlnA
    pwc_adj(i_fstart: i_fend) = pwc_adj(i_fstart: i_fend) &
        * dcmplx(dtau_w(i_fstart: i_fend),0.) &
        * dcmplx(wp_taper(i_fstart: i_fend),0.)
    qwc_adj(i_fstart: i_fend) = qwc_adj(i_fstart: i_fend) &
        * dcmplx(dlnA_w(i_fstart: i_fend),0.) &
        * dcmplx(wq_taper(i_fstart: i_fend),0.)

    ! IFFT into the time domain
    call fftinv(LNPT,pwc_adj,dble(REVERSE_FFT),dble(deltat),dtau_pj_t)
    call fftinv(LNPT,qwc_adj,dble(REVERSE_FFT),dble(deltat),dlnA_qj_t)

    ! create adjoint source
    ! applies taper to time signal
    fp(1:nlen) = fp(1:nlen) + real(tas(1:nlen,ictaper) * dtau_pj_t(1:nlen),CUSTOM_REAL)
    fq(1:nlen) = fq(1:nlen) + real(tas(1:nlen,ictaper) * dlnA_qj_t(1:nlen),CUSTOM_REAL)

    enddo  ! ictaper
    deallocate(tas)
    deallocate(syn_dtw_ho_all)
    deallocate(syn_vtw_ho_all)

end subroutine mtm_adj
! ------------------------------------------------------------------------------------
subroutine mtm_DD_adj(s1,s2,npts,deltat,nlen,df,i_fstart,i_fend,ddtau_w,ddlnA_w,&
        err_dt_cc_obs,err_dt_cc_syn,err_dlnA_cc_obs,err_dlnA_cc_syn, &
        err_dtau_mt_obs,err_dtau_mt_syn,err_dlnA_mt_obs,err_dlnA_mt_syn, &
        fp1,fp2,fq1,fq2)
    use constants
    implicit none

    integer, intent(in) :: npts,nlen
    real(kind=CUSTOM_REAL), dimension(*), intent(in) :: s1,s2
    real(kind=CUSTOM_REAL), intent(in) :: deltat,df
    real(kind=CUSTOM_REAL), intent(in) :: err_dt_cc_obs,err_dt_cc_syn,err_dlnA_cc_obs,err_dlnA_cc_syn
    integer, intent(in) :: i_fstart,i_fend
    real(kind=CUSTOM_REAL), dimension(*), intent(in) :: ddtau_w,ddlnA_w
    real(kind=CUSTOM_REAL), dimension(*), intent(in) :: err_dtau_mt_obs,err_dtau_mt_syn,err_dlnA_mt_obs,err_dlnA_mt_syn
    real(kind=CUSTOM_REAL), dimension(*), intent(out) :: fp1,fp2,fq1,fq2

    ! frequency-domain taper 
    integer :: nfrange, window_type =3
    real(kind=CUSTOM_REAL), dimension(NPT) :: w_taper, wp_taper, wq_taper
    real(kind=CUSTOM_REAL) :: ffac, ddtau_wtr, ddlnA_wtr, err_t, err_A

    !! fft -- double size
    ! FFT parameters 
    real(kind=SIZE_DOUBLE) :: ampmax,wtr_mtm_use,wtr_use
    integer :: i,ictaper,iom

    !! multitaper
    real(kind=SIZE_DOUBLE), dimension(NPT) :: ey1,ey2
    real(kind=SIZE_DOUBLE), dimension(:,:),allocatable :: tas

    !! fft to evaluate adj
    ! tapered data in time domain
    real(kind=CUSTOM_REAL), dimension(NPT) :: s1_dtw_h, s1_vtw_h
    real(kind=CUSTOM_REAL), dimension(NPT) :: s2_dtw_h, s2_vtw_h
    ! fft and ifft (forward and backward)
    complex (SIZE_DOUBLE), dimension(:,:),allocatable :: s1_dtw_ho_all,s1_vtw_ho_all
    complex (SIZE_DOUBLE), dimension(:,:),allocatable :: s2_dtw_ho_all,s2_vtw_ho_all
    complex (SIZE_DOUBLE), dimension(NPT) :: Mtr1_mtm,Mtr2_mtm,Mtr3_mtm,Mtr4_mtm,Mtr5_mtm
    complex (SIZE_DOUBLE), dimension(NPT) :: p1_adj,pw1_adj,p2_adj,pw2_adj 
    complex (SIZE_DOUBLE), dimension(NPT) :: q1_adj,qw1_adj,q2_adj,qw2_adj
    complex (SIZE_DOUBLE), dimension(NPT) :: fp1_adj,fp2_adj,fq1_adj,fq2_adj
    real(kind=SIZE_DOUBLE), dimension(NPT) :: fp1_adj_t,fp2_adj_t,fq1_adj_t,fq2_adj_t

    ! frequency-domain taper 
    !    integer :: nfrange, window_type =3
    !    real(kind=CUSTOM_REAL), dimension(NPT) :: w_taper, wp_taper, wq_taper, ey1, ey2
    !    real(kind=CUSTOM_REAL) :: ffac, ddtau_wtr, ddlnA_wtr, err_t, err_A

    ! mtm adj
    !    real(kind=CUSTOM_REAL), dimension(NPT) :: s1_dtw_h, s1_vtw_h,s2_dtw_h, s2_vtw_h
    !    complex (CUSTOM_REAL), dimension(:,:),allocatable ::  s1_dtw_ho_all,s1_vtw_ho_all
    !    complex (CUSTOM_REAL), dimension(:,:),allocatable ::  s2_dtw_ho_all,s2_vtw_ho_all
    !    comple (CUSTOM_REAL), dimension(NPT) :: Mtr1_mtm,Mtr2_mtm,Mtr3_mtm,Mtr4_mtm,Mtr5_mtm
    !    complex (CUSTOM_REAL), dimension(NPT) :: p1_adj,pw1_adj,p2_adj,pw2_adj
    !    complex (CUSTOM_REAL), dimension(NPT) :: q1_adj,qw1_adj,q2_adj,qw2_adj
    !    complex (CUSTOM_REAL), dimension(NPT) :: fp1_adj,fp2_adj,fq1_adj,fq2_adj
    !    real(kind=CUSTOM_REAL), dimension(NPT) :: fp1_adj_t,fp2_adj_t,fq1_adj_t,fq2_adj_t
    !    real(kind=CUSTOM_REAL), dimension(:,:),allocatable :: tas

    ! frequency-domain tapers
    ! THIS CHOICE WILL HAVE AN EFFECT ON THE ADJOINT SOURCES
    nfrange = i_fend -i_fstart
    w_taper(:) = 0.
    do i = i_fstart, i_fend
    ! boxcar 
    ! if(window_type ==1) w_taper(i) = 1. ! boxcar window
    ! welch 
    ! if(window_type ==2) w_taper(i) = 1. - (2.0/nw)**2 * ((i-1) - nw/2.0)**2
    ! cosine
    if(window_type ==3) w_taper(i) = &
        1. - cos(PI*(i-i_fstart)/(nfrange))**ipwr_w
    enddo

    ! compute normalization factor for w_taper
    ! note: 2 is needed for the integration from -inf to inf
    ffac = 2.0 * df * sum(w_taper(i_fstart:i_fend) )
    if (DISPLAY_DETAILS) then
        print*
        print*, ' normalize frequency-domain taper'
        print*, ' df = ', df
        print *, 'Taper normalization factor, ffac = ', ffac
    endif
    ! no error estimate
    ! only adds normalization factor
    wp_taper = w_taper / ffac
    wq_taper = w_taper / ffac

    ! add error estimated 
    if (USE_ERROR_CC) then
        !          ! CC error estimate
        do i = i_fstart,i_fend
        wp_taper(i) = wp_taper(i) / (err_dt_cc_obs * err_dt_cc_syn)
        wq_taper(i) = wq_taper(i) / (err_dlnA_cc_obs * err_dlnA_cc_syn)
        enddo
    endif
    if (USE_ERROR_MT) then
        ! MT jack-knife error estimate
        ddtau_wtr = (WTR * sum(abs(ddtau_w(i_fstart:i_fend)))/ nfrange)**2
        ddlnA_wtr = (WTR * sum(abs(ddlnA_w(i_fstart:i_fend)))/ nfrange)**2
        print*
        print*,'to add mt error :'
        print*,'water level for err_dtau is ',ddtau_wtr
        print*,'water level for err_dlnA is ',ddlnA_wtr

        do i = i_fstart,i_fend
        err_t = err_dtau_mt_obs(i) * err_dtau_mt_syn(i)
        err_A = err_dlnA_mt_obs(i) * err_dlnA_mt_syn(i)
        if (err_t < ddtau_wtr)  err_t = err_t + ddtau_wtr
        if (err_A < ddlnA_wtr)  err_A = err_A + ddlnA_wtr
        wp_taper(i) = wp_taper(i) / err_t 
        wq_taper(i) = wq_taper(i) / err_A 
        enddo
    endif

    if( DISPLAY_DETAILS) then
        open(1,file=trim(output_dir)//'/frequency_taper',status='unknown')
        do  i =  i_fstart,i_fend
        write(1,'(I5,3e15.5)') i,w_taper(i),wp_taper(i),wq_taper(i)
        enddo
        close(1)
    endif

    ! allocate MT variables
    allocate(s1_dtw_ho_all(NPT,ntaper))
    allocate(s1_vtw_ho_all(NPT,ntaper))
    allocate(s2_dtw_ho_all(NPT,ntaper))
    allocate(s2_vtw_ho_all(NPT,ntaper))


    ! calculate the tapers
    allocate(tas(NPT,ntaper))
    call staper(nlen, NW, ntaper, tas, NPT, ey1, ey2)

    !! constant terms for adj 
    Mtr1_mtm = 0.
    Mtr2_mtm = 0.
    Mtr3_mtm = 0.
    Mtr4_mtm = 0.
    Mtr5_mtm = 0.

    do ictaper = 1,ntaper

    ! tapered synthetic displacement
    s1_dtw_h(1:nlen) = s1(1:nlen) * real(tas(1:nlen,ictaper),CUSTOM_REAL)
    s2_dtw_h(1:nlen) = s2(1:nlen) * real(tas(1:nlen,ictaper),CUSTOM_REAL)

    ! compute velocity of tapered syn
    call compute_vel(s1_dtw_h,npts,deltat,nlen,s1_vtw_h)
    call compute_vel(s2_dtw_h,npts,deltat,nlen,s2_vtw_h)


    ! single-tapered complex synthetic displacement and velocity
    s1_dtw_ho_all(:,ictaper) = 0.
    s1_vtw_ho_all(:,ictaper) = 0.
    s2_dtw_ho_all(:,ictaper) = 0.
    s2_vtw_ho_all(:,ictaper) = 0.
    s1_dtw_ho_all(1:nlen,ictaper) = dcmplx(s1_dtw_h(1:nlen),0.)
    s1_vtw_ho_all(1:nlen,ictaper) = dcmplx(s1_vtw_h(1:nlen),0.)
    s2_dtw_ho_all(1:nlen,ictaper) = dcmplx(s2_dtw_h(1:nlen),0.)
    s2_vtw_ho_all(1:nlen,ictaper) = dcmplx(s2_vtw_h(1:nlen),0.)

    ! apply FFT get complex spectra
    call fft(LNPT,s1_dtw_ho_all(:,ictaper),dble(FORWARD_FFT),dble(deltat))
    call fft(LNPT,s1_vtw_ho_all(:,ictaper),dble(FORWARD_FFT),dble(deltat))
    call fft(LNPT,s2_dtw_ho_all(:,ictaper),dble(FORWARD_FFT),dble(deltat))
    call fft(LNPT,s2_vtw_ho_all(:,ictaper),dble(FORWARD_FFT),dble(deltat))

    Mtr1_mtm(:) = Mtr1_mtm(:) + s1_vtw_ho_all(:,ictaper) &
        * conjg(s2_vtw_ho_all(:,ictaper))
    Mtr2_mtm(:) = Mtr2_mtm(:) + s2_vtw_ho_all(:,ictaper) &
        * conjg(s1_vtw_ho_all(:,ictaper))

    Mtr3_mtm(:) = Mtr3_mtm(:) + s1_dtw_ho_all(:,ictaper) &
        * conjg(s2_dtw_ho_all(:,ictaper))
    Mtr4_mtm(:) = Mtr4_mtm(:) + s2_dtw_ho_all(:,ictaper) &
        * conjg(s1_dtw_ho_all(:,ictaper))
    Mtr5_mtm(:) = Mtr5_mtm(:) + s2_dtw_ho_all(:,ictaper) &
        * conjg(s2_dtw_ho_all(:,ictaper))
    enddo ! ictaper


    ! compute p_j, q_j, P_j, Q_j and adjoint source fp, fq
    fp1(1:npts) = 0.
    fq1(1:npts) = 0.
    fp2(1:npts) = 0.
    fq2(1:npts) = 0.

    do ictaper = 1,ntaper

    ! compute p_j(w) and q_j(w)
    p1_adj(:) = cmplx(0.,0.)
    pw1_adj(:) = cmplx(0.,0.)
    p2_adj(:) = cmplx(0.,0.)
    pw2_adj(:) = cmplx(0.,0.)
    q1_adj(:) = cmplx(0.,0.)
    qw1_adj(:) = cmplx(0.,0.)
    q2_adj(:) = cmplx(0.,0.)
    qw2_adj(:) = cmplx(0.,0.)

    fp1_adj(:) = cmplx(0.,0.)
    fp2_adj(:) = cmplx(0.,0.)
    fq1_adj(:) = cmplx(0.,0.)
    fq2_adj(:) = cmplx(0.,0.)

    do i = i_fstart,i_fend
    !p1_adj(i)  =  -0.5* conjg(s2_vtw_ho_all(i,ictaper)) / Mtr1_mtm(i)
    pw1_adj(i) = - 0.5 * s2_vtw_ho_all(i,ictaper) / Mtr2_mtm(i)

    !p2_adj(i)  =  0.5 * conjg(s1_vtw_ho_all(i,ictaper)) / Mtr2_mtm(i)
    pw2_adj(i) =  0.5 * s1_vtw_ho_all(i,ictaper) / Mtr1_mtm(i)

    q1_adj(i)  =  0.5 * conjg(s2_dtw_ho_all(i,ictaper)) / Mtr3_mtm(i)
    !   qw1_adj(i) =  0.5 * s2_dtw_ho_all(i,ictaper) / Mtr4_mtm(i)
    ! 
    q2_adj(i)  =  0.5 * conjg(s1_dtw_ho_all(i,ictaper)) / Mtr4_mtm(i) &
        - conjg(s2_dtw_ho_all(i,ictaper)) / Mtr5_mtm(i)
    !   qw2_adj(i) =  0.5 * s1_dtw_ho_all(i,ictaper) / Mtr3_mtm(i) &
    !                 - s2_dtw_ho_all(i,ictaper) / Mtr5_mtm(i)
    enddo

    ! compute P_j(w) and Q_j(w)
    ! NOTE: the MT measurement is incorporated here
    !             also note that wp_taper and wq_taper can contain
    !             uncertainty estimations
    ! adds misfit measurement dtau, dlnA
    fp1_adj(i_fstart: i_fend) &
        ! = (conjg(p1_adj(i_fstart: i_fend)) +
    =  2.0*(pw1_adj(i_fstart:i_fend)) &
        * dcmplx(ddtau_w(i_fstart: i_fend),0.) &
        * dcmplx(wp_taper(i_fstart: i_fend),0.)
    fp2_adj(i_fstart: i_fend) &
        ! = (conjg(p2_adj(i_fstart: i_fend)) +
    = 2.0*(pw2_adj(i_fstart: i_fend)) &
        * dcmplx(ddtau_w(i_fstart: i_fend),0.) &
        * dcmplx(wp_taper(i_fstart: i_fend),0.)
    fq1_adj(i_fstart: i_fend) &
        = (q1_adj(i_fstart: i_fend) + conjg(q1_adj(i_fstart: i_fend)))&
        * dcmplx(ddlnA_w(i_fstart: i_fend),0.) &
        * dcmplx(wq_taper(i_fstart: i_fend),0.)
    fq2_adj(i_fstart: i_fend) &
        = (q2_adj(i_fstart: i_fend) + conjg(q2_adj(i_fstart: i_fend)))&
        * dcmplx(ddlnA_w(i_fstart: i_fend),0.) &
        * dcmplx(wq_taper(i_fstart: i_fend),0.)


    ! IFFT into the time domain
    call fftinv(LNPT,fp1_adj,dble(REVERSE_FFT),dble(deltat),fp1_adj_t)
    call fftinv(LNPT,fp2_adj,dble(REVERSE_FFT),dble(deltat),fp2_adj_t)
    call fftinv(LNPT,fq1_adj,dble(REVERSE_FFT),dble(deltat),fq1_adj_t)
    call fftinv(LNPT,fq2_adj,dble(REVERSE_FFT),dble(deltat),fq2_adj_t)

    ! create adjoint source
    ! applies taper to time signal
    fp1(1:nlen) = fp1(1:nlen) + real(tas(1:nlen,ictaper) * &
        fp1_adj_t(1:nlen),CUSTOM_REAL)
    fp2(1:nlen) = fp2(1:nlen) + real(tas(1:nlen,ictaper) * &
        fp2_adj_t(1:nlen),CUSTOM_REAL)
    fq1(1:nlen) = fq1(1:nlen) + real(tas(1:nlen,ictaper) * &
        fq1_adj_t(1:nlen),CUSTOM_REAL)
    fq2(1:nlen) = fq2(1:nlen) + real(tas(1:nlen,ictaper) * &
        fq2_adj_t(1:nlen),CUSTOM_REAL)
    enddo  ! ictaper
    deallocate(tas)
    deallocate(s1_dtw_ho_all)
    deallocate(s1_vtw_ho_all)
    deallocate(s2_dtw_ho_all)
    deallocate(s2_vtw_ho_all)

end subroutine mtm_DD_adj

!==============================================================================
!        subroutines used in mtm_measure() and mtm_adj()
!==============================================================================

!-----------------------------------------------------------------------
subroutine CC_similarity(d1,d2,npts,&
        i_tstart1,i_tend1,i_tstart2,i_tend2,&
        window_type,&
        cc_max)
    !! measure the similarity of two waveforms based on cross-correlatiobs

    use constants
    implicit none

    ! inputs & outputs 
    real(kind=CUSTOM_REAL), dimension(*), intent(in) :: d1,d2
    integer, intent(in) :: i_tstart1,i_tend1,i_tstart2,i_tend2
    integer, intent(in) :: npts,window_type
    real(kind=CUSTOM_REAL), intent(out) :: cc_max

    ! window
    integer :: nlen1,nlen2,nlen
    real(kind=CUSTOM_REAL), dimension(npts) :: d1_tw,d2_tw
    ! cc 
    integer :: ishift
    real(kind=CUSTOM_REAL) :: dlnA

    !! initialization
    cc_max=0.0

    !! window
    call cc_window(d1,npts,window_type,i_tstart1,i_tend1,0,0.d0,nlen1,d1_tw)
    call cc_window(d2,npts,window_type,i_tstart2,i_tend2,0,0.d0,nlen2,d2_tw)
    if(nlen1<1 .or. nlen1>npts) print*,'check nlen1 ',nlen1
    if(nlen2<1 .or. nlen2>npts) print*,'check nlen2 ',nlen2
    nlen = max(nlen1,nlen2)
    !! cc
    call xcorr_calc(d1_tw,d2_tw,npts,1,nlen,ishift,dlnA,cc_max) !T(d1-d2)

end subroutine cc_similarity

! ---------------------------------------------------------------------------
subroutine frequency_limit (dat,nlen,deltat,i_fstart,i_fend)
    use constants
    implicit none

    integer, intent(in) :: nlen
    real(kind=CUSTOM_REAL), dimension(*), intent(in) :: dat
    real(kind=CUSTOM_REAL), intent(in) ::  deltat
    integer, intent(out) :: i_fstart,i_fend

    integer :: i,fnum 
    real(kind=CUSTOM_REAL) :: df,ampmax, wtr_use
    integer :: i_ampmax, i_fstart_stop, i_fend_stop
    complex (SIZE_DOUBLE), dimension(NPT) :: dat_dtwo     


    !! find reliable frequency limits 
    ! fft of untapered inputs 
    df = 1./(NPT*deltat)
    dat_dtwo = cmplx(0.,0.)
    dat_dtwo(1:nlen) = cmplx(real(dat(1:nlen),SIZE_DOUBLE),0.0_SIZE_DOUBLE)
    call fft(LNPT,dat_dtwo,dble(FORWARD_FFT),dble(deltat))

    fnum = NPT/2+1
    ! water level based untapered dat1
    ampmax = 0.
    i_ampmax = 1
    do i = 1, fnum   ! loop over frequencies
    if( abs(dat_dtwo(i)) > ampmax) then
        ampmax =  abs(dat_dtwo(i))
        i_ampmax = i
    endif
    enddo
    wtr_use = ampmax * WTR

    ! i_fend 
    i_fend = fnum
    i_fend_stop = 0
    do i = 1,fnum
    if( abs(dat_dtwo(i)) <= abs(wtr_use) .and. i_fend_stop==0 .and. i > i_ampmax ) then
        i_fend_stop = 1
        i_fend = i
    endif
    if( abs(dat_dtwo(i)) >= 10.*abs(wtr_use) .and. i_fend_stop==1 .and. i > i_ampmax) then
        i_fend_stop = 0
        i_fend = i
    endif
    enddo

    ! i_fstart 
    i_fstart = 1
    i_fstart_stop = 0
    do i = fnum,1,-1
    if( abs(dat_dtwo(i)) <= abs(wtr_use) .and. i_fstart_stop==0 .and. i < i_ampmax ) then
        i_fstart_stop = 1
        i_fstart = i
    endif
    if( abs(dat_dtwo(i)) >= 10.*abs(wtr_use) .and. i_fstart_stop==1 .and. i < i_ampmax) then
        i_fstart_stop = 0
        i_fstart = i
    endif
    enddo

    if( DISPLAY_DETAILS) then
        open(1,file=trim(output_dir)//'/spectrum',status='unknown')
        do  i =  1,2*i_fend
        write(1,'(f15.5,2e15.5)') (i-1)*df,abs(dat_dtwo(i)),wtr_use
        enddo
        close(1)
    endif

end subroutine frequency_limit
! ---------------------------------------------------------------------------
subroutine cc_error(data_dtw,syn_dtw,npts,deltat,nlen,ishift,dlnA,sigma_dt,sigma_dlnA)
    ! CHT: Estimate the uncertainty in the CC measurement
    !      based on the integrated waveform difference between the data
    !      and the reconstructed synthetics.
    ! NOTE: We implement the exact equations that are in the Latex notes.

    use constants
    implicit none
    real(kind=CUSTOM_REAL), dimension(*), intent(in) :: data_dtw, syn_dtw
    integer, intent(in) :: npts,nlen, ishift
    real(kind=CUSTOM_REAL), intent(in) :: deltat, dlnA
    real(kind=CUSTOM_REAL), intent(inout) :: sigma_dt, sigma_dlnA

    real(kind=CUSTOM_REAL), dimension(npts) :: syn_dtw_cc, syn_vtw_cc
    real(kind=CUSTOM_REAL) :: sigma_dt_top, sigma_dlnA_top, sigma_dt_bot, sigma_dlnA_bot
    integer :: i,j

    ! make corrections to syn based on cc measurement
    syn_dtw_cc(:)=0.0
    do i=1,nlen
    j=i-ishift 
    if(j>=1 .and. j <=nlen) then
        syn_dtw_cc(i) = syn_dtw(j) * exp(dlnA)
    endif
    enddo

    syn_vtw_cc(:)=0.0
    ! compute cc-corrected synthetic velocity
    call compute_vel(syn_dtw_cc,npts,deltat,nlen,syn_vtw_cc)

    ! estimated uncertainty in cross-correlation travltime and amplitude
    sigma_dt_top   = sum( (data_dtw(1:nlen) - syn_dtw_cc(1:nlen) )**2 )
    sigma_dt_bot   = sum( syn_vtw_cc(1:nlen)**2 )
    sigma_dlnA_top = sigma_dt_top
    sigma_dlnA_bot = sum( (syn_dtw_cc(1:nlen))**2 )/(dlnA * dlnA)
    sigma_dt       = sqrt( sigma_dt_top / sigma_dt_bot )
    sigma_dlnA     = sqrt( sigma_dlnA_top / sigma_dlnA_bot )

    ! make sure that the uncertainty estimates are not below the water level;
    ! otherwise, the adjoint sources will blow up unreasonably
    if( sigma_dt < DT_SIGMA_MIN) sigma_dt = DT_SIGMA_MIN
    if( sigma_dlnA < DLNA_SIGMA_MIN) sigma_dlnA = DLNA_SIGMA_MIN

    if(DISPLAY_DETAILS) then
        print*
        print*,'error estimation based on cc'
        print*,'sigma_dt top = ',sigma_dt_top, 'bot =', sigma_dt_bot, &
            ' ratio =', sqrt( sigma_dt_top / sigma_dt_bot )
        print*,'sigma_dlnA top = ',sigma_dlnA_top, 'bot =', sigma_dlnA_bot,&
            ' ratio =', sqrt( sigma_dlnA_top / sigma_dlnA_bot )
        print *, 'estimate sigma_dt   : ', sigma_dt
        print *, 'estimate sigma_dlnA : ', sigma_dlnA
        open(1,file=trim(output_dir)//'/cc_error.dat',status='unknown')
        do i = 1,nlen
        write(1,'(I5,4e15.5)') i, data_dtw(i), syn_dtw(i), syn_dtw_cc(i), syn_vtw_cc(i)
        enddo
        close(1)
    endif

end subroutine cc_error
! ---------------------------------------------------------------------------
subroutine write_phase(phase_func, wvec, i_right, tshift, phi_wt, dtau_wt)
    use constants
    implicit none

    ! The transfer function maps the data2 to the CC-corrected data1;

    complex (CUSTOM_REAL), dimension(*), intent(in) :: phase_func
    real(kind=CUSTOM_REAL), dimension(*), intent(in) :: wvec
    real(kind=CUSTOM_REAL), intent(in) ::  tshift
    integer, intent(in) ::  i_right
    real(kind=CUSTOM_REAL), dimension(*), intent(out) :: phi_wt, dtau_wt

    integer :: i, j
    real(kind=CUSTOM_REAL) :: smth, smth1, smth2 ! f0

    phi_wt(1:NPT) = 0.

    ! loop to calculate phase and amplitude
    do i = 1, i_right
    phi_wt(i) = atan2( aimag(phase_func(i)) , real(phase_func(i)) )
    enddo

    ! NOTE: the CC measurements dT (tshift) and dlnA are BOTH included
    dtau_wt(1) = tshift
    do i = 1, i_right

    if (i > 1 .and. i < i_right) then
        ! check the smoothness (2nd-order derivative) by 2*pi changes
        smth  =  phi_wt(i+1) + phi_wt(i-1) - 2.0 * phi_wt(i)
        smth1 = (phi_wt(i+1) + TWOPI) + phi_wt(i-1) - 2.0 * phi_wt(i)
        smth2 = (phi_wt(i+1) - TWOPI) + phi_wt(i-1) - 2.0 * phi_wt(i)
        if(abs(smth1).lt.abs(smth).and.abs(smth1).lt.abs(smth2).and. abs(phi_wt(i) - phi_wt(i+1)) > PHASE_STEP)then
            if (DISPLAY_DETAILS) print *, 'phase correction : 2 pi', i, phi_wt(i) - phi_wt(i+1)
            do j = i+1, i_right
            phi_wt(j) = phi_wt(j) + TWOPI
            enddo
        endif
        if(abs(smth2).lt.abs(smth).and.abs(smth2).lt.abs(smth1).and. abs(phi_wt(i) - phi_wt(i+1)) > PHASE_STEP)then
            if (DISPLAY_DETAILS) print *, 'phase correction : - 2 pi', i, phi_wt(i) - phi_wt(i+1)
            do j = i+1, i_right
            phi_wt(j) = phi_wt(j) - TWOPI
            enddo
        endif
    endif

    ! add the CC measurements to the transfer function
    if (i > 1) dtau_wt(i) = (0.5/wvec(i)) * phi_wt(i) + tshift

    enddo

end subroutine write_phase
! ---------------------------------------------------------------------------
subroutine write_trans(trans, wvec,i_left, i_right, tshift, dlnA,phi_wt, abs_wt, dtau_wt, dlnA_wt)
    ! The transfer function maps the synthetics to the CC-deconstructed data;
    ! the CC measurements then need to be applied to match the original data.

    use constants
    implicit none
    complex (CUSTOM_REAL), dimension(*),intent(in) :: trans
    real(kind=CUSTOM_REAL), dimension(*),intent(in) :: wvec
    real(kind=CUSTOM_REAL), intent(in) :: tshift, dlnA
    integer, intent(in) :: i_left, i_right
    real(kind=CUSTOM_REAL), dimension(*), intent(out) :: phi_wt, abs_wt,dtau_wt,dlnA_wt

    integer :: i

    ! initialization 
    phi_wt(1:NPT) = 0.d0 
    abs_wt(1:NPT) = 0.d0
    dtau_wt(1:NPT) = 0.d0
    dlnA_wt(1:NPT) = 0.d0 

    ! loop to calculate phase and amplitude
    do i = i_left, i_right
    phi_wt(i) = atan2( aimag(trans(i)) , real(trans(i)) )
    abs_wt(i) = abs(trans(i))
    enddo

    ! check the smoothenss of phi 
    call smooth_phase(phi_wt,i_left,i_right)

    ! NOTE: the CC measurements dT (tshift) and dlnA are BOTH included 
    ! print*,' add the CC measurements to the transfer function'
    if (i_left > 1) dtau_wt(i_left:i_right) = (-1./wvec(i_left:i_right)) * phi_wt(i_left:i_right) + tshift
    if (i_left == 1) dtau_wt(i_left+1:i_right) = (-1./wvec(i_left+1:i_right)) * phi_wt(i_left+1:i_right) + tshift
    dtau_wt(1) = tshift 
    dlnA_wt(i_left:i_right) = log(abs_wt(i_left:i_right)) + dlnA

end subroutine write_trans
! --------------------------------------------------------------------
subroutine smooth_phase(phi_wt,i_left,i_right)
    ! check the smoothness of phi
    use constants
    implicit none

    real(kind=CUSTOM_REAL), dimension(*),intent(out) :: phi_wt
    integer, intent(in) :: i_left, i_right

    integer :: i, j
    real(kind=CUSTOM_REAL) :: smth, smth1, smth2 ! f0

    do i = i_left, i_right
    if (i > i_left .and. i < i_right) then
        !  print*,' check the smoothness (2nd-order derivative) by 2*pi changes'
        smth  =  phi_wt(i+1) + phi_wt(i-1) - 2.0 * phi_wt(i)
        smth1 = (phi_wt(i+1) + TWOPI) + phi_wt(i-1) - 2.0 * phi_wt(i)
        smth2 = (phi_wt(i+1) - TWOPI) + phi_wt(i-1) - 2.0 * phi_wt(i)
        if(abs(smth1).lt.abs(smth).and.abs(smth1).lt.abs(smth2).and. abs(phi_wt(i) - phi_wt(i+1)) > PHASE_STEP)then
            if (DISPLAY_DETAILS) print *, 'phase correction : 2pi',i,phi_wt(i) - phi_wt(i+1)
            do j = i+1, i_right
            phi_wt(j) = phi_wt(j) + TWOPI
            enddo
        endif

        if(abs(smth2).lt.abs(smth).and.abs(smth2).lt.abs(smth1).and. abs(phi_wt(i) - phi_wt(i+1)) > PHASE_STEP)then
            if (DISPLAY_DETAILS) print *, 'phase correction : - 2 pi ', i, phi_wt(i) - phi_wt(i+1)
            do j = i+1, i_right
            phi_wt(j) = phi_wt(j) - TWOPI
            enddo
        endif

    endif 
    enddo

end subroutine smooth_phase
! --------------------------------------------------------------------
subroutine unwrap_phase(phi_wt,i_left,i_right)
    ! unwrap phi
    use constants
    implicit none

    real(kind=CUSTOM_REAL), dimension(*),intent(out) :: phi_wt
    integer, intent(in) :: i_left, i_right

    integer :: i
    real(kind=CUSTOM_REAL) :: dp, dps, dp_corr

    do i = i_left, i_right
    dp_corr = 0.0
    if (i > i_left .and. i <= i_right) then
        !% Incremental phase variations 
        dp = phi_wt(i)-phi_wt(i-1)
        !% Equivalent phase variations in [-pi,pi)
        dps = mod(dp+PI,TWOPI)-PI
        !% Preserve variation sign for pi vs. -pi
        if(dps==-PI .and. dp>0.0) dps=PI
        !% Incremental phase corrections
        dp_corr = dps-dp
        !% Ignore correction when incr. variation is < CUTOFF
        if(abs(dp)<CUTOFF) dp_corr=0.0
        !% Integrate corrections and add to P to produce smoothed phase values
        phi_wt(i:i_right) = phi_wt(i:i_right) + dp_corr
    endif
    enddo

end subroutine unwrap_phase
! --------------------------------------------------------------------
subroutine time_costaper(dat,syn,npts,istart,iend)
    use constants
    implicit none

    ! input parameters
    integer, intent(in) :: istart,iend,npts
    real(kind=CUSTOM_REAL), dimension(*), intent(out) :: dat,syn

    integer :: B1,B2,B3,B4,nlen_taper
    integer i, n_width_left,n_width_right,nlen
    real(kind=CUSTOM_REAL) :: fac

    ! set the number of points corresponding to the taper width
    if(istart<iend) then
        nlen_taper=floor((iend-istart)/4.0)
        B1=max(istart-nlen_taper,1)
        B4=min(iend+nlen_taper,npts)
        B2=istart
        B3=iend
    endif
    n_width_left=B2-B1
    n_width_right=B4-B3

    ! left side
    nlen=2*n_width_left
    do i = 1,n_width_left
    !fac = 1.                                         ! boxcar window
    !fac = 1 - sfac1*((i-1) - real(nlen)/2.)**2       ! welch window
    fac = 1. - cos(PI*(i-1)/(nlen-1))**ipwr_t
    dat(B1+i-1)=dat(B1-1+i)*fac
    syn(B1+i-1)=syn(B1-1+i)*fac
    enddo

    ! right side
    nlen=2*n_width_right
    do i = 1,n_width_right
    !fac = 1.                                         ! boxcar window
    !fac = 1 - sfac1*((i-1) - real(nlen)/2.)**2       ! welch window
    fac = 1. - cos(PI*(i-1)/(nlen-1))**ipwr_t
    dat(B4-i+1)=dat(B4-i+1)*fac
    syn(B4-i+1)=syn(B4-i+1)*fac
    enddo

    ! beyond the window is zero
    dat(1:B1)=0.0;
    dat(B4:npts)=0.0
    syn(1:B1)=0.0;
    syn(B4:npts)=0.0

end subroutine time_costaper
! ---------------------------------------------------------------------------
subroutine cc_window(dat,npts,window_type,istart,iend,ishift,dlnA,nlen,dat_win)
    ! delay by ishift and scale by exp(dlnA)
    use constants
    implicit none

    real(kind=CUSTOM_REAL), dimension(*), intent(in) :: dat
    real(kind=CUSTOM_REAL), intent(in) :: dlnA
    integer, intent(in) :: npts, istart,iend, ishift,window_type
    integer, intent(out) :: nlen
    real(kind=CUSTOM_REAL), dimension(*), intent(out) :: dat_win

    integer ::  i, j
    real(kind=CUSTOM_REAL) :: sfac1,  fac

    !print*
    !print*,'window with corrections: isfhit = ', ishift, ' dlnA = ',dlnA,&
    !'window type -- ',window_type

    nlen = iend - istart+1

    ! some constants
    sfac1 = (2./real(nlen))**2   ! for Welch window

    ! initialization
    dat_win(1:npts) = 0.d0

    do i = 1, nlen
    if(window_type ==2) then 
        fac = 1 - sfac1*((i-1) - real(nlen)/2.)**2   
    elseif(window_type ==3) then 
        fac = 1. - cos(PI*(i-1)/(nlen-1))**ipwr_t    
    elseif(window_type ==4) then 
        fac = 0.5 - 0.5*cos(TWOPI*(i-1)/(nlen-1))    
    else 
        fac = 1. ! boxcar window
    endif

    ! index 
    j = i + (istart-1) - ishift

    if(j>=1 .and. j <=npts) then
        dat_win(i) = dat(j) * exp(dlnA) * fac ! shift and scale
    endif

    enddo

end subroutine cc_window
! ---------------------------------------------------------------------------
subroutine cc_window_inverse(dat_win,npts,window_type,istart,iend,ishift,dlnA,dat)
    use constants
    implicit none

    real(kind=CUSTOM_REAL), dimension(*), intent(in) :: dat_win
    real(kind=CUSTOM_REAL), intent(in) :: dlnA
    integer, intent(in) :: npts, istart,iend, ishift,window_type
    real(kind=CUSTOM_REAL), dimension(*), intent(out) :: dat

    integer ::  i, j, nlen
    real :: sfac1, fac

    !  print*
    !  print*,'window with corrections: isfhit = ', ishift, ' dlnA = ',dlnA

    nlen = iend - istart

    ! some constants
    sfac1 = (2./real(nlen))**2   ! for Welch window

    ! initialization
    dat(1:npts) = 0.d0

    do i = 1, nlen
    if(window_type ==2) then 
        fac = 1 - sfac1*((i-1) - real(nlen)/2.)**2
    else if(window_type ==3) then  
        fac = 1. - cos(PI*(i-1)/(nlen-1))**ipwr_t 
    else if(window_type ==4) then 
        fac = 0.5 - 0.5*cos(TWOPI*(i-1)/(nlen-1))
    else 
        fac = 1. ! boxcar window
    endif

    j = i + (istart-1) - ishift

    ! window dat2 using exact window info 
    if(j>=1 .and. j <=npts) then
        dat(j) = dat_win(i) * exp(-dlnA) * fac
    endif
    enddo
end subroutine cc_window_inverse
! ---------------------------------------------------------------------------
subroutine compute_cc(dat1, dat2, nlen, dt, ishift, tshift, dlnA, cc_max)
    use constants
    implicit none

    ! time shift MEASUREMENT between data1 and data2
    ! CHT: modified the subroutine to resemble the one used in FLEXWIN
    ! tshift = T(dat1) - T(dat2)
    ! dlnA = 0.5 * log(A(dat1)**2 / A(dat2)**2)

    real(kind=CUSTOM_REAL), dimension(*), intent(in) :: dat1, dat2
    integer, intent(in) :: nlen
    real(kind=CUSTOM_REAL), intent(in) :: dt
    real(kind=CUSTOM_REAL), intent(out) :: tshift, dlnA, cc_max
    integer, intent(out) :: ishift

    real(kind=CUSTOM_REAL) :: cc, norm_s, norm ! cr_shift
    integer i1, i2, i, j, i_left, i_right, id_left, id_right
    ! initialise shift and cross correlation to zero
    ishift = 0
    cc_max = 0.0

    ! index of window limits
    i1 = 1
    i2 = nlen

    ! length of window (number of points, including ends)
    !nlen = i2 - i1 + 1

    ! power of synthetic signal in window
    norm_s = sqrt(sum(dat2(i1:i2)*dat2(i1:i2)))

    ! left and right limits of index (time) shift search
    ! NOTE: This looks OUTSIDE the time window of interest to compute TSHIFT and
    ! CC.
    !       How far to look outside, in theory, should be another parameter.
    !       However, it does not matter as much if the data and synthetics are
    !          zeroed outside the windows.
    i_left = - nlen
    i_right = nlen

    ! i is the index to shift to be applied to DATA (data)
    do i = i_left, i_right

    ! normalization factor varies as you take different windows of data
    id_left = max(1,i1+i)      ! left index for data window
    id_right = min(nlen,i2+i)  ! right index for data window
    norm = norm_s * sqrt(sum(dat1(id_left:id_right)*(dat1(id_left:id_right))))

    ! cc as a function of i
    cc = 0.
    do j = i1, i2   ! loop over full window length
    if((j+i).ge.1 .and. (j+i).le.nlen) cc = cc + dat2(j)*dat1(j+i)  ! d ishifted by i
    enddo
    cc = cc/norm

    if (cc > cc_max) then
        ! CHT: do not allow time shifts larger than the specified input range
        ! This is an important criterion, since it may pick TSHIFT_MIN or
        ! TSHIFT_MAX
        ! if cc_max within the interval occurs on the boundary.
        !         if( (i*dt >= TSHIFT_MIN).and.(i*dt <= TSHIFT_MAX) ) then
        cc_max = cc
        ishift = i
        !         endif
    endif

    enddo
    tshift = ishift*dt
    ! The previously used expression for dlnA of Dahlen and Baig (2002),
    ! is a first-order perturbation of ln(A1/A2) = (A1-A2)/A2 .
    ! The new expression is better suited to getting Gaussian-distributed
    ! values between -1 and 1, with dlnA = 0 indicating perfect fit, as before.    
    dlnA = 0.5 * log( sum(dat1(i1:i2) * dat1(i1:i2)) / sum(dat2(i1:i2) * dat2(i1:i2)) )

end subroutine compute_cc
!-----------------------------------------------------------------------------
subroutine compute_vel(syn,npts,deltat,nlen,syn_vel)
    use constants
    implicit none 
    real(kind=CUSTOM_REAL), dimension(*),intent(in) :: syn 
    real(kind=CUSTOM_REAL), intent(in) :: deltat 
    integer, intent(in) :: npts, nlen
    real(kind=CUSTOM_REAL), dimension(*),intent(out) :: syn_vel

    ! index
    integer :: itime 

    ! initialization 
    syn_vel(1:npts)=0.0

    do itime=2,nlen-1
    syn_vel(itime)=(syn(itime+1)-syn(itime-1))/(2*deltat)
    enddo
    ! boundaries
    syn_vel(1)=(syn(2)-syn(1))/deltat
    syn_vel(nlen)=(syn(nlen)-syn(nlen-1))/deltat

end subroutine compute_vel
!-----------------------------------------------------------------------------

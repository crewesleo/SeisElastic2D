program misfit_adjoint
    ! to calculate adjoint source and misfit

#ifdef USE_MPI
    use mpi
#endif

    use seismo_parameters
    implicit none

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer, parameter :: NARGS = 7
    integer ::  ier,i,icomp,irec,send_data_tag,receive_data_tag
    real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: adj_master 
    character(len=MAX_STRING_LEN) :: data_names(MAX_DATA_NUM)
    character(len=MAX_STRING_LEN) :: data_names_comma_delimited
    character(len=MAX_STRING_LEN) :: arg(NARGS)
    character(len=MAX_STRING_LEN) :: input_dir
    character(len=MAX_FILENAME_LEN) :: filename
    real t1,t2

#ifdef USE_MPI
    call MPI_INIT(ier)
    call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ier)
    call MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ier)
#else
    nproc = 1
    myrank = 0
#endif

    if (DISPLAY_DETAILS .and. myrank==0)&
        print *,"Running misfit_adjoint.exe on",NPROC,"processors"
    call cpu_time(t1)

    ! parse command line arguments
    if (command_argument_count() /= NARGS) then
        if (myrank == 0) then
            print *, 'USAGE:  mpirun -np NPROC bin/misfit_adjoint.exe ... ' 
            stop ' Please check command line arguments'
        endif
    endif

#ifdef USE_MPI
    call MPI_BARRIER(MPI_COMM_WORLD,ier)
#endif

    ! get input parameters
    do i = 1, NARGS
    call get_command_argument(i,arg(i), status=ier)
    enddo
    read(arg(1),*) compute_adjoint
    data_names_comma_delimited = arg(2)
    measurement_list=arg(3)
    misfit_type_list=arg(4)
    input_dir=arg(5)
    read(arg(6),*) VISCOELASTIC
    measurement_attenuation=arg(7)
    call split_string(data_names_comma_delimited,',',data_names,ndata)

    !! gloabl initialization
    misfit=0.0_CUSTOM_REAL

    ! loop over comp
    do icomp=1,ndata

    ! step 1 -- load data and syn
    call initialize(input_dir,adjustl(data_names(icomp)))

    ! exist data/syn
    if(nrec_proc > 0 ) then  ! non-zero file

        ! step 2 -- preprocessing data and syn 
        if(maxval(abs(seism_obs(:,:)))<SMALL_VAL) cycle
        call process_obs(seism_obs)
        if(maxval(abs(seism_syn(:,:)))<SMALL_VAL) cycle
        call process_syn(seism_syn)

        !!! Normalization by PWY 02-08-2018
!        seism_obs=seism_obs/maxval(abs(seism_obs(:,:)))
!        seism_syn=seism_syn/maxval(abs(seism_syn(:,:)))

        if(compute_adjoint .and. DISPLAY_DETAILS) then
            write(filename,'(a)') &
                trim(input_dir)//'/DATA_obs/'//trim(adjustl(data_names(icomp)))//'_processed.bin'
            print*,'SAVE processed seismic_obs -- ',trim(filename)
            open(unit=IOUT,file=trim(filename),status='unknown',form='unformatted',iostat=ier)
            if (ier /= 0) then
                print*, 'Error: could not open model file: ',trim(filename)
                stop 'Error reading neighbors external mesh file'
            endif
            write(IOUT) seism_obs
            close(IOUT)
            write(filename,'(a)') &
                trim(input_dir)//'/DATA_syn/'//trim(adjustl(data_names(icomp)))//'_processed.bin'
            print*,'SAVE processed seismic_syn -- ',trim(filename)
            open(unit=IOUT,file=trim(filename),status='unknown',form='unformatted',iostat=ier)
            if (ier /= 0) then
                print*, 'Error: could not open model file: ',trim(filename)
                stop 'Error reading neighbors external mesh file'
            endif
            write(IOUT) seism_syn
            close(IOUT)
        endif

        ! step 3 -- misfit-adjoint
        num_AD=0
        num_DD=0
        misfit_AD=0.0_CUSTOM_REAL
        misfit_DD=0.0_CUSTOM_REAL

        ! absolute mwasurement  
        if(index(misfit_type_list,'AD')>0) then 
            call Absolute_diff()

            if(DISPLAY_DETAILS .and. compute_adjoint) then
                print*
                print*, 'summed squared misfit_',trim(measurement_list),&
                    '_AD=',misfit_AD  
                print*, trim(adjustl(data_names(icomp))), &
                    ': total number of AD measurements :', &
                    num_AD, num_AD*100.0/nrec_proc,'%'
            endif
        endif

        ! differential measurement
        if (nrec_proc>1 .and. index(misfit_type_list,'DD')>0) then
            call Relative_diff(input_dir,adjustl(data_names(icomp)))

            if(DISPLAY_DETAILS .and. compute_adjoint) then
                print*
                print*, 'summed squared misfit_',trim(measurement_list),&
                    '_DD=',misfit_DD
                print*, trim(adjustl(data_names(icomp))), &
                    ': total number of DD measurements :', &
                    num_DD, num_DD*100.0/(nrec_proc*(nrec_proc-1)/2),'%'
            endif
        endif

        ! misfit (AD+DD+...)
        misfit=misfit+ misfit_AD/max(num_AD,1) + misfit_DD/max(num_DD,1)

        ! adjoint
        if(compute_adjoint) then 
            seism_adj = seism_adj_AD/max(num_AD,1) + seism_adj_DD/max(num_DD,1)
        endif

        call finalize(input_dir,adjustl(data_names(icomp)))

    endif ! nrec_proc>0
    enddo ! icomp

    ! step 5 -- save misfit
    write(filename,'(a,i6.6,a)') &
        trim(input_dir)//'/proc',myrank,'_misfit.dat'
    OPEN (IOUT, FILE=trim(filename),status='unknown',iostat = ier)
    if(ier/=0) then
        print*,'Error opening data misfit file: ',trim(filename)
        stop
    else
        write(IOUT,*) misfit
    endif
    close(IOUT)
    if(DISPLAY_DETAILS .and. compute_adjoint .and. nrec_proc>0) then
        print*
        print*,'SAVE misfit -- ',trim(filename)
        print*,'myrank=',myrank,' final misfit_',trim(measurement_list), &
            '_',trim(misfit_type_list), '= ',misfit
    endif

    call cpu_time(t2)
    if(DISPLAY_DETAILS .and. compute_adjoint .and. myrank==0) &
        print *,'Computation time with CPU:',t2-t1

#ifdef USE_MPI
    ! stop all the processes and exit
    call MPI_FINALIZE(ier)
#endif

end program misfit_adjoint
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine initialize(directory,data_name)
    use seismo_parameters
    implicit none
    integer :: ier,iker,imod,irec
    integer :: filesize, nrec_obs, nrec_syn
    logical :: ex_obs, ex_syn
    character(len=MAX_FILENAME_LEN) :: char_i, filen, filename_obs, filename_syn
    character(len=MAX_STRING_LEN) :: directory
    character(len=MAX_STRING_LEN) :: data_name
    character(len=MAX_FILENAME_LEN) :: stf_filename,filename
    logical :: ex_stf
    real(kind=CUSTOM_REAL), dimension(:,:),allocatable :: temp
    integer :: irow, ncolum
    real(kind=CUSTOM_REAL) :: dis_sr

    !! data file format
    filen='empty'
    select case(solver)
    case('specfem2D')   ! single file
        if(myrank==0 .and. data_name == 'x') filen='Ux_file_single.su'
        if(myrank==0 .and. data_name == 'y') filen='Uy_file_single.su'
        if(myrank==0 .and. data_name == 'z') filen='Uz_file_single.su'
        if(myrank==0 .and. data_name == 'p') filen='Up_file_single.su'

    case('specfem3D')  
        write(char_i, '(I5)') myrank        ! convert integer to char
        if(data_name == 'x') write(filen,'(a)') trim(adjustl(char_i))//'_dx_SU'
        if(data_name == 'y') write(filen,'(a)') trim(adjustl(char_i))//'_dy_SU'
        if(data_name == 'z') write(filen,'(a)') trim(adjustl(char_i))//'_dz_SU'
        if(data_name == 'p') write(filen,'(a)') trim(adjustl(char_i))//'_dp_SU'

    case default
        print*,'Currently only work for specfem2D and specfem3D solver'
        stop
    end select

    !! load data & syn
    nrec_proc=0
    write(filename_obs,'(a)') &
        trim(directory)//'/DATA_obs/'//trim(filen)
    write(filename_syn,'(a)') &
        trim(directory)//'/DATA_syn/'//trim(filen)
    !exist and size
    inquire (file=trim(filename_obs), exist=ex_obs)
    inquire (file=trim(filename_syn), exist=ex_syn)

    if(ex_obs .and. ex_syn) then ! exist
        inquire(file=trim(filename_obs), size=filesize)
        nrec_obs=filesize/(240+4*NSTEP) 
        inquire(file=trim(filename_syn), size=filesize)
        nrec_syn=filesize/(240+4*NSTEP)
        if(nrec_obs == nrec_syn) then
            nrec_proc=nrec_obs
            if(DISPLAY_DETAILS) then
                print*,'myrank=',myrank,' LOAD nrec_proc traces : ',nrec_proc
                print*,'seism_obs -- ',trim(filename_obs)
                print*,'seism_syn -- ',trim(filename_syn)
            endif
        else
            print*,'size for obs and syn file is not the same : ',nrec_obs, nrec_syn
            stop
        endif

        !! allocate 
        allocate(seism_obs(NSTEP,nrec_proc))
        allocate(seism_syn(NSTEP,nrec_proc))
        allocate(seism_adj(NSTEP,nrec_proc))
        allocate(seism_adj_AD(NSTEP,nrec_proc))
        allocate(seism_adj_DD(NSTEP,nrec_proc))
        allocate(st_xval(nrec_proc))
        allocate(st_yval(nrec_proc))
        allocate(st_zval(nrec_proc))
        allocate(win_start(nrec_proc))
        allocate(win_end(nrec_proc))
        !allocate(which_proc_receiver(NPROC))

        ! initialization
        seism_obs = 0.0_CUSTOM_REAL
        seism_syn = 0.0_CUSTOM_REAL
        seism_adj = 0.0_CUSTOM_REAL
        seism_adj_AD = 0.0_CUSTOM_REAL
        seism_adj_DD = 0.0_CUSTOM_REAL
        win_start=NSTEP
        win_end=1
        ! which_proc_receiver=0

        ! allocate obs and syn
        call readSU(filename_obs,seism_obs)
        call readSU(filename_syn,seism_syn)

        ! normalization modified by PWY
!        seism_obs=seism_obs/maxval(abs(seism_obs(:,:)))
!        seism_syn=seism_syn/maxval(abs(seism_syn(:,:)))

!        call writeSU(filename_obs,seism_obs)
!        call writeSU(filename_syn,seism_syn)

        if(DISPLAY_DETAILS) then
            print*,'Min / Max of seism_obs : ',&
                minval(seism_obs(:,:)),maxval(seism_obs(:,:))
            print*,'Min / Max of seism_syn : ',&
                minval(seism_syn(:,:)),maxval(seism_syn(:,:))
        endif
    endif ! exist

end subroutine initialize
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine process_obs(seism)
    use seismo_parameters
    integer :: ier,irec
    real(kind=CUSTOM_REAL) :: seism(NSTEP,nrec_proc)
    real(kind=CUSTOM_REAL) :: trace(NSTEP)
    real(kind=CUSTOM_REAL) :: wtr
    integer :: itstart,itend
    integer :: itime
    integer :: NA
    real(kind=CUSTOM_REAL) :: time(NSTEP),tas(NSTEP)
    real(kind=CUSTOM_REAL) :: dis_sr

    do itime=1,NSTEP
    time=(itime-1)*deltat+t0
    enddo

    do irec = 1,nrec_proc
    trace(:) = seism(:,irec)
    ! initialization
    istart=1
    iend=NSTEP

    wtr=maxval(abs(trace(istart:iend)))
    if(wtr>SMALL_VAL) then  ! non-zero trace
        dis_sr=sqrt((st_xval(irec)-x_source)**2 &
            +(st_yval(irec)-y_source)**2 &
            +(st_zval(irec)-z_source)**2)

        !! mute 
        if(mute_near .eq. 1 .and. dis_sr<=offset_near ) then
            trace(1:NSTEP) = 0.0
            istart=NSTEP
            iend=0
        endif
        if(mute_far .eq. 1 .and. dis_sr>=offset_far) then
            trace(1:NSTEP) =0.0
            istart=NSTEP
            iend=0
        endif
        !! laplace damping spatially and temporally 
        if (is_laplace .eq. 1 .and. istart<iend) then
            ! spatial
            trace=trace*exp(-(dis_sr*lambda_x)**2)
            ! temporal
            trace=trace*exp(-(time*lambda_t)**2)
        endif
        !! dip-window using slopes 
        if(is_window .eq. 1 .and. istart<iend) then
            istart=max(int((t0+dis_sr/Vmax-taper_len)/deltat),istart)
            iend=min(int((t0+dis_sr/Vmin+wavelet_len+taper_len)/deltat),iend)
            tas(1:NSTEP)=0.d0
            if(iend>istart) then
                call window(NSTEP,istart,iend,window_type,tas)
            else
                istart=NSTEP
                iend=1
            endif
            trace=trace*tas
        endif ! window
        !! WT filtering
        if( Wscale .gt. 0.and. istart<iend) then
            call WT(trace,NSTEP,Wscale,NA)
            !istart=max(istart-NA,1)
            !iend=min(iend+NA,NSTEP)
        endif

        win_start(irec)=istart
        win_end(irec)=iend

    endif
    seism(:,irec)=trace(:)
    enddo
!!!! normalization
!    seism=seism/maxval(abs(seism(:,:)))
end subroutine process_obs
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine process_syn(seism)
    use seismo_parameters
    integer :: ier,irec
    real(kind=CUSTOM_REAL) :: seism(NSTEP,nrec_proc)
    real(kind=CUSTOM_REAL) :: trace(NSTEP)
    real(kind=CUSTOM_REAL) :: wtr
    integer :: itstart,itend
    integer :: itime
    integer :: NA
    real(kind=CUSTOM_REAL) :: time(NSTEP),tas(NSTEP)
    real(kind=CUSTOM_REAL) :: dis_sr

    do itime=1,NSTEP
    time=(itime-1)*deltat+t0
    enddo

    do irec = 1,nrec_proc 
    trace(:) = seism(:,irec)
    ! initialization
    istart=1
    iend=NSTEP

    wtr=maxval(abs(trace(istart:iend)))
    if(wtr>SMALL_VAL) then  ! non-zero trace
        dis_sr=sqrt((st_xval(irec)-x_source)**2 &
            +(st_yval(irec)-y_source)**2 &
            +(st_zval(irec)-z_source)**2)

        !! mute 
        if(mute_near .eq. 1 .and. dis_sr<=offset_near ) then
            trace(1:NSTEP) = 0.0
            istart=NSTEP
            iend=0
        endif
        if(mute_far .eq. 1 .and. dis_sr>=offset_far) then
            trace(1:NSTEP) =0.0
            istart=NSTEP
            iend=0
        endif
        !! laplace damping spatially and temporally 
        if (is_laplace .eq. 1 .and. istart<iend) then
            ! spatial
            trace=trace*exp(-(dis_sr*lambda_x)**2)
            ! temporal
            trace=trace*exp(-(time*lambda_t)**2)
        endif
        ! dip-window using slopes 
        if(is_window .eq. 1 .and. istart<iend) then
            istart=max(int((t0+dis_sr/Vmax-taper_len)/deltat),istart)
            iend=min(int((t0+dis_sr/Vmin+wavelet_len+taper_len)/deltat),iend)
            tas(1:NSTEP)=0.d0
            if(iend>istart) then
                call window(NSTEP,istart,iend,window_type,tas)
            else
                istart=NSTEP
                iend=1
            endif
            trace=trace*tas 
        endif ! window
        ! WT filtering
        if( Wscale .gt. 0.and. istart<iend) then
            call WT(trace,NSTEP,Wscale,NA)
            !istart=max(istart-NA,1)
            !iend=min(iend+NA,NSTEP)
        endif

        win_start(irec)=istart
        win_end(irec)=iend

    endif ! non-zero trace
    seism(:,irec)=trace(:)
    enddo ! irec
!!!! normalization
!    seism=seism/maxval(abs(seism(:,:)))
end subroutine process_syn
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine process_adj(trace,istart,iend,dis_sr)
    use seismo_parameters
    real(kind=CUSTOM_REAL) :: trace(NSTEP)
    integer :: itstart,itend
    real(kind=CUSTOM_REAL) :: dis_sr
    real(kind=CUSTOM_REAL) :: wtr
    integer :: NA
    integer :: itime
    real(kind=CUSTOM_REAL) :: time(NSTEP),tas(NSTEP)
    real(kind=CUSTOM_REAL), dimension(:),allocatable :: stf_reverse

    do itime=1,NSTEP
    time=(itime-1)*deltat+t0
    enddo

    wtr=maxval(abs(trace(istart:iend)))
    if(wtr>SMALL_VAL) then  ! non-zero trace
        ! WT filtering
        if( Wscale .gt. 0) then
            call WT(trace,NSTEP,Wscale,NA)
        endif
        ! dip-window using slopes 
        if(is_window .eq. 1) then
            ! window
            tas(1:NSTEP)=0.d0
            call window(NSTEP,istart,iend,window_type,tas)
            trace=trace*tas
        endif
        !! laplace damping spatially and temporally 
        if (is_laplace .eq. 1 .and. istart<iend) then
            ! spatial
            trace=trace*exp(-(dis_sr*lambda_x)**2)
            ! temporal
            trace=trace*exp(-(time*lambda_t)**2)
        endif
        !! mute 
        if(mute_near .eq. 1 .and. dis_sr<=offset_near ) then
            trace(1:NSTEP) = 0.0
        endif
        if(mute_far .eq. 1 .and. dis_sr>=offset_far) then                                                                        
            trace(1:NSTEP) = 0.0
        endif

    endif ! non-zero trace

end subroutine process_adj
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  `  
subroutine Absolute_diff()
    use seismo_parameters
    integer :: irec, itype, ntype
    real(kind=CUSTOM_REAL) :: seism(NSTEP,nrec_proc)
    real(kind=CUSTOM_REAL) :: d(NSTEP),s(NSTEP),adj_trace(NSTEP),adj(NSTEP)
    real(kind=CUSTOM_REAL) :: wtr
    real(kind=CUSTOM_REAL) :: misfit_value, misfit_trace
    integer :: ntstart,ntend,nlen
    integer :: ishift
    real(kind=CUSTOM_REAL) :: dlnA,cc_max
    character(len=MAX_STRING_LEN) :: measurement_types(MAX_MISFIT_TYPE)
    character(len=2) :: measurement_type
    character(len=3) :: measurement_attenuation_type
    character, parameter :: delimiter='+'

    ! initialization
    d = 0.0_CUSTOM_REAL
    s = 0.0_CUSTOM_REAL

    call split_string(measurement_list,delimiter,measurement_types,ntype)

    do irec=1,nrec_proc

    ! get data 
    d(:)=seism_obs(:,irec)
    s(:)=seism_syn(:,irec)
    !s(:)=0.d0
    !! misfit and adjoint evaluation
    ! initialization
    adj_trace = 0.0_CUSTOM_REAL
    misfit_trace=0.0

    dis_sr=sqrt((st_xval(irec)-x_source)**2 &
        +(st_yval(irec)-y_source)**2 &
        +(st_zval(irec)-z_source)**2)

    ! window info
    ntstart = win_start(irec)
    ntend= win_end(irec)

    !! quality control (window, norm, cc similarity)
    if (ntstart>=ntend ) cycle
    if (norm2(d(ntstart:ntend),1)<SMALL_VAL .or. &
        norm2(s(ntstart:ntend),1)<SMALL_VAL .or. &
        (ntend-ntstart)*deltat<window_len) &
        cycle

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
    if(DISPLAY_DETAILS .and. compute_adjoint) then
        print*
        print*,'irec=',irec
        print*,'window -- ',ntstart,ntend
    endif

    do itype=1,ntype
    ! initialization
    misfit_value=0.0
    adj(:)=0.d0

    measurement_type=trim(measurement_types(itype))
    measurement_attenuation_type=trim(measurement_attenuation)
    if(VISCOELASTIC) then
        print*
        print*,'Wenyong Pan is testing viscoelastic not'
        print*,'misfit type for attenuation: ',measurement_attenuation_type
    endif
    call misfit_adj_AD(measurement_type,d,s,NSTEP,&
        VISCOELASTIC,measurement_attenuation_type,& 
        JOINT_MISFIT,misfit_lambda,&
        deltat,f0,ntstart,ntend,&
        window_type,compute_adjoint, &
        misfit_value,adj) 

    ! sum over itype of misfit and adjoint
    misfit_trace = misfit_trace + misfit_value**2
    if(DISPLAY_DETAILS .and. compute_adjoint) then
        print*,'misfit_',trim(measurement_type),'_AD=',misfit_value
        print*,'squared misfit_',trim(measurement_type),'_AD=',misfit_value**2
    endif

    if(compute_adjoint) adj_trace = adj_trace + adj ! FWI
!    if(compute_adjoint) adj_trace = adj_trace + d ! RTM
    enddo !itype

    if (ntype>=1) then 
        ! window sum 
        misfit_AD = misfit_AD + misfit_trace / ntype
        if(compute_adjoint) then 
            !  call process_adj(adj_trace,ntstart,ntend,dis_sr)
            if(DISPLAY_DETAILS) then
                print*, 'Min/Max of adj :',minval(adj_trace(:)),maxval(adj_trace(:))
            endif

            seism_adj_AD(:,irec) = seism_adj_AD(:,irec) + adj_trace(:) / ntype
        endif
        num_AD = num_AD +1
    endif

    enddo ! irec

end subroutine Absolute_diff

subroutine Yunsong_diff()
    use seismo_parameters; use corr_focus_mod
    integer :: irec, itype, ntype
    real(kind=CUSTOM_REAL) :: seism(NSTEP,nrec_proc)
    real(kind=CUSTOM_REAL) :: d(NSTEP),s(NSTEP),adj_trace(NSTEP),adj(NSTEP)
    real(kind=CUSTOM_REAL) :: wtr
    real(kind=CUSTOM_REAL) :: misfit_value, misfit_trace
    integer :: ntstart,ntend,nlen
    integer :: ishift
    real(kind=CUSTOM_REAL) :: dlnA,cc_max,Terr_max,f0_temp,deltat_temp
    character(len=MAX_STRING_LEN) :: measurement_types(MAX_MISFIT_TYPE)
    character(len=2) :: measurement_type
    character(len=3) :: measurement_attenuation_type
    character, parameter :: delimiter='+'

    ! initialization
    d = 0.0_CUSTOM_REAL
    s = 0.0_CUSTOM_REAL

    call split_string(measurement_list,delimiter,measurement_types,ntype)

    do irec=1,nrec_proc

    ! get data 
    d(:)=seism_obs(:,irec)
    s(:)=seism_syn(:,irec)
    !s(:)=0.d0
    !! misfit and adjoint evaluation
    ! initialization
    adj_trace = 0.0_CUSTOM_REAL
    misfit_trace=0.0
    dis_sr=sqrt((st_xval(irec)-x_source)**2 &
        +(st_yval(irec)-y_source)**2 &
        +(st_zval(irec)-z_source)**2)

    ! window info
    ntstart = win_start(irec)
    ntend= win_end(irec)

    !! quality control (window, norm, cc similarity)
    if (ntstart>=ntend ) cycle
    if (norm2(d(ntstart:ntend),1)<SMALL_VAL .or. &
        norm2(s(ntstart:ntend),1)<SMALL_VAL .or. &
        (ntend-ntstart)*deltat<window_len) &
        cycle

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
    if(DISPLAY_DETAILS .and. compute_adjoint) then
        print*
        print*,'irec=',irec
        print*,'window -- ',ntstart,ntend
    endif

    do itype=1,ntype
    ! initialization
    misfit_value=0.0
    adj(:)=0.d0

    measurement_type=trim(measurement_types(itype))
    measurement_attenuation_type=trim(measurement_attenuation)
    call misfit_adj_AD(measurement_type,d,s,NSTEP,&
        VISCOELASTIC,measurement_attenuation_type,&
        deltat,f0,ntstart,ntend,&
        window_type,compute_adjoint, &
        misfit_value,adj)

    ! sum over itype of misfit and adjoint
    misfit_trace = misfit_trace + misfit_value**2
    if(DISPLAY_DETAILS .and. compute_adjoint) then
        print*,'misfit_',trim(measurement_type),'_AD=',misfit_value
        print*,'squared misfit_',trim(measurement_type),'_AD=',misfit_value**2
    endif

    if(compute_adjoint) adj_trace = adj_trace + adj ! FWI
!    if(compute_adjoint) adj_trace = adj_trace + d ! RTM
    enddo !itype

    if (ntype>=1) then
        ! window sum 
        misfit_AD = misfit_AD + misfit_trace / ntype
        if(compute_adjoint) then
            !  call process_adj(adj_trace,ntstart,ntend,dis_sr)
            if(DISPLAY_DETAILS) then
                print*, 'Min/Max of adj :',minval(adj_trace(:)),maxval(adj_trace(:))
            endif

            !seism_adj_AD(:,irec) = seism_adj_AD(:,irec) + adj_trace(:) / ntype
        endif
        num_AD = num_AD +1
    endif

    enddo ! irec

    Terr_max=0.25
    f0_temp=8
    deltat_temp=5e-4
    call corr_focusing_Adj_Src(seism,misfit_value,&
         seism_syn,seism_obs,f0_temp,deltat_temp,Terr_max)

    if(compute_adjoint) seism_adj_AD=seism_adj_AD+seism

end subroutine Yunsong_diff

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Relative_diff(input_dir,data_name)
    use seismo_parameters
    integer :: irec,jrec,itype, ntype
    character(len=MAX_STRING_LEN) :: input_dir
    character(len=MAX_STRING_LEN) :: data_name
    character(len=MAX_FILENAME_LEN) :: filename
    logical :: ex
    real(kind=CUSTOM_REAL) :: seism(NSTEP,nrec_proc)
    real(kind=CUSTOM_REAL) :: d(NSTEP),s(NSTEP),adj(NSTEP),adj_pair(NSTEP)
    real(kind=CUSTOM_REAL) :: d_ref(NSTEP),s_ref(NSTEP),adj_ref(NSTEP),adj_ref_pair(NSTEP)
    real(kind=CUSTOM_REAL) :: dis_sr1, dis_sr2, dis_rr
    real(kind=CUSTOM_REAL) :: cc_max_obs
    real(kind=CUSTOM_REAL) :: misfit_value,misfit_pair
    integer :: ntstart,ntend,nlen
    integer :: ntstart_ref,ntend_ref,nlen_ref

    character(len=MAX_STRING_LEN) :: measurement_types(MAX_MISFIT_TYPE)
    character(len=2) :: measurement_type
    character, parameter :: delimiter='+'

    call split_string(measurement_list,delimiter,measurement_types,ntype)

    ! initialization
    d = 0.0_CUSTOM_REAL
    s = 0.0_CUSTOM_REAL
    d_ref = 0.0_CUSTOM_REAL
    s_ref = 0.0_CUSTOM_REAL
    cc_max_obs = 0.0

    allocate(is_pair(nrec_proc,nrec_proc))
    is_pair=0

    write(filename,'(a)') trim(input_dir)//'DATA_obs/'//trim(data_name)//'.similarity.dat'
    ! if existence, read, otherwise, write 
    inquire (file=trim(filename), exist=ex)
    OPEN (UNIT=IIN, FILE=filename,iostat=ier)
    do while(ier==0)
    read(IIN,*,iostat=ier) irec,jrec,cc_max_obs,is_pair(irec,jrec)
    enddo

    ! loop over master trace
    do irec=1,nrec_proc

    ! get data 
    d(:)=seism_obs(:,irec)
    s(:)=seism_syn(:,irec)

    dis_sr1=sqrt((st_xval(irec)-x_source)**2 &
        +(st_yval(irec)-y_source)**2 &
        +(st_zval(irec)-z_source)**2)

    ! window info
    ntstart = win_start(irec)
    ntend = win_end(irec)
    nlen=ntend-ntstart+1
    if(nlen<1 .or. nlen>NSTEP) cycle

    if(norm2(d(ntstart:ntend),1)<SMALL_VAL .or. &
        norm2(s(ntstart:ntend),1)<SMALL_VAL) cycle ! zero  master trace
    ! loop over reference trace
    do jrec=irec+1,nrec_proc
    ! get data 
    d_ref(:)=seism_obs(:,jrec)
    s_ref(:)=seism_syn(:,jrec)

    ! window info
    ntstart_ref = win_start(jrec)
    ntend_ref= win_end(jrec)
    nlen_ref=ntend_ref-ntstart_ref+1
    if(nlen_ref<1 .or. nlen_ref>NSTEP) cycle

    if(norm2(d_ref(ntstart_ref:ntend_ref),1)<SMALL_VAL .or. &
        norm2(s_ref(ntstart_ref:ntend_ref),1)<SMALL_VAL) cycle  ! zero reference trace

    if(.not. ex) then
        cc_max_obs=0.0
        dis_rr=sqrt((st_xval(jrec)-st_xval(irec))**2 &
            +(st_yval(jrec)-st_yval(irec))**2 &
            +(st_zval(jrec)-st_zval(irec))**2)
        if(dis_rr<=DD_max .and. dis_rr>=DD_min)&
            call CC_similarity(d,d_ref,NSTEP,&
            ntstart,ntend,ntstart_ref,ntend_ref,window_type,cc_max_obs)
        if(cc_max_obs>cc_threshold)  is_pair(irec,jrec) = 1
    endif !! ex

    if(DISPLAY_DETAILS .and. compute_adjoint) then
        print*      
        print*,'pair irec=',irec, 'jrec=',jrec           
        print*,'window -- ',ntstart,ntend,ntstart_ref,ntend_ref
        print*, 'rr distance dis_rr, DD_min/DD_max: ',dis_rr, DD_min, DD_max 
        print*, 'cc similarity -- ', cc_max_obs 
        print*,'is_pair : ',is_pair(irec,jrec)
    endif 

    if(is_pair(irec,jrec)==1) then
        ! initialization
        misfit_pair = 0.0
        adj_pair = 0.0
        adj_ref_pair = 0.0

        dis_sr2=sqrt((st_xval(jrec)-x_source)**2 &
            +(st_yval(jrec)-y_source)**2 &
            +(st_zval(jrec)-z_source)**2)

        ! number of double difference measurements
        num_DD = num_DD+1

        do itype=1,ntype
        ! initialization
        misfit_value = 0.0
        adj = 0.0
        adj_ref = 0.0

        measurement_type=trim(measurement_types(itype))

        call misfit_adj_DD(measurement_type,d,d_ref,s,s_ref,NSTEP,deltat,f0,&
            ntstart,ntend,ntstart_ref,ntend_ref,window_type,compute_adjoint,&
            misfit_value,adj,adj_ref)

        ! sum over itype of misfit and adjoint
        misfit_pair=misfit_pair+misfit_value**2
        if(DISPLAY_DETAILS .and. compute_adjoint) then
            print*,'misfit_',trim(measurement_type),'_DD=',misfit_value
            print*,'squared misfit_',trim(measurement_type),'_DD=',misfit_value**2
        endif

        if(compute_adjoint) then
            adj_pair=adj_pair+adj
            adj_ref_pair=adj_ref_pair+adj_ref
        endif

        enddo ! itype

        if (ntype>=1) then 

            ! sum of misfit over all stations
            misfit_DD=misfit_DD+ misfit_pair/ntype

            if(compute_adjoint) then
                call process_adj(adj_pair,ntstart,ntend,dis_sr1)
                call process_adj(adj_ref_pair,ntstart_ref,ntend_ref,dis_sr2)

                if(DISPLAY_DETAILS) then
                    print*, 'Min/Max of adj :',minval(adj_pair(:)),maxval(adj_pair(:))
                    print*, 'Min/Max of adj_ref:',minval(adj_ref_pair(:)),maxval(adj_ref_pair(:))
                endif

                ! sum of adjoint source over pair
                seism_adj_DD(:,irec) = seism_adj_DD(:,irec) +  adj_pair(:)/ntype
                seism_adj_DD(:,jrec) = seism_adj_DD(:,jrec) + adj_ref_pair(:)/ntype
            endif
        endif

        !! save waveform similarity
        if(.not. ex) write(IIN,'(2I5,1e15.5,I5)') irec,jrec,cc_max_obs,is_pair(irec,jrec)

    endif  ! is_pair

    enddo  ! jrec 
    enddo ! irec 

    close(IIN) ! close similarity file

    deallocate(is_pair)

end subroutine Relative_diff
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine finalize(directory,data_name)
    use seismo_parameters
    implicit none
    integer :: irec
    character(len=MAX_STRING_LEN) :: data_name
    character(len=MAX_FILENAME_LEN) :: char_i,filen,filename
    character(len=MAX_STRING_LEN) :: directory
    real(kind=CUSTOM_REAL) :: dis_sr

    !! write SU adj
    if(compute_adjoint) then

        !! adjoint file format
        select case(solver)
        case('specfem2D')
            if(myrank==0 .and. data_name == 'x') filen='Ux_file_single.su.adj'
            if(myrank==0 .and. data_name == 'y') filen='Uy_file_single.su.adj'
            if(myrank==0 .and. data_name == 'z') filen='Uz_file_single.su.adj'
            if(myrank==0 .and. data_name == 'p') filen='Up_file_single.su.adj'

        case('specfem3D')
            write(char_i, '(I5)') myrank        ! convert integer to char
            if(data_name == 'x') write(filen,'(a)') trim(adjustl(char_i))//'_dx_SU.adj'
            if(data_name == 'y') write(filen,'(a)') trim(adjustl(char_i))//'_dy_SU.adj'
            if(data_name == 'z') write(filen,'(a)') trim(adjustl(char_i))//'_dz_SU.adj'
            if(data_name == 'p') write(filen,'(a)') trim(adjustl(char_i))//'_dp_SU.adj'

        case default
            print*,'Currently only work for specfem2D and specfem3D solver'
            stop
        end select

        !! write adjoint source
        write(filename,'(a)') &
            trim(directory)//'/SEM/'//trim(filen)
!!              trim(directory)//'/SEM/'//trim(filen)
        if(DISPLAY_DETAILS) then
            print*
            print*,'SAVE seism_adj -- ',trim(filename)
            if(DISPLAY_DETAILS) print*,'Min / Max of final seism_adj : ',&
                minval(seism_adj(:,:)),maxval(seism_adj(:,:))
        endif

        call writeSU(filename,seism_adj)
    endif

    deallocate(seism_obs)
    deallocate(seism_syn)
    deallocate(seism_adj)
    deallocate(seism_adj_AD)
    deallocate(seism_adj_DD)
    deallocate(st_xval)
    deallocate(st_yval)
    deallocate(st_zval)
    deallocate(win_start)
    deallocate(win_end)

end subroutine finalize
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine readSU(filename,seism)
    use seismo_parameters
    integer :: ier,irec
    real(kind=CUSTOM_REAL) :: seism(NSTEP,nrec_proc)
    character(len=MAX_FILENAME_LEN) :: filename

    seism = 0.0_CUSTOM_REAL
    ! open(IIN,file=trim(filename),access='direct',recl=240+4*NSTEP,iostat = ier)
    open(IIN,file=trim(filename),status='old',form='unformatted',&
        access='direct',recl=240+4*NSTEP,iostat = ier)

    if (ier /= 0) then
        print*, 'Error: could not open data file: ',trim(filename)
        stop
    endif

    do irec = 1, nrec_proc
    read(IIN,rec=irec,iostat=ier) r4head, seism(:,irec)
    ! header info
    z_source=r4head(13) ! Source depth below surface (sdepth) 
    x_source=r4head(19) ! Source x coord (sx)
    y_source=r4head(20) ! Source y coord  (sy)
    st_zval(irec)=r4head(11) ! Receiver group elevation (gelev)
    st_xval(irec)=r4head(21) ! Receiver x coord (gx)
    st_yval(irec)=r4head(22) ! Receiver y coord (gy)
    ! header2=int(r4head(29), kind=2)
    if (DISPLAY_DETAILS .and. irec==1) print *, 'xs,ys,zs',&
        x_source,y_source,z_source
    if (DISPLAY_DETAILS .and. irec==1) print *, 'xr,yr,zr',&
        st_xval(irec),st_yval(irec),st_zval(irec)

    enddo

    close(IIN)

end subroutine readSU
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine writeSU(filename,seism)
    use seismo_parameters
    integer :: ier,irec,itime
    real(kind=CUSTOM_REAL) :: seism(NSTEP,nrec_proc)
    character(len=MAX_FILENAME_LEN) :: filename
    integer :: deltat_int2

    open(IOUT,file=trim(filename),status='unknown',access='direct',recl=4,iostat=ier)
    if (ier /= 0) then
        print*, 'Error: could not open data file: ',trim(filename)
        stop
    endif

    do irec = 1, nrec_proc
    ! write header
    write(IOUT,rec=(irec-1)*60+(irec-1)*NSTEP+1)  irec !receiver ID
    write(IOUT,rec=(irec-1)*60+(irec-1)*NSTEP+10) NINT(st_xval(irec)-x_source)  ! offset
    write(IOUT,rec=(irec-1)*60+(irec-1)*NSTEP+19) NINT(x_source)                ! source location xs
    write(IOUT,rec=(irec-1)*60+(irec-1)*NSTEP+20) NINT(y_source)                ! source location ys
    write(IOUT,rec=(irec-1)*60+(irec-1)*NSTEP+21) NINT(st_xval(irec))           ! receiver location xr
    write(IOUT,rec=(irec-1)*60+(irec-1)*NSTEP+22) NINT(st_yval(irec))           ! receiver location zr
    if (nrec_proc>1) write(IOUT,rec=(irec-1)*60+(irec-1)*NSTEP+48) SNGL(st_xval(2)-st_xval(1)) ! receiver interval
    header2(1)=0  ! dummy
    header2(2)=int(NSTEP, kind=2)
    write(IOUT,rec=(irec-1)*60+(irec-1)*NSTEP+29) header2
    header2(1)=NINT(deltat*1.0d6, kind=2) ! deltat (unit: 10^{-6} s)
    header2(2)=0  ! dummy
    write(IOUT,rec=(irec-1)*60+(irec-1)*NSTEP+30) header2

    ! the "60" in the following corresponds to 240 bytes header (note the
    ! reclength is 4 bytes)
    do itime = 1, NSTEP
    write(IOUT,rec=irec*60+(irec-1)*NSTEP+itime) sngl(seism_adj(itime,irec))
    enddo
    enddo ! irec

    close(IOUT)

end subroutine writeSU

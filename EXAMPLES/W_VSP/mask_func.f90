program mask_func
    ! To mask source

#ifdef USE_MPI
    use mpi
#endif

    use seismo_parameters
    implicit none

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
    integer, parameter :: NARGS = 4
    character(len=MAX_STRING_LEN) :: arg(NARGS)
    character(len=MAX_STRING_LEN) :: input_dir,output_dir
    INTEGER :: i, ier
    real t1,t2

#ifdef USE_MPI
    call MPI_INIT(ier)
    call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ier)
    call MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ier)
#else
    nproc = 1
    myrank = 0
#endif

    call cpu_time(t1)

    ! parse command line arguments
    if (command_argument_count() /= NARGS) then
        if (myrank == 0) then
            print *, 'USAGE:  mpirun -np NPROC bin/sum_kernel.exe ...'
            stop ' Please check command line arguments'
        endif
    endif

#ifdef USE_MPI
    call MPI_BARRIER(MPI_COMM_WORLD,ier)
#endif

    do i = 1, NARGS
    call get_command_argument(i,arg(i), status=ier)
    enddo

    read(arg(1),*) x_source
    read(arg(2),*) z_source
    input_dir=arg(3)
    output_dir=arg(4)

    !! initialization  -- get number of spectral elements
    call initialize(input_dir)

    !! generate mask file
    if(MASK_SOURCE) call add_source_mask()
    if (MASK_STATION) call add_station_mask(input_dir)
    if (MASK_MODEL) call add_model_mask()
    !! save mask file
    call finalize(output_dir)

    call cpu_time(t2)

#ifdef USE_MPI
    ! stop all the processes and exit
    call MPI_FINALIZE(ier)
#endif

end program mask_func
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine initialize(directory)
    use seismo_parameters
    implicit none
    integer :: ier
    character(len=MAX_FILENAME_LEN) :: filename
    character(len=MAX_STRING_LEN) :: directory
    logical :: ex=.false.

    write(filename,'(a,i6.6,a)') trim(directory)//'/proc',myrank,'_'//trim(IBOOL_NAME)
    open(IIN,file=trim(filename),status='old',action='read',form='unformatted',iostat=ier)
    if (ier /= 0) then
        print *,'Error: could not open database file: ',trim(filename)
        stop 'Error opening _NSPEC_IBOOL file'
    endif
    read(IIN) nspec
    close(IIN)

    allocate(mask(NGLLX,NGLLY,NGLLZ,NSPEC))
    mask = 1.0
    allocate(xstore(NGLLX,NGLLY,NGLLZ,NSPEC))
    allocate(ystore(NGLLX,NGLLY,NGLLZ,NSPEC))
    allocate(zstore(NGLLX,NGLLY,NGLLZ,NSPEC))
    xstore =  0.0
    ystore = 0.0 
    zstore = 0.0

    write(filename, '(a,i6.6,a)') trim(directory)//'/proc',myrank,'_x.bin'
    ! gets the coordinate x of the points located in my slice
    open(unit=IIN,file=trim(filename),status='old',action='read',form='unformatted',iostat=ier)
    if (ier /= 0) then
        print *,'Error: could not open database file: ',trim(filename)
        stop 'Error reading neighbors external mesh file'
    endif
    ! global point arrays
    read(IIN) xstore
    close(IIN)

    !exist
    write(filename, '(a,i6.6,a)') trim(directory)//'/proc',myrank,'_y.bin'
    inquire (file=trim(filename), exist=ex)
    if(ex) then
        ! gets the coordinate y of the points located in my slice
        open(unit=IIN,file=trim(filename),status='old',action='read',form='unformatted',iostat=ier)
        if (ier /= 0) then
            print *,'Error: could not open database file: ',trim(filename)
            stop 'Error reading neighbors external mesh file'
        endif
        ! global point arrays
        read(IIN) ystore
        close(IIN)
    endif

    write(filename, '(a,i6.6,a)') trim(directory)//'/proc',myrank,'_z.bin'
    ! gets the coordinate z of the points located in my slice
    open(unit=IIN,file=trim(filename),status='old',action='read',form='unformatted',iostat=ier)
    if (ier /= 0) then
        print *,'Error: could not open database file: ',trim(filename)
        stop 'Error reading neighbors external mesh file'
    endif
    ! global point arrays
    read(IIN) zstore
    close(IIN)

end subroutine initialize
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine add_source_mask()
    use seismo_parameters
    implicit none
    integer :: ix,iy,iz,ispec
    real :: dis

    if(myrank==0) then
        print* 
        print*, 'mask source at:', x_source, y_source, z_source
        print* 
    endif

    do ix = 1,NGLLX
    do iy = 1,NGLLY
    do iz = 1,NGLLZ
    do ispec = 1,NSPEC
    dis = sqrt((xstore(ix,iy,iz,ispec)-x_source)**2+&
        (ystore(ix,iy,iz,ispec)-y_source)**2+ &
        (zstore(ix,iy,iz,ispec)-z_source)**2)
    if(dis<=source_radius) then 
        ! gaussian
        if (MASK_DAMP) then
          mask(ix,iy,iz,ispec)=mask(ix,iy,iz,ispec) * &
              (1.0 - exp(-4.0*(dis/source_radius)**2))**4
        endif
        if (.not. MASK_DAMP) then
          ! mute 
          mask(ix,iy,iz,ispec)=0.0
        endif
    endif
    enddo
    enddo
    enddo
    enddo

end subroutine add_source_mask

subroutine add_model_mask()
    use seismo_parameters
    implicit none
    integer :: ix,iy,iz,ispec
    real :: dis

    if(myrank==0) then
        print*
        print*, 'mask source at:', x_source, y_source, z_source
        print*
    endif

    do ix = 1,NGLLX
    do iy = 1,NGLLY
    do iz = 1,NGLLZ
    do ispec = 1,NSPEC
    if(zstore(ix,iy,iz,ispec)<=mask_z) then
        ! gaussian
        !mask(ix,iy,iz,ispec)=mask(ix,iy,iz,ispec) * &
        !    (1.0 - exp(-4.0*(dis/source_radius)**2))**4
        ! mute 
        mask(ix,iy,iz,ispec)=0.0
    endif
    if(zstore(ix,iy,iz,ispec)>mask_zend) then
        ! gaussian
        !mask(ix,iy,iz,ispec)=mask(ix,iy,iz,ispec) * &
        !    (1.0 - exp(-4.0*(dis/source_radius)**2))**4
        ! mute 
        mask(ix,iy,iz,ispec)=0.0
    endif
    enddo
    enddo
    enddo
    enddo
end subroutine add_model_mask
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine add_station_mask(directory)
    use seismo_parameters
    implicit none
    character(len=MAX_STRING_LEN) :: directory
    integer :: ier
    character(len=MAX_FILENAME_LEN) :: filename
    integer :: ix,iy,iz,ispec,irec
    real :: dis
    CHARACTER (LEN=MAX_STRING_LEN) :: station_name, network_name
    real(kind=CUSTOM_REAL) :: stele,stbur

    allocate(st_xval(NREC))
    allocate(st_yval(NREC))
    allocate(st_zval(NREC))
    st_xval = 0.0 
    st_yval = 0.0 
    st_zval = 0.0

    if(myrank==0) then
        print*
        print*, 'mask stations'
        print*
    endif

    !! read STATION file 
    write(filename,'(a)') trim(directory)//'/STATIONS'
    OPEN (UNIT=IIN,FILE=trim(filename),STATUS='OLD',action='read',iostat=ier)
    if(ier /=0) then
        print*,'Error opening file. File load station file'
        stop
    else
        do irec=1,NREC ! trace loop
        read(IIN,*) station_name,network_name,st_xval(irec),st_zval(irec),stele,stbur
        enddo
    end if
    close(IIN)

    do ix = 1,NGLLX
    do iy = 1,NGLLY
    do iz = 1,NGLLZ
    do ispec = 1,NSPEC
    do irec=1,NREC
    dis = sqrt((xstore(ix,iy,iz,ispec)-st_xval(irec))**2+&
        (ystore(ix,iy,iz,ispec)-st_yval(irec))**2+ &
        (zstore(ix,iy,iz,ispec)-st_zval(irec))**2)
    if(dis<=station_radius) then
        ! gaussian damping
        mask(ix,iy,iz,ispec)=mask(ix,iy,iz,ispec) * & 
            (1.0 - exp(-4.0*(dis/station_radius)**2))**4
        ! zero mute
        !  mask(ix,iy,iz,ispec)= 0.0
    endif
    enddo
    enddo
    enddo
    enddo
    enddo

    deallocate(st_xval)
    deallocate(st_yval)
    deallocate(st_zval)

end subroutine add_station_mask
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine finalize(directory)
    use seismo_parameters
    implicit none
    integer :: ier
    character(len=MAX_FILENAME_LEN) :: filename
    character(len=MAX_STRING_LEN) :: directory

    write(filename,'(a,i6.6,a)') trim(directory)//'/proc',myrank,'_mask.bin'
    open(IOUT,file=trim(filename),status='unknown',action='write',form='unformatted',iostat=ier)
    if (ier /= 0) then
        print *,'Error: could not open source mask file: ',trim(filename)
        stop 'Error opening source mask file'
    endif
    write(IOUT) mask
    close(IOUT)

    deallocate(mask)
    deallocate(xstore)
    deallocate(ystore)
    deallocate(zstore)

end subroutine finalize
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

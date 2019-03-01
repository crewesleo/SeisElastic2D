program model_update
    ! To update model along search direction
    ! yanhuay@princeton.edu

    use seismo_parameters
    implicit none

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
    integer, parameter :: NARGS = 4
    INTEGER :: itime, ier, isrc,i,j
    character(len=MAX_STRING_LEN) :: model_names(MAX_KERNEL_NUM)
    character(len=MAX_STRING_LEN) :: model_names_comma_delimited
    character(len=MAX_STRING_LEN) :: kernel_names(MAX_KERNEL_NUM)
    character(len=MAX_STRING_LEN) :: kernel_names_comma_delimited
    character(len=MAX_STRING_LEN) :: arg(NARGS)
    character(len=MAX_STRING_LEN) :: directory
    real t1,t2
    character, parameter :: delimiter=','

    call cpu_time(t1)

    ! parse command line arguments
    if (command_argument_count() /= NARGS) then
        print *, 'USAGE:  mpirun -np NPROC bin/model_update.exe nproc directory MODEL_NAME, step_length'
        stop ' Please check command line arguments'
    endif

    do i = 1, NARGS
    call get_command_argument(i,arg(i), status=ier)
    enddo

    read(arg(1),*) nproc
    directory=arg(2)
    model_names_comma_delimited = arg(3)
    read(arg(4),*) step_length
!    read(arg(5),*) scaling_factor
!    kernel_names_comma_delimited = arg(5)
    if (myrank == 0) write(*,'(a,f15.2,a)') 'try step length -- ',step_length*100,'%'

    call split_string(model_names_comma_delimited,delimiter,model_names,nmod)
!    call split_string(kernel_names_comma_delimited,delimiter,kernel_names,nmod)
    !! initialization  -- get number of spectral elements
    call initialize(directory)

    !! update model
    !! check parameterization
    !if (model_names(1)=='vp') then
    call update(directory)
    !    print *, 'Velocity-density Parameterization !!!'
    !endif
!    if (kernel_names(1)=='kappa_kernel') then
!        call update_modulus(directory)
!        print *, 'Modulus-density Parameterization !!!'
!    endif

    !! save updated model  
    call finalize(directory,adjustl(model_names(1:nmod)))

end program model_update
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine initialize(directory)
    use seismo_parameters
    implicit none
    integer :: ier,filesize
    character(len=MAX_FILENAME_LEN) :: filename
    character(len=MAX_STRING_LEN) :: directory
    character(len=MAX_STRING_LEN) :: model_name

    ! slice database file
    allocate(nspec_proc(nproc))
    nspec_proc=0
    do myrank=0,nproc-1

    ! nspec
    write(filename,'(a,i6.6,a)') &
        trim(directory)//'/misfit_kernel/proc',myrank,'_'//trim(IBOOL_NAME) 
    write(filename,'(a,i6.6,a)') &
        trim(directory)//'/misfit_kernel/proc',myrank,'_'//trim(IBOOL_NAME)
    open(IIN,file=trim(filename),status='old',action='read',form='unformatted',iostat=ier)            
    if (ier /= 0) then                          
        print *,'Error: could not open database file:',trim(filename)
        stop 'Error opening NSPEC_IBOOL file'       
    endif                                                
    read(IIN) nspec_proc(myrank+1)
    close(IIN)                                                          

    enddo

    nspec=sum(nspec_proc)
    if(DISPLAY_DETAILS) print*,'NGLLX*NGLLY*NGLLZ*NSPEC*nmod:',NGLLX,NGLLY,NGLLZ,NSPEC,nmod

    allocate(m_new(NGLLX*NGLLY*NGLLZ*NSPEC*nmod))
    m_new = 0.0_CUSTOM_REAL
    allocate(p_new(NGLLX*NGLLY*NGLLZ*NSPEC*nmod))
    p_new = 0.0_CUSTOM_REAL
    allocate(m_try(NGLLX*NGLLY*NGLLZ*NSPEC*nmod))
    m_try = 0.0_CUSTOM_REAL

end subroutine initialize
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine update(directory)
    use seismo_parameters
    implicit none
    integer :: ier,imod
    character(len=MAX_FILENAME_LEN) :: filename
    character(len=MAX_STRING_LEN) :: directory
    real(kind=CUSTOM_REAL), dimension(:,:,:,:,:), allocatable :: temp_store
    !! LOAD p_new
    write(filename,'(a)') &
        trim(directory)//'/optimizer/p_new.bin'
    print*,'LOAD p_new -- ', trim(filename)
    open(unit=IIN,file=trim(filename),status='old',action='read',form='unformatted',iostat=ier)
    if (ier /= 0) then
        print*, 'Error: could not open model file: ',trim(filename)
        stop 'Error reading neighbors external mesh file'
    endif
    read(IIN) p_new
    close(IIN)

    !! LOAD m_new
    write(filename,'(a)') &
        trim(directory)//'/optimizer/m_new.bin'
    print*,'LOAD m_new -- ', trim(filename)
    open(unit=IIN,file=trim(filename),status='old',action='read',form='unformatted',iostat=ier)
    if (ier /= 0) then
        print*, 'Error: could not open model file: ',trim(filename)
        stop 'Error reading neighbors external mesh file'
    endif
    read(IIN) m_new
    close(IIN)
    if(DISPLAY_DETAILS .and. myrank==0) print *,'Min / Max m_new = ', &
        minval(m_new(:)),maxval(m_new(:))
 !!   step_length=100
    !! update 
!    p_new=p_new*100
    m_try = m_new * (1+step_length*p_new) 
    if(DISPLAY_DETAILS .and. myrank==0) print *,'Min / Max m_try = ', &
        minval(m_try(:)),maxval(m_try(:))

end subroutine update
!! HTI elastic constant parameterization
subroutine update_hti_ec(directory)
    use seismo_parameters
    implicit none
    integer :: ier,imod
    character(len=MAX_FILENAME_LEN) :: filename
    character(len=MAX_STRING_LEN) :: directory
    real(kind=CUSTOM_REAL), dimension(:,:,:,:,:), allocatable :: temp_store
    !! LOAD p_new
    write(filename,'(a)') &
        trim(directory)//'/optimizer/p_new.bin'
    print*,'LOAD p_new -- ', trim(filename)
    open(unit=IIN,file=trim(filename),status='old',action='read',form='unformatted',iostat=ier)
    if (ier /= 0) then
        print*, 'Error: could not open model file: ',trim(filename)
        stop 'Error reading neighbors external mesh file'
    endif
    read(IIN) p_new
    close(IIN)

    !! LOAD m_new
    write(filename,'(a)') &
        trim(directory)//'/optimizer/m_new.bin'
    print*,'LOAD m_new -- ', trim(filename)
    open(unit=IIN,file=trim(filename),status='old',action='read',form='unformatted',iostat=ier)
    if (ier /= 0) then
        print*, 'Error: could not open model file: ',trim(filename)
        stop 'Error reading neighbors external mesh file'
    endif
    read(IIN) m_new
    close(IIN)
    if(DISPLAY_DETAILS .and. myrank==0) print *,'Min / Max m_new = ', &
        minval(m_new(:)),maxval(m_new(:))
 !!   step_length=100
    !! update 

    m_try = m_new * (1+step_length*p_new)
    if(DISPLAY_DETAILS .and. myrank==0) print *,'Min / Max m_try = ', &
        minval(m_new(:)),maxval(m_new(:))

end subroutine update_hti_ec

!! Modulus-density update
subroutine update_modulus(directory)
    use seismo_parameters
    implicit none
    integer :: ier,imod
    character(len=MAX_FILENAME_LEN) :: filename
    character(len=MAX_STRING_LEN) :: directory
    real(kind=CUSTOM_REAL), dimension(:,:,:,:,:), allocatable :: temp_store
    real(kind=CUSTOM_REAL), dimension(:), allocatable :: temp_store_modulus
    real(kind=CUSTOM_REAL), dimension(:,:,:,:,:), allocatable :: m_new_modulus
    real(kind=CUSTOM_REAL), dimension(:), allocatable :: m_try_modulus
    real(kind=CUSTOM_REAL), dimension(:,:,:,:,:), allocatable :: m_new_velocity
    !! LOAD p_new
    write(filename,'(a)') &
        trim(directory)//'/optimizer/p_new.bin'
    print*,'LOAD p_new -- ', trim(filename)
    open(unit=IIN,file=trim(filename),status='old',action='read',form='unformatted',iostat=ier)
    if (ier /= 0) then
        print*, 'Error: could not open model file: ',trim(filename)
        stop 'Error reading neighbors external mesh file'
    endif
    read(IIN) p_new
    close(IIN)

    !! LOAD m_new
    write(filename,'(a)') &
        trim(directory)//'/optimizer/m_new.bin'
    print*,'LOAD m_new -- ', trim(filename)
    open(unit=IIN,file=trim(filename),status='old',action='read',form='unformatted',iostat=ier)
    if (ier /= 0) then
        print*, 'Error: could not open model file: ',trim(filename)
        stop 'Error reading neighbors external mesh file'
    endif
    read(IIN) m_new
    close(IIN)
    if(DISPLAY_DETAILS .and. myrank==0) print *,'Min / Max m_new = ', &
        minval(m_new(:)),maxval(m_new(:))
    !!   step_length=100
    !! update
    allocate(temp_store(NGLLX,NGLLY,NGLLZ,NSPEC,nmod))
    temp_store=0.0_CUSTOM_REAL
    allocate(temp_store_modulus(NGLLX*NGLLY*NGLLZ*NSPEC*nmod))
    temp_store_modulus=0.0_CUSTOM_REAL
    allocate(m_new_modulus(NGLLX,NGLLY,NGLLZ,NSPEC,nmod))
    m_new_modulus=0.0_CUSTOM_REAL
    allocate(m_try_modulus(NGLLX*NGLLY*NGLLZ*NSPEC*nmod))
    m_try_modulus=0.0_CUSTOM_REAL
    allocate(m_new_velocity(NGLLX,NGLLY,NGLLZ,NSPEC,nmod))
    m_new_velocity=0.0_CUSTOM_REAL

    ! p_new_modulus=reshape(p_new,shape(p_new_modulus))
    temp_store=reshape(m_new,shape(temp_store))
    !! Transform velocity models into modulus models
    ! kappa
    m_new_modulus(:,:,:,:,1) = temp_store(:,:,:,:,3)*temp_store(:,:,:,:,1)**2-4*temp_store(:,:,:,:,3)*temp_store(:,:,:,:,2)**2/3
    ! mu
    m_new_modulus(:,:,:,:,2) = temp_store(:,:,:,:,3)*temp_store(:,:,:,:,2)**2
    ! rho
    m_new_modulus(:,:,:,:,3) = temp_store(:,:,:,:,3)
    !! update modulus density parameters
    temp_store_modulus = reshape(m_new_modulus,shape(temp_store_modulus))
    m_try_modulus = temp_store_modulus * (1+step_length*p_new)
    !! transform updated budulus models into velocity models
    temp_store=reshape(m_try_modulus,shape(temp_store))    
    m_new_velocity(:,:,:,:,2) = SQRT(temp_store(:,:,:,:,2)/temp_store(:,:,:,:,3))
    m_new_velocity(:,:,:,:,1) = SQRT((temp_store(:,:,:,:,1)+4*temp_store(:,:,:,:,2)/3)/temp_store(:,:,:,:,3))
    m_new_velocity(:,:,:,:,3) = temp_store(:,:,:,:,3)
    m_try = reshape(m_new_velocity,shape(m_try))
!    m_try = m_new * (1+step_length*p_new)
    if(DISPLAY_DETAILS .and. myrank==0) print *,'Min / Max m_try = ', &
        minval(m_try(:)),maxval(m_try(:))

end subroutine update_modulus

!! Lame parameterization: Lambda, Mup, Rhopp
subroutine update_lame(directory)
    use seismo_parameters
    implicit none
    integer :: ier,imod
    character(len=MAX_FILENAME_LEN) :: filename
    character(len=MAX_STRING_LEN) :: directory
    real(kind=CUSTOM_REAL), dimension(:,:,:,:,:), allocatable :: temp_store
    real(kind=CUSTOM_REAL), dimension(:), allocatable :: temp_store_lame
    real(kind=CUSTOM_REAL), dimension(:,:,:,:,:), allocatable :: m_new_lame
    real(kind=CUSTOM_REAL), dimension(:), allocatable :: m_try_lame
    real(kind=CUSTOM_REAL), dimension(:,:,:,:,:), allocatable :: m_new_velocity
    real(kind=CUSTOM_REAL), dimension(:,:,:,:,:), allocatable :: p_new_lame
    !! LOAD p_new
    write(filename,'(a)') &
        trim(directory)//'/optimizer/p_new.bin'
    print*,'LOAD p_new -- ', trim(filename)
    open(unit=IIN,file=trim(filename),status='old',action='read',form='unformatted',iostat=ier)
    if (ier /= 0) then
        print*, 'Error: could not open model file: ',trim(filename)
        stop 'Error reading neighbors external mesh file'
    endif
    read(IIN) p_new
    close(IIN)

    !! LOAD m_new
    write(filename,'(a)') &
        trim(directory)//'/optimizer/m_new.bin'
    print*,'LOAD m_new -- ', trim(filename)
    open(unit=IIN,file=trim(filename),status='old',action='read',form='unformatted',iostat=ier)
    if (ier /= 0) then
        print*, 'Error: could not open model file: ',trim(filename)
        stop 'Error reading neighbors external mesh file'
    endif
    read(IIN) m_new
    close(IIN)
    if(DISPLAY_DETAILS .and. myrank==0) print *,'Min / Max m_new = ', &
        minval(m_new(:)),maxval(m_new(:))
    !!   step_length=100
    !! update
    allocate(temp_store(NGLLX,NGLLY,NGLLZ,NSPEC,nmod))
    temp_store=0.0_CUSTOM_REAL
    allocate(temp_store_lame(NGLLX*NGLLY*NGLLZ*NSPEC*nmod))
    temp_store_lame=0.0_CUSTOM_REAL
    allocate(m_new_lame(NGLLX,NGLLY,NGLLZ,NSPEC,nmod))
    m_new_lame=0.0_CUSTOM_REAL
    allocate(m_try_lame(NGLLX*NGLLY*NGLLZ*NSPEC*nmod))
    m_try_lame=0.0_CUSTOM_REAL
    allocate(m_new_velocity(NGLLX,NGLLY,NGLLZ,NSPEC,nmod))
    m_new_velocity=0.0_CUSTOM_REAL
    allocate(p_new_lame(NGLLX,NGLLY,NGLLZ,NSPEC,nmod))
    p_new_lame=0.0_CUSTOM_REAL

    temp_store=reshape(m_new,shape(temp_store))
    !! Transform velocity models into lame models
    ! lambda=rho*vp^2-2*rho*vs^2
    m_new_lame(:,:,:,:,1) = temp_store(:,:,:,:,3)*temp_store(:,:,:,:,1)**2-2*temp_store(:,:,:,:,3)*temp_store(:,:,:,:,2)**2
    ! mup
    m_new_lame(:,:,:,:,2) = temp_store(:,:,:,:,3)*temp_store(:,:,:,:,2)**2
    ! rhoppp
    m_new_lame(:,:,:,:,3) = temp_store(:,:,:,:,3)
    !! update lame density parameters
    !! transform velocity sensitivity kernels into lame sensitivity kernels
    !! K_lambda=K_alpha(lambda/(2*(lambda+2*mup))), K_mup=K_alpha*mup/(lambda+2*mup)+K_beta/2, K_rhoppp=-K_alpha/2-k_beta/2+K_rhop
    temp_store=reshape(p_new,shape(temp_store))
    ! K_lambda
    p_new_lame(:,:,:,:,1) = temp_store(:,:,:,:,1) * m_new_lame(:,:,:,:,1)/(2*(m_new_lame(:,:,:,:,1)+2*m_new_lame(:,:,:,:,2)))
    ! K_mup
    p_new_lame(:,:,:,:,2) = temp_store(:,:,:,:,1) * m_new_lame(:,:,:,:,2)/((m_new_lame(:,:,:,:,1)+2*m_new_lame(:,:,:,:,2))) + temp_store(:,:,:,:,2)/2
    ! K_rhoppp
    p_new_lame(:,:,:,:,3) = - temp_store(:,:,:,:,1)/2 - temp_store(:,:,:,:,2)/2 + temp_store(:,:,:,:,3) 

    p_new=reshape(p_new_lame,shape(p_new))
    temp_store_lame = reshape(m_new_lame,shape(temp_store_lame))
    m_try_lame = temp_store_lame * (1+step_length*p_new)
    !! transform updated lame models into velocity models
    temp_store=reshape(m_try_lame,shape(temp_store))
    m_new_velocity(:,:,:,:,2) = SQRT(temp_store(:,:,:,:,2)/temp_store(:,:,:,:,3))
    m_new_velocity(:,:,:,:,1) = SQRT((temp_store(:,:,:,:,1)+2*temp_store(:,:,:,:,2))/temp_store(:,:,:,:,3))
    m_new_velocity(:,:,:,:,3) = temp_store(:,:,:,:,3)
    m_try = reshape(m_new_velocity,shape(m_try))
!    m_try = m_new * (1+step_length*p_new)
    if(DISPLAY_DETAILS .and. myrank==0) print *,'Min / Max m_try = ', &
        minval(m_try(:)),maxval(m_try(:))

end subroutine update_lame
! Impedance-I parameterization
subroutine update_impedance_I(directory)
    use seismo_parameters
    implicit none
    integer :: ier,imod
    character(len=MAX_FILENAME_LEN) :: filename
    character(len=MAX_STRING_LEN) :: directory
    real(kind=CUSTOM_REAL), dimension(:,:,:,:,:), allocatable :: temp_store
    real(kind=CUSTOM_REAL), dimension(:), allocatable :: temp_store_impedance
    real(kind=CUSTOM_REAL), dimension(:,:,:,:,:), allocatable :: m_new_impedance
    real(kind=CUSTOM_REAL), dimension(:), allocatable :: m_try_impedance
    real(kind=CUSTOM_REAL), dimension(:,:,:,:,:), allocatable :: m_new_velocity
    !! LOAD p_new
    write(filename,'(a)') &
        trim(directory)//'/optimizer/p_new.bin'
    print*,'LOAD p_new -- ', trim(filename)
    open(unit=IIN,file=trim(filename),status='old',action='read',form='unformatted',iostat=ier)
    if (ier /= 0) then
        print*, 'Error: could not open model file: ',trim(filename)
        stop 'Error reading neighbors external mesh file'
    endif
    read(IIN) p_new
    close(IIN)

    !! LOAD m_new
    write(filename,'(a)') &
        trim(directory)//'/optimizer/m_new.bin'
    print*,'LOAD m_new -- ', trim(filename)
    open(unit=IIN,file=trim(filename),status='old',action='read',form='unformatted',iostat=ier)
    if (ier /= 0) then
        print*, 'Error: could not open model file: ',trim(filename)
        stop 'Error reading neighbors external mesh file'
    endif
    read(IIN) m_new
    close(IIN)
    if(DISPLAY_DETAILS .and. myrank==0) print *,'Min / Max m_new = ', &
        minval(m_new(:)),maxval(m_new(:))
    !!   step_length=100
    !! update
    allocate(temp_store(NGLLX,NGLLY,NGLLZ,NSPEC,nmod))
    temp_store=0.0_CUSTOM_REAL
    allocate(temp_store_impedance(NGLLX*NGLLY*NGLLZ*NSPEC*nmod))
    temp_store_impedance=0.0_CUSTOM_REAL
    allocate(m_new_impedance(NGLLX,NGLLY,NGLLZ,NSPEC,nmod))
    m_new_impedance=0.0_CUSTOM_REAL
    allocate(m_try_impedance(NGLLX*NGLLY*NGLLZ*NSPEC*nmod))
    m_try_impedance=0.0_CUSTOM_REAL
    allocate(m_new_velocity(NGLLX,NGLLY,NGLLZ,NSPEC,nmod))
    m_new_velocity=0.0_CUSTOM_REAL

    ! p_new_modulus=reshape(p_new,shape(p_new_modulus))
    temp_store=reshape(m_new,shape(temp_store))
    !! Transform velocity models into impedance_I models
    ! IP
    m_new_impedance(:,:,:,:,1) = temp_store(:,:,:,:,3)*temp_store(:,:,:,:,1)
    ! IS
    m_new_impedance(:,:,:,:,2) = temp_store(:,:,:,:,3)*temp_store(:,:,:,:,2)
    ! rhopp
    m_new_impedance(:,:,:,:,3) = temp_store(:,:,:,:,3)
    !! update impedance density parameters
    temp_store_impedance = reshape(m_new_impedance,shape(temp_store_impedance))
    !! transform into impedance sensitivity kernels
    !! K_IP=K_alpha, K_IS=K_beta, K_rhopp=-K_alpha-k_beta+K_rhop
    temp_store=reshape(p_new,shape(temp_store))
    temp_store(:,:,:,:,3) = temp_store(:,:,:,:,3)- temp_store(:,:,:,:,2) - temp_store(:,:,:,:,1)
    p_new=reshape(temp_store,shape(p_new))
    m_try_impedance = temp_store_impedance * (1+step_length*p_new)
    !! transform updated impedance models into velocity models
    temp_store=reshape(m_try_impedance,shape(temp_store))
    m_new_velocity(:,:,:,:,1) = temp_store(:,:,:,:,1)/temp_store(:,:,:,:,3)
    m_new_velocity(:,:,:,:,2) = temp_store(:,:,:,:,2)/temp_store(:,:,:,:,3) 
    m_new_velocity(:,:,:,:,3) = temp_store(:,:,:,:,3)
    m_try = reshape(m_new_velocity,shape(m_try))
!    m_try = m_new * (1+step_length*p_new)
    if(DISPLAY_DETAILS .and. myrank==0) print *,'Min / Max m_try = ', &
        minval(m_try(:)),maxval(m_try(:))

end subroutine update_impedance_I

subroutine update_impedance_II(directory)
    use seismo_parameters
    implicit none
    integer :: ier,imod
    character(len=MAX_FILENAME_LEN) :: filename
    character(len=MAX_STRING_LEN) :: directory
    real(kind=CUSTOM_REAL), dimension(:,:,:,:,:), allocatable :: temp_store
    real(kind=CUSTOM_REAL), dimension(:), allocatable :: temp_store_impedance
    real(kind=CUSTOM_REAL), dimension(:,:,:,:,:), allocatable :: m_new_impedance
    real(kind=CUSTOM_REAL), dimension(:), allocatable :: m_try_impedance
    real(kind=CUSTOM_REAL), dimension(:,:,:,:,:), allocatable :: m_new_velocity
    !! LOAD p_new
    write(filename,'(a)') &
        trim(directory)//'/optimizer/p_new.bin'
    print*,'LOAD p_new -- ', trim(filename)
    open(unit=IIN,file=trim(filename),status='old',action='read',form='unformatted',iostat=ier)
    if (ier /= 0) then
        print*, 'Error: could not open model file: ',trim(filename)
        stop 'Error reading neighbors external mesh file'
    endif
    read(IIN) p_new
    close(IIN)

    !! LOAD m_new
    write(filename,'(a)') &
        trim(directory)//'/optimizer/m_new.bin'
    print*,'LOAD m_new -- ', trim(filename)
    open(unit=IIN,file=trim(filename),status='old',action='read',form='unformatted',iostat=ier)
    if (ier /= 0) then
        print*, 'Error: could not open model file: ',trim(filename)
        stop 'Error reading neighbors external mesh file'
    endif
    read(IIN) m_new
    close(IIN)
    if(DISPLAY_DETAILS .and. myrank==0) print *,'Min / Max m_new = ', &
        minval(m_new(:)),maxval(m_new(:))
    !!   step_length=100
    !! update
    allocate(temp_store(NGLLX,NGLLY,NGLLZ,NSPEC,nmod))
    temp_store=0.0_CUSTOM_REAL
    allocate(temp_store_impedance(NGLLX*NGLLY*NGLLZ*NSPEC*nmod))
    temp_store_impedance=0.0_CUSTOM_REAL
    allocate(m_new_impedance(NGLLX,NGLLY,NGLLZ,NSPEC,nmod))
    m_new_impedance=0.0_CUSTOM_REAL
    allocate(m_try_impedance(NGLLX*NGLLY*NGLLZ*NSPEC*nmod))
    m_try_impedance=0.0_CUSTOM_REAL
    allocate(m_new_velocity(NGLLX,NGLLY,NGLLZ,NSPEC,nmod))
    m_new_velocity=0.0_CUSTOM_REAL

    ! p_new_modulus=reshape(p_new,shape(p_new_modulus))
    temp_store=reshape(m_new,shape(temp_store))
    !! Transform velocity models into impedance-II models
    ! alphap
    m_new_impedance(:,:,:,:,1) = temp_store(:,:,:,:,1)
    ! betap
    m_new_impedance(:,:,:,:,2) = temp_store(:,:,:,:,2)
    ! IPp
    m_new_impedance(:,:,:,:,3) = temp_store(:,:,:,:,1)*temp_store(:,:,:,:,3)
    !! update impedance density parameters
    temp_store_impedance = reshape(m_new_impedance,shape(temp_store_impedance))
    !! transform into impedance sensitivity kernels
    !! K_alphap=K_alpha-K_rhop, K_betap=K_beta, K_IPp=K_rhop
    temp_store=reshape(p_new,shape(temp_store))
    temp_store(:,:,:,:,1) = temp_store(:,:,:,:,1)- temp_store(:,:,:,:,3)
    p_new=reshape(temp_store,shape(p_new))
    m_try_impedance = temp_store_impedance * (1+step_length*p_new)
    !! transform updated impedance-II models into velocity models
    temp_store=reshape(m_try_impedance,shape(temp_store))
    m_new_velocity(:,:,:,:,1) = temp_store(:,:,:,:,1)
    m_new_velocity(:,:,:,:,2) = temp_store(:,:,:,:,2)
    m_new_velocity(:,:,:,:,3) = temp_store(:,:,:,:,3)/temp_store(:,:,:,:,1)
    m_try = reshape(m_new_velocity,shape(m_try))
!    m_try = m_new * (1+step_length*p_new)
    if(DISPLAY_DETAILS .and. myrank==0) print *,'Min / Max m_try = ', &
        minval(m_try(:)),maxval(m_try(:))

end subroutine update_impedance_II

subroutine update_impedance_III(directory)
    use seismo_parameters
    implicit none
    integer :: ier,imod
    character(len=MAX_FILENAME_LEN) :: filename
    character(len=MAX_STRING_LEN) :: directory
    real(kind=CUSTOM_REAL), dimension(:,:,:,:,:), allocatable :: temp_store
    real(kind=CUSTOM_REAL), dimension(:), allocatable :: temp_store_impedance
    real(kind=CUSTOM_REAL), dimension(:,:,:,:,:), allocatable :: m_new_impedance
    real(kind=CUSTOM_REAL), dimension(:), allocatable :: m_try_impedance
    real(kind=CUSTOM_REAL), dimension(:,:,:,:,:), allocatable :: m_new_velocity
    !! LOAD p_new
    write(filename,'(a)') &
        trim(directory)//'/optimizer/p_new.bin'
    print*,'LOAD p_new -- ', trim(filename)
    open(unit=IIN,file=trim(filename),status='old',action='read',form='unformatted',iostat=ier)
    if (ier /= 0) then
        print*, 'Error: could not open model file: ',trim(filename)
        stop 'Error reading neighbors external mesh file'
    endif
    read(IIN) p_new
    close(IIN)

    !! LOAD m_new
    write(filename,'(a)') &
        trim(directory)//'/optimizer/m_new.bin'
    print*,'LOAD m_new -- ', trim(filename)
    open(unit=IIN,file=trim(filename),status='old',action='read',form='unformatted',iostat=ier)
    if (ier /= 0) then
        print*, 'Error: could not open model file: ',trim(filename)
        stop 'Error reading neighbors external mesh file'
    endif
    read(IIN) m_new
    close(IIN)
    if(DISPLAY_DETAILS .and. myrank==0) print *,'Min / Max m_new = ', &
        minval(m_new(:)),maxval(m_new(:))
    !!   step_length=100
    !! update
    allocate(temp_store(NGLLX,NGLLY,NGLLZ,NSPEC,nmod))
    temp_store=0.0_CUSTOM_REAL
    allocate(temp_store_impedance(NGLLX*NGLLY*NGLLZ*NSPEC*nmod))
    temp_store_impedance=0.0_CUSTOM_REAL
    allocate(m_new_impedance(NGLLX,NGLLY,NGLLZ,NSPEC,nmod))
    m_new_impedance=0.0_CUSTOM_REAL
    allocate(m_try_impedance(NGLLX*NGLLY*NGLLZ*NSPEC*nmod))
    m_try_impedance=0.0_CUSTOM_REAL
    allocate(m_new_velocity(NGLLX,NGLLY,NGLLZ,NSPEC,nmod))
    m_new_velocity=0.0_CUSTOM_REAL

    ! p_new_modulus=reshape(p_new,shape(p_new_modulus))
    temp_store=reshape(m_new,shape(temp_store))
    !! Transform velocity models into impedance-II models
    ! alphap
    m_new_impedance(:,:,:,:,1) = temp_store(:,:,:,:,1)
    ! betap
    m_new_impedance(:,:,:,:,2) = temp_store(:,:,:,:,2)
    ! IPp
    m_new_impedance(:,:,:,:,3) = temp_store(:,:,:,:,2)*temp_store(:,:,:,:,3)
    !! update impedance density parameters
    temp_store_impedance = reshape(m_new_impedance,shape(temp_store_impedance))
    !! transform into impedance sensitivity kernels
    !! K_alphap=K_alpha-K_rhop, K_betap=K_beta, K_IPp=K_rhop
    temp_store=reshape(p_new,shape(temp_store))
    temp_store(:,:,:,:,2) = temp_store(:,:,:,:,2)- temp_store(:,:,:,:,3)
    p_new=reshape(temp_store,shape(p_new))
    m_try_impedance = temp_store_impedance * (1+step_length*p_new)
    !! transform updated impedance-II models into velocity models
    temp_store=reshape(m_try_impedance,shape(temp_store))
    m_new_velocity(:,:,:,:,1) = temp_store(:,:,:,:,1)
    m_new_velocity(:,:,:,:,2) = temp_store(:,:,:,:,2)
    m_new_velocity(:,:,:,:,3) = temp_store(:,:,:,:,3)/temp_store(:,:,:,:,2)
    m_try = reshape(m_new_velocity,shape(m_try))
!    m_try = m_new * (1+step_length*p_new)
    if(DISPLAY_DETAILS .and. myrank==0) print *,'Min / Max m_try = ', &
        minval(m_try(:)),maxval(m_try(:))

end subroutine update_impedance_III
!! Impedance_new IP, IS, VS
subroutine update_impedance_new(directory)
    use seismo_parameters
    implicit none
    integer :: ier,imod
    character(len=MAX_FILENAME_LEN) :: filename
    character(len=MAX_STRING_LEN) :: directory
    real(kind=CUSTOM_REAL), dimension(:,:,:,:,:), allocatable :: temp_store
    real(kind=CUSTOM_REAL), dimension(:), allocatable :: temp_store_impedance
    real(kind=CUSTOM_REAL), dimension(:,:,:,:,:), allocatable :: m_new_impedance
    real(kind=CUSTOM_REAL), dimension(:), allocatable :: m_try_impedance
    real(kind=CUSTOM_REAL), dimension(:,:,:,:,:), allocatable :: m_new_velocity
    real(kind=CUSTOM_REAL), dimension(:,:,:,:,:), allocatable :: p_new_impedance_new
    !! LOAD p_new
    write(filename,'(a)') &
        trim(directory)//'/optimizer/p_new.bin'
    print*,'LOAD p_new -- ', trim(filename)
    open(unit=IIN,file=trim(filename),status='old',action='read',form='unformatted',iostat=ier)
    if (ier /= 0) then
        print*, 'Error: could not open model file: ',trim(filename)
        stop 'Error reading neighbors external mesh file'
    endif
    read(IIN) p_new
    close(IIN)

    !! LOAD m_new
    write(filename,'(a)') &
        trim(directory)//'/optimizer/m_new.bin'
    print*,'LOAD m_new -- ', trim(filename)
    open(unit=IIN,file=trim(filename),status='old',action='read',form='unformatted',iostat=ier)
    if (ier /= 0) then
        print*, 'Error: could not open model file: ',trim(filename)
        stop 'Error reading neighbors external mesh file'
    endif
    read(IIN) m_new
    close(IIN)
    if(DISPLAY_DETAILS .and. myrank==0) print *,'Min / Max m_new = ', &
        minval(m_new(:)),maxval(m_new(:))
    !!   step_length=100
    !! update
    allocate(temp_store(NGLLX,NGLLY,NGLLZ,NSPEC,nmod))
    temp_store=0.0_CUSTOM_REAL
    allocate(temp_store_impedance(NGLLX*NGLLY*NGLLZ*NSPEC*nmod))
    temp_store_impedance=0.0_CUSTOM_REAL
    allocate(m_new_impedance(NGLLX,NGLLY,NGLLZ,NSPEC,nmod))
    m_new_impedance=0.0_CUSTOM_REAL
    allocate(m_try_impedance(NGLLX*NGLLY*NGLLZ*NSPEC*nmod))
    m_try_impedance=0.0_CUSTOM_REAL
    allocate(m_new_velocity(NGLLX,NGLLY,NGLLZ,NSPEC,nmod))
    m_new_velocity=0.0_CUSTOM_REAL
    allocate(p_new_impedance_new(NGLLX,NGLLY,NGLLZ,NSPEC,nmod))
    p_new_impedance_new=0.0_CUSTOM_REAL

    ! p_new_modulus=reshape(p_new,shape(p_new_modulus))
    temp_store=reshape(m_new,shape(temp_store))
    !! Transform velocity models into impedance_I models
    ! IP
    m_new_impedance(:,:,:,:,1) = temp_store(:,:,:,:,3)*temp_store(:,:,:,:,1)
    ! IS
    m_new_impedance(:,:,:,:,2) = temp_store(:,:,:,:,3)*temp_store(:,:,:,:,2)
    ! betappp
    m_new_impedance(:,:,:,:,3) = temp_store(:,:,:,:,2)
    !! update impedance density parameters
    temp_store_impedance = reshape(m_new_impedance,shape(temp_store_impedance))
    !! transform into impedance sensitivity kernels
    !! K_IP=K_alpha, K_IS=-K_alpha+K_rhop, K_betappp=K_alpha+k_beta-K_rhop
    temp_store=reshape(p_new,shape(temp_store))
    p_new_impedance_new(:,:,:,:,1) = temp_store(:,:,:,:,1)
    p_new_impedance_new(:,:,:,:,2) = -temp_store(:,:,:,:,1) + temp_store(:,:,:,:,3)
    p_new_impedance_new(:,:,:,:,3) = temp_store(:,:,:,:,1) + temp_store(:,:,:,:,2) - temp_store(:,:,:,:,3)
    p_new=reshape(p_new_impedance_new,shape(p_new))
    m_try_impedance = temp_store_impedance * (1+step_length*p_new)
    !! transform updated impedance models into velocity models
    temp_store=reshape(m_try_impedance,shape(temp_store))
    m_new_velocity(:,:,:,:,1) = temp_store(:,:,:,:,1)*temp_store(:,:,:,:,3)/temp_store(:,:,:,:,2)
    m_new_velocity(:,:,:,:,2) = temp_store(:,:,:,:,3)
    m_new_velocity(:,:,:,:,3) = temp_store(:,:,:,:,2)/temp_store(:,:,:,:,3)
    m_try = reshape(m_new_velocity,shape(m_try))
!    m_try = m_new * (1+step_length*p_new)
    if(DISPLAY_DETAILS .and. myrank==0) print *,'Min / Max m_try = ', &
        minval(m_try(:)),maxval(m_try(:))

end subroutine update_impedance_new
! Update HTI thomsen parameterization, alpha, beta, epsilon, delta, rho
subroutine update_hti_thom(directory)
    use seismo_parameters
    implicit none
    integer :: ier,imod
    character(len=MAX_FILENAME_LEN) :: filename
    character(len=MAX_STRING_LEN) :: directory
    real(kind=CUSTOM_REAL), dimension(:,:,:,:,:), allocatable :: temp_store
    real(kind=CUSTOM_REAL), dimension(:), allocatable :: temp_store_hti_thom
    real(kind=CUSTOM_REAL), dimension(:,:,:,:,:), allocatable :: m_new_hti_thom
    real(kind=CUSTOM_REAL), dimension(:), allocatable :: m_try_hti_thom
    real(kind=CUSTOM_REAL), dimension(:,:,:,:,:), allocatable :: m_new_hti_ec
    real(kind=CUSTOM_REAL), dimension(:,:,:,:,:), allocatable :: p_new_hti_thom
    !! LOAD p_new
    write(filename,'(a)') &
        trim(directory)//'/optimizer/p_new.bin'
    print*,'LOAD p_new -- ', trim(filename)
    open(unit=IIN,file=trim(filename),status='old',action='read',form='unformatted',iostat=ier)
    if (ier /= 0) then
        print*, 'Error: could not open model file: ',trim(filename)
        stop 'Error reading neighbors external mesh file'
    endif
    read(IIN) p_new
    close(IIN)

    !! LOAD m_new
    write(filename,'(a)') &
        trim(directory)//'/optimizer/m_new.bin'
    print*,'LOAD m_new -- ', trim(filename)
    open(unit=IIN,file=trim(filename),status='old',action='read',form='unformatted',iostat=ier)
    if (ier /= 0) then
        print*, 'Error: could not open model file: ',trim(filename)
        stop 'Error reading neighbors external mesh file'
    endif
    read(IIN) m_new
    close(IIN)
    
    if(DISPLAY_DETAILS .and. myrank==0) print *,'Min / Max m_new = ', &
        minval(m_new(:)),maxval(m_new(:))
    !!   step_length=100
    !! update
    allocate(temp_store(NGLLX,NGLLY,NGLLZ,NSPEC,nmod))
    temp_store=0.0_CUSTOM_REAL
    allocate(temp_store_hti_thom(NGLLX*NGLLY*NGLLZ*NSPEC*nmod))
    temp_store_hti_thom=0.0_CUSTOM_REAL
    allocate(m_new_hti_thom(NGLLX,NGLLY,NGLLZ,NSPEC,nmod))
    m_new_hti_thom=0.0_CUSTOM_REAL
    allocate(m_try_hti_thom(NGLLX*NGLLY*NGLLZ*NSPEC*nmod))
    m_try_hti_thom=0.0_CUSTOM_REAL
    allocate(m_new_hti_ec(NGLLX,NGLLY,NGLLZ,NSPEC,nmod))
    m_new_hti_ec=0.0_CUSTOM_REAL
    allocate(p_new_hti_thom(NGLLX,NGLLY,NGLLZ,NSPEC,nmod))
    p_new_hti_thom=0.0_CUSTOM_REAL

    temp_store=reshape(m_new,shape(temp_store))
    !! Transform elastic constants models into thomsen parameterization models
    ! alpha
    m_new_hti_thom(:,:,:,:,1) = SQRT(temp_store(:,:,:,:,3)/temp_store(:,:,:,:,5))
    ! beta
    m_new_hti_thom(:,:,:,:,2) = SQRT(temp_store(:,:,:,:,4)/temp_store(:,:,:,:,5))
    ! epsilon
    m_new_hti_thom(:,:,:,:,3) = (temp_store(:,:,:,:,1)-temp_store(:,:,:,:,3))/&
                                (2.0_CUSTOM_REAL*temp_store(:,:,:,:,3))
    ! delta
    m_new_hti_thom(:,:,:,:,4) = ((temp_store(:,:,:,:,2)+temp_store(:,:,:,:,4))**2)/(2.0_CUSTOM_REAL*temp_store(:,:,:,:,3)*&
      (temp_store(:,:,:,:,3)-temp_store(:,:,:,:,4)))-(temp_store(:,:,:,:,3)-temp_store(:,:,:,:,4))/&
      (2.0_CUSTOM_REAL*temp_store(:,:,:,:,3))
    ! rhop
    m_new_hti_thom(:,:,:,:,5)=temp_store(:,:,:,:,5)
    !! update thomsen parameters
    temp_store_hti_thom = reshape(m_new_hti_thom,shape(temp_store_hti_thom))
    !!
    !p_new=reshape(p_new_hti_thom,shape(p_new))
    m_try_hti_thom = temp_store_hti_thom * (1+step_length*p_new)
    !! transform updated thomsen models into elastic constants models
    temp_store=reshape(m_try_hti_thom,shape(temp_store))
    ! c11=2*rho*alpha**2(2*epsilon+1)
    m_new_hti_ec(:,:,:,:,1) = (temp_store(:,:,:,:,5)*temp_store(:,:,:,:,1)*temp_store(:,:,:,:,1))*&
                              (2.0_CUSTOM_REAL*temp_store(:,:,:,:,3)+1.0_CUSTOM_REAL)
    ! c13=sqrt(2rho*alpha**2*delta*(rho*alpha**2-rho*beta**2)+(rho*alpha**2-rho*beta**2)**2)-rho*beta**2
    m_new_hti_ec(:,:,:,:,2) = SQRT(2.0_CUSTOM_REAL*temp_store(:,:,:,:,5)*temp_store(:,:,:,:,1)*&
         temp_store(:,:,:,:,1)*(temp_store(:,:,:,:,5)*temp_store(:,:,:,:,1)*temp_store(:,:,:,:,1)-&
         temp_store(:,:,:,:,5)*temp_store(:,:,:,:,2)*temp_store(:,:,:,:,2))*temp_store(:,:,:,:,4)+(&
         temp_store(:,:,:,:,5)*temp_store(:,:,:,:,1)*temp_store(:,:,:,:,1)-&
         temp_store(:,:,:,:,5)*temp_store(:,:,:,:,2)*temp_store(:,:,:,:,2))*&
         (temp_store(:,:,:,:,5)*temp_store(:,:,:,:,1)*temp_store(:,:,:,:,1)-&
         temp_store(:,:,:,:,5)*temp_store(:,:,:,:,2)*temp_store(:,:,:,:,2)))-&
         temp_store(:,:,:,:,5)*temp_store(:,:,:,:,2)*temp_store(:,:,:,:,2)
    ! c33
    m_new_hti_ec(:,:,:,:,3) = temp_store(:,:,:,:,5)*temp_store(:,:,:,:,1)*temp_store(:,:,:,:,1)
    ! c55
    m_new_hti_ec(:,:,:,:,4) = temp_store(:,:,:,:,5)*temp_store(:,:,:,:,2)*temp_store(:,:,:,:,2)
    ! rho
    m_new_hti_ec(:,:,:,:,5) = temp_store(:,:,:,:,5)

    m_try = reshape(m_new_hti_ec,shape(m_try))
!    m_try = m_new * (1+step_length*p_new)
    if(DISPLAY_DETAILS .and. myrank==0) print *,'Min / Max m_try = ', &
        minval(m_try(:)),maxval(m_try(:))

end subroutine update_hti_thom

! vp, sigma=vp/vs, rho parameterization
subroutine update_sigma(directory)
    use seismo_parameters
    implicit none
    integer :: ier,imod
    character(len=MAX_FILENAME_LEN) :: filename
    character(len=MAX_STRING_LEN) :: directory
    real(kind=CUSTOM_REAL), dimension(:,:,:,:,:), allocatable :: temp_store
    real(kind=CUSTOM_REAL), dimension(:), allocatable :: temp_store_sigma
    real(kind=CUSTOM_REAL), dimension(:,:,:,:,:), allocatable :: m_new_sigma
    real(kind=CUSTOM_REAL), dimension(:), allocatable :: m_try_sigma
    real(kind=CUSTOM_REAL), dimension(:,:,:,:,:), allocatable :: m_new_velocity
    !! LOAD p_new
    write(filename,'(a)') &
        trim(directory)//'/optimizer/p_new.bin'
    print*,'LOAD p_new -- ', trim(filename)
    open(unit=IIN,file=trim(filename),status='old',action='read',form='unformatted',iostat=ier)
    if (ier /= 0) then
        print*, 'Error: could not open model file: ',trim(filename)
        stop 'Error reading neighbors external mesh file'
    endif
    read(IIN) p_new
    close(IIN)

    !! LOAD m_new
    write(filename,'(a)') &
        trim(directory)//'/optimizer/m_new.bin'
    print*,'LOAD m_new -- ', trim(filename)
    open(unit=IIN,file=trim(filename),status='old',action='read',form='unformatted',iostat=ier)
    if (ier /= 0) then
        print*, 'Error: could not open model file: ',trim(filename)
        stop 'Error reading neighbors external mesh file'
    endif
    read(IIN) m_new
    close(IIN)
    if(DISPLAY_DETAILS .and. myrank==0) print *,'Min / Max m_new = ', &
        minval(m_new(:)),maxval(m_new(:))

    allocate(temp_store(NGLLX,NGLLY,NGLLZ,NSPEC,nmod))
    temp_store=0.0_CUSTOM_REAL
    allocate(temp_store_sigma(NGLLX*NGLLY*NGLLZ*NSPEC*nmod))
    temp_store_sigma=0.0_CUSTOM_REAL
    allocate(m_new_sigma(NGLLX,NGLLY,NGLLZ,NSPEC,nmod))
    m_new_sigma=0.0_CUSTOM_REAL
    allocate(m_try_sigma(NGLLX*NGLLY*NGLLZ*NSPEC*nmod))
    m_try_sigma=0.0_CUSTOM_REAL
    allocate(m_new_velocity(NGLLX,NGLLY,NGLLZ,NSPEC,nmod))
    m_new_velocity=0.0_CUSTOM_REAL

    ! p_new_modulus=reshape(p_new,shape(p_new_modulus))
    temp_store=reshape(m_new,shape(temp_store))
    !! Transform velocity models into sigma models
    ! alphap
    m_new_sigma(:,:,:,:,1) = temp_store(:,:,:,:,1)
    ! betap
    m_new_sigma(:,:,:,:,2) = temp_store(:,:,:,:,1)/temp_store(:,:,:,:,2)
    ! IPp
    m_new_sigma(:,:,:,:,3) = temp_store(:,:,:,:,3)
    !! update sigma model parameters
    temp_store_sigma = reshape(m_new_sigma,shape(temp_store_sigma))
    !! transform into sigma sensitivity kernels
    !! K_{\hat{alpha}}=K_alpha+K_beta, K_{sigma}=-K_beta, K_{\hat{rho}}=K_rhop
    temp_store=reshape(p_new,shape(temp_store))
    temp_store(:,:,:,:,1) = temp_store(:,:,:,:,1)+ temp_store(:,:,:,:,2)
    temp_store(:,:,:,:,2) = -temp_store(:,:,:,:,2)
    p_new=reshape(temp_store,shape(p_new))
    m_try_sigma = temp_store_sigma * (1+step_length*p_new)
    !! transform updated sigma models into velocity models
    temp_store=reshape(m_try_sigma,shape(temp_store))
    m_new_velocity(:,:,:,:,1) = temp_store(:,:,:,:,1)
    m_new_velocity(:,:,:,:,2) = temp_store(:,:,:,:,1)/temp_store(:,:,:,:,2)
    m_new_velocity(:,:,:,:,3) = temp_store(:,:,:,:,3)
    m_try = reshape(m_new_velocity,shape(m_try))
!    m_try = m_new * (1+step_length*p_new)
    if(DISPLAY_DETAILS .and. myrank==0) print *,'Min / Max m_try = ', &
        minval(m_try(:)),maxval(m_try(:))

end subroutine update_sigma

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine finalize(directory,model_names)
    use seismo_parameters
    implicit none
    integer :: ier,imod
    integer :: nspec_start,nspec_end
    character(len=MAX_STRING_LEN) :: model_names(nmod)
    character(len=MAX_FILENAME_LEN) :: filename
    character(len=MAX_STRING_LEN) :: directory
    real(kind=CUSTOM_REAL), dimension(:,:,:,:,:),allocatable :: temp_store

    allocate(temp_store(NGLLX,NGLLY,NGLLZ,NSPEC,nmod))
    temp_store = 0.0_CUSTOM_REAL

    temp_store=reshape(m_try,shape(temp_store))

    do myrank=0,nproc-1
    nspec_start=sum(nspec_proc(1:myrank))+1
    nspec_end=sum(nspec_proc(1:myrank))+nspec_proc(myrank+1)
    do imod=1,nmod
    write(filename,'(a,i6.6,a)') &
        trim(directory)//'/m_try/proc',myrank,'_'//&
        trim(model_names(imod))//'.bin'
    if (myrank == 0) print*,'SAVE m_try -- ', trim(filename)
    open(unit=IOUT,file=trim(filename),status='unknown',form='unformatted',iostat=ier)
    if (ier /= 0) then
        print*, 'Error: could not open gradient file: ',trim(filename)
        stop 'Error reading neighbors external mesh file'
    endif
    write(IOUT) temp_store(:,:,:,nspec_start:nspec_end,imod)
    close(IOUT)
    enddo ! imod
    enddo ! myrank

    deallocate(temp_store)
    deallocate(m_new)
    deallocate(m_try)
    deallocate(p_new)
    deallocate(nspec_proc)
end subroutine finalize 

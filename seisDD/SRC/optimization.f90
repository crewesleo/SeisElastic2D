program optimization
    ! To calculate gradient and update direction
    ! yanhuay@princeton.edu
    ! modified by PWY for diagonal pseudo Hessian preconditioning in anisotropic
    ! media
    use seismo_parameters
    implicit none

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
    integer, parameter :: NARGS = 6
    INTEGER :: itime, ier, isrc,i,j,iter
    character(len=MAX_STRING_LEN) :: kernel_names(MAX_KERNEL_NUM)
    character(len=MAX_STRING_LEN) :: kernel_names_comma_delimited
    character(len=MAX_STRING_LEN) :: model_names(MAX_KERNEL_NUM)
    character(len=MAX_STRING_LEN) :: model_names_comma_delimited
    character(len=MAX_STRING_LEN) :: precond_name
    character(len=MAX_STRING_LEN) :: pdh_names(MAX_KERNEL_NUM)! diagonal pseudo Hessian names
    character(len=MAX_STRING_LEN) :: arg(NARGS)
    character(len=MAX_STRING_LEN) :: directory
    real t1,t2
    character, parameter :: delimiter=','

    call cpu_time(t1)

    ! parse command line arguments
    if (command_argument_count() /= NARGS) then
        if (DISPLAY_DETAILS) then
            print *, 'USAGE:  mpirun -np NPROC bin/gradient.exe ...'
            stop ' Please check command line arguments'
        endif
    endif

    do i = 1, NARGS
    call get_command_argument(i,arg(i), status=ier)
    enddo

    read(arg(1),*) nproc
    directory=arg(2) 
    kernel_names_comma_delimited = arg(3)
    precond_name=arg(4)
    model_names_comma_delimited = arg(5)
    read(arg(6),*) iter

    call split_string(kernel_names_comma_delimited,delimiter,kernel_names,nker)
    call split_string(model_names_comma_delimited,delimiter,model_names,nmod)
    call split_string(precond_name,delimiter,pdh_names,npre) ! PWY

    if (nker .ne. nmod) then
        print*, 'number of kernel ',nker,' is not equal to number of model ',nmod
        stop
    endif
    if(precond) then
        print*,'optimization with preconditioning'
        print*,'preconditioner -- ', trim(adjustl(precond_name))
        print*,'preconditioner -- ', trim(adjustl(pdh_names(1)))
    else
        print*,'optimization without preconditioning'
    endif

    !! initialization  -- get number of spectral elements
!    call initialize(directory,adjustl(kernel_names(1:NKER)),&
!        adjustl(precond_name),adjustl(model_names(1:nmod))) 
    ! PWY
    call initialize(directory,adjustl(kernel_names(1:NKER)),adjustl(precond_name),adjustl(model_names(1:nmod)))

    !! optimization(update) direction
    call update_direction(directory,adjustl(kernel_names(1:NKER)),iter)

    !! save update direction 
    call finalize(directory)   

    call cpu_time(t2)

    if (DISPLAY_DETAILS .and. myrank==0) print *,'Computation time with CPU:',t2-t1

end program optimization

subroutine initialize(directory,kernel_names,precond_name,model_names)
!subroutine initialize(directory,kernel_names,pdh_names,model_names)

    use seismo_parameters
    implicit none
    integer :: ier,iker,imod
    integer :: filesize,nspec_start,nspec_end
    real(kind=CUSTOM_REAL) :: wtr
    real(kind=CUSTOM_REAL) :: max_kernel_before_pre
    real(kind=CUSTOM_REAL) :: max_kernel_after_pre
    character(len=MAX_FILENAME_LEN) :: filename
    character(len=MAX_STRING_LEN) :: directory
    character(len=MAX_STRING_LEN) :: kernel_names(nker)
    character(len=MAX_STRING_LEN) :: model_names(nmod)
    character(len=MAX_STRING_LEN) :: precond_name
    character(len=MAX_STRING_LEN) :: pdh_names(nker)
    character(len=MAX_FILENAME_LEN) :: filename_pre
    real(kind=CUSTOM_REAL), dimension(:,:,:,:,:),allocatable :: temp_store
    !real(kind=CUSTOM_REAL), dimension(:,:,:,:),allocatable :: preconditioner
    real(kind=CUSTOM_REAL), dimension(:,:,:,:,:),allocatable :: preconditioner
    real(kind=CUSTOM_REAL), dimension(:,:,:,:,:),allocatable :: preconditioner_norm

    !! Modulue-density parameterization
    real(kind=CUSTOM_REAL), dimension(:,:,:,:),allocatable :: temp_store_kappa
    real(kind=CUSTOM_REAL), dimension(:,:,:,:),allocatable :: temp_store_mu
    real(kind=CUSTOM_REAL), dimension(:,:,:,:),allocatable :: temp_store_rho
    real(kind=CUSTOM_REAL), dimension(:,:,:,:,:),allocatable :: temp_store_modulus

    character, parameter :: delimiter_pre=','

    call split_string(precond_name,delimiter_pre,pdh_names,npre) ! PWY
!    pdh_names(1)='pdh_hti_thom_alpha'
!    pdh_names(2)='pdh_hti_thom_beta'
!    pdh_names(3)='pdh_hti_thom_epsilon'
!    pdh_names(4)='pdh_hti_thom_delta'
!    pdh_names(5)='pdh_hti_thom_rhop'
!    print*,'Diagonal pseudo Hessian names ',trim(pdh_names(1))
!    print*,'Diagonal pseudo Hessian names ',trim(pdh_names(2))
!    print*,'Diagonal pseudo Hessian names ',trim(pdh_names(3))
!    print*,'Diagonal pseudo Hessian names ',trim(pdh_names(4))
!    print*,'Diagonal pseudo Hessian names ',trim(pdh_names(5))
    ! slice database file
    allocate(nspec_proc(nproc))
    nspec_proc=0  

    do myrank=0,nproc-1

    ! nspec
    write(filename,'(a,i6.6,a)') &
        trim(directory)//'/misfit_kernel/proc',myrank,'_'//trim(IBOOL_NAME) 
    open(IIN,file=trim(filename),status='old',action='read',form='unformatted',iostat=ier)
    if (ier /= 0) then
        print *,'Error: could not open database file:',trim(filename)
        stop 'Error opening _NSPEC_IBOOL file'
    endif
    read(IIN) nspec_proc(myrank+1)
    close(IIN)

    if(DISPLAY_DETAILS .and. myrank==0) print*, 'nspec_proc=',nspec_proc(myrank+1)

    enddo

    nspec=sum(nspec_proc)
    if(DISPLAY_DETAILS) print*,'NGLLX*NGLLY*NGLLZ*NSPEC*nmod:',NGLLX,NGLLY,NGLLZ,NSPEC,nmod

    !! gloabl
    allocate(g_new(NGLLX*NGLLY*NGLLZ*NSPEC*nmod))
    g_new = 0.0_CUSTOM_REAL
    allocate(p_new(NGLLX*NGLLY*NGLLZ*NSPEC*nmod))
    p_new = 0.0_CUSTOM_REAL
    allocate(m_new(NGLLX*NGLLY*NGLLZ*NSPEC*nmod))
    m_new = 0.0_CUSTOM_REAL

    !! local
    allocate(temp_store(NGLLX,NGLLY,NGLLZ,NSPEC,nmod))
    allocate(preconditioner(NGLLX,NGLLY,NGLLZ,NSPEC,nker))
    allocate(preconditioner_norm(NGLLX,NGLLY,NGLLZ,NSPEC,nker))

    temp_store = 0.0_CUSTOM_REAL
    preconditioner = 0.0_CUSTOM_REAL
    preconditioner_norm = 0.0_CUSTOM_REAL

    !! prepare g_new
    do myrank=0,nproc-1

    nspec_start=sum(nspec_proc(1:myrank))+1
    nspec_end=sum(nspec_proc(1:myrank))+nspec_proc(myrank+1)

    do iker=1,nker
    if(smooth) then
        write(filename,'(a,i6.6,a)') &
            trim(directory)//'/misfit_kernel/proc',myrank,&
            '_'//trim(kernel_names(iker))//'_smooth.bin'
    else
        write(filename,'(a,i6.6,a)') &
            trim(directory)//'/misfit_kernel/proc',myrank,&
            '_'//trim(kernel_names(iker))//'.bin'
    endif
    if(myrank==0) print*,'LOAD misfit_kernel -- ',trim(filename)
    open(IIN,file=trim(filename),status='old',action='read',form='unformatted',iostat=ier)
    if (ier /= 0) then
        print*, 'Error: could not open gradient file: ',trim(filename)
        stop 'Error: could not open gradient file: '
    endif
    read(IIN) temp_store(:,:,:,nspec_start:nspec_end,iker)
    close(IIN)
!    enddo  ! iker

    ! preconditioner 
    if (precond) then 
        if(smooth) then
            print*,'write Hessian names : PWY ',trim(pdh_names(iker))
            write(filename_pre,'(a,i6.6,a)') &
                trim(directory)//'/misfit_kernel/proc',myrank,&
                '_'//trim(pdh_names(iker))//'_smooth.bin'
        else
            write(filename_pre,'(a,i6.6,a)') &
                trim(directory)//'/misfit_kernel/proc',myrank,&
                '_'//trim(pdh_names(iker))//'.bin'
        endif
        if(myrank==0) print*,'LOAD hessian_kernel -- ',trim(filename_pre)
        open(IIN,file=trim(filename_pre),status='old',action='read',form='unformatted',iostat=ier)
        if (ier /= 0) then
            print*, 'Error: could not open hessian file: ',trim(filename_pre)
            stop 'Error: could not open hessian file: '
        endif
        read(IIN) preconditioner(:,:,:,nspec_start:nspec_end,iker)
        close(IIN)

    endif ! non-empty preconditioner
    enddo ! iker
    enddo  ! myrank

    !! preconditioning [misfit_kernel --> g_new]
    wtr = maxval(abs(preconditioner(:,:,:,:,:)))
    max_kernel_before_pre=maxval(abs(temp_store(:,:,:,:,:)))
    !wtr_precond = 1.0d-3
    if(wtr>SMALL_VAL .and. precond)then
        do iker=1,nker
        !wtr = maxval(abs(preconditioner(:,:,:,:,iker)))
        ! normalization
        preconditioner_norm(:,:,:,:,iker)=abs(preconditioner(:,:,:,:,iker))/maxval(abs(preconditioner(:,:,:,:,:)))
        wtr = maxval(abs(preconditioner_norm(:,:,:,:,iker)))
        print*,'Number of preconditioner', iker
        print*,'Maximum values of preconditioner after normalization',maxval(preconditioner_norm(:,:,:,:,iker))
        print*,'Maximum values of kernel before preconditioning',maxval(temp_store(:,:,:,:,iker))
        ! preconditioning
        temp_store(:,:,:,:,iker) = &
            temp_store(:,:,:,:,iker) / (preconditioner_norm(:,:,:,:,iker) + wtr_precond * wtr)

        print*,'Maximum values of kernel after preconditioning',maxval(temp_store(:,:,:,:,iker))
        enddo
    endif
    max_kernel_after_pre=maxval(abs(temp_store(:,:,:,:,:)))
    !! recover the magnitudes
    temp_store(:,:,:,:,:)=temp_store(:,:,:,:,:)*max_kernel_before_pre/max_kernel_after_pre
    print*,'Maximum values of kernel after balancing',maxval(temp_store(:,:,:,:,:))

    !! convert to 1D vector
    g_new=reshape(temp_store,shape(g_new))
    if(DISPLAY_DETAILS .and. myrank==0) then
        print *,'myrank=',myrank,' Min / Max g_new = ',&
            minval(g_new(:)),maxval(g_new(:))
    endif

    !! prepare m_new
    temp_store = 0.0_CUSTOM_REAL
    !! prepare g_new
    do myrank=0,nproc-1

    nspec_start=sum(nspec_proc(1:myrank))+1
    nspec_end=sum(nspec_proc(1:myrank))+nspec_proc(myrank+1)
    do imod=1,nmod
    write(filename,'(a,i6.6,a)') &
        trim(directory)//'/m_current/proc',myrank,&
        '_'//trim(model_names(imod))//'.bin'
    if(myrank==0) print*,'LOAD m_current -- ',trim(filename)
    open(IIN,file=trim(filename),status='old',action='read',form='unformatted',iostat=ier)
    if (ier /= 0) then
        print*, 'Error: could not open m_new file: ',trim(filename)
        stop 'Error: could not open m_new file: '
    endif
    read(IIN) temp_store(:,:,:,nspec_start:nspec_end,imod)
    close(IIN)
    enddo ! imod
    enddo ! myrank
    m_new=reshape(temp_store,shape(m_new))

    if(DISPLAY_DETAILS .and. myrank==0) then
        print *,'myrank=',myrank,' Min / Max m_new = ',&
            minval(m_new(:)),maxval(m_new(:))
    endif

    deallocate(temp_store)
    deallocate(preconditioner)
    deallocate(preconditioner_norm)
    deallocate(nspec_proc)
end subroutine initialize
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine update_direction(directory,kernel_names,iter)
    use seismo_parameters
    implicit none
    integer :: iker,ier,iter
    character(len=MAX_STRING_LEN) :: kernel_names(nker)
    real(kind=CUSTOM_REAL), dimension(:,:,:,:),allocatable :: temp_store
    character(len=MAX_FILENAME_LEN) :: filename
    character(len=MAX_STRING_LEN) :: directory
    real(kind=CUSTOM_REAL) :: pmax
    !! CG
    integer :: cgstep
    real(kind=CUSTOM_REAL), dimension(:),allocatable :: g_old
    real(kind=CUSTOM_REAL), dimension(:),allocatable :: p_old
    real(kind=CUSTOM_REAL) :: temp
    !! BFGS
    integer :: BFGS_step,m
    real(kind=CUSTOM_REAL), dimension(:),allocatable :: m_old
    real(kind=CUSTOM_REAL), dimension(:,:),allocatable :: Deltam
    real(kind=CUSTOM_REAL), dimension(:,:),allocatable :: Deltag

    allocate(g_old(size(g_new)))
    allocate(p_old(size(g_new)))
    allocate(m_old(size(g_new)))
    allocate(Deltam(size(g_new),BFGS_stepmax))
    allocate(Deltag(size(g_new),BFGS_stepmax))
    g_old = 0.0_CUSTOM_REAL
    p_old = 0.0_CUSTOM_REAL
    m_old = 0.0_CUSTOM_REAL
    Deltam = 0.0_CUSTOM_REAL
    Deltag = 0.0_CUSTOM_REAL

    select case(opt_scheme)
    case("SD") !! steepest descent method  
        if(myrank==0) print*, 'steepest descent for iter  ',iter
        !! search direction
        call SD(g_new, size(g_new),p_new)
    case ("CG") 
        if(iter==1) then   !! first iter step, do SD
            !! search direction 
            if(myrank==0) print*, 'steepest descent for iter 1 '
            call SD(g_new, size(g_new), p_new)
            cgstep = 1

        else !! not the first iter step, try CG
            if(myrank==0) print*, 'CG for iter ',iter

            ! additional file needed: cgstep
            write(filename,'(a)') trim(directory)//'/optimizer/cgstep.dat'
            OPEN(IIN,FILE=trim(filename),STATUS='old',action='read',iostat=ier)
            if(ier>0) then
                print*, 'Error: could not open cgstep file:',trim(filename)

                stop 'Error: could not open cgstep file:'
            else
                read(IIN,*) cgstep
            end if
            close(IIN)
            if(myrank==0) print*,'LOAD cgstep -- ',trim(filename)

            !! second if   
            if( cgstep > cgstepmax ) then ! exceed max cg steps, do SD
                print*, 'restarting NLCG ... [periodic restart]'
                cgstep = 1
                !! search direction 
                print*, 'steepest descent for restarting iter=',&
                    iter, ' cgstep=', cgstep
                call SD(g_new, size(g_new), p_new)

            elseif(cgstep>=1 .and. cgstep<=cgstepmax) then !! not exceed maxcg steps, try CG 
                ! additional file needed: g_old, p_old                                 
                write(filename,'(a)')  trim(directory)//'/optimizer/g_old.bin'

                print*,'LOAD g_old -- ', trim(filename)
                open(unit=IIN,file=trim(filename),status='old',action='read',form='unformatted',iostat=ier)
                if (ier /= 0) then
                    print*, 'Error: could not open g_old file: ',trim(filename)
                    stop 'Error: could not open g_old file: '
                endif
                read(IIN) g_old
                close(IIN) 

                !! p_old
                write(filename,'(a)')  trim(directory)//'/optimizer/p_old.bin'
                print*,'LOAD p_old -- ', trim(filename)
                open(unit=IIN,file=trim(filename),status='old',action='read',form='unformatted',iostat=ier)
                if (ier /= 0) then
                    print*, 'Error: could not open p_old file: ',trim(filename)
                    stop 'Error: could not open p_old file: '
                endif
                read(IIN) p_old
                close(IIN)
                !! search direction 
                print*, 'conjugate gradient direction for iter=',&
                    iter, ' cgstep=', cgstep
                call NLCG(g_new, g_old, p_old, size(g_new), CG_scheme, cgstep, p_new)

            endif !! cgstep
        endif !! iter==1

        !! save cgstep
        write(filename,'(a)') trim(directory)//'/optimizer/cgstep.dat'
        OPEN(IOUT,FILE=trim(filename),STATUS='unknown',iostat=ier)
        if(ier>0) then
            print*, 'Error: could not open cgstep file:',trim(filename)

            stop 'Error: could not open cgstep file:'
        else
            write(IOUT,*) cgstep
        end if
        close(IOUT)
        print*,'SAVE cgstep -- ',trim(filename)


    case("QN") !! Qausi-Newton (L-BFGS) method   
        !! first if 
        if(iter==1) then   !! first iter step, do SD
            !! search direction 
            print*, 'steepest descent for iter 1 '
            call SD(g_new, size(g_new), p_new)
            BFGS_step = 1

        else !! not the first iter step, try L_BFGS
            print*, 'L-BFGS for iter ',iter

            ! additional file needed: BFGS_step, m_old, g_old
            write(filename,'(a)') trim(directory)//'/optimizer/BFGS_step.dat'
            OPEN(IIN,FILE=filename,STATUS='old',action='read',iostat=ier)
            if(ier>0) then
                print*,'Error opening BFGS_step.dat file : ',trim(filename)
                stop 'Error opening BFGS_step.dat file : '
            else
                read(IIN,*) BFGS_step
            end if
            close(IIN)
            print*,'LOAD old BFGS_step -- ',trim(filename)

            !! m_old -->  m_new -m_old
            write(filename,'(a)')  trim(directory)//'/optimizer/m_old.bin'
            print*,'LOAD m_old -- ', trim(filename)
            open(unit=IIN,file=trim(filename),status='old',action='read',form='unformatted',iostat=ier)
            if (ier /= 0) then
                print*, 'Error: could not open m_old file: ',trim(filename)
                stop 'Error: could not open m_old file: '
            endif
            read(IIN) m_old
            close(IIN)
            if(DISPLAY_DETAILS .and. myrank==0) then
                print *,' Min / Max m_old = ', &
                    minval(m_old(:)),maxval(m_old(:))
            endif

            !! g_old --> g_new - g_old
            write(filename,'(a)')  trim(directory)//'/optimizer/g_old.bin'
            print*,'LOAD g_old -- ', trim(filename)
            open(unit=IIN,file=trim(filename),status='old',action='read',form='unformatted',iostat=ier)
            if (ier /= 0) then
                print*, 'Error: could not open g_old file: ',trim(filename)
                stop 'Error: could not open g_old file: '
            endif
            read(IIN) g_old
            close(IIN)

            !! m- step L-BFGS (accumulative steps and max steps)
            m = min(BFGS_step,BFGS_stepmax)

            !! Deltam, Deltag: renew
            Deltam(:,1) = m_new(:)-m_old(:)
            Deltag(:,1) = g_new(:)-g_old(:)

            if(m>1 .and. m<=BFGS_stepmax) then !! consider multiple previous steps
                ! additonal files: old Deltam, Deltag
                !! LOAD Deltam
                write(filename,'(a)')  trim(directory)//'/optimizer/Deltam.bin'
                print*,'LOAD Deltam -- ', trim(filename)
                open(unit=IIN,file=trim(filename),status='old',action='read',form='unformatted',iostat=ier)
                if (ier /= 0) then
                    print*, 'Error: could not open Deltam file: ',trim(filename)
                    stop 'Error: could not open Deltam file: '
                endif
                read(IIN) Deltam(:,2:BFGS_stepmax)
                close(IIN)

                !! LOAD Deltag
                write(filename,'(a)')  trim(directory)//'/optimizer/Deltag.bin'
                print*,'LOAD Deltag -- ', trim(filename)
                open(unit=IIN,file=trim(filename),status='old',action='read',form='unformatted',iostat=ier)
                if (ier /= 0) then
                    print*, 'Error: could not open Deltam file: ',trim(filename)
                    stop 'Error: could not open Deltam file: '
                endif
                read(IIN) Deltag(:,2:BFGS_stepmax)
                close(IIN)
            endif ! m>1

            ! BFGS direction
            print*, 'L-BFGS direction for iter=',iter, &
                ' BFGS_step=', BFGS_step, ' m=',m

            ! B-BFGS
            call LBFGS(Deltam, Deltag, g_new, size(g_new), m, BFGS_step, p_new)

            !! check restarting or not 
            if(BFGS_step ==1) then
                Deltag(:,:)=0.0
                Deltam(:,:)=0.0
            endif
        endif !iter==1

        if(BFGS_step > 1) then
            !! SAVE Deltam
            write(filename,'(a)')  trim(directory)//'/optimizer/Deltam.bin'
            print*,'SAVE Deltam -- ', trim(filename)
            open(unit=IOUT,file=trim(filename),status='unknown',form='unformatted',iostat=ier)
            if (ier /= 0) then
                print*, 'Error: could not open Deltam file: ',trim(filename)
                stop 'Error: could not open Deltam file: '
            endif
            write(IOUT) Deltam
            close(IOUT)

            !! SAVE Deltag
            write(filename,'(a)')  trim(directory)//'/optimizer/Deltag.bin'
            if(DISPLAY_DETAILS .and. myrank==0) print*,'SAVE Deltag -- ', trim(filename)
            open(unit=IOUT,file=trim(filename),status='unknown',form='unformatted',iostat=ier)
            if (ier /= 0) then
                print*, 'Error: could not open Deltam file: ',trim(filename)
                stop 'Error: could not open Deltam file: '
            endif
            write(IOUT) Deltag
            close(IOUT)
        endif

        !! save BFGS_step
        write(filename,'(a)') trim(directory)//'/optimizer/BFGS_step.dat'
        OPEN(IOUT,FILE=trim(filename),STATUS='unknown',iostat=ier)
        if(ier>0) then
            print*, 'Error: could not open cgstep file:',trim(filename)
            stop 'Error: could not open cgstep file:'
        else
            write(IOUT,*) BFGS_step
        end if
        close(IOUT)
        print*,'SAVE BFGS_step -- ',trim(filename)

    case default
        print*, 'opt_scheme must be among "SD"/"CG"/"QN" ...'
        stop 'opt_scheme must be among "SD"/"CG"/"QN" ...'

    end select      

    if(DISPLAY_DETAILS) then
        print*,'Min / Max direction = ',&
            minval(p_new(:)),maxval(p_new(:))
    endif

    !! noramlize p_new
    pmax=maxval(abs(p_new(:)))
    if(pmax>SMALL_VAL)  p_new = p_new / pmax

    deallocate(g_old)
    deallocate(p_old)
    deallocate(Deltam)
    deallocate(Deltag)

end subroutine update_direction
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine finalize(directory)
    use seismo_parameters
    implicit none
    integer :: ier
    character(len=MAX_FILENAME_LEN) :: filename
    character(len=MAX_STRING_LEN) :: directory

    !! SAVE gradient
    write(filename,'(a)')  trim(directory)//'/optimizer/g_new.bin'
    print*,'SAVE g_new -- ', trim(filename)
    open(unit=IOUT,file=trim(filename),status='unknown',form='unformatted',iostat=ier)
    if (ier /= 0) then
        print*, 'Error: could not open gradient file: ',trim(filename)
        stop 'Error: could not open gradient file: '
    endif
    write(IOUT) g_new
    close(IOUT) 

    if(DISPLAY_DETAILS) then
        print *,'Min / Max gradient = ',&
            minval(g_new(:)),maxval(g_new(:))
    endif

    !! SAVE direction
    write(filename,'(a)')  trim(directory)//'/optimizer/p_new.bin'
    print*,'SAVE p_new -- ',trim(filename) 
    open(unit=IOUT,file=trim(filename),status='unknown',form='unformatted',iostat=ier)
    if (ier /= 0) then
        print*, 'Error: could not open direction file: ',trim(filename)
        stop 'Error: could not open direction file: '
    endif
    write(IOUT) p_new
    close(IOUT) 

    if(DISPLAY_DETAILS) then
        print *,'Min / Max direction = ', &
            minval(p_new(:)),maxval(p_new(:))
    endif

    !! SAVE model 
    write(filename,'(a)')  trim(directory)//'/optimizer/m_new.bin'
    print*,'SAVE m_new -- ',trim(filename)
    open(unit=IOUT,file=trim(filename),status='unknown',form='unformatted',iostat=ier)
    if (ier /= 0) then
        print*, 'Error: could not open model file: ',trim(filename)
        stop 'Error: could not open model file: '
    endif
    write(IOUT) m_new
    close(IOUT)

    if(DISPLAY_DETAILS) then
        print *,'Min / Max model = ', &
            minval(m_new(:)),maxval(m_new(:))
    endif

    ! deallocate
    deallocate(g_new)
    deallocate(p_new)
    deallocate(m_new)

end subroutine finalize

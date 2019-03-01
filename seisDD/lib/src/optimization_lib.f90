!! main subroutines for optimization scheme
!! created by Yanhua O. Yuan ( yanhuay@princeton.edu)

!!!!!!!!!!!!!!!! OPTIMIZATION !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine SD(g_new, NSTEP, p_new)
    use constants
    implicit none

    !! steepest descent method
    real(kind=CUSTOM_REAL), dimension(*), intent(in) :: g_new
    integer, intent(in) :: NSTEP
    real(kind=CUSTOM_REAL), dimension(*), intent(out) :: p_new

    !! initialization 
    p_new (1:NSTEP) = 0.0_CUSTOM_REAL 

    !! SD 
    p_new(1:NSTEP) = - g_new(1:NSTEP)

end subroutine SD
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine NLCG(g_new, g_old, p_old, NSTEP, CG_method, cgstep, p_new)
    !! non-linear conjugate method
    !! NOTES :: The biggest threat to the CG method is a loss of conjugacy in the search vectors 
    !! one way is to periodically restart the CG method
    use constants
    implicit none

    real(kind=CUSTOM_REAL), dimension(*), intent(in) :: g_new, g_old, p_old
    integer, intent(in) :: NSTEP
    character(len=2), intent(in) :: CG_method
    integer, intent(out) :: cgstep
    real(kind=CUSTOM_REAL), dimension(*), intent(out) :: p_new

    real(kind=CUSTOM_REAL) :: beta
    real(kind=CUSTOM_REAL) :: top, bot  

    !! initialization 
    p_new (1:NSTEP) = 0.0
    cgstep = cgstep+1

    !! beta formula  
    select case (CG_method)
    case ("FR") ! Fletcher-Reeves 
        top = sum(g_new(1:NSTEP) * g_new(1:NSTEP))
        bot = sum(g_old(1:NSTEP) * g_old(1:NSTEP))
        beta = top / bot
    case ("PR") ! Polak-Ribiere 
        top = sum(g_new(1:NSTEP) * (g_new(1:NSTEP)-g_old(1:NSTEP)))
        bot = sum(g_old(1:NSTEP) * g_old(1:NSTEP)) 
        beta = top / bot
    case ("HS") ! Hestenes-Stiefel 
        top = sum(g_new(1:NSTEP) * (g_new(1:NSTEP)-g_old(1:NSTEP)))
        bot = -1.0 * sum(p_old(1:NSTEP) * (g_new(1:NSTEP)-g_old(1:NSTEP)))
        beta = top / bot
    case ("DY") ! Dai-Yuan: 
        top = sum(g_new(1:NSTEP) * g_new(1:NSTEP))
        bot = -1.0 * sum(p_old(1:NSTEP) * (g_new(1:NSTEP)-g_old(1:NSTEP)))
        beta = top / bot
    case default
        print*, 'CG_method must be among "FR"/"PR"/"HS"/"DY" ...';
        stop
    end select

    ! handle loss of conjugacy that results from the non-quadratic terms 
    if(beta<=0.0) then 
        print*, 'restarting NLCG ... [negative beta]'
        beta = 0.0  
        cgstep = 1
    elseif (dot_product(p_new(1:NSTEP), g_new(1:NSTEP)) > 0.0) then 
        print*, 'restarting NLCG ... [not a descent direction]' 
        beta = 0.0 
        cgstep = 1
    elseif (abs(dot_product(g_new(1:NSTEP), g_old(1:NSTEP)) &
            / dot_product(g_new(1:NSTEP), g_new(1:NSTEP))) > CG_threshold ) then 
        print*, 'restarting NLCG ... [loss of conjugacy]'
        beta = 0.0
        cgstep = 1
    endif  !! restart NLCG  

    !! search direction 
    p_new(1:NSTEP) = - g_new(1:NSTEP) + beta * p_old(1:NSTEP)


end subroutine NLCG
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine LBFGS(Dm, Dg, g_new, NSTEP, m, BFGS_step, p_new)
    !! Quasi-Newton method (L-BFGS)
    use constants
    implicit none

    integer, intent(in) :: NSTEP, m
    real(kind=CUSTOM_REAL), intent(in) :: Dm(NSTEP,m),Dg(NSTEP,m)
    real(kind=CUSTOM_REAL), dimension(*), intent(in) :: g_new
    integer, intent(out) :: BFGS_step ! for next iter
    real(kind=CUSTOM_REAL), dimension(*), intent(out) :: p_new

    real(kind=CUSTOM_REAL), dimension(m) :: alpha, beta, rho 
    integer :: i 
    real(kind=CUSTOM_REAL), dimension(NSTEP) :: q, z  
    real(kind=CUSTOM_REAL) :: invH

    !! initialization 
    p_new (1:NSTEP) = 0.0

    !!  initialization 
    q(1:NSTEP) = g_new(1:NSTEP) 

    !! the first loop
    do i = 1, m 
    rho(i) = 1.0/dot_product(Dg(1:NSTEP,i),Dm(1:NSTEP,i))
    alpha(i) = rho(i) * dot_product(Dm(1:NSTEP,i), q(1:NSTEP))
    q(1:NSTEP) = q(1:NSTEP) - alpha(i) * Dg(1:NSTEP,i)
    enddo

    invH = dot_product(Dg(1:NSTEP,1), Dm(1:NSTEP,1)) / &
        dot_product(Dg(1:NSTEP,1), Dg(1:NSTEP,1))
    z(1:NSTEP)  =  invH * q(1:NSTEP)

    !! the second loop 
    do i =  m, 1, -1 
    beta(i) = rho(i) * dot_product(Dg(1:NSTEP,i), z(1:NSTEP))
    z(1:NSTEP) = z(1:NSTEP) + (alpha(i) -beta(i)) * Dm(1:NSTEP,i)
    enddo

    !! restart L-BFGS
    if ( dot_product(g_new(1:NSTEP), -z(1:NSTEP)) &
        / dot_product(g_new(1:NSTEP), g_new(1:NSTEP)) > 0.0 ) then
        print*, 'restarting L-BFGS ... [not the descent direction]'
        p_new(1:NSTEP) = - g_new(1:NSTEP)
        BFGS_step = 1
    else
        !! L-BFGS search direction 
        p_new(1:NSTEP) = - z(1:NSTEP)
        BFGS_step = BFGS_step+1
    endif

end subroutine LBFGS

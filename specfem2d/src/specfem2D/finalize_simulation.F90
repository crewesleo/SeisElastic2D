
subroutine finalize_simulation()

#ifdef USE_MPI
  use mpi
#endif

  use specfem_par

  implicit none

integer i,ispec,j,iglob

#ifdef USE_MPI
  include "precision.h"
#endif


  real(kind=4),dimension(:,:,:),allocatable :: rhop_save, vp_save, vs_save, x_save, z_save
  real(kind=4),dimension(:,:,:),allocatable :: QKappa_save, Qmu_save !YY
  real(kind=4),dimension(:,:,:),allocatable :: rho_save, kappa_save, mu_save
  real(kind=4),dimension(:,:,:),allocatable :: c11_save, c13_save, c15_save, c33_save, c35_save, c55_save ! PWY
  real(kind=4),dimension(:,:,:),allocatable :: c12_save, c23_save, c25_save ! PWY
  real(kind=4),dimension(:,:,:),allocatable :: hti_thom_vp_save, hti_thom_vs_save, hti_thom_epsilon_save
  real(kind=4),dimension(:,:,:),allocatable :: hti_thom_delta_save, hti_thom_rhop_save

  real(kind=4),dimension(:,:,:),allocatable :: hti_thom2_vp_save, hti_thom2_vs_save, hti_thom2_epsilon_save
  real(kind=4),dimension(:,:,:),allocatable :: hti_thom2_eta_save, hti_thom2_rhop_save, hti_thom2_delta_save

  real(kind=4),dimension(:,:,:),allocatable :: hti_vel_vp_save, hti_vel_vs_save, hti_vel_vph_save
  real(kind=4),dimension(:,:,:),allocatable :: hti_vel_vpn_save, hti_vel_rhop_save

  real(kind=4),dimension(:,:,:),allocatable :: lambda_save,sigma_save
  real(kind=4),dimension(:,:,:),allocatable :: ip_save, is_save
  real(kind=4),dimension(:,:,:),allocatable :: theta_save !!! TTI titled angle

  real(kind=4),dimension(:,:,:),allocatable :: tti_thom_vp_save, tti_thom_vs_save, tti_thom_epsilon_save
  real(kind=4),dimension(:,:,:),allocatable :: tti_thom_delta_save, tti_thom_theta_save, tti_thom_rhop_save


if ( trim(SAVE_MODEL) /= 'default' ) then
   allocate(rhop_save(NGLLX,NGLLZ,nspec))
   allocate(vp_save(NGLLX,NGLLZ,nspec))
   allocate(vs_save(NGLLX,NGLLZ,nspec))
   allocate(rho_save(NGLLX,NGLLZ,nspec))
   allocate(kappa_save(NGLLX,NGLLZ,nspec))
   allocate(mu_save(NGLLX,NGLLZ,nspec))
   allocate(QKappa_save(NGLLX,NGLLZ,nspec)) !YY
   allocate(Qmu_save(NGLLX,NGLLZ,nspec)) !YY
   allocate(x_save(NGLLX,NGLLZ,nspec))
   allocate(z_save(NGLLX,NGLLZ,nspec))
   allocate(ip_save(NGLLX,NGLLZ,nspec))
   allocate(is_save(NGLLX,NGLLZ,nspec))
   allocate(lambda_save(NGLLX,NGLLZ,nspec))
   allocate(sigma_save(NGLLX,NGLLZ,nspec))
   allocate(c11_save(NGLLX,NGLLZ,nspec))
   allocate(c13_save(NGLLX,NGLLZ,nspec))
   allocate(c15_save(NGLLX,NGLLZ,nspec))
   allocate(c33_save(NGLLX,NGLLZ,nspec))
   allocate(c35_save(NGLLX,NGLLZ,nspec))
   allocate(c55_save(NGLLX,NGLLZ,nspec))
   allocate(c12_save(NGLLX,NGLLZ,nspec))
   allocate(c23_save(NGLLX,NGLLZ,nspec))
   allocate(c25_save(NGLLX,NGLLZ,nspec))
   allocate(theta_save(NGLLX,NGLLZ,nspec))
   allocate(hti_thom_rhop_save(NGLLX,NGLLZ,nspec))
   allocate(hti_thom_vp_save(NGLLX,NGLLZ,nspec))
   allocate(hti_thom_vs_save(NGLLX,NGLLZ,nspec))
   allocate(hti_thom_epsilon_save(NGLLX,NGLLZ,nspec))
   allocate(hti_thom_delta_save(NGLLX,NGLLZ,nspec))
   allocate(hti_thom2_rhop_save(NGLLX,NGLLZ,nspec))
   allocate(hti_thom2_vp_save(NGLLX,NGLLZ,nspec))
   allocate(hti_thom2_vs_save(NGLLX,NGLLZ,nspec))
   allocate(hti_thom2_epsilon_save(NGLLX,NGLLZ,nspec))
   allocate(hti_thom2_eta_save(NGLLX,NGLLZ,nspec))
   allocate(hti_thom2_delta_save(NGLLX,NGLLZ,nspec))
   allocate(hti_vel_rhop_save(NGLLX,NGLLZ,nspec))
   allocate(hti_vel_vp_save(NGLLX,NGLLZ,nspec))
   allocate(hti_vel_vs_save(NGLLX,NGLLZ,nspec))
   allocate(hti_vel_vph_save(NGLLX,NGLLZ,nspec))
   allocate(hti_vel_vpn_save(NGLLX,NGLLZ,nspec))
   allocate(tti_thom_rhop_save(NGLLX,NGLLZ,nspec))
   allocate(tti_thom_vp_save(NGLLX,NGLLZ,nspec))
   allocate(tti_thom_vs_save(NGLLX,NGLLZ,nspec))
   allocate(tti_thom_epsilon_save(NGLLX,NGLLZ,nspec))
   allocate(tti_thom_delta_save(NGLLX,NGLLZ,nspec))
   allocate(tti_thom_theta_save(NGLLX,NGLLZ,nspec))

do ispec=1,nspec
          do j = 1,NGLLZ
              do i = 1,NGLLX

              QKappa_save(i,j,ispec) = QKappastore(i,j,ispec) !YY
              Qmu_save(i,j,ispec) = Qmustore(i,j,ispec) !YY
              !rho_save(i,j,ispec)            = density(1,kmato(ispec)) !YY
              rhop_save(i,j,ispec)            = rhostore(i,j,ispec)  !YY
             ! lambdal_unrelaxed_elastic      = poroelastcoef(1,1,kmato(ispec))
             ! mul_unrelaxed_elastic          = poroelastcoef(2,1,kmato(ispec))
             ! kappa_save(i,j,ispec)          = lambdal_unrelaxed_elastic + TWO*mul_unrelaxed_elastic/3._CUSTOM_REAL
             ! vp_save(i,j,ispec)             = sqrt((kappa_save(i,j,ispec) + &
             !                                   4._CUSTOM_REAL*mul_unrelaxed_elastic/ &
             !                                   3._CUSTOM_REAL)/density(1,kmato(ispec)))
             ! vs_save(i,j,ispec)             = sqrt(mul_unrelaxed_elastic/density(1,kmato(ispec)))
              vp_save(i,j,ispec)             = sqrt((kappastore(i,j,ispec) + &
                  4._CUSTOM_REAL*mustore(i,j,ispec)/ &
                  3._CUSTOM_REAL)/rhostore(i,j,ispec))
              vs_save(i,j,ispec)             = sqrt(mustore(i,j,ispec)/rhostore(i,j,ispec))
              rho_save(i,j,ispec)= rhostore(i,j,ispec)
              kappa_save(i,j,ispec)= kappastore(i,j,ispec)
              mu_save(i,j,ispec)= mustore(i,j,ispec)

              c11_save(i,j,ispec) = c11store(i,j,ispec) ! PWY
              c13_save(i,j,ispec) = c13store(i,j,ispec) ! PWY
              c15_save(i,j,ispec) = c15store(i,j,ispec) ! PWY
              c33_save(i,j,ispec) = c33store(i,j,ispec) ! PWY
              c35_save(i,j,ispec) = c35store(i,j,ispec) ! PWY
              c55_save(i,j,ispec) = c55store(i,j,ispec) ! PWY
              c12_save(i,j,ispec) = c12store(i,j,ispec) ! PWY
              c23_save(i,j,ispec) = c23store(i,j,ispec) ! PWY
              c25_save(i,j,ispec) = c25store(i,j,ispec) ! PWY
              
              if (ANISO .and. (trim(M_PAR)=='ttiecu')) then
                !! when using ttiecu model parameterizations, the input elastic
                !constant parameters must be unrotataed elastic constants
                theta_save(i,j,ispec) = 0.5 * ATAN(2*(c15_save(i,j,ispec)+c35_save(i,j,ispec))/&
                   (c11_save(i,j,ispec)-c33_save(i,j,ispec)))
              endif

              if (ANISO .and. (trim(M_PAR)=='ttithom')) then
                tti_thom_rhop_save(i,j,ispec) = rhostore(i,j,ispec)
                tti_thom_vp_save(i,j,ispec) = SQRT(c33store(i,j,ispec)/rhostore(i,j,ispec))
                tti_thom_vs_save(i,j,ispec) = SQRT(c55store(i,j,ispec)/rhostore(i,j,ispec))
                tti_thom_epsilon_save(i,j,ispec) = (c11store(i,j,ispec)-c33store(i,j,ispec))/&
                  (2._CUSTOM_REAL*c33store(i,j,ispec))

                tti_thom_delta_save(i,j,ispec) = (c13store(i,j,ispec)+c55store(i,j,ispec))**2/&
                  (2._CUSTOM_REAL*c33store(i,j,ispec)*(c33store(i,j,ispec)-c55store(i,j,ispec)))-&
                  (c33store(i,j,ispec)-c55store(i,j,ispec))/(2._CUSTOM_REAL*c33store(i,j,ispec))

                tti_thom_theta_save(i,j,ispec) = 0.5 * ATAN(2*(c15_save(i,j,ispec)+c35_save(i,j,ispec))/&
                   (c11_save(i,j,ispec)-c33_save(i,j,ispec)))
              endif
              lambda_save(i,j,ispec) = rhostore(i,j,ispec)*vp_save(i,j,ispec)**2 - &
                  2._CUSTOM_REAL*mustore(i,j,ispec)
              ip_save(i,j,ispec) = rhostore(i,j,ispec)*vp_save(i,j,ispec)
              is_save(i,j,ispec) = rhostore(i,j,ispec)*vs_save(i,j,ispec)
              sigma_save(i,j,ispec) = vp_save(i,j,ispec)/vs_save(i,j,ispec)

              if (ANISO .and. (trim(M_PAR)=='htithom' .OR. trim(M_PAR)=='vtithom')) then
                hti_thom_rhop_save(i,j,ispec) = rhostore(i,j,ispec)
                hti_thom_vp_save(i,j,ispec) = SQRT(c33store(i,j,ispec)/rhostore(i,j,ispec))
                hti_thom_vs_save(i,j,ispec) = SQRT(c55store(i,j,ispec)/rhostore(i,j,ispec))
                hti_thom_epsilon_save(i,j,ispec) = (c11store(i,j,ispec)-c33store(i,j,ispec))/&
                  (2._CUSTOM_REAL*c33store(i,j,ispec))
                hti_thom_delta_save(i,j,ispec) = (c13store(i,j,ispec)+c55store(i,j,ispec))**2/&
                  (2._CUSTOM_REAL*c33store(i,j,ispec)*(c33store(i,j,ispec)-c55store(i,j,ispec)))-&
                  (c33store(i,j,ispec)-c55store(i,j,ispec))/(2._CUSTOM_REAL*c33store(i,j,ispec))
              endif
              if (ANISO .and. (trim(M_PAR)=='htithom2' .OR. trim(M_PAR)=='vtithom2')) then
                hti_thom2_rhop_save(i,j,ispec) = rhostore(i,j,ispec)
                hti_thom2_vp_save(i,j,ispec) = SQRT(c11store(i,j,ispec)/rhostore(i,j,ispec))
                hti_thom2_vs_save(i,j,ispec) = SQRT(c55store(i,j,ispec)/rhostore(i,j,ispec))
                hti_thom2_epsilon_save(i,j,ispec) = (c11store(i,j,ispec)-c33store(i,j,ispec))/&
                  (2._CUSTOM_REAL*c33store(i,j,ispec))
                hti_thom2_delta_save(i,j,ispec) = (c13store(i,j,ispec)+c55store(i,j,ispec))**2/&
                  (2._CUSTOM_REAL*c33store(i,j,ispec)*(c33store(i,j,ispec)-c55store(i,j,ispec)))-&
                  (c33store(i,j,ispec)-c55store(i,j,ispec))/(2._CUSTOM_REAL*c33store(i,j,ispec))
                hti_thom2_eta_save(i,j,ispec) = (hti_thom2_epsilon_save(i,j,ispec)-hti_thom2_delta_save(i,j,ispec))/&
                  (1._CUSTOM_REAL + 2._CUSTOM_REAL * hti_thom2_delta_save(i,j,ispec))
              endif
              if (ANISO .and. (trim(M_PAR)=='htivel' .OR. trim(M_PAR)=='vtivel')) then
                hti_vel_rhop_save(i,j,ispec) = rhostore(i,j,ispec)
                hti_vel_vp_save(i,j,ispec) = SQRT(c33store(i,j,ispec)/rhostore(i,j,ispec))
                hti_vel_vs_save(i,j,ispec) = SQRT(c55store(i,j,ispec)/rhostore(i,j,ispec))
                hti_thom_epsilon_save(i,j,ispec) = (c11store(i,j,ispec)-c33store(i,j,ispec))/&
                  (2._CUSTOM_REAL*c33store(i,j,ispec))
                hti_thom_delta_save(i,j,ispec) = (c13store(i,j,ispec)+c55store(i,j,ispec))**2/&
                  (2._CUSTOM_REAL*c33store(i,j,ispec)*(c33store(i,j,ispec)-c55store(i,j,ispec)))-&
                  (c33store(i,j,ispec)-c55store(i,j,ispec))/(2._CUSTOM_REAL*c33store(i,j,ispec))
                hti_vel_vph_save(i,j,ispec) = hti_vel_vp_save(i,j,ispec) * &
                  SQRT(1._CUSTOM_REAL + 2._CUSTOM_REAL * hti_thom_epsilon_save(i,j,ispec))
                hti_vel_vpn_save(i,j,ispec) = hti_vel_vp_save(i,j,ispec) * &
                  SQRT(1._CUSTOM_REAL + 2._CUSTOM_REAL * hti_thom_delta_save(i,j,ispec))
              endif
              
              iglob = ibool(i,j,ispec)
              x_save(i,j,ispec)              = coord(1,iglob)
              z_save(i,j,ispec)              = coord(2,iglob)
              enddo
        enddo
enddo

!print*, 'YY, min/max of rho_save -- ',minval(rho_save(:,:,:)),maxval(rho_save(:,:,:))
!print*, 'YY, min/max of vp_save --',minval(vp_save(:,:,:)),maxval(vp_save(:,:,:))
!print*, 'YY, min/max of vs_save --',minval(vs_save(:,:,:)),maxval(vs_save(:,:,:))

  if(trim(SAVE_MODEL) == 'legacy') then

    write(inputname,'(a,i6.6,a)') 'DATA/proc',myrank,'_model_velocity.dat_input'
    open(unit=1001,file=inputname,status='unknown')
    do ispec = 1,nspec
      do j = 1,NGLLZ
        do i = 1,NGLLX
          iglob = ibool(i,j,ispec)
          write(1001,'(6e15.5e4)') x_save(i,j,ispec), x_save(i,j,ispec),z_save(i,j,ispec),rho_save(i,j,ispec),&
                                   vp_save(i,j,ispec),vs_save(i,j,ispec)
        enddo
      enddo
    enddo
    close(1001)


  else if(trim(SAVE_MODEL)=='ascii') then

    write(inputname,'(a,i6.6,a)') 'DATA/proc',myrank,'_rho_vp_vs.dat'
    open(unit=1001,file= inputname,status='unknown')
    do ispec = 1,nspec
      do j = 1,NGLLZ
        do i = 1,NGLLX
          iglob = ibool(i,j,ispec)
          write(1001,'(5e15.5e4)') x_save(i,j,ispec),z_save(i,j,ispec),rho_save(i,j,ispec),vp_save(i,j,ispec),vs_save(i,j,ispec)
        enddo
      enddo
    enddo
    close(1001)

  else if((trim(SAVE_MODEL) == 'binary') .or. (trim(SAVE_MODEL) == 'gll')) then

          write(outputname,'(a,i6.6,a)') 'DATA/proc',myrank,'_rho.bin'
          open(unit=172,file=outputname,status='unknown',form='unformatted')
          write(172) rho_save
          close(172)
        if (.not. ANISO .and. (trim(M_PAR)=='isodv')) then
          write(outputname,'(a,i6.6,a)') 'DATA/proc',myrank,'_vp.bin'
          open(unit=172,file=outputname,status='unknown',form='unformatted')
          write(172) vp_save
          close(172)
          write(outputname,'(a,i6.6,a)') 'DATA/proc',myrank,'_vs.bin'
          open(unit=172,file=outputname,status='unknown',form='unformatted')
          write(172) vs_save
          close(172)
          !YY 
          write(outputname,'(a,i6.6,a)') 'DATA/proc',myrank,'_QKappa.bin'
          open(unit=172,file=outputname,status='unknown',form='unformatted')
          write(172) QKappa_save
          close(172)
          !YY 
          write(outputname,'(a,i6.6,a)') 'DATA/proc',myrank,'_Qmu.bin'
          open(unit=172,file=outputname,status='unknown',form='unformatted')
          write(172) Qmu_save
          close(172)
          ! PWY
        endif
        if (.not. ANISO .and. (trim(M_PAR)=='isodm')) then
          write(outputname,'(a,i6.6,a)') 'DATA/proc',myrank,'_kappa.bin'
          open(unit=172,file=outputname,status='unknown',form='unformatted')
          write(172) kappa_save
          close(172)
          write(outputname,'(a,i6.6,a)') 'DATA/proc',myrank,'_mu.bin'
          open(unit=172,file=outputname,status='unknown',form='unformatted')
          write(172) mu_save
          close(172)
          !YY 
          write(outputname,'(a,i6.6,a)') 'DATA/proc',myrank,'_QKappa.bin'
          open(unit=172,file=outputname,status='unknown',form='unformatted')
          write(172) QKappa_save
          close(172)
          !YY 
          write(outputname,'(a,i6.6,a)') 'DATA/proc',myrank,'_Qmu.bin'
          open(unit=172,file=outputname,status='unknown',form='unformatted')
          write(172) Qmu_save
          close(172)
          ! PWY
        endif

        if (.not. ANISO .and. (trim(M_PAR)=='isodl')) then
          write(outputname,'(a,i6.6,a)') 'DATA/proc',myrank,'_lambda.bin'
          open(unit=172,file=outputname,status='unknown',form='unformatted')
          write(172) lambda_save
          close(172)
          write(outputname,'(a,i6.6,a)') 'DATA/proc',myrank,'_mu.bin'
          open(unit=172,file=outputname,status='unknown',form='unformatted')
          write(172) mu_save
          close(172)
          !YY 
          write(outputname,'(a,i6.6,a)') 'DATA/proc',myrank,'_QKappa.bin'
          open(unit=172,file=outputname,status='unknown',form='unformatted')
          write(172) QKappa_save
          close(172)
          !YY 
          write(outputname,'(a,i6.6,a)') 'DATA/proc',myrank,'_Qmu.bin'
          open(unit=172,file=outputname,status='unknown',form='unformatted')
          write(172) Qmu_save
          close(172)
          ! PWY
        endif

        if (.not. ANISO .and. (trim(M_PAR)=='isodip')) then
          write(outputname,'(a,i6.6,a)') 'DATA/proc',myrank,'_ip.bin'
          open(unit=172,file=outputname,status='unknown',form='unformatted')
          write(172) ip_save
          close(172)
          write(outputname,'(a,i6.6,a)') 'DATA/proc',myrank,'_is.bin'
          open(unit=172,file=outputname,status='unknown',form='unformatted')
          write(172) is_save
          close(172)
          !YY 
          write(outputname,'(a,i6.6,a)') 'DATA/proc',myrank,'_QKappa.bin'
          open(unit=172,file=outputname,status='unknown',form='unformatted')
          write(172) QKappa_save
          close(172)
          !YY 
          write(outputname,'(a,i6.6,a)') 'DATA/proc',myrank,'_Qmu.bin'
          open(unit=172,file=outputname,status='unknown',form='unformatted')
          write(172) Qmu_save
          close(172)
          ! PWY
        endif

        if (.not. ANISO .and. (trim(M_PAR)=='isovipi')) then
          write(outputname,'(a,i6.6,a)') 'DATA/proc',myrank,'_vp.bin'
          open(unit=172,file=outputname,status='unknown',form='unformatted')
          write(172) vp_save
          close(172)

          write(outputname,'(a,i6.6,a)') 'DATA/proc',myrank,'_vs.bin'
          open(unit=172,file=outputname,status='unknown',form='unformatted')
          write(172) vs_save
          close(172)

          write(outputname,'(a,i6.6,a)') 'DATA/proc',myrank,'_ip.bin'
          open(unit=172,file=outputname,status='unknown',form='unformatted')
          write(172) ip_save
          close(172)
          ! 
          write(outputname,'(a,i6.6,a)') 'DATA/proc',myrank,'_QKappa.bin'
          open(unit=172,file=outputname,status='unknown',form='unformatted')
          write(172) QKappa_save
          close(172)
          ! 
          write(outputname,'(a,i6.6,a)') 'DATA/proc',myrank,'_Qmu.bin'
          open(unit=172,file=outputname,status='unknown',form='unformatted')
          write(172) Qmu_save
          close(172)
          ! PWY
        endif

        if (.not. ANISO .and. (trim(M_PAR)=='isovipii')) then
          write(outputname,'(a,i6.6,a)') 'DATA/proc',myrank,'_vp.bin'
          open(unit=172,file=outputname,status='unknown',form='unformatted')
          write(172) vp_save
          close(172)

          write(outputname,'(a,i6.6,a)') 'DATA/proc',myrank,'_vs.bin'
          open(unit=172,file=outputname,status='unknown',form='unformatted')
          write(172) vs_save
          close(172)

          write(outputname,'(a,i6.6,a)') 'DATA/proc',myrank,'_is.bin'
          open(unit=172,file=outputname,status='unknown',form='unformatted')
          write(172) is_save
          close(172)
          ! 
          write(outputname,'(a,i6.6,a)') 'DATA/proc',myrank,'_QKappa.bin'
          open(unit=172,file=outputname,status='unknown',form='unformatted')
          write(172) QKappa_save
          close(172)
          ! 
          write(outputname,'(a,i6.6,a)') 'DATA/proc',myrank,'_Qmu.bin'
          open(unit=172,file=outputname,status='unknown',form='unformatted')
          write(172) Qmu_save
          close(172)
          ! PWY
        endif

        if (.not. ANISO .and. (trim(M_PAR)=='isosigma')) then
          write(outputname,'(a,i6.6,a)') 'DATA/proc',myrank,'_vp.bin'
          open(unit=172,file=outputname,status='unknown',form='unformatted')
          write(172) vp_save
          close(172)

          write(outputname,'(a,i6.6,a)') 'DATA/proc',myrank,'_sigma.bin'
          open(unit=172,file=outputname,status='unknown',form='unformatted')
          write(172) sigma_save
          close(172)
          ! 
          write(outputname,'(a,i6.6,a)') 'DATA/proc',myrank,'_QKappa.bin'
          open(unit=172,file=outputname,status='unknown',form='unformatted')
          write(172) QKappa_save
          close(172)
          ! 
          write(outputname,'(a,i6.6,a)') 'DATA/proc',myrank,'_Qmu.bin'
          open(unit=172,file=outputname,status='unknown',form='unformatted')
          write(172) Qmu_save
          close(172)
          ! PWY
        endif

        if (ANISO .and. (trim(M_PAR)=='ttiec')) then
          write(outputname,'(a,i6.6,a)') 'DATA/proc',myrank,'_c11.bin'
          open(unit=172,file=outputname,status='unknown',form='unformatted')
          write(172) c11_save
          close(172)
          ! PWY
          write(outputname,'(a,i6.6,a)') 'DATA/proc',myrank,'_c13.bin'
          open(unit=172,file=outputname,status='unknown',form='unformatted')
          write(172) c13_save
          close(172)
          ! PWY
          write(outputname,'(a,i6.6,a)') 'DATA/proc',myrank,'_c15.bin'
          open(unit=172,file=outputname,status='unknown',form='unformatted')
          write(172) c15_save
          close(172)
          ! PWY
          write(outputname,'(a,i6.6,a)') 'DATA/proc',myrank,'_c33.bin'
          open(unit=172,file=outputname,status='unknown',form='unformatted')
          write(172) c33_save
          close(172)
          ! PWY
          write(outputname,'(a,i6.6,a)') 'DATA/proc',myrank,'_c35.bin'
          open(unit=172,file=outputname,status='unknown',form='unformatted')
          write(172) c35_save
          close(172)
          ! PWY
          write(outputname,'(a,i6.6,a)') 'DATA/proc',myrank,'_c55.bin'
          open(unit=172,file=outputname,status='unknown',form='unformatted')
          write(172) c55_save
          close(172)
          ! PWY
          write(outputname,'(a,i6.6,a)') 'DATA/proc',myrank,'_c12.bin'
          open(unit=172,file=outputname,status='unknown',form='unformatted')
          write(172) c12_save
          close(172)
          ! PWY
          write(outputname,'(a,i6.6,a)') 'DATA/proc',myrank,'_c23.bin'
          open(unit=172,file=outputname,status='unknown',form='unformatted')
          write(172) c23_save
          close(172)
          ! PWY
          write(outputname,'(a,i6.6,a)') 'DATA/proc',myrank,'_c25.bin'
          open(unit=172,file=outputname,status='unknown',form='unformatted')
          write(172) c25_save
          close(172)
        endif
        if (ANISO .and. (trim(M_PAR)=='ttiecu')) then
          write(outputname,'(a,i6.6,a)') 'DATA/proc',myrank,'_c11.bin'
          open(unit=172,file=outputname,status='unknown',form='unformatted')
          write(172) c11_save
          close(172)
          ! PWY
          write(outputname,'(a,i6.6,a)') 'DATA/proc',myrank,'_c13.bin'
          open(unit=172,file=outputname,status='unknown',form='unformatted')
          write(172) c13_save
          close(172)
          ! PWY
          write(outputname,'(a,i6.6,a)') 'DATA/proc',myrank,'_c15.bin'
          open(unit=172,file=outputname,status='unknown',form='unformatted')
          write(172) c15_save
          close(172)
          ! PWY
          write(outputname,'(a,i6.6,a)') 'DATA/proc',myrank,'_c33.bin'
          open(unit=172,file=outputname,status='unknown',form='unformatted')
          write(172) c33_save
          close(172)
          ! PWY
          write(outputname,'(a,i6.6,a)') 'DATA/proc',myrank,'_c35.bin'
          open(unit=172,file=outputname,status='unknown',form='unformatted')
          write(172) c35_save
          close(172)
          ! PWY
          write(outputname,'(a,i6.6,a)') 'DATA/proc',myrank,'_c55.bin'
          open(unit=172,file=outputname,status='unknown',form='unformatted')
          write(172) c55_save
          close(172)
          ! PWY
          write(outputname,'(a,i6.6,a)') 'DATA/proc',myrank,'_c12.bin'
          open(unit=172,file=outputname,status='unknown',form='unformatted')
          write(172) c12_save
          close(172)
          ! PWY
          write(outputname,'(a,i6.6,a)') 'DATA/proc',myrank,'_c23.bin'
          open(unit=172,file=outputname,status='unknown',form='unformatted')
          write(172) c23_save
          close(172)
          ! PWY
          write(outputname,'(a,i6.6,a)') 'DATA/proc',myrank,'_c25.bin'
          open(unit=172,file=outputname,status='unknown',form='unformatted')
          write(172) c25_save
          close(172)

          ! PWY
          write(outputname,'(a,i6.6,a)') 'DATA/proc',myrank,'_theta.bin'
          open(unit=172,file=outputname,status='unknown',form='unformatted')
          write(172) theta_save
          close(172)
        endif
        if (ANISO .and. (trim(M_PAR)=='ttithom')) then
          write(outputname,'(a,i6.6,a)') 'DATA/proc',myrank,'_tti_thom_rhop.bin'
          open(unit=172,file=outputname,status='unknown',form='unformatted')
          write(172) tti_thom_rhop_save
          close(172)

          write(outputname,'(a,i6.6,a)') 'DATA/proc',myrank,'_tti_thom_vp.bin'
          open(unit=172,file=outputname,status='unknown',form='unformatted')
          write(172) tti_thom_vp_save
          close(172)

          write(outputname,'(a,i6.6,a)') 'DATA/proc',myrank,'_tti_thom_vs.bin'
          open(unit=172,file=outputname,status='unknown',form='unformatted')
          write(172) tti_thom_vs_save
          close(172)

          write(outputname,'(a,i6.6,a)') 'DATA/proc',myrank,'_tti_thom_epsilon.bin'
          open(unit=172,file=outputname,status='unknown',form='unformatted')
          write(172) tti_thom_epsilon_save
          close(172)

          write(outputname,'(a,i6.6,a)') 'DATA/proc',myrank,'_tti_thom_delta.bin'
          open(unit=172,file=outputname,status='unknown',form='unformatted')
          write(172) tti_thom_delta_save
          close(172)

          write(outputname,'(a,i6.6,a)') 'DATA/proc',myrank,'_theta.bin'
          open(unit=172,file=outputname,status='unknown',form='unformatted')
          write(172) tti_thom_theta_save
          close(172)
        endif
        if (ANISO .and. (trim(M_PAR)=='htiec' .OR. trim(M_PAR)=='vtiec')) then  
          write(outputname,'(a,i6.6,a)') 'DATA/proc',myrank,'_c11.bin'
          open(unit=172,file=outputname,status='unknown',form='unformatted')
          write(172) c11_save
          close(172)
          ! PWY
          write(outputname,'(a,i6.6,a)') 'DATA/proc',myrank,'_c13.bin'
          open(unit=172,file=outputname,status='unknown',form='unformatted')
          write(172) c13_save
          close(172)
          ! PWY
          write(outputname,'(a,i6.6,a)') 'DATA/proc',myrank,'_c15.bin'
          open(unit=172,file=outputname,status='unknown',form='unformatted')
          write(172) c15_save
          close(172)
          ! PWY
          write(outputname,'(a,i6.6,a)') 'DATA/proc',myrank,'_c33.bin'
          open(unit=172,file=outputname,status='unknown',form='unformatted')
          write(172) c33_save
          close(172)
          ! PWY
          write(outputname,'(a,i6.6,a)') 'DATA/proc',myrank,'_c35.bin'
          open(unit=172,file=outputname,status='unknown',form='unformatted')
          write(172) c35_save
          close(172)
          ! PWY
          write(outputname,'(a,i6.6,a)') 'DATA/proc',myrank,'_c55.bin'
          open(unit=172,file=outputname,status='unknown',form='unformatted')
          write(172) c55_save
          close(172)
          ! PWY
          write(outputname,'(a,i6.6,a)') 'DATA/proc',myrank,'_c12.bin'
          open(unit=172,file=outputname,status='unknown',form='unformatted')
          write(172) c12_save
          close(172)
          ! PWY
          write(outputname,'(a,i6.6,a)') 'DATA/proc',myrank,'_c23.bin'
          open(unit=172,file=outputname,status='unknown',form='unformatted')
          write(172) c23_save
          close(172)
          ! PWY
          write(outputname,'(a,i6.6,a)') 'DATA/proc',myrank,'_c25.bin'
          open(unit=172,file=outputname,status='unknown',form='unformatted')
          write(172) c25_save
          close(172)
        endif
        if (ANISO .and. (trim(M_PAR)=='htithom' .OR. trim(M_PAR)=='vtithom')) then
          write(outputname,'(a,i6.6,a)') 'DATA/proc',myrank,'_hti_thom_rhop.bin'
          open(unit=172,file=outputname,status='unknown',form='unformatted')
          write(172) hti_thom_rhop_save
          close(172)

          write(outputname,'(a,i6.6,a)') 'DATA/proc',myrank,'_hti_thom_vp.bin'
          open(unit=172,file=outputname,status='unknown',form='unformatted')
          write(172) hti_thom_vp_save
          close(172)

          write(outputname,'(a,i6.6,a)') 'DATA/proc',myrank,'_hti_thom_vs.bin'
          open(unit=172,file=outputname,status='unknown',form='unformatted')
          write(172) hti_thom_vs_save
          close(172)

          write(outputname,'(a,i6.6,a)') 'DATA/proc',myrank,'_hti_thom_epsilon.bin'
          open(unit=172,file=outputname,status='unknown',form='unformatted')
          write(172) hti_thom_epsilon_save
          close(172)

          write(outputname,'(a,i6.6,a)') 'DATA/proc',myrank,'_hti_thom_delta.bin'
          open(unit=172,file=outputname,status='unknown',form='unformatted')
          write(172) hti_thom_delta_save
          close(172)
        endif
        if (ANISO .and. (trim(M_PAR)=='htithom2' .OR. trim(M_PAR)=='vtithom2')) then
          write(outputname,'(a,i6.6,a)') 'DATA/proc',myrank,'_hti_thom2_rhop.bin'
          open(unit=172,file=outputname,status='unknown',form='unformatted')
          write(172) hti_thom2_rhop_save
          close(172)

          write(outputname,'(a,i6.6,a)') 'DATA/proc',myrank,'_hti_thom2_vp.bin'
          open(unit=172,file=outputname,status='unknown',form='unformatted')
          write(172) hti_thom2_vp_save
          close(172)

          write(outputname,'(a,i6.6,a)') 'DATA/proc',myrank,'_hti_thom2_vs.bin'
          open(unit=172,file=outputname,status='unknown',form='unformatted')
          write(172) hti_thom2_vs_save
          close(172)

          write(outputname,'(a,i6.6,a)') 'DATA/proc',myrank,'_hti_thom2_epsilon.bin'
          open(unit=172,file=outputname,status='unknown',form='unformatted')
          write(172) hti_thom2_epsilon_save
          close(172)

          write(outputname,'(a,i6.6,a)') 'DATA/proc',myrank,'_hti_thom2_delta.bin'
          open(unit=172,file=outputname,status='unknown',form='unformatted')
          write(172) hti_thom2_delta_save
          close(172)

          write(outputname,'(a,i6.6,a)') 'DATA/proc',myrank,'_hti_thom2_eta.bin'
          open(unit=172,file=outputname,status='unknown',form='unformatted')
          write(172) hti_thom2_eta_save
          close(172)
        endif
        if (ANISO .and. (trim(M_PAR)=='htivel' .OR. trim(M_PAR)=='vtivel')) then
          write(outputname,'(a,i6.6,a)') 'DATA/proc',myrank,'_hti_vel_rhop.bin'
          open(unit=172,file=outputname,status='unknown',form='unformatted')
          write(172) hti_vel_rhop_save
          close(172)

          write(outputname,'(a,i6.6,a)') 'DATA/proc',myrank,'_hti_vel_vp.bin'
          open(unit=172,file=outputname,status='unknown',form='unformatted')
          write(172) hti_vel_vp_save
          close(172)

          write(outputname,'(a,i6.6,a)') 'DATA/proc',myrank,'_hti_vel_vs.bin'
          open(unit=172,file=outputname,status='unknown',form='unformatted')
          write(172) hti_vel_vs_save
          close(172)

          write(outputname,'(a,i6.6,a)') 'DATA/proc',myrank,'_hti_vel_vph.bin'
          open(unit=172,file=outputname,status='unknown',form='unformatted')
          write(172) hti_vel_vph_save
          close(172)

          write(outputname,'(a,i6.6,a)') 'DATA/proc',myrank,'_hti_vel_vpn.bin'
          open(unit=172,file=outputname,status='unknown',form='unformatted')
          write(172) hti_vel_vpn_save
          close(172)
        endif        
          write(outputname,'(a,i6.6,a)') 'DATA/proc',myrank,'_x.bin'
          open(unit=172,file=outputname,status='unknown',form='unformatted')
          write(172) x_save
          close(172)
          write(outputname,'(a,i6.6,a)') 'DATA/proc',myrank,'_z.bin'
          open(unit=172,file=outputname,status='unknown',form='unformatted')
          write(172) z_save
          close(172)

          write(outputname,'(a,i6.6,a)') 'DATA/proc',myrank,'_jacobian.bin'
          open(unit=172,file=outputname,status='unknown',form='unformatted')
          write(172) jacobian
          close(172)


  else
       stop 'Save Model not implemented for external and tomo'

  endif !Type of model


endif !save model




if (GPU_MODE) call prepare_cleanup_device(Mesh_pointer, &
                              any_acoustic,any_elastic, &
                              STACEY_BOUNDARY_CONDITIONS, &
                              ANISOTROPY, &
                              APPROXIMATE_HESS_KL)


  if(output_wavefield_dumps) deallocate(mask_ibool)


!!!! Displacement Etienne GPU

! stores absorbing boundary contributions into files
      if(anyabs .and. SAVE_FORWARD .and. SIMULATION_TYPE == 1 .and. (.not. PML_BOUNDARY_CONDITIONS)) then

      if (any_acoustic) then

        !--- left absorbing boundary
        if(nspec_left >0) write(65) b_absorb_acoustic_left
        !--- right absorbing boundary
        if(nspec_right >0) write(66) b_absorb_acoustic_right
        !--- bottom absorbing boundary
        if(nspec_bottom >0) write(67) b_absorb_acoustic_bottom
        !--- top absorbing boundary
        if(nspec_top >0) write(68) b_absorb_acoustic_top

      endif !any acoustic

      close(65)
      close(66)
      close(67)
      close(68)
      close(72)

 if(any_elastic) then

        !--- left absorbing boundary
        if(nspec_left >0) write(35) b_absorb_elastic_left
        !--- right absorbing boundary
        if(nspec_right >0)  write(36) b_absorb_elastic_right
        !--- bottom absorbing boundary
        if(nspec_bottom >0)  write(37) b_absorb_elastic_bottom
        !--- top absorbing boundary
        if(nspec_top >0) write(38) b_absorb_elastic_top

   endif !any elastic

      close(35)
      close(36)
      close(37)
      close(38)
      close(71)


    if(any_poroelastic) then
      close(25)
      close(45)
      close(26)
      close(46)
      close(29)
      close(47)
      close(28)
      close(48)
    endif

  endif

!
!--- save last frame
!
  if(SAVE_FORWARD .and. SIMULATION_TYPE ==1 .and. any_elastic) then
    if ( myrank == 0 ) then
      write(IOUT,*)
      write(IOUT,*) 'Saving elastic last frame...'
      write(IOUT,*)
    endif
    write(outputname,'(a,i6.6,a)') 'lastframe_elastic',myrank,'.bin'
    open(unit=55,file='OUTPUT_FILES/'//outputname,status='unknown',form='unformatted')

        write(55) displ_elastic
        write(55) veloc_elastic
        write(55) accel_elastic

    close(55)
  endif

  if(SAVE_FORWARD .and. SIMULATION_TYPE ==1 .and. any_poroelastic) then
    if ( myrank == 0 ) then
      write(IOUT,*)
      write(IOUT,*) 'Saving poroelastic last frame...'
      write(IOUT,*)
    endif
    write(outputname,'(a,i6.6,a)') 'lastframe_poroelastic_s',myrank,'.bin'
    open(unit=55,file='OUTPUT_FILES/'//outputname,status='unknown',form='unformatted')
    write(outputname,'(a,i6.6,a)') 'lastframe_poroelastic_w',myrank,'.bin'
    open(unit=56,file='OUTPUT_FILES/'//outputname,status='unknown',form='unformatted')

      write(55) displs_poroelastic
      write(55) velocs_poroelastic
      write(55) accels_poroelastic
      write(56) displw_poroelastic
      write(56) velocw_poroelastic
      write(56) accelw_poroelastic

    close(55)
    close(56)
  endif

  if(SAVE_FORWARD .and. SIMULATION_TYPE ==1 .and. any_acoustic) then
    if ( myrank == 0 ) then
      write(IOUT,*)
      write(IOUT,*) 'Saving acoustic last frame...'
      write(IOUT,*)
    endif
    write(outputname,'(a,i6.6,a)') 'lastframe_acoustic',myrank,'.bin'
    open(unit=55,file='OUTPUT_FILES/'//outputname,status='unknown',form='unformatted')
      write(55) potential_acoustic
      write(55) potential_dot_acoustic
      write(55) potential_dot_dot_acoustic
    close(55)
  endif


  deallocate(v0x_left)
  deallocate(v0z_left)
  deallocate(t0x_left)
  deallocate(t0z_left)

  deallocate(v0x_right)
  deallocate(v0z_right)
  deallocate(t0x_right)
  deallocate(t0z_right)

  deallocate(v0x_bot)
  deallocate(v0z_bot)
  deallocate(t0x_bot)
  deallocate(t0z_bot)

!----  close energy file
  if(output_energy .and. myrank == 0) close(IOUT_ENERGY)


! print exit banner
  if (myrank == 0) call datim(simulation_title)

!
!----  close output file
!
  if(IOUT /= ISTANDARD_OUTPUT) close(IOUT)

!
!----  end MPI
!
#ifdef USE_MPI
  call MPI_FINALIZE(ier)
#endif


end subroutine finalize_simulation

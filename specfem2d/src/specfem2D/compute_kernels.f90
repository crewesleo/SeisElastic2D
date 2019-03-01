
!========================================================================
!
!                   S P E C F E M 2 D  Version 7 . 0
!                   --------------------------------
!
!     Main historical authors: Dimitri Komatitsch and Jeroen Tromp
!                        Princeton University, USA
!                and CNRS / University of Marseille, France
!                 (there are currently many more authors!)
! (c) Princeton University and CNRS / University of Marseille, April 2014
!
! This software is a computer program whose purpose is to solve
! the two-dimensional viscoelastic anisotropic or poroelastic wave equation
! using a spectral-element method (SEM).
!
! This software is governed by the CeCILL license under French law and
! abiding by the rules of distribution of free software. You can use,
! modify and/or redistribute the software under the terms of the CeCILL
! license as circulated by CEA, CNRS and Inria at the following URL
! "http://www.cecill.info".
!
! As a counterpart to the access to the source code and rights to copy,
! modify and redistribute granted by the license, users are provided only
! with a limited warranty and the software's author, the holder of the
! economic rights, and the successive licensors have only limited
! liability.
!
! In this respect, the user's attention is drawn to the risks associated
! with loading, using, modifying and/or developing or reproducing the
! software by the user in light of its specific status of free software,
! that may mean that it is complicated to manipulate, and that also
! therefore means that it is reserved for developers and experienced
! professionals having in-depth computer knowledge. Users are therefore
! encouraged to load and test the software's suitability as regards their
! requirements in conditions enabling the security of their systems and/or
! data to be ensured and, more generally, to use and operate it in the
! same conditions as regards security.
!
! The full text of the license is available in file "LICENSE".
!
!=====================================================================
!
!-------------------------------------------------------------------------------------------------
!

  subroutine compute_kernels_el()

! elastic kernel calculations
! see e.g. Tromp et al. (2005)

  use specfem_par, only: elastic,rho_k,rhorho_el_hessian_temp1,rhorho_el_hessian_temp2, &
                         rho_kl,mu_kl,kappa_kl,rhop_kl,beta_kl,alpha_kl,bulk_c_kl,bulk_beta_kl, &
                         rhorho_el_hessian_final1,rhorho_el_hessian_final2, &
                         nglob,nspec,ibool,accel_elastic,b_displ_elastic,b_accel_elastic, &
                         rhol_global,mul_global,kappal_global, &
                         hti_c11l_global,hti_c13l_global,hti_c33l_global,hti_c55l_global,hti_rhol_global,&
                         hti_alphal_global,hti_betal_global,hti_epsilonl_global,hti_deltal_global,hti_rhopl_global,&
                         tti_alphal_global,tti_betal_global,tti_epsilonl_global,tti_deltal_global,tti_rhopl_global,&
                         hti_etal_global,hti_alphahl_global,hti_alphanl_global, &
                         ipl_global, isl_global, lambdal_global, sigmal_global, vpl_global, vsl_global, &
                         density,poroelastcoef,kmato,assign_external_model,rhoext,vsext,vpext,&
                         c11ext,c13ext,c33ext,c55ext,c15ext,c35ext,&
                         c11uext,c13uext,c33uext,c55uext,thetaext,&
                         m1ext,m2ext,m3ext,m4ext,m5ext,m6ext,&
                         deltat,p_sv,displ_elastic,&
                         mu_k,kappa_k,elastic,ibool,hprime_xx,hprime_zz,xix,xiz,gammax,gammaz,&
                         hti_ec_c11_k,hti_ec_c13_k,hti_ec_c33_k,hti_ec_c55_k,hti_ec_rho_k, &
                         hti_ec_c11_kl,hti_ec_c13_kl,hti_ec_c33_kl,hti_ec_c55_kl,hti_ec_rho_kl, &
                         hti_coeffl_global, tti_coeffl_global, &
                         hti_thom_alpha_kl,hti_thom_beta_kl,hti_thom_epsilon_kl,hti_thom_delta_kl, &
                         hti_thom_rhop_kl, ANISO,&
                         pdh_hti_ec_c11,pdh_hti_ec_c13,pdh_hti_ec_c33,pdh_hti_ec_c55,pdh_hti_ec_rho, &
                         pdh_hti_thom_alpha,pdh_hti_thom_beta,pdh_hti_thom_epsilon,pdh_hti_thom_delta,pdh_hti_thom_rhop, &
                         pdh_hti_ec_c11_temp,pdh_hti_ec_c13_temp,pdh_hti_ec_c33_temp,pdh_hti_ec_c55_temp,pdh_hti_ec_rho_temp, &
                         pdh_hti_ec_c11_temp1,pdh_hti_ec_c13_temp1,pdh_hti_ec_c33_temp1,pdh_hti_ec_c55_temp1,pdh_hti_ec_rho_temp1, &
                         pdh_hti_thom_alpha_temp,pdh_hti_thom_beta_temp,pdh_hti_thom_epsilon_temp,pdh_hti_thom_delta_temp, &
                         pdh_hti_thom_rhop_temp,&
                         Qalpha_k,Qbeta_k,Qkappa_k,Qmu_k,Qalpha_kl,Qbeta_kl,Qkappa_kl,Qmu_kl,&
                         Qalphal_global,Qbetal_global,Qkappal_global,Qmul_global,&
                         hti_thom2_alpha_kl,hti_thom2_beta_kl,hti_thom2_epsilon_kl,hti_thom2_delta_kl,hti_thom2_eta_kl,&
                         hti_thom2_rhop_kl,&
                         pdh_hti_thom2_alpha,pdh_hti_thom2_beta,pdh_hti_thom2_epsilon,pdh_hti_thom2_delta,pdh_hti_thom2_rhop,&
                         pdh_hti_thom2_eta,&
                         hti_vel_alpha_kl,hti_vel_beta_kl,hti_vel_alphah_kl,hti_vel_alphan_kl,hti_vel_rhop_kl,&
                         pdh_hti_vel_alpha,pdh_hti_vel_beta,pdh_hti_vel_alphah,pdh_hti_vel_alphan,pdh_hti_vel_rhop,&
                         tti_ec_c11_k,tti_ec_c13_k,tti_ec_c15_k,tti_ec_c33_k,tti_ec_c35_k,tti_ec_c55_k,tti_ec_rho_k,&
                         tti_ec_c11_kl,tti_ec_c13_kl,tti_ec_c15_kl,tti_ec_c33_kl,tti_ec_c35_kl,tti_ec_c55_kl,tti_ec_rho_kl,&
                         tti_ecu_c11_kl,tti_ecu_c13_kl,tti_ecu_c33_kl,tti_ecu_c55_kl,tti_ecu_theta_kl,tti_ecu_rho_kl,&
                         pdh_tti_ec_c11_temp,pdh_tti_ec_c13_temp,pdh_tti_ec_c15_temp,pdh_tti_ec_c33_temp,&
                         pdh_tti_ec_c35_temp,pdh_tti_ec_c55_temp,pdh_tti_ec_rho_temp,&
                         pdh_tti_ec_c11_temp1,pdh_tti_ec_c13_temp1,pdh_tti_ec_c15_temp1,pdh_tti_ec_c33_temp1,&
                         pdh_tti_ec_c35_temp1,pdh_tti_ec_c55_temp1,pdh_tti_ec_rho_temp1,&
                         pdh_tti_ec_c11,pdh_tti_ec_c13,pdh_tti_ec_c15,pdh_tti_ec_c33,&
                         pdh_tti_ec_c35,pdh_tti_ec_c55,pdh_tti_ec_rho,&
                         pdh_tti_ecu_c11_temp,pdh_tti_ecu_c13_temp,pdh_tti_ecu_c33_temp,pdh_tti_ecu_c55_temp, &
                         pdh_tti_ecu_theta_temp,pdh_tti_ecu_rho_temp,&
                         pdh_tti_ecu_c11,pdh_tti_ecu_c13,pdh_tti_ecu_c33,pdh_tti_ecu_c55,pdh_tti_ecu_theta,pdh_tti_ecu_rho, &
                         tti_c11l_global,tti_c13l_global,tti_c15l_global,tti_c33l_global,&
                         tti_c35l_global,tti_c55l_global,tti_rhol_global, &
                         tti_c11ul_global,tti_c13ul_global,tti_c33ul_global,tti_c55ul_global,tti_thetaul_global,tti_rhoul_global,&
                         dl_lambda_kl, dl_mu_kl, dl_rho_kl, dip_ip_kl, dip_is_kl, dip_rho_kl, vipi_vp_kl, vipi_vs_kl, vipi_ip_kl, &
                         vipii_vp_kl, vipii_vs_kl, vipii_is_kl, sigma_vp_kl, sigma_sigma_kl, sigma_rho_kl, & 
                         dl_lambda_k, dl_mu_k, dl_rho_k, dip_ip_k, dip_is_k, dip_rho_k, vipi_vp_k, vipi_vs_k, vipi_ip_k, &
                         vipii_vp_k, vipii_vs_k, vipii_is_k, sigma_vp_k, sigma_sigma_k, sigma_rho_k, &
                         pdh_vp, pdh_vs, pdh_rhop, pdh_kappa, pdh_mu, pdh_rho, &
                         pdh_dl_lambda, pdh_dl_mu, pdh_dl_rho, pdh_dip_ip,pdh_dip_is,pdh_dip_rho,&
                         pdh_vipi_vp,pdh_vipi_vs,pdh_vipi_ip, pdh_vipii_vp,pdh_vipii_vs,pdh_vipii_is, &
                         pdh_sigma_vp,pdh_sigma_sigma,pdh_sigma_rho, &
                         pdh_vp_temp, pdh_vs_temp, pdh_rhop_temp, pdh_kappa_temp, pdh_mu_temp, pdh_rho_temp, &
                         pdh_dl_lambda_temp, pdh_dl_mu_temp, pdh_dl_rho_temp, pdh_dip_ip_temp,pdh_dip_is_temp,pdh_dip_rho_temp,&
                         pdh_vipi_vp_temp,pdh_vipi_vs_temp,pdh_vipi_ip_temp, pdh_vipii_vp_temp,pdh_vipii_vs_temp,pdh_vipii_is_temp, &
                         pdh_sigma_vp_temp,pdh_sigma_sigma_temp,pdh_sigma_rho_temp, &
                         pdh_vp_temp1, pdh_vs_temp1, pdh_rhop_temp1, pdh_kappa_temp1, pdh_mu_temp1, pdh_rho_temp1, &
                         pdh_dl_lambda_temp1, pdh_dl_mu_temp1, pdh_dl_rho_temp1, pdh_dip_ip_temp1,pdh_dip_is_temp1,pdh_dip_rho_temp1,&
                         pdh_vipi_vp_temp1,pdh_vipi_vs_temp1,pdh_vipi_ip_temp1, &
                         pdh_vipii_vp_temp1,pdh_vipii_vs_temp1,pdh_vipii_is_temp1, &
                         pdh_sigma_vp_temp1,pdh_sigma_sigma_temp1,pdh_sigma_rho_temp1,&
                         tti_ecu_c11_k,tti_ecu_c13_k,tti_ecu_c33_k,tti_ec_c55_k,tti_ecu_theta_k,tti_ecu_rho_k,&
                         m1l_global,m2l_global,m3l_global,m4l_global,m5l_global,m6l_global,p1l_global,p2l_global,M_PAR,&
                         tti_thom_alpha_kl,tti_thom_beta_kl,tti_thom_epsilon_kl,tti_thom_delta_kl,tti_thom_theta_kl,&
                         tti_thom_rhop_kl,pdh_tti_thom_alpha,pdh_tti_thom_beta,pdh_tti_thom_epsilon,pdh_tti_thom_delta,&
                         pdh_tti_thom_rhop,pdh_tti_thom_theta,pdh_tti_thom_alpha_temp,pdh_tti_thom_beta_temp,pdh_tti_thom_epsilon_temp,&
                         pdh_tti_thom_delta_temp,pdh_tti_thom_rhop_temp,pdh_tti_thom_theta_temp

  implicit none
  include "constants.h"

  !local variables
  integer :: i,j,k,ispec,iglob
  real(kind=CUSTOM_REAL) :: dux_dxi,dux_dgamma,duy_dxi,duy_dgamma,duz_dxi,duz_dgamma
  real(kind=CUSTOM_REAL) :: dux_dxl,duy_dxl,duz_dxl,dux_dzl,duy_dzl,duz_dzl
  real(kind=CUSTOM_REAL) :: b_dux_dxi,b_dux_dgamma,b_duy_dxi,b_duy_dgamma,b_duz_dxi,b_duz_dgamma
  real(kind=CUSTOM_REAL) :: b_dux_dxl,b_duy_dxl,b_duz_dxl,b_dux_dzl,b_duy_dzl,b_duz_dzl
  real(kind=CUSTOM_REAL) :: dsxx,dsxz,dszz
  real(kind=CUSTOM_REAL) :: b_dsxx,b_dsxz,b_dszz

  ! Jacobian matrix and determinant
  double precision :: xixl,xizl,gammaxl,gammazl

  do ispec = 1,nspec
    if( elastic(ispec) ) then
      do j=1,NGLLZ; do i=1,NGLLX
        ! derivative along x and along z
        dux_dxi = 0._CUSTOM_REAL; duy_dxi = 0._CUSTOM_REAL; duz_dxi = 0._CUSTOM_REAL
        dux_dgamma = 0._CUSTOM_REAL; duy_dgamma = 0._CUSTOM_REAL; duz_dgamma = 0._CUSTOM_REAL
        b_dux_dxi = 0._CUSTOM_REAL; b_duy_dxi = 0._CUSTOM_REAL; b_duz_dxi = 0._CUSTOM_REAL
        b_dux_dgamma = 0._CUSTOM_REAL; b_duy_dgamma = 0._CUSTOM_REAL; b_duz_dgamma = 0._CUSTOM_REAL

        ! first double loop over GLL points to compute and store gradients
        ! we can merge the two loops because NGLLX == NGLLZ
        do k = 1,NGLLX
          dux_dxi = dux_dxi + displ_elastic(1,ibool(k,j,ispec))*hprime_xx(i,k)
          duy_dxi = duy_dxi + displ_elastic(2,ibool(k,j,ispec))*hprime_xx(i,k)
          duz_dxi = duz_dxi + displ_elastic(3,ibool(k,j,ispec))*hprime_xx(i,k)
          dux_dgamma = dux_dgamma + displ_elastic(1,ibool(i,k,ispec))*hprime_zz(j,k)
          duy_dgamma = duy_dgamma + displ_elastic(2,ibool(i,k,ispec))*hprime_zz(j,k)
          duz_dgamma = duz_dgamma + displ_elastic(3,ibool(i,k,ispec))*hprime_zz(j,k)


          b_dux_dxi = b_dux_dxi + b_displ_elastic(1,ibool(k,j,ispec))*hprime_xx(i,k)
          b_duy_dxi = b_duy_dxi + b_displ_elastic(2,ibool(k,j,ispec))*hprime_xx(i,k)
          b_duz_dxi = b_duz_dxi + b_displ_elastic(3,ibool(k,j,ispec))*hprime_xx(i,k)
          b_dux_dgamma = b_dux_dgamma + b_displ_elastic(1,ibool(i,k,ispec))*hprime_zz(j,k)
          b_duy_dgamma = b_duy_dgamma + b_displ_elastic(2,ibool(i,k,ispec))*hprime_zz(j,k)
          b_duz_dgamma = b_duz_dgamma + b_displ_elastic(3,ibool(i,k,ispec))*hprime_zz(j,k)

        enddo

        xixl = xix(i,j,ispec)
        xizl = xiz(i,j,ispec)
        gammaxl = gammax(i,j,ispec)
        gammazl = gammaz(i,j,ispec)

        ! derivatives of displacement
        dux_dxl = dux_dxi*xixl + dux_dgamma*gammaxl
        dux_dzl = dux_dxi*xizl + dux_dgamma*gammazl
        duy_dxl = duy_dxi*xixl + duy_dgamma*gammaxl
        duy_dzl = duy_dxi*xizl + duy_dgamma*gammazl
        duz_dxl = duz_dxi*xixl + duz_dgamma*gammaxl
        duz_dzl = duz_dxi*xizl + duz_dgamma*gammazl


        b_dux_dxl = b_dux_dxi*xixl + b_dux_dgamma*gammaxl
        b_dux_dzl = b_dux_dxi*xizl + b_dux_dgamma*gammazl

        b_duy_dxl = b_duy_dxi*xixl + b_duy_dgamma*gammaxl
        b_duy_dzl = b_duy_dxi*xizl + b_duy_dgamma*gammazl

        b_duz_dxl = b_duz_dxi*xixl + b_duz_dgamma*gammaxl
        b_duz_dzl = b_duz_dxi*xizl + b_duz_dgamma*gammazl

        iglob = ibool(i,j,ispec)
        if( p_sv ) then !P-SV waves
          dsxx =  dux_dxl
          dsxz = HALF * (duz_dxl + dux_dzl)
          dszz =  duz_dzl

          b_dsxx =  b_dux_dxl
          b_dsxz = HALF * (b_duz_dxl + b_dux_dzl)
          b_dszz =  b_duz_dzl

          kappa_k(iglob) = (dux_dxl + duz_dzl) *  (b_dux_dxl + b_duz_dzl)
          mu_k(iglob) = dsxx * b_dsxx + dszz * b_dszz + &
                        2._CUSTOM_REAL * dsxz * b_dsxz - 1._CUSTOM_REAL/3._CUSTOM_REAL * kappa_k(iglob)
          !!! ATTENUATION
          Qkappa_k(iglob) = -kappa_k(iglob)
          Qmu_k(iglob) = -mu_k(iglob)
          Qalpha_k(iglob) = Qkappa_k(iglob)
          Qbeta_k(iglob) = Qmu_k(iglob)
          
          ! diagonal pseudo-Hessian
          pdh_kappa_temp(iglob) = (b_dux_dxl + b_duz_dzl) * (b_dux_dxl + b_duz_dzl)
          pdh_mu_temp(iglob) = b_dsxx * b_dsxx + b_dszz * b_dszz + SQRT(2._CUSTOM_REAL) * b_dsxz * b_dsxz - &
                           SQRT(1._CUSTOM_REAL/3._CUSTOM_REAL) * pdh_kappa_temp(iglob)
          pdh_rho_temp(iglob) = b_displ_elastic(1,iglob)**2 + b_displ_elastic(2,iglob)**2 + b_displ_elastic(3,iglob)**2

          !pdh_vp_temp(iglob) = 2 * (1._CUSTOM_REAL + 4._CUSTOM_REAL*mul_global(iglob)/(3._CUSTOM_REAL*kappal_global(iglob))) * &
          !                 pdh_kappa_temp(iglob)
          !pdh_vs_temp(iglob) = pdh_mu_temp(iglob) - 4._CUSTOM_REAL*mul_global(iglob)*pdh_kappa_temp(iglob)/ &
          !                 (3._CUSTOM_REAL*kappal_global(iglob))
          !pdh_rhop_temp(iglob) = pdh_rho_temp(iglob) + pdh_kappa_temp(iglob) + pdh_mu_temp(iglob)

          !!! HTI elastic constants parameterization PWY
          if (ANISO) then
            hti_ec_c11_k(iglob) = dux_dxl * b_dux_dxl
            hti_ec_c13_k(iglob) = (duz_dzl * b_dux_dxl) + (dux_dxl * b_duz_dzl) 
            hti_ec_c33_k(iglob) = duz_dzl * b_duz_dzl
            hti_ec_c55_k(iglob) = (duz_dxl + dux_dzl) * (b_duz_dxl + b_dux_dzl)

            tti_ec_c11_k(iglob) = dux_dxl * b_dux_dxl
            tti_ec_c13_k(iglob) = (duz_dzl * b_dux_dxl) + (dux_dxl * b_duz_dzl)
            tti_ec_c15_k(iglob) = 2._CUSTOM_REAL * (duz_dxl * b_dux_dxl) + 2._CUSTOM_REAL * (dux_dxl * b_duz_dxl)
            tti_ec_c33_k(iglob) = duz_dzl * b_duz_dzl
            tti_ec_c35_k(iglob) = 2._CUSTOM_REAL * (duz_dxl * b_duz_dzl) + 2._CUSTOM_REAL * (duz_dzl * b_duz_dxl)
            tti_ec_c55_k(iglob) = 4._CUSTOM_REAL * (duz_dxl * b_duz_dxl)
            !! diagonal pseudo Hessian
            pdh_hti_ec_c11_temp(iglob) = b_dux_dxl * b_dux_dxl
            pdh_hti_ec_c13_temp(iglob) = (b_duz_dzl + b_dux_dxl) * (b_duz_dzl + b_dux_dxl)
            pdh_hti_ec_c33_temp(iglob) = b_duz_dzl * b_duz_dzl
            pdh_hti_ec_c55_temp(iglob) = (b_duz_dxl + b_dux_dzl) * (b_duz_dxl + b_dux_dzl)

            pdh_hti_ec_c11_temp1(iglob) = b_dux_dxl
            pdh_hti_ec_c13_temp1(iglob) = (b_duz_dzl + b_dux_dxl)
            pdh_hti_ec_c33_temp1(iglob) = b_duz_dzl
            pdh_hti_ec_c55_temp1(iglob) = (b_duz_dxl + b_dux_dzl)

            pdh_tti_ec_c11_temp(iglob) = b_dux_dxl * b_dux_dxl
            pdh_tti_ec_c13_temp(iglob) = (b_duz_dzl + b_dux_dxl) * (b_duz_dzl + b_dux_dxl)
            pdh_tti_ec_c15_temp(iglob) = 4._CUSTOM_REAL * (b_duz_dxl + b_dux_dxl) * (b_duz_dxl + b_dux_dxl)
            pdh_tti_ec_c33_temp(iglob) = b_duz_dzl * b_duz_dzl
            pdh_tti_ec_c35_temp(iglob) = 4._CUSTOM_REAL * (b_duz_dxl + b_duz_dzl) * (b_duz_dxl + b_duz_dzl)
            pdh_tti_ec_c55_temp(iglob) = 4._CUSTOM_REAL * (b_duz_dxl * b_duz_dxl)
          endif
        else !SH (membrane) waves
          mu_k(iglob) = duy_dxl * b_duy_dxl + duy_dzl * b_duy_dzl
        endif
      enddo; enddo
    endif
  enddo

  do iglob = 1,nglob
    rho_k(iglob) =  accel_elastic(1,iglob)*b_displ_elastic(1,iglob) + &
                    accel_elastic(2,iglob)*b_displ_elastic(2,iglob) + &
                    accel_elastic(3,iglob)*b_displ_elastic(3,iglob)
    hti_ec_rho_k(iglob) = rho_k(iglob) ! density kl hti elastic constant parameterization PWY
    tti_ec_rho_k(iglob) = rho_k(iglob)
    rhorho_el_hessian_temp1(iglob) = b_accel_elastic(1,iglob)*b_accel_elastic(1,iglob) + &
                                     b_accel_elastic(2,iglob)*b_accel_elastic(2,iglob) + &
                                     b_accel_elastic(3,iglob)*b_accel_elastic(3,iglob)
    pdh_hti_ec_rho_temp(iglob)     = rhorho_el_hessian_temp1(iglob)
    pdh_tti_ec_rho_temp(iglob)     = rhorho_el_hessian_temp1(iglob)
!    pdh_rho_temp(iglob) = rhorho_el_hessian_temp1(iglob)
!    pdh_vp_temp(iglob) = 2._CUSTOM_REAL * (1._CUSTOM_REAL + 4._CUSTOM_REAL*mul_global(iglob)/(3._CUSTOM_REAL*kappal_global(iglob))) * &
!         pdh_kappa_temp(iglob)
!    pdh_vs_temp(iglob) = pdh_mu_temp(iglob) - 4._CUSTOM_REAL*mul_global(iglob)*pdh_kappa_temp(iglob)/ &
!        (3._CUSTOM_REAL*kappal_global(iglob))
!    pdh_rhop_temp(iglob) = pdh_rho_temp(iglob) + pdh_kappa_temp(iglob) + pdh_mu_temp(iglob)
    rhorho_el_hessian_temp2(iglob) = accel_elastic(1,iglob)*b_accel_elastic(1,iglob) + &
                                     accel_elastic(2,iglob)*b_accel_elastic(2,iglob) + &
                                     accel_elastic(3,iglob)*b_accel_elastic(3,iglob)

 !   pdh_rho_temp(iglob) = b_displ_elastic(1,iglob)**2 + b_displ_elastic(2,iglob)**2 + b_displ_elastic(3,iglob)**2
  enddo

  do ispec = 1, nspec
    if( elastic(ispec) ) then
      do j = 1, NGLLZ
        do i = 1, NGLLX
          iglob = ibool(i,j,ispec)
          if( .not. assign_external_model ) then
            rhol_global(iglob) = density(1,kmato(ispec))
            mul_global(iglob) = poroelastcoef(2,1,kmato(ispec))
            kappal_global(iglob) = poroelastcoef(3,1,kmato(ispec)) - &
                                   4._CUSTOM_REAL*mul_global(iglob) / 3._CUSTOM_REAL
          else
            rhol_global(iglob)   = rhoext(i,j,ispec)
            mul_global(iglob)    = rhoext(i,j,ispec)*vsext(i,j,ispec)*vsext(i,j,ispec)
            kappal_global(iglob) = rhoext(i,j,ispec)*vpext(i,j,ispec)*vpext(i,j,ispec) - &
                                   4._CUSTOM_REAL*mul_global(iglob) / 3._CUSTOM_REAL
            
            lambdal_global(iglob) = kappal_global(iglob) - 2._CUSTOM_REAL*mul_global(iglob) / 3._CUSTOM_REAL
            ipl_global(iglob) = rhoext(i,j,ispec)*vpext(i,j,ispec)
            isl_global(iglob) = rhoext(i,j,ispec)*vsext(i,j,ispec)
            sigmal_global(iglob) = vpext(i,j,ispec) / vsext(i,j,ispec)
            vpl_global(iglob) = vpext(i,j,ispec)
            vsl_global(iglob) = vsext(i,j,ispec)
            !! Attenuation
            Qkappal_global(iglob) = -kappal_global(iglob)
            Qmul_global(iglob) = -mul_global(iglob)
            Qalphal_global(iglob) = Qkappal_global(iglob)
            Qbetal_global(iglob) = Qmul_global(iglob)

            !tti_c11ul_global(iglob) = c11uext(i,j,ispec)
            !tti_c13ul_global(iglob) = c13uext(i,j,ispec)
            !tti_c33ul_global(iglob) = c33uext(i,j,ispec)
            !tti_c55ul_global(iglob) = c55uext(i,j,ispec)
            !tti_thetaul_global(iglob) = thetaext(i,j,ispec)

            pdh_vp_temp(iglob) = SQRT(2 * (1._CUSTOM_REAL + 4._CUSTOM_REAL*mul_global(iglob)/(3._CUSTOM_REAL*kappal_global(iglob)))) * &
                           pdh_kappa_temp(iglob)

            pdh_vs_temp(iglob) = SQRT(2._CUSTOM_REAL) * pdh_mu_temp(iglob) - SQRT(8._CUSTOM_REAL*mul_global(iglob)/&
                           (3._CUSTOM_REAL*kappal_global(iglob)))*pdh_kappa_temp(iglob)

            pdh_rhop_temp(iglob) = pdh_rho_temp(iglob) + pdh_kappa_temp(iglob) + pdh_mu_temp(iglob)

            pdh_dl_lambda_temp(iglob) = SQRT(2*(1 + 2*mul_global(iglob)/(3*kappal_global(iglob)))) * &
                           pdh_kappa_temp(iglob)

            pdh_dl_mu_temp(iglob) = pdh_mu_temp(iglob) + SQRT(2*mul_global(iglob)/(3*kappal_global(iglob)))*&
                           pdh_kappa_temp(iglob)

            pdh_dl_rho_temp(iglob) = pdh_rho_temp(iglob)

            pdh_dip_ip_temp(iglob) = SQRT(2 * (1._CUSTOM_REAL + 4._CUSTOM_REAL*mul_global(iglob)/&
                           (3._CUSTOM_REAL*kappal_global(iglob)))) * pdh_kappa_temp(iglob)

            pdh_dip_is_temp(iglob) = SQRT(2._CUSTOM_REAL) * pdh_mu_temp(iglob) - SQRT(8._CUSTOM_REAL*mul_global(iglob)/&
                           (3._CUSTOM_REAL*kappal_global(iglob)))*pdh_kappa_temp(iglob)

            pdh_dip_rho_temp(iglob) = pdh_rho_temp(iglob) - pdh_kappa_temp(iglob) - pdh_mu_temp(iglob)

            pdh_vipi_vp_temp(iglob) = SQRT(2 * (1._CUSTOM_REAL + 8._CUSTOM_REAL*mul_global(iglob)/&
                           (3._CUSTOM_REAL*kappal_global(iglob)))) * pdh_kappa_temp(iglob) - pdh_mu_temp(iglob) -pdh_rho_temp(iglob)

            pdh_vipi_vs_temp(iglob) = SQRT(2._CUSTOM_REAL) * pdh_mu_temp(iglob) - SQRT(8._CUSTOM_REAL*mul_global(iglob)/&
                           (3._CUSTOM_REAL*kappal_global(iglob)))*pdh_kappa_temp(iglob)

            pdh_vipi_ip_temp(iglob) = pdh_rho_temp(iglob) + pdh_kappa_temp(iglob) + pdh_mu_temp(iglob)

            pdh_vipii_vp_temp(iglob) = SQRT(2 * (1._CUSTOM_REAL + 8._CUSTOM_REAL*mul_global(iglob)/&
                           (3._CUSTOM_REAL*kappal_global(iglob)))) * pdh_kappa_temp(iglob)

            pdh_vipii_vs_temp(iglob) = pdh_mu_temp(iglob) - SQRT(11._CUSTOM_REAL*mul_global(iglob)/&
                           (3._CUSTOM_REAL*kappal_global(iglob)))*pdh_kappa_temp(iglob)-pdh_rho_temp(iglob)

            pdh_vipii_is_temp(iglob) = pdh_rho_temp(iglob) + pdh_kappa_temp(iglob) + pdh_mu_temp(iglob)

            pdh_sigma_vp_temp(iglob) = SQRT(2._CUSTOM_REAL) * (pdh_mu_temp(iglob) + pdh_kappa_temp(iglob))
            
            pdh_sigma_sigma_temp(iglob) = SQRT(2._CUSTOM_REAL) * pdh_mu_temp(iglob) - SQRT(8._CUSTOM_REAL*mul_global(iglob)/&
                           (3._CUSTOM_REAL*kappal_global(iglob)))*pdh_kappa_temp(iglob)
            
            pdh_sigma_rho_temp(iglob) = pdh_rho_temp(iglob) + pdh_kappa_temp(iglob) + pdh_mu_temp(iglob)

            !! HTI parameters global
            if (ANISO) then
              hti_rhol_global(iglob)   = rhoext(i,j,ispec)
              hti_c11l_global(iglob)   = c11ext(i,j,ispec)
              hti_c13l_global(iglob)   = c13ext(i,j,ispec)
              hti_c33l_global(iglob)   = c33ext(i,j,ispec)
              hti_c55l_global(iglob)   = c55ext(i,j,ispec)

              tti_rhol_global(iglob)   = rhoext(i,j,ispec)
              tti_c11l_global(iglob)   = c11ext(i,j,ispec)
              tti_c13l_global(iglob)   = c13ext(i,j,ispec)
              tti_c33l_global(iglob)   = c33ext(i,j,ispec)
              tti_c55l_global(iglob)   = c55ext(i,j,ispec)
              tti_c15l_global(iglob)   = c15ext(i,j,ispec)
              tti_c35l_global(iglob)   = c35ext(i,j,ispec)

              tti_c11ul_global(iglob) = c11uext(i,j,ispec)
              tti_c13ul_global(iglob) = c13uext(i,j,ispec)
              tti_c33ul_global(iglob) = c33uext(i,j,ispec)
              tti_c55ul_global(iglob) = c55uext(i,j,ispec)
              tti_thetaul_global(iglob) = thetaext(i,j,ispec)
              tti_rhoul_global(iglob) = rhoext(i,j,ispec)
              
              m1l_global(iglob) = m1ext(i,j,ispec)
              m2l_global(iglob) = m2ext(i,j,ispec)
              m3l_global(iglob) = m3ext(i,j,ispec)
              m4l_global(iglob) = m4ext(i,j,ispec)
              m5l_global(iglob) = m5ext(i,j,ispec)
              m6l_global(iglob) = m6ext(i,j,ispec)

              p1l_global(iglob) = c11uext(i,j,ispec) - 2._CUSTOM_REAL * c13uext(i,j,ispec) + &
                 c33uext(i,j,ispec) - 4._CUSTOM_REAL * c55uext(i,j,ispec)

              p2l_global(iglob) = c11uext(i,j,ispec) - c13uext(i,j,ispec)
              
              hti_rhopl_global(iglob)    = rhoext(i,j,ispec)
              hti_alphal_global(iglob)   = SQRT(c33ext(i,j,ispec)/rhoext(i,j,ispec))
              hti_betal_global(iglob)    = SQRT(c55ext(i,j,ispec)/rhoext(i,j,ispec))
              hti_deltal_global(iglob) = (hti_c13l_global(iglob) + hti_c55l_global(iglob))**2/(2._CUSTOM_REAL*&
                            hti_c33l_global(iglob)*(hti_c33l_global(iglob)-hti_c55l_global(iglob)))-&
                            (hti_c33l_global(iglob) - hti_c55l_global(iglob))/(2._CUSTOM_REAL * hti_c33l_global(iglob))

              hti_epsilonl_global(iglob)   = (hti_c11l_global(iglob)-hti_c33l_global(iglob))/&
                            (2._CUSTOM_REAL * hti_c33l_global(iglob))

              hti_etal_global(iglob) = (hti_epsilonl_global(iglob) - hti_deltal_global(iglob))/&
                            (1._CUSTOM_REAL + 2._CUSTOM_REAL * hti_deltal_global(iglob))
              
              hti_alphahl_global(iglob) = hti_alphal_global(iglob) * SQRT(1._CUSTOM_REAL+2*hti_epsilonl_global(iglob))
              hti_alphanl_global(iglob) = hti_alphal_global(iglob) * SQRT(1._CUSTOM_REAL+2*hti_deltal_global(iglob))

              hti_coeffl_global(nglob) = SQRT(2._CUSTOM_REAL*hti_deltal_global(iglob)*hti_alphal_global(iglob)*&
                            hti_alphal_global(iglob)*(hti_alphal_global(iglob)*hti_alphal_global(iglob)-&
                            hti_betal_global(iglob)*hti_betal_global(iglob))+&
                            (hti_alphal_global(iglob)*hti_alphal_global(iglob)-&
                            hti_betal_global(iglob)*hti_betal_global(iglob))**2)

              tti_rhopl_global(iglob)    = rhoext(i,j,ispec)
              tti_alphal_global(iglob)   = SQRT(c33uext(i,j,ispec)/rhoext(i,j,ispec))
              tti_betal_global(iglob)    = SQRT(c55uext(i,j,ispec)/rhoext(i,j,ispec))
              tti_deltal_global(iglob) = (tti_c13ul_global(iglob) + tti_c55ul_global(iglob))**2/(2._CUSTOM_REAL*&
                            tti_c33ul_global(iglob)*(tti_c33ul_global(iglob)-tti_c55ul_global(iglob)))-&
                            (tti_c33ul_global(iglob) - tti_c55ul_global(iglob))/(2._CUSTOM_REAL * tti_c33ul_global(iglob))

              tti_epsilonl_global(iglob)   = (tti_c11ul_global(iglob)-tti_c33ul_global(iglob))/&
                            (2._CUSTOM_REAL * tti_c33ul_global(iglob))
              
              tti_coeffl_global(nglob) = SQRT(2._CUSTOM_REAL*tti_deltal_global(iglob)*tti_alphal_global(iglob)*&
                            tti_alphal_global(iglob)*(tti_alphal_global(iglob)*tti_alphal_global(iglob)-&
                            tti_betal_global(iglob)*tti_betal_global(iglob))+&
                            (tti_alphal_global(iglob)*tti_alphal_global(iglob)-&
                            tti_betal_global(iglob)*tti_betal_global(iglob))**2)

              !!! Kernels in TTI unrotated elastic-constant parameterization
              !tti_ecu_c11_k(iglob) = tti_ec_c11_k(iglob)*tti_c11ul_global(iglob)*m1ext(i,j,ispec)/(tti_c11l_global(iglob) * 8) + &
              !                       tti_ec_c13_k(iglob)*tti_c11ul_global(iglob)*m2ext(i,j,ispec)/(tti_c13l_global(iglob)*8) + &
              !                       tti_ec_c15_k(iglob)*tti_c11ul_global(iglob)*m5ext(i,j,ispec)/(tti_c15l_global(iglob)*4) + &
              !                       tti_ec_c33_k(iglob)*tti_c11ul_global(iglob)*m3ext(i,j,ispec)/(tti_c33l_global(iglob)*8) + &
              !                       tti_ec_c35_k(iglob)*tti_c11ul_global(iglob)*m5ext(i,j,ispec)/(tti_c35l_global(iglob)*4) - &
              !                       tti_ec_c35_k(iglob)*tti_c11ul_global(iglob)*2*m4ext(i,j,ispec)/(tti_c35l_global(iglob)*4) + &
              !                       tti_ec_c55_k(iglob)*tti_c11ul_global(iglob)*m2ext(i,j,ispec)/(tti_c55l_global(iglob)*8)

              !!! diagonal pseudo Hessian
              pdh_hti_thom_alpha_temp(iglob) = 2._CUSTOM_REAL * hti_rhopl_global(iglob) * hti_alphal_global(iglob) * &
                pdh_hti_ec_c33_temp1(iglob) + 2._CUSTOM_REAL * hti_rhopl_global(iglob) * hti_alphal_global(iglob) * &
                (2._CUSTOM_REAL * hti_epsilonl_global(iglob) + 1._CUSTOM_REAL) * pdh_hti_ec_c11_temp1(iglob) + &
                4._CUSTOM_REAL*hti_rhopl_global(iglob)*hti_deltal_global(iglob)*hti_alphal_global(iglob)* &
                hti_alphal_global(iglob)*hti_alphal_global(iglob)*pdh_hti_ec_c13_temp1(iglob)/hti_coeffl_global(nglob)-&
                2._CUSTOM_REAL*hti_rhopl_global(iglob)*hti_deltal_global(iglob)*hti_alphal_global(iglob)*&
                hti_betal_global(iglob)*hti_betal_global(iglob)*pdh_hti_ec_c13_temp1(iglob)/hti_coeffl_global(nglob)+&
                2._CUSTOM_REAL*hti_rhopl_global(iglob)*hti_alphal_global(iglob)*&
                (hti_alphal_global(iglob)**2-hti_betal_global(iglob)**2)*pdh_hti_ec_c13_temp1(iglob)/hti_coeffl_global(nglob)

              pdh_hti_thom_beta_temp(iglob) = 2._CUSTOM_REAL*hti_rhopl_global(iglob)*hti_betal_global(iglob)*&
                pdh_hti_ec_c55_temp1(iglob)-2._CUSTOM_REAL*hti_rhopl_global(iglob)*hti_deltal_global(iglob)*hti_alphal_global(iglob)*&
                hti_alphal_global(iglob)*hti_betal_global(iglob)*pdh_hti_ec_c13_temp1(iglob)/hti_coeffl_global(nglob) - &
                2._CUSTOM_REAL*hti_rhopl_global(iglob)*hti_betal_global(iglob)*&
                (hti_alphal_global(iglob)**2-hti_betal_global(iglob)**2)*pdh_hti_ec_c13_temp1(iglob)/hti_coeffl_global(nglob)-& 
                2._CUSTOM_REAL*hti_rhopl_global(iglob)*hti_betal_global(iglob)*pdh_hti_ec_c13_temp1(iglob)
              
               pdh_hti_thom_epsilon_temp(iglob) = 2._CUSTOM_REAL*hti_rhopl_global(iglob)*pdh_hti_ec_c13_temp1(iglob)*&
                hti_alphal_global(iglob)**2

               pdh_hti_thom_delta_temp(iglob) = hti_rhopl_global(iglob)*pdh_hti_ec_c13_temp1(iglob)*&
                (hti_alphal_global(iglob)**2-hti_betal_global(iglob)**2)*hti_alphal_global(iglob)**2/hti_coeffl_global(nglob)

               pdh_hti_thom_rhop_temp(iglob) = pdh_hti_ec_rho_temp(iglob)+pdh_hti_ec_c33_temp1(iglob)*hti_alphal_global(iglob)**2 +&
                pdh_hti_ec_c55_temp1(iglob)*hti_betal_global(iglob)**2+pdh_hti_ec_c13_temp1(iglob)*hti_coeffl_global(nglob)-&
                pdh_hti_ec_c13_temp1(iglob)*hti_betal_global(iglob)**2+pdh_hti_ec_c11_temp1(iglob)*&
                (2._CUSTOM_REAL * hti_epsilonl_global(iglob)+1._CUSTOM_REAL)*hti_alphal_global(iglob)**2

            endif
          endif

          rho_kl(i,j,ispec) = rho_kl(i,j,ispec) - rhol_global(iglob)  * rho_k(iglob) * deltat
          mu_kl(i,j,ispec) =  mu_kl(i,j,ispec) - TWO * mul_global(iglob) * mu_k(iglob) * deltat
          kappa_kl(i,j,ispec) = kappa_kl(i,j,ispec) - kappal_global(iglob) * kappa_k(iglob) * deltat
          ! attenuation
          Qkappa_kl(i,j,ispec) = -1._CUSTOM_REAL * kappa_kl(i,j,ispec)
          Qmu_kl(i,j,ispec) = -1._CUSTOM_REAL * mu_kl(i,j,ispec)
          !Qmu_kl(i,j,ispec) = 1._CUSTOM_REAL * kappa_kl(i,j,ispec)
          Qalpha_kl(i,j,ispec) = -1._CUSTOM_REAL * kappa_kl(i,j,ispec)
          Qbeta_kl(i,j,ispec) = -1._CUSTOM_REAL * mu_kl(i,j,ispec)
          
          Qalpha_kl(i,j,ispec) = -1._CUSTOM_REAL * kappa_kl(i,j,ispec)*(3.0*kappal_global(iglob)+4.0*mul_global(iglob))/ &
                               (3.0*kappal_global(iglob)+mul_global(iglob))
          Qalpha_kl(i,j,ispec) = 0.5 * Qalpha_kl(i,j,ispec) - 0.5 * kappa_kl(i,j,ispec)*(3.0*mul_global(iglob)*Qkappal_global(iglob))/ &
                               ((3.0*kappal_global(iglob)+4.0*mul_global(iglob))*Qmul_global(iglob) - &
                                3.0*mul_global(iglob)*Qkappal_global(iglob)) + 0.5 * Qmu_kl(i,j,ispec)
          ! 2018-06-27 by PWY
          !Qalpha_kl(i,j,ispec) = -1*kappa_kl(i,j,ispec)*Qmul_global(iglob)*(kappal_global(iglob) + 4*mul_global(iglob)/3)/&
          !              (Qmul_global(iglob)*(kappal_global(iglob)+4*mul_global(iglob)/3)-mul_global(iglob)*Qkappal_global(iglob))
          
          Qkappa_kl(i,j,ispec) = 0.5*Qalpha_kl(i,j,ispec) + 0.5* Qmu_kl(i,j,ispec) 
          
          !
          rhop_kl(i,j,ispec) = rho_kl(i,j,ispec) + kappa_kl(i,j,ispec) + mu_kl(i,j,ispec)
          beta_kl(i,j,ispec) = TWO * (mu_kl(i,j,ispec) - 4._CUSTOM_REAL * mul_global(iglob)/&
                        (3._CUSTOM_REAL * kappal_global(iglob)) * kappa_kl(i,j,ispec))
          alpha_kl(i,j,ispec) = TWO * (1._CUSTOM_REAL + 4._CUSTOM_REAL * mul_global(iglob)/&
                         (3._CUSTOM_REAL * kappal_global(iglob))) * kappa_kl(i,j,ispec)

          
          dl_lambda_kl(i,j,ispec) = (1._CUSTOM_REAL - 2._CUSTOM_REAL*mul_global(iglob)/(3._CUSTOM_REAL*kappal_global(iglob)))*&
                     kappa_kl(i,j,ispec) 

          dl_mu_kl(i,j,ispec) = mu_kl(i,j,ispec) + 2._CUSTOM_REAL*mul_global(iglob)*kappa_kl(i,j,ispec)/&
                     (3._CUSTOM_REAL*kappal_global(iglob))

          dl_rho_kl(i,j,ispec) = rho_kl(i,j,ispec)
          
          dip_ip_kl(i,j,ispec) = alpha_kl(i,j,ispec)
          dip_is_kl(i,j,ispec) = beta_kl(i,j,ispec)
          dip_rho_kl(i,j,ispec) = rho_kl(i,j,ispec) - kappa_kl(i,j,ispec) - mu_kl(i,j,ispec)

          vipi_vp_kl(i,j,ispec) = alpha_kl(i,j,ispec) - rhop_kl(i,j,ispec)
          vipi_vs_kl(i,j,ispec) = beta_kl(i,j,ispec)
          vipi_ip_kl(i,j,ispec) = rhop_kl(i,j,ispec)

          vipii_vp_kl(i,j,ispec) = alpha_kl(i,j,ispec)
          vipii_vs_kl(i,j,ispec) = beta_kl(i,j,ispec) - rhop_kl(i,j,ispec)
          vipii_is_kl(i,j,ispec) = rhop_kl(i,j,ispec)
                                                        
          sigma_vp_kl(i,j,ispec) = alpha_kl(i,j,ispec) + beta_kl(i,j,ispec)
          sigma_sigma_kl(i,j,ispec) = - beta_kl(i,j,ispec)
          sigma_rho_kl(i,j,ispec) = rhop_kl(i,j,ispec)

          !! Diagonal pseudo-Hessian in isotropic-elastic media
          pdh_vp(i,j,ispec) = pdh_vp(i,j,ispec) + pdh_vp_temp(iglob)*pdh_vp_temp(iglob)*deltat
          pdh_vs(i,j,ispec) = pdh_vs(i,j,ispec) + pdh_vs_temp(iglob)*pdh_vs_temp(iglob)*deltat
          pdh_rhop(i,j,ispec) = pdh_rhop(i,j,ispec) + 0.2 * pdh_vs_temp(iglob)*pdh_vs_temp(iglob)*deltat

          pdh_kappa(i,j,ispec) = pdh_kappa(i,j,ispec) +kappal_global(iglob)* pdh_kappa_temp(iglob)*&
                             pdh_kappa_temp(iglob)*deltat
          pdh_mu(i,j,ispec) = pdh_mu(i,j,ispec) + mul_global(iglob)*pdh_mu_temp(iglob)*&
                             pdh_mu_temp(iglob)*deltat
          pdh_rho(i,j,ispec) = pdh_rho(i,j,ispec) + rhol_global(iglob)*pdh_rho_temp(iglob)*&
                             pdh_rho_temp(iglob)*deltat

          pdh_dl_lambda(i,j,ispec) = pdh_dl_lambda(i,j,ispec) + lambdal_global(iglob)*pdh_dl_lambda_temp(iglob)*&
                             pdh_dl_lambda_temp(iglob)*deltat

          pdh_dl_mu(i,j,ispec) = pdh_dl_mu(i,j,ispec) + mul_global(iglob)*pdh_dl_mu_temp(iglob)*&
                             pdh_dl_mu_temp(iglob)*deltat

          pdh_dl_rho(i,j,ispec) = pdh_dl_rho(i,j,ispec) + rhol_global(iglob)*pdh_dl_rho_temp(iglob)*&
                             pdh_dl_rho_temp(iglob)*deltat

          pdh_dip_ip(i,j,ispec) = pdh_dip_ip(i,j,ispec) + ipl_global(iglob)*pdh_dip_ip_temp(iglob)*&
                             pdh_dip_ip_temp(iglob)*deltat

          pdh_dip_is(i,j,ispec) = pdh_dip_is(i,j,ispec) + isl_global(iglob)*pdh_dip_is_temp(iglob)*&
                             pdh_dip_is_temp(iglob)*deltat

          pdh_dip_rho(i,j,ispec) = pdh_dip_rho(i,j,ispec) + pdh_dip_rho_temp(iglob)*&
                             pdh_dip_rho_temp(iglob)*deltat

          pdh_vipi_vp(i,j,ispec) = pdh_vipi_vp(i,j,ispec) + vsl_global(iglob)*pdh_vipi_vs_temp(iglob)*&
                             pdh_vipi_vs_temp(iglob)*deltat

          pdh_vipi_vs(i,j,ispec) = pdh_vipi_vs(i,j,ispec) + vsl_global(iglob)*pdh_vipi_vs_temp(iglob)*&
                             pdh_vipi_vs_temp(iglob)*deltat

          pdh_vipi_ip(i,j,ispec) = pdh_vipi_ip(i,j,ispec) + vsl_global(iglob)*pdh_vipi_vs_temp(iglob)*&
                             pdh_vipi_vs_temp(iglob)*deltat

          pdh_vipii_vp(i,j,ispec) = pdh_vipii_vp(i,j,ispec) + vsl_global(iglob)*pdh_vipii_vs_temp(iglob)*&
                             pdh_vipii_vs_temp(iglob)*deltat

          pdh_vipii_vs(i,j,ispec) = pdh_vipii_vs(i,j,ispec) + vsl_global(iglob)*pdh_vipii_vs_temp(iglob)*&
                             pdh_vipii_vs_temp(iglob)*deltat

          pdh_vipii_is(i,j,ispec) = pdh_vipii_is(i,j,ispec) + vsl_global(iglob)*pdh_vipii_vs_temp(iglob)*&
                             pdh_vipii_vs_temp(iglob)*deltat

          pdh_sigma_vp(i,j,ispec) = pdh_sigma_vp(i,j,ispec) + rhol_global(iglob)*pdh_sigma_rho_temp(iglob)*&
                             pdh_sigma_rho_temp(iglob)*deltat

          pdh_sigma_sigma(i,j,ispec) = pdh_sigma_sigma(i,j,ispec) + rhol_global(iglob)*pdh_sigma_rho_temp(iglob)*&
                             pdh_sigma_rho_temp(iglob)*deltat

          pdh_sigma_rho(i,j,ispec) = pdh_sigma_rho(i,j,ispec) + rhol_global(iglob)*pdh_sigma_rho_temp(iglob)*&
                             pdh_sigma_rho_temp(iglob)*deltat

          !pdh_vp(i,j,ispec) = 2 * (1._CUSTOM_REAL + 4._CUSTOM_REAL*mul_global(iglob)/(3._CUSTOM_REAL*kappal_global(iglob)))*&
          !                  pdh_kappa(i,j,ispec)

          !pdh_vs(i,j,ispec) = pdh_mu(i,j,ispec) - (4._CUSTOM_REAL*mul_global(iglob)/(3._CUSTOM_REAL*kappal_global(iglob)))*&
          !                 pdh_kappa(i,j,ispec)

          !pdh_rhop(i,j,ispec) = pdh_kappa(i,j,ispec) + pdh_mu(i,j,ispec) + pdh_rho(i,j,ispec)
         
          !
          bulk_c_kl(i,j,ispec) =  TWO * kappa_kl(i,j,ispec)
          bulk_beta_kl(i,j,ispec) =  TWO * mu_kl(i,j,ispec)

          rhorho_el_hessian_final1(i,j,ispec) = rhorho_el_hessian_final1(i,j,ispec) + &
                                  rhorho_el_hessian_temp1(iglob) * deltat
          rhorho_el_hessian_final2(i,j,ispec) = rhorho_el_hessian_final2(i,j,ispec) + &
                                  rhorho_el_hessian_temp2(iglob) * deltat
          
          ! Kernels in HTI media with elastic constants parameterization
          if (ANISO) then
            hti_ec_rho_kl(i,j,ispec) = hti_ec_rho_kl(i,j,ispec) - hti_rhol_global(iglob) * &
                                  hti_ec_rho_k(iglob) * deltat
            hti_ec_c11_kl(i,j,ispec) = hti_ec_c11_kl(i,j,ispec) - hti_c11l_global(iglob) * &
                                  hti_ec_c11_k(iglob) * deltat
            hti_ec_c13_kl(i,j,ispec) = hti_ec_c13_kl(i,j,ispec) - hti_c13l_global(iglob) * &
                                  hti_ec_c13_k(iglob) * deltat
            hti_ec_c33_kl(i,j,ispec) = hti_ec_c33_kl(i,j,ispec) - hti_c33l_global(iglob) * &
                                  hti_ec_c33_k(iglob) * deltat
            hti_ec_c55_kl(i,j,ispec) = hti_ec_c55_kl(i,j,ispec) - hti_c55l_global(iglob) * &
                                  hti_ec_c55_k(iglob) * deltat

            tti_ec_rho_kl(i,j,ispec) = tti_ec_rho_kl(i,j,ispec) - tti_rhol_global(iglob) * &
                                  tti_ec_rho_k(iglob) * deltat
            tti_ec_c11_kl(i,j,ispec) = tti_ec_c11_kl(i,j,ispec) - tti_c11l_global(iglob) * &
                                  tti_ec_c11_k(iglob) * deltat
            tti_ec_c13_kl(i,j,ispec) = tti_ec_c13_kl(i,j,ispec) - tti_c13l_global(iglob) * &
                                  tti_ec_c13_k(iglob) * deltat
            tti_ec_c33_kl(i,j,ispec) = tti_ec_c33_kl(i,j,ispec) - tti_c33l_global(iglob) * &
                                  tti_ec_c33_k(iglob) * deltat
            tti_ec_c55_kl(i,j,ispec) = tti_ec_c55_kl(i,j,ispec) - tti_c55l_global(iglob) * &
                                  tti_ec_c55_k(iglob) * deltat
            tti_ec_c15_kl(i,j,ispec) = tti_ec_c15_kl(i,j,ispec) - tti_c15l_global(iglob) * &
                                  tti_ec_c15_k(iglob) * deltat
            tti_ec_c35_kl(i,j,ispec) = tti_ec_c35_kl(i,j,ispec) - tti_c35l_global(iglob) * &
                                  tti_ec_c35_k(iglob) * deltat

            tti_ecu_c11_kl(i,j,ispec) = tti_ecu_c11_kl(i,j,ispec) - tti_c11ul_global(iglob) * &
                                  tti_ecu_c11_k(iglob) * deltat

          !! Diagonal Hessian
            pdh_hti_ec_c11(i,j,ispec) = pdh_hti_ec_c11(i,j,ispec) + &
                                  pdh_hti_ec_c11_temp(iglob) * deltat
            pdh_hti_ec_c13(i,j,ispec) = pdh_hti_ec_c13(i,j,ispec) + &
                                  pdh_hti_ec_c13_temp(iglob) * deltat
            pdh_hti_ec_c33(i,j,ispec) = pdh_hti_ec_c33(i,j,ispec) + &
                                  pdh_hti_ec_c33_temp(iglob) * deltat
            pdh_hti_ec_c55(i,j,ispec) = pdh_hti_ec_c55(i,j,ispec) + &
                                  pdh_hti_ec_c55_temp(iglob) * deltat
            pdh_hti_ec_rho(i,j,ispec) = pdh_hti_ec_rho(i,j,ispec) + &
                                  pdh_hti_ec_rho_temp(iglob) * deltat

          !! Diagonal Hessian
            pdh_tti_ec_c11(i,j,ispec) = pdh_tti_ec_c11(i,j,ispec) + &
                                  pdh_tti_ec_c11_temp(iglob) * deltat
            pdh_tti_ec_c13(i,j,ispec) = pdh_tti_ec_c13(i,j,ispec) + &
                                  pdh_tti_ec_c13_temp(iglob) * deltat
            pdh_tti_ec_c33(i,j,ispec) = pdh_tti_ec_c33(i,j,ispec) + &
                                  pdh_tti_ec_c33_temp(iglob) * deltat
            pdh_tti_ec_c55(i,j,ispec) = pdh_tti_ec_c55(i,j,ispec) + &
                                  pdh_tti_ec_c55_temp(iglob) * deltat
            pdh_tti_ec_c15(i,j,ispec) = pdh_tti_ec_c15(i,j,ispec) + &
                                  pdh_tti_ec_c15_temp(iglob) * deltat
            pdh_tti_ec_c35(i,j,ispec) = pdh_tti_ec_c35(i,j,ispec) + &
                                  pdh_tti_ec_c35_temp(iglob) * deltat
            pdh_tti_ec_rho(i,j,ispec) = pdh_tti_ec_rho(i,j,ispec) + &
                                  pdh_hti_ec_rho_temp(iglob) * deltat

            pdh_hti_thom_alpha(i,j,ispec) = pdh_hti_thom_alpha(i,j,ispec)+&
                            pdh_hti_thom_alpha_temp(iglob) * pdh_hti_thom_alpha_temp(iglob) * deltat *&
                            hti_alphal_global(iglob)

            pdh_hti_thom_beta(i,j,ispec) = pdh_hti_thom_beta(i,j,ispec)+&
                            pdh_hti_thom_beta_temp(iglob) * pdh_hti_thom_beta_temp(iglob) * deltat * &
                            hti_betal_global(iglob)

            pdh_hti_thom_epsilon(i,j,ispec) = pdh_hti_thom_epsilon(i,j,ispec)+&
                            pdh_hti_thom_epsilon_temp(iglob) * pdh_hti_thom_epsilon_temp(iglob) * deltat * &
                            hti_epsilonl_global(iglob)

            pdh_hti_thom_delta(i,j,ispec) = pdh_hti_thom_delta(i,j,ispec)+&
                            pdh_hti_thom_delta_temp(iglob) * pdh_hti_thom_delta_temp(iglob) * deltat * &
                            hti_deltal_global(iglob)

            pdh_hti_thom_rhop(i,j,ispec) = pdh_hti_thom_rhop(i,j,ispec)+&
                            pdh_hti_thom_rhop_temp(iglob) * pdh_hti_thom_rhop_temp(iglob) * deltat * &
                            hti_rhopl_global(iglob)

            pdh_hti_thom2_alpha(i,j,ispec) = pdh_hti_thom2_alpha(i,j,ispec)+&
                            pdh_hti_thom_alpha_temp(iglob) * pdh_hti_thom_alpha_temp(iglob) * deltat *&
                            hti_alphal_global(iglob)

            pdh_hti_thom2_beta(i,j,ispec) = pdh_hti_thom2_beta(i,j,ispec)+&
                            pdh_hti_thom_beta_temp(iglob) * pdh_hti_thom_beta_temp(iglob) * deltat * &
                            hti_betal_global(iglob)

            pdh_hti_thom2_epsilon(i,j,ispec) = pdh_hti_thom2_epsilon(i,j,ispec)+&
                            pdh_hti_thom_epsilon_temp(iglob) * pdh_hti_thom_epsilon_temp(iglob) * deltat * &
                            hti_epsilonl_global(iglob)

            pdh_hti_thom2_delta(i,j,ispec) = pdh_hti_thom2_delta(i,j,ispec)+&
                            pdh_hti_thom_delta_temp(iglob) * pdh_hti_thom_delta_temp(iglob) * deltat * &
                            hti_deltal_global(iglob)

            pdh_hti_thom2_rhop(i,j,ispec) = pdh_hti_thom2_rhop(i,j,ispec)+&
                            pdh_hti_thom_rhop_temp(iglob) * pdh_hti_thom_rhop_temp(iglob) * deltat * &
                            hti_rhopl_global(iglob)

            pdh_hti_thom2_eta(i,j,ispec) = pdh_hti_thom2_eta(i,j,ispec)+&
                            pdh_hti_thom_delta_temp(iglob) * pdh_hti_thom_delta_temp(iglob) * deltat * &
                            hti_deltal_global(iglob)

            pdh_hti_vel_alpha(i,j,ispec) = pdh_hti_vel_alpha(i,j,ispec)+&
                            pdh_hti_thom_alpha_temp(iglob) * pdh_hti_thom_alpha_temp(iglob) * deltat *&
                            hti_alphal_global(iglob)

            pdh_hti_vel_beta(i,j,ispec) = pdh_hti_vel_beta(i,j,ispec)+&
                            pdh_hti_thom_beta_temp(iglob) * pdh_hti_thom_beta_temp(iglob) * deltat * &
                            hti_betal_global(iglob)

            pdh_hti_vel_alphah(i,j,ispec) = pdh_hti_vel_alphah(i,j,ispec)+&
                            pdh_hti_thom_alpha_temp(iglob) * pdh_hti_thom_alpha_temp(iglob) * deltat *&
                            hti_alphal_global(iglob)

            pdh_hti_vel_alphan(i,j,ispec) = pdh_hti_vel_alphan(i,j,ispec)+&
                            pdh_hti_thom_alpha_temp(iglob) * pdh_hti_thom_alpha_temp(iglob) * deltat *&
                            hti_alphal_global(iglob)

            pdh_hti_vel_rhop(i,j,ispec) = pdh_hti_vel_rhop(i,j,ispec)+&
                            pdh_hti_thom_rhop_temp(iglob) * pdh_hti_thom_rhop_temp(iglob) * deltat * &
                            hti_rhopl_global(iglob)
          ! Kernels in HTI media with Thomsen parameterization
            hti_thom_rhop_kl(i,j,ispec) = hti_thom_rhop_kl(i,j,ispec) - hti_rhopl_global(iglob) * &
             hti_ec_rho_k(iglob)*deltat - hti_rhopl_global(iglob) * hti_alphal_global(iglob)*hti_alphal_global(iglob)*&
             hti_ec_c33_k(iglob)*deltat - hti_rhopl_global(iglob) * hti_betal_global(iglob)*hti_betal_global(iglob)*&
             hti_ec_c55_k(iglob)*deltat - hti_rhopl_global(iglob) * (hti_coeffl_global(nglob)-hti_betal_global(iglob)*&
             hti_betal_global(iglob))*hti_ec_c13_k(iglob)*deltat+hti_alphal_global(iglob)*hti_alphal_global(iglob)*&
             (2._CUSTOM_REAL*hti_epsilonl_global(iglob)+1._CUSTOM_REAL)*hti_ec_c11_k(iglob)*deltat

            hti_thom_alpha_kl(i,j,ispec) = hti_thom_alpha_kl(i,j,ispec) - 2._CUSTOM_REAL*hti_alphal_global(iglob)*&
             hti_alphal_global(iglob)*hti_rhopl_global(iglob)*hti_ec_c33_k(iglob)*deltat-&
             2._CUSTOM_REAL*hti_alphal_global(iglob)*hti_alphal_global(iglob)*hti_rhopl_global(iglob)*&
             (2._CUSTOM_REAL*hti_epsilonl_global(iglob)+1._CUSTOM_REAL)*hti_ec_c11_k(iglob)*deltat-&
             4._CUSTOM_REAL*hti_rhopl_global(iglob)*hti_deltal_global(iglob)*hti_alphal_global(iglob)*&
             hti_alphal_global(iglob)*hti_alphal_global(iglob)*hti_alphal_global(iglob)*hti_ec_c13_k(iglob)*deltat/&
             hti_coeffl_global(nglob)+2._CUSTOM_REAL*hti_rhopl_global(iglob)*hti_deltal_global(iglob)*hti_alphal_global(iglob)*&
             hti_alphal_global(iglob)*hti_betal_global(iglob)*hti_betal_global(iglob)*hti_ec_c13_k(iglob)*deltat/&
             hti_coeffl_global(nglob)-2._CUSTOM_REAL*hti_rhopl_global(iglob)*hti_alphal_global(iglob)*&
             hti_alphal_global(iglob)*(hti_alphal_global(iglob)**2-hti_betal_global(iglob)**2)*&
             hti_ec_c13_k(iglob)*deltat/hti_coeffl_global(nglob)

            hti_thom_beta_kl(i,j,ispec) = hti_thom_beta_kl(i,j,ispec) - 2._CUSTOM_REAL*hti_rhopl_global(iglob)*&
             hti_betal_global(iglob)*hti_betal_global(iglob)*hti_ec_c55_k(iglob)*deltat+&
             2._CUSTOM_REAL*hti_rhopl_global(iglob)*hti_deltal_global(iglob)*hti_alphal_global(iglob)*&
             hti_alphal_global(iglob)*hti_betal_global(iglob)*hti_betal_global(iglob)*hti_ec_c13_k(iglob)*deltat/&
             hti_coeffl_global(nglob) + 2._CUSTOM_REAL*hti_rhopl_global(iglob)*hti_betal_global(iglob)*&
             hti_betal_global(iglob)*(hti_alphal_global(iglob)**2-hti_betal_global(iglob)**2)*hti_ec_c13_k(iglob)*deltat/&
             hti_coeffl_global(nglob) + 2._CUSTOM_REAL*hti_rhopl_global(iglob)*hti_betal_global(iglob)*&
             hti_betal_global(iglob)*hti_ec_c13_k(iglob)*deltat

            hti_thom_epsilon_kl(i,j,ispec) = hti_thom_epsilon_kl(i,j,ispec) - 2._CUSTOM_REAL*hti_rhopl_global(iglob)*&
             hti_alphal_global(iglob)*hti_alphal_global(iglob)*hti_epsilonl_global(iglob)*hti_ec_c11_k(iglob)*deltat

            hti_thom_delta_kl(i,j,ispec) = hti_thom_delta_kl(i,j,ispec) - hti_rhopl_global(iglob)*&
             hti_alphal_global(iglob)*hti_alphal_global(iglob)*(hti_alphal_global(iglob)**2-hti_betal_global(iglob)**2)*&
             hti_deltal_global(iglob)*hti_ec_c13_k(iglob)*deltat/hti_coeffl_global(nglob)

            ! Horizontal velocity and thomsen parameter parameterization
            hti_thom2_rhop_kl(i,j,ispec) = hti_thom_rhop_kl(i,j,ispec)

            hti_thom2_alpha_kl(i,j,ispec) = hti_thom_alpha_kl(i,j,ispec)

            hti_thom2_beta_kl(i,j,ispec) = hti_thom_beta_kl(i,j,ispec)

            hti_thom2_epsilon_kl(i,j,ispec)= hti_thom_epsilon_kl(i,j,ispec)+ &
             hti_thom_delta_kl(i,j,ispec)*hti_epsilonl_global(iglob)/&
             (hti_epsilonl_global(iglob)-hti_etal_global(iglob)) - hti_thom_alpha_kl(i,j,ispec)*hti_epsilonl_global(iglob)/&
             (1._CUSTOM_REAL + 2._CUSTOM_REAL * hti_epsilonl_global(iglob))

            hti_thom2_delta_kl(i,j,ispec) = hti_thom_delta_kl(i,j,ispec)

            hti_thom2_eta_kl(i,j,ispec) = - 1._CUSTOM_REAL * hti_thom_delta_kl(i,j,ispec)*2._CUSTOM_REAL*hti_etal_global(iglob)/&
             (1._CUSTOM_REAL + 2._CUSTOM_REAL * hti_etal_global(iglob))-&
             hti_thom_delta_kl(i,j,ispec) * hti_etal_global(iglob)/(hti_epsilonl_global(iglob) - hti_etal_global(iglob))

            ! Velocity parameter parameterization
            hti_vel_rhop_kl(i,j,ispec) = hti_thom_rhop_kl(i,j,ispec)
            hti_vel_beta_kl(i,j,ispec) = hti_thom_beta_kl(i,j,ispec)
!            hti_vel_alpha_kl(i,j,ispec) = hti_thom_alpha_kl(i,j,ispec)
!            hti_vel_alpha_kl(i,j,ispec) = alpha_kl(i,j,ispec)
!            hti_vel_alphah_kl(i,j,ispec) = hti_thom_alpha_kl(i,j,ispec)
!            hti_vel_alphan_kl(i,j,ispec) = hti_thom_alpha_kl(i,j,ispec)
            hti_vel_alpha_kl(i,j,ispec) = hti_thom_alpha_kl(i,j,ispec)-2._CUSTOM_REAL*hti_thom_epsilon_kl(i,j,ispec)* &
             hti_alphahl_global(iglob)**2/(hti_alphahl_global(iglob)**2-hti_alphal_global(iglob)**2) - &
             2._CUSTOM_REAL*hti_thom_delta_kl(i,j,ispec)*hti_alphanl_global(iglob)**2/&
             (hti_alphanl_global(iglob)**2-hti_alphal_global(iglob)**2)

            hti_vel_alphah_kl(i,j,ispec) = 2._CUSTOM_REAL*hti_thom_epsilon_kl(i,j,ispec)*hti_alphahl_global(iglob)**2/&
             (hti_alphahl_global(iglob)**2-hti_alphal_global(iglob)**2)

            hti_vel_alphan_kl(i,j,ispec) = 2._CUSTOM_REAL*hti_thom_delta_kl(i,j,ispec)*hti_alphanl_global(iglob)**2/&
             (hti_alphanl_global(iglob)**2-hti_alphal_global(iglob)**2)

            !hti_thom_epsilon_kl(i,j,ispec) = hti_thom_epsilon_kl(i,j,ispec)-2._CUSTOM_REAL*hti_rhopl_global(iglob)*&
             !hti_alphal_global(iglob)*hti_alphal_global(iglob)*hti_ec_c11_k(iglob)*deltat
            !hti_thom_delta_kl(i,j,ispec) = hti_thom_delta_kl(i,j,ispec) - hti_rhopl_global(iglob)*&
             !hti_alphal_global(iglob)*hti_alphal_global(iglob)*(hti_alphal_global(iglob)**2-hti_betal_global(iglob)**2)*&
             !hti_ec_c13_k(iglob)*deltat/hti_coeffl_global(nglob)
            if ((trim(M_PAR) == 'ttiecu') .OR. (trim(M_PAR) == 'ttithom')) then
              tti_ecu_c11_kl(i,j,ispec) = tti_ec_c11_kl(i,j,ispec) * tti_c11ul_global(iglob) * m1l_global(iglob) / &
                 (tti_c11l_global(iglob) * 8._CUSTOM_REAL) + &
                 tti_ec_c13_kl(i,j,ispec) * tti_c11ul_global(iglob) * m2l_global(iglob) / (tti_c13l_global(iglob) * 8._CUSTOM_REAL) + &
                 tti_ec_c15_kl(i,j,ispec) * tti_c11ul_global(iglob) * m5l_global(iglob) / (tti_c15l_global(iglob) * 4._CUSTOM_REAL) + &
                 tti_ec_c33_kl(i,j,ispec) * tti_c11ul_global(iglob) * m3l_global(iglob) / (tti_c33l_global(iglob) * 8._CUSTOM_REAL) + &
                 tti_ec_c35_kl(i,j,ispec) * tti_c11ul_global(iglob) * m5l_global(iglob) / (tti_c35l_global(iglob) * 4._CUSTOM_REAL) - &
                 tti_ec_c35_kl(i,j,ispec) * tti_c11ul_global(iglob) * 2._CUSTOM_REAL * m4l_global(iglob) / &
                 (tti_c35l_global(iglob) * 4._CUSTOM_REAL) + &
                 tti_ec_c55_kl(i,j,ispec) * tti_c11ul_global(iglob) * m2l_global(iglob)/(tti_c55l_global(iglob) * 8._CUSTOM_REAL)

              tti_ecu_c13_kl(i,j,ispec) = tti_ec_c11_kl(i,j,ispec) * tti_c13ul_global(iglob) * m2l_global(iglob) / &
                 (tti_c11l_global(iglob) * 4._CUSTOM_REAL) + &
                 tti_ec_c13_kl(i,j,ispec) * tti_c13ul_global(iglob) * m6l_global(iglob) / (tti_c13l_global(iglob) * 8._CUSTOM_REAL) - &
                 tti_ec_c15_kl(i,j,ispec) * tti_c13ul_global(iglob) * m4l_global(iglob) / (tti_c15l_global(iglob) * 2._CUSTOM_REAL) + &
                 tti_ec_c33_kl(i,j,ispec) * tti_c13ul_global(iglob) * m2l_global(iglob) / (tti_c33l_global(iglob) * 4._CUSTOM_REAL) + &
                 tti_ec_c35_kl(i,j,ispec) * tti_c13ul_global(iglob) * m4l_global(iglob) / (tti_c35l_global(iglob) * 2._CUSTOM_REAL) - &
                 tti_ec_c55_kl(i,j,ispec) * tti_c13ul_global(iglob) * m2l_global(iglob)/(tti_c55l_global(iglob) * 4._CUSTOM_REAL)

              tti_ecu_c33_kl(i,j,ispec) = tti_ec_c11_kl(i,j,ispec) * tti_c33ul_global(iglob) * m3l_global(iglob) / &
                 (tti_c11l_global(iglob) * 8._CUSTOM_REAL) + &
                 tti_ec_c13_kl(i,j,ispec) * tti_c33ul_global(iglob) * m2l_global(iglob) / (tti_c13l_global(iglob) * 8._CUSTOM_REAL) + &
                 tti_ec_c15_kl(i,j,ispec) * tti_c33ul_global(iglob) * m4l_global(iglob) / (tti_c15l_global(iglob) * 2._CUSTOM_REAL) - &
                 tti_ec_c15_kl(i,j,ispec) * tti_c33ul_global(iglob) * m5l_global(iglob) / (tti_c15l_global(iglob) * 4._CUSTOM_REAL) + &
                 tti_ec_c33_kl(i,j,ispec) * tti_c33ul_global(iglob) * m1l_global(iglob) / (tti_c33l_global(iglob) * 8._CUSTOM_REAL) - &
                 tti_ec_c35_kl(i,j,ispec) * tti_c33ul_global(iglob) * m5l_global(iglob) / (tti_c35l_global(iglob) * 4._CUSTOM_REAL) + &
                 tti_ec_c55_kl(i,j,ispec) * tti_c33ul_global(iglob) * m2l_global(iglob)/(tti_c55l_global(iglob) * 8._CUSTOM_REAL)

              tti_ecu_c55_kl(i,j,ispec) = tti_ec_c11_kl(i,j,ispec) * tti_c55ul_global(iglob) * m2l_global(iglob) / &
                 (tti_c11l_global(iglob) * 2._CUSTOM_REAL) - &
                 tti_ec_c13_kl(i,j,ispec) * tti_c55ul_global(iglob) * m2l_global(iglob) / (tti_c13l_global(iglob) * 2._CUSTOM_REAL) - &
                 tti_ec_c15_kl(i,j,ispec) * tti_c55ul_global(iglob) * m4l_global(iglob) / (tti_c15l_global(iglob) * 1._CUSTOM_REAL) + &
                 tti_ec_c33_kl(i,j,ispec) * tti_c55ul_global(iglob) * m2l_global(iglob) / (tti_c33l_global(iglob) * 2._CUSTOM_REAL) + &
                 tti_ec_c35_kl(i,j,ispec) * tti_c55ul_global(iglob) * m4l_global(iglob) / (tti_c35l_global(iglob) * 1._CUSTOM_REAL) + &
                 tti_ec_c55_kl(i,j,ispec) * tti_c55ul_global(iglob) /(tti_c55l_global(iglob) * 1._CUSTOM_REAL) - &
                 tti_ec_c55_kl(i,j,ispec) * tti_c55ul_global(iglob) * m2l_global(iglob)/(tti_c55l_global(iglob) * 2._CUSTOM_REAL)

              tti_ecu_theta_kl(i,j,ispec) = - tti_ec_c11_kl(i,j,ispec) * tti_thetaul_global(iglob) * p2l_global(iglob) * &
                 SIN(2._CUSTOM_REAL * tti_thetaul_global(iglob)) / (tti_c11l_global(iglob) * 1._CUSTOM_REAL) - &
                 tti_ec_c11_kl(i,j,ispec) * tti_thetaul_global(iglob) * p1l_global(iglob) * &
                 SIN(4._CUSTOM_REAL * tti_thetaul_global(iglob))/(tti_c11l_global(iglob) * 2._CUSTOM_REAL) + &
                 tti_ec_c13_kl(i,j,ispec) * tti_thetaul_global(iglob) * p1l_global(iglob) * &
                 SIN(4._CUSTOM_REAL * tti_thetaul_global(iglob))/(tti_c13l_global(iglob) * 2._CUSTOM_REAL) + &
                 tti_ec_c15_kl(i,j,ispec) * tti_thetaul_global(iglob) * p2l_global(iglob) * &
                 COS(2._CUSTOM_REAL * tti_thetaul_global(iglob))/(tti_c15l_global(iglob) * 1._CUSTOM_REAL) + &
                 tti_ec_c15_kl(i,j,ispec) * tti_thetaul_global(iglob) * p1l_global(iglob) * &
                 COS(4._CUSTOM_REAL * tti_thetaul_global(iglob))/(tti_c15l_global(iglob) * 2._CUSTOM_REAL) + &
                 tti_ec_c33_kl(i,j,ispec) * tti_thetaul_global(iglob) * p2l_global(iglob) * &
                 SIN(2._CUSTOM_REAL * tti_thetaul_global(iglob))/(tti_c33l_global(iglob) * 1._CUSTOM_REAL) - &
                 tti_ec_c33_kl(i,j,ispec) * tti_thetaul_global(iglob) * p1l_global(iglob) * &
                 SIN(4._CUSTOM_REAL * tti_thetaul_global(iglob))/(tti_c33l_global(iglob) * 2._CUSTOM_REAL) + &
                 tti_ec_c35_kl(i,j,ispec) * tti_thetaul_global(iglob) * p2l_global(iglob) * &
                 COS(2._CUSTOM_REAL * tti_thetaul_global(iglob))/(tti_c35l_global(iglob) * 1._CUSTOM_REAL) - &
                 tti_ec_c35_kl(i,j,ispec) * tti_thetaul_global(iglob) * p1l_global(iglob) * &
                 COS(4._CUSTOM_REAL * tti_thetaul_global(iglob))/(tti_c35l_global(iglob) * 2._CUSTOM_REAL) + &
                 tti_ec_c55_kl(i,j,ispec) * tti_thetaul_global(iglob) * p1l_global(iglob) * &
                 SIN(4._CUSTOM_REAL * tti_thetaul_global(iglob))/(tti_c55l_global(iglob) * 2._CUSTOM_REAL)


              tti_ecu_rho_kl(i,j,ispec) = tti_ec_rho_kl(i,j,ispec)
              !tti_ecu_c11_kl(i,j,ispec) = tti_ecu_c11_kl(i,j,ispec) - tti_c11l_global(iglob) * &
              !                      tti_ec_c11_k(iglob) * deltat

              !tti_ecu_c13_kl(i,j,ispec) = tti_ec_c13_kl(i,j,ispec)
              !tti_ecu_c33_kl(i,j,ispec) = tti_ec_c33_kl(i,j,ispec)
              !tti_ecu_c55_kl(i,j,ispec) = tti_ec_c55_kl(i,j,ispec)
              !tti_ecu_theta_kl(i,j,ispec) = tti_ec_c15_kl(i,j,ispec)
              !tti_ecu_rho_kl(i,j,ispec) = tti_ec_rho_kl(i,j,ispec)

              pdh_tti_ecu_c11(i,j,ispec) = pdh_tti_ecu_c11(i,j,ispec) + &
                                  pdh_tti_ec_c11_temp(iglob) * deltat
              pdh_tti_ecu_c13(i,j,ispec) = pdh_tti_ecu_c13(i,j,ispec) + &
                                  pdh_tti_ec_c13_temp(iglob) * deltat
              pdh_tti_ecu_c33(i,j,ispec) = pdh_tti_ecu_c33(i,j,ispec) + &
                                  pdh_tti_ec_c33_temp(iglob) * deltat
              pdh_tti_ecu_c55(i,j,ispec) = pdh_tti_ecu_c55(i,j,ispec) + &
                                  pdh_tti_ec_c55_temp(iglob) * deltat
              pdh_tti_ecu_theta(i,j,ispec) = pdh_tti_ecu_theta(i,j,ispec) + &
                                  pdh_tti_ec_c15_temp(iglob) * deltat
              pdh_tti_ecu_rho(i,j,ispec) = pdh_tti_ecu_rho(i,j,ispec) + &
                                  pdh_hti_ec_rho_temp(iglob) * deltat

             !!! TTITH model parameterization kernels
             ! tti_thom_rhop_kl(i,j,ispec) = tti_ecu_rho_kl(i,j,ispec)
             tti_thom_rhop_kl(i,j,ispec) = tti_ecu_rho_kl(i,j,ispec) + tti_alphal_global(iglob)*tti_alphal_global(iglob)*&
               tti_ecu_c33_kl(i,j,ispec) + tti_betal_global(iglob)*tti_betal_global(iglob)*tti_ecu_c55_kl(i,j,ispec) + &
               tti_coeffl_global(nglob)*tti_ecu_c13_kl(i,j,ispec)

             tti_thom_alpha_kl(i,j,ispec) = 2._CUSTOM_REAL*tti_ecu_c33_kl(i,j,ispec) + 2._CUSTOM_REAL*tti_ecu_c11_kl(i,j,ispec)+&
               2._CUSTOM_REAL*tti_alphal_global(iglob)*tti_alphal_global(iglob)*tti_alphal_global(iglob)*tti_alphal_global(iglob)*&
               (2._CUSTOM_REAL*tti_deltal_global(iglob)+1._CUSTOM_REAL)*tti_ecu_c13_kl(i,j,ispec)/(tti_coeffl_global(iglob)*&
               (tti_coeffl_global(iglob)-tti_betal_global(iglob)**2))-2._CUSTOM_REAL*tti_alphal_global(iglob)*tti_alphal_global(iglob)*&
               tti_betal_global(iglob)*tti_betal_global(iglob)*(tti_deltal_global(iglob)+1._CUSTOM_REAL)*tti_ecu_c13_kl(i,j,ispec)/&
               (tti_coeffl_global(iglob)*(tti_coeffl_global(iglob)-tti_betal_global(iglob)**2))

             tti_thom_beta_kl(i,j,ispec) = 2._CUSTOM_REAL*tti_ecu_c55_kl(i,j,ispec) - 2._CUSTOM_REAL * tti_betal_global(iglob) * &
               tti_betal_global(iglob)*(tti_alphal_global(iglob)*tti_alphal_global(iglob)*tti_deltal_global(iglob)-&
               tti_betal_global(iglob)**2 + tti_alphal_global(iglob)**2)*tti_ecu_c13_kl(i,j,ispec)/&
               (tti_coeffl_global(iglob)*(tti_coeffl_global(iglob)-tti_betal_global(iglob)**2)) - &
               2._CUSTOM_REAL*tti_betal_global(iglob)*tti_betal_global(iglob)*tti_ecu_c13_kl(i,j,ispec)/&
               (tti_coeffl_global(iglob)-tti_betal_global(iglob)**2)

             tti_thom_epsilon_kl(i,j,ispec) = 2._CUSTOM_REAL*tti_epsilonl_global(iglob) * tti_ecu_c11_kl(i,j,ispec)/&
               (2._CUSTOM_REAL*tti_epsilonl_global(iglob) + 1._CUSTOM_REAL)

             tti_thom_delta_kl(i,j,ispec) = tti_deltal_global(iglob) * tti_alphal_global(iglob)*tti_alphal_global(iglob) * &
               (tti_alphal_global(iglob)**2 - tti_betal_global(iglob)**2) * tti_ecu_c13_kl(i,j,ispec) / &
               (tti_coeffl_global(iglob)*(tti_coeffl_global(iglob)-tti_betal_global(iglob)**2))

             tti_thom_theta_kl(i,j,ispec) = tti_ecu_theta_kl(i,j,ispec)

             pdh_tti_thom_alpha(i,j,ispec) = pdh_hti_thom_alpha(i,j,ispec)
             pdh_tti_thom_beta(i,j,ispec) = pdh_hti_thom_beta(i,j,ispec)
             pdh_tti_thom_epsilon(i,j,ispec) = pdh_hti_thom_epsilon(i,j,ispec)
             pdh_tti_thom_delta(i,j,ispec) = pdh_hti_thom_delta(i,j,ispec)
             pdh_tti_thom_rhop(i,j,ispec) = pdh_hti_thom_rhop(i,j,ispec)
             pdh_tti_thom_theta(i,j,ispec) = pdh_tti_ecu_theta(i,j,ispec)
            endif 
            !tti_ecu_c11_kl(i,j,ispec) = &
            !   tti_ec_c11_kl(i,j,ispec) * tti_c11ul_global(iglob) * m1l_global(iglob)/(tti_c11l_global(iglob) * 8._CUSTOM_REAL) + &
            !   tti_ec_c13_kl(i,j,ispec) * tti_c11ul_global(iglob) * m2l_global(iglob)/(tti_c13l_global(iglob) * 8._CUSTOM_REAL) + &
            !   tti_ec_c15_kl(i,j,ispec) * tti_c11ul_global(iglob) * m5l_global(iglob)/(tti_c15l_global(iglob) * 4._CUSTOM_REAL) + &
            !   tti_ec_c33_kl(i,j,ispec) * tti_c11ul_global(iglob) * m3l_global(iglob)/(tti_c33l_global(iglob) * 8._CUSTOM_REAL) + &
            !   tti_ec_c35_kl(i,j,ispec) * tti_c11ul_global(iglob) * m5l_global(iglob)/(tti_c35l_global(iglob) * 4._CUSTOM_REAL) - &
            !   tti_ec_c35_kl(i,j,ispec) * tti_c11ul_global(iglob) * 2._CUSTOM_REAL * m4l_global(iglob) / &
            !   (tti_c35l_global(iglob) * 4._CUSTOM_REAL) + &
            !   tti_ec_c55_kl(i,j,ispec) * tti_c11ul_global(iglob) * m2l_global(iglob)/(tti_c55l_global(iglob) * 8._CUSTOM_REAL)
           endif
        enddo
      enddo
    endif
  enddo

  end subroutine compute_kernels_el

!
!-------------------------------------------------------------------------------------------------
!

  subroutine compute_kernels_ac()

! acoustic kernel calculations
! see e.g. Tromp et al. (2005)

  use specfem_par, only: nspec,acoustic,ibool,kappal_ac_global,rhol_ac_global,&
                         poroelastcoef,density,kmato,assign_external_model,rhoext,vpext,deltat,&
                         hprime_xx,hprime_zz,xix,xiz,gammax,gammaz,&
                         potential_acoustic,b_potential_acoustic,b_potential_dot_dot_acoustic,&
                         accel_ac,b_accel_ac,b_displ_ac,&
                         rho_ac_kl,kappa_ac_kl,rhop_ac_kl,alpha_ac_kl,rhorho_ac_hessian_final1,&
                         rhorho_ac_hessian_final2
  implicit none
  include "constants.h"

  !local variables
  integer :: i,j,k,ispec,iglob
  real(kind=CUSTOM_REAL) :: tempx1l,tempx2l,b_tempx1l,b_tempx2l,bb_tempx1l,bb_tempx2l
  double precision :: xixl,xizl,gammaxl,gammazl

  do ispec = 1, nspec
    if( acoustic(ispec) ) then
      do j = 1, NGLLZ
        do i = 1, NGLLX
          iglob = ibool(i,j,ispec)
          if( .not. assign_external_model ) then
            kappal_ac_global(iglob) = poroelastcoef(3,1,kmato(ispec))
            rhol_ac_global(iglob) = density(1,kmato(ispec))
          else
            kappal_ac_global(iglob) = rhoext(i,j,ispec)*vpext(i,j,ispec)*vpext(i,j,ispec)
            rhol_ac_global(iglob)   = rhoext(i,j,ispec)
          endif

! calcul the displacement by computing the gradient of potential / rho
! and calcul the acceleration by computing the gradient of potential_dot_dot / rho
          tempx1l = ZERO
          tempx2l = ZERO
          b_tempx1l = ZERO
          b_tempx2l = ZERO
          bb_tempx1l = ZERO
          bb_tempx2l = ZERO
          do k = 1,NGLLX
            ! derivative along x
            !tempx1l = tempx1l + potential_dot_dot_acoustic(ibool(k,j,ispec))*hprime_xx(i,k)
            tempx1l = tempx1l + potential_acoustic(ibool(k,j,ispec))*hprime_xx(i,k) !!! YANGL
            b_tempx1l = b_tempx1l + b_potential_acoustic(ibool(k,j,ispec))*hprime_xx(i,k)
            bb_tempx1l = bb_tempx1l + b_potential_dot_dot_acoustic(ibool(k,j,ispec))*hprime_xx(i,k)
            ! derivative along z
            !tempx2l = tempx2l + potential_dot_dot_acoustic(ibool(i,k,ispec))*hprime_zz(j,k)
            tempx2l = tempx2l + potential_acoustic(ibool(i,k,ispec))*hprime_zz(j,k) !!! YANGL
            b_tempx2l = b_tempx2l + b_potential_acoustic(ibool(i,k,ispec))*hprime_zz(j,k)
            bb_tempx2l = bb_tempx2l + b_potential_dot_dot_acoustic(ibool(i,k,ispec))*hprime_zz(j,k)
          enddo

          xixl = xix(i,j,ispec)
          xizl = xiz(i,j,ispec)
          gammaxl = gammax(i,j,ispec)
          gammazl = gammaz(i,j,ispec)

          ! derivatives of potential
          accel_ac(1,iglob) = (tempx1l*xixl + tempx2l*gammaxl) / rhol_ac_global(iglob)
          accel_ac(2,iglob) = (tempx1l*xizl + tempx2l*gammazl) / rhol_ac_global(iglob)
          b_displ_ac(1,iglob) = (b_tempx1l*xixl + b_tempx2l*gammaxl) / rhol_ac_global(iglob)
          b_displ_ac(2,iglob) = (b_tempx1l*xizl + b_tempx2l*gammazl) / rhol_ac_global(iglob)
          b_accel_ac(1,iglob) = (bb_tempx1l*xixl + bb_tempx2l*gammaxl) / rhol_ac_global(iglob)
          b_accel_ac(2,iglob) = (bb_tempx1l*xizl + bb_tempx2l*gammazl) / rhol_ac_global(iglob)
        enddo !i = 1, NGLLX
      enddo !j = 1, NGLLZ
    endif
  enddo

  do ispec = 1,nspec
    if( acoustic(ispec) ) then
      do j = 1, NGLLZ
        do i = 1, NGLLX
          iglob = ibool(i,j,ispec)
          !<YANGL
          !!!! old expression (from elastic kernels)
          !!!rho_ac_kl(i,j,ispec) = rho_ac_kl(i,j,ispec) - rhol_ac_global(iglob)  * &
          !!!      dot_product(accel_ac(:,iglob),b_displ_ac(:,iglob)) * deltat
          !!!kappa_ac_kl(i,j,ispec) = kappa_ac_kl(i,j,ispec) - kappal_ac_global(iglob) * &
          !!!      potential_dot_dot_acoustic(iglob)/kappal_ac_global(iglob) * &
          !!!      b_potential_dot_dot_acoustic(iglob)/kappal_ac_global(iglob)&
          !!!      * deltat
          !!!! new expression (from PDE-constrained optimization, coupling terms changed as well)
          rho_ac_kl(i,j,ispec) = rho_ac_kl(i,j,ispec) + rhol_ac_global(iglob) * &
                                 dot_product(accel_ac(:,iglob),b_displ_ac(:,iglob)) * deltat
          kappa_ac_kl(i,j,ispec) = kappa_ac_kl(i,j,ispec) + kappal_ac_global(iglob) * &
                                   potential_acoustic(iglob)/kappal_ac_global(iglob) * &
                                   b_potential_dot_dot_acoustic(iglob)/kappal_ac_global(iglob) * deltat
          !>YANGL
          rhop_ac_kl(i,j,ispec) = rho_ac_kl(i,j,ispec) + kappa_ac_kl(i,j,ispec)
          alpha_ac_kl(i,j,ispec) = TWO *  kappa_ac_kl(i,j,ispec)
          rhorho_ac_hessian_final1(i,j,ispec) =  rhorho_ac_hessian_final1(i,j,ispec) + &
                                                 dot_product(accel_ac(:,iglob),accel_ac(:,iglob)) * deltat
          rhorho_ac_hessian_final2(i,j,ispec) =  rhorho_ac_hessian_final2(i,j,ispec) + &
                                                 dot_product(accel_ac(:,iglob),b_accel_ac(:,iglob)) * deltat
        enddo
      enddo
    endif
  enddo

  end subroutine compute_kernels_ac

!
!-------------------------------------------------------------------------------------------------
!

  subroutine compute_kernels_po()

! kernel calculations
! see e.g. Morency et al. (2009)

  use specfem_par, only: nglob,nspec,poroelastic,ibool,deltat, &
                         kmato,porosity,tortuosity,density,permeability,poroelastcoef, &
                         ratio,B_biot,M_biot,C_biot,cpIsquare,cpIIsquare,cssquare, &
                         accels_poroelastic,accelw_poroelastic,velocw_poroelastic, &
                         b_displs_poroelastic,b_displw_poroelastic, &
                         rhot_k,rhof_k,sm_k,eta_k,B_k,C_k, &
                         phil_global,tortl_global,rhol_s_global,rhol_f_global,rhol_bar_global, &
                         etal_f_global,permlxx_global,permlxz_global,permlzz_global,mulfr_global, &
                         rhot_kl,rhof_kl,sm_kl,eta_kl,B_kl,C_kl,M_kl,M_k, &
                         mufr_kl,mufr_k,rhol_bar_global,rhob_kl,rhofb_kl,Bb_kl,Cb_kl,Mb_kl, &
                         mufrb_kl,phi_kl,rhobb_kl,rhofbb_kl,phib_kl,cpI_kl,cpII_kl,cs_kl,ratio_kl
  implicit none
  include "constants.h"

  !local variables
  integer :: i,j,ispec,iglob
  real(kind=CUSTOM_REAL) :: rholb,dd1

  do iglob =1,nglob
    rhot_k(iglob) = accels_poroelastic(1,iglob) * b_displs_poroelastic(1,iglob) + &
                    accels_poroelastic(2,iglob) * b_displs_poroelastic(2,iglob)
    rhof_k(iglob) = accelw_poroelastic(1,iglob) * b_displs_poroelastic(1,iglob) + &
                    accelw_poroelastic(2,iglob) * b_displs_poroelastic(2,iglob) + &
                    accels_poroelastic(1,iglob) * b_displw_poroelastic(1,iglob) + &
                    accels_poroelastic(2,iglob) * b_displw_poroelastic(2,iglob)
     sm_k(iglob)  = accelw_poroelastic(1,iglob) * b_displw_poroelastic(1,iglob) + &
                    accelw_poroelastic(2,iglob) * b_displw_poroelastic(2,iglob)
     eta_k(iglob) = velocw_poroelastic(1,iglob) * b_displw_poroelastic(1,iglob) + &
                    velocw_poroelastic(2,iglob) * b_displw_poroelastic(2,iglob)
  enddo

  do ispec = 1, nspec
    if( poroelastic(ispec) ) then
      do j = 1, NGLLZ
        do i = 1, NGLLX
          iglob = ibool(i,j,ispec)
          phil_global(iglob) = porosity(kmato(ispec))
          tortl_global(iglob) = tortuosity(kmato(ispec))
          rhol_s_global(iglob) = density(1,kmato(ispec))
          rhol_f_global(iglob) = density(2,kmato(ispec))
          rhol_bar_global(iglob) =  (1._CUSTOM_REAL - phil_global(iglob))*rhol_s_global(iglob) + &
                                    phil_global(iglob)*rhol_f_global(iglob)
          etal_f_global(iglob) = poroelastcoef(2,2,kmato(ispec))
          permlxx_global(iglob) = permeability(1,kmato(ispec))
          permlxz_global(iglob) = permeability(2,kmato(ispec))
          permlzz_global(iglob) = permeability(3,kmato(ispec))
          mulfr_global(iglob) = poroelastcoef(2,3,kmato(ispec))

          rhot_kl(i,j,ispec) = rhot_kl(i,j,ispec) - deltat * rhol_bar_global(iglob) * rhot_k(iglob)
          rhof_kl(i,j,ispec) = rhof_kl(i,j,ispec) - deltat * rhol_f_global(iglob) * rhof_k(iglob)
          sm_kl(i,j,ispec) = sm_kl(i,j,ispec) - &
                             deltat * rhol_f_global(iglob)*tortl_global(iglob)/phil_global(iglob) * sm_k(iglob)
          !at the moment works with constant permeability
          eta_kl(i,j,ispec) = eta_kl(i,j,ispec) - deltat * etal_f_global(iglob)/permlxx_global(iglob) * eta_k(iglob)
          B_kl(i,j,ispec) = B_kl(i,j,ispec) - deltat * B_k(iglob)
          C_kl(i,j,ispec) = C_kl(i,j,ispec) - deltat * C_k(iglob)
          M_kl(i,j,ispec) = M_kl(i,j,ispec) - deltat * M_k(iglob)
          mufr_kl(i,j,ispec) = mufr_kl(i,j,ispec) - TWO * deltat * mufr_k(iglob)
          ! density kernels
          rholb = rhol_bar_global(iglob) - phil_global(iglob)*rhol_f_global(iglob)/tortl_global(iglob)
          rhob_kl(i,j,ispec) = rhot_kl(i,j,ispec) + B_kl(i,j,ispec) + mufr_kl(i,j,ispec)
          rhofb_kl(i,j,ispec) = rhof_kl(i,j,ispec) + C_kl(i,j,ispec) + M_kl(i,j,ispec) + sm_kl(i,j,ispec)
          Bb_kl(i,j,ispec) = B_kl(i,j,ispec)
          Cb_kl(i,j,ispec) = C_kl(i,j,ispec)
          Mb_kl(i,j,ispec) = M_kl(i,j,ispec)
          mufrb_kl(i,j,ispec) = mufr_kl(i,j,ispec)
          phi_kl(i,j,ispec) = - sm_kl(i,j,ispec) - M_kl(i,j,ispec)
          ! wave speed kernels
          dd1 = (1._CUSTOM_REAL+rholb/rhol_f_global(iglob))*ratio**2 + &
                2._CUSTOM_REAL*ratio + tortl_global(iglob)/phil_global(iglob)

          rhobb_kl(i,j,ispec) = rhob_kl(i,j,ispec) - &
                phil_global(iglob)*rhol_f_global(iglob)/(tortl_global(iglob)*B_biot) * &
                (cpIIsquare + (cpIsquare - cpIIsquare)*( (phil_global(iglob) / &
                tortl_global(iglob)*ratio +1._CUSTOM_REAL)/dd1 + &
                (rhol_bar_global(iglob)**2*ratio**2/rhol_f_global(iglob)**2*(phil_global(iglob) / &
                tortl_global(iglob)*ratio+1._CUSTOM_REAL)*(phil_global(iglob)/tortl_global(iglob)*ratio + &
                phil_global(iglob)/tortl_global(iglob) * &
                (1._CUSTOM_REAL+rhol_f_global(iglob)/rhol_bar_global(iglob))-1._CUSTOM_REAL) )/dd1**2 ) - &
                FOUR_THIRDS*cssquare )*Bb_kl(i,j,ispec) - &
                rhol_bar_global(iglob)*ratio**2/M_biot * (cpIsquare - cpIIsquare)* &
                (phil_global(iglob)/tortl_global(iglob)*ratio + &
                1._CUSTOM_REAL)**2/dd1**2*Mb_kl(i,j,ispec) + &
                rhol_bar_global(iglob)*ratio/C_biot * (cpIsquare - cpIIsquare)* (&
                (phil_global(iglob)/tortl_global(iglob)*ratio+1._CUSTOM_REAL)/dd1 - &
                phil_global(iglob)*ratio/tortl_global(iglob)*(phil_global(iglob) / &
                tortl_global(iglob)*ratio+1._CUSTOM_REAL)*&
                (1._CUSTOM_REAL+rhol_bar_global(iglob)*ratio/rhol_f_global(iglob))/dd1**2)*Cb_kl(i,j,ispec)+ &
                phil_global(iglob)*rhol_f_global(iglob)*cssquare / &
                (tortl_global(iglob)*mulfr_global(iglob))*mufrb_kl(i,j,ispec)
          rhofbb_kl(i,j,ispec) = rhofb_kl(i,j,ispec) + &
                 phil_global(iglob)*rhol_f_global(iglob)/(tortl_global(iglob)*B_biot) * &
                 (cpIIsquare + (cpIsquare - cpIIsquare)*( (phil_global(iglob)/ &
                 tortl_global(iglob)*ratio +1._CUSTOM_REAL)/dd1+&
                 (rhol_bar_global(iglob)**2*ratio**2/rhol_f_global(iglob)**2*(phil_global(iglob)/ &
                 tortl_global(iglob)*ratio+1)*(phil_global(iglob)/tortl_global(iglob)*ratio+ &
                 phil_global(iglob)/tortl_global(iglob)*&
                 (1._CUSTOM_REAL+rhol_f_global(iglob)/rhol_bar_global(iglob))-1._CUSTOM_REAL) )/dd1**2 )- &
                 FOUR_THIRDS*cssquare )*Bb_kl(i,j,ispec) + &
                 rhol_bar_global(iglob)*ratio**2/M_biot * (cpIsquare - cpIIsquare)* &
                 (phil_global(iglob)/tortl_global(iglob)*ratio + &
                 1._CUSTOM_REAL)**2/dd1**2*Mb_kl(i,j,ispec) - &
                 rhol_bar_global(iglob)*ratio/C_biot * (cpIsquare - cpIIsquare)* (&
                 (phil_global(iglob)/tortl_global(iglob)*ratio+1._CUSTOM_REAL)/dd1 - &
                 phil_global(iglob)*ratio/tortl_global(iglob)*(phil_global(iglob)/ &
                 tortl_global(iglob)*ratio+1._CUSTOM_REAL)*&
                 (1._CUSTOM_REAL+rhol_bar_global(iglob)*ratio/rhol_f_global(iglob))/dd1**2)*Cb_kl(i,j,ispec)- &
                 phil_global(iglob)*rhol_f_global(iglob)*cssquare/ &
                 (tortl_global(iglob)*mulfr_global(iglob))*mufrb_kl(i,j,ispec)
          phib_kl(i,j,ispec) = phi_kl(i,j,ispec) - &
                 phil_global(iglob)*rhol_bar_global(iglob)/(tortl_global(iglob)*B_biot) &
                 * ( cpIsquare - rhol_f_global(iglob)/rhol_bar_global(iglob)*cpIIsquare- &
                 (cpIsquare-cpIIsquare)*( (TWO*ratio**2*phil_global(iglob)/ &
                 tortl_global(iglob) + (1._CUSTOM_REAL+&
                 rhol_f_global(iglob)/rhol_bar_global(iglob))* &
                 (TWO*ratio*phil_global(iglob)/tortl_global(iglob)+&
                 1._CUSTOM_REAL))/dd1 + (phil_global(iglob)/tortl_global(iglob)*ratio+ &
                 1._CUSTOM_REAL)*(phil_global(iglob)*&
                 ratio/tortl_global(iglob)+phil_global(iglob)/tortl_global(iglob)* &
                 (1._CUSTOM_REAL+rhol_f_global(iglob)/&
                 rhol_bar_global(iglob))-1._CUSTOM_REAL)*((1._CUSTOM_REAL+ &
                 rhol_bar_global(iglob)/rhol_f_global(iglob)-&
                 TWO*phil_global(iglob)/tortl_global(iglob))*ratio**2+TWO*ratio)/dd1**2 ) - &
                 FOUR_THIRDS*rhol_f_global(iglob)*cssquare/rhol_bar_global(iglob) )*Bb_kl(i,j,ispec) + &
                 rhol_f_global(iglob)/M_biot * (cpIsquare-cpIIsquare)*(&
                 TWO*ratio*(phil_global(iglob)/tortl_global(iglob)*ratio+1._CUSTOM_REAL)/dd1 - &
                 (phil_global(iglob)/tortl_global(iglob)*ratio+1._CUSTOM_REAL)**2*( &
                 (1._CUSTOM_REAL+rhol_bar_global(iglob)/&
                 rhol_f_global(iglob)-TWO*phil_global(iglob)/tortl_global(iglob))*ratio**2+TWO*ratio)/dd1**2 &
                 )*Mb_kl(i,j,ispec) + &
                 phil_global(iglob)*rhol_f_global(iglob)/(tortl_global(iglob)*C_biot)* &
                 (cpIsquare-cpIIsquare)*ratio* (&
                 (1._CUSTOM_REAL+rhol_f_global(iglob)/rhol_bar_global(iglob)*ratio)/dd1 - &
                 (phil_global(iglob)/tortl_global(iglob)*ratio+1._CUSTOM_REAL)* &
                 (1._CUSTOM_REAL+rhol_bar_global(iglob)/&
                 rhol_f_global(iglob)*ratio)*((1._CUSTOM_REAL+rhol_bar_global(iglob)/rhol_f_global(iglob)-TWO*&
                 phil_global(iglob)/tortl_global(iglob))*ratio+TWO)/dd1**2&
                 )*Cb_kl(i,j,ispec) -&
                 phil_global(iglob)*rhol_f_global(iglob)*cssquare &
                 /(tortl_global(iglob)*mulfr_global(iglob))*mufrb_kl(i,j,ispec)
          cpI_kl(i,j,ispec) = 2._CUSTOM_REAL*cpIsquare/B_biot*rhol_bar_global(iglob)*( &
                 1._CUSTOM_REAL-phil_global(iglob)/tortl_global(iglob) + &
                 (phil_global(iglob)/tortl_global(iglob)*ratio+ &
                 1._CUSTOM_REAL)*(phil_global(iglob)/tortl_global(iglob)*&
                 ratio+phil_global(iglob)/tortl_global(iglob)* &
                 (1._CUSTOM_REAL+rhol_f_global(iglob)/rhol_bar_global(iglob))-&
                 1._CUSTOM_REAL)/dd1 &
                 )* Bb_kl(i,j,ispec) +&
                 2._CUSTOM_REAL*cpIsquare*rhol_f_global(iglob)*tortl_global(iglob)/(phil_global(iglob)*M_biot) *&
                 (phil_global(iglob)/tortl_global(iglob)*ratio+1._CUSTOM_REAL)**2/dd1*Mb_kl(i,j,ispec)+&
                 2._CUSTOM_REAL*cpIsquare*rhol_f_global(iglob)/C_biot * &
                 (phil_global(iglob)/tortl_global(iglob)*ratio+1._CUSTOM_REAL)* &
                 (1._CUSTOM_REAL+rhol_bar_global(iglob)/&
                 rhol_f_global(iglob)*ratio)/dd1*Cb_kl(i,j,ispec)
          cpII_kl(i,j,ispec) = 2._CUSTOM_REAL*cpIIsquare*rhol_bar_global(iglob)/B_biot * (&
                 phil_global(iglob)*rhol_f_global(iglob)/(tortl_global(iglob)*rhol_bar_global(iglob)) - &
                 (phil_global(iglob)/tortl_global(iglob)*ratio+ &
                 1._CUSTOM_REAL)*(phil_global(iglob)/tortl_global(iglob)*&
                 ratio+phil_global(iglob)/tortl_global(iglob)* &
                 (1._CUSTOM_REAL+rhol_f_global(iglob)/rhol_bar_global(iglob))-&
                 1._CUSTOM_REAL)/dd1  ) * Bb_kl(i,j,ispec) +&
                 2._CUSTOM_REAL*cpIIsquare*rhol_f_global(iglob)*tortl_global(iglob)/(phil_global(iglob)*M_biot) * (&
                 1._CUSTOM_REAL - (phil_global(iglob)/tortl_global(iglob)*ratio+ &
                 1._CUSTOM_REAL)**2/dd1  )*Mb_kl(i,j,ispec) + &
                 2._CUSTOM_REAL*cpIIsquare*rhol_f_global(iglob)/C_biot * (&
                 1._CUSTOM_REAL - (phil_global(iglob)/tortl_global(iglob)*ratio+ &
                 1._CUSTOM_REAL)*(1._CUSTOM_REAL+&
                 rhol_bar_global(iglob)/rhol_f_global(iglob)*ratio)/dd1  )*Cb_kl(i,j,ispec)
          cs_kl(i,j,ispec) = - 8._CUSTOM_REAL/3._CUSTOM_REAL*cssquare* &
                 rhol_bar_global(iglob)/B_biot*(1._CUSTOM_REAL-&
                 phil_global(iglob)*rhol_f_global(iglob)/(tortl_global(iglob)* &
                 rhol_bar_global(iglob)))*Bb_kl(i,j,ispec) + &
                 2._CUSTOM_REAL*(rhol_bar_global(iglob)-rhol_f_global(iglob)*&
                 phil_global(iglob)/tortl_global(iglob))/&
                 mulfr_global(iglob)*cssquare*mufrb_kl(i,j,ispec)
          ratio_kl(i,j,ispec) = ratio*rhol_bar_global(iglob)*phil_global(iglob)/(tortl_global(iglob)*B_biot) * &
                 (cpIsquare-cpIIsquare) * ( &
                 phil_global(iglob)/tortl_global(iglob)*(2._CUSTOM_REAL*ratio+1._CUSTOM_REAL+rhol_f_global(iglob)/ &
                 rhol_bar_global(iglob))/dd1 - (phil_global(iglob)/tortl_global(iglob)*ratio+1._CUSTOM_REAL)*&
                 (phil_global(iglob)/tortl_global(iglob)*ratio+phil_global(iglob)/tortl_global(iglob)*(&
                 1._CUSTOM_REAL+rhol_f_global(iglob)/rhol_bar_global(iglob))-1._CUSTOM_REAL)*(2._CUSTOM_REAL*ratio*(&
                 1._CUSTOM_REAL+rhol_bar_global(iglob)/rhol_f_global(iglob)-phil_global(iglob)/tortl_global(iglob)) +&
                 2._CUSTOM_REAL)/dd1**2  )*Bb_kl(i,j,ispec) + &
                 ratio*rhol_f_global(iglob)*tortl_global(iglob)/(phil_global(iglob)*M_biot)*(cpIsquare-cpIIsquare) * &
                 2._CUSTOM_REAL*phil_global(iglob)/tortl_global(iglob) * (&
                 (phil_global(iglob)/tortl_global(iglob)*ratio+1._CUSTOM_REAL)/dd1 - &
                 (phil_global(iglob)/tortl_global(iglob)*ratio+1._CUSTOM_REAL)**2*( &
                 (1._CUSTOM_REAL+rhol_bar_global(iglob)/&
                 rhol_f_global(iglob)-phil_global(iglob)/tortl_global(iglob))*ratio+ &
                 1._CUSTOM_REAL)/dd1**2 )*Mb_kl(i,j,ispec) +&
                 ratio*rhol_f_global(iglob)/C_biot*(cpIsquare-cpIIsquare) * (&
                 (2._CUSTOM_REAL*phil_global(iglob)*rhol_bar_global(iglob)* &
                 ratio/(tortl_global(iglob)*rhol_f_global(iglob))+&
                 phil_global(iglob)/tortl_global(iglob)+rhol_bar_global(iglob)/rhol_f_global(iglob))/dd1 - &
                 2._CUSTOM_REAL*phil_global(iglob)/tortl_global(iglob)*(phil_global(iglob)/tortl_global(iglob)*ratio+&
                 1._CUSTOM_REAL)*(1._CUSTOM_REAL+rhol_bar_global(iglob)/rhol_f_global(iglob)*ratio)*((1._CUSTOM_REAL+&
                 rhol_bar_global(iglob)/rhol_f_global(iglob)- &
                 phil_global(iglob)/tortl_global(iglob))*ratio+1._CUSTOM_REAL)/&
                 dd1**2 )*Cb_kl(i,j,ispec)
        enddo
      enddo
    endif
  enddo

  end subroutine compute_kernels_po


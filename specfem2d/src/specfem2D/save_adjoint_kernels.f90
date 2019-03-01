
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
!========================================================================

subroutine save_adjoint_kernels()

  use specfem_par, only : myrank, nspec, ibool, coord, save_ASCII_kernels, &
                          any_acoustic, any_elastic, any_poroelastic, &
                          rho_ac_kl, kappa_ac_kl, alpha_ac_kl, rhop_ac_kl, &
                          rho_kl, kappa_kl, mu_kl, rhop_kl, alpha_kl, beta_kl, &
                          bulk_c_kl, bulk_beta_kl, &
                          rhorho_ac_hessian_final1, rhorho_ac_hessian_final2, &
                          rhorho_el_hessian_final1, rhorho_el_hessian_final2, &
                          rhot_kl, rhof_kl, sm_kl, eta_kl, mufr_kl, B_kl, &
                          C_kl, M_kl, rhob_kl, rhofb_kl, phi_kl, Bb_kl, Cb_kl, Mb_kl, mufrb_kl, &
                          rhobb_kl, rhofbb_kl, phib_kl, cpI_kl, cpII_kl, cs_kl, ratio_kl, GPU_MODE, &
                          hti_ec_rho_kl,hti_ec_c11_kl,hti_ec_c13_kl,hti_ec_c33_kl,hti_ec_c55_kl, &
                          hti_thom_rhop_kl,hti_thom_alpha_kl,hti_thom_beta_kl,&
                          hti_thom_epsilon_kl,hti_thom_delta_kl,ANISO, M_PAR, &
                          pdh_hti_ec_c11, pdh_hti_ec_c13, pdh_hti_ec_c33, pdh_hti_ec_c55, pdh_hti_ec_rho, &
                          pdh_hti_thom_alpha, pdh_hti_thom_beta, pdh_hti_thom_epsilon, pdh_hti_thom_delta, &
                          pdh_hti_thom_rhop,&
                          Qalpha_kl,Qbeta_kl,Qkappa_kl,Qmu_kl, &
                          hti_thom2_rhop_kl,hti_thom2_alpha_kl,hti_thom2_beta_kl,&
                          hti_thom2_epsilon_kl,hti_thom2_delta_kl,hti_thom2_eta_kl,&
                          pdh_hti_thom2_alpha, pdh_hti_thom2_beta, pdh_hti_thom2_epsilon, pdh_hti_thom2_delta, &
                          pdh_hti_thom2_rhop,pdh_hti_thom2_eta, &
                          hti_vel_rhop_kl,hti_vel_alpha_kl,hti_vel_beta_kl,hti_vel_alphah_kl,hti_vel_alphan_kl,&
                          pdh_hti_vel_alpha, pdh_hti_vel_beta, pdh_hti_vel_alphah, pdh_hti_vel_alphan,&
                          pdh_hti_vel_rhop,&
                          tti_ec_rho_kl,tti_ec_c11_kl,tti_ec_c13_kl,tti_ec_c15_kl,tti_ec_c33_kl,tti_ec_c35_kl,tti_ec_c55_kl,&
                          pdh_tti_ec_c11, pdh_tti_ec_c13, pdh_tti_ec_c33,&
                          pdh_tti_ec_c15, pdh_tti_ec_c55, pdh_tti_ec_c35, pdh_tti_ec_rho, &
                          dl_lambda_kl, dl_mu_kl, dl_rho_kl, dip_ip_kl, dip_is_kl, dip_rho_kl, &
                          vipi_vp_kl, vipi_vs_kl, vipi_ip_kl, vipii_vp_kl, vipii_vs_kl, vipii_is_kl, &
                          sigma_vp_kl, sigma_sigma_kl, sigma_rho_kl, &
                          pdh_vp, pdh_vs, pdh_rhop, pdh_kappa, pdh_mu, pdh_rho, pdh_dl_lambda, pdh_dl_mu, pdh_dl_rho,&
                          pdh_dip_ip, pdh_dip_is, pdh_dip_rho, pdh_vipi_vp, pdh_vipi_vs, pdh_vipi_ip, &
                          pdh_vipii_vp, pdh_vipii_vs, pdh_vipii_is, pdh_sigma_vp, pdh_sigma_sigma, pdh_sigma_rho, &
                          tti_ecu_c11_kl,tti_ecu_c13_kl,tti_ecu_c33_kl,tti_ecu_c55_kl,tti_ecu_theta_kl,tti_ecu_rho_kl,&
                          pdh_tti_ecu_c11, pdh_tti_ecu_c13, pdh_tti_ecu_c33, pdh_tti_ecu_c55,pdh_tti_ecu_theta,pdh_tti_ecu_rho,&
                          tti_thom_alpha_kl,tti_thom_beta_kl,tti_thom_epsilon_kl,tti_thom_delta_kl,tti_thom_rhop_kl,tti_thom_theta_kl,&
                          pdh_tti_thom_alpha,pdh_tti_thom_beta,pdh_tti_thom_epsilon,pdh_tti_thom_delta,&
                          pdh_tti_thom_theta,pdh_tti_thom_rhop

  include "constants.h"

  integer :: i, j, ispec, iglob
  double precision :: xx, zz

  if ( myrank == 0 ) then
    write(IOUT,*) 'Writing Kernels file'
  endif

  if(any_acoustic) then
    if(save_ASCII_kernels) then ! ascii format
      do ispec = 1, nspec
        do j = 1, NGLLZ
          do i = 1, NGLLX
            iglob = ibool(i,j,ispec)
            xx = coord(1,iglob)
            zz = coord(2,iglob)
            write(95,'(4e15.5e4)')xx,zz,rho_ac_kl(i,j,ispec),kappa_ac_kl(i,j,ispec)
            write(96,'(4e15.5e4)')xx,zz,rhop_ac_kl(i,j,ispec),alpha_ac_kl(i,j,ispec)
            !write(96,'(4e15.5e4)')rhorho_ac_hessian_final1(i,j,ispec),
            !rhorho_ac_hessian_final2(i,j,ispec),&
            !                rhop_ac_kl(i,j,ispec),alpha_ac_kl(i,j,ispec)
          enddo
        enddo
      enddo
      close(95)
      close(96)

    else ! binary format
       write(200)rho_ac_kl
       write(201)kappa_ac_kl
       write(202)rhop_ac_kl
       write(203)alpha_ac_kl
       close(200)
       close(201)
       close(202)
       close(203)

      if (SAVE_DIAGONAL_HESSIAN) then
        write(212)rhorho_ac_hessian_final1
        write(213)rhorho_ac_hessian_final2
        close(212)
        close(213)
      endif

    endif
  endif

  if(any_elastic) then
    if(save_ASCII_kernels)then ! ascii format
    do ispec = 1, nspec
        do j = 1, NGLLZ
          do i = 1, NGLLX
            iglob = ibool(i,j,ispec)
            xx = coord(1,iglob)
            zz = coord(2,iglob)
            write(97,'(5e15.5e4)')xx,zz,rho_kl(i,j,ispec),kappa_kl(i,j,ispec),mu_kl(i,j,ispec)
            write(98,'(5e15.5e4)')xx,zz,rhop_kl(i,j,ispec),alpha_kl(i,j,ispec),beta_kl(i,j,ispec)
            !write(98,'(5e15.5e4)')rhorho_el_hessian_final1(i,j,ispec),
            !rhorho_el_hessian_final2(i,j,ispec),&
            !rhop_kl(i,j,ispec),alpha_kl(i,j,ispec),beta_kl(i,j,ispec)
          enddo
        enddo
      enddo
      close(97)
      close(98)
    else if (ANISO) then
     if ((trim(M_PAR) == 'htiec' .OR. trim(M_PAR) == 'vtiec')) then
      write(216)hti_ec_rho_kl
      write(217)hti_ec_c11_kl
      write(218)hti_ec_c13_kl
      write(219)hti_ec_c33_kl
      write(220)hti_ec_c55_kl
      close(216)
      close(217)
      close(218)
      close(219)
      close(220)
      write(226)pdh_hti_ec_c11
      write(227)pdh_hti_ec_c13
      write(228)pdh_hti_ec_c33
      write(229)pdh_hti_ec_c55
      write(230)pdh_hti_ec_rho
      close(226)
      close(227)
      close(228)
      close(229)
      close(230)
     endif
     if ((trim(M_PAR) == 'ttiec')) then
      write(258)tti_ec_rho_kl
      write(259)tti_ec_c11_kl
      write(260)tti_ec_c13_kl
      write(261)tti_ec_c15_kl
      write(262)tti_ec_c33_kl
      write(263)tti_ec_c35_kl
      write(264)tti_ec_c55_kl
      close(258)
      close(259)
      close(260)
      close(261)
      close(262)
      close(263)
      close(264)
      write(265)pdh_tti_ec_rho
      write(266)pdh_tti_ec_c11
      write(267)pdh_tti_ec_c13
      write(268)pdh_tti_ec_c15
      write(269)pdh_tti_ec_c33
      write(270)pdh_tti_ec_c35
      write(271)pdh_tti_ec_c55
      close(265)
      close(266)
      close(267)
      close(268)
      close(269)
      close(270)
      close(271)
     endif
     if ((trim(M_PAR) == 'ttiecu')) then
      write(308)tti_ecu_c11_kl
      write(309)tti_ecu_c13_kl
      write(310)tti_ecu_c33_kl
      write(311)tti_ecu_c55_kl
      write(312)tti_ecu_theta_kl
      write(313)tti_ecu_rho_kl
      write(314)pdh_tti_ecu_c11
      write(315)pdh_tti_ecu_c13
      write(316)pdh_tti_ecu_c33
      write(317)pdh_tti_ecu_c55
      write(318)pdh_tti_ecu_theta
      write(319)pdh_tti_ecu_rho
      close(308)
      close(309)
      close(310)
      close(311)
      close(312)
      close(313)
      close(314)
      close(315)
      close(316)
      close(317)
      close(318)
      close(319)
     endif
     if ((trim(M_PAR) == 'ttithom')) then
      write(320)tti_thom_alpha_kl
      write(321)tti_thom_beta_kl
      write(322)tti_thom_epsilon_kl
      write(323)tti_thom_delta_kl
      write(324)tti_thom_theta_kl
      write(325)tti_thom_rhop_kl
      write(326)pdh_tti_thom_alpha
      write(327)pdh_tti_thom_beta
      write(328)pdh_tti_thom_epsilon
      write(329)pdh_tti_thom_delta
      write(330)pdh_tti_thom_theta
      write(331)pdh_tti_thom_rhop
      close(320)
      close(321)
      close(322)
      close(323)
      close(324)
      close(325)
      close(326)
      close(327)
      close(328)
      close(329)
      close(330)
      close(331)
     endif
     if ((trim(M_PAR) == 'htithom' .OR. trim(M_PAR) == 'vtithom')) then
      write(221)hti_thom_rhop_kl
      write(222)hti_thom_alpha_kl
      write(223)hti_thom_beta_kl
      write(224)hti_thom_epsilon_kl
      write(225)hti_thom_delta_kl
      close(221)
      close(222)
      close(223)
      close(224)
      close(225)
      write(231)pdh_hti_thom_alpha
      write(232)pdh_hti_thom_beta
      write(233)pdh_hti_thom_epsilon
      write(234)pdh_hti_thom_delta
      write(235)pdh_hti_thom_rhop
      close(231)
      close(232)
      close(233)
      close(234)
      close(235)
     endif
     if ((trim(M_PAR) == 'htithom2' .OR. trim(M_PAR) == 'vtithom2')) then
      write(236)hti_thom2_rhop_kl
      write(237)hti_thom2_alpha_kl
      write(238)hti_thom2_beta_kl
      write(239)hti_thom2_epsilon_kl
      write(240)hti_thom2_delta_kl
      write(241)hti_thom2_eta_kl
      close(236)
      close(237)
      close(238)
      close(239)
      close(240)
      close(241)
      write(242)pdh_hti_thom2_alpha
      write(243)pdh_hti_thom2_beta
      write(244)pdh_hti_thom2_epsilon
      write(245)pdh_hti_thom2_delta
      write(246)pdh_hti_thom2_rhop
      write(247)pdh_hti_thom2_eta
      close(242)
      close(243)
      close(244)
      close(245)
      close(246)
      close(247)
     endif
     if ((trim(M_PAR) == 'htivel' .OR. trim(M_PAR) == 'vtivel')) then
      write(248)hti_vel_rhop_kl
      write(249)hti_vel_alpha_kl
      write(250)hti_vel_beta_kl
      write(251)hti_vel_alphah_kl
      write(252)hti_vel_alphan_kl
      close(248)
      close(249)
      close(250)
      close(251)
      close(252)
      write(253)pdh_hti_vel_alpha
      write(254)pdh_hti_vel_beta
      write(255)pdh_hti_vel_alphah
      write(256)pdh_hti_vel_alphan
      write(257)pdh_hti_vel_rhop
      close(253)
      close(254)
      close(255)
      close(256)
      close(257)
     endif
    else if (.not. ANISO) then! binary format
     if ((trim(M_PAR) == 'isodv' .OR. trim(M_PAR) == 'isodm')) then
      write(207)rhop_kl
      write(208)alpha_kl
      write(209)beta_kl
      close(207)
      close(208)
      close(209)
      write(204)rho_kl
      write(205)kappa_kl
      write(206)mu_kl
      close(204)
      close(205)
      close(206)
      write(287)pdh_vp
      write(288)pdh_vs
      write(289)pdh_rhop
      close(287)
      close(288)
      close(289)
      write(290)pdh_kappa
      write(291)pdh_mu
      write(292)pdh_rho
      close(290)
      close(291)
      close(292)
     endif
     if ((trim(M_PAR) == 'isodl')) then
      write(272)dl_lambda_kl
      write(273)dl_mu_kl
      write(274)dl_rho_kl
      close(272)
      close(273)
      close(274)
      write(293)pdh_dl_lambda
      write(294)pdh_dl_mu
      write(295)pdh_dl_rho
      close(293)
      close(294)
      close(295)
     endif
     if ((trim(M_PAR) == 'isodip')) then
      write(275)dip_ip_kl
      write(276)dip_is_kl
      write(277)dip_rho_kl
      close(275)
      close(276)
      close(277)
      write(296)pdh_dip_ip
      write(297)pdh_dip_is
      write(298)pdh_dip_rho
      close(296)
      close(297)
      close(298)
     endif
     if ((trim(M_PAR) == 'isovipi')) then
      write(278)vipi_vp_kl
      write(279)vipi_vs_kl
      write(280)vipi_ip_kl
      close(278)
      close(279)
      close(280)
      write(299)pdh_vipi_vp
      write(300)pdh_vipi_vs
      write(301)pdh_vipi_ip
      close(299)
      close(300)
      close(301)
     endif
     if ((trim(M_PAR) == 'isovipii')) then
      write(281)vipii_vp_kl
      write(282)vipii_vs_kl
      write(283)vipii_is_kl
      close(281)
      close(282)
      close(283)
      write(302)pdh_vipii_vp
      write(303)pdh_vipii_vs
      write(304)pdh_vipii_is
      close(302)
      close(303)
      close(304)
     endif
     if ((trim(M_PAR) == 'isosigma')) then
      write(284)sigma_vp_kl
      write(285)sigma_sigma_kl
      write(286)sigma_rho_kl
      close(284)
      close(285)
      close(286)
      write(305)pdh_sigma_vp
      write(306)pdh_sigma_sigma
      write(307)pdh_sigma_rho
      close(305)
      close(306)
      close(307)
     endif
      write(210)Qkappa_kl
      write(211)Qmu_kl
      close(210)
      close(211)

      if (SAVE_DIAGONAL_HESSIAN) then
        write(214)rhorho_el_hessian_final1
        write(215)rhorho_el_hessian_final2
        close(214)
        close(215)
      endif

    endif
  endif

if (.not. GPU_MODE )  then

  if(any_poroelastic) then

      if (.not. SAVE_ASCII_KERNELS) stop 'poroelastic simulations must use SAVE_ASCII_KERNELS'

    do ispec = 1, nspec
      do j = 1, NGLLZ
        do i = 1, NGLLX
          iglob = ibool(i,j,ispec)
          xx = coord(1,iglob)
          zz = coord(2,iglob)
          write(144,'(5e11.3)')xx,zz,mufr_kl(i,j,ispec),B_kl(i,j,ispec),C_kl(i,j,ispec)
          write(155,'(5e11.3)')xx,zz,M_kl(i,j,ispec),rhot_kl(i,j,ispec),rhof_kl(i,j,ispec)
          write(16,'(5e11.3)')xx,zz,sm_kl(i,j,ispec),eta_kl(i,j,ispec)
          write(17,'(5e11.3)')xx,zz,mufrb_kl(i,j,ispec),Bb_kl(i,j,ispec),Cb_kl(i,j,ispec)
          write(18,'(5e11.3)')xx,zz,Mb_kl(i,j,ispec),rhob_kl(i,j,ispec),rhofb_kl(i,j,ispec)
          write(19,'(5e11.3)')xx,zz,phi_kl(i,j,ispec),eta_kl(i,j,ispec)
          write(20,'(5e11.3)')xx,zz,cpI_kl(i,j,ispec),cpII_kl(i,j,ispec),cs_kl(i,j,ispec)
          write(21,'(5e11.3)')xx,zz,rhobb_kl(i,j,ispec),rhofbb_kl(i,j,ispec),ratio_kl(i,j,ispec)
          write(22,'(5e11.3)')xx,zz,phib_kl(i,j,ispec),eta_kl(i,j,ispec)
        enddo
      enddo
    enddo
    close(144)
    close(155)
    close(16)
    close(17)
    close(18)
    close(19)
    close(20)
    close(21)
    close(22)
  endif

endif

end subroutine save_adjoint_kernels


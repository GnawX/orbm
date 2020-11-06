!
! Copyright (C) 2001-2013 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
SUBROUTINE orbm
  !-----------------------------------------------------------------------
  !
  ! This routine calculates the orbital magnetization using the modern theory
  ! of orbital magnetization.
  !  m = \alpha/2 Im \sum_nk <d u_nk| \times (H+E-2\mu) |d u_nk>
  ! References:
  !   [1] Phys. Rev. B 63, 245101 (2001)       (norm-conserving GIPAW)
  !   [2] Phys. Rev. B 76, 024401 (2007)       (ultrasoft)
  !   [3] Phys. Rev. B 76, 165122 (2007)       (metallic systems)
  !   [4] Phys. Rev. Lett. 88, 086403 (2002)   (EPR g-tensor)
  ! Contributors:
  !   D. Ceresoli                        bare susceptibility and current
  !   A. P. Seitsonen and U. Gerstmann   GIPAW contributions
  !   E. Kucukbenli                      Ultrasoft and PAW
  !
  USE kinds,                  ONLY : dp
  USE io_global,              ONLY : stdout, ionode
  USE io_files,               ONLY : nwordwfc, iunwfc
  USE cell_base,              ONLY : omega, tpiba, tpiba2
  USE wavefunctions,   ONLY : evc
  USE klist,                  ONLY : nks, wk, xk, igk_k, ngk
  USE wvfct,                  ONLY : nbnd, npwx, wg, g2kin, current_k
  USE gvecw,                  ONLY : gcutw
  USE lsda_mod,               ONLY : current_spin, isk
  USE becmod,                 ONLY : calbec, becp, allocate_bec_type, deallocate_bec_type
  USE symme,                  ONLY : symmatrix
  USE constants,              ONLY : pi
  USE gvect,                  ONLY : ngm, g
  USE fft_base,               ONLY : dfftp, dffts
  USE uspp,                   ONLY : vkb, okvan
  USE lsda_mod,               ONLY : nspin
  USE gipaw_module,           ONLY : tens_fmt, q_gipaw, iverbosity, alpha, evq, &
                                     avogadro, g_e, gprime, filcurr, filfield, filnics, &
                                     nbnd_occ, a0_to_cm, isolve, &
                                     conv_threshold, job, restart_mode
  USE paw_gipaw,              ONLY : paw_vkb, paw_becp
  USE ions_base,              ONLY : nat
  USE buffers,                ONLY : get_buffer
  USE mp_pools,               ONLY : my_pool_id, me_pool, root_pool,  &
                                     inter_pool_comm, intra_pool_comm
  USE mp_images,              ONLY : my_image_id, inter_image_comm, nimage
  USE mp,                     ONLY : mp_sum, mp_barrier
#ifdef __BANDS
  USE gipaw_module,           ONLY : ibnd_start, ibnd_end
  USE mp_bands,               ONLY : intra_bgrp_comm, inter_bgrp_comm
#endif
  !-- local variables ----------------------------------------------------
  IMPLICIT NONE

  ! the following three quantities are for norm-conserving PPs
  complex(dp), allocatable, dimension(:,:,:) :: p_evc         ! p_k|evc>
  complex(dp), allocatable, dimension(:,:,:) :: vel_evc       ! v_{k,k}|evc>
  complex(dp), allocatable, dimension(:,:,:) :: evc1          ! du/dk
  ! temporary working array, same size as evc/evq
  complex(dp), allocatable :: aux(:,:)
  complex(dp), allocatable :: hpsi(:)
  real(dp) :: berry(3), mlc(3), mic(3)

  integer :: ik, iq, ik0, iq0
  integer :: ipol, jpol, i, ibnd, isign, ispin, ii, jj
  real(dp) :: tmp(3,3), q(3), k_plus_q(3), braket
  integer :: s_min, s_maj, s_weight
  complex(dp), external :: zdotc
  real(dp), external :: get_clock
  integer :: npw
  integer :: ind(2,3)
  
  ind(:,1) = (/ 2, 3 /)
  ind(:,2) = (/ 3, 1 /)
  ind(:,3) = (/ 1, 2 /)
  mlc = 0.d0
  mic = 0.d0
  berry=0.d0
 
  call start_clock('orbm')
  !-----------------------------------------------------------------------
  ! allocate memory
  !-----------------------------------------------------------------------
  allocate ( p_evc(npwx,nbnd,3), vel_evc(npwx,nbnd,3) )
  allocate ( aux(npwx,nbnd), evc1(npwx,nbnd,3), hpsi(npwx) )

  
  


  ! print memory estimate
  call gipaw_memory_report



  write(stdout, '(5X,''Computing the orbital magnetization'',$)')
  write(stdout, '(5X,''isolve='',I1,4X,''ethr='',E12.4)') isolve, conv_threshold

  ! check for recover file
  call check_for_restart_info ( ik0, iq0 )

  !====================================================================
  ! loop over k-points
  !====================================================================
  q(:) = 0.d0
  do ik = 1, nks
    ! skip if already done in a previous run
    if ( ik < ik0 ) cycle 

#ifdef __MPI
    if (me_pool == root_pool) &
    write(*, '(5X,''k-point #'',I5,'' of '',I5,6X,''pool #'',I3,4X,''cpu time:'',F10.1)') &
      ik, nks, my_pool_id+1, get_clock('GIPAW')
#else
    write(stdout, '(5X,''k-point #'',I5,'' of '',I5,4X,''cpu time:'',F10.1)') &
      ik, nks, get_clock('GIPAW')
#endif

    ! read wfcs from file and compute becp
    call get_buffer (evc, nwordwfc, iunwfc, ik)

    ! calculate du/dk    
    vel_evc(:,:,:) = (0.d0,0.d0)
    do ipol = 1,3
       call apply_vel(evc, vel_evc(1,1,ipol), ik, ipol, q)
       aux(:,:) = vel_evc(:,:,ipol)
       call greenfunction(ik, aux, evc1(1,1,ipol), q)
    enddo
    
    ! set up hamiltonian
    call allocate_bec_type(nkb, nbnd, becp)
    current_k = ik
    current_spin = isk(ik)
    npw = ngk(ik)

    call gk_sort(xk(1,ik), ngm, g, gcutw, npw, igk_k(1,ik), g2kin)
    g2kin(:) = g2kin(:) * tpiba2
    call init_us_2(npw, igk_k(1,ik), xk(1,ik), vkb)
    
#ifdef __BANDS
    call calbec_bands (npwx, npw, nkb, vkb, evc, becp%k, nbnd, ibnd_start, ibnd_end)
#else
    call calbec (npw, vkb, evc, becp, nbnd)
#endif
    
    ! loop over the bands
#ifdef __BANDS
    do ibnd = ibnd_start, ibnd_end
#else
    do ibnd = 1, nbnd
#endif
       do ipol = 1,3
          ii = ind(1, ipol)
          jj = ind(2, ipol)
          ! IC term
          braket = zdotc(npw, evc1(1,ibnd,ii), 1, evc1(1,ibnd,jj), 1)
          mic(ipol) = mic(ipol) + 2.d0*wg(ibnd,ik)*et(ibnd,ik)*imag(braket)
          berry(ipol) = berry(ipol) + 2.d0*wg(ibnd,ik)*et(ibnd,ik)*imag(braket)
          ! LC term
          call h_psi(npwx, npw, 1, evc1(1:npwx,ibnd,jj), hpsi)
          braket = zdotc(npw, evc1(1,ibnd,ii), 1, hpsi, 1)
          mlc(ipol) = mlc(ipol) + wg(ibnd,ik)*imag(braket)
          
          call h_psi(npwx, npw, 1, evc1(1:npwx,ibnd,ii), hpsi)
          braket = zdotc(npw, evc1(1,ibnd,jj), 1, hpsi, 1)
          mlc(ipol) = mlc(ipol) - wg(ibnd,ik)*imag(braket)
       enddo ! ipol
    enddo ! ibnd
  enddo ! ik
    
#ifdef __MPI
#ifdef __BANDS
  ! reduce over G-vectors
  call mp_sum( mlc, intra_bgrp_comm )
  call mp_sum( mic, intra_bgrp_comm )
  call mp_sum( berry, intra_bgrp_comm )
  ! reduce over band groups
  call mp_sum( mlc, inter_bgrp_comm )
  call mp_sum( mic, inter_bgrp_comm )
  call mp_sum( berry, inter_bgrp_comm )
#else
  ! reduce over G-vectors
  call mp_sum( mlc, intra_pool_comm )
  call mp_sum( mic, intra_pool_comm )
  call mp_sum( berry, intra_pool_comm )
#endif
  ! reduce over k-point pools
  call mp_sum( mlc, inter_pool_comm )
  call mp_sum( mic, inter_pool_comm )
  call mp_sum( berry, inter_pool_comm )
    
  morb = mlc + mic - 2.d0*ef*berry

  write(stdout,'(5X,''M_orb              = '',3(F14.6))') morb
    



  enddo  ! end of loop over k-point

  ! Parallel reductions
#ifdef __MPI
#ifdef __BANDS
  ! reduce over G-vectors
  call mp_sum( f_sum, intra_bgrp_comm )
  call mp_sum( f_sum_occ, intra_bgrp_comm )
  call mp_sum( f_sum_nelec, intra_bgrp_comm )
  call mp_sum( q_pGv, intra_bgrp_comm )
  call mp_sum( q_vGv, intra_bgrp_comm )
  call mp_sum( delta_g_rmc, intra_bgrp_comm)
  ! reduce over band groups
  call mp_sum( f_sum, inter_bgrp_comm )
  call mp_sum( f_sum_occ, inter_bgrp_comm )
  call mp_sum( f_sum_nelec, inter_bgrp_comm )
  call mp_sum( q_pGv, inter_bgrp_comm )
  call mp_sum( q_vGv, inter_bgrp_comm )
  call mp_sum( delta_g_rmc, inter_bgrp_comm) ! TODO: check this
#else
  ! reduce over G-vectors
  call mp_sum( f_sum, intra_pool_comm )
  call mp_sum( f_sum_occ, intra_pool_comm )
  call mp_sum( f_sum_nelec, intra_pool_comm )
  call mp_sum( q_pGv, intra_pool_comm )
  call mp_sum( q_vGv, intra_pool_comm )
  call mp_sum( delta_g_rmc, intra_pool_comm)
#endif

  ! reduce over k-point pools
  call mp_sum( f_sum, inter_pool_comm )
  call mp_sum( f_sum_occ, inter_pool_comm )
  call mp_sum( f_sum_nelec, inter_pool_comm )
  call mp_sum( q_pGv, inter_pool_comm )
  call mp_sum( q_vGv, inter_pool_comm )
  call mp_sum( delta_g_rmc, inter_pool_comm)

  call mp_sum( j_bare_s, inter_pool_comm )
  call mp_sum( sigma_diamagnetic, inter_pool_comm )
  call mp_sum( sigma_paramagnetic, inter_pool_comm )
  call mp_sum( sigma_paramagnetic_us, inter_pool_comm )
  call mp_sum( sigma_paramagnetic_aug, inter_pool_comm )
  call mp_sum( delta_g_rmc_gipaw, inter_pool_comm)
  call mp_sum( delta_g_so_dia, inter_pool_comm )
  call mp_sum( delta_g_so_para, inter_pool_comm )
  call mp_sum( delta_g_so_para_us, inter_pool_comm )
  call mp_sum( delta_g_so_para_aug, inter_pool_comm )

  ! reduce over images (q-star)
  if (nimage > 1) then
    call mp_sum( f_sum, inter_image_comm )
    call mp_sum( f_sum_occ, inter_image_comm )
    call mp_sum( f_sum_nelec, inter_image_comm )
    call mp_sum( q_pGv, inter_image_comm )
    call mp_sum( q_vGv, inter_image_comm )
    call mp_sum( j_bare_s, inter_image_comm )
    call mp_sum( sigma_diamagnetic, inter_image_comm )
    call mp_sum( sigma_paramagnetic, inter_image_comm )
    call mp_sum( sigma_paramagnetic_us, inter_image_comm )
    call mp_sum( sigma_paramagnetic_aug, inter_image_comm )
    call mp_sum( delta_g_rmc, inter_image_comm)
    call mp_sum( delta_g_rmc_gipaw, inter_image_comm)
    call mp_sum( delta_g_so_dia, inter_image_comm )
    call mp_sum( delta_g_so_para, inter_image_comm )
    call mp_sum( delta_g_so_para_us, inter_image_comm )
    call mp_sum( delta_g_so_para_aug, inter_image_comm )
  endif
#endif

  
  !====================================================================
  ! print out results
  !====================================================================
  write(stdout,'(5X,''End of magnetic susceptibility calculation'')')
  write(stdout,*)

  ! free memory as soon as possible
  deallocate( p_evc, vel_evc, aux, G_vel_evc )
  deallocate( svel_evc, u_svel_evc )
  
  ! f-sum rule
  call symmatrix (f_sum)
  call symmatrix (f_sum_occ)
  if (iverbosity > 0) then
    write(stdout, '(5X,''f-sum rule (1st term):'')')
    write(stdout, tens_fmt) f_sum
    write(stdout, '(5X,''f-sum rule (2nd term):'')')
    write(stdout, tens_fmt) f_sum_occ
  endif
  write(stdout, '(5X,''f-sum rule (should be '',F10.4''):'')') f_sum_nelec
  write(stdout, tens_fmt) f_sum + f_sum_occ
  if (job == 'f-sum') return

  ! F_{ij} = (2 - \delta_{ij}) Q_{ij}
  do ipol = 1, 3
    do jpol = 1, 3
      f_pGv(ipol,jpol,:) = 2.0_dp*q_pGv(ipol,jpol,:)
      if (ipol == jpol) f_pGv(ipol,jpol,:) = q_pGv(ipol,jpol,:)

      f_vGv(ipol,jpol,:) = 2.0_dp*q_vGv(ipol,jpol,:)
      if (ipol == jpol) f_vGv(ipol,jpol,:) = q_vGv(ipol,jpol,:)
    enddo
  enddo
  
  ! compute chi_bare both pGv and vGv terms
  chi_bare_pGv(:,:) = f_pGv(:,:,1) - 2.0_dp*f_pGv(:,:,0) + f_pGv(:,:,-1)
  chi_bare_pGv(:,:) = -0.5_dp * chi_bare_pGv(:,:) * alpha ** 2 &
       / ( q_gipaw * tpiba)**2
  if (iverbosity > 0) then
    write(stdout, '(5X,''chi_bare pGv (HH) in paratec units:'')')
    write(stdout, '(3(5X,3(F12.6,2X)/))') chi_bare_pGv(:,:) / alpha ** 2
  endif
  call symmatrix (chi_bare_pGv)
  if (iverbosity > 0) then
    write(stdout, '(3(5X,3(F12.6,2X)/))') chi_bare_pGv(:,:) / alpha ** 2
  endif
  
  chi_bare_vGv(:,:) = f_vGv(:,:,1) - 2.0_dp*f_vGv(:,:,0) + f_vGv(:,:,-1)
  chi_bare_vGv(:,:) = -0.5_dp * chi_bare_vGv(:,:) * alpha ** 2 &
       / ( q_gipaw * tpiba)**2
  if (iverbosity > 0) then
    write(stdout, '(5X,''chi_bare vGv (VV) in paratec units:'')')
    write(stdout, '(3(5X,3(F12.6,2X)/))') chi_bare_vGv(:,:) / alpha ** 2
  endif
  call symmatrix(chi_bare_vGv)
  if (iverbosity > 0) then
    write(stdout, '(3(5X,3(F12.6,2X)/))') chi_bare_vGv(:,:) / alpha ** 2
  endif

  ! convert from atomic units to 10^{-6} cm^3 / mol
  tmp(:,:) = chi_bare_pGv(:,:) * 1e6_dp * a0_to_cm**3.0_dp * avogadro
  write(stdout, '(5X,''chi_bare pGv (HH) in 10^{-6} cm^3/mol:'')')
  write(stdout, tens_fmt) tmp(:,:)

  tmp(:,:) = chi_bare_vGv(:,:) * 1e6_dp * a0_to_cm**3.0_dp * avogadro
  write(stdout, '(5X,''chi_bare vGv (VV) in 10^{-6} cm^3/mol:'')')
  write(stdout, tens_fmt) tmp(:,:)

  !--------------------------------------------------------------------
  ! now get the current and induced field
  !--------------------------------------------------------------------
  chi_bare_pGv(:,:) = chi_bare_pGv(:,:) / omega
  j_bare_s(:,:,:,:) = j_bare_s(:,:,:,:) * alpha / ( 2.d0 * q_gipaw * tpiba * omega )

  ! interpolate induced current
  allocate( j_bare(dfftp%nnr,3,3,nspin) )
  call interpolate_current(j_bare_s, j_bare)
  deallocate( j_bare_s )

  ! symmetrize the current
  do ispin = 1, nspin
#ifdef __MPI
    call psymmetrize_field(j_bare(:,:,:,ispin), 1)
#else
    call symmetrize_field(j_bare(:,:,:,ispin), 1)
#endif
  enddo

  ! compute induced field
  allocate( B_ind_r(dfftp%nnr,3,3,nspin), B_ind(ngm,3,3,nspin) )
  call biot_savart(j_bare, B_ind, B_ind_r)

  ! write fields to disk
  do i = 1, nspin
    if (trim(filcurr) /= '') call write_tensor_field(filcurr, i, j_bare(1,1,1,i))
    if (trim(filfield) /= '') call write_tensor_field(filfield, i, B_ind_r(1,1,1,i))
  enddo
  if (trim(filnics) /= '') call write_nics(filnics, B_ind_r(1,1,1,1))

  !--------------------------------------------------------------------
  ! calculate the chemical shifts or g-tensor
  !--------------------------------------------------------------------
  if (job == 'nmr') then
    ! compute bare chemical shift and print all results
    call compute_sigma_bare(B_ind, chi_bare_pGv, sigma_bare, sigma_shape)
    call print_chemical_shifts(sigma_shape, sigma_bare, sigma_diamagnetic, sigma_paramagnetic, &
                               sigma_paramagnetic_us, sigma_paramagnetic_aug, sigma_tot)
    if (ionode) then
       call output_magres_begin('nmr')
       call output_magres_nmr(chi_bare_pGv, chi_bare_vGv, sigma_tot)
       call output_magres_end
    endif
  elseif (job == 'g_tensor') then
    call get_rho_up_down
    call compute_delta_g_so(j_bare, s_maj, s_min, delta_g_so)
    call compute_delta_g_soo(j_bare, B_ind_r, s_maj, s_min, delta_g_soo, delta_g_soo2)
    call print_g_tensor(delta_g_rmc, delta_g_rmc_gipaw, delta_g_so, &
                        delta_g_soo, delta_g_soo2, delta_g_so_para, &
                        delta_g_so_para_us, delta_g_so_para_aug, &
                        delta_g_so_dia)
  endif

  ! deallocate fields
  deallocate( j_bare, B_ind_r, B_ind )
  
  call restart_cleanup ( )
  call stop_clock('suscept_crystal')
  
CONTAINS

  !====================================================================
  ! This routine evaluates the terms at q == 0
  !====================================================================
  SUBROUTINE suscept_crystal_inner_qzero
    IMPLICIT NONE

    ! check if my_image_id has to process this term
    if (mod(iq,nimage) /= my_image_id) return

    !if (job /= 'f-sum') then
       call compute_u_kq(ik, q)
    !else
    !   etq(:,ik) = et(:,ik)
    !endif

    if (ik == ik0 .and. iq < iq0 ) return
    call save_info_for_restart (ik, iq)

    ! 1. the diamagnetic contribution to the field: Eq.(58) of [1]/Eq.(9) of [4]
    if (job == 'nmr') then
      diamagnetic_corr_tensor = 0.0d0
      call diamagnetic_correction (diamagnetic_corr_tensor)
      sigma_diamagnetic = sigma_diamagnetic + diamagnetic_corr_tensor
    elseif (job == 'g_tensor') then
      diamagnetic_corr_tensor_so = 0.0d0
      call diamagnetic_correction_so (diamagnetic_corr_tensor_so)
      delta_g_so_dia = delta_g_so_dia + s_weight * diamagnetic_corr_tensor_so
    endif

    ! 2. the paramagnetic US augmentation: Eq.(30) of [2]
    if (okvan) then
      if (job == 'nmr') then
        paramagnetic_corr_tensor_aug = 0.d0
        call paramagnetic_correction_aug (paramagnetic_corr_tensor_aug, j_bare_s)
        sigma_paramagnetic_aug = sigma_paramagnetic_aug + paramagnetic_corr_tensor_aug
      elseif (job == 'g_tensor') then
        paramagnetic_corr_tensor_aug_so = 0.d0
        call paramagnetic_correction_aug_so (paramagnetic_corr_tensor_aug_so, j_bare_s)
        delta_g_so_para_aug = delta_g_so_para_aug + s_weight * paramagnetic_corr_tensor_aug_so
      endif
    endif

    ! 2. relativistic mass corrections for EPR
    if (job == 'g_tensor') call rmc(s_weight, delta_g_rmc, delta_g_rmc_gipaw)

    ! compute p_k|evc>, v_{k,k}|evc>, G_k v_{k,k}|evc> and s_{k,k}|evc>
    call apply_operators
    if (okvan) then 
        evq(:,:) = evc(:,:)
        call apply_occ_occ_us
    endif

    ! 3. f-sum rule
    do ipol = 1, 3 
      do jpol = 1, 3
#ifdef __BANDS
        do ibnd = ibnd_start, ibnd_end
#else
        do ibnd = 1, nbnd_occ(ik)
#endif
          ! count number of electrons
          if (ipol == 1 .and. jpol == 1) then
              braket = -real(zdotc(npw, evc(1,ibnd), 1, evc(1,ibnd), 1), dp)
              f_sum_nelec = f_sum_nelec + wg(ibnd,ik) * braket
          endif

          ! this is the "p-G-v" term
          braket = 2.d0*real(zdotc(npw, p_evc(1,ibnd,ipol), 1, G_vel_evc(1,ibnd,jpol), 1), dp)
!DEBUG if (ipol==1 .and. jpol==1) print*, 'AAA:', ibnd, braket*wg(ibnd,ik), ef-et(ibnd,ik)
          f_sum(ipol,jpol) = f_sum(ipol,jpol) + wk(ik) * braket

          ! this is the "occ-occ" term
          if (okvan) then
              braket = zdotc(npw, p_evc(1,ibnd,ipol), 1, u_svel_evc(1,ibnd,jpol), 1)
              f_sum_occ(ipol,jpol) = f_sum_occ(ipol,jpol) + wk(ik) * braket
         endif
        enddo
      enddo
    enddo

    ! 5. pGv and vGv contribution to chi_{bare}
    if (job /= 'f-sum') then
      do i = 1, 3
        call add_to_tensor(i,q_pGv(:,:,0), p_evc, G_vel_evc)
        call add_to_tensor(i,q_vGv(:,:,0), vel_evc, G_vel_evc)
      enddo
    endif
  END SUBROUTINE suscept_crystal_inner_qzero


  !====================================================================
  ! This routine evaluates the terms at q /= 0
  !====================================================================
  SUBROUTINE suscept_crystal_inner_q
    IMPLICIT NONE

    ! check if my_image_id has to process this term
    if (mod(iq,nimage) /= my_image_id) return
    if (job == 'f-sum') return
      
    if (ik == ik0 .and. iq < iq0) return
    call save_info_for_restart(ik, iq)
        
    ! compute the wfcs at k+q
    call compute_u_kq(ik, q)
        
    ! compute p_k|evc>, v_k|evc> and G_{k+q} v_{k+q,k}|evc>
    call apply_operators
      
    k_plus_q(1:3) = xk(1:3,ik) + q(1:3)
    call init_gipaw_2_no_phase(npw, igk_k(1,ik), k_plus_q, paw_vkb)

    ! pGv and vGv contribution to chi_bare
    call add_to_tensor(i,q_pGv(:,:,isign), p_evc, G_vel_evc)
    call add_to_tensor(i,q_vGv(:,:,isign), vel_evc, G_vel_evc)
        
    ! now the j_bare term 
    call add_to_current(j_bare_s(:,:,:,current_spin), evc, G_vel_evc)
    if (okvan) then
      call apply_occ_occ_us
      call add_to_current(j_bare_s(:,:,:,current_spin), evc, u_svel_evc)
    endif

    ! paramagnetic terms
    if (job == 'nmr') then
      call paramagnetic_correction(paramagnetic_corr_tensor, paramagnetic_corr_tensor_us, &
           G_vel_evc, u_svel_evc, i)
      call add_to_sigma_para(paramagnetic_corr_tensor, sigma_paramagnetic)
      if (okvan) call add_to_sigma_para(paramagnetic_corr_tensor_us, sigma_paramagnetic_us)
    elseif (job == 'g_tensor') then
      call paramagnetic_correction_so(paramagnetic_corr_tensor_so, paramagnetic_corr_tensor_us_so, &
           G_vel_evc, u_svel_evc, i)
      paramagnetic_corr_tensor_so = paramagnetic_corr_tensor_so * s_weight
      paramagnetic_corr_tensor_us_so = paramagnetic_corr_tensor_us_so * s_weight
      call add_to_sigma_para_so(paramagnetic_corr_tensor_so, delta_g_so_para)
      if (okvan) call add_to_sigma_para_so(paramagnetic_corr_tensor_us_so, delta_g_so_para_us)
    endif
  END SUBROUTINE suscept_crystal_inner_q
 

  !====================================================================
  ! compute p_k|evc>, v_{k+q,k}|evc>, G_{k+q} v_{k+q,k}|evc> and s_{k+q,k}|evc>
  !====================================================================
  SUBROUTINE apply_operators
    IMPLICIT NONE
    integer ipol

    p_evc(:,:,:) = (0.d0,0.d0)
    vel_evc(:,:,:) = (0.d0,0.d0)
    if(okvan) svel_evc(:,:,:) = (0.d0,0.d0)
    
    do ipol = 1, 3
      call apply_p(evc, p_evc(1,1,ipol), ik, ipol, q)
      call apply_vel(evc, vel_evc(1,1,ipol), ik, ipol, q)
      if (okvan) call apply_vel_NL('S', evc, svel_evc(1,1,ipol), ik, ipol, q)
      ! necessary because aux is overwritten by subroutine greenfunction
      aux(:,:) = vel_evc(:,:,ipol)
      call greenfunction(ik, aux, G_vel_evc(1,1,ipol), q)
    enddo
  END SUBROUTINE apply_operators


  !====================================================================
  ! compute |evq><evq|s_{k+q,k}|evc>
  !====================================================================
  SUBROUTINE apply_occ_occ_us
    IMPLICIT NONE
    integer ipol
    complex(dp), allocatable :: ps(:,:)

    allocate( ps(nbnd,nbnd) )
    do ipol = 1, 3
      ps = (0.d0,0.d0)
      aux(:,:) = svel_evc(:,:,ipol)
      CALL ZGEMM('C', 'N', nbnd_occ(ik), nbnd_occ(ik), npw, &
                (1.d0,0.d0), evq(1,1), npwx, aux(1,1), npwx, (0.d0,0.d0), &
                ps(1,1), nbnd)
#ifdef __MPI
#  ifdef __BANDS
      call mp_sum(ps, intra_bgrp_comm)
#  else
      call mp_sum(ps, intra_pool_comm)
#  endif
#endif
      aux = (0.d0,0.d0)
      CALL ZGEMM('N', 'N', npw, nbnd_occ(ik), nbnd_occ(ik), &
                (1.d0,0.d0), evq(1,1), npwx, ps(1,1), nbnd, (0.d0,0.d0), &
                aux(1,1), npwx)
      u_svel_evc(:,:,ipol) = -aux(:,:)
    enddo
    deallocate(ps)
  END SUBROUTINE apply_occ_occ_us


  !====================================================================
  ! add contribution the Q tensors
  ! Q_{\alpha,\beta} += <(e_i \times ul)_\alpha | (e_i \times ur)_\beta>
  !====================================================================
  SUBROUTINE add_to_tensor(i, qt, ul, ur)
    IMPLICIT NONE
    integer, intent(in) :: i
    real(dp), intent(inout) :: qt(3,3)
    complex(dp), intent(in) :: ul(npwx,nbnd,3), ur(npwx,nbnd,3)
    real(dp) :: braket
    integer :: ibnd, ia, ib, comp_ia, comp_ib, ind(3,3), mult(3,3)

    ! index for the cross product
    ind(:,1) = (/ 1, 3, 2/);  mult(:,1) = (/ 0,-1, 1 /)
    ind(:,2) = (/ 3, 2, 1/);  mult(:,2) = (/ 1, 0,-1 /)
    ind(:,3) = (/ 2, 1, 3/);  mult(:,3) = (/-1, 1, 0 /)

    do ia = 1, 3    ! ia = alpha
      comp_ia = ind(ia,i)
      if (mult(ia,i) == 0) cycle

      do ib = 1, 3    ! ib = beta
        comp_ib = ind(ib,i)
        if (mult(ib,i) == 0) cycle

#ifdef __BANDS
        do ibnd = ibnd_start, ibnd_end
#else
        do ibnd = 1, nbnd_occ(ik)
#endif
          braket = real(zdotc(npw, ul(1,ibnd,comp_ia), 1, &
                                   ur(1,ibnd,comp_ib), 1), dp)
          !WAS: qt(ia,ib) = qt(ia,ib) + wg(ibnd,ik) * braket * mult(ia,i) * mult(ib,i)
          qt(ia,ib) = qt(ia,ib) + wk(ik) * braket * mult(ia,i) * mult(ib,i)
        enddo  ! ibnd

      enddo  ! ib
    enddo  ! ia
  END SUBROUTINE add_to_tensor


  !====================================================================
  ! add contribution the the current
  ! j(r)_{\alpha,\beta} += <ul|J(r)|(B\times e_i \cdot ur)>
  !====================================================================
  SUBROUTINE add_to_current(j, ul, ur)
    IMPLICIT NONE
    real(dp), intent(inout) :: j(dffts%nnr,3,3)
    complex(dp), intent(in) :: ul(npwx,nbnd), ur(npwx,nbnd,3)
    real(dp) :: fact
    integer :: ibdir, icomp, ind(3,3), mult(3,3)

    ! index for the cross product
    ind(:,1) = (/ 1, 3, 2/);  mult(:,1) = (/ 0,-1, 1 /)
    ind(:,2) = (/ 3, 2, 1/);  mult(:,2) = (/ 1, 0,-1 /)
    ind(:,3) = (/ 2, 1, 3/);  mult(:,3) = (/-1, 1, 0 /)

    ! loop over B direction
    do ibdir = 1, 3
      if (i == ibdir) cycle
      icomp = ind(ibdir, i)
      fact = real(mult(ibdir,i)*isign)
      call j_para(fact, ul(1,1), ur(1,1,icomp), ik, q, j(1,1,ibdir))
    enddo
  END SUBROUTINE add_to_current


  !====================================================================
  ! Add contribution to current: Eq.(46) of [1]
  !====================================================================
  SUBROUTINE add_to_sigma_para(paramagnetic_correction, sigma_paramagnetic)
    IMPLICIT NONE
    real(dp), intent(in) :: paramagnetic_correction(3,3,nat)
    real(dp), intent(inout) :: sigma_paramagnetic(3,3,nat)
    real(dp) :: fact
    integer :: ibdir, icomp, ipol, ind(3,3), mult(3,3)
    
    ! index for the cross product
    ind(:,1) = (/ 1, 3, 2/);  mult(:,1) = (/ 0,-1, 1 /)
    ind(:,2) = (/ 3, 2, 1/);  mult(:,2) = (/ 1, 0,-1 /)
    ind(:,3) = (/ 2, 1, 3/);  mult(:,3) = (/-1, 1, 0 /)
    
    ! loop over B direction
    do ibdir = 1, 3
      if (i == ibdir) cycle
      icomp = ind(ibdir,i)
      fact = real(mult(ibdir,i)*isign)
      
      do ipol = 1, 3
         sigma_paramagnetic ( ipol, icomp, : ) &
              = sigma_paramagnetic ( ipol, icomp, : ) &
              + fact * paramagnetic_correction ( ipol, ibdir, : ) &
              / ( 2 * q_gipaw * tpiba )
      end do
    enddo
  END SUBROUTINE add_to_sigma_para


  SUBROUTINE add_to_sigma_para_so(paramagnetic_correction, sigma_paramagnetic)
    IMPLICIT NONE
    real(dp), intent(in) :: paramagnetic_correction(3,3)
    real(dp), intent(inout) :: sigma_paramagnetic(3,3)
    real(dp) :: fact
    integer :: ibdir, icomp, ipol, ind(3,3), mult(3,3)
    
    ! index for the cross product
    ind(:,1) = (/ 1, 3, 2/);  mult(:,1) = (/ 0,-1, 1 /)
    ind(:,2) = (/ 3, 2, 1/);  mult(:,2) = (/ 1, 0,-1 /)
    ind(:,3) = (/ 2, 1, 3/);  mult(:,3) = (/-1, 1, 0 /)
    
    ! loop over B direction
    do ibdir = 1, 3
      if (i == ibdir) cycle
      icomp = ind(ibdir,i)
      fact = real(mult(ibdir,i)*isign)
      
      do ipol = 1, 3
         sigma_paramagnetic ( ipol, icomp ) &
              = sigma_paramagnetic ( ipol, icomp ) &
              + fact * paramagnetic_correction ( ipol, ibdir ) &
              / ( 2 * q_gipaw * tpiba )
      end do
    enddo
  END SUBROUTINE add_to_sigma_para_so


  !====================================================================
  ! Restart and checkpointing section: read restart, if present
  !====================================================================
  SUBROUTINE check_for_restart_info(ik0_, iq0_)
    USE io_files, ONLY: seqopn
    IMPLICIT NONE

    integer, intent (out) :: ik0_, iq0_
    logical :: exst
    integer :: iunrec
    integer, external :: find_free_unit 
    character(80) :: job_

    ik0_ = 0
    iq0_ = -1

    if (restart_mode == 'from_scratch') then
       write(stdout, '(5X,''Starting from scratch'')')
       return
    endif

    iunrec = find_free_unit()
    call seqopn (iunrec, 'gipaw_recover', 'unformatted', exst)
    if ( exst ) then
       read(iunrec) job_
       if (job_ /= job) exst = .false.
    endif

    if (exst) then
       write(stdout, '(5X,''Restarting from a previous run'')')
       read (iunrec) ik0_, iq0_
       read (iunrec) f_sum(:,:), f_sum_occ(:,:), f_sum_nelec 
       read (iunrec) q_pgv(:,:,:), q_vgv(:,:,:)
       read (iunrec) j_bare_s(:,:,:,:)
       read (iunrec) sigma_bare, sigma_diamagnetic, sigma_paramagnetic, &
                     sigma_paramagnetic_us, sigma_paramagnetic_aug 
       read (iunrec) delta_g_rmc, delta_g_rmc_gipaw, delta_g_so, delta_g_soo, &
                     delta_g_soo2, delta_g_so_para, delta_g_so_para_aug, &
                     delta_g_so_para_us, delta_g_so_dia 
       write(stdout,'(5X,''Resuming from k-point #'',I5,3X,''and q #'',I3)') ik0_, iq0_
    endif
    close( unit = iunrec, status = 'keep' )
    call save_info_for_restart(ik0_,iq0_)
    return
  END SUBROUTINE check_for_restart_info


  !====================================================================
  ! Restart and checkpointing section: write restart
  !====================================================================
  SUBROUTINE save_info_for_restart(ik0_, iq0_)
    USE io_files,   ONLY : seqopn
    USE check_stop, ONLY : check_stop_now
    IMPLICIT NONE
    integer, intent (in) :: ik0_, iq0_
    integer :: iunrec
    integer, external :: find_free_unit 
    logical :: exst

    iunrec = find_free_unit()
    call seqopn(iunrec, 'gipaw_recover', 'unformatted', exst)

    write(iunrec) job
    write(iunrec) ik0_, iq0_
    write(iunrec) f_sum(:,:), f_sum_occ(:,:), f_sum_nelec 
    write(iunrec) q_pgv(:,:,:), q_vgv(:,:,:)
    write(iunrec) j_bare_s(:,:,:,:)
    write(iunrec) sigma_bare, sigma_diamagnetic, sigma_paramagnetic, &
                  sigma_paramagnetic_us, sigma_paramagnetic_aug 
    write(iunrec) delta_g_rmc, delta_g_rmc_gipaw, delta_g_so, delta_g_soo, &
                  delta_g_soo2, delta_g_so_para, delta_g_so_para_aug, &
                  delta_g_so_para_us, delta_g_so_dia 
    close( unit = iunrec, status = 'keep' )

    IF (check_stop_now()) THEN
       write(stdout,'(5X,''Stopping at k-point #'',I5,3X,''and q #'',I3)') ik0_, iq0_
       ! close files, print timings and stop the code
       call gipaw_closefil
       call print_clock_gipaw()
       call stop_code( .false. )
    end if
    return
  END SUBROUTINE save_info_for_restart


  !====================================================================
  ! Restart and checkpointing section: delete the restart file
  !====================================================================
  SUBROUTINE restart_cleanup ( )
    USE io_files, ONLY: seqopn
    IMPLICIT NONE
    integer :: iunrec
    integer, external :: find_free_unit 
    logical :: exst

    iunrec = find_free_unit()
    call seqopn(iunrec, 'gipaw_recover', 'unformatted', exst)
    close(unit = iunrec, status = 'delete')
    return
  END SUBROUTINE restart_cleanup


END SUBROUTINE orbm

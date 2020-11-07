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
  USE wavefunctions,          ONLY : evc
  USE noncollin_module,       ONLY : noncolin, npol
  USE klist,                  ONLY : nks, wk, xk, igk_k, ngk
  USE wvfct,                  ONLY : nbnd, npwx, wg, g2kin, current_k
  USE ener,                   ONLY : ef
  USE uspp,                   ONLY : nkb, vkb
  USE gvecw,                  ONLY : gcutw
  USE lsda_mod,               ONLY : lsda, current_spin, isk
  USE becmod,                 ONLY : calbec, becp, allocate_bec_type, deallocate_bec_type
  USE constants,              ONLY : pi
  USE gvect,                  ONLY : ngm, g
  USE fft_base,               ONLY : dfftp, dffts
  USE uspp,                   ONLY : vkb, okvan
  USE lsda_mod,               ONLY : nspin
  USE gipaw_module,           ONLY : tens_fmt, q_gipaw, iverbosity, alpha, evq, &
                                     avogadro, g_e, gprime, filcurr, filfield, filnics, &
                                     nbnd_occ, a0_to_cm, isolve, &
                                     conv_threshold, job, restart_mode
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
  complex(dp), allocatable, dimension(:,:,:) :: vel_evc       ! v_{k,k}|evc>
  complex(dp), allocatable, dimension(:,:,:) :: evc1          ! du/dk
  ! temporary working array, same size as evc/evq
  complex(dp), allocatable :: aux(:,:)
  complex(dp), allocatable :: hpsi(:)
  real(dp) :: berry(3), mlc(3), mic(3), morb

  integer :: ik
  integer :: i, ibnd,  ii, jj
  real(dp) :: q(3), braket
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
  allocate ( vel_evc(npwx*npol,nbnd,3), evc1(npwx*npol,nbnd,3) )
  allocate ( aux(npwx*npol,nbnd),  hpsi(npwx*npol) )
  call allocate_bec_type(nkb, nbnd, becp)

  ! print memory estimate
  call gipaw_memory_report

  write(stdout, '(5X,''Computing the orbital magnetization'',$)')
  write(stdout, '(5X,''isolve='',I1,4X,''ethr='',E12.4)') isolve, conv_threshold


  !====================================================================
  ! loop over k-points on the pool
  !====================================================================
  q(:) = 0.d0
  do ik = 1, nks

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
    
    ! initialize k, spin, g2kin used in h_psi    
    current_k = ik
    if (lsda) current_spin = isk(ik)
    npw = ngk(ik)
    call g2_kin( ik )
    if ( nkb >0 ) call init_us_2(npw, igk_k(1,ik), xk(1,ik), vkb)
    
    ! calculate du/dk    
    vel_evc(:,:,:) = (0.d0,0.d0)
    do i = 1,3
       call apply_vel(evc, vel_evc(1,1,i), ik, i, q)
       aux(:,:) = vel_evc(:,:,i)
       call greenfunction(ik, aux, evc1(1,1,i), q)
    enddo
    
    ! calculate orbital magnetization
    ! loop over the bands
#ifdef __BANDS
    do ibnd = ibnd_start, ibnd_end
#else
    do ibnd = 1, nbnd_occ(ik)
#endif
       do i = 1,3
          ii = ind(1, i)
          jj = ind(2, i)
          ! IC term
          braket = zdotc(npwx*npol, evc1(1,ibnd,ii), 1, evc1(1,ibnd,jj), 1)
          mic(i) = mic(i) + 2.d0*wg(ibnd,ik)*et(ibnd,ik)*imag(braket)
          berry(i) = berry(i) + 2.d0*wg(ibnd,ik)*et(ibnd,ik)*imag(braket)
          ! LC term
          call h_psi(npwx, npw, 1, evc1(:,ibnd,jj), hpsi)
          braket = zdotc(npwx*npol, evc1(1,ibnd,ii), 1, hpsi, 1)
          mlc(i) = mlc(i) + wg(ibnd,ik)*imag(braket)
          
          call h_psi(npwx, npw, 1, evc1(:,ibnd,ii), hpsi)
          braket = zdotc(npwx*npol, evc1(1,ibnd,jj), 1, hpsi, 1)
          mlc(i) = mlc(i) - wg(ibnd,ik)*imag(braket)
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
#endif
    
  morb = mlc + mic - 2.d0*ef*berry

  
  
  !====================================================================
  ! print out results
  !====================================================================
  write(stdout,'(5X,''End of orbital magnetization calculation'')')
  write(stdout,*)
  write(stdout,'(5X,''M_orb              = '',3(F14.6))') morb

  ! free memory as soon as possible
  deallocate( vel_evc, aux, evc1, hpsi )
  call deallocate_bec_type (becp)

  
  !call restart_cleanup ( )
  call stop_clock('orbm')

END SUBROUTINE orbm

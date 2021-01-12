! first and second derivative of Hamiltonian
!-----------------------------------------------------------------------
SUBROUTINE apply_vel2(psi, vel_psi, vel2_psi, ik)
  !-----------------------------------------------------------------------
  !
  ! ... Apply the velocity operator
  ! ...   v = p + dV^{NL}_{k,k}/dk = dH_k/dk
  ! ...
  ! ... Here we use Hartree atomic units, so that:
  ! ...   V^{NL} => V^{NL} * ryd_to_hartree
  !-----------------------------------------------------------------------
  USE kinds,                ONLY : DP
  USE klist,                ONLY : xk, igk_k, ngk
  USE wvfct,                ONLY : nbnd, npwx, npw, current_k, et
  USE lsda_mod,             ONLY : current_spin, isk, lsda
  USE becmod,               ONLY : becp, allocate_bec_type, deallocate_bec_type
  USE cell_base,            ONLY : tpiba
  USE uspp,                 ONLY : nkb, vkb
  USE orbm_module,          ONLY : q_orbm
  USE noncollin_module,     ONLY : npol

  !-- paramters ----------------------------------------------------------
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: ik         ! k-point
  COMPLEX(DP), INTENT(IN) :: psi(npwx*npol,nbnd)
  COMPLEX(DP), INTENT(OUT) :: vel_psi(npwx*npol,nbnd,3)
  COMPLEX(DP), INTENT(OUT) :: vel2_psi(npwx*npol,nbnd,3,3)

  !-- local variables ----------------------------------------------------
  real(dp), parameter :: ryd_to_hartree = 0.5d0
  complex(dp), allocatable :: vkb_save(:,:), vq0(:,:)
  complex(dp), allocatable :: vq(:,:,:,:)
  complex(dp), allocatable :: vq2(:,:,:,:)
  real(dp) :: dk, oldk(3)
  integer :: iq(3,3), iq2(3,3), i, ibnd

  call start_clock('apply_vel2')
  
  vel_psi = (0.d0,0.d0)
  vel2_psi = (0.d0,0.d0)
  
  iq(:,1) = (/1,0,0/)
  iq(:,2) = (/0,1,0/)
  iq(:,3) = (/0,0,1/)

  iq2(:,1) = (/0,1,1/)
  iq2(:,2) = (/1,0,1/)
  iq2(:,3) = (/1,1,0/)

  
  allocate(vq(npwx*npol,nbnd,3,2))
  allocate(vq2(npwx*npol,nbnd,3,2))
  allocate(vq0(npwx*npol,nbnd))
  
  do ibnd = 1, nbnd
     vq0(:,ibnd) = psi(:,ibnd)*et(ibnd,ik)
  enddo


  dk = q_orbm/2.d0
  
  ! allocate temporary arrays, save old NL-potential
  allocate(vkb_save(npwx,nkb))
  call allocate_bec_type(nkb, nbnd, becp)
  
  vkb_save = vkb
  oldk(:) = xk(:,ik)

  current_k = ik
  if (lsda) current_spin = isk(ik)
  npw = ngk(ik)
  
  
  do i = 1,3
     xk(:,ik) = oldk(:,ik) + iq(:,1)*dk
     call g2_kin(ik)
     call init_us_2(npw, igk_k(1,ik), xk(1,ik), vkb)
     call h_psi(npwx, npw, nbnd, psi, vq(:,:,i,1))
  enddo
  do i = 1,3
     xk(:,ik) = oldk(:,ik) - iq(:,1)*dk
     call g2_kin(ik)
     call init_us_2(npw, igk_k(1,ik), xk(1,ik), vkb)
     call h_psi(npwx, npw, nbnd, psi, vq(:,:,i,2))
  enddo
  do i = 1,3
     xk(:,ik) = oldk(:,ik) + iq2(:,1)*dk
     call g2_kin(ik)
     call init_us_2(npw, igk_k(1,ik), xk(1,ik), vkb)
     call h_psi(npwx, npw, nbnd, psi, vq2(:,:,i,1))
  enddo
  do i = 1,3
     xk(:,ik) = oldk(:,ik) - iq2(:,1)*dk
     call g2_kin(ik)
     call init_us_2(npw, igk_k(1,ik), xk(1,ik), vkb)
     call h_psi(npwx, npw, nbnd, psi, vq2(:,:,i,2))
  enddo
  
  call deallocate_bec_type (becp)
  
  do i = 1,3
     vel_psi(:,:,i) = (vq(:,:,i,1) + vq(:,:,i,2))*ryd_to_hartree/(2.d0*dk*tpiba)
  enddo

  do i = 1,3
     vel2_psi(:,:,i,i) = (vq(:,:,i,1) + vq(:,:,i,2) - 2.d0*vq0(:,:))*ryd_to_hartree/(dk*dk*tpiba*tpiba)
  enddo
  vel2_psi(:,:,2,3) = (vq2(:,:,1,1) + vq2(:,:,1,2) + 2.d0*vq0(:,:) - vq(:,:,2,1) - &
                       vq(:,:,2,2) - vq(:,:,3,1) - vq(:,:,3,2))*ryd_to_hartree/(2.d0*dk*dk*tpiba*tpiba)
  vel2_psi(:,:,3,1) = (vq2(:,:,2,1) + vq2(:,:,2,2) + 2.d0*vq0(:,:) - vq(:,:,3,1) - &
                       vq(:,:,3,2) - vq(:,:,1,1) - vq(:,:,1,2))*ryd_to_hartree/(2.d0*dk*dk*tpiba*tpiba)
  vel2_psi(:,:,1,2) = (vq2(:,:,3,1) + vq2(:,:,3,2) + 2.d0*vq0(:,:) - vq(:,:,1,1) - &
                       vq(:,:,1,2) - vq(:,:,2,1) - vq(:,:,2,2))*ryd_to_hartree/(2.d0*dk*dk*tpiba*tpiba)
  vel2_psi(:,:,3,2) = vel2_psi(:,:,2,3)
  vel2_psi(:,:,1,3) = vel2_psi(:,:,3,1)
  vel2_psi(:,:,2,1) = vel2_psi(:,:,1,2)


  ! restore NL-potential at k
  xk(:,ik) = oldk(:)
  vkb = vkb_save
  
  ! free memory
  deallocate(vkb_save, vq0, vq, vq2)

  call stop_clock('apply_vel2')

END SUBROUTINE apply_vel2

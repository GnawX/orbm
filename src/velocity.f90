! Copyright (C) 2001-2013 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

! The velocity operator is composed of three terms:
! v(\epsilon) = (1/i)[r, H - \epsilon S] = p + (1/i)[r, V_NL] - \epsilon (1/i)[r, S]
!
! the three terms are calculated in the following routines 
! apply_p      => apply p to the wave functions
! apply_vel_NL => apply (1/i)[r,V_NL] or (1/i)[r,S] to the wave functions
!
! Finally, the apply_vel subroutine is a driver that applies the velocity operator
!
!-----------------------------------------------------------------------
SUBROUTINE apply_p(psi, p_psi, ik, ipol)
  !-----------------------------------------------------------------------
  !
  ! ... Apply the kinetic part of the velocity operator
  ! ... |p_psi> = (G+k+q/2)_{ipol} |psi>
  !  
  USE kinds,                ONLY : DP
  USE klist,                ONLY : xk, igk_k, ngk
  USE wvfct,                ONLY : nbnd, npwx
  USE orbm_module,         ONLY : nbnd_occ
  USE gvect,                ONLY : g
  USE cell_base,            ONLY : tpiba
  USE noncollin_module,     ONLY : noncolin, npol

  !-- parameters ---------------------------------------------------------
  IMPLICIT none
  INTEGER, INTENT(IN) :: ik               ! k-point
  INTEGER, INTENT(IN) :: ipol             ! cartesian direction (1..3)
  COMPLEX(DP), INTENT(IN) :: psi(npwx*npol,nbnd)
  COMPLEX(DP), INTENT(OUT) :: p_psi(npwx*npol,nbnd)

  !-- local variables ----------------------------------------------------
  REAL(DP) :: gk(npwx)
  INTEGER :: ibnd
  INTEGER :: npw
  
  npw = ngk(ik)
  gk(1:npw) = (xk(ipol,ik)+g(ipol,igk_k(1:npw,ik)))*tpiba
  
  do ibnd = 1, nbnd_occ(ik)
    p_psi(1:npw,ibnd) = p_psi(1:npw,ibnd) + gk(1:npw)*psi(1:npw,ibnd)
    if (noncolin) then
      p_psi(npwx+1:npwx+npw,ibnd) = p_psi(npwx+1:npwx+npw,ibnd) + &
                              gk(1:npw)*psi(npwx+1:npwx+npw,ibnd)
    endif
  enddo

END SUBROUTINE apply_p

 
!-----------------------------------------------------------------------
SUBROUTINE apply_vel_NL(psi, vel_psi, ik, ipol)
  !-----------------------------------------------------------------------
  !
  ! ... Apply the non-local part of the velocity operator:
  ! ...   (1/i)[r,V_NL] => dV^{NL}_{k+q,k}/dk
  ! ... here we use Hartree atomic units, so that:
  ! ...   V^{NL} => V^{NL} * ryd_to_hartree
  ! ...
  !-----------------------------------------------------------------------
  USE kinds,                ONLY : DP
  USE klist,                ONLY : xk, igk_k, ngk
  USE wvfct,                ONLY : nbnd, npwx, current_k
  USE becmod,               ONLY : becp, calbec, allocate_bec_type, deallocate_bec_type
  USE uspp,                 ONLY : nkb, vkb
  USE cell_base,            ONLY : tpiba
  USE orbm_module,         ONLY : q_gipaw, nbnd_occ
  USE lsda_mod,             ONLY : current_spin, lsda, isk, nspin
  USE noncollin_module,     ONLY : noncolin, npol

  !-- paramters ----------------------------------------------------------
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: ipol       ! cartesian direction (1..3)
  INTEGER, INTENT(IN) :: ik         ! k-point
  COMPLEX(DP), INTENT(IN) :: psi(npwx*npol,nbnd)
  COMPLEX(DP), INTENT(INOUT) :: vel_psi(npwx*npol,nbnd)

  !-- local variables ----------------------------------------------------
  real(dp), parameter :: ryd_to_hartree = 0.5d0
  complex(dp), allocatable :: aux(:,:), vkb_save(:,:)
  real(dp) :: dk, dxk(3)
  integer :: isign
  integer :: npw

  ! if no projectors, return
  if (nkb == 0) return

  ! set dk
  dk = q_gipaw/2.d0


  ! initialization
  npw = ngk(ik)
  current_k = ik
  if (lsda) current_spin = isk(ik)

  ! allocate temporary arrays, save old NL-potential
  allocate(aux(npwx*npol,nbnd), vkb_save(npwx,nkb))
  vkb_save = vkb


  !====================================================================
  ! compute (1/2|dk|) ( V^{NL}_{k+dk,k+dk} |psi> - 
  !                     V^{NL}_{k-dk,k-dk} |psi> )
  ! or the same, with S.
  !====================================================================
  do isign = -1,1,2
      dxk(:) = xk(:,ik)
      dxk(ipol) = dxk(ipol) + isign * dk     ! k \pm dk

      ! compute <\beta(k \pm dk)| and project on |psi>
      call init_us_2_no_phase(npw, igk_k(1,ik), dxk, vkb)
      call allocate_bec_type(nkb, nbnd_occ(ik), becp)

      call calbec (npw, vkb, psi, becp, nbnd_occ(ik))

      aux = (0.d0,0.d0)
      
      ! apply |\beta(k \pm dk+q)>D<\beta(k \pm dk)| to |psi>

      call add_vuspsi(npwx, npw, nbnd_occ(ik), aux)
      call deallocate_bec_type(becp)

      vel_psi = vel_psi + dble(isign) * ryd_to_hartree * aux/(2.d0*dk*tpiba)

  enddo


  ! restore NL-potential at k
  vkb = vkb_save
  
  ! free memory
  deallocate(aux, vkb_save)

END SUBROUTINE apply_vel_NL


!-----------------------------------------------------------------------
SUBROUTINE apply_vel(psi, vel_psi, ik, ipol)
  !-----------------------------------------------------------------------
  !
  ! ... Apply the velocity operator
  !-----------------------------------------------------------------------
  USE kinds,                ONLY : DP
  USE wvfct,                ONLY : nbnd, npwx 
  USE noncollin_module,     ONLY : npol


  !-- paramters ----------------------------------------------------------
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: ipol       ! cartesian direction (1..3)
  INTEGER, INTENT(IN) :: ik         ! k-point
  COMPLEX(DP), INTENT(IN) :: psi(npwx*npol,nbnd)
  COMPLEX(DP), INTENT(OUT) :: vel_psi(npwx*npol,nbnd)


  call start_clock('apply_vel')

  call apply_vel_NL(psi, vel_psi, ik, ipol)

  call apply_p(psi, vel_psi, ik, ipol)

  call stop_clock('apply_vel')

END SUBROUTINE apply_vel


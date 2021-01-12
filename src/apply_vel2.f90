! first and second derivative of Hamiltonian
!-----------------------------------------------------------------------
SUBROUTINE apply_vel(psi, vel_psi, ik, ipol)
  !-----------------------------------------------------------------------
  !
  ! ... Apply the velocity operator
  ! ...   v = p + dV^{NL}_{k,k}/dk = dH_k/dk
  ! ...
  ! ... Here we use Hartree atomic units, so that:
  ! ...   V^{NL} => V^{NL} * ryd_to_hartree
  !-----------------------------------------------------------------------
  USE kinds,                ONLY : DP
  USE ions_base,            ONLY : nat, ityp, ntyp => nsp
  USE klist,                ONLY : xk
  USE wvfct,                ONLY : nbnd, npwx, npw, igk, current_k
  USE lsda_mod,             ONLY : current_spin, isk
  USE becmod,               ONLY : becp
  USE cell_base,            ONLY : tpiba
  USE g_tensor_module,      ONLY : init_us_2_no_phase
  USE g_tensor_module,      ONLY : q_gipaw
  USE pwcom

  !-- paramters ----------------------------------------------------------
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: ipol       ! cartesian direction (1..3)
  INTEGER, INTENT(IN) :: ik         ! k-point
  COMPLEX(DP), INTENT(IN) :: psi(npwx,nbnd)
  COMPLEX(DP), INTENT(OUT) :: vel_psi(npwx,nbnd)

  !-- local variables ----------------------------------------------------
  real(dp), parameter :: ryd_to_hartree = 0.5d0
  !real(dp), parameter :: q_gipaw = 0.02d0
  complex(dp), allocatable :: aux(:,:), vkb_save(:,:)
  real(dp) :: dk, oldk(3)
  integer :: isign

  call start_clock('apply_vel')
  vel_psi = (0.d0,0.d0)

  ! set dk (= delta_k ?)
  !dk = q_gipaw/2.d0/tpiba
  dk = q_gipaw/tpiba
  
  ! allocate temporary arrays, save old NL-potential
  allocate(aux(npwx,nbnd), vkb_save(npwx,nkb))
  vkb_save = vkb
  oldk(:) = xk(:,ik)

  current_k = ik
  if (lsda) current_spin = isk(ik)
  call gk_sort(xk(1,ik), ngm, g, ecutwfc/tpiba2, npw, igk, g2kin)

  !====================================================================
  ! compute (1/2|dk|) ( H_{k+dk} |psi> - H_{k-dk}|psi> )
  !====================================================================
  allocate(becp(nkb,nbnd))
  do isign = -1,1,2
    xk(ipol,ik) = oldk(ipol) + isign * dk     ! k + dk

    ! compute <\beta(k \pm dk)| and project on |psi>
    call g2_kin(ik)
    call init_us_2_no_phase(npw, igk, xk(1,ik), vkb)
    aux = (0.d0,0.d0)
    call h_psi(npwx, npw, nbnd, psi, aux)
    vel_psi = vel_psi + dble(isign) * ryd_to_hartree * aux/(2.d0*dk*tpiba)
  enddo
  deallocate(becp)

  ! restore NL-potential at k
  xk(:,ik) = oldk(:)
  vkb = vkb_save
  
  ! free memory
  deallocate(aux, vkb_save)

  call stop_clock('apply_vel')

END SUBROUTINE apply_vel

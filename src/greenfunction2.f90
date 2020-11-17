! Copyright (C) 2001-2005 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
SUBROUTINE greenfunction2(ik, psi, g_psi)
  !-----------------------------------------------------------------------
  !
  ! ... Apply the Green function operator
  ! ...
  ! ... [H^(0) - alpha P_i - E_i^(0)] |Gpsi_i> = - Q_i H^(1)|psi_i^(0)>
  ! ...
  ! ... We use Hartree atomic units; since G is the inverse of an
  ! ... energy: G => G / ryd_to_hartree
  !
  USE kinds,                       ONLY : DP
  USE io_global,                   ONLY : stdout  
  USE becmod,                      ONLY : bec_type, becp, calbec, &
                                          allocate_bec_type, deallocate_bec_type
  USE wavefunctions,               ONLY : evc
  USE pwcom,                       ONLY : ef
  USE wvfct,                       ONLY : nbnd, et, npwx, g2kin
  USE gvect,                       ONLY : g
  USE uspp,                        ONLY : nkb, vkb
  USE mp_pools,                    ONLY : intra_pool_comm
  USE mp,                          ONLY : mp_sum
  USE ldaU,                        ONLY : lda_plus_u, wfcU
  USE io_files,                    ONLY : iunhub, nwordwfcU
  USE buffers,                     ONLY : get_buffer
  USE cell_base,                   ONLY : tpiba
  USE klist,                       ONLY : lgauss, xk, degauss, ngauss, igk_k, ngk
  USE noncollin_module,            ONLY : noncolin, npol
  USE orbm_module

  !-- parameters ---------------------------------------------------------
  IMPLICIT none
  INTEGER, INTENT(IN) :: ik
  COMPLEX(DP), INTENT(INOUT) :: psi(npwx*npol,nbnd)  ! psi is H1*psi and is changed on output!!!
  COMPLEX(DP), INTENT(OUT) :: g_psi(npwx*npol,nbnd)


  !-- local variables ----------------------------------------------------
  real(dp), parameter :: ryd_to_hartree = 0.5d0
  complex(dp), allocatable :: ps(:,:), work (:)
  real(dp), allocatable :: h_diag (:,:), eprec (:)
  real(dp) :: anorm, thresh, gk(3)
  integer :: ibnd, jbnd, ig, lter
  logical :: conv_root
  complex(dp), external :: zdotc
  real(dp), external :: wgauss, w0gauss
  real(dp) :: wg1, w0g, wgp, wwg, deltae, theta
  external ch_psi_all, cg_psi
  integer :: npw
 
  ! start clock
  call start_clock ('greenf')
  npw = ngk(ik)

  ! allocate memory
  allocate (work(npwx*npol), ps(nbnd,nbnd), h_diag(npwx*npol,nbnd), eprec(nbnd))
  call allocate_bec_type(nkb, nbnd, becp)

  !====================================================================
  ! apply -Q to the r.h.s.
  !====================================================================
  ! project on <evq|: ps(i,j) = <evq(i)|psi(j)>
  ps = (0.d0,0.d0)
  
  do ibnd = 1, nbnd
     ps(ibnd) = zdotc(npwx*npol, evc(1,ibnd), 1, psi(1,ibnd), 1)
  enddo
  
#ifdef __MPI
  call mp_sum(ps, intra_pool_comm)
#endif

  do ibnd = 1, nbnd
     psi(:,ibnd) = evc(:,ibnd)*ps(ibnd) - psi(:,ibnd)
  enddo
  

  !====================================================================
  ! solve the linear system (apply G_{k+q})
  !====================================================================
  ! convergence treshold
  thresh = sqrt(conv_threshold)   ! sqrt(of that of PARATEC)

  ! use the hamiltonian at k
  !do ig = 1, npw
  !  gk(1) = (xk(1,ik) + g(1,igk_k(ig,ik)) ) * tpiba
  !  gk(2) = (xk(2,ik) + g(2,igk_k(ig,ik)) ) * tpiba
  !  gk(3) = (xk(3,ik) + g(3,igk_k(ig,ik)) ) * tpiba
  !  g2kin (ig) = gk(1)**2 + gk(2)**2 + gk(3)**2
  !enddo

  ! preconditioning of the linear system
  work = (0.d0,0.d0)

  do ibnd = 1, nbnd
     work(1:npw) = g2kin(1:npw)*evc(1:npw, ibnd)
     if (noncolin) then
         work(npwx+1:npwx+npw) = g2kin(1:npw)*evc(npwx+1:npwx+npw,ibnd)
     endif
     eprec (ibnd) = 1.35d0 * zdotc (npwx*npol, evc (1, ibnd), 1, work, 1)
  enddo

#ifdef __MPI
  call mp_sum ( eprec, intra_pool_comm )
#endif

  h_diag = 0.d0

  do ibnd = 1, nbnd
     do ig = 1, npw
        h_diag (ig, ibnd) = 1.d0 / max (1.0d0, g2kin (ig) / eprec (ibnd) )
        if (noncolin) then
           h_diag (ig+npwx, ibnd) = h_diag(ig, ibnd)
        endif
     enddo
  enddo


  !call calbec (npw, vkb, psi, becp, nbnd)


  ! initial guess
  g_psi(:,:) = (0.d0, 0.d0)

  ! solve linear system  
  conv_root = .true.
  call cgsolve_all (ch_psi_all2, cg_psi, et(1,ik), psi, g_psi, &
       h_diag, npwx, npw, thresh, ik, lter, conv_root, anorm, &
       nbnd_occ(ik), npol )


  if (iverbosity > 20) &
    write(stdout, '(5X,''cgsolve_all iterations '',I3,4X,''anorm='',E12.2)')  lter, anorm

  ! convert to Hartree
  g_psi(:,:) = g_psi(:,:) / ryd_to_hartree

  FLUSH(stdout) 
  call stop_clock('greenf')
 
  ! free memory
  deallocate (work, h_diag, eprec, ps)
  call deallocate_bec_type (becp)

END SUBROUTINE greenfunction2


subroutine ch_psi_all2 (n, h, ah, e, ik, m)
  !-----------------------------------------------------------------------
  !
  ! This routine applies the operator ( H - \epsilon S + alpha_pv P_v)
  ! to a vector h. The result is given in Ah.
  !
  USE kinds,        ONLY : dp
  USE wvfct,        ONLY : npwx, nbnd
  USE uspp,         ONLY : vkb
!  USE becmod,       ONLY : becp, calbec
  USE orbm_module, ONLY : nbnd_occ, alpha_pv
  USE mp_pools,     ONLY : intra_pool_comm
  USE mp,           ONLY : mp_sum
  USE noncollin_module,      ONLY : noncolin, npol
  USE wavefunctions,         ONLY : evc
  
  implicit none

  integer :: n, m, ik
  ! input: the dimension of h
  ! input: the number of bands
  ! input: the k point
  real(DP) :: e (m)
  ! input: the eigenvalue

  complex(DP) :: h (npwx*npol, m), ah (npwx*npol, m)
  ! input: the vector
  ! output: the operator applied to the vector
  !
  !   local variables
  !
  integer :: ibnd, ikq, ig
  ! counter on bands
  ! the point k+q
  ! counter on G vetors

  complex(DP), allocatable :: ps (:), hpsi (:,:), spsi (:,:)
  ! scalar products
  ! the product of the Hamiltonian and h
  ! the product of the S matrix and h

  call start_clock ('ch_psi')
  allocate (ps  ( m ))    
  allocate (hpsi( npwx*npol , m))    
  allocate (spsi( npwx*npol , m))    
  hpsi (:,:) = (0.d0, 0.d0)
  spsi (:,:) = (0.d0, 0.d0)
  !
  !   compute the product of the hamiltonian with the h vector
  !

  call h_psi ( npwx, n, m, h, hpsi )
  CALL s_psi ( npwx, n, m, h, spsi )
  
  call start_clock ('last')
  !
  !   then we compute the operator H-epsilon S
  !
  ah = (0.d0,0.d0)

  do ibnd = 1, m
     do ig = 1, n
        ah (ig, ibnd) = hpsi (ig, ibnd) - e (ibnd) * spsi (ig, ibnd)
     enddo
  enddo
  IF (noncolin) THEN
     DO ibnd = 1, m
        DO ig = 1, n
           ah(ig+npwx, ibnd) = hpsi(ig+npwx, ibnd) - e(ibnd) * spsi(ig+npwx, ibnd)
        ENDDO
     ENDDO
  ENDIF
  !
  !   Here we compute the projector in the valence band
  !
  ikq = ik
  ps (:) = (0.d0, 0.d0)
  
  do ibnd = 1, m
     ps(ibnd) = zdotc(npwx*npol, evc(1,ibnd), 1, spsi(1,ibnd), 1)
  enddo
  
  ps(:) = ps(:)*alpha_pv
  
#ifdef __MPI
  call mp_sum(ps, intra_pool_comm)
#endif

  do ibnd = 1, m
     hpsi(:,ibnd) = evc(:,ibnd)*ps(ibnd) + hpsi(:,ibnd)
  enddo
  
  spsi(:,:) = hpsi(:,:)

  !
  !    And apply S again
  !

  !call calbec(n, vkb, hpsi, becp, m)
  !call s_psi (npwx, n, m, hpsi, spsi)


  do ibnd = 1, m
     do ig = 1, n
        ah (ig, ibnd) = ah (ig, ibnd) + spsi (ig, ibnd)
     enddo
  enddo
  IF (noncolin) THEN
       DO ibnd = 1, m
          DO ig = 1, n
             ah (ig+npwx, ibnd) = ah (ig+npwx, ibnd) + spsi (ig+npwx, ibnd)
          ENDDO
       ENDDO
   END IF

  deallocate (spsi)
  deallocate (hpsi)
  deallocate (ps)
  call stop_clock ('last')
  call stop_clock ('ch_psi')
  return
end subroutine ch_psi_all2


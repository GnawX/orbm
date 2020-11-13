! Copyright (C) 2001-2005 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
SUBROUTINE greenfunction(ik, psi, g_psi)
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
  USE optic_module

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
  ! apply -Q_{k+q} to the r.h.s.
  !====================================================================
  ! project on <evq|: ps(i,j) = <evq(i)|psi(j)>
  ps = (0.d0,0.d0)

  if (lgauss) then
     ! metallic case
     if (noncolin) then
         CALL zgemm('C', 'N', nbnd, nbnd_occ (ik), npwx*npol, (1.d0,0.d0), evc(1,1), &
                     npwx*npol, psi(1,1), npwx*npol, (0.d0,0.d0), ps(1,1), nbnd)
     else
         CALL zgemm('C', 'N', nbnd, nbnd_occ (ik), npw, (1.d0,0.d0), evc(1,1), &
                     npwx, psi(1,1), npwx, (0.d0,0.d0), ps(1,1), nbnd)
     endif
  
     do ibnd = 1, nbnd_occ(ik)
        wg1 = wgauss ((ef-et(ibnd,ik)) / degauss, ngauss)
        w0g = w0gauss((ef-et(ibnd,ik)) / degauss, ngauss) / degauss
        do jbnd = 1, nbnd
           wgp = wgauss ( (ef - et(jbnd,ik)) / degauss, ngauss)
           deltae = et(jbnd,ik) - et(ibnd,ik)
           theta = wgauss (deltae / degauss, 0)
           wwg = wg1 * (1.d0 - theta) + wgp * theta
           if (jbnd <= nbnd_occ(ik)) then
              if (abs (deltae) > 1d-5) then
                 wwg = wwg + alpha_pv * theta * (wgp - wg1) / deltae
              else
                 ! if the two energies are too close takes the limit of the 0/0 ratio
                 wwg = wwg - alpha_pv * theta * w0g
              endif
           endif
           ps(jbnd,ibnd) = wwg * ps(jbnd,ibnd)
        enddo
        call dscal(2*npwx*npol, wg1, psi(1,ibnd), 1)
     enddo

  else
     ! insulators
     if (noncolin) then
     
         CALL zgemm('C', 'N', nbnd_occ (ik), nbnd_occ (ik), npwx*npol, (1.d0,0.d0), evc(1,1), &
                    npwx*npol, psi(1,1), npwx*npol, (0.d0,0.d0), ps(1,1), nbnd)
     else
         CALL zgemm('C', 'N', nbnd_occ (ik), nbnd_occ (ik), npw, (1.d0,0.d0), evc(1,1), &
                    npwx, psi(1,1), npwx, (0.d0,0.d0), ps(1,1), nbnd)
     endif
 
  endif

#ifdef __MPI
  call mp_sum(ps, intra_pool_comm)
#endif

  if (lgauss) then
     ! metallic case
     if (noncolin) then
          CALL zgemm( 'N', 'N', npwx*npol, nbnd_occ(ik), nbnd, (1.d0,0.d0), &
          evc(1,1), npwx*npol, ps(1,1), nbnd, (-1.d0,0.d0), psi(1,1), npwx*npol )
     else
          CALL zgemm( 'N', 'N', npw, nbnd_occ(ik), nbnd, (1.d0,0.d0), &
          evc(1,1), npwx, ps(1,1), nbnd, (-1.d0,0.d0), psi(1,1), npwx )
     endif

  else
     ! insulators
     if (noncolin) then
         CALL zgemm( 'N', 'N', npwx*npol, nbnd_occ(ik), nbnd_occ(ik), (1.d0,0.d0), &
         evc(1,1), npwx*npol, ps(1,1), nbnd, (-1.d0,0.d0), psi(1,1), npwx*npol )
     else
         CALL zgemm( 'N', 'N', npw, nbnd_occ(ik), nbnd_occ(ik), (1.d0,0.d0), &
         evc(1,1), npwx, ps(1,1), nbnd, (-1.d0,0.d0), psi(1,1), npwx )
     endif

  endif

  !! this is the old code for norm-conserving:
  !! |psi> = -(1 - |evq><evq|) |psi>
  !!CALL zgemm('N', 'N', npw, nbnd_occ(ik), nbnd_occ(ik), &
  !!           (1.d0,0.d0), evq(1,1), npwx, ps(1,1), nbnd, (-1.d0,0.d0), &
  !!           psi(1,1), npwx)


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
  call cgsolve_all (ch_psi_all, cg_psi, et(1,ik), psi, g_psi, &
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

END SUBROUTINE greenfunction

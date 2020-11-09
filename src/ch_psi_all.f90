!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

!-----------------------------------------------------------------------
subroutine ch_psi_all (n, h, ah, e, ik, m)
  !-----------------------------------------------------------------------
  !
  ! This routine applies the operator ( H - \epsilon S + alpha_pv P_v)
  ! to a vector h. The result is given in Ah.
  !
  USE kinds,        ONLY : dp
  USE wvfct,        ONLY : npwx, nbnd
  USE uspp,         ONLY : vkb
  USE becmod,       ONLY : becp, calbec
  USE gipaw_module, ONLY : nbnd_occ, alpha_pv
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

  complex(DP), allocatable :: ps (:,:), hpsi (:,:), spsi (:,:)
  ! scalar products
  ! the product of the Hamiltonian and h
  ! the product of the S matrix and h

  call start_clock ('ch_psi')
  allocate (ps  ( nbnd , m))    
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
  ps (:,:) = (0.d0, 0.d0)
  if (noncolin) then
     CALL zgemm ('C', 'N', nbnd_occ (ikq) , m, npwx*npol, (1.d0, 0.d0) , evc, &
            npwx*npol, spsi, npwx*npol, (0.d0, 0.d0) , ps, nbnd)
  else
     call zgemm ('C', 'N', nbnd_occ (ikq) , m, n, (1.d0, 0.d0) , evc, &
            npwx, spsi, npwx, (0.d0, 0.d0) , ps, nbnd)
  endif
  ps (:,:) = ps(:,:) * alpha_pv

#ifdef __MPI
  call mp_sum ( ps, intra_pool_comm )
#endif

  hpsi (:,:) = (0.d0, 0.d0)
  if (noncolin) then
      CALL zgemm ('N', 'N', npwx*npol, m, nbnd_occ (ikq) , (1.d0, 0.d0) , evc, &
            npwx*npol, ps, nbnd, (1.d0, 0.d0) , hpsi, npwx*npol)
  else
      call zgemm ('N', 'N', n, m, nbnd_occ (ikq) , (1.d0, 0.d0) , evc, &
            npwx, ps, nbnd, (1.d0, 0.d0) , hpsi, npwx)
  endif
  spsi(:,:) = hpsi(:,:)

  !
  !    And apply S again
  !

  call calbec(n, vkb, hpsi, becp, m)
  call s_psi (npwx, n, m, hpsi, spsi)


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
end subroutine ch_psi_all

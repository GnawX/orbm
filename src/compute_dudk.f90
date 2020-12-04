!-----------------------------------------------------------------------
! Compute the k derivative of evc -> evc1
!-----------------------------------------------------------------------
SUBROUTINE compute_dudk(ik)
  
  USE kinds,                  ONLY : dp
  USE cell_base,              ONLY : tpiba2
  USE uspp,                   ONLY : nkb, vkb
  USE gvect,                  ONLY : ngm, g
  USE gvecw,                  ONLY : gcutw
  USE klist,                  ONLY : nks, nkstot, ngk, xk, igk_k 
  USE wvfct,                  ONLY : nbnd, npwx, et, current_k, g2kin
  USE lsda_mod,               ONLY : nspin, lsda, isk, current_spin
  USE orbm_module,            ONLY : evc1, vel_evc, dudk_method
  USE noncollin_module,       ONLY : npol
  USE io_global,              ONLY : stdout
  
  implicit none

  COMPLEX(dp), allocatable :: aux(:,:)
  integer :: ik, ipol, npw


  write(stdout,*)
  write(stdout,'(5X,''Computing du/dk '',$)')


  if (trim(dudk_method) == 'covariant') then
  
     write(stdout,'(''(covariant derivative)'')') 
     call dudk_covariant(ik)
    
  elseif (trim(dudk_method) == 'kdotp') then
  
     write(stdout,'(''(k \dot p perturbation)'')') 
     allocate ( aux(npwx*npol, nbnd) )
     npw = ngk(ik)
     ! initialize k, spin, g2kin used in h_psi    
     current_k = ik
     if (lsda) current_spin = isk(ik)
     call gk_sort(xk(1,ik), ngm, g, gcutw, npw, igk_k(1,ik), g2kin)
     g2kin(:) = g2kin(:) * tpiba2
     if (nkb > 0) call init_us_2(npw, igk_k(1,ik), xk(1,ik), vkb)
    
     do ipol = 1,3
        aux(:,:) = vel_evc(:,:,ipol)
        call greenfunction(ik, aux, evc1(1,1,ipol))
     enddo
       

  elseif (trim(dudk_method) == 'sos') then
  
     write(stdout,'(''(sum over states)'')')
     call dudk_sos(ik)

  else
    write(stdout,*)
    call errore('compute_dudk', 'unknown du/dk method: '//trim(dudk_method), 1)
  endif


END SUBROUTINE compute_dudk
  
  
  
  SUBROUTINE dudk_covariant(ik)
  USE kinds,                ONLY : dp  
  USE cell_base,            ONLY : tpiba
  USE wvfct,                ONLY : nbnd, npwx
  USE wavefunctions,        ONLY : evc
  USE klist,                ONLY : xk, ngk
  USE noncollin_module,     ONLY : npol, noncolin
  USE mp_pools,             ONLY : intra_pool_comm
  USE mp,                   ONLY : mp_sum
  !USE io_global,            ONLY : stdout
  !USE buffers,              ONLY : get_buffer, save_buffer
  USE matrix_inversion
  USE orbm_module         
  IMPLICIT NONE
  complex(dp), external :: zdotc 
  complex(dp), allocatable :: overlap(:,:), evc0(:,:)
  real(dp) :: q(3), delta_k
  integer :: ik, i, sig, ibnd, jbnd
  integer :: ipol, npw
  
  npw = ngk(ik)
  ! allocate overlap matrix 
  allocate(overlap(nbnd_occ(ik), nbnd_occ(ik)))
  !allocate(evc0(npwx*npol,nbnd))

  delta_k = q_orbm/2.d0
  !q=0.d0
  
  !call compute_u_kq(ik,q)
  !evc0 = evq

  ! loop over crystal directions
  do ipol = 1, 3
    evc1(:,:,ipol) = (0.d0,0.d0)
    
    ! loop over +/-1
    do sig = -1, 1, 2

      ! set the k-point
      q(:) = 0.d0
      q(ipol) = delta_k * sig
      !!write(*,'(5X,''overlap: ik='',I4,'' ipol='',I1,'' sig='',I2)') &
      !!      ik, ipol, sig
      call compute_u_kq(ik, q)
      ! compute overlaps
      if (noncolin) then
     
         CALL zgemm('C', 'N', nbnd_occ (ik), nbnd_occ (ik), npwx*npol, (1.d0,0.d0), evc(1,1), &
                    npwx*npol, evq(1,1), npwx*npol, (0.d0,0.d0), overlap(1,1), nbnd_occ(ik))
      else

         CALL zgemm('C', 'N', nbnd_occ (ik), nbnd_occ (ik), npw, (1.d0,0.d0), evc(1,1), &
                    npwx, evq(1,1), npwx, (0.d0,0.d0), overlap(1,1), nbnd_occ(ik))
      endif
      !do ibnd = 1, nbnd_occ(ik)
      !  do jbnd = 1, nbnd_occ(ik)
      !    overlap(ibnd,jbnd) = zdotc(npwx*npol, evc(1,ibnd), 1, evq(1,jbnd), 1)
      ! enddo
      !enddo
    
#ifdef __MPI
  call mp_sum(overlap, intra_pool_comm)
#endif
      !!write(*,'(5X,''inverse: ik='',I4,'' ipol='',I1,'' sig='',I2)') &
      !!      ik, ipol, sig
      call invmat(nbnd_occ(ik), overlap)

      ! compute the covariant derivative
      !!write(*,'(5X,''dual   : ik='',I4,'' ipol='',I1,'' sig='',I2)') &
      !!      ik, ipol, sig
      if (noncolin) then
         CALL zgemm( 'N', 'N', npwx*npol, nbnd_occ(ik), nbnd_occ(ik), (1.d0,0.d0)*sig*0.5d0/(delta_k*tpiba), &
         evq(1,1), npwx*npol, overlap(1,1), nbnd_occ(ik), (1.d0,0.d0), evc1(1,1,ipol), npwx*npol )
      else
         CALL zgemm( 'N', 'N', npw, nbnd_occ(ik), nbnd_occ(ik), (1.d0,0.d0)*sig*0.5d0/(delta_k*tpiba), &
         evq(1,1), npwx, overlap(1,1), nbnd_occ(ik), (1.d0,0.d0), evc1(1,1,ipol), npwx )
      endif
      !do ibnd = 1, nbnd_occ(ik)
      !  do jbnd = 1, nbnd_occ(ik)
      !    dudk(1:npwx*npol,ibnd,ipol) = dudk(1:npwx*npol,ibnd,ipol) + &
      !                  sig * 0.5d0/(delta_k*tpiba) * &
      !                  overlap(jbnd,ibnd) * evq(1:npwx*npol,jbnd)
      !  enddo
      !enddo
    enddo ! sig
  enddo ! ipol
  deallocate(overlap)
END SUBROUTINE dudk_covariant

SUBROUTINE dudk_sos(ik)

   USE kinds,                       ONLY : DP
   USE wvfct,                       ONLY : nbnd, et, npwx
   USE wavefunctions,               ONLY : evc
   USE noncollin_module,            ONLY : noncolin, npol
   USE orbm_module,                 ONLY : vel_evc, evc1, ry2ha
   USE mp_pools,                    ONLY : intra_pool_comm
   USE mp,                          ONLY : mp_sum
   USE klist,                       ONLY : ngk

   IMPLICIT none
   INTEGER :: ik, ipol
   ! local
   INTEGER :: npw, ibnd, jbnd
   COMPLEX(DP), ALLOCATABLE :: ps(:,:)
   REAL(DP), PARAMETER :: delta = 1.d-5
   
   ALLOCATE( ps(nbnd, nbnd) )
   
   npw = ngk(ik)

   DO ipol = 1,3
   
   if (noncolin) then
     
      CALL zgemm('C', 'N', nbnd, nbnd, npwx*npol, (1.d0,0.d0), evc(1,1), &
                    npwx*npol, vel_evc(1,1,ipol), npwx*npol, (0.d0,0.d0), ps(1,1), nbnd)
   else
      CALL zgemm('C', 'N', nbnd, nbnd, npw, (1.d0,0.d0), evc(1,1), &
                    npwx, vel_evc(1,1,ipol), npwx, (0.d0,0.d0), ps(1,1), nbnd)
   endif
 
#ifdef __MPI
   call mp_sum(ps, intra_pool_comm)
#endif

   DO ibnd =  1, nbnd
      DO jbnd = 1, nbnd
         IF ( abs(et(ibnd,ik)-et(jbnd,ik)) < delta )  THEN
            ps(jbnd, ibnd) = (0.d0,0.d0)
         ELSE
            ps(jbnd, ibnd) = ps(jbnd, ibnd)/(et(ibnd,ik)-et(jbnd,ik))/ry2ha
         ENDIF
      ENDDO
   ENDDO
   
   if (noncolin) then
   
       CALL zgemm( 'N', 'N', npwx*npol, nbnd, nbnd, (1.d0,0.d0), &
         evc(1,1), npwx*npol, ps(1,1), nbnd, (0.d0,0.d0), evc1(1,1,ipol), npwx*npol )
   else
       CALL zgemm( 'N', 'N', npw, nbnd, nbnd, (1.d0,0.d0), &
         evc(1,1), npwx, ps(1,1), nbnd, (0.d0,0.d0), evc1(1,1,ipol), npwx )
         
   endif

   ENDDO
   
   DEALLOCATE(ps)
   
END SUBROUTINE dudk_sos

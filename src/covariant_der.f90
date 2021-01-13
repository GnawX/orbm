  SUBROUTINE covariant_der(ik)
  ! calculate first order wave function evc1
  ! first order velocity or inverse effective mass invmass
  USE kinds,                ONLY : dp  
  USE cell_base,            ONLY : tpiba
  USE wvfct,                ONLY : nbnd, npwx
  USE wavefunctions,        ONLY : evc
  USE klist,                ONLY : xk, ngk
  USE noncollin_module,     ONLY : npol, noncolin
  USE mp_pools,             ONLY : intra_pool_comm
  USE mp,                   ONLY : mp_sum
  USE matrix_inversion
  USE orbm_module         
  IMPLICIT NONE
  complex(dp), external :: zdotc 
  complex(dp), allocatable :: overlap(:,:), vel(:,:,:)
  real(dp) :: q(3), dk, xkold(3)
  integer :: ik, i, sig, ibnd, jbnd
  integer :: ipol, jpol, npw
  
  npw = ngk(ik)
  ! allocate overlap matrix 
  allocate(overlap(nbnd_occ(ik), nbnd_occ(ik)))
  allocate(vel(npwx*npol,nbnd,3))
  allocate(invmass(nbnd,nbnd,3,3))
  !allocate(evc0(npwx*npol,nbnd))

  dk = q_orbm/2.d0
  !q=0.d0
  
  !call compute_u_kq(ik,q)
  !evc0 = evq

  evc1(:,:,:) = (0.d0,0.d0)
  invmass(:,:,:,:) = (0.d0,0.d0)
  ! loop over crystal directions
  do ipol = 1, 3
    !evc1(:,:,ipol) = (0.d0,0.d0)
    
    ! loop over +/-1
    do sig = -1, 1, 2

      ! set the k-point
      q(:) = 0.d0
      q(ipol) = dk * sig

      call compute_u_kq(ik, q)
      ! compute overlaps
      if (noncolin) then
     
         CALL zgemm('C', 'N', nbnd_occ (ik), nbnd_occ (ik), npwx*npol, (1.d0,0.d0), evc(1,1), &
                    npwx*npol, evq(1,1), npwx*npol, (0.d0,0.d0), overlap(1,1), nbnd_occ(ik))
      else

         CALL zgemm('C', 'N', nbnd_occ (ik), nbnd_occ (ik), npw, (1.d0,0.d0), evc(1,1), &
                    npwx, evq(1,1), npwx, (0.d0,0.d0), overlap(1,1), nbnd_occ(ik))
      endif

    
#ifdef __MPI
  call mp_sum(overlap, intra_pool_comm)
#endif

      call invmat(nbnd_occ(ik), overlap)

      ! compute the covariant derivative

      if (noncolin) then
         CALL zgemm( 'N', 'N', npwx*npol, nbnd_occ(ik), nbnd_occ(ik), (1.d0,0.d0)*sig*0.5d0/(dk*tpiba), &
         evq(1,1), npwx*npol, overlap(1,1), nbnd_occ(ik), (1.d0,0.d0), evc1(1,1,ipol), npwx*npol )
      else
         CALL zgemm( 'N', 'N', npw, nbnd_occ(ik), nbnd_occ(ik), (1.d0,0.d0)*sig*0.5d0/(dk*tpiba), &
         evq(1,1), npwx, overlap(1,1), nbnd_occ(ik), (1.d0,0.d0), evc1(1,1,ipol), npwx )
      endif
      
      ! derivative of velocity matrix elements
      xkold(:) = xk(:,ik)
      xk(:,ik) = xk(:,ik) + q(:)
      
      vel = 0.d0
      do jpol = 1,3
      
         call apply_vel(evq, vel(1,1,jpol), ik, jpol)
      
         if (noncolin) then
            CALL zgemm('C', 'N', nbnd, nbnd, npwx*npol, (1.d0,0.d0)*sig*0.5d0/(dk*tpiba), evq(1,1), &
                    npwx*npol, vel(1,1,jpol), npwx*npol, (1.d0,0.d0), invmass(1,1,ipol,jpol), nbnd)
         else

            CALL zgemm('C', 'N', nbnd, nbnd, npw, (1.d0,0.d0)*sig*0.5d0/(dk*tpiba), evq(1,1), &
                    npwx, vel(1,1,jpol), npwx, (1.d0,0.d0), invmass(1,1,ipol,jpol), nbnd)
         endif
         
      enddo
      
      xk(:,ik) = xkold(:)

    enddo ! sig
  enddo ! ipol
  
  deallocate(overlap,vel)
  
END SUBROUTINE covariant_der

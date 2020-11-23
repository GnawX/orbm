  SUBROUTINE dudk_covariant(ik, dudk)
  USE kinds,                ONLY : dp  
  USE cell_base,            ONLY : tpiba
  USE wvfct,                ONLY : npw, nbnd, npwx
  USE wavefunction,         ONLY : evc
  USE klist,                ONLY : xk
  USE noncollin_module,     ONLY : npol, noncolin
  !USE io_global,            ONLY : stdout
  !USE buffers,              ONLY : get_buffer, save_buffer
  USE matrix_inversion
  USE orbm_module         
  IMPLICIT NONE
  
  complex(dp), allocatable :: overlap(:,:)
  complex(dp) :: dudk(npwx*npol,nbnd,3)
  real(dp) :: q(3), delta_k
  integer :: ik, i, sig, 
  integer :: ipol
  
  ! allocate overlap matrix 
  allocate(overlap(nbnd_occ(ik), nbnd_occ(ik)))

  delta_k = q_orbm/2.d0

  ! loop over crystal directions
  do ipol = 1, 3
    dudk(:,:,:) = (0.d0,0.d0)
    
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
         evq(1,1), npwx*npol, overlap(1,1), nbnd_occ(ik), (1.d0,0.d0), dudk(1,1,ipol), npwx*npol )
      else
         CALL zgemm( 'N', 'N', npw, nbnd_occ(ik), nbnd_occ(ik), (1.d0,0.d0)*sig*0.5d0/(delta_k*tpiba), &
         evq(1,1), npwx, overlap(1,1), nbnd_occ(ik), (1.d0,0.d0), dudk(1,1,ipol), npwx )
      endif
      
    enddo ! sig
  enddo ! ipol
  deallocate(overlap)
END SUBROUTINE dudk_covariant

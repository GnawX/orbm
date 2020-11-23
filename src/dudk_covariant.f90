  SUBROUTINE dudk_covariant(ik, occ, dudk)
  USE kinds,                ONLY : dp  
  USE cell_base,            ONLY : bg, tpiba2, tpiba
  USE wvfct,                ONLY : npw, nbnd, npwx
  USE wavefunctions_module, ONLY : evc
  USE klist,                ONLY : xk
  USE io_files,             ONLY : nwordwfc, iunwfc
  USE io_global,            ONLY : stdout
  USE buffers,              ONLY : get_buffer, save_buffer
  USE g_tensor_module,      ONLY : q_conv_thr, q_gipaw
  IMPLICIT NONE
  !real(dp), parameter :: delta_k = 0.05_dp
  complex(dp), external :: zdotc
  complex(dp), allocatable :: overlap(:,:), evc0(:,:)
  complex(dp) :: dudk(npwx,nbnd,3)
  real(dp) :: q_gipaw3(3), ethr_q, delta_k
  integer :: ik, i, occ, sig, iter
  integer :: ibnd, jbnd, ipol

  if (occ == 0) return

  ! allocate overlap matrix and evc0
  allocate(overlap(occ,occ), evc0(npwx,nbnd))

  ! read the wavefunction
  call get_buffer(evc, nwordwfc, iunwfc, ik)
  q_gipaw3(:) = 0.d0
  ethr_q = q_conv_thr

  call compute_u_kq(ik, q_gipaw3, iter, ethr_q)
  call save_buffer(evc, nwordwfc, iunwfc, ik)
  evc0(1:npwx,1:nbnd) = evc(1:npwx,1:nbnd)
  !delta_k = q_gipaw/2.d0/tpiba
  delta_k = q_gipaw/tpiba

  ! loop over crystal directions
  do ipol = 1, 3
    dudk(1:npwx,1:occ,ipol) = (0.d0,0.d0)
    
    ! loop over +/-1
    do sig = -1, 1, 2

      ! set the k-point
      q_gipaw3(:) = 0.d0
      q_gipaw3(ipol) = delta_k * sig
      !!write(*,'(5X,''overlap: ik='',I4,'' ipol='',I1,'' sig='',I2)') &
      !!      ik, ipol, sig
      evc(1:npwx,1:nbnd) = evc0(1:npwx,1:nbnd)
      ethr_q = q_conv_thr
      call compute_u_kq(ik, q_gipaw3, iter, ethr_q)

      ! compute overlaps
      do ibnd = 1, occ
        do jbnd = 1, occ
          overlap(ibnd,jbnd) = zdotc(npw, evc0(1,ibnd), 1, evc(1,jbnd), 1)
       enddo
      enddo
#ifdef __PARA
      call reduce(2*occ*occ, overlap)
#endif

      !!write(*,'(5X,''inverse: ik='',I4,'' ipol='',I1,'' sig='',I2)') &
      !!      ik, ipol, sig
      call invert_matrix(occ, overlap)

      ! compute the covariant derivative
      !!write(*,'(5X,''dual   : ik='',I4,'' ipol='',I1,'' sig='',I2)') &
      !!      ik, ipol, sig
      do ibnd = 1, occ
        do jbnd = 1, occ
          dudk(1:npw,ibnd,ipol) = dudk(1:npw,ibnd,ipol) + &
                        sig * 0.5d0/(delta_k*tpiba) * &
                        overlap(jbnd,ibnd) * evc(1:npw,jbnd)
        enddo
      enddo

    enddo ! sig
  enddo ! ipol
  deallocate(overlap, evc0)
END SUBROUTINE dudk_covariant

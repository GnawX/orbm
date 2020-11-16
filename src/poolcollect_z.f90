SUBROUTINE poolcollect_z( length, nks, f_in, nkstot, f_out )
  !----------------------------------------------------------------------------
  !! Collects a real array f_in, distributed across pools, from all pools,
  !! into a real array f_out.
  !
  !! * On input: f_in(length,nks) contains data for the "nks" k-points
  !!   of the current pool, on all pools;
  !! * On output: f_out(length,nkstot) contains data for all "nkstot" k-points
  !!   on all pools.
  !
  !! f_in and f_out must differ! Honors "kunit"
  !
  USE kinds,     ONLY : DP
  USE mp_pools,  ONLY : my_pool_id, npool, kunit, &
                        inter_pool_comm, intra_pool_comm
  USE mp,        ONLY : mp_sum
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: length
  !! first dimension of arrays
  INTEGER, INTENT(IN) :: nks
  !! number of k-points per pool
  INTEGER, INTENT(IN) :: nkstot
  !! total number of k-points
  COMPLEX(DP), INTENT(IN) :: f_in(length,length,nks,3)
  !! pool-distributed function
  COMPLEX(DP), INTENT(OUT) :: f_out(length,length,nkstot,3)
  !! pool-collected function
  !
  ! ... local variables
  !
  INTEGER :: nbase, rest, nks1
  !
  nks1 = kunit * ( nkstot / kunit / npool )
  !
  rest = ( nkstot - nks1 * npool ) / kunit
  !
  IF ( ( my_pool_id + 1 ) <= rest ) nks1 = nks1 + kunit
  !
  IF (nks1 /= nks) CALL errore( 'xk_collect', 'inconsistent number of k-points', 1 )
  !
  ! ... calculates nbase = the position in the list of the first point that
  ! ...                    belong to this npool - 1
  !
  nbase = nks * my_pool_id
  !
  IF ( ( my_pool_id + 1 ) > rest ) nbase = nbase + rest * kunit
  !
  ! copy the original points in the correct position of the list
  !
  f_out = (0.0_DP, 0.0_DP)
  f_out(:,:,nbase+1:nbase+nks,:) = f_in(:,:,1:nks,:)
  !
  CALL mp_sum( f_out, inter_pool_comm )
  !
  !
  RETURN
  !
END SUBROUTINE poolcollect_z

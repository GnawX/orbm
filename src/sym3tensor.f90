   SUBROUTINE symmatrix3( mat3 )
     !-----------------------------------------------------------------------
     !! Symmetrize a function \(f(i,j,k)\) (e.g. nonlinear susceptibility),
     !! where \(i,j,k\) are the cartesian components.  
     !! BEWARE: input in crystal axis, output in cartesian axis.
     !
     IMPLICIT NONE
     !
     REAL(DP), INTENT(INOUT) :: mat3(3,3,3)
     !! function f(i,j,k) to symmetrize
     !
     ! ... local variables
     !
     INTEGER :: isym, i,j,k,l,m,n
     REAL(DP) :: work(3,3,3)
     !
     IF (nsym > 1) THEN
        !
        work (:,:,:) = 0.0_dp
        DO isym = 1, nsym
           DO i = 1, 3
              DO j = 1, 3
                 DO k = 1, 3
                    DO l = 1, 3
                       DO m = 1, 3
                          DO n = 1, 3
                             work (i, j, k) = work (i, j, k) + &
                                s (i, l, isym) * s (j, m, isym) * &
                                s (k, n, isym) * mat3 (l, m, n)
                          END DO
                       END DO
                    END DO
                 END DO
              END DO
           END DO
        END DO
        mat3 = work/ DBLE(nsym)
        !
     END IF
     !
     ! Bring to cartesian axis
     !
     CALL crys_to_cart_mat3 ( mat3 ) 
     !
   END SUBROUTINE symmatrix3

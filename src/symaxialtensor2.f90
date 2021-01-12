   SUBROUTINE symaxialtensor2( matr )
     !-----------------------------------------------------------------------
     USE kinds,       ONLY : DP
     USE symm_base,   ONLY : sr, nsym 
     !
     IMPLICIT NONE
     !
     REAL(DP), INTENT(INOUT) :: matr(3,3)
     !! the function \(f(i,j)\) to symmetrize
     !
     ! ... local variables
     !
     INTEGER :: isym, i,j,k,l
     REAL(DP) :: work (3,3), det
     !
     IF (nsym == 1) RETURN
     !
     ! bring matrix to crystal axis
     !
     !CALL cart_to_crys ( matr )
     !
     ! symmetrize in crystal axis
     !
     work (:,:) = 0.0_dp
     DO isym = 1, nsym
        det = sr(1,1,isym)*sr(2,2,isym)*sr(3,3,isym) + &
              sr(1,2,isym)*sr(2,3,isym)*sr(3,1,isym) + &
              sr(1,3,isym)*sr(2,1,isym)*sr(3,2,isym) - &
              sr(1,3,isym)*sr(2,2,isym)*sr(3,1,isym) - &  
              sr(1,1,isym)*sr(2,3,isym)*sr(3,2,isym) - &  
              sr(1,2,isym)*sr(2,1,isym)*sr(3,3,isym)   
        DO i = 1, 3
           DO j = 1, 3
              DO k = 1, 3
                 DO l = 1, 3
                    work (i,j) = work (i,j) + &
                       det * sr (i,k,isym) * sr (j,l,isym) * matr (k,l)
                 END DO
              END DO
           END DO
        END DO
     END DO
     matr (:,:) = work (:,:) / DBLE(nsym)
     !
     ! bring matrix back to cartesian axis
     !
     !CALL crys_to_cart ( matr )
     !
   END SUBROUTINE symaxialtensor2

!   SUBROUTINE sym3tensor( mat3 )
!     !-----------------------------------------------------------------------
!     !! Symmetrize a function \(f(i,j,k)\) (e.g. nonlinear susceptibility),
!     !! where \(i,j,k\) are the cartesian components.  
!     !
!     USE kinds,       ONLY : DP
!     USE symm_base,   ONLY : s, nsym 
!     IMPLICIT NONE
!     !
!     REAL(DP), INTENT(INOUT) :: mat3(3,3,3)
!     !! function f(i,j,k) to symmetrize
!     !
!     ! ... local variables
!     !
!     INTEGER :: isym, i,j,k,l,m,n
!     REAL(DP) :: work(3,3,3)
!     !
!     IF (nsym == 1) RETURN
!     
!     CALL cart_to_crys_m3 ( mat3 )
!        !
!        work (:,:,:) = 0.0_dp
!        DO isym = 1, nsym
!           DO i = 1, 3
!              DO j = 1, 3
!                 DO k = 1, 3
!                    DO l = 1, 3
!                       DO m = 1, 3
!                          DO n = 1, 3
!                             work (i, j, k) = work (i, j, k) + &
!                                s (i, l, isym) * s (j, m, isym) * &
!                                s (k, n, isym) * mat3 (l, m, n)
!                          END DO
!                       END DO
!                    END DO
!                 END DO
!              END DO
!           END DO
!        END DO
!        mat3 = work/ DBLE(nsym)
!        !
!     !
!     ! Bring to cartesian axis
!     !
!     CALL crys_to_cart_m3 ( mat3 ) 
!     !
!   END SUBROUTINE sym3tensor
!   
!   SUBROUTINE cart_to_crys_m3( mat3 )
!     !-----------------------------------------------------------------------
!     !! Crystal to cartesian axis conversion for \(f(i,j,k)\) objects.
!     !
!     USE kinds,       ONLY : DP
!     USE cell_base,   ONLY : at
!     IMPLICIT NONE
!     !
!     REAL(DP), INTENT(INOUT) :: mat3(3,3,3)
!     !! Axis conversion tensor
!     !
!     REAL(DP) :: work(3,3,3)
!     INTEGER :: i,j,k,l,m,n
!     !
!     ! ... local variables
!     !
!     work(:,:,:) = 0.0_dp
!     DO i = 1, 3
!        DO j = 1, 3
!           DO k = 1, 3
!              DO l = 1, 3
!                 DO m = 1, 3
!                    DO n = 1, 3
!                       work (i, j, k) = work (i, j, k) +  &
!                          mat3 (l, m, n) * at (l, i) * at (m, j) * at (n, k)
!                    END DO
!                 END DO
!              END DO
!           END DO
!        END DO
!     END DO
!     mat3(:,:,:) = work (:,:,:)
!     !
!   END SUBROUTINE cart_to_crys_m3
!   
!   SUBROUTINE crys_to_cart_m3( mat3 )
!     !-----------------------------------------------------------------------
!     !! Crystal to cartesian axis conversion for \(f(i,j,k)\) objects.
!     !
!     USE kinds,       ONLY : DP
!     USE cell_base,   ONLY : bg
!     IMPLICIT NONE
!     !
!     REAL(DP), INTENT(INOUT) :: mat3(3,3,3)
!     !! Axis conversion tensor
!     !
!     REAL(DP) :: work(3,3,3)
!     INTEGER :: i,j,k,l,m,n
!     !
!     ! ... local variables
!     !
!     work(:,:,:) = 0.0_dp
!     DO i = 1, 3
!        DO j = 1, 3
!           DO k = 1, 3
!              DO l = 1, 3
!                 DO m = 1, 3
!                    DO n = 1, 3
!                       work (i, j, k) = work (i, j, k) +  &
!                          mat3 (l, m, n) * bg (i, l) * bg (j, m) * bg (k, n)
!                    END DO
!                 END DO
!              END DO
!           END DO
!        END DO
!     END DO
!     mat3(:,:,:) = work (:,:,:)
!     !
!   END SUBROUTINE crys_to_cart_m3

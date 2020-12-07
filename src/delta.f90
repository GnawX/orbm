  FUNCTION delta(x) RESULT(y)
     USE kinds,       ONLY : DP
     USE orbm_module, ONLY : sigma
     USE constants,   ONLY : pi
     IMPLICIT NONE
     REAL(DP), INTENT(IN) :: x
     REAL(DP) :: y
     y = EXP(-(x/sigma)**2/2)/sigma/SQRT(2*pi)
  END FUNCTION

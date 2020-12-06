  FUNCTION delta(x) RESULT(y)
     USE orbm_module
     IMPLICIT NONE
     REAL(DP), INTENT(IN) :: x
     REAL(DP) :: y
     y = EXP(-(x/sigma)**2/2)/sigma/SQRT(2*pi)
  END FUNCTION

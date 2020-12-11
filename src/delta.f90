  FUNCTION delta(x) RESULT(y)
     USE kinds,       ONLY : DP
     USE orbm_module, ONLY : sigma, smear
     USE constants,   ONLY : pi
     IMPLICIT NONE
     REAL(DP), INTENT(IN) :: x
     REAL(DP) :: y
     if (trim(smear) == 'gauss') then
        y = EXP(-(x/sigma)**2)/sigma/SQRT(pi)
     elseif (trim(smear) == 'lorenz') then
        y = sigma/pi/(x*x + sigma*sigma)
     endif
  END FUNCTION

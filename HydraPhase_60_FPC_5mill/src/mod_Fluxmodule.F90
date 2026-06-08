MODULE FluxModule
  USE GlobalVariables
  IMPLICIT NONE

CONTAINS

  SUBROUTINE CallFluxSolver()
    IF     (myflux == 1) THEN
      CALL hllczien()
    ELSEIF (myflux == 2) THEN
      CALL hllcbatten()
    ELSEIF (myflux == 3) THEN
      CALL jameson()
    ELSE   
      PRINT*, "Invalid flux method!"
      STOP
    END IF
#ifdef FIVE_PHASE
  CALL jamesonqrs() 
#endif
  END SUBROUTINE CallFluxSolver

END MODULE FluxModule


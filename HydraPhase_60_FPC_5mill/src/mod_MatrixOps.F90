MODULE MatrixOps
  !$acc routine(det) seq
  CONTAINS
  FUNCTION det(A11, A21, A31, A12, A22, A32, A13, A23, A33) RESULT(res)
    DOUBLE PRECISION, INTENT(IN) :: A11, A21, A31, A12, A22, A32, A13, A23, A33
    DOUBLE PRECISION :: res

    res = A11 * (A22 * A33 - A23 * A32) - A21 * (A12 * A33 - A13 * A32) + A31 * (A12 * A23 - A13 * A22)
  END FUNCTION det
  
 !$acc routine(fvenkat) seq
  FUNCTION fvenkat(delp, delm, epssq) RESULT(res)
    DOUBLE PRECISION, INTENT(IN) :: delp, delm, epssq
    DOUBLE PRECISION :: res

    res = (((delp * delp + epssq) * delm + 2.0d0 * delm * delm * delp) / &
          (delp * delp + 2.0d0 * delm * delm + delp * delm + epssq)) / delm
  END FUNCTION fvenkat

END MODULE MatrixOps



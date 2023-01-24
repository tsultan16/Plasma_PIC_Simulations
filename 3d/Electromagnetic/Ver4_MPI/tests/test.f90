PROGRAM array_test

IMPLICIT NONE


INTEGER, PARAMETER :: nx = 20 
INTEGER :: i


DO i = 3, nx-2, nx-5
    PRINT*,'i =',i
END DO






END PROGRAM array_test
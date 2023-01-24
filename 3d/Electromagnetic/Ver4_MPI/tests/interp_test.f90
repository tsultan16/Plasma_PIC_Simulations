PROGRAM interpolation_test

IMPLICIT NONE


REAL(8) :: x, xj0, xj, xj1, xj2
INTEGER :: ix

ix = 1
x = 1.9d0

xj = S1_1D(ix-x)
xj1 = S1_1D((ix+1)-x)

PRINT*,''
PRINT*,'1st order Interpolation: x_j, x_j+1 = ',xj, xj1
PRINT*,''

xj0 = S2_1D((ix-1)-x)
xj = S2_1D(ix-x)
xj1 = S2_1D((ix+1)-x)
xj2 = S2_1D((ix+2)-x)

PRINT*,''
PRINT*,'2nd order Interpolation: x_j-1, x_j, x_j+1, x_j+2 = ',xj0,xj,xj1,xj2
PRINT*,'Sum = ', xj0+xj+xj1+xj2
PRINT*,''


CONTAINS




! 1D particle 1st order shape function (piecewise linear)
FUNCTION S1_1D(x) RESULT (sx)

    REAL*8, INTENT(IN) :: x
    REAL*8 :: sx

    IF(ABS(x) .LT. 1.d0) THEN
        sx = 1.d0 - ABS(x)     
    ELSE    
        sx = 0.d0
    END IF

END FUNCTION S1_1D


! 1D particle 2nd order shape function (quadratic spline)
FUNCTION S2_1D(x) RESULT (sx)

    REAL*8, INTENT(IN) :: x
    REAL*8 :: sx

    IF(ABS(x) .LT. 0.5d0) THEN
        sx = (3.d0/4.d0) - x**2     
    ELSE IF(ABS(x) .GE. 0.5d0 .AND. ABS(x) .LE. 1.5d0)THEN
        sx = 0.5d0 * ( 1.5d0 - ABS(x))**2
    ELSE    
        sx = 0.d0
    END IF

END FUNCTION S2_1D





END PROGRAM interpolation_test

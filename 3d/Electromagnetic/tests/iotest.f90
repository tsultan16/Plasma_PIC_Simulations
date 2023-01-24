PROGRAM iotest

IMPLICIT NONE

REAL*8 :: A(2,2,2), B(2,2,2)


A(1,1,1) = 1
A(2,1,1) = 2
A(1,2,1) = 3
A(2,2,1) = 4

A(1,1,2) = 5
A(2,1,2) = 6
A(1,2,2) = 7
A(2,2,2) = 8


!CALL writefile()


CALL readfile()

PRINT*,' B = '
PRINT*,B(1:2,1,1)
PRINT*,B(1:2,2,1)
PRINT*,B(1:2,1,2)
PRINT*,B(1:2,2,2)
PRINT*,''


CONTAINS



SUBROUTINE writefile()

    OPEN(UNIT=1, FILE='output.dat', FORM = 'UNFORMATTED')
    
    write(1) A

    CLOSE(UNIT=1)    


END SUBROUTINE writefile




SUBROUTINE readfile()

    OPEN(UNIT=2, FILE='output.dat', FORM = 'UNFORMATTED', STATUS='OLD')
    
    READ(2) B

    CLOSE(UNIT=2) 


END SUBROUTINE readfile




END PROGRAM iotest
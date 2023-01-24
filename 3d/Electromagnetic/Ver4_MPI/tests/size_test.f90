PROGRAM size_test
IMPLICIT NONE


REAL(8) :: A(8)


A =0.d0

PRINT*,'Memory(bytes) = ',SIZEOF(A)


END PROGRAM size_test
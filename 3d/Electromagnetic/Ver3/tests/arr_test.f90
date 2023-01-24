PROGRAM arr_test

IMPLICIT NONE


REAL*8, ALLOCATABLE :: a(:,:,:)
INTEGER :: i

ALLOCATE(a(2,2,2))

a(1,1,1) = 1
a(2,1,1) = 2
a(1,2,1) = 3
a(2,2,1) = 4
a(1,1,2) = 5
a(2,1,2) = 6
a(1,2,2) = 7
a(2,2,2) = 8



CALL see_1d(a)


DEALLOCATE(a)

CONTAINS


SUBROUTINE see_1d(b)

    REAL*8, INTENT(INOUT) :: b(1)
    INTEGER :: i
    
    DO i = 1, 2*2*2
       WRITE(*,FMT='( "b(",i3, " ) = ",f5.1)') i, b(i)
    END DO    

END SUBROUTINE see_1d



END PROGRAM arr_test
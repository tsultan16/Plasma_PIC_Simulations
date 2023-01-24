MODULE domain_mod

USE constants_mod
USE data_mod
USE particleMover_mod
USE ISO_C_BINDING

IMPLICIT NONE

INCLUDE 'mpif.h'



! MPI variables
INTEGER :: status(MPI_STATUS_SIZE)


! RMA variables
INTEGER :: sizedouble, winsize
INTEGER :: win 
TYPE(C_PTR) :: rma_cmem_fields  ! C pointer to start of RMA window memory block
REAL(8), POINTER :: mem_fields(:) ! Fortran pointer that will be associated with the RMA window's C pointer



! Notes: Need to be careful when handling 1d and 2d domain decompositions with the same subroutines. 
! Better to just make separate subroutines for that...
!


CONTAINS

!################
! General setup ! 
!################

! top-level MPI initialization subroutine
SUBROUTINE initialize_mpi()
    
    INTEGER :: i, nwin
        
    ! Sanity checks
    IF((nranks_y .GT. 1 .AND. ndims .EQ. 1) .OR. (nranks_y .EQ. 1 .AND. ndims .EQ. 2) ) THEN
        PRINT*,'ERROR. Need ndims = 2 for nranks_y > 1.'
        STOP
    END IF
    
    PRINT*,''
    
    ! initialize
    CALL MPI_INIT(ierr)

    ! Get total number of processes
    CALL MPI_COMM_SIZE(MPI_COMM_WORLD, numprocs(1), ierr)

    ! set up domain decomposition
    CALL setup_domain()

    
    ! allocate memory for MPI buffers and RMA memory windows
    CALL setup_rma()

   

PRINT*,''
PRINT*,'My rank, coordinate =',myrank, mycoord
PRINT*,''


END SUBROUTINE initialize_mpi



SUBROUTINE setup_rma()

    INTEGER :: i, nwin

 
    ! Compute buffer size needed for field data, two-layers deep (choose from max of x,y boundary sizes).
    ! and particle data.
    
    field_buffer_size_x = 2 * nvars_fields * (ny+5) * (nz+5)
    field_buffer_size_y = 2 * nvars_fields * (nx+5) * (nz+5) 
    
    ! allocate field data buffer
    ALLOCATE(field_buffer_xm(1:field_buffer_size_x),field_buffer_xp(1:field_buffer_size_x)) 
    IF(ndims .EQ. 2) ALLOCATE(field_buffer_ym(1:field_buffer_size_y),field_buffer_yp(1:field_buffer_size_y)) 
    
   
   ! determine window size in bytes    
   CALL MPI_SIZEOF(DBLE(1), sizedouble, ierr) ! size in bytes per slot
   nwin = 2*field_buffer_size_x  ! total # of required slots
   
   IF(ndims .EQ. 2) nwin = nwin + 2*field_buffer_size_y  ! in 2D, need extra space to store ghost cells coming from neighbors in y direction
   
   winsize = nwin * sizedouble  ! total # of bytes required for the field data RMA window

   IF(myrank .EQ. 0) PRINT*,'sizedouble, RMA winsize = ',sizedouble, winsize 

   ! rma memory window allocation 
   CALL MPI_WIN_ALLOCATE(winsize, sizedouble, MPI_INFO_NULL, comm2d, rma_cmem_fields, win, ierr)
   CALL C_F_POINTER(rma_cmem_fields, mem_fields, (/ nwin /)) ! associate C pointer returned by win_allocate to our fortran pointer
   
   ! clear my rma window
   mem_fields = 0.d0


END SUBROUTINE setup_rma


! domain decomposition subroutine
SUBROUTINE setup_domain()

    ! set up cart communicator properties
    IF(bndry .EQ. 1) THEN
        isperiodic = .TRUE.          
    ELSE IF(bndry .EQ. 2) THEN
        isperiodic = .FALSE.   
    END IF
    
	reorder = .TRUE. ! allow automatic rank reordering 
	dims = (/ nranks_x, nranks_y /) 
	
    ! create cart communicator handle
	CALL MPI_CART_CREATE(MPI_COMM_WORLD, ndims, dims, isperiodic, reorder, comm2d, ierr)

! get my rank 
    CALL MPI_COMM_RANK(comm2d, myrank, ierr)

    ! get my rank co-ordinates
    CALL MPI_CART_COORDS(comm2d, myrank, ndims, mycoord, ierr)

    ! get neighbor ranks in x and y directions
    !neighbor_rank = -2
    CALL MPI_CART_SHIFT(comm2d, 0, 1, neighbor_rank(1), neighbor_rank(2), ierr) 
    !IF(ndims .EQ. 2 ) CALL MPI_CART_SHIFT(comm2d, 1, 1, neighbor_rank(3), neighbor_rank(4), ierr)
    CALL MPI_CART_SHIFT(comm2d, 1, 1, neighbor_rank(3), neighbor_rank(4), ierr)

       
    ! find diagonal neighbor ranks
    IF(ndims .EQ. 2) THEN
  
        IF(mycoord(1) .GT. 0 .AND. mycoord(1) .LT. nranks_x-1 .AND. &
        mycoord(2) .GT. 0 .AND. mycoord(2) .LT. nranks_y-1 ) THEN
    
            neighbor_rank(5) = neighbor_rank(1) - 1 ! bottom left
            neighbor_rank(7) = neighbor_rank(1) + 1 ! top left
            neighbor_rank(6) = neighbor_rank(2) - 1 ! bottom right
            neighbor_rank(8) = neighbor_rank(2) + 1 ! top right
    
    
         ELSE IF(mycoord(1) .EQ. 0 .AND.  &
         mycoord(2) .GT. 0 .AND. mycoord(2) .LT. nranks_y-1 ) THEN
         
            neighbor_rank(5) = -2 ! bottom left
            neighbor_rank(7) = -2 ! top left
            neighbor_rank(6) = neighbor_rank(2) - 1 ! bottom right
            neighbor_rank(8) = neighbor_rank(2) + 1 ! top right
    
         ELSE IF(mycoord(1) .EQ. nranks_x-1 .AND.  &
         mycoord(2) .GT. 0 .AND. mycoord(2) .LT. nranks_y-1 ) THEN
        
    
            neighbor_rank(5) = neighbor_rank(1) - 1 ! bottom left
            neighbor_rank(7) = neighbor_rank(1) + 1 ! top left
            neighbor_rank(6) = -2 ! bottom right
            neighbor_rank(8) = -2 ! top right
    
         ELSE IF(mycoord(1) .GT. 0 .AND. mycoord(1) .LT. nranks_x-1 .AND. &
         mycoord(2) .EQ. 0 ) THEN
    
            neighbor_rank(5) = -2 ! bottom left
            neighbor_rank(7) = neighbor_rank(1) + 1 ! top left
            neighbor_rank(6) = -2 ! bottom right
            neighbor_rank(8) = neighbor_rank(2) + 1 ! top right
    
        
         ELSE IF(mycoord(1) .GT. 0 .AND. mycoord(1) .LT. nranks_x-1 .AND. &
         mycoord(2) .EQ. nranks_y-1 ) THEN
    
            neighbor_rank(5) = neighbor_rank(1) - 1 ! bottom left
            neighbor_rank(7) = -2 ! top left
            neighbor_rank(6) = neighbor_rank(2) - 1 ! bottom right
            neighbor_rank(8) = -2 ! top right
    
       ELSE IF(mycoord(1) .EQ. 0 .AND. mycoord(2) .EQ. 0 ) THEN
         
            neighbor_rank(5) = -2 ! bottom left
            neighbor_rank(7) = -2 ! top left
            neighbor_rank(6) = -2 ! bottom right
            neighbor_rank(8) = neighbor_rank(2) + 1 ! top right
    
         ELSE IF(mycoord(1) .EQ. nranks_x-1 .AND. mycoord(2) .EQ. nranks_y-1 ) THEN
        
    
            neighbor_rank(5) = neighbor_rank(1) - 1 ! bottom left
            neighbor_rank(7) = -2 ! top left
            neighbor_rank(6) = -2 ! bottom right
            neighbor_rank(8) = -2 ! top right
    
         ELSE IF(mycoord(1) .EQ. nranks_x-1 .AND. mycoord(2) .EQ. 0 ) THEN
    
            neighbor_rank(5) = -2 ! bottom left
            neighbor_rank(7) = neighbor_rank(1) + 1 ! top left
            neighbor_rank(6) = -2 ! bottom right
            neighbor_rank(8) = -2 ! top right
    
        
         ELSE IF(mycoord(1) .EQ. 0  .AND. mycoord(2) .EQ. nranks_y-1 ) THEN
    
            neighbor_rank(5) = -2 ! bottom left
            neighbor_rank(7) = -2 ! top left
            neighbor_rank(6) = neighbor_rank(2) - 1 ! bottom right
            neighbor_rank(8) = -2 ! top right
   
        END IF
    
    ELSE
        neighbor_rank(5:8) = -2
    END IF
    

    !IF(myrank .EQ. 0) THEN
        PRINT*,'####################################################################'
        PRINT*,'Myrank, coords = ',myrank, mycoord
        WRITE(*,FMT='("Neighbors: left, right, bottom, top, bl, br, tl, tr = ",I2,2X,I2,2X,I2,2X,I2,2X,I2,2X,I2,2X,I2,2X,I2)') &
        neighbor_rank(1), neighbor_rank(2), neighbor_rank(3), neighbor_rank(4), &
        neighbor_rank(5), neighbor_rank(6), neighbor_rank(7), neighbor_rank(8)        
        PRINT*,'####################################################################'
    !END IF

    ! compute spatial cell indices of domain bottom and number of ghost cells on the boundaries 
    xlow = 1 + mycoord(1) * nx  
    ylow = 1 + mycoord(2) * ny
    
    ! bottom end ghost cells
    IF(mycoord(1)  .EQ. 0) THEN    
        nb_xlow = 2  ! two ghost cells 
        
        bxlow = -1
        exlow = 0
        
    ELSE 
        nb_xlow = 1 ! one ghost cell
        
        bxlow = 1
        exlow = 1         
        
    END IF    
    IF(mycoord(2)  .EQ. 0) THEN   
        nb_ylow = 2  ! two ghost cells 
        
        bylow = -1
        eylow = 0
        
        
    ELSE 
        nb_ylow = 1 ! one ghost cell
        
        bylow = 1
        eylow = 1
        
    END IF    
    
    ! top-end ghost cells
    IF(mycoord(1)  .EQ. nranks_x-1) THEN   
        nb_xhi = 3  ! three ghost cells 
        
        bxhi = nx + 2
        exhi = nx + 3
        
    ELSE 
        nb_xhi = 1 ! one ghost cell

        bxhi = nx 
        exhi = nx 
        
    END IF
    
    IF(mycoord(2)  .EQ. nranks_y-1) THEN   
        nb_yhi = 3  ! three ghost cells 
        
        byhi = ny + 2
        eyhi = ny + 3
        
    ELSE 
        nb_yhi = 1 ! one ghost cell
        
        byhi = ny
        eyhi = ny
        
    END IF    
    
    bzlow = -1
    ezlow = 0
        
    bzhi = nz+2
    ezhi = nz+3
    

    ! periodic boundaries
    IF(bndry .EQ. 1) THEN

        bxlow = 1
        bxhi = nx
        bylow = 1
        byhi = ny
        bzlow = 1
        bzhi = nz

        exlow = 1
        exhi = nx
        eylow = 1
        eyhi = ny
        ezlow = 1
        ezhi = nz
        
    END IF    

    
    PRINT*,''    
    PRINT*,'bx,y,z low = ', bxlow, bylow, bzlow    
    PRINT*,'bx,y,z high = ', bxhi, byhi, bzhi
    PRINT*,'ex,y,z low = ', exlow, eylow, ezlow    
    PRINT*,'ex,y,z hi = ', exhi, eyhi, ezhi    
    PRINT*,''    
    

END SUBROUTINE setup_domain



SUBROUTINE terminate_mpi()

    DEALLOCATE(field_buffer_xm,field_buffer_xp)
    IF(ndims .EQ. 2) DEALLOCATE(field_buffer_ym,field_buffer_yp)
   
    CALL MPI_FINALIZE(ierr)

END SUBROUTINE terminate_mpi


!########################################
! MPI data buffering and communications ! 
!########################################


! Top-level routine for magnetic Field exchange (outermost interior layers)
SUBROUTINE exchange_B(layers, offset)
      
    INTEGER, INTENT(IN) :: layers, offset ! specifies number of boundary cell layers we want to exchange and which ones

    ! clear all buffers
    field_buffer_xm = 0.d0
    field_buffer_xp = 0.d0
    IF(ndims .EQ. 2) THEN 
        field_buffer_ym = 0.d0
        field_buffer_yp = 0.d0
    END IF
    
    CALL exchange_fields(Bx, By, Bz, layers, offset, .FALSE.)
   

END SUBROUTINE exchange_B



! Top-level routine for electric field exchange (outermost interior layers)
SUBROUTINE exchange_E(layers, offset)
      
    INTEGER, INTENT(IN) :: layers, offset ! specifies number of boundary cell layers we want to exchange and which ones
    
    ! clear all buffers
    field_buffer_xm = 0.d0
    field_buffer_xp = 0.d0
    IF(ndims .EQ. 2) THEN 
        field_buffer_ym = 0.d0
        field_buffer_yp = 0.d0
    END IF
    
    CALL exchange_fields(Ex, Ey, Ez, layers, offset, .FALSE.)  
       

END SUBROUTINE exchange_E



! Top-level routine for current exchange (outermost interior layer and first boundary exterior layer)
SUBROUTINE exchange_current(layers, offset)

    INTEGER, INTENT(IN) :: layers, offset ! specifies number of boundary cell layers we want to exchange and which ones


     ! clear all buffers
    field_buffer_xm = 0.d0
    field_buffer_xp = 0.d0
    IF(ndims .EQ. 2) THEN 
        field_buffer_ym = 0.d0
        field_buffer_yp = 0.d0
    END IF
    
    CALL exchange_fields(Jx, Jy, Jz, layers, offset, .TRUE.)  
    

END SUBROUTINE exchange_current



! low-level field exchange routine with RMA puts
SUBROUTINE exchange_fields(fx, fy, fz, layers, offset, append)

    REAL(8), INTENT(INOUT) :: fx(-1:nx+3,-1:ny+3,-1:nz+3), fy(-1:nx+3,-1:ny+3,-1:nz+3), fz(-1:nx+3,-1:ny+3,-1:nz+3)
    INTEGER, INTENT(IN) :: layers ! specifies number of boundary cell layers we want to exchange
    INTEGER, INTENT(IN) :: offset ! specifies which layers we want to exchange
    LOGICAL, INTENT(IN) :: append ! specifies whether to append to or replace with incoming ghost cells in field arrays during unpack (e.g. J needs to be appended, while E/B are replaced)
    INTEGER :: ierr, i, j, k
    INTEGER(KIND=MPI_ADDRESS_KIND) :: disp(4) 
   
   
    IF(nranks_x .GT. 1) THEN
   
    ! clear my RMA window (this is important!!!)
    mem_fields = 0.d0

    ! put fence around rma put operations
    CALL MPI_WIN_FENCE(0, win, ierr)


        ! Set displacement units/offsets. 
        ! RMA window broken into 4 chunks organized as follows: 
        ! (i) The first (ny+5)*(nz+5) slots are for target x- ghost cells
        ! (ii) The next (ny+5)*(nz+5) slots are for target x+ ghost cells
        ! (iii) The next (nx+5)*(nz+5) slots are for target y- ghost cells
        ! (iii) The next (nx+5)*(nz+5) slots are for target y- ghost cells
        !
        ! Calculate the displacement units/offsets accordingly.     
        
        disp(1) = 0
        disp(2) = disp(1) + layers*field_buffer_size_x/2
        disp(3) = disp(2) + layers*field_buffer_size_x/2
        disp(4) = disp(3) + layers*field_buffer_size_y/2


        ! put my x- interior layer into x+ ghost cells of left neighbor
        
        IF(neighbor_rank(1) .NE. -2) THEN
        
            CALL field_packup(field_buffer_xm, 1, fx, fy, fz, layers, offset)  
            CALL MPI_PUT(field_buffer_xm, layers*field_buffer_size_x/2, MPI_DOUBLE_PRECISION, neighbor_rank(1) , &
                     disp(2), layers*field_buffer_size_x/2, MPI_DOUBLE_PRECISION, win, ierr)  

        END IF
        
        ! put my x+ interior layer into x- ghost cells of right neighbor
        IF(neighbor_rank(2) .NE. -2) THEN
        
            CALL field_packup(field_buffer_xp, 2, fx, fy, fz, layers, offset)        
            CALL MPI_PUT(field_buffer_xp, layers*field_buffer_size_x/2, MPI_DOUBLE_PRECISION, neighbor_rank(2) , &
                     disp(1), layers*field_buffer_size_x/2, MPI_DOUBLE_PRECISION, win, ierr)  
        
        END IF

        IF(ndims .EQ. 2) THEN

            ! put my y- interior layer into y+ ghost cells of bottom neighbor
            IF(neighbor_rank(3) .NE. -2) THEN

                CALL field_packup(field_buffer_ym, 3, fx, fy, fz, layers, offset)
                CALL MPI_PUT(field_buffer_ym, layers*field_buffer_size_y/2, MPI_DOUBLE_PRECISION, neighbor_rank(3) , &
                            disp(4), layers*field_buffer_size_y/2, MPI_DOUBLE_PRECISION, win, ierr)  

            END IF
            
            ! put my y+ interior layer into y- ghost cells of top neighbor
            IF(neighbor_rank(4) .NE. -2) THEN

                CALL field_packup(field_buffer_yp, 4, fx, fy, fz, layers, offset)
                CALL MPI_PUT(field_buffer_yp, layers*field_buffer_size_y/2, MPI_DOUBLE_PRECISION, neighbor_rank(4) , &
                            disp(3), layers*field_buffer_size_y/2, MPI_DOUBLE_PRECISION, win, ierr)  

            END IF
            
        END IF

    CALL MPI_WIN_FENCE(0, win, ierr)

    !PRINT*,'RMA puts completed..'

    ! retrieve data from my RMA window
    CALL field_unpack(fx, fy, fz, layers, offset, append)    


    END IF
    

END SUBROUTINE exchange_fields



! load boundary field (B, E, J) data onto buffer
! Note: offset = 0 means last interior layer, offset = 1 means first exterior layer (cells 
! starting from offset and moving outwards into the exterior will be packed up. E.g. If offset = 0 and layers =2, then for the x- boundary
! the routine will packup the layers x=1 and x=0, in that order. Similarly, if offset = 1 and layers =2, then for the x- boundary
! the routine will packup the layers x=0 and x=-1)

SUBROUTINE field_packup(buffer_out, bndry, fx, fy, fz, layers, offset)

    REAL(8), INTENT(INOUT) :: buffer_out(:)
    INTEGER, INTENT(IN) :: bndry
    REAL(8), INTENT(IN) :: fx(-1:nx+3,-1:ny+3,-1:nz+3), fy(-1:nx+3,-1:ny+3,-1:nz+3), fz(-1:nx+3,-1:ny+3,-1:nz+3)
    INTEGER, INTENT(IN) :: layers, offset  
    INTEGER :: i, j, k, ix


  
    SELECT CASE(bndry) 
    
        !x-
        CASE(1) 
                                      
             ix = 1
             DO k = -1, nz+3
                DO j = -1, ny+3     
                   DO i = 1 - offset , 1 - offset  - (layers - 1) , -1                   
                
                        buffer_out(ix)   = fx(i,j,k)   
                        buffer_out(ix+1) = fy(i,j,k)   
                        buffer_out(ix+2) = fz(i,j,k)
                        ix = ix + 3
                        
                    END DO
                END DO
            END DO    
           
        !x+
        CASE(2) 
     
            ix = 1
            DO k = -1, nz+3
                DO j = -1, ny+3         
                    DO i = offset , offset + layers - 1                    
            
                        buffer_out(ix)   = fx(nx+i,j,k)   
                        buffer_out(ix+1) = fy(nx+i,j,k)   
                        buffer_out(ix+2) = fz(nx+i,j,k)   
                        ix = ix + 3
    
                    END DO
                END DO
            END DO    
    
        !y-
        CASE(3) 
     
             ix = 1
             DO k = -1, nz+3
                DO j = 1 - offset , 1 - offset  - (layers - 1), -1    
                    DO i = -1, nx+3         
                  
                        buffer_out(ix)   = fx(i,j,k)   
                        buffer_out(ix+1) = fy(i,j,k)   
                        buffer_out(ix+2) = fz(i,j,k)   
                        ix = ix + 3
                    
                        END DO
                    END DO
                END DO    
    
        ! y+   
        CASE(4) 
       
             ix = 1
             DO k = -1, nz+3
                DO j =  offset , offset + layers - 1    
                    DO i = -1, nx+3         
                    
                        buffer_out(ix)   = fx(i,ny+j,k)   
                        buffer_out(ix+1) = fy(i,ny+j,k)   
                        buffer_out(ix+2) = fz(i,ny+j,k)   
                        ix = ix + 3
                    
                    END DO
                END DO
            END DO    
      
    END SELECT
    

END SUBROUTINE field_packup



! Retrieve boundary field data put in my RMA window by my neighbors (this routine essentially does the inverse of the field_packup routine)
SUBROUTINE field_unpack(fx, fy, fz, layers, offset, append)

    REAL(8), INTENT(INOUT) :: fx(-1:nx+3,-1:ny+3,-1:nz+3), fy(-1:nx+3,-1:ny+3,-1:nz+3), fz(-1:nx+3,-1:ny+3,-1:nz+3)
    INTEGER, INTENT(IN) :: layers, offset
    LOGICAL, INTENT(IN) :: append
    INTEGER :: i, j, k, ix
   
   
    ! For E and B field arrays, make sure that the boundary ghost cells have been cleared out before this subroutine is called.
    ! (On the other hand, currents need to be appended to the existing values for the pre-filter exchange. For post-filter exchanges, we do not appended
    ! but replace the existing values. )


    !x-  

    IF(neighbor_rank(1) .NE. -2 .AND. .NOT. append) THEN

    ix = 1    
    DO k = -1, nz+3
        DO j = -1, ny+3         
            DO i = offset, offset + (layers-1)                
                       
               fx(i,j,k) = 0.d0 
               fy(i,j,k) = 0.d0    
               fz(i,j,k) = 0.d0      
               ix = ix + 3
                        
            END DO
        END DO
    END DO    

    END IF 
     
    !x+
    
    IF(neighbor_rank(2) .NE. -2 .AND. .NOT. append) THEN

    ix = 1 + 3*(nz+5)*(ny+5)*layers 
    DO k = -1, nz+3
        DO j = -1, ny+3         
            DO i = offset, offset + (layers-1)                    
                        
                fx(nx+1-i,j,k) = 0.d0       
                fy(nx+1-i,j,k) = 0.d0    
                fz(nx+1-i,j,k) = 0.d0     
                ix = ix + 3
                    
            END DO
        END DO
    END DO   
 
    END IF

    IF(ndims .EQ. 2) THEN    

        !y-
        
        IF(neighbor_rank(3) .NE. -2 .AND. .NOT. append) THEN

        ix = 1 + 2*3*(nz+5)*(ny+5)*layers 
        DO k = -1, nz+3
            DO j = offset, offset + (layers-1)
                DO i = -1, nx+3         
                    
                    fx(i,j,k) = 0.d0    
                    fy(i,j,k) = 0.d0   
                    fz(i,j,k) = 0.d0    
                    ix = ix + 3
                    
                END DO
            END DO
        END DO           

        END IF
        
        ! y+

        IF(neighbor_rank(4) .NE. -2 .AND. .NOT. append) THEN
        
        ix = 1 + 2*3*(nz+5)*(ny+5)*layers + 3*(nz+5)*(nx+5)*layers  
        DO k = -1, nz+3
            DO j = offset, offset + (layers-1) 
                DO i = -1, nx+3         
                    
                    fx(i,ny+1-j,k) = 0.d0    
                    fy(i,ny+1-j,k) = 0.d0    
                    fz(i,ny+1-j,k) = 0.d0    
                    ix = ix + 3
                    
                END DO
            END DO
        END DO  
               
        END IF
        
    END IF
    
    
    ! now collect boundary data from RMA window and store them on the grid arrays

    !x-                   
    IF(neighbor_rank(1) .NE. -2) THEN

    ix = 1    
    DO k = -1, nz+3
        DO j = -1, ny+3         
            DO i = offset, offset + (layers-1)                      
                       
               fx(i,j,k) = fx(i,j,k) + mem_fields(ix)       
               fy(i,j,k) = fy(i,j,k) + mem_fields(ix+1)    
               fz(i,j,k) = fz(i,j,k) + mem_fields(ix+2)     
               ix = ix + 3
                        
            END DO
        END DO
    END DO    

    END IF 
    
    !x+
    IF(neighbor_rank(2) .NE. -2) THEN

    ix = 1 + 3*(nz+5)*(ny+5)*layers 
    DO k = -1, nz+3
        DO j = -1, ny+3         
            DO i = offset, offset + (layers-1)                      
                        
                fx(nx+1-i,j,k) = fx(nx+1-i,j,k) + mem_fields(ix)       
                fy(nx+1-i,j,k) = fy(nx+1-i,j,k) + mem_fields(ix+1)    
                fz(nx+1-i,j,k) = fz(nx+1-i,j,k) + mem_fields(ix+2)     
                ix = ix + 3
                    
            END DO
        END DO
    END DO    
    
    END IF

    IF(ndims .EQ. 2) THEN

        !y-
        IF(neighbor_rank(3) .NE. -2) THEN

        ix = 1 + 2*3*(nz+5)*(ny+5)*layers 
        DO k = -1, nz+3
            DO j = offset, offset + (layers-1)   
                DO i = -1, nx+3         
                    
                    fx(i,j,k) = fx(i,j,k) + mem_fields(ix)    
                    fy(i,j,k) = fy(i,j,k) + mem_fields(ix+1)   
                    fz(i,j,k) = fz(i,j,k) + mem_fields(ix+2)    
                    ix = ix + 3
                    
                END DO
            END DO
        END DO    
        
        END IF
    
        ! y+ 
        IF(neighbor_rank(4) .NE. -2) THEN
        
        ix = 1 + 2*3*(nz+5)*(ny+5)*layers + 3*(nz+5)*(nx+5)*layers  
        DO k = -1, nz+3
            DO j = offset, offset + (layers-1)    
                DO i = -1, nx+3        
                    
                    fx(i,ny+1-j,k) = fx(i,ny+1+j,k) + mem_fields(ix)    
                    fy(i,ny+1-j,k) = fy(i,ny+1+j,k) + mem_fields(ix+1)    
                    fz(i,ny+1-j,k) = fz(i,ny+1+j,k) + mem_fields(ix+2)    
                    ix = ix + 3
                    
                END DO
            END DO
        END DO    

        END IF
        
    END IF


END SUBROUTINE field_unpack



! Subroutine for particle data exchange across domain faces and edges 
! send and recvs may need to be subcycled due to variable number of boundary particle 
SUBROUTINE exchange_particles()              
    
    INTEGER :: bndry, tag, bsize(8), barrier_request
	INTEGER :: requests1d(2), requests2d(8), status(MPI_STATUS_SIZE), ierr, msize
    LOGICAL :: barrier_done, barrier_active, flag
    
    
    
    IF(nranks_x .GT. 1) THEN
    
    tag = 1
    
    ! clear in-buffer
    particle_buffer_in = 0.d0        
       
    ! particles already packed and ready to go  (first slot contains particle number in buffer)           
    bparticles_xm(0) = Np_esc_xm 
    bparticles_xp(0) = Np_esc_xp
    IF(ndims .EQ. 2) THEN    
        bparticles_ym(0) = Np_esc_ym
        bparticles_yp(0) = Np_esc_yp
        bparticles_bl(0) = Np_esc_bl 
        bparticles_br(0) = Np_esc_br
        bparticles_tl(0) = Np_esc_tl
        bparticles_tr(0) = Np_esc_tr
    END IF
    
    bsize(1) = 1 + Np_esc_xm * nvars_particles
    bsize(2) = 1 + Np_esc_xp * nvars_particles
    IF(ndims .EQ. 2) THEN    
        bsize(3) = 1 + Np_esc_ym * nvars_particles
        bsize(4) = 1 + Np_esc_yp * nvars_particles
        bsize(5) = 1 + Np_esc_bl * nvars_particles
        bsize(6) = 1 + Np_esc_br * nvars_particles
        bsize(7) = 1 + Np_esc_tl * nvars_particles
        bsize(8) = 1 + Np_esc_tr * nvars_particles
    END IF
    
    ! non-blocking synchronous sends to all neighbors     
    
    IF (ndims .EQ. 1) THEN  
        ! send across faces
        CALL MPI_ISSEND(bparticles_xm(0),        &  ! send buffer
                   bsize(1),                 &  ! send count
                   MPI_DOUBLE_PRECISION,    &  ! data type
                   neighbor_rank(1),        &  ! dest
                   tag,                     &  ! send tag
                   comm2d,                  &
                   requests1d(1),             &
                   ierr)
    
        CALL MPI_ISEND(bparticles_xp(0),        &  ! send buffer
                   bsize(2),                   &  ! send count
                   MPI_DOUBLE_PRECISION,     &  ! data type
                   neighbor_rank(2),        &  ! dest
                   tag,                     &  ! send tag
                   comm2d,                  &
                   requests1d(2),             &
                   ierr)

     
    ELSE IF(ndims .EQ. 2) THEN    
    
        ! send across faces
        CALL MPI_ISSEND(bparticles_xm(0),        &  ! send buffer
                   bsize(1),                 &  ! send count
                   MPI_DOUBLE_PRECISION,    &  ! data type
                   neighbor_rank(1),        &  ! dest
                   tag,                     &  ! send tag
                   comm2d,                  &
                   requests2d(1),             &
                   ierr)
    
        CALL MPI_ISEND(bparticles_xp(0),        &  ! send buffer
                   bsize(2),                   &  ! send count
                   MPI_DOUBLE_PRECISION,    &  ! data type
                   neighbor_rank(2),        &  ! dest
                   tag,                     &  ! send tag
                   comm2d,                  &
                   requests2d(2),             &
                   ierr)
    
        CALL MPI_ISEND(bparticles_ym(0),        &  ! send buffer
                    bsize(3),                &  ! send count
                   MPI_DOUBLE_PRECISION,    &  ! data type
                   neighbor_rank(3),        &  ! dest
                   tag,                     &  ! send tag
                   comm2d,                  &
                   requests2d(3),             &
                   ierr)
         
    
        CALL MPI_ISEND(bparticles_yp(0),        &  ! send buffer
                   bsize(4),                   &  ! send count
                   MPI_DOUBLE_PRECISION,    &  ! data type
                   neighbor_rank(4),        &  ! dest
                   tag,                     &  ! send tag
                   comm2d,                  &
                   requests2d(4),             &
                   ierr)
         
        ! send across edges
        CALL MPI_ISSEND(bparticles_bl(0),        &  ! send buffer
                   bsize(5),                 &  ! send count
                   MPI_DOUBLE_PRECISION,    &  ! data type
                   neighbor_rank(5),        &  ! dest
                   tag,                     &  ! send tag
                   comm2d,                  &
                   requests2d(5),             &
                   ierr)
    
        CALL MPI_ISEND(bparticles_br(0),        &  ! send buffer
                   bsize(6),                   &  ! send count
                   MPI_DOUBLE_PRECISION,    &  ! data type
                   neighbor_rank(6),        &  ! dest
                   tag,                     &  ! send tag
                   comm2d,                  &
                   requests2d(6),             &
                   ierr)
         
    
        CALL MPI_ISEND(bparticles_tl(0),        &  ! send buffer
                   bsize(7),                &  ! send count
                   MPI_DOUBLE_PRECISION,    &  ! data type
                   neighbor_rank(7),        &  ! dest
                   tag,                     &  ! send tag
                   comm2d,                  &
                   requests2d(7),             &
                   ierr)
         
    
        CALL MPI_ISEND(bparticles_tr(0),        &  ! send buffer
                   bsize(8),                &  ! send count
                   MPI_DOUBLE_PRECISION,    &  ! data type
                   neighbor_rank(8),        &  ! dest
                   tag,                     &  ! send tag
                   comm2d,                  &
                   requests2d(8),             &
                   ierr)
         
    END IF
    
    
    barrier_done = .FALSE.
    barrier_active = .FALSE.
    
    ! start a non-blocking barrier and check for completion of sends and issue
    ! receives accordingly
    DO WHILE(.NOT. barrier_done)       
    
        ! probe for send completions
        CALL MPI_IPROBE(MPI_ANY_SOURCE, tag, comm2d, flag, status, ierr)

        ! get message size from IPROBE status
        !CALL MPI_GET_COUNT(status, MPI_DOUBLE_PRECISION, msize, ierr)
 
        ! if a send is complete, receive a message
        IF(flag) THEN
        
            
            msize = particle_buffer_size

            CALL MPI_RECV(particle_buffer_in,  &  ! recv buffer
                   msize,                      &  ! recv count
                   MPI_DOUBLE_PRECISION,       &  ! data type
                   MPI_ANY_SOURCE,             &  ! source
                   tag,                        &  ! recv tag
                   comm2d,                     &
                   status, ierr)
  
  
            ! unpack incoming  data
            !PRINT*,'Unpacking in-buffer...'
            CALL unpack_particle_data(particle_buffer_in)             
  
  
        END IF
    
    
        IF(.NOT. barrier_active) THEN  
        
            ! start a non-blocking barrier
            IF(ndims .EQ. 1) THEN
            
                CALL MPI_TESTALL(2, requests1d, flag, MPI_STATUSES_IGNORE, ierr)
                IF(flag) THEN
                    CALL MPI_IBARRIER(comm2d, barrier_request, ierr) ! non-blocking barrier
                    barrier_active = .TRUE.                
                END IF   
                
            ELSE IF(ndims .EQ. 2) THEN    

                CALL MPI_TESTALL(8, requests2d, flag, MPI_STATUSES_IGNORE, ierr)
                IF(flag) THEN
                    CALL MPI_IBARRIER(comm2d, barrier_request, ierr) ! non-blocking barrier
                    barrier_active = .TRUE.                
                END IF     
                
            END IF
        
        ELSE
         
            ! test for completion of sends
            CALL MPI_TEST(barrier_request, barrier_done, MPI_STATUS_IGNORE, ierr)

        END IF
    
        
    END DO

    END IF        
         
    !PRINT*,'Particle Exchange complete...My_rank = ',myrank


END SUBROUTINE exchange_particles




SUBROUTINE unpack_particle_data(buffer_in)

    REAL*8, INTENT(INOUT) :: buffer_in(:)
    INTEGER :: i, ix, Nin
    REAL*8 :: x, y, z, ux, uy, uz
    INTEGER :: sp
    
                 
                 
    ! find out how many particles are in incoming buffer             
    Nin = buffer_in(1)
  
    ! take particles from buffer and insert into main particle array
    ix = 2
    DO i = 1, Nin
                   
        sp = buffer_in(ix)             
        x  = buffer_in(ix+1) 
        y  = buffer_in(ix+2) 
        z  = buffer_in(ix+3) 
        ux = buffer_in(ix+4) 
        uy = buffer_in(ix+5) 
        uz = buffer_in(ix+6) 
                
        CALL create_particle(sp,x,y,z,ux,uy,uz)

        ix = ix + 7
                
    END DO
     
    ! reset particle count in buffer (otherwise the particles might get re-buffered during other boundary exchanges)
    buffer_in(1) = 0.d0


END SUBROUTINE unpack_particle_data


END MODULE domain_mod
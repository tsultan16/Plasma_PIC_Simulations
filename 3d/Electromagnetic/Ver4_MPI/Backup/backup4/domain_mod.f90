MODULE domain_mod

USE constants_mod
USE data_mod
USE particleMover_mod
USE MPI

IMPLICIT NONE



CONTAINS

!################
! General setup ! 
!################

! top-level MPI initialization subroutine
SUBROUTINE initialize_mpi()
    
    ! initialize
    CALL MPI_INIT(ierr)

    ! Get total number of processes
    CALL MPI_COMM_SIZE(MPI_COMM_WORLD, numprocs(1), ierr)

    ! set up domain decomposition
    CALL setup_domain()

    ! get my rank 
    CALL MPI_COMM_RANK(comm2d, myrank, ierr)

    ! get my rank co-ordinates
    CALL MPI_CART_COORDS(comm2d, myrank, ndims, mycoord, ierr)

    ! get neighbor ranks
    CALL MPI_CART_SHIFT(comm2d, 0, 1, neighbor_rank(1), neighbor_rank(2), ierr) 
    CALL MPI_CART_SHIFT(comm2d, 1, 1, neighbor_rank(3), neighbor_rank(4), ierr)

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
        
    PRINT*,''    
    PRINT*,'bx,y,z low = ', bxlow, bylow, bzlow    
    PRINT*,'bx,y,z high = ', bxhi, byhi, bzhi
    PRINT*,'ex,y,z low = ', exlow, eylow, ezlow    
    PRINT*,'ex,y,z hi = ', exhi, eyhi, ezhi    
    PRINT*,''    

    !*******************************
    !allocate memory for MPI buffers
    !*******************************
    
    ! Compute buffer size needed for field data, one-layer deep (choose from max of x,y boundary sizes).
    ! and particle data.
    
    field_buffer_size = nvars_fields * MAX((nx+5)*(nz+5), (ny+5)*(nz+5)) 
    ALLOCATE(field_buffer_in(1:field_buffer_size),field_buffer_out(1:field_buffer_size))
    !ALLOCATE(field_buffer_in_2(1:field_buffer_size),field_buffer_out_2(1:field_buffer_size))
    !ALLOCATE(field_buffer_in_3(1:field_buffer_size),field_buffer_out_3(1:field_buffer_size))
    !ALLOCATE(field_buffer_in_4(1:field_buffer_size),field_buffer_out_4(1:field_buffer_size))
    
    particle_buffer_size = 1 + nvars_particles*max_particles    
    ALLOCATE(particle_buffer_in(1:particle_buffer_size))
   

PRINT*,''
PRINT*,'My rank, coordinate =',myrank, mycoord
PRINT*,''


END SUBROUTINE initialize_mpi



! domain decomposition subroutine
SUBROUTINE setup_domain()

    ! set up cart communicator properties
	isperiodic(1) = .FALSE. ! periodic boundary setting in 1st dimension
	isperiodic(2) = .FALSE. ! periodic boundary setting in 2nd dimension
	reorder = .TRUE. ! allow automatic rank reordering 
	dims = (/ nranks_x, nranks_y /) 
	
    ! create cart communicator handle
	CALL MPI_CART_CREATE(MPI_COMM_WORLD, ndims, dims, isperiodic, reorder, comm2d, ierr)


END SUBROUTINE setup_domain



SUBROUTINE terminate_mpi()

    DEALLOCATE(field_buffer_in, field_buffer_out)
    DEALLOCATE(particle_buffer_in)

    CALL MPI_FINALIZE(ierr)

END SUBROUTINE terminate_mpi


!########################################
! MPI data buffering and communications ! 
!########################################


! Magnetic Field exchange
SUBROUTINE exchange_B()

    INTEGER :: bndry, s_tag(4), r_tag(4), bsize
	INTEGER :: requests(2), status(MPI_STATUS_SIZE)
    
    
    s_tag = (/2, 1, 4, 3/) 
    r_tag = (/1, 2, 3, 4/)
    
    ! Loop over boundaries
    DO bndry = 1, 4
 
        ! clear buffers
        field_buffer_out = 0.d0
        field_buffer_in = 0.d0
 
        !PRINT*,'My_rank = ',myrank
        !PRINT*,'bndry, neighbor_rank = ',bndry, neighbor_rank(bndry)
        
        ! pack outgoing field data
        !PRINT*,'Packing out buffer...'
        CALL packup_B(field_buffer_out, bndry, bsize)        
                
        ! mpi sendrecv
        !PRINT*,'Begin communication...'

        CALL MPI_ISEND(field_buffer_out,        &  ! send buffer
                       bsize,                   &  ! send count
                       MPI_DOUBLE_PRECISION,    &  ! data type
                       neighbor_rank(bndry),    &  ! dest
                       s_tag(bndry),            &  ! send tag
                       comm2d,                  &
                       requests(1),             &
                       ierr)
         
        CALL MPI_IRECV(field_buffer_in,         &  ! send buffer
                       bsize,       &  ! send count
                       MPI_DOUBLE_PRECISION,    &  ! data type
                       neighbor_rank(bndry),    &  ! source
                       r_tag(bndry),            &  ! recv tag
                       comm2d,                  &
                       requests(2),             &
                       ierr)
                              
                
        ! now wait until data exchange has completed
        CALL MPI_WAITALL(2, requests, MPI_STATUSES_IGNORE, ierr)    

        ! unpack incoming field data
        !PRINT*,'Unpacking in-buffer...'
        CALL unpack_B(field_buffer_in, bndry)               
  
  
  END DO
    

END SUBROUTINE exchange_B



! Electric field exchange
SUBROUTINE exchange_E()

    INTEGER :: bndry, s_tag(4), r_tag(4), bsize
	INTEGER :: requests(2), status(MPI_STATUS_SIZE)
    
    
    s_tag = (/2, 1, 4, 3/) 
    r_tag = (/1, 2, 3, 4/)
    
    
    ! Loop over boundaries
    DO bndry = 1, 4
 
        ! clear buffers
        field_buffer_out = 0.d0
        field_buffer_in = 0.d0
 
        !PRINT*,'My_rank = ',myrank
        !PRINT*,'bndry, neighbor_rank = ',bndry, neighbor_rank(bndry)
        
        ! pack outgoing field data
        !PRINT*,'Packing out buffer...'
        CALL packup_E(field_buffer_out, bndry, bsize)        
                
        ! mpi sendrecv
        !PRINT*,'Begin communication...'
      
        CALL MPI_ISEND(field_buffer_out,        &  ! send buffer
                       bsize,                   &  ! send count
                       MPI_DOUBLE_PRECISION,    &  ! data type
                       neighbor_rank(bndry),    &  ! dest
                       s_tag(bndry),                     &  ! send tag
                       comm2d,                  &
                       requests(1),             &
                       ierr)
         
        CALL MPI_IRECV(field_buffer_in,         &  ! send buffer
                       bsize,                   &  ! send count
                       MPI_DOUBLE_PRECISION,    &  ! data type
                       neighbor_rank(bndry),    &  ! source
                       r_tag(bndry),                     &  ! recv tag
                       comm2d,                  &
                       requests(2),             &
                       ierr)
                              
                
        ! now wait until data exchange has completed
        CALL MPI_WAITALL(2, requests, MPI_STATUSES_IGNORE, ierr)    

        ! unpack incoming field data
        !PRINT*,'Unpacking in-buffer...'
        CALL unpack_E(field_buffer_in, bndry)               
  
  
  END DO
    

END SUBROUTINE exchange_E


! Electric field exchange
SUBROUTINE exchange_current()

    INTEGER :: bndry, s_tag(4), r_tag(4), bsize
	INTEGER :: requests(2), status(MPI_STATUS_SIZE)
    REAL*8 :: maxj
    
    
    s_tag = (/2, 1, 4, 3/) 
    r_tag = (/1, 2, 3, 4/)
    
    
    ! Loop over boundaries
    DO bndry = 1, 4
 
        ! clear buffers
        field_buffer_out = 0.d0
        field_buffer_in = 0.d0
 
        !PRINT*,'My_rank = ',myrank
        !PRINT*,'bndry, neighbor_rank = ',bndry, neighbor_rank(bndry)
        
        ! pack outgoing field data
        !PRINT*,'Packing out buffer...'
        CALL packup_J(field_buffer_out, bndry, bsize)        
                
        ! mpi sendrecv
        !PRINT*,'Begin communication...'
      
        CALL MPI_ISEND(field_buffer_out,        &  ! send buffer
                       bsize,                   &  ! send count
                       MPI_DOUBLE_PRECISION,    &  ! data type
                       neighbor_rank(bndry),    &  ! dest
                       s_tag(bndry),                     &  ! send tag
                       comm2d,                  &
                       requests(1),             &
                       ierr)
         
        CALL MPI_IRECV(field_buffer_in,         &  ! send buffer
                       bsize,                   &  ! send count
                       MPI_DOUBLE_PRECISION,    &  ! data type
                       neighbor_rank(bndry),    &  ! source
                       r_tag(bndry),                     &  ! recv tag
                       comm2d,                  &
                       requests(2),             &
                       ierr)
                              
                
        ! now wait until data exchange has completed
        CALL MPI_WAITALL(2, requests, MPI_STATUSES_IGNORE, ierr)    

        ! unpack incoming field data
        !PRINT*,'Unpacking in-buffer...'
        CALL unpack_J(field_buffer_in, bndry)               
  
  
  END DO
    

END SUBROUTINE exchange_current


! Magnetic field data pack up into outgoing buffer
SUBROUTINE packup_B(buffer_out, bndry, bsize)

    REAL*8, INTENT(INOUT) :: buffer_out(:)
    INTEGER, INTENT(IN) :: bndry
    INTEGER, INTENT(INOUT) :: bsize
    INTEGER :: i, j, k, ix


    SELECT CASE(bndry) 
    
        !x-
        CASE(1) 
                
             IF(nb_xlow .EQ. 1) THEN   
     
                 ix = 1
                 DO k = -1, nz+3
                    DO j = 0, ny+1         
                    
                        buffer_out(ix)   = Bx(1,j,k)   
                        buffer_out(ix+1) = By(1,j,k)   
                        buffer_out(ix+2) = Bz(1,j,k)
                        ix = ix + 3
                    
                    END DO
                END DO    
    
             END IF
             
             bsize = 3 * (ny+2) * (nz+5)
     
        !x+
        CASE(2) 

             IF(nb_xhi .EQ. 1) THEN   
     
                 ix = 1
                 DO k = -1, nz+3
                    DO j = 0, ny+1         
                    
                        buffer_out(ix)   = Bx(nx,j,k)   
                        buffer_out(ix+1) = By(nx,j,k)   
                        buffer_out(ix+2) = Bz(nx,j,k)   
                        ix = ix + 3
                    
                    END DO
                END DO    
    
             END IF

            bsize = 3 * (ny+2) * (nz+5)

        !y-
        CASE(3) 

                IF(nb_ylow .EQ. 1) THEN   
     
                 ix = 1
                 DO k = -1, nz+3
                    DO i = 0, nx+1         
                    
                        buffer_out(ix)   = Bx(i,1,k)   
                        buffer_out(ix+1) = By(i,1,k)   
                        buffer_out(ix+2) = Bz(i,1,k)   
                        ix = ix + 3
                    
                    END DO
                END DO    
    
             END IF

            bsize = 3 * (nx+2) * (nz+5)

        ! y+   
        CASE(4) 
  
            IF(nb_yhi .EQ. 1) THEN   
     
                 ix = 1
                 DO k = -1, nz+3
                    DO i = 0, nx+1         
                    
                        buffer_out(ix)   = Bx(i,ny,k)   
                        buffer_out(ix+1) = By(i,ny,k)   
                        buffer_out(ix+2) = Bz(i,ny,k)   
                        ix = ix + 3
                    
                    END DO
                END DO    
    
             END IF

             bsize = 3 * (nx+2) * (nz+5)
  
    END SELECT


END SUBROUTINE packup_B


! Electric field data pack up into outgoing buffer
SUBROUTINE packup_E(buffer_out, bndry, bsize)

    REAL*8, INTENT(INOUT) :: buffer_out(:)
    INTEGER, INTENT(IN) :: bndry
    INTEGER, INTENT(INOUT) :: bsize
    INTEGER :: i, j, k, ix


    SELECT CASE(bndry) 
    
        !x-
        CASE(1) 
                
             IF(nb_xlow .EQ. 1) THEN   
     
                 ix = 1
                 DO k = -1, nz+3
                    DO j = -1, ny+3         
                    
                        buffer_out(ix)   = Ex(1,j,k)   
                        buffer_out(ix+1) = Ey(1,j,k)   
                        buffer_out(ix+2) = Ez(1,j,k)
                        ix = ix + 3
                    
                    END DO
                END DO    
    
             END IF
             
             bsize = 3 * (ny+5) * (nz+5)
     
        !x+
        CASE(2) 

             IF(nb_xhi .EQ. 1) THEN   
     
                 ix = 1
                 DO k = -1, nz+3
                    DO j = -1, ny+3         
                    
                        buffer_out(ix)   = Ex(nx,j,k)   
                        buffer_out(ix+1) = Ey(nx,j,k)   
                        buffer_out(ix+2) = Ez(nx,j,k)   
                        ix = ix + 3
                    
                    END DO
                END DO    
    
             END IF

            bsize = 3 * (ny+5) * (nz+5)

        !y-
        CASE(3) 

                IF(nb_ylow .EQ. 1) THEN   
     
                 ix = 1
                 DO k = -1, nz+3
                    DO i = -1, nx+3         
                    
                        buffer_out(ix)   = Ex(i,1,k)   
                        buffer_out(ix+1) = Ey(i,1,k)   
                        buffer_out(ix+2) = Ez(i,1,k)   
                        ix = ix + 3
                    
                    END DO
                END DO    
    
             END IF

            bsize = 3 * (nx+5) * (nz+5)

        ! y+   
        CASE(4) 
  
            IF(nb_yhi .EQ. 1) THEN   
     
                 ix = 1
                 DO k = -1, nz+3
                    DO i = -1, nx+3         
                    
                        buffer_out(ix)   = Ex(i,ny,k)   
                        buffer_out(ix+1) = Ey(i,ny,k)   
                        buffer_out(ix+2) = Ez(i,ny,k)   
                        ix = ix + 3
                    
                    END DO
                END DO    
    
             END IF

             bsize = 3 * (nx+5) * (nz+5)
  
    END SELECT


END SUBROUTINE packup_E



! Electric field data pack up into outgoing buffer
SUBROUTINE packup_J(buffer_out, bndry, bsize)

    REAL*8, INTENT(INOUT) :: buffer_out(:)
    INTEGER, INTENT(IN) :: bndry
    INTEGER, INTENT(INOUT) :: bsize
    INTEGER :: i, j, k, ix


    SELECT CASE(bndry) 
    
        !x-
        CASE(1) 
                
             IF(nb_xlow .EQ. 1) THEN   
     
                 ix = 1
                 DO k = 0, nz+1
                    DO j = 0, ny+1         
                    
                        buffer_out(ix)   = Jx(0,j,k)   
                        buffer_out(ix+1) = Jy(0,j,k)   
                        buffer_out(ix+2) = Jz(0,j,k)
                        ix = ix + 3
                    
                    END DO
                END DO    
    
             END IF
             
             bsize = 3 * (ny+2) * (nz+2)
     
        !x+
        CASE(2) 

             IF(nb_xhi .EQ. 1) THEN   
     
                 ix = 1
                 DO k = 0, nz+1
                    DO j = 0, ny+1         
                    
                        buffer_out(ix)   = Jx(nx+1,j,k)   
                        buffer_out(ix+1) = Jy(nx+1,j,k)   
                        buffer_out(ix+2) = Jz(nx+1,j,k)   
                        ix = ix + 3
                    
                    END DO
                END DO    
    
             END IF

            bsize = 3 * (ny+2) * (nz+2)

        !y-
        CASE(3) 

                IF(nb_ylow .EQ. 1) THEN   
     
                 ix = 1
                 DO k = 0, nz+1
                    DO i = 0, nx+1         
                    
                        buffer_out(ix)   = Jx(i,0,k)   
                        buffer_out(ix+1) = Jy(i,0,k)   
                        buffer_out(ix+2) = Jz(i,0,k)   
                        ix = ix + 3
                    
                    END DO
                END DO    
    
             END IF

            bsize = 3 * (nx+2) * (nz+2)

        ! y+   
        CASE(4) 
  
            IF(nb_yhi .EQ. 1) THEN   
     
                 ix = 1
                 DO k = 0, nz+1
                    DO i = 0, nx+1         
                    
                        buffer_out(ix)   = Jx(i,ny+1,k)   
                        buffer_out(ix+1) = Jy(i,ny+1,k)   
                        buffer_out(ix+2) = Jz(i,ny+1,k)   
                        ix = ix + 3
                    
                    END DO
                END DO    
    
             END IF

             bsize = 3 * (nx+2) * (nz+2)
  
    END SELECT


END SUBROUTINE packup_J



! retrieve magnetic field data from incoming mpi buffer
SUBROUTINE unpack_B(buffer_in, bndry)

    REAL*8, INTENT(IN) :: buffer_in(:)
    INTEGER, INTENT(IN) :: bndry
    INTEGER :: i, j, k, ix

    SELECT CASE(bndry) 

        !x-
        CASE(1) 
                
             IF(nb_xlow .EQ. 1) THEN   
     
                 ix = 1
                 DO k = -1, nz+3
                    DO j = 0, ny+1         
                    
                        Bx(0,j,k) = buffer_in(ix)       
                        By(0,j,k) = buffer_in(ix+1)    
                        Bz(0,j,k) = buffer_in(ix+2)     
                        ix = ix + 3
                    
                    END DO
                END DO    
    
             END IF
     
        !x+
        CASE(2) 

             IF(nb_xhi .EQ. 1) THEN   
     
                 ix = 1
                 DO k = -1, nz+3
                    DO j = 0, ny+1         
                    
                        Bx(nx+1,j,k) = buffer_in(ix)       
                        By(nx+1,j,k) = buffer_in(ix+1)    
                        Bz(nx+1,j,k) = buffer_in(ix+2)     
                        ix = ix + 3
                    
                    END DO
                END DO    
    
             END IF


        !y-
        CASE(3) 

                IF(nb_ylow .EQ. 1) THEN   
     
                 ix = 1
                 DO k = -1, nz+3
                    DO i = 0, nx+1         
                    
                        Bx(i,0,k) = buffer_in(ix)    
                        By(i,0,k) = buffer_in(ix+1)   
                        Bz(i,0,k) = buffer_in(ix+2)    
                        ix = ix + 3
                    
                    END DO
                END DO    
    
             END IF

        ! y+   
        CASE(4) 
  
            IF(nb_yhi .EQ. 1) THEN   
     
                 ix = 1
                 DO k = -1, nz+3
                    DO i = 0, nx+1         
                    
                        Bx(i,ny+1,k) = buffer_in(ix)    
                        By(i,ny+1,k) = buffer_in(ix+1)    
                        Bz(i,ny+1,k) = buffer_in(ix+2)    
                        ix = ix + 3
                    
                    END DO
                END DO    
    
             END IF

  
    END SELECT



END SUBROUTINE unpack_B


! retrieve magnetic field data from incoming mpi buffer
SUBROUTINE unpack_E(buffer_in, bndry)

    REAL*8, INTENT(IN) :: buffer_in(:)
    INTEGER, INTENT(IN) :: bndry
    INTEGER :: i, j, k, ix

    SELECT CASE(bndry) 

        !x-
        CASE(1) 
                
             IF(nb_xlow .EQ. 1) THEN   
     
                 ix = 1
                 DO k = -1, nz+3
                    DO j = -1, ny+3         
                    
                        Ex(0,j,k) = buffer_in(ix)       
                        Ey(0,j,k) = buffer_in(ix+1)    
                        Ez(0,j,k) = buffer_in(ix+2)     
                        ix = ix + 3
                    
                    END DO
                END DO    
    
             END IF
     
        !x+
        CASE(2) 

             IF(nb_xhi .EQ. 1) THEN   
     
                 ix = 1
                 DO k = -1, nz+3
                    DO j = -1, ny+3         
                    
                        Ex(nx+1,j,k) = buffer_in(ix)       
                        Ey(nx+1,j,k) = buffer_in(ix+1)    
                        Ez(nx+1,j,k) = buffer_in(ix+2)     
                        ix = ix + 3
                    
                    END DO
                END DO    
    
             END IF


        !y-
        CASE(3) 

                IF(nb_ylow .EQ. 1) THEN   
     
                 ix = 1
                 DO k = -1, nz+3
                    DO i = -1, nx+3         
                    
                        Ex(i,0,k) = buffer_in(ix)    
                        Ey(i,0,k) = buffer_in(ix+1)   
                        Ez(i,0,k) = buffer_in(ix+2)    
                        ix = ix + 3
                    
                    END DO
                END DO    
    
             END IF

        ! y+   
        CASE(4) 
  
            IF(nb_yhi .EQ. 1) THEN   
     
                 ix = 1
                 DO k = -1, nz+3
                    DO i = -1, nx+3         
                    
                        Ex(i,ny+1,k) = buffer_in(ix)    
                        Ey(i,ny+1,k) = buffer_in(ix+1)    
                        Ez(i,ny+1,k) = buffer_in(ix+2)    
                        ix = ix + 3
                    
                    END DO
                END DO    
    
             END IF

  
    END SELECT



END SUBROUTINE unpack_E



! retrieve magnetic field data from incoming mpi buffer
SUBROUTINE unpack_J(buffer_in, bndry)

    REAL*8, INTENT(IN) :: buffer_in(:)
    INTEGER, INTENT(IN) :: bndry
    INTEGER :: i, j, k, ix

    SELECT CASE(bndry) 

        !x-
        CASE(1) 
                
             IF(nb_xlow .EQ. 1) THEN   
     
                 ix = 1
                 DO k = 0, nz+1
                    DO j = 0, ny+1         
                    
                        Jx(1,j,k) = Jx(1,j,k) + buffer_in(ix)       
                        Jy(1,j,k) = Jy(1,j,k) + buffer_in(ix+1)    
                        Jz(1,j,k) = Jz(1,j,k) + buffer_in(ix+2)     
                        ix = ix + 3
                    
                    END DO
                END DO    
    
             END IF
     
        !x+
        CASE(2) 

             IF(nb_xhi .EQ. 1) THEN   
     
                 ix = 1
                 DO k = 0, nz+1
                    DO j = 0, ny+1         
                    
                        Jx(nx,j,k) = Jx(nx,j,k) + buffer_in(ix)       
                        Jy(nx,j,k) = Jy(nx,j,k) + buffer_in(ix+1)    
                        Jz(nx,j,k) = Jz(nx,j,k) + buffer_in(ix+2)     
                        ix = ix + 3
                    
                    END DO
                END DO    
    
             END IF


        !y-
        CASE(3) 

                IF(nb_ylow .EQ. 1) THEN   
     
                 ix = 1
                 DO k = 0, nz+1
                    DO i = 0, nx+1         
                    
                        Jx(i,1,k) = Jx(i,1,k) + buffer_in(ix)    
                        Jy(i,1,k) = Jy(i,1,k) + buffer_in(ix+1)   
                        Jz(i,1,k) = Jz(i,1,k) + buffer_in(ix+2)    
                        ix = ix + 3
                    
                    END DO
                END DO    
    
             END IF

        ! y+   
        CASE(4) 
  
            IF(nb_yhi .EQ. 1) THEN   
     
                 ix = 1
                 DO k = 0, nz+1
                    DO i = 0, nx+1         
                    
                        Jx(i,ny,k) = Jx(i,ny,k) + buffer_in(ix)    
                        Jy(i,ny,k) = Jy(i,ny,k) + buffer_in(ix+1)    
                        Jz(i,ny,k) = Jz(i,ny,k) + buffer_in(ix+2)    
                        ix = ix + 3
                    
                    END DO
                END DO    
    
             END IF

  
    END SELECT



END SUBROUTINE unpack_J




! Subroutine for particle data exchange    
! send and recvs may need to be subcycled due to variable number of boundary particle 
SUBROUTINE exchange_particles()              
    
    INTEGER :: bndry, tag, bsize
	INTEGER :: requests(2), status(MPI_STATUS_SIZE)
    
    tag = 1
    
    ! clear in-buffer
    particle_buffer_in = 0.d0        
       
    ! particles already packed and ready to go             
                
    !*************
    ! x- boundary
    !*************
    bparticles_xm(0) = Np_esc_xm ! first slot contains particle number in buffer
    bsize = 1 + Np_esc_xm * nvars_particles
    
    CALL MPI_ISEND(bparticles_xm(0),        &  ! send buffer
                   bsize,                   &  ! send count
                   MPI_DOUBLE_PRECISION,    &  ! data type
                   neighbor_rank(1),    &  ! dest
                   tag,                     &  ! send tag
                   comm2d,                  &
                   requests(1),             &
                   ierr)
         
    CALL MPI_IRECV(particle_buffer_in,      &  ! recv buffer
                   particle_buffer_size,    &  ! recv count
                   MPI_DOUBLE_PRECISION,    &  ! data type
                   neighbor_rank(1),    &  ! source
                   tag,                     &  ! recv tag
                   comm2d,                  &
                   requests(2),             &
                   ierr)
                              
                
    ! now wait until data exchange has completed
    CALL MPI_WAITALL(2, requests, MPI_STATUSES_IGNORE, ierr)    
    ! unpack incoming field data
    !PRINT*,'Unpacking in-buffer...'
    CALL unpack_particle_data(particle_buffer_in)               
  
  
  
    !*************
    ! x+ boundary
    !*************
    bparticles_xp(0) = Np_esc_xp ! first slot contains particle number in buffer
    bsize = 1 + Np_esc_xp * nvars_particles
    
    CALL MPI_ISEND(bparticles_xp(0),        &  ! send buffer
                   bsize,                   &  ! send count
                   MPI_DOUBLE_PRECISION,    &  ! data type
                   neighbor_rank(2),        &  ! dest
                   tag,                     &  ! send tag
                   comm2d,                  &
                   requests(1),             &
                   ierr)
         
    CALL MPI_IRECV(particle_buffer_in,      &  ! recv buffer
                   particle_buffer_size,    &  ! recv count
                   MPI_DOUBLE_PRECISION,    &  ! data type
                   neighbor_rank(2),    &  ! source
                   tag,                     &  ! recv tag
                   comm2d,                  &
                   requests(2),             &
                   ierr)
                              
                
    ! now wait until data exchange has completed
    CALL MPI_WAITALL(2, requests, MPI_STATUSES_IGNORE, ierr)    
     ! unpack incoming field data
    !PRINT*,'Unpacking in-buffer...'
    CALL unpack_particle_data(particle_buffer_in)               
    
    
    !*************
    ! y- boundary
    !*************
    bparticles_ym(0) = Np_esc_ym ! first slot contains particle number in buffer
    bsize = 1 + Np_esc_ym * nvars_particles
    
    CALL MPI_ISEND(bparticles_ym(0),        &  ! send buffer
                   bsize,                   &  ! send count
                   MPI_DOUBLE_PRECISION,    &  ! data type
                   neighbor_rank(3),    &  ! dest
                   tag,                     &  ! send tag
                   comm2d,                  &
                   requests(1),             &
                   ierr)
         
    CALL MPI_IRECV(particle_buffer_in,      &  ! recv buffer
                   particle_buffer_size,    &  ! recv count
                   MPI_DOUBLE_PRECISION,    &  ! data type
                   neighbor_rank(3),    &  ! source
                   tag,                     &  ! recv tag
                   comm2d,                  &
                   requests(2),             &
                   ierr)
                              
                
    ! now wait until data exchange has completed
    CALL MPI_WAITALL(2, requests, MPI_STATUSES_IGNORE, ierr)    
     ! unpack incoming field data
    !PRINT*,'Unpacking in-buffer...'
    CALL unpack_particle_data(particle_buffer_in)               
  

    !*************
    ! y+ boundary
    !*************
    bparticles_yp(0) = Np_esc_yp ! first slot contains particle number in buffer
    bsize = 1 + Np_esc_yp * nvars_particles
    
    CALL MPI_ISEND(bparticles_yp(0),        &  ! send buffer
                   bsize,                   &  ! send count
                   MPI_DOUBLE_PRECISION,    &  ! data type
                   neighbor_rank(4),    &  ! dest
                   tag,                     &  ! send tag
                   comm2d,                  &
                   requests(1),             &
                   ierr)
         
    CALL MPI_IRECV(particle_buffer_in,      &  ! recv buffer
                   particle_buffer_size,    &  ! recv count
                   MPI_DOUBLE_PRECISION,    &  ! data type
                   neighbor_rank(4),    &  ! source
                   tag,                     &  ! recv tag
                   comm2d,                  &
                   requests(2),             &
                   ierr)
                              
                
    ! now wait until data exchange has completed
    CALL MPI_WAITALL(2, requests, MPI_STATUSES_IGNORE, ierr)    
    ! unpack incoming field data
    !PRINT*,'Unpacking in-buffer...'
    CALL unpack_particle_data(particle_buffer_in)            

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
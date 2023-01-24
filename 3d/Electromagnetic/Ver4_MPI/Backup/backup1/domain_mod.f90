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
    ELSE 
        nb_xlow = 1 ! one ghost cell
    END IF    
    IF(mycoord(2)  .EQ. 0) THEN   
        nb_ylow = 2  ! two ghost cells 
    ELSE 
        nb_ylow = 1 ! one ghost cell
    END IF    
    
    ! top-end ghost cells
    IF(mycoord(1)  .EQ. nranks_x-1) THEN   
        nb_xhi = 3  ! three ghost cells 
    ELSE 
        nb_xhi = 1 ! one ghost cell
    END IF    
    IF(mycoord(2)  .EQ. nranks_y-1) THEN   
        nb_yhi = 3  ! three ghost cells 
    ELSE 
        nb_yhi = 1 ! one ghost cell
    END IF    
    

    !*******************************
    !allocate memory for MPI buffers
    !*******************************
    
    ! Compute buffer size needed for field data, one-layer deep (choose from max of x,y boundary sizes).
    ! and particle data.
    
    field_buffer_size = nvars_fields * MAX((nx+2)*(nz+2), (ny+2)*(nz+2)) 
    ALLOCATE(field_buffer_in(1:field_buffer_size),field_buffer_out(1:field_buffer_size))
    
    particle_buffer_size = 2 + nvars_particles*max_particles    
    ALLOCATE(particle_buffer_in(1:particle_buffer_size),particle_buffer_out(1:particle_buffer_size))
   

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
    DEALLOCATE(particle_buffer_in, particle_buffer_out)

    CALL MPI_FINALIZE(ierr)

END SUBROUTINE terminate_mpi


!########################################
! MPI data buffering and communications ! 
!########################################


! Subroutine for field data exchange     
SUBROUTINE boundary_exchange_fields()

    INTEGER :: bndry, tag
	INTEGER :: requests(2), send_status(MPI_STATUS_SIZE), recv_status(MPI_STATUS_SIZE)   
    LOGICAL :: send_cmplt, recv_cmplt      
    
    
    ! Loop over boundaries
    DO bndry = 1, 4
 
        !PRINT*,'My_rank = ',myrank
        !PRINT*,'bndry, neighbor_rank = ',bndry, neighbor_rank(bndry)

        send_cmplt = .FALSE.
        recv_cmplt = .FALSE.
        
        ! pack outgoing field data
        !PRINT*,'Packing out buffer...'
        CALL packup_field_data(field_buffer_out, bndry)        
                
        ! mpi sendrecv
        !PRINT*,'Begin communication...'
        tag = 1
      
        CALL MPI_ISEND(field_buffer_out,        &  ! send buffer
                       field_buffer_size,       &  ! send count
                       MPI_DOUBLE_PRECISION,    &  ! data type
                       neighbor_rank(bndry),    &  ! dest
                       tag,                     &  ! send tag
                       comm2d,                  &
                       requests(1),             &
                       ierr)
         
        CALL MPI_IRECV(field_buffer_in,         &  ! send buffer
                       field_buffer_size,       &  ! send count
                       MPI_DOUBLE_PRECISION,    &  ! data type
                       neighbor_rank(bndry),    &  ! source
                       tag,                     &  ! recv tag
                       comm2d,                  &
                       requests(2),             &
                       ierr)
                
        ! now wait until data exchange has completed
        DO WHILE(.NOT. send_cmplt .OR. .NOT. recv_cmplt) 
            CALL MPI_TEST(requests(1), send_cmplt, send_status, ierr)
            CALL MPI_TEST(requests(2), recv_cmplt, recv_status, ierr)

            !PRINT*,'Myrank, STATUS: send, recv = ',myrank, send_cmplt, recv_cmplt

        END DO

        ! unpack incoming field data
        !PRINT*,'Unpacking in buffer...'
        CALL unpack_field_data(field_buffer_in, bndry)               
  
  END DO

  PRINT*,'Field data exchange completed...'

END SUBROUTINE boundary_exchange_fields



! Subroutine for particle data exchange    
! send and recvs may need to be subcycled due to variable number of boundary particle 
SUBROUTINE boundary_exchange_particles()

    INTEGER :: bndry
    LOGICAL :: buffer_out_cmplt, buffer_in_cmplt, &
               send_cmplt, recv_cmplt
	INTEGER :: tag, requests(2), &
               send_status(MPI_STATUS_SIZE), recv_status(MPI_STATUS_SIZE), &
               send_buff_size               
    
    
    
    ! Loop over boundaries
    DO bndry = 1, 4
    
         !PRINT*,''
         !PRINT*,'My_rank = ',myrank
         !PRINT*,'bndry, neighbor_rank = ',bndry, neighbor_rank(bndry)
         !PRINT*,''

         ! clear buffer (really only need to reset particle numbers slot in header)
         particle_buffer_out(1) = 0.d0
         particle_buffer_in(1) = 0.d0
         
         buffer_out_cmplt = .FALSE.
         buffer_in_cmplt = .FALSE.

            ! subcycle if more particles than will fit in buffer or if more particles are expected to arrive
            DO WHILE(.NOT. buffer_out_cmplt .OR. .NOT. buffer_in_cmplt)
                
                send_cmplt = .FALSE.
                recv_cmplt = .FALSE.
              
              
                IF(.NOT. buffer_out_cmplt) THEN
                
                    ! pack outgoing particle data 
                    !PRINT*,'myrank, Packing out buffer...',myrank
                    CALL packup_particle_data(particle_buffer_out, bndry, buffer_out_cmplt)
                    !PRINT*,'myrank, bndry, Nout=',myrank, bndry, particle_buffer_out(1)
                    !IF(INT(particle_buffer_out(1)) .GT. 0) PRINT*,'myrank, bndry, Buffer_out contents = ',&
                    !                                      myrank,bndry, particle_buffer_out

                    ! compute send buffer size
                    send_buff_size = 3 + INT(particle_buffer_out(1))*6  ! first three slots are header info
                    
                    ! mpi send
                    tag = 0
                    CALL MPI_ISEND(particle_buffer_out,     &  ! send buffer
                                   send_buff_size,          &  ! send count
                                   MPI_DOUBLE_PRECISION,    &  ! data type
                                   neighbor_rank(bndry),    &  ! dest
                                   tag,                     &  ! send tag
                                   comm2d,                  &
                                   requests(1),             &
                                   ierr)
                ELSE
                    send_cmplt = .TRUE. 
                END IF

   
                IF(.NOT. buffer_in_cmplt) THEN

                    ! mpi recv
                    tag = 0
                    CALL MPI_IRECV(particle_buffer_in,      &  ! recv buffer
                                   particle_buffer_size,    &  ! recv max count
                                   MPI_DOUBLE_PRECISION,    &  ! data type
                                   neighbor_rank(bndry),    &  ! source
                                   tag,                     &  ! recv tag
                                   comm2d,                  &
                                   requests(2),             &
                                   ierr)
                ELSE
                    recv_cmplt = .TRUE.
                END IF
                
                
                ! wait until data exchange has completed
                DO WHILE(.NOT. send_cmplt .OR. .NOT. recv_cmplt) 
                    IF(.NOT. send_cmplt) CALL MPI_TEST(requests(1), send_cmplt, send_status, ierr)
                    IF(.NOT. recv_cmplt) CALL MPI_TEST(requests(2), recv_cmplt, recv_status, ierr)
                    
                    !PRINT*,'Myrank, COM STATUS: send, recv = ',myrank, send_cmplt, recv_cmplt

                END DO
                
                ! now unpack incoming particle data
                !PRINT*,'Unpacking in buffer...'
                !PRINT*,'Nin, Nleft=',INT(particle_buffer_in(1)),INT(particle_buffer_in(2))
                !IF(INT(particle_buffer_in(1)) .GT. 0) PRINT*,'myrank, bndry, Buffer_in contents =', myrank, &
                !                                      bndry, particle_buffer_in

                CALL unpack_particle_data(particle_buffer_in)              

                ! check to see if more particles are inbound
                IF(INT(particle_buffer_in(2)) .EQ. 0) THEN
                    buffer_in_cmplt = .TRUE.
                END IF
                
                !PRINT*,'Myrank, Buffer STATUS: out, in = ',myrank, buffer_out_cmplt, buffer_in_cmplt
                
            END DO
              
  END DO


  PRINT*,'Particle Exchange complete...My_rank = ',myrank


END SUBROUTINE boundary_exchange_particles



! Field data pack up into outgoing buffer (outflow boundary only, for now...)
SUBROUTINE packup_field_data(buffer_out, bndry)

    REAL*8, INTENT(INOUT) :: buffer_out(:)
    INTEGER, INTENT(IN) :: bndry
    INTEGER :: i, j, k, ix


    SELECT CASE(bndry) 
    
        !x-
        CASE(1) 
                
             IF(nb_xlow .EQ. 1) THEN   
     
                 ix = 1
                 DO k = 0, nz+1
                    DO j = 0, ny+1         
                    
                        ! field data ordering: Ex, Ey, Ez, Bx, By, Bz, Jx, Jy, Jz                        
                        buffer_out(ix)   = Ex(1,j,k)   
                        buffer_out(ix+1) = Ey(1,j,k)   
                        buffer_out(ix+2) = Ez(1,j,k)   
                        buffer_out(ix+3) = Bx(1,j,k)   
                        buffer_out(ix+4) = By(1,j,k)   
                        buffer_out(ix+5) = Bz(1,j,k)   
                        buffer_out(ix+6) = Jx(1,j,k)   
                        buffer_out(ix+7) = Jy(1,j,k)   
                        buffer_out(ix+8) = Jz(1,j,k)   
                        
                        ix = ix + 9
                    
                    END DO
                END DO    
    
             END IF
     
        !x+
        CASE(2) 

             IF(nb_xhi .EQ. 1) THEN   
     
                 ix = 1
                 DO k = 0, nz+1
                    DO j = 0, ny+1         
                    
                        ! field data ordering: Ex, Ey, Ez, Bx, By, Bz, Jx, Jy, Jz                        
                        buffer_out(ix)   = Ex(nx,j,k)   
                        buffer_out(ix+1) = Ey(nx,j,k)   
                        buffer_out(ix+2) = Ez(nx,j,k)   
                        buffer_out(ix+3) = Bx(nx,j,k)   
                        buffer_out(ix+4) = By(nx,j,k)   
                        buffer_out(ix+5) = Bz(nx,j,k)   
                        buffer_out(ix+6) = Jx(nx,j,k)   
                        buffer_out(ix+7) = Jy(nx,j,k)   
                        buffer_out(ix+8) = Jz(nx,j,k)   
                        
                        ix = ix + 9
                    
                    END DO
                END DO    
    
             END IF


        !y-
        CASE(3) 

                IF(nb_ylow .EQ. 1) THEN   
     
                 ix = 1
                 DO k = 0, nz+1
                    DO i = 0, nx+1         
                    
                        ! field data ordering: Ex, Ey, Ez, Bx, By, Bz, Jx, Jy, Jz                        
                        buffer_out(ix)   = Ex(i,1,k)   
                        buffer_out(ix+1) = Ey(i,1,k)   
                        buffer_out(ix+2) = Ez(i,1,k)   
                        buffer_out(ix+3) = Bx(i,1,k)   
                        buffer_out(ix+4) = By(i,1,k)   
                        buffer_out(ix+5) = Bz(i,1,k)   
                        buffer_out(ix+6) = Jx(i,1,k)   
                        buffer_out(ix+7) = Jy(i,1,k)   
                        buffer_out(ix+8) = Jz(i,1,k)   
                        
                        ix = ix + 9
                    
                    END DO
                END DO    
    
             END IF

        ! y+   
        CASE(4) 
  
            IF(nb_yhi .EQ. 1) THEN   
     
                 ix = 1
                 DO k = 0, nz+1
                    DO i = 0, nx+1         
                    
                        ! field data ordering: Ex, Ey, Ez, Bx, By, Bz, Jx, Jy, Jz                        
                        buffer_out(ix)   = Ex(i,ny,k)   
                        buffer_out(ix+1) = Ey(i,ny,k)   
                        buffer_out(ix+2) = Ez(i,ny,k)   
                        buffer_out(ix+3) = Bx(i,ny,k)   
                        buffer_out(ix+4) = By(i,ny,k)   
                        buffer_out(ix+5) = Bz(i,ny,k)   
                        buffer_out(ix+6) = Jx(i,ny,k)   
                        buffer_out(ix+7) = Jy(i,ny,k)   
                        buffer_out(ix+8) = Jz(i,ny,k)   
                        
                        ix = ix + 9
                    
                    END DO
                END DO    
    
             END IF

  
    END SELECT


END SUBROUTINE packup_field_data



! retrieve field data from incoming mpi buffer
SUBROUTINE unpack_field_data(buffer_in, bndry)

    REAL*8, INTENT(IN) :: buffer_in(:)
    INTEGER, INTENT(IN) :: bndry
    INTEGER :: i, j, k

    SELECT CASE(bndry) 

        !x-
        CASE(1) 
                
             IF(nb_xlow .EQ. 1) THEN   
     
                 ix = 1
                 DO k = 0, nz+1
                    DO j = 0, ny+1         
                    
                        ! field data ordering: Ex, Ey, Ez, Bx, By, Bz, Jx, Jy, Jz                        
                        Ex(0,j,k) = buffer_in(ix)       
                        Ey(0,j,k) = buffer_in(ix+1)    
                        Ez(0,j,k) = buffer_in(ix+2)     
                        Bx(0,j,k) = buffer_in(ix+3)   
                        By(0,j,k) = buffer_in(ix+4)   
                        Bz(0,j,k) = buffer_in(ix+5)     
                        Jx(0,j,k) = buffer_in(ix+6)    
                        Jy(0,j,k) = buffer_in(ix+7)    
                        Jz(0,j,k) = buffer_in(ix+8)   
                        
                        ix = ix + 9
                    
                    END DO
                END DO    
    
             END IF
     
        !x+
        CASE(2) 

             IF(nb_xhi .EQ. 1) THEN   
     
                 ix = 1
                 DO k = 0, nz+1
                    DO j = 0, ny+1         
                    
                        ! field data ordering: Ex, Ey, Ez, Bx, By, Bz, Jx, Jy, Jz                        
                        Ex(nx+1,j,k) = buffer_in(ix)       
                        Ey(nx+1,j,k) = buffer_in(ix+1)    
                        Ez(nx+1,j,k) = buffer_in(ix+2)     
                        Bx(nx+1,j,k) = buffer_in(ix+3)   
                        By(nx+1,j,k) = buffer_in(ix+4)   
                        Bz(nx+1,j,k) = buffer_in(ix+5)     
                        Jx(nx+1,j,k) = buffer_in(ix+6)    
                        Jy(nx+1,j,k) = buffer_in(ix+7)    
                        Jz(nx+1,j,k) = buffer_in(ix+8)  
                        
                        ix = ix + 9
                    
                    END DO
                END DO    
    
             END IF


        !y-
        CASE(3) 

                IF(nb_ylow .EQ. 1) THEN   
     
                 ix = 1
                 DO k = 0, nz+1
                    DO i = 0, nx+1         
                    
                        ! field data ordering: Ex, Ey, Ez, Bx, By, Bz, Jx, Jy, Jz                        
                        Ex(i,0,k) = buffer_in(ix)    
                        Ey(i,0,k) = buffer_in(ix+1)   
                        Ez(i,0,k) = buffer_in(ix+2)    
                        Bx(i,0,k) = buffer_in(ix+3)    
                        By(i,0,k) = buffer_in(ix+4)    
                        Bz(i,0,k) = buffer_in(ix+5)    
                        Jx(i,0,k) = buffer_in(ix+6)    
                        Jy(i,0,k) = buffer_in(ix+7)    
                        Jz(i,0,k) = buffer_in(ix+8)    
                        
                        ix = ix + 9
                    
                    END DO
                END DO    
    
             END IF

        ! y+   
        CASE(4) 
  
            IF(nb_yhi .EQ. 1) THEN   
     
                 ix = 1
                 DO k = 0, nz+1
                    DO i = 0, nx+1         
                    
                        ! field data ordering: Ex, Ey, Ez, Bx, By, Bz, Jx, Jy, Jz                        
                        Ex(i,ny+1,k) = buffer_in(ix)    
                        Ey(i,ny+1,k) = buffer_in(ix+1)    
                        Ez(i,ny+1,k) = buffer_in(ix+2)    
                        Bx(i,ny+1,k) = buffer_in(ix+3)    
                        By(i,ny+1,k) = buffer_in(ix+4)   
                        Bz(i,ny+1,k) = buffer_in(ix+5)    
                        Jx(i,ny+1,k) = buffer_in(ix+6)    
                        Jy(i,ny+1,k) = buffer_in(ix+7)    
                        Jz(i,ny+1,k) = buffer_in(ix+8)   
                        
                        ix = ix + 9
                    
                    END DO
                END DO    
    
             END IF

  
    END SELECT



END SUBROUTINE unpack_field_data



! Particle data pack up into outgoing buffer (outflow boundary only, for now...)
! Particles are loaded onto buffer in batches of 'max_particles'. The 'particles_buffered'
! flag gets set when all particles have been buffered.
SUBROUTINE packup_particle_data(buffer_out, bndry, buffer_out_cmplt)

    REAL*8, INTENT(INOUT) :: buffer_out(:)
    INTEGER, INTENT(IN) :: bndry
    LOGICAL, INTENT(INOUT) :: buffer_out_cmplt
    INTEGER :: i, j, k, ix, Nout


    SELECT CASE(bndry) 
    
        !x-
        CASE(1) 
                         
            ! find out how many particles to put in out buffer             
            Nout  = MIN(max_particles, Np_esc_xm)            
            
            ! first two slots are header info: (number of particles, number of particles leftover)
            buffer_out(1) = Nout
            buffer_out(2) = Np_esc_xm - Nout
            
            ! loop over all escaped electrons through the boundary and load them onto buffer 
            ! (starting from the back of the array)
            ix = 3
            
            DO i = Np_esc_xm, Np_esc_xm - Nout + 1, -1
               
                buffer_out(ix)   = bparticles_xm(i)%species
                buffer_out(ix+1) = bparticles_xm(i)%x + Lx
                buffer_out(ix+2) = bparticles_xm(i)%y
                buffer_out(ix+3) = bparticles_xm(i)%z
                buffer_out(ix+4) = bparticles_xm(i)%ux
                buffer_out(ix+5) = bparticles_xm(i)%uy
                buffer_out(ix+6) = bparticles_xm(i)%uz
                
                ix = ix + 7
                
            END DO
     
            Np_esc_xm = Np_esc_xm - Nout
            
            IF(Np_esc_xm .EQ. 0) THEN 
                buffer_out_cmplt = .TRUE.
            END IF
     
        !x+
        CASE(2) 
        
            ! find out how many particles to put in out buffer             
            Nout  = MIN(max_particles, Np_esc_xp)            
            
            ! first two slots are header info: (number of particles, number of particles leftover)
            buffer_out(1) = Nout
            buffer_out(2) = Np_esc_xp - Nout
            
            ! loop over all escaped electrons through the boundary and load them onto buffer 
            ! (starting from the back of the array)
            ix = 3
            
            DO i = Np_esc_xp, Np_esc_xp - Nout + 1, -1
               
                buffer_out(ix)   = bparticles_xp(i)%species
                buffer_out(ix+1) = bparticles_xp(i)%x - Lx
                buffer_out(ix+2) = bparticles_xp(i)%y
                buffer_out(ix+3) = bparticles_xp(i)%z
                buffer_out(ix+4) = bparticles_xp(i)%ux
                buffer_out(ix+5) = bparticles_xp(i)%uy
                buffer_out(ix+6) = bparticles_xp(i)%uz
                
                ix = ix + 7
                
            END DO
     
            Np_esc_xp = Np_esc_xp - Nout
            
            IF(Np_esc_xp .EQ. 0) THEN 
                buffer_out_cmplt = .TRUE.
            END IF
            
        !y-
        CASE(3) 

            ! find out how many particles to put in out buffer             
            Nout  = MIN(max_particles, Np_esc_ym)            
            
            ! first two slots are header info: (number of particles, number of particles leftover)
            buffer_out(1) = Nout
            buffer_out(2) = Np_esc_ym - Nout
            
            ! loop over all escaped electrons through the boundary and load them onto buffer 
            ! (starting from the back of the array)
            ix = 3
            
            DO i = Np_esc_ym, Np_esc_ym - Nout + 1, -1
               
                buffer_out(ix)   = bparticles_ym(i)%species
                buffer_out(ix+1) = bparticles_ym(i)%x 
                buffer_out(ix+2) = bparticles_ym(i)%y + Ly
                buffer_out(ix+3) = bparticles_ym(i)%z
                buffer_out(ix+4) = bparticles_ym(i)%ux
                buffer_out(ix+5) = bparticles_ym(i)%uy
                buffer_out(ix+6) = bparticles_ym(i)%uz
                
                ix = ix + 7
                
            END DO
     
            Np_esc_ym = Np_esc_ym - Nout
            
            IF(Np_esc_ym .EQ. 0) THEN 
                buffer_out_cmplt = .TRUE.
            END IF

        ! y+   
        CASE(4) 
  
            ! find out how many particles to put in out buffer             
            Nout  = MIN(max_particles, Np_esc_yp)            
            
            ! first two slots are header info: (number of particles, number of particles leftover)
            buffer_out(1) = Nout
            buffer_out(2) = Np_esc_yp - Nout
            
            ! loop over all escaped electrons through the boundary and load them onto buffer 
            ! (starting from the back of the array)
            ix = 3
            
            DO i = Np_esc_yp, Np_esc_yp - Nout + 1, -1
               
                buffer_out(ix)   = bparticles_yp(i)%species
                buffer_out(ix+1) = bparticles_yp(i)%x 
                buffer_out(ix+2) = bparticles_yp(i)%y - Ly
                buffer_out(ix+3) = bparticles_yp(i)%z
                buffer_out(ix+4) = bparticles_yp(i)%ux
                buffer_out(ix+5) = bparticles_yp(i)%uy
                buffer_out(ix+6) = bparticles_yp(i)%uz
                
                ix = ix + 7
                
            END DO
     
            Np_esc_yp = Np_esc_yp - Nout
            
            IF(Np_esc_yp .EQ. 0) THEN 
                buffer_out_cmplt = .TRUE.
            END IF

  
    END SELECT



END SUBROUTINE packup_particle_data



SUBROUTINE unpack_particle_data(buffer_in)

    REAL*8, INTENT(IN) :: buffer_in(:)
    INTEGER :: i, ix, Nin
    REAL*8 :: x, y, z, ux, uy, uz
    INTEGER :: sp
    
                 
                 
    ! find out how many particles are in incoming buffer             
    Nin = buffer_in(1)
      
    
    ! make sure they will fit in main partile array (... need to do this correctly...)
    IF(Nin .GT. Np_max-Np_in(1) .OR. Nin .GT. Np_max-Np_in(2)) THEN
        PRINT*,'Particle overflow error...too many particles in MPI &
                incoming buffer. Terminating program.'
        STOP
    END IF 
  
    ! take particles from buffer and insert into main particle array
    ix = 3
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
     

END SUBROUTINE unpack_particle_data


END MODULE domain_mod
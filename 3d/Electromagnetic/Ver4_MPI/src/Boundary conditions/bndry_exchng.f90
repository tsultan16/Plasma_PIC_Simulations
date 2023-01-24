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
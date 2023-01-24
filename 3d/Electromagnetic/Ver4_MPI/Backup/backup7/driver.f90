PROGRAM driver_plasma_3d

USE constants_mod
USE data_mod
USE domain_mod
USE init_mod
USE fieldSolver_mod
USE particleMover_mod
USE solar_wind_mod
USE OMP_LIB ! openmp runtime library

IMPLICIT NONE


INTEGER :: ii, i, j, k, tstart 
REAL*8 :: x1, y1, x2, y2, Jx_max = 0.d0, Jy_max = 0.d0, Jx_min = 0.d0, Jy_min = 0.d0, J_max = 0.d0, &
          E_max = 0.d0, B_max = 0.d0, t1, t2, t3, t4, t5, t6, tavg = 0.d0, tsim = 0.d0, tfile = 0.d0, &
          S_avg = 0.d0, S_exact = 0.d0


! check that there are more particles than cells
IF(.NOT. override_checks) THEN

    IF(Np_max .LT. nx*ny) THEN
        PRINT*,'Need N > nx*ny.'
        PRINT*,'Terminating program..'
        STOP
    END IF

    IF(MOD(Np_max,nx*ny) .NE. 0) THEN
        PRINT*,'Need N to be a multiple of nx*ny.'
        PRINT*,'Terminating program..'
        STOP
    END If

    IF(ABS(wc_e) .GT. 1.d-15 .AND. interpolation_type .NE. 2) THEN
        PRINT*,'For non-zero uniform magnetic field, need to use interpolation_type = 2 (momentum conserving).'
        PRINT*,'Terminating program...'
        STOP
    END IF

END IF


! set up computational domain
CALL set_domain()


! set up MPI
CALL initialize_mpi()


! allocate memory for fields
CALL create_grid_arrays()



IF(file_restart) THEN
     CALL read_restart()
     tstart = tstart + 1
ELSE

    tstart = 1

    ! initialize particle distribution and state 
    CALL init_particles()


END IF

! add sinusoidal perturbation to electrons
!CALL add_perturbations()

! compute initial charge distribution
!CALL setrho()

! intialize uniform fields (Ey and Bz)
!IF(myrank .NE. 0) 
!Bz = 0.01

!Interplanetary fields
!Bz = -1.
!Ey = -Bz*v0x



CALL add_ring_current()

CALL grid_avg_fields()

!CALL output_bin(0)


CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)

t1 = MPI_Wtime()



! enter simulation loop
DO i_sim = tstart, tstart+maxsteps
	
    IF(myrank .EQ. 0) THEN
        PRINT*,''
        WRITE(*,FMT='("Time step = ",i4,", Completion Status = ",f6.1," %")') &
            i_sim, 100.*REAL(i_sim-tstart)/REAL(maxsteps)
        PRINT*,''
	END IF
    
    t5 = MPI_Wtime()

    ! advance by one time-step 
    CALL solve() 
    
    
    !************************************
 	!PRINT*,''
	!PRINT*,''
	PRINT*,''
    PRINT*,'Myrank, Np_in =',myrank, Np_in(1),Np_in(2) 
    !DO i = 1, Np_in(1)
        !PRINT*,'Particle position: x, y = ',particles(1,1)%x,particles(1,1)%y
        !PRINT*,'Particle velocity: ux, uy = ',particles(1,1)%ux,particles(1,1)%uy
    !END DO
 	PRINT*,''
	!PRINT*,''
	!PRINT*,''

    !************************************
    
    
    ! compute poynting flux through computational boundary
    !CALL compute_poynting_flux() 
    
	! save to file
	IF(MOD(i_sim,tskip) .EQ. 0)THEN
        !PRINT*,'Saving to file...'   
        t3 = MPI_Wtime()
        CALL output_bin(i_sim)
        t4 = MPI_Wtime()
        tfile = tfile + t4 - t3        
    END IF
    
	J_max = MAX(J_max,MAXVAL(Jx**2+Jy**2+Jz**2))
	E_max = MAX(E_max,MAXVAL(Ex**2+Ey**2+Ez**2))
	B_max = MAX(B_max,MAXVAL(Bx**2+By**2+Bz**2))
    

	
    t6 = MPI_Wtime()
    tavg = tavg + (t6 - t5)

    IF(myrank .EQ. 0) THEN
        PRINT*,''
        WRITE(*,FMT='(" Total Kinetic Energy = ",e6.1)') K_E
        PRINT*,''
	
        PRINT*,'Completed time-step.'
        PRINT*,''
        WRITE(*,FMT='(" Estimated time to completion (minutes) = ",f8.2)') &
          (tavg/DBLE(i_sim-tstart))*DBLE(tstart+maxsteps-i_sim)/60.d0
        PRINT*,''
    END IF
    
    !CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
    
END DO


CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)


! save to restart file
CALL output_restart()


t2 = MPI_Wtime() 
tsim = t2 - t1

PRINT*,''
!PRINT*,'Jx_max, Jy_max=',Jx_max, Jy_max
!PRINT*,'Jy_min, Jy_max=',Jy_min, Jy_max

!IF(myrank .EQ. 0) THEN

PRINT*,'B_max, maxloc = ',SQRT(B_max)/c, MAXLOC(Bz)     
PRINT*,'E_max, maxloc = ',SQRT(E_max), MAXLOC(Ey)  
PRINT*,'J_max = ',SQRT(J_max)   
PRINT*,''
PRINT*,''
PRINT*,'Total number of particle injections: electrons, ions = ',Ntot_inj(1),Ntot_inj(2)
PRINT*,'Total number of escapes: electrons, ions = ',Ntot_esc(1), Ntot_esc(2)
PRINT*,''
PRINT*,''
PRINT*,'Simulation complete!'
PRINT*,'Total Time Elapsed (seconds) = ',tsim
PRINT*,'File I/O time(seconds) = ',tfile
PRINT*,'Fractional File I/O time = ',tfile/tsim
PRINT*,''

!END IF

! deallocate memory for fields and particles
CALL destroy_grid_arrays()

CALL terminate_mpi()



CONTAINS


! initialize grid and other important quantities
SUBROUTINE set_domain()

    INTEGER :: i, j, k, n

    IF(.NOT. current_filter_on .AND. N_filter .GT. 0) THEN
        PRINT*,'Error. Need N_filter = 0  if current filter is switched off...'
        STOP
    END IF

    ! set grid dimensions
    Lx = DBLE(nx)
    Ly = DBLE(ny)
    Lz = DBLE(nz)

    ! set cell size
    dx = 1.d0
    dy = 1.d0
    dz = 1.d0
    
    ! set grid bounds
    xmin = 0.5 * dx
    xmax = xmin + Lx
    ymin = 0.5 * dy
    ymax = ymin + Ly
    zmin = 0.5 * dz
    zmax = zmin + Lz

    ! set some useful grid variables
    mx = nx + 5
    my = ny + 5
    mz = nz + 5
    ix = 1
    iy = mx
    iz = iy * my
    lot = iz * mz
    
    ! set time-step size according to Courant condition (In 3D, COUR = c dt/dx < 1/SQRT(3) ~ 0.5 => c dt < 0.5 => dt < 0.5/c) 
    dt = 1.d0

    ! set mass and charge
    q = (/ q_e, q_i /)
    m = (/ m_e, m_i /)
    
    
    ! ring current initialization
    o3 = 0.d0!-0.25/10.d0 !-0.25
    o2 = -31 * o3 ! set o2 = -n*o3, where n must be an integer (2n = time steps it takes for ring current to reach it's final value)
    o1 = o2
    o = o1
    
    IF(file_restart) THEN
        o1 = 0.d0
        o = -o3*31*(31+1)*(2*31+1)/3.d0     
    END IF
    
    ! ring center
    ie = 4 !nx/2 !18 
    je = ny/2 !ny/2
    ke = nz/2
    
    ! set weights for current filtering
    n = 1
    DO k = -1,1
        DO j = -1,1
            DO i = -1,1
                sm(n) = 0.015625*(2-i*i)*(2-j*j)*(2-k*k)
                ms(n) = i*ix + j*iy + k*iz
                n = n + 1
            END DO
        END DO
    END DO

    ! pre-calculated table of values for a half-Maxwellian velocity distribution (can be used for the flux of particles through the lateral faces)
    DO n =1,64
       table(n) = 0.5d0 * SQRT(LOG(16384.) - 2.D0 * LOG(DBLE(2*n-1)))
    END DO

    ! set taus88 PRNG non-zero seed values (hardcoded for now.)
    s1 = 1234
    s2 = 5678
    s3 = 9123


    PRINT*,''
    PRINT*,'Computational Domain:'
    WRITE(*,FMT='("nx, ny, nz = ",i4,2x,i4,2x,i4)') nx,ny,nz
    WRITE(*,FMT='("Lx = ",f5.1)') Lx
    WRITE(*,FMT='("Ly = ",f5.1)') Ly
    WRITE(*,FMT='("Lz = ",f5.1)') Lz 
    WRITE(*,FMT='("xmin, xmax = ",f5.1,2x,f5.1)') xmin,xmax
    WRITE(*,FMT='("ymin, ymax = ",f5.1,2x,f5.1)') ymin,ymax
    WRITE(*,FMT='("zmin, zmax = ",f5.1,2x,f5.1)') zmin,zmax
    WRITE(*,FMT='("dt = ", f5.1)') dt
   
    PRINT*,''
    WRITE(*,FMT='("mx, my, mz, lot = ",i4,2x,i4,2x,i4,2x,i8)') mx,my,mz,lot
    WRITE(*,FMT='("ix, iy, iz = ",i4,2x,i4,2x,i4)') ix,iy,iz
   
    PRINT*,'q/m : electrons, ions =', q(1)/m(1), q(2)/m(2)

   
    PRINT*,'' 
    
    
END SUBROUTINE set_domain


! Top-level solve routine
SUBROUTINE solve()

    INTEGER :: k
    REAL*8 :: p

    !********************
    ! Half-update B field
    !********************
    PRINT*,'Half updating B field...'
    CALL half_update_B()
    

    !***************
	! Move particles
    !***************
       
    ! first compute grid average fields  
    CALL grid_avg_fields()
      
    PRINT*,'Moving particles...'
    CALL move_particles(i_sim)
	
    !*****************************
    ! Deposit charges and currents
    !*****************************
    PRINT*,'Depositing charges and currents...'
    IF(current_type .EQ. 1) CALL deposit_currents_zigzag()
    IF(current_type .EQ. 2) CALL deposit_currents_esirkepov()
    IF(current_type .EQ. 3) CALL deposit_currents_buneman()
    IF(deposit_charges_On ) CALL deposit_charges()


    ! add ring current
    IF(mycoord(1) .EQ. 1 .AND. mycoord(2) .EQ. 0) CALL add_ring_current()


    ! pre-filter exchange
    CALL exchange_current(2)  ! Filtering needs to be done after exchange. For multiple filtering cycles, boundary currents 
                              ! have to be exchanged before each cycle. Otherwise, current conservation will be violated.           

    ! repeat 'N_filter' times. (Note: In the limit N_filter -> infinity, the binomial filter becomes a gaussian filter)
    DO k = 1 , N_filter      
                     
                                 
        ! Filtering the current makes the spatial distriution smoother and attenuates short wavelength (high k) components of the fourier transform.
        ! Note: Filtering is absolutely crucial. The electromagnetic waves generated by our finite difference equations are dispersive, the wave speed
        ! decreases monotonically for higher k. This results in numerical Cerenkov radiation from particles moving at relativistic speeds. 
        ! (e.g. try moving a single point charge at relativistic speed. It will generate what looks like a bow shock leaving a triangular region
        ! in it's wake containing the outward going e.m. waves. Of course, the particle is moving faster than those waves hence the triangular wake..)
        ! Filtering will attentuate those (aliased) waves at higher frequencies and so mitagate this unphysical effect to a significant extent.
        IF(current_filter_on)THEN
            CALL filter(Jx)
            CALL filter(Jy)
            CALL filter(Jz)
        END IF
        

        ! do a post filter-exchange if more filter cycles need to be done
        IF(k .LT. N_filter)  CALL exchange_current(1) ! only need one boundary layer for next filtering cycle     
        
    END DO

    

    !**************
    ! Evolve fields
    !**************
	PRINT*,'Full updating E and B fields...'	
	CALL full_update_EB() 
    
 
	! compute fv
	!CALL compute_velocity_distribution()


    ! inject solar wind particles
    PRINT*,'Injecting particles...'
    !CALL inject_solarwind()   
    !CALL inject_solarwind()   
    !CALL inject_solarwind()  

    
    !************************************
    ! Apply particle boundary conditions
    !************************************
    PRINT*,'Applying particle boundary conditions...'
    CALL particle_boundary()
    PRINT*,'Exchanging particles...'
    CALL exchange_particles()
    
 
    
    IF(Np_in(1) .GT. Np_max .OR. Np_in(2) .GT. Np_max) THEN
        PRINT*,'Particle array full...terminating program.'
        STOP
    END IF 


END SUBROUTINE solve



SUBROUTINE compute_poynting_flux()
    
    INTEGER :: i,j,k
    REAL*8 :: S, w_d

    S = 0.d0
    
    ! flux through z+, z- faces
    !$OMP PARALLEL DO REDUCTION(+:S)
    DO j = 1, ny
        DO i = 1, nx
            S = S + Ex_grid(i,j,nz)*By_grid(i,j,nz) - Ey_grid(i,j,nz)*Bx_grid(i,j,nz) &
                  - (Ex_grid(i,j,1)*By_grid(i,j,1) - Ey_grid(i,j,1)*Bx_grid(i,j,1))
        END DO
    END DO    
    !$OMP END PARALLEL DO
 
    ! flux through y+, y- faces
    !$OMP PARALLEL DO REDUCTION(+:S)
    DO k = 1, nz
        DO i = 1, nx
            S = S + Ez_grid(i,ny,k)*Bx_grid(i,ny,k) - Ex_grid(i,ny,k)*Bz_grid(i,ny,k) &
                  - (Ez_grid(i,1,k)*Bx_grid(i,1,k) - Ex_grid(i,1,k)*Bz_grid(i,1,k))
        END DO
    END DO    
    !$OMP END PARALLEL DO
    
    ! flux through x+, x- faces
    !$OMP PARALLEL DO REDUCTION(+:S)
    DO k = 1, nz
        DO j = 1, ny
            S = S + Ey_grid(nx,j,k)*Bz_grid(nx,j,k) - Ez_grid(nx,j,k)*By_grid(nx,j,k)
            S = S - (Ey_grid(1,j,k)*Bz_grid(1,j,k) - Ez_grid(1,j,k)*By_grid(1,j,k))
        END DO
    END DO  
    !$OMP END PARALLEL DO

    
    ! poynting flux
    S = S
    
    ! exact time-averaged dipole power radiated (Larmor formula)
    w_d = twopi/1000.d0 
    S_exact = (2.d0/3.d0)*(q_i**2/c**3) * 0.5d0*(10.d0*(w_d)**2)**2 
    
    IF(i_sim .GE. 601 .AND. i_sim .LE. 2600) THEN
       S_avg = S_avg + S
    END IF
    
    WRITE(20,*) i_sim,S,S_exact
    
END SUBROUTINE compute_poynting_flux



SUBROUTINE output_bin(ts)

    INTEGER, INTENT(IN) :: ts
	
    INTEGER :: i, j, k
    REAL*8 :: dv
    REAL*8 :: x, y, Bmag, Emag
    CHARACTER(LEN=180) :: filename1, filename2, filename3, filename4, filename_e_xz, filename_i_xz, &
                         filename_e_yz, filename_i_yz, filename_e_xy, filename_i_xy, filename_cur
    CHARACTER(LEN=6) :: uniti, coordx, coordy
    !INTEGER, PARAMETER :: N_out = 100

    IF(ts<10) THEN
        WRITE(uniti,'(I1.1)') ts
    ELSE IF(ts>=10 .and. ts<100) THEN
        WRITE(uniti,'(I2.2)') ts
    ELSE IF(ts>=100 .and. ts<1000) THEN
        WRITE (uniti,'(I3.3)') ts
    ELSE IF(ts>=1000 .and. ts<10000) THEN
        WRITE (uniti,'(I4.3)') ts
    ELSE IF(ts>=10000 .and. ts<100000) THEN
        WRITE (uniti,'(I5.3)') ts  
    END IF
  
    WRITE(coordx,'(I1.1)') mycoord(1)
    WRITE(coordy,'(I1.1)') mycoord(2)
  
  
    !filename1 = trim('Frames/Snapshots/Bvector_xz_t=')//TRIM(uniti)//TRIM('.dat')    
    !filename2 = trim('Frames/Snapshots/Bvector_xy_t=')//TRIM(uniti)//TRIM('.dat')
    
    filename1 = trim('Frames/Snapshots/Bvector_xz_rank=')//TRIM(coordx)//TRIM(coordy)//TRIM('_t=')//TRIM(uniti)//TRIM('.dat')
    filename2 = trim('Frames/Snapshots/Bvector_xy_rank=')//TRIM(coordx)//TRIM(coordy)//TRIM('_t=')//TRIM(uniti)//TRIM('.dat')
    
    
    filename3 = trim('Frames/Snapshots/density_t=')//TRIM(uniti)//TRIM('.dat')
    filename4 = trim('Frames/Snapshots/density_i_xz_t=')//TRIM(uniti)//TRIM('.dat')
    
    !filename_e_xy = trim('Frames/Snapshots/electrons_xy_t=')//TRIM(uniti)//TRIM('.dat')
    !filename_i_xy = trim('Frames/Snapshots/ions_xy_t=')//TRIM(uniti)//TRIM('.dat')   
    
    filename_e_xy = trim('Frames/Snapshots/electrons_xy_rank=')//TRIM(coordx)//TRIM(coordy)//TRIM('_t=')//TRIM(uniti)//TRIM('.dat')
    filename_i_xy = trim('Frames/Snapshots/ions_xy_rank=')//TRIM(coordx)//TRIM(coordy)//TRIM('_t=')//TRIM(uniti)//TRIM('.dat')
    
    filename_e_xz = trim('Frames/Snapshots/electrons_xz_rank=')//TRIM(coordx)//TRIM(coordy)//TRIM('_t=')//TRIM(uniti)//TRIM('.dat')
    filename_i_xz = trim('Frames/Snapshots/ions_xz_rank=')//TRIM(coordx)//TRIM(coordy)//TRIM('_t=')//TRIM(uniti)//TRIM('.dat')
    
    !filename_cur = trim('Frames/Snapshots/Jvector_xz_t=')//TRIM(uniti)//TRIM('.dat')
    filename_cur = trim('Frames/Snapshots/Jvector_xz_rank=')//TRIM(coordx)//TRIM(coordy)//TRIM('_t=')//TRIM(uniti)//TRIM('.dat')

    !print*,'filename=',filename

    !$OMP PARALLEL
    
    !$OMP SECTIONS PRIVATE(i,j)
    
    !$OMP SECTION
    IF(fields_out) THEN
    
        OPEN(UNIT=11,FILE=filename1, FORM = 'UNFORMATTED')
        OPEN(UNIT=12,FILE=filename2, FORM = 'UNFORMATTED')
        OPEN(UNIT=13,FILE=filename3, FORM = 'UNFORMATTED')
        !OPEN(UNIT=14,FILE=filename4, FORM = 'UNFORMATTED')
       
        ! xy plane
        !DO j = 1, ny
        !    DO i = 1, nx
        !        WRITE(11) i*dx,j*dy,Ex_grid(i,j,nz/2),Ey_grid(i,j,nz/2) 
        !        WRITE(11) i*dx,j*dy,Jx(i,j,nz/2),Jy(i,j,nz/2) 
        !    END DO
        !END DO
 
 
        ! apply rescaling to vectors, i.e. make all vector's have unit length
 
         ! xz plane
         DO j = 1, nz
            DO i = 1, nx
                Bmag = Ex_grid(i,je,j)**2 + Ez_grid(i,je,j)**2 
                Bmag = SQRT(Bmag)
                IF(MOD(j,1) .EQ. 0 .AND. MOD(i,1) .EQ. 0) THEN
                    
                    !IF(Bmag .LT. 200.d0) THEN
                        !WRITE(11)  i*dx+mycoord(1)*Lx ,j*dz,Ex_grid(i,je,j)/Bmag,Ez_grid(i,je,j)/Bmag 
                        WRITE(11)  i*dx+mycoord(1)*Lx ,j*dz,Ex_grid(i,je,j)/Bmag,Ez_grid(i,je,j)/Bmag 
                    !ELSE
                       ! WRITE(11) i*dx,j*dz,0.d0,0.d0
                    !END IF           
                END IF
                
                !Bmag = (Bx_grid(i,ny/2,j)/c)**2 + (By_grid(i,ny/2,j)/c)**2 + (Bz_grid(i,ny/2,j)/c)**2 
                !Bmag = SQRT(Bmag)
             
                WRITE(13) i*dx,j*dz,Bmag
                !WRITE(14) i*dx,j*dz,ni(i,ny/2,j)
                  
            END DO
        END DO
       
        ! xy plane
        
        
        DO j = 1, ny
            DO i = 1, nx
                Bmag = Ex_grid(i,j,ke)**2 + Ey_grid(i,j,ke)**2 
                Bmag = SQRT(Bmag)
                 IF(MOD(j,2) .EQ. 0 .AND. MOD(i,2) .EQ. 0) THEN
                  !   IF(Bmag .LT. 200.d0) THEN
                          WRITE(12) i*dx+mycoord(1)*nx ,j*dy+mycoord(2)*ny,Ex_grid(i,j,ke)/Bmag,Ey_grid(i,j,ke)/Bmag
                         !WRITE(12) i*dx+mycoord(1)*nx ,j*dy+mycoord(2)*ny,Jx(i,j,ke),Jy(i,j,ke)
                   !  ELSE
                 !       WRITE(12) i*dx+mycoord(1)*Lx ,j*dy+mycoord(2)*Ly,0.d0,0.d0 
                  !   END IF
                 END IF      
            END DO
        END DO

        

        CLOSE(UNIT=11)
        CLOSE(UNIT=12)
        CLOSE(UNIT=13)
        !CLOSE(UNIT=14)
                
	END IF	    

    !$OMP SECTION 
    IF(particles_out) THEN
        
        OPEN(UNIT=9,FILE=filename_e_xz, FORM = 'UNFORMATTED')
        OPEN(UNIT=7,FILE=filename_e_xy, FORM = 'UNFORMATTED') 
        IF(ns .EQ. 2) THEN
            OPEN(UNIT=10,FILE=filename_i_xz, FORM = 'UNFORMATTED')
            OPEN(UNIT=8,FILE=filename_i_xy, FORM = 'UNFORMATTED')
        END IF
        
        
        !OPEN(UNIT=9,FILE=filename_e_xz)
        !OPEN(UNIT=7,FILE=filename_e_xy) 
        !IF(ns .EQ. 2) THEN
        !    OPEN(UNIT=10,FILE=filename_i_xz)
        !    OPEN(UNIT=8,FILE=filename_i_xy)
        !END IF
        
        
        ! save particle states
        DO i = 1, Np_in(1)
            !IF( ABS(particles(1,i)%y -ny/2) .LT. 1)
            !IF( ABS(particles(1,i)%z -nz/2) .LT. 1)
            IF( ABS(particles(1,i)%y -je) .LT. 1) WRITE(9) particles(1,i)%x+(mycoord(1))*Lx,particles(1,i)%z            
            IF( ABS(particles(1,i)%z -ke) .LT. 1) WRITE(7) particles(1,i)%x+(mycoord(1))*Lx,particles(1,i)%y+(mycoord(2))*Ly
        END DO
       
        IF(ns .EQ. 2) THEN       
            DO i = 1, Np_in(2) 
                !IF( ABS(particles(2,i)%y -ny/2) .LT. 1) WRITE(10) particles(2,i)%x,particles(2,i)%z
                !IF( ABS(particles(2,i)%z -nz/2) .LT. 1) WRITE(8+15*myrank,*) particles(2,i)%x,particles(2,i)%y
                IF( ABS(particles(2,i)%y -je) .LT. 1) WRITE(10) particles(2,i)%x+(mycoord(1))*Lx,particles(2,i)%z
                IF( ABS(particles(2,i)%z -ke) .LT. 1) WRITE(8) particles(2,i)%x+(mycoord(1))*Lx,particles(2,i)%y+(mycoord(2))*Ly
            END DO
        END IF
        
        CLOSE(UNIT=9)
        CLOSE(UNIT=7)
        IF(ns .EQ. 2) THEN
            CLOSE(UNIT=10)
            CLOSE(UNIT=8)
        END IF
        
	END IF
	
    !$OMP SECTION
    IF(current_out)THEN
        OPEN(UNIT=15,FILE=filename_cur, FORM = 'UNFORMATTED')

        DO j = 1, nz
            DO i = 1, nx
                WRITE(15) i*dx+mycoord(1)*Lx,j*dz,Jx(i,ny/2,j),Jz(i,ny/2,j) 
            END DO
        END DO
   
        CLOSE(UNIT=15)
    END IF
	
	!IF(energy_out) THEN
	!	IF(i_sim .GT. 0) WRITE(13) i_sim, KE, ESE2, KE+ESE2, Px, vd1, vspread1, vd2, vspread2 
	!END IF

    !$OMP END SECTIONS

    !$OMP END PARALLEL 
	

    
    
	!IF(fv_out)THEN
		!dv = (vmax-vmin)/vbins 
		!DO i = 1, vbins
		!	WRITE(17,*) vmin+(i-0.5)*dv, fv(1,i), fv(2,i)
		!END DO
	!END IF
	
END SUBROUTINE output_bin



SUBROUTINE output_restart()

    INTEGER :: i, j, k
    CHARACTER(LEN=180) :: filename_bx, filename_by, filename_bz, &
                          filename_ex, filename_ey, filename_ez, &
                          filename_particles, filename_np, filename_tlast
    CHARACTER(LEN=6) :: uniti, coordx, coordy


    PRINT*,''
    PRINT*,'Writing to restart file...'
    PRINT*,''
  
    WRITE(coordx,'(I1.1)') mycoord(1)
    WRITE(coordy,'(I1.1)') mycoord(2)
  
 
    filename_bx = trim('Frames/Restart/restart_bx_rank=')//TRIM(coordx)//TRIM(coordy)//TRIM('.dat')
    filename_by = trim('Frames/Restart/restart_by_rank=')//TRIM(coordx)//TRIM(coordy)//TRIM('.dat')
    filename_bz = trim('Frames/Restart/restart_bz_rank=')//TRIM(coordx)//TRIM(coordy)//TRIM('.dat')
    filename_ex = trim('Frames/Restart/restart_ex_rank=')//TRIM(coordx)//TRIM(coordy)//TRIM('.dat')
    filename_ey = trim('Frames/Restart/restart_ey_rank=')//TRIM(coordx)//TRIM(coordy)//TRIM('.dat')
    filename_ez = trim('Frames/Restart/restart_ez_rank=')//TRIM(coordx)//TRIM(coordy)//TRIM('.dat')
    filename_particles  = trim('Frames/Restart/restart_particles_rank=')//TRIM(coordx)//TRIM(coordy)//TRIM('.dat') 
    filename_np  = trim('Frames/Restart/restart_np_rank=')//TRIM(coordx)//TRIM(coordy)//TRIM('.dat') 
    filename_tlast  = trim('Frames/Restart/restart_tlast_rank=')//TRIM(coordx)//TRIM(coordy)//TRIM('.dat') 


    OPEN(UNIT=20,FILE=filename_bx, FORM = 'UNFORMATTED')
    OPEN(UNIT=21,FILE=filename_by, FORM = 'UNFORMATTED')
    OPEN(UNIT=22,FILE=filename_bz, FORM = 'UNFORMATTED')
    OPEN(UNIT=23,FILE=filename_ex, FORM = 'UNFORMATTED')
    OPEN(UNIT=24,FILE=filename_ey, FORM = 'UNFORMATTED')
    OPEN(UNIT=25,FILE=filename_ez, FORM = 'UNFORMATTED')
    OPEN(UNIT=26,FILE=filename_particles, FORM = 'UNFORMATTED')
    OPEN(UNIT=27,FILE=filename_np, FORM = 'UNFORMATTED')
    OPEN(UNIT=28,FILE=filename_tlast, FORM = 'UNFORMATTED')


    !$OMP PARALLEL
    
    !$OMP SECTIONS 
    
    !$OMP SECTION    
    WRITE(20) Bx

    !$OMP SECTION
    WRITE(21) By

    !$OMP SECTION
    WRITE(22) Bz
        
    !$OMP SECTION    
    WRITE(23) Ex

    !$OMP SECTION
    WRITE(24) Ey

    !$OMP SECTION
    WRITE(25) Ez
    
    !$OMP SECTION
    WRITE(26) particles

    !$OMP SECTION
    WRITE(27) Np_in
    WRITE(28) i_sim


    !$OMP END SECTIONS

    !$OMP END PARALLEL 


    CLOSE(UNIT=20)
    CLOSE(UNIT=21)
    CLOSE(UNIT=22)
    CLOSE(UNIT=23)
    CLOSE(UNIT=24)
    CLOSE(UNIT=25)
    CLOSE(UNIT=26)
    CLOSE(UNIT=27)
    CLOSE(UNIT=28)



END SUBROUTINE output_restart





SUBROUTINE read_restart()

    INTEGER :: i, j, k
    CHARACTER(LEN=180) :: filename_bx, filename_by, filename_bz, &
                          filename_ex, filename_ey, filename_ez, &
                          filename_particles, filename_np, filename_tlast
    CHARACTER(LEN=6) :: uniti, coordx, coordy


  
    PRINT*,''
    PRINT*,'Reading from restart file...'
    PRINT*,''
    
    WRITE(coordx,'(I1.1)') mycoord(1)
    WRITE(coordy,'(I1.1)') mycoord(2)
  
 
    filename_bx = trim('Frames/Restart/restart_bx_rank=')//TRIM(coordx)//TRIM(coordy)//TRIM('.dat')
    filename_by = trim('Frames/Restart/restart_by_rank=')//TRIM(coordx)//TRIM(coordy)//TRIM('.dat')
    filename_bz = trim('Frames/Restart/restart_bz_rank=')//TRIM(coordx)//TRIM(coordy)//TRIM('.dat')
    filename_ex = trim('Frames/Restart/restart_ex_rank=')//TRIM(coordx)//TRIM(coordy)//TRIM('.dat')
    filename_ey = trim('Frames/Restart/restart_ey_rank=')//TRIM(coordx)//TRIM(coordy)//TRIM('.dat')
    filename_ez = trim('Frames/Restart/restart_ez_rank=')//TRIM(coordx)//TRIM(coordy)//TRIM('.dat')
    filename_particles  = trim('Frames/Restart/restart_particles_rank=')//TRIM(coordx)//TRIM(coordy)//TRIM('.dat') 
    filename_np  = trim('Frames/Restart/restart_np_rank=')//TRIM(coordx)//TRIM(coordy)//TRIM('.dat') 
    filename_tlast  = trim('Frames/Restart/restart_tlast_rank=')//TRIM(coordx)//TRIM(coordy)//TRIM('.dat') 


    OPEN(UNIT=20,FILE=filename_bx, FORM = 'UNFORMATTED', STATUS='OLD')
    OPEN(UNIT=21,FILE=filename_by, FORM = 'UNFORMATTED', STATUS='OLD')
    OPEN(UNIT=22,FILE=filename_bz, FORM = 'UNFORMATTED', STATUS='OLD')
    OPEN(UNIT=23,FILE=filename_ex, FORM = 'UNFORMATTED', STATUS='OLD')
    OPEN(UNIT=24,FILE=filename_ey, FORM = 'UNFORMATTED', STATUS='OLD')
    OPEN(UNIT=25,FILE=filename_ez, FORM = 'UNFORMATTED', STATUS='OLD')
    OPEN(UNIT=26,FILE=filename_particles, FORM = 'UNFORMATTED', STATUS='OLD')
    OPEN(UNIT=27,FILE=filename_np, FORM = 'UNFORMATTED', STATUS='OLD')
    OPEN(UNIT=28,FILE=filename_tlast, FORM = 'UNFORMATTED', STATUS='OLD')

    !$OMP PARALLEL
    
    !$OMP SECTIONS 
    
    !$OMP SECTION    
    READ(20) Bx

    !$OMP SECTION
    READ(21) By

    !$OMP SECTION
    READ(22) Bz
        
    !$OMP SECTION    
    READ(23) Ex

    !$OMP SECTION
    READ(24) Ey

    !$OMP SECTION
    READ(25) Ez

    !$OMP SECTION
    READ(26) particles

    !$OMP SECTION
    READ(27) Np_in
    READ(28) tstart
    

    !$OMP END SECTIONS

    !$OMP END PARALLEL 


    CLOSE(UNIT=20)
    CLOSE(UNIT=21)
    CLOSE(UNIT=22)
    CLOSE(UNIT=23)
    CLOSE(UNIT=24)
    CLOSE(UNIT=25)
    CLOSE(UNIT=26)
    CLOSE(UNIT=27)
    CLOSE(UNIT=28)


    PRINT*,''
    PRINT*,'Completed reading from restart file.'
    PRINT*,''
    

END SUBROUTINE read_restart




END PROGRAM driver_plasma_3d
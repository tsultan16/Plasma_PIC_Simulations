PROGRAM driver_plasma_3d

USE constants_mod
USE data_mod
USE init_mod
USE fieldSolver_mod
USE particleMover_mod
USE OMP_LIB ! openmp runtime library

IMPLICIT NONE


INTEGER :: ii, i, j 
REAL*8 :: x1, y1, x2, y2, Jx_max = 0.d0, Jy_max = 0.d0, Jx_min = 0.d0, Jy_min = 0.d0, J_max = 0.d0, &
          E_max = 0.d0, B_max = 0.d0, t1, t2, t3, t4, t5, t6, tavg = 0.d0, tsim = 0.d0, tfile = 0.d0, &
          S_avg = 0.d0, S_exact = 0.d0
CHARACTER*60 :: file_particles_e,file_particles_i, file_horcut, file_vercut, file_energy, file_fv 


!***********************I/O Files*****************************
file_particles_e = 'Frames/particles_e.dat'
file_particles_i = 'Frames/particles_i.dat'
file_horcut = 'Frames/fields_horcut.dat'
file_vercut = 'Frames/fields_vercut.dat'
file_energy = 'Frames/energy.dat'
file_fv = 'Frames/fv.dat'

IF(particles_out) THEN
    OPEN(UNIT=9, FILE=file_particles_e, FORM = 'UNFORMATTED')
    IF(ns .EQ. 2) OPEN(UNIT=10, FILE=file_particles_i, FORM = 'UNFORMATTED')
END IF
IF(energy_out) OPEN(UNIT=13, FILE=file_energy, FORM = 'UNFORMATTED')
IF(fv_out) OPEN(UNIT=17, FILE=file_fv, FORM = 'UNFORMATTED')

OPEN(UNIT=20, FILE='Frames/poynting.txt')
!*************************************************************


! check that there are more particles than cells
IF(.NOT. override_checks) THEN

    IF(N .LT. nx*ny) THEN
        PRINT*,'Need N > nx*ny.'
        PRINT*,'Terminating program..'
        STOP
    END IF

    IF(MOD(N,nx*ny) .NE. 0) THEN
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

! allocate memory for fields
CALL create_grid_arrays()


! initialize particle distribution and state
CALL init_particles()


! add sinusoidal perturbation to electrons
CALL add_perturbations()

! compute initial charge distribution
CALL setrho()


! intialize uniform fields (Ey and Bz)
Ex = 0.d0
Ey = 0.d0
Ez = 0.d0
Bx = 0.d0
By = 0.d0
Bz = 0.d0

!DO i = 0,nx+1
!    DO j = 0, ny+1
!        Bz(i,j) = 1.d0 + (i*dx)*1.0d0  ! linearly increasing in y-direction
!    END DO
!END DO
CALL grid_avg_fields()

! set v(-dt/2)
CALL setv()



CALL output_bin(0)

t1 = OMP_get_wtime() 

! enter simulation loop
DO i_sim = 1, maxsteps
	
	PRINT*,''
	WRITE(*,FMT='("Time step = ",i4,", dt = ",f4.1,", Completion = ",f5.2," %")') &
         i_sim, dt, 100.*REAL(i_sim)/REAL(maxsteps)
    PRINT*,''
	
    t5 = OMP_get_wtime() 

    ! advance by one time-step 
    CALL solve() 
	
    ! compute poynting flux through computational boundary
    !CALL compute_poynting_flux() 
    
	! save to file
	IF(MOD(i_sim-1,tskip) .EQ. 0)THEN
        t3 = OMP_get_wtime() 
        CALL output_bin(i_sim-1)
        t4 = OMP_get_wtime()
        tfile = tfile + t4 - t3        
    END IF
    !J_max = MAX(J_max,MAXVAL(Jx**2+Jy**2+Jz**2))
	!E_max = MAX(E_max,MAXVAL(Ex_grid**2+Ey_grid**2+Ez_grid**2))
	B_max = MAX(B_max,MAXVAL(Bx_grid**2+By_grid**2+Bz_grid**2))


	PRINT*,''
	WRITE(*,FMT='(" Total Kinetic Energy = ",e6.1)') KE
	PRINT*,''
	
    t6 = OMP_get_wtime() 
    tavg = tavg + (t6 - t5)

    PRINT*,'Completed time-step.'
    PRINT*,''
	WRITE(*,FMT='(" Expected time to completion (minutes) = ",f8.2)') (tavg/DBLE(i_sim))*DBLE(maxsteps-i_sim)/60.d0
	PRINT*,''
    
    
END DO

t2 = OMP_get_wtime() 
tsim = t2 - t1

PRINT*,''
!PRINT*,'Jx_max, Jy_max=',Jx_max, Jy_max
!PRINT*,'Jy_min, Jy_max=',Jy_min, Jy_max
!PRINT*,'E_max = ',SQRT(E_max)
PRINT*,'B_max = ',SQRT(B_max)
!PRINT*,'J_max = ',SQRT(J_max)
PRINT*,''

! deallocate memory for fields and particles
CALL destroy_grid_arrays()

PRINT*,''
PRINT*,'Simulation complete!'
PRINT*,'Time Elapsed (seconds) = ',tsim
PRINT*,'Fractional File I/O time(seconds) = ',tfile/tsim
PRINT*,''


CLOSE(UNIT=10)
CLOSE(UNIT=11)
CLOSE(UNIT=12)
CLOSE(UNIT=13)
CLOSE(UNIT=14)
CLOSE(UNIT=15)
CLOSE(UNIT=16)
CLOSE(UNIT=17)
CLOSE(UNIT=111)
CLOSE(UNIT=112)
CLOSE(UNIT=20)

CONTAINS


SUBROUTINE set_domain()

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

    ! ring current initialization
    o3 = -0.25
    o2 = -31 * o3 ! set o2 = -n*o3, where n must be an integer (2n = time steps it takes for ring current to reach it's final value)
    o1 = o2
    o = o1

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
    
    
END SUBROUTINE set_domain


! Top-level solve routine
SUBROUTINE solve()

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
    IF(deposit_charges_On ) CALL deposit_charges()
  
    !************************************
    ! Apply particle boundary conditions
    !************************************
    CALL particle_boundary()

    !**************
    ! Evolve fields
    !**************
	PRINT*,'Full updating E and B fields...'	
	CALL full_update_EB()
    
    
	! compute fv
	!CALL compute_velocity_distribution()
	

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
    S_exact = (2.d0/3.d0)*(ee**2/c**3) * 0.5d0*(10.d0*(w_d)**2)**2 
    
    IF(i_sim .GE. 601 .AND. i_sim .LE. 2600) THEN
       S_avg = S_avg + S
    END IF
    
    WRITE(20,*) i_sim,S,S_exact
    
END SUBROUTINE compute_poynting_flux


SUBROUTINE create_grid_arrays()


	ALLOCATE(particles(ns,N))
    ALLOCATE(rho(-1:nx+3,-1:ny+3,-1:nz+3), phi(-1:nx+3,-1:ny+3,-1:nz+3))
    ALLOCATE(Ex(-1:nx+3,-1:ny+3,-1:nz+3), Ey(-1:nx+3,-1:ny+3,-1:nz+3),Ez(-1:nx+3,-1:ny+3,-1:nz+3))  ! -2,-1,0,n+1,n+1,n+3 are ghost cells
    ALLOCATE(Bx(-1:nx+3,-1:ny+3,-1:nz+3), By(-1:nx+3,-1:ny+3,-1:nz+3),Bz(-1:nx+3,-1:ny+3,-1:nz+3))  
    ALLOCATE(Ex_grid(-1:nx+1,-1:ny+1,-1:nz+1), Ey_grid(-1:nx+1,-1:ny+1,-1:nz+1), Ez_grid(-1:nx+1,-1:ny+1,-1:nz+1))
    ALLOCATE(Bx_grid(-1:nx+1,-1:ny+1,-1:nz+1), By_grid(-1:nx+1,-1:ny+1,-1:nz+1), Bz_grid(-1:nx+1,-1:ny+1,-1:nz+1))  
    ALLOCATE(Jx(-1:nx+3,-1:ny+3,-1:nz+3), Jy(-1:nx+3,-1:ny+3,-1:nz+3), Jz(-1:nx+3,-1:ny+3,-1:nz+3))
    ALLOCATE(ne(0:nx+1,0:ny+1,0:nz+1), ni(0:nx+1,0:ny+1,0:nz+1))
    ALLOCATE(buffer(-1:nx+1,-1:ny+1,-1:nz+1))
    ALLOCATE(fv(1:ns,1:vbins))
    ALLOCATE(pos_buffer(nx,N,3))

    rho = 0.d0
    phi =0.d0
    Jx = 0.d0
    Jy = 0.d0
    Jz = 0.d0

END SUBROUTINE create_grid_arrays


SUBROUTINE destroy_grid_arrays()

    DEALLOCATE(rho, phi, Ex, Ey, Ez, Bx, By, Bz, ne, ni, particles,fv)
    DEALLOCATE(Jx, Jy, Jz)
    DEALLOCATE(buffer)
    DEALLOCATE(pos_buffer)

END SUBROUTINE destroy_grid_arrays



! computes space-averaged velocity distribution function
! and drift and velocity spreads for each species
SUBROUTINE compute_velocity_distribution()

	INTEGER :: i, j, i_v
    REAL*8 :: dv, vv
	
	! clear fv, vd vt
	fv = 0.d0
	vd1 = 0.d0
	vd2 = 0.d0
	vspread1 = 0.d0
	vspread2 = 0.d0
	
	dv = (vmax-vmin)/vbins 
    ! loop over all particles
    DO i = 1, ns
		DO j = 1, N
            ! compute velocity bin index
			vv = (particles(i,j)%ux-vmin)/dv
		    i_v = vv 
			
			! place particle in bin (NOTE: this is a first order interpolation)
			IF(i_v .GE. 1 .AND. i_v .LE. vbins) THEN
				fv(i,i_v)   = fv(i,i_v) + (1.d0 - (vv-i_v))
				IF(i_v .LE. vbins-1) fv(i,i_v+1) = fv(i,i_v+1) + (vv-i_v)
			END IF

            IF(i .EQ. 1) THEN
				vd1 = vd1 + particles(i,j)%ux
			ELSE IF(i .EQ. 2) THEN
				vd2 = vd2 + particles(i,j)%ux
			END IF
			
		END DO
	END DO
	
	vd1 = vd1/N
	vd2 = vd2/N
	
	! compute velocity spreads 
	DO j = 1,N
		vspread1 = vspread1 +  particles(1,j)%ux**2
		IF(ns .EQ.  2) vspread2 = vspread2 +  particles(2,j)%ux**2		
	END DO
	vspread1 = vspread1/N
	vspread2 = vspread2/N

	vspread1 = vspread1-vd1**2
	vspread2 = vspread2-vd1**2
	
	
	
END SUBROUTINE compute_velocity_distribution


SUBROUTINE output_bin(ts)

    INTEGER, INTENT(IN) :: ts
	
    INTEGER :: i, j, k
    REAL*8 :: dv
    REAL*8 :: x, y, Bmag
    CHARACTER(LEN=40) :: filename1, filename2, filename3, filename4
    CHARACTER(LEN=6) :: uniti
    INTEGER, PARAMETER :: N_out = 100

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
  
    filename1 = trim('Frames/Fields/vector_xy_t=')//TRIM(uniti)//TRIM('.dat')
    filename2 = trim('Frames/Fields/fields_vercut_t=')//TRIM(uniti)//TRIM('.dat')
    filename3 = trim('Frames/Fields/fields_xz_t=')//TRIM(uniti)//TRIM('.dat')
    filename4 = trim('Frames/Fields/B_xz_t=')//TRIM(uniti)//TRIM('.dat')
    
    !print*,'filename=',filename

    !$OMP PARALLEL
    
    !$OMP SECTIONS PRIVATE(i,j)
    
    !$OMP SECTION
    IF(fields_out) THEN
    
        OPEN(UNIT=11,FILE=filename1, FORM = 'UNFORMATTED')
        OPEN(UNIT=12,FILE=filename2, FORM = 'UNFORMATTED')
        OPEN(UNIT=13,FILE=filename3, FORM = 'UNFORMATTED')
        OPEN(UNIT=14,FILE=filename4, FORM = 'UNFORMATTED')
       
        ! xy plane
        !DO j = 1, ny
        !    DO i = 1, nx
        !        WRITE(11) i*dx,j*dy,Ex_grid(i,j,nz/2),Ey_grid(i,j,nz/2) 
        !        WRITE(11) i*dx,j*dy,Jx(i,j,nz/2),Jy(i,j,nz/2) 
        !    END DO
        !END DO
 
         ! xz plane
         DO j = 1, nz
            DO i = 1, nx
                WRITE(11) i*dx,j*dz,Bx_grid(i,ny/2,j)/c,Bz_grid(i,ny/2,j)/c 
                !WRITE(11) i*dx,j*dz,Ex_grid(i,ny/2,j),Ez_grid(i,ny/2,j) 
                
                Bmag = (Bx_grid(i,ny/2,j)/c)**2 + (By_grid(i,ny/2,j)/c)**2 + (Bz_grid(i,ny/2,j)/c)**2 
                Bmag = SQRT(Bmag)
                WRITE(13) i*dx,j*dz,ne(i,ny/2,j),ni(i,ny/2,j)
                WRITE(14) i*dx,j*dz,Bmag  
            END DO
        END DO
       
        !DO j = 1, ny
        !    DO i = 1, nx
        !        IF(MOD(j,10) .EQ. 0 .AND. MOD(i,10) .EQ. 0) WRITE(12) i*dx,j*dy,Ex(i,j),Ey(i,j) 
        !    END DO
        !END DO
        
        DO i = 1, nz
            WRITE(12) DBLE(i)*dz,Bz(INT(0.5*nx),ny/2,i)/c 
        END DO
        
        
        CLOSE(UNIT=11)
        CLOSE(UNIT=12)
        CLOSE(UNIT=13)
        CLOSE(UNIT=14)
                
	END IF	    

    !$OMP SECTION 
    IF(particles_out) THEN
        ! save particle states
        DO i = 1, N_out
          
           j = i * (N / N_out)
        
            IF(.NOT. particles(1,j)%oob)THEN !.AND. INT(particles(1,i)%y) .EQ. ny/2) THEN
                WRITE(9) particles(1,j)%x,particles(1,j)%z
            ELSE   
                WRITE(9) -1.d0,-1.d0
            END IF            
                IF(ns .EQ. 2) THEN
                    IF(.NOT. particles(2,j)%oob)THEN ! .AND. INT(particles(2,i)%y) .EQ. ny/2) THEN
                        WRITE(10) particles(2,j)%x,particles(2,j)%z
                    ELSE   
                        WRITE(10) -1.d0,-1.d0
                    END IF            
            END IF
        END DO
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


END PROGRAM driver_plasma_3d
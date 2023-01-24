PROGRAM driver_plasma_2d

USE constants_mod
USE data_mod
USE init_mod
USE fieldSolver_mod
USE particleMover_mod
IMPLICIT NONE


INTEGER :: i_sim, ii, i, j 
REAL*8 :: x1, y1, x2, y2, Jx_max = 0.d0, Jy_max = 0.d0, Jx_min = 0.d0, Jy_min = 0.d0, J_max = 0.d0, E_max = 0.d0

IF(particles_out) THEN
    OPEN(UNIT=9, FILE='Frames/particles_e.txt')
    IF(ns .EQ. 2) OPEN(UNIT=10, FILE='Frames/particles_i.txt')
END IF
IF(fields_out) THEN
    OPEN(UNIT=111, FILE='Frames/fields_horcut.txt')
    OPEN(UNIT=112, FILE='Frames/fields_vercut.txt')  
END IF  
IF(fieldsk_out) OPEN(UNIT=12, FILE='Frames/fieldsk.txt')
IF(energy_out) OPEN(UNIT=13, FILE='Frames/energy.txt')
IF(modes_out) OPEN(UNIT=14, FILE='Frames/modes.txt')
!OPEN(UNIT=15, FILE='frequencies.txt' )
!OPEN(UNIT=16, FILE='single_particle.txt' )
IF(fv_out) OPEN(UNIT=17, FILE='fv.txt')

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



CALL output(0)

! enter simulation loop
DO i_sim = 1, maxsteps
	
	PRINT*,''
	PRINT*,'Time step, dt, % complete = ',i_sim,dt,100.*REAL(i_sim)/REAL(maxsteps)
	PRINT*,''
	
    ! advance by one time-step 
    CALL solve() 
	
	! save to file
	IF(MOD(i_sim-1,tskip) .EQ. 0) CALL output(i_sim-1)
	
    !Jx_max = MAX(Jx_max,MAXVAL(Jx))
    !Jy_max = MAX(Jy_max,MAXVAL(Jy))
	!Jx_min = MIN(Jx_min,MINVAL(Jx))
    !Jy_min = MIN(Jy_min,MINVAL(Jy))
	!J_max = MAX(J_max,MAXVAL(Jx**2+Jy**2))
	E_max = MAX(E_max,MAXVAL(Ex_grid**2+Ey_grid**2+Ez_grid**2))

	PRINT*,''
	PRINT*,'Total KE, ESE, Px, Py = ',KE, ESE, Px, Py
	PRINT*,''
	
    PRINT*,'Completed time-step.'

END DO

PRINT*,''
!PRINT*,'Jx_max, Jy_max=',Jx_max, Jy_max
!PRINT*,'Jy_min, Jy_max=',Jy_min, Jy_max
PRINT*,'E_max = ',SQRT(E_max)
!PRINT*,''


! deallocate memory for fields and particles
CALL destroy_grid_arrays()

PRINT*,'Simulation complete!'

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



CONTAINS


SUBROUTINE set_domain()


    ! set grid dimensions
    Lx = 1.d0
    Ly = 1.d0
    Lz = 1.d0

    ! set cell size
    dx = Lx/REAL(nx)
    dy = Ly/REAL(ny)
    dz = Ly/REAL(nz)
    
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
    dt = COUR*MIN(dx,dy,dz)/c 

    PRINT*,''
    PRINT*,'Computational Domain:'
    PRINT*,'nx, ny, nz =',nx,ny,nz
    PRINT*,'Lx = ',Lx 
    PRINT*,'Ly = ',Ly
    PRINT*,'Lz = ',Lz 
    PRINT*,'xmin, xmax =',xmin,xmax
    PRINT*,'ymin, ymax =',ymin,ymax
    PRINT*,'zmin, zmax =',zmin,zmax
    PRINT*,'dt = ',dt
    PRINT*,''
    
    
END SUBROUTINE set_domain


! Top-level solve routine
SUBROUTINE solve()

    !****************
	! Move particles
    !***************
    PRINT*,'Moving particles...'
    CALL move_particles(i_sim)
	
    !*****************************
    ! Deposit charges and currents
    !*****************************
    PRINT*,'Depositing charges and currents...'
    CALL deposit_currents()
    CALL deposit_charges()
  
    !************************************
    ! Apply particle boundary conditions
    !************************************
    CALL particle_boundary()

    !**************
    ! Evolve fields
    !**************
	PRINT*,'Updating fields...'	
	CALL compute_fields()
	
	
	! compute fv
	!CALL compute_velocity_distribution()
	

END SUBROUTINE solve



SUBROUTINE create_grid_arrays()


	ALLOCATE(particles(ns,N))
    ALLOCATE(rho(-1:nx+3,-1:ny+3,-1:nz+3), phi(-1:nx+3,-1:ny+3,-1:nz+3))
    ALLOCATE(Ex(-2:nx+3,-2:ny+3,-2:nz+3), Ey(-2:nx+3,-2:ny+3,-1:nz+3),Ez(-2:nx+3,-2:ny+3,-2:nz+3))  ! -2,-1,0,n+1,n+1,n+3 are ghost cells
    ALLOCATE(Bx(-2:nx+3,-2:ny+3,-2:nz+3), By(-2:nx+3,-2:ny+3,-1:nz+3),Bz(-2:nx+3,-2:ny+3,-2:nz+3))  
    ALLOCATE(Ex_grid(-1:nx+1,-1:ny+1,-1:nz+1), Ey_grid(-1:nx+1,-1:ny+1,-1:nz+1), Ez_grid(-1:nx+1,-1:ny+1,-1:nz+1))
    ALLOCATE(Bx_grid(-1:nx+1,-1:ny+1,-1:nz+1), By_grid(-1:nx+1,-1:ny+1,-1:nz+1), Bz_grid(-1:nx+1,-1:ny+1,-1:nz+1))  
    ALLOCATE(Jx(-2:nx+3,-2:ny+3,-2:nz+3), Jy(-2:nx+3,-2:ny+3,-2:nz+3), Jz(-2:nx+3,-2:ny+3,-2:nz+3))
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



SUBROUTINE output(ts)

    INTEGER, INTENT(IN) :: ts
	
    INTEGER :: i, j, k
    REAL*8 :: dv
    REAL*8 :: x, y
    CHARACTER(LEN=40) :: filename1, filename2
    CHARACTER(LEN=6) :: uniti

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
  
    filename1 = trim('Frames/Fields/fields_xy_t=')//TRIM(uniti)//TRIM('.txt')
    filename2 = trim('Frames/Fields/E_t=')//TRIM(uniti)//TRIM('.txt')
    
    !print*,'filename=',filename

    IF(fields_out) THEN
       OPEN(UNIT=11,FILE=filename1)
       !OPEN(UNIT=12,FILE=filename2)
    
       ! save rho, phi and E
        DO j = 1, nz
            DO i = 1, ny
                WRITE(11,*) i*dy,j*dz,Ey_grid(nx/2,i,j),Ez_grid(nx/2,i,j) 
            END DO
        END DO
        
        !DO j = 1, ny
        !    DO i = 1, nx
        !        IF(MOD(j,10) .EQ. 0 .AND. MOD(i,10) .EQ. 0) WRITE(12,*) i*dx,j*dy,Ex(i,j),Ey(i,j) 
        !    END DO
        !END DO
        
        CLOSE(UNIT=11)
    
        !CLOSE(UNIT=12)
    
        
        !DO i = 1, nx
        !    WRITE(111,*) i*dx,rho(i,ny/2),ne(i,ny/2),Ex(i,ny/2),Ey(i,ny/2)!, ne(i),ni(i) 
        !END DO
    
        !DO i = 1, ny
        !    WRITE(112,*) i*dy,Bz(nx/2,i) 
        !END DO
        
	END IF	
    
    IF(particles_out) THEN
	! save particle states
	DO j = 1, N
		WRITE(9,*) particles(1,j)%y,particles(1,j)%z, particles(1,j)%ux,particles(1,j)%uy
		IF(ns .EQ. 2) WRITE(10,*) particles(2,j)%y,particles(2,j)%z, particles(2,j)%ux,particles(2,j)%uy
	END DO
	
	!WRITE(16,*) particles(1,N/2)%x,particles(1,N/2)%y, particles(2,N/2)%x,particles(2,N/2)%y

	END IF
	
	
	IF(energy_out) THEN
		IF(i_sim .GT. 0) WRITE(13,*) i_sim, KE, ESE2, KE+ESE2, Px, vd1, vspread1, vd2, vspread2 
	END IF

	
	!IF(fv_out)THEN
		!dv = (vmax-vmin)/vbins 
		!DO i = 1, vbins
		!	WRITE(17,*) vmin+(i-0.5)*dv, fv(1,i), fv(2,i)
		!END DO
	!END IF
	
END SUBROUTINE output

END PROGRAM driver_plasma_2d
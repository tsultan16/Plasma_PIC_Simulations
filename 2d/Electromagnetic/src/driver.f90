PROGRAM driver_plasma_2d

USE constants_mod
USE data_mod
USE init_mod
USE fieldSolver_mod
USE particleMover_mod
IMPLICIT NONE


INTEGER :: i_sim, ii, i, j 
REAL*8 :: w_mode, w_exact
REAL*8 :: x1, y1, x2, y2, Jx_max = 0.d0, Jy_max = 0.d0, Jx_min = 0.d0, Jy_min = 0.d0, J_max = 0.d0, E_max = 0.d0
REAL*8 :: rho_old(0:nx+1,0:ny+1)

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


! allocate memory for fields
ALLOCATE(rho(-1:nx+3,-1:ny+3), phi(-1:nx+3,-1:ny+3))
ALLOCATE(Ex(-1:nx+1,-1:ny+1), Ey(-1:nx+3,-1:ny+3),Bz(-1:nx+3,-1:ny+3))  ! -1,0,nx+1,nx+1,nx+3 (ny+1,ny+1,ny+3) are ghost cells
ALLOCATE(Ex_grid(0:nx+1,0:ny+1), Ey_grid(0:nx+1,0:ny+1),Bz_grid(0:nx+1,0:ny+1))
ALLOCATE(Jx(-1:nx+3,-1:ny+3), Jy(-1:nx+3,-1:ny+3))
ALLOCATE(ne(0:nx+1,0:ny+1), ni(0:nx+1,0:ny+1))
ALLOCATE(buffer(-1:nx+1,-1:ny+1))

ALLOCATE(rhok(-nx/2:nx/2 - 1,-ny/2:ny/2-1 ,2), phik(-nx/2:nx/2 - 1,-ny/2:ny/2-1 ,2), &
         Ksqr(-nx/2:nx/2 - 1,-ny/2:ny/2-1 ,2),ESEk(-nx/2:nx/2 - 1,-ny/2:ny/2-1),SMk(-nx/2:nx/2 - 1,-ny/2:ny/2-1))
ALLOCATE(ESEkt(1:maxsteps,-nx/2:nx/2 - 1,-ny/2:ny/2-1))

ALLOCATE(fv(1:ns,1:vbins))

Jx = 0.d0
Jy = 0.d0

! compute dx, dy
dx = Lx/REAL(nx)
dy = Ly/REAL(nx)


! set grid bounds
xmin = 0.5*dx
xmax = xmin+Lx
ymin = 0.5*dy
ymax = ymin+Ly


! set dt 
dt = COUR*dx/c 
!0.01*twopi/wp_e

! initialize particle distribution and state
CALL init_particles()


! add sinusoidal perturbation to electrons
CALL add_perturbations()

! compute initial charge distribution
CALL setrho()

! compute initial electrostatic field

! intialize uniform fields (Ey and Bz)
Ex = 0.d0
Ey = 0.d0
Bz = 0.d0
!DO i = 0,nx+1
!    DO j = 0, ny+1
!        Bz(i,j) = 1.d0 + (i*dx)*1.0d0  ! linearly increasing in y-direction
!    END DO
!END DO
!CALL poisson_solve()
CALL grid_avg_fields()

! set v(-dt/2)
CALL setv()



CALL output(0)

PRINT*,'xmin,xmax,dx=',xmin,xmax,dx
PRINT*,'ymin,ymax,dy=',ymin,ymax,dy

! enter simulation loop
DO i_sim = 1, maxsteps
	
	PRINT*,''
	PRINT*,'Time step, dt, % complete = ',i_sim,dt,100.*REAL(i_sim)/REAL(maxsteps)
	PRINT*,''
	
    !rho_old = rho
    
	PRINT*,'Moving particles...'
	! move particles
    CALL move_particles(i_sim)
	
    
	PRINT*,'Updating fields...'
	! evolve fields
	!CALL compute_fields()
	
	ESEkt(i_sim,:,:) = ESEk(:,:)
	
	! compute fv
	!CALL compute_velocity_distribution()
	
	! save to file
	IF(MOD(i_sim-1,tskip) .EQ. 0) CALL output(i_sim-1)
	
    !Jx_max = MAX(Jx_max,MAXVAL(Jx))
    !Jy_max = MAX(Jy_max,MAXVAL(Jy))
	!Jx_min = MIN(Jx_min,MINVAL(Jx))
    !Jy_min = MIN(Jy_min,MINVAL(Jy))
	!J_max = MAX(J_max,MAXVAL(Jx**2+Jy**2))
	E_max = MAX(E_max,MAXVAL(Ex**2+Ey**2))

	PRINT*,''
	PRINT*,'Total KE, ESE, Px, Py = ',KE, ESE, Px, Py
	PRINT*,''
	
    PRINT*,'Completed time-step.'

END DO

PRINT*,''
!PRINT*,'Jx_min, Jx_max=',Jx_min, Jx_max
!PRINT*,'Jy_min, Jy_max=',Jy_min, Jy_max
PRINT*,'E_max = ',SQRT(E_max)
PRINT*,''

! compute mode frequencies
!w_exact = SQRT(wp_e**2+wc_e**2)
!PRINT*,''
!PRINT*,'SQRT*(wpe^2 + wce^2)=',w_exact
!PRINT*,''
!DO ii = 1, 13
!	CALL compute_frequency(ii, w_mode)
!	WRITE(15,*) ii*twopi/nx, w_mode, w_exact
!END DO

! deallocate memory for fields and particles
DEALLOCATE(rho, Ex, Ey, Bz, phi, ne, ni, particles,fv)
DEALLOCATE(Jx,Jy)
DEALLOCATE(buffer)
DEALLOCATE(rhok, phik, Ksqr,ESEk,SMk, ESEkt) 

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


SUBROUTINE charge_conservation()

    INTEGER :: i,j
    REAL*8 :: f(1:nx, 1:ny)
    
    ! compute net charge flux out of each cell
    DO i = 1, nx
        DO j = 1, ny
            f(i,j) =(rho(i,j)-rho_old(i,j))+ (Jx(i,j)-Jx(i-1,j)) + (Jy(i,j)-Jy(i,j-1))   
        END DO
    END DO    


    PRINT*,''
	PRINT*,'rho(t+dt) - rho(t) + div(J) = '
	DO j = ny, 1, -1
        DO i = 1, nx
            WRITE(*,FMT='(1f8.2)',ADVANCE='NO') f(i,j)
        END DO
        PRINT*,''
    END DO
    

END SUBROUTINE charge_conservation



SUBROUTINE compute_frequency(modex, modey, w_mode)

	INTEGER, INTENT(IN)   ::  modex, modey
	REAL*8, INTENT(INOUT) ::  w_mode
	INTEGER :: i,j

    INTEGER :: maxpos1, maxpos2	
	REAL*8  :: maxval1, period
	LOGICAL :: maxfound
	
	! find first local max
	maxfound = .FALSE.
	maxval1= ESEkt(1,modex,modey)	
    maxpos1 = 1	
	DO i = 2, maxSteps
	
	    IF(maxfound) EXIT
	 		
		IF(ESEkt(i,modex,modey) .GT. ESEkt(i-1,modex,modey) .AND. &
           ESEkt(i,modex,modey) .GT. ESEkt(i+1,modex,modey) ) THEN
		    maxval1 = ESEkt(i,modex,modey)
			maxpos1 = i 	
			maxfound = .TRUE.
		END IF			
				
	END DO

	PRINT*,'Mode#',modex,modey,', MAXPOS1 =',maxpos1
	
    ! find second local max
	maxfound = .FALSE.
	maxval1= ESEkt(maxpos1+1,modex,modey)		
	maxpos2 = maxpos1+1
	DO i = maxpos1+1, maxSteps

	    IF(maxfound) EXIT
	 		
		IF(ESEkt(i,modex,modey) .GT. ESEkt(i-1,modex,modey) .AND. &
           ESEkt(i,modex,modey) .GT. ESEkt(i+1,modex,modey) ) THEN
		    maxval1 = ESEkt(i,modex,modey)
			maxpos2 = i 	
			maxfound = .TRUE.
		END IF			
		
	END DO

	PRINT*,'Mode#',modex,modey,', MAXPOS2 =',maxpos2

    period = 2.d0 *(maxpos2-maxpos1)*dt ! multiply by 2 because ESEk(k)~ phik*rhok ~ exp(iwt)*exp(iwt) ~ cos(2wt) 
    w_mode = twopi/period
    PRINT*,'Mode#,',modex, modey,' w = ',w_mode
	
	

END SUBROUTINE compute_frequency


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
  
    filename1 = trim('Frames/Fields/fields_t=')//TRIM(uniti)//TRIM('.txt')
    filename2 = trim('Frames/Fields/E_t=')//TRIM(uniti)//TRIM('.txt')
    
    !print*,'filename=',filename

    IF(fields_out) THEN
       OPEN(UNIT=11,FILE=filename1)
       OPEN(UNIT=12,FILE=filename2)
    
       ! save rho, phi and E
        DO j = 1, ny
            DO i = 1, nx
                WRITE(11,*) i*dx,j*dy,rho(i,j),Bz(i,j),Ex(i,j), &
                Ey(i,j),Jx(i,j),Jy(i,j) 
            END DO
        END DO
        
        !DO j = 1, ny
        !    DO i = 1, nx
        !        IF(MOD(j,10) .EQ. 0 .AND. MOD(i,10) .EQ. 0) WRITE(12,*) i*dx,j*dy,Ex(i,j),Ey(i,j) 
        !    END DO
        !END DO
        
        CLOSE(UNIT=11)
    
        CLOSE(UNIT=12)
    
        
        DO i = 1, nx
            WRITE(111,*) i*dx,rho(i,ny/2),ne(i,ny/2),Ex(i,ny/2),Ey(i,ny/2)!, ne(i),ni(i) 
        END DO
    
        DO i = 1, ny
            WRITE(112,*) i*dy,Bz(nx/2,i) 
        END DO
        
	END IF	
    
    IF(particles_out) THEN
	! save particle states
	DO j = 1, N
		WRITE(9,*) particles(1,j)%x,particles(1,j)%y, particles(1,j)%ux,particles(1,j)%uy
		IF(ns .EQ. 2) WRITE(10,*) particles(2,j)%x,particles(2,j)%y, particles(2,j)%ux,particles(2,j)%uy
	END DO
	
	!WRITE(16,*) particles(1,N/2)%x,particles(1,N/2)%y, particles(2,N/2)%x,particles(2,N/2)%y

	END IF
	
	
	
	
	IF(fieldsk_out) THEN
	! save Re[rhok] and Re[phik]
    !DO k = 1, nx-1
		!WRITE(12,*) k*twopi/L, rhok(k,1),phik(k,1),ESEk(k) 
	!END DO
	END IF
	
	IF(energy_out) THEN
		IF(i_sim .GT. 0) WRITE(13,*) i_sim, KE, ESE2, KE+ESE2, Px, vd1, vspread1, vd2, vspread2 
	END IF

	IF(modes_out) THEN
	!IF(i_sim .GT. 0) WRITE(14,*) i_sim, ESEk(1),ESEk(2), ESEk(4), ESEK(6),ESEk(8), ESEk(10), ESEk(12), ESEk(14) ,ESEk(16) 
	END IF
	
	IF(fv_out)THEN
		dv = (vmax-vmin)/vbins 
		!DO i = 1, vbins
		!	WRITE(17,*) vmin+(i-0.5)*dv, fv(1,i), fv(2,i)
		!END DO
	END IF
	
END SUBROUTINE output

END PROGRAM driver_plasma_2d
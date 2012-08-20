!###########################################################################
!
!     ##################################################################
!     ######                                                      ######
!     ######                SUBROUTINE READ_CM1_RAW               ######
!     ######                                                      ######
!     ##################################################################
!
!     PURPOSE:
!
!     Read raw binary cm1 files
!
!############################################################################
!
!     Author:  Lou Wicker
!
!     Creation Date:  10 April 2012
!
!############################################################################

  SUBROUTINE READ_CM1_RAW(time, u, v, w, nx, ny, nz)

    implicit none
    
! Passed variables
  
!   character(LEN=11), INTENT(IN) :: prefix         ! number of input locations to be interpolated to
    integer, INTENT(IN)          :: time
    integer, INTENT(IN)          :: nx, ny, nz     ! grid dimensions
    real(kind=4),    INTENT(OUT)         :: u(nz,ny,nx+1)
    real(kind=4),    INTENT(OUT)         :: v(nz,ny+1,nx)
    real(kind=4),    INTENT(OUT)         :: w(nz+1,ny,nx)

! Local variables

    character*32 :: filename
    character*4  :: stime
    integer      :: i, j, k, lprefix, irec
    real(kind=4), allocatable, dimension(:,:,:) :: u_loc, v_loc, w_loc

! Read raw binary u-field

    print *, "Allocating arrays of dimensions:  ", nx, ny, nz
    allocate(u_loc(nx+1,ny,nz))
    allocate(v_loc(nx,ny+1,nz))
    allocate(w_loc(nx,ny,nx+1))

!   write(filename,'(a32)') "################################"
    write(stime,'(i4.4)') time

!   lprefix = len(prefix)
!   print *, prefix(1:lprefix)
    lprefix = 20

    filename(1:lprefix) = "raw_data/cm1out_"//stime

    filename(lprefix+1:lprefix+6) = "_u.dat"
    print *, "Now reading in U-file.......", filename(1:lprefix+6)
    open(unit=20,file=filename(1:lprefix+6),form='unformatted',access='direct',recl=((nx+1)*ny),status='old',convert='BIG_ENDIAN')
    irec = 1
    DO k = 1,nz
     READ(20, rec=irec) ((u_loc(i,j,k),i=1,nx+1),j=1,ny)
     irec = irec + 1
    ENDDO
    close(20)

    filename(lprefix+1:lprefix+6) = "_v.dat"
    print *, "Now reading in V-file.......", filename(1:lprefix+6)
    open(unit=20,file=filename(1:lprefix+6),form='unformatted',access='direct',recl=((ny+1)*nx),status='old',convert='BIG_ENDIAN')
    irec = 1
    DO k = 1,nz
     READ(20, rec=irec) ((v_loc(i,j,k),i=1,nx),j=1,ny+1)
     irec = irec + 1
    ENDDO
    close(20)

    filename(lprefix+1:lprefix+6) = "_w.dat"
    print *, "Now reading in W-file.......", filename(1:lprefix+6)
    open(unit=20,file=filename(1:lprefix+6),form='unformatted',access='direct',recl=(ny*nx),status='old',convert='BIG_ENDIAN')
    irec = 1
    DO k = 1,nz+1
     READ(20, rec=irec) ((w_loc(i,j,k),i=1,nx),j=1,ny)
     irec = irec + 1
    ENDDO
    close(20)

    print *, "ALL DATA READ IN"

    DO i = 1,nx
     DO j = 1,ny
      DO k = 1,nz

      if( i == nx ) u(k,j,nx+1) = u_loc(nx+1,j,k)
      if( j == ny ) v(k,ny+1,i) = v_loc(i,ny+1,k)
      if( k == nz ) w(nz+1,j,i) = w_loc(i,j,nz+1)
    
      u(k,j,i) = u_loc(i,j,k)
      v(k,j,i) = v_loc(i,j,k)
      w(k,j,i) = w_loc(i,j,k)

      ENDDO
     ENDDO
    ENDDO
    
   print *, "DATA SWAPPED, returning...."

   deallocate(u_loc, v_loc, w_loc)

  RETURN
  END
!###########################################################################
!
!     ##################################################################
!     ######                                                      ######
!     ######                REAL FUNCTION TLINT                   ######
!     ######                                                      ######
!     ##################################################################
!
!     PURPOSE:
!
!     This function returns the tri-linearly interpolated value of a scalar
!     quantity at the given (y,z,x) location.
!
!############################################################################
!
!     Author:  David Dowell
!
!     Creation Date:  16 March 2002, updated 29 June 2004
!
!     30 Aug. 2007: (erm) Had to split up IF statement to avoid out-of-bounds
!                   references
!     Updated by Mike Coniglio, 07/2008
!
!     Adapted by Lou Wicker in March 2012 for trajectory calculations
!
!############################################################################

  SUBROUTINE TRILINEAR(x, y, z, nloc, v, xc, yc, zc, nx, ny, nz, xmin, ymin, missing, w)

    implicit none
    
! Passed variables
  
    integer, INTENT(IN)    :: nx, ny, nz     ! grid dimensions
    integer, INTENT(IN)    :: nloc           ! number of input locations to be interpolated to
    real,    INTENT(IN)    :: v(nz,ny,nx)    ! data array
    real,    INTENT(IN)    :: x(nloc)        ! location at which to compute interpolated value
    real,    INTENT(IN)    :: y(nloc)        ! location at which to compute interpolated value
    real,    INTENT(IN)    :: z(nloc)        ! location at which to compute interpolated value
    real,    INTENT(IN)    :: xmin, ymin     ! coordinates of southwest corner of grid
    real,    INTENT(IN)    :: xc(nx)         ! coordinates corresponding to grid indices
    real,    INTENT(IN)    :: yc(ny)         ! coordinates corresponding to grid indices
    real,    INTENT(IN)    :: zc(nz)         ! coordinates corresponding to grid indices
    real,    INTENT(IN)    :: missing        ! bad/missing data flag
    real,    INTENT(OUT)   :: w(nloc)        ! field v interpolated to the x/y/z locs
    
! Local variables

    integer i, j, k, n
    integer ip1, jp1, kp1
    real q1, q2, vb, vt
    real dx, dx1, dx2, dy, dy1, dy2, dz1, dz2, dz
    integer find_index
    integer i1, i2, j1, j2, k1, k2
  
! Main loop over the number of input locations

    DO n = 1, nloc
    
      IF ( x(n) == missing .or. y(n) == missing .or. z(n) == missing ) THEN
        w(n) = missing
        return
      ENDIF

      i   = find_index(x(n)-xmin, xc, nx)
      ip1 = min(nx, i + 1)
      
      j = find_index(y(n)-ymin, yc, ny)
      jp1 = min(ny, j + 1)
    
      IF ( nz .gt. 1 ) THEN
      
        k   = find_index(z(n), zc, nz)
        kp1 = min(nz, k + 1)
        
! Special boundary conditions for lower boundary

        IF( k == 0 .and. zc(1) .ne. 0.0 ) THEN
          k   = 1                               ! This is the boundary condition for freeslip tangential velocity at bottom
          kp1 = 1
        ENDIF

! Special boundary conditions for top boundary

        IF( k == nz+1 .and. zc(1) .ne. 0.0 ) THEN
          k   = nz                              ! This is the boundary condition for freeslip tangential velocity at top
          kp1 = nz
        ENDIF
      
      ELSE
        k   = 1
        kp1 = 1
      ENDIF
  
! Check to see if the locations are out of bounds first, set return value to missing

      IF ( (i < 1) .or. (j < 1) .or. (k < 1) .or. (i > nx) .or. (j > ny) .or. (k > nz) ) THEN
        
        w(n) = missing
        
      ELSEIF ( (v(k,j,i)      .eq. missing) .or. (v(k,j,ip1)  .eq. missing) .or. &
               (v(k,jp1,ip1)  .eq. missing) .or. (v(k,jp1,i)  .eq. missing) .or. &
               (v(kp1,j,i)    .eq. missing) .or. (v(kp1,j,ip1).eq. missing) .or. &
               (v(kp1,jp1,ip1).eq. missing) .or. (v(kp1,jp1,i).eq. missing)) THEN

        w(n) = missing
  
      ELSE
    
        IF( i == nx ) THEN            ! This boundary conditions is ZERO GRADIENT
         dx  = 1.0
         dx1 = 1.0
!        i1 = i-1
         i1 = i 
         i2 = i
        ELSE
         dx = xc(i+1)-xc(i)
         dx1 = Min( dx, xc(i+1) - (x(n) - xmin) )
         i1 = i
         i2 = i+1
        ENDIF
        dx2 = dx-dx1
    
        IF( j == ny ) then            ! This boundary conditions is ZERO GRADIENT
         dy  = 1.0
         dy1 = 1.0
!        j1 = j-1
         j1 = j 
         j2 = j
        ELSE
         dy = yc(j+1)-yc(j)
         dy1 = Min( dy, yc(j+1) - ( y(n) - ymin ) )
         j1 = j
         j2 = j+1
        ENDIF
        dy2 = dy-dy1
    
        IF ( nz .gt. 1 ) THEN
        
          IF( k == nz ) THEN          ! This boundary conditions is ZERO GRADIENT
           dz = zc(k)-zc(k-1)
           dz1 = dz
!          k1 = k-1
           k1 = k 
           k2 = k
          ELSE
           dz = zc(k+1)-zc(k)
           dz1 = Min( dz, zc(k+1) - z(n) )
           k1 = k
           k2 = k+1
          ENDIF
        ELSE
          dz  = 1.0
          dz1 = 1.0
          k1  = 1
          k2  = 1
        ENDIF
        dz2 = dz-dz1
          
        q1   = dx1*v(k1,j1,i1) + dx2*v(k1,j1,i2)
        q2   = dx1*v(k1,j2,i1) + dx2*v(k1,j2,i2)
        vb   = (dy1*q1 + dy2*q2) / ( dx*dy )
        q1   = dx1*v(k2,j1,i1) + dx2*v(k2,j1,i2)
        q2   = dx1*v(k2,j2,i1) + dx2*v(k2,j2,i2)
        vt   = (dy1*q1 + dy2*q2) / ( dx*dy )
        w(n) = (dz1*vb + dz2*vt) / dz
    
      ENDIF
    
    ENDDO

  RETURN
  END

!###########################################################################
!
!     ##################################################################
!     ######                                                      ######
!     ######          SUBROUTINE TRAJ_INTEGRATE                   ######
!     ######                                                      ######
!     ##################################################################
!
!     PURPOSE:
!
!
!############################################################################
!
!     Adapted by Lou Wicker in March 2012 for trajectory calculations
!
!############################################################################

  SUBROUTINE TRAJ_INTEGRATE(x, y, z, u0, v0, w0, u1, v1, w1, xc, xe, yc, ye, zc, ze, &
                            dt, t0, t1, xmin, ymin, xnew, ynew, znew, unew, vnew, wnew, nx, ny, nz, nstep, nloc)

    implicit none
    
! Passed variables
  
    integer, INTENT(IN)    :: nx, ny, nz          ! grid dimensions for centered values
    integer, INTENT(IN)    :: nstep               ! number of time steps to take
    integer, INTENT(IN)    :: nloc                ! number of trajectories
    real,    INTENT(IN)    :: dt                  ! trajectory time step
    real,    INTENT(IN)    :: t0                  ! time of u0/v0/w0 variables
    real,    INTENT(IN)    :: t1                  ! time of u1/v1/w1 variables
    real,    INTENT(IN)    :: x(nloc)             ! starting location
    real,    INTENT(IN)    :: y(nloc)             ! starting location
    real,    INTENT(IN)    :: z(nloc)             ! starting location
    real,    INTENT(IN)    :: xmin, ymin          ! coordinates of southwest corner of grid
    
    real,    INTENT(IN)    :: u0(nz,ny,nx+1)      ! data array
    real,    INTENT(IN)    :: u1(nz,ny,nx+1)      ! data array
    real,    INTENT(IN)    :: v0(nz,ny+1,nx)      ! data array
    real,    INTENT(IN)    :: v1(nz,ny+1,nx)      ! data array
    real,    INTENT(IN)    :: w0(nz+1,ny,nx)      ! data array
    real,    INTENT(IN)    :: w1(nz+1,ny,nx)      ! data array

    real,    INTENT(IN)    :: xc(nx)              ! coordinates corresponding to grid indices
    real,    INTENT(IN)    :: yc(ny)              ! coordinates corresponding to grid indices
    real,    INTENT(IN)    :: zc(nz)              ! coordinates corresponding to grid indices
    real,    INTENT(IN)    :: xe(nx+1)            ! coordinates corresponding to grid indices
    real,    INTENT(IN)    :: ye(ny+1)            ! coordinates corresponding to grid indices
    real,    INTENT(IN)    :: ze(nz+1)            ! coordinates corresponding to grid indices
    real,    INTENT(OUT)   :: xnew(nstep,nloc)    ! x-paths of trajectory
    real,    INTENT(OUT)   :: ynew(nstep,nloc)    ! y-paths of trajectory
    real,    INTENT(OUT)   :: znew(nstep,nloc)    ! z-paths of trajectory
    real,    INTENT(OUT)   :: unew(nstep,nloc)    ! x-paths of trajectory
    real,    INTENT(OUT)   :: vnew(nstep,nloc)    ! y-paths of trajectory
    real,    INTENT(OUT)   :: wnew(nstep,nloc)    ! z-paths of trajectory
    
! Local variables

    integer n
    real    xtmp(nloc), ytmp(nloc), ztmp(nloc)
    real    us0(nloc), us1(nloc), us2(nloc), us3(nloc)
    real    vs0(nloc), vs1(nloc), vs2(nloc), vs3(nloc)
    real    ws0(nloc), ws1(nloc), ws2(nloc), ws3(nloc)
    real    velo0(nloc), velo1(nloc)
    real    dt_half, t10, t11, t20, t21, traj_time
    
    real, parameter :: missing = -999.
    
    real, parameter :: z_min = 0.0
     
    xnew(1,:) = x(:)
    ynew(1,:) = y(:)
    znew(1,:) = z(:)
    
    dt_half = 0.5 * dt
    
    DO n = 1,nstep-1
    
      IF( mod(n,10) .eq. 0 ) THEN
        print *, 'Traj_integrate on the ',n,'th step'
      ENDIF
    
      xtmp(:) = xnew(n,:)
      ytmp(:) = ynew(n,:)
      ztmp(:) = znew(n,:)
      
      traj_time = t0 + n * dt   ! t0 is the START of trajectory, for backward integration, this means t0 > t1 & dt < 0.
      
      IF (t0 == t1 ) THEN
        t10 = 1.0
        t11 = 0.0
        t20 = 1.0
        t21 = 0.0
      ELSE
        t10 = (t1 - traj_time - dt_half) / (t1 - t0)
        t11 = 1.0 - t10
        t20 = (t1 - traj_time - dt) / (t1 - t0)
        t21 = 1.0 - t20
      ENDIF
 
! Now find velocities at x/y/z and t = n

      CALL TRILINEAR(xtmp, ytmp, ztmp, nloc, u0, xe, yc, zc, nx+1, ny, nz, xmin, ymin, missing, velo0)    
      us0(:) = velo0
      
      CALL TRILINEAR(xtmp, ytmp, ztmp, nloc, v0, xc, ye, zc, nx, ny+1, nz, xmin, ymin, missing, velo0)    
      vs0(:) = velo0

      CALL TRILINEAR(xtmp, ytmp, ztmp, nloc, w0, xc, yc, ze, nx, ny, nz+1, xmin, ymin, missing, velo0)    
      ws0(:) = velo0
      
      IF( n == 1 ) THEN
        unew(n,:) = us0
        vnew(n,:) = vs0
        wnew(n,:) = ws0
      ENDIF

! Update to n+1/2

      xtmp = xnew(n,:) + dt_half*us0
      ytmp = ynew(n,:) + dt_half*vs0
      
      WHERE( ztmp .lt. ze(2) )                    ! Here we use the special expodential decay toward the ground if zpos < ze(2)
                                                  ! This boundary condition will NEVER allow a particle to go below ground                                 
        ztmp = znew(n,:) * exp(dt_half*ws0/znew(n,:))
      ELSEWHERE
        ztmp = znew(n,:) + dt_half*ws0
      END WHERE
  
! Now find velocities at x/y/z and t = n+1/2
! Now we have to interpolate in time as well

      CALL TRILINEAR(xtmp, ytmp, ztmp, nloc, u0, xe, yc, zc, nx+1, ny, nz, xmin, ymin, missing, velo0)    
      CALL TRILINEAR(xtmp, ytmp, ztmp, nloc, u1, xe, yc, zc, nx+1, ny, nz, xmin, ymin, missing, velo1)    
      us1(:) = t10*velo0 + t11*velo1
      
      CALL TRILINEAR(xtmp, ytmp, ztmp, nloc, v0, xc, ye, zc, nx, ny+1, nz, xmin, ymin, missing, velo0)    
      CALL TRILINEAR(xtmp, ytmp, ztmp, nloc, v1, xc, ye, zc, nx, ny+1, nz, xmin, ymin, missing, velo1)    
      vs1(:) = t10*velo0 + t11*velo1

      CALL TRILINEAR(xtmp, ytmp, ztmp, nloc, w0, xc, yc, ze, nx, ny, nz+1, xmin, ymin, missing, velo0)    
      CALL TRILINEAR(xtmp, ytmp, ztmp, nloc, w1, xc, yc, ze, nx, ny, nz+1, xmin, ymin, missing, velo1)    
      ws1(:) = t10*velo0 + t11*velo1
      
! Update to n+1/2*

      xtmp = xnew(n,:) + dt_half*us1
      ytmp = ynew(n,:) + dt_half*vs1
      
      WHERE( ztmp .lt. ze(2) )                    ! Here we use the special expodential decay toward the ground if zpos < ze(2)
                                                  ! This boundary condition will NEVER allow a particle to go below ground                                 
        ztmp = znew(n,:) * exp(dt_half*ws1/znew(n,:))
      ELSEWHERE
        ztmp = znew(n,:) + dt_half*ws1
      END WHERE

! Now find velocities at x/y/z and t = n+1/2*
! Now we have to interpolate in time as well

      CALL TRILINEAR(xtmp, ytmp, ztmp, nloc, u0, xe, yc, zc, nx+1, ny, nz, xmin, ymin, missing, velo0)    
      CALL TRILINEAR(xtmp, ytmp, ztmp, nloc, u1, xe, yc, zc, nx+1, ny, nz, xmin, ymin, missing, velo1)    
      us2(:) = t10*velo0 + t11*velo1
      
      CALL TRILINEAR(xtmp, ytmp, ztmp, nloc, v0, xc, ye, zc, nx, ny+1, nz, xmin, ymin, missing, velo0)    
      CALL TRILINEAR(xtmp, ytmp, ztmp, nloc, v1, xc, ye, zc, nx, ny+1, nz, xmin, ymin, missing, velo1)    
      vs2(:) = t10*velo0 + t11*velo1

      CALL TRILINEAR(xtmp, ytmp, ztmp, nloc, w0, xc, yc, ze, nx, ny, nz+1, xmin, ymin, missing, velo0)    
      CALL TRILINEAR(xtmp, ytmp, ztmp, nloc, w1, xc, yc, ze, nx, ny, nz+1, xmin, ymin, missing, velo1)    
      ws2(:) = t10*velo0 + t11*velo1
      
! Update to n+1*

      xtmp = xnew(n,:) + dt*us2
      ytmp = ynew(n,:) + dt*vs2
      
      WHERE( ztmp .lt. ze(2) )                    ! Here we use the special expodential decay toward the ground if zpos < ze(2)
                                                  ! This boundary condition will NEVER allow a particle to go below ground                                 
        ztmp = znew(n,:) * exp(dt*ws2/znew(n,:))
      ELSEWHERE
        ztmp = znew(n,:) + dt*ws2
      END WHERE

! Now find velocities at x/y/z and t = n+1
! Now we have to interpolate in time as well

      CALL TRILINEAR(xtmp, ytmp, ztmp, nloc, u0, xe, yc, zc, nx+1, ny, nz, xmin, ymin, missing, velo0)    
      CALL TRILINEAR(xtmp, ytmp, ztmp, nloc, u1, xe, yc, zc, nx+1, ny, nz, xmin, ymin, missing, velo1)    
      us3(:) = t20*velo0 + t21*velo1
      
      CALL TRILINEAR(xtmp, ytmp, ztmp, nloc, v0, xc, ye, zc, nx, ny+1, nz, xmin, ymin, missing, velo0)    
      CALL TRILINEAR(xtmp, ytmp, ztmp, nloc, v1, xc, ye, zc, nx, ny+1, nz, xmin, ymin, missing, velo1)    
      vs3(:) = t20*velo0 + t21*velo1

      CALL TRILINEAR(xtmp, ytmp, ztmp, nloc, w0, xc, yc, ze, nx, ny, nz+1, xmin, ymin, missing, velo0)    
      CALL TRILINEAR(xtmp, ytmp, ztmp, nloc, w1, xc, yc, ze, nx, ny, nz+1, xmin, ymin, missing, velo1)    
      ws3(:) = t20*velo0 + t21*velo1   
      
! Update to n+1 using class RK3 formula

      xtmp = xnew(n,:) + dt*( us0 + 2.0*(us1+us2) + us3 ) / 6.0
      ytmp = ynew(n,:) + dt*( vs0 + 2.0*(vs1+vs2) + vs3 ) / 6.0
      
      ws3 = ( ws0 + 2.0*(ws1+ws2) + ws3 ) / 6.0
      
      WHERE( ztmp .lt. ze(2) )                    ! Here we use the special expodential decay toward the ground if zpos < ze(2)
                                                  ! This boundary condition will NEVER allow a particle to go below ground                                 
        ztmp = znew(n,:) * exp(dt*ws3/znew(n,:))
      ELSEWHERE
        ztmp = znew(n,:) + dt*ws3
      END WHERE

! Update location arrays and start over!
      
      xnew(n+1,:) = xtmp
      ynew(n+1,:) = ytmp
      znew(n+1,:) = max(ztmp, z_min)

      unew(n+1,:) = us3
      vnew(n+1,:) = vs3
      wnew(n+1,:) = ws3
      
    ENDDO
      
  RETURN
  END

!###########################################################################
!
!     ##################################################################
!     ######                                                      ######
!     ######             INTEGER FUNCTION FIND_INDEX              ######
!     ######                                                      ######
!     ##################################################################
!
!     PURPOSE:
!
!     This function returns the array index (here, the value returned by
!     find_index is designated as i) such that x is between xa(i) and xa(i+1).
!     If x is less than xa(1), then i=0 is returned.  If x is greater than
!     xa(n), then i=n+1 is returned.  It is assumed that the values of
!     xa increase monotonically with increasing i.
!
!############################################################################
!
!     Author:  David Dowell (based on "locate" algorithm in Numerical Recipes)
!
!     Creation Date:  17 November 2004
!
!     Fixed up a bit by LJW
!
!############################################################################

  INTEGER FUNCTION FIND_INDEX(x, xa, np)

    implicit none

    integer np                          ! array size
    real xa(np)                         ! array of locations
    real x                              ! location of interest
    integer il, im, iu                  ! lower and upper limits, and midpoint

    il = 0
    iu = np+1
    
! Check for coordinate array increasing....

    IF( xa(np) < xa(1) ) THEN
      print *, "FIND_INDEX ERROR, COORDINATE ARRAY MUST BE INCREASING, RET=0"
      print *, "FIND_INDEX:  XC(1):   ", xa(1)
      print *, "FIND_INDEX:  XC(NX):  ", xa(np)
      find_index = 0
      RETURN
    ENDIF
    
    IF ( x > xa(np) ) THEN
      il = np+1
      
    ELSEIF ( x < xa(1) ) THEN
      il = 0
      
    ELSE
    
      DO 
      
        IF( (iu - il) < 2 ) EXIT      ! This tests to see if we have found the point
        
        im = (il + iu)/2
        
        IF( x >= xa(im) ) THEN
          il = im
        ELSE
           iu=im
        ENDIF
        
      ENDDO
  
    ENDIF

    find_index = il

  RETURN
  END


      program T9
      use constants
      use vectors
      use particles
      use omp_lib
      implicit none

! unit system of the simulation
      real, parameter ::        &
        m_unit = M_sun,         &
        l_unit = AU,            &
        t_unit = year,          &
        v_unit = l_unit/t_unit, &
        a_unit = v_unit/t_unit, &
        G = G_cgs * m_unit * t_unit**2 / l_unit**3

      real, parameter :: dtmin = 1e-10 * year/t_unit
      real :: t,dt,dt0,dtout,tend,tout
      real :: rad,rho,mass,dx,tff,grid
      real :: rs,rs2,accmax,tcpu0,tcpu1
      real, dimension(3) :: pos
      integer :: i,j,k,l,it,io,np,num
      type(particle), dimension(:), allocatable :: sys
      character(len=80) :: filename
      logical :: var_dt

! statement function for the initial grid positions
      grid(i) = real(2*i-num-1)*dx

! input job parameters
      read(*,*) tend
      read(*,*) dtout
      read(*,*) dt0
      read(*,*) var_dt
      read(*,*) rs
      read(*,*) rad
      read(*,*) rho
      read(*,*) num

      tend  = tend  * year  / t_unit
      dtout = dtout * year  / t_unit
      rad   = rad   * AU    / l_unit
      rho   = rho   * M_sun/AU**3 / (m_unit/l_unit**3)
      if(var_dt) then
        dt0 = dt0   * year  / t_unit * AU/year**2 / a_unit
      else
        dt  = dt0   * year  / t_unit
      endif

      tff = sqrt(3.*pi/(32.*g*rho))
      write(0,'("Free-fall time = ",es12.4," yr")') tff * t_unit/year
      dx = rad/real(num)
      rs2 = (rs*dx)**2

! count the particles and set up the grid
      do l=0,1
        np = 0
        do i=1,num
          do j=1,num
            do k=1,num
              pos = [grid(i),grid(j),grid(k)]
              if(veclen(pos).lt.rad) then
                if(l.eq.0) then
                  np = np + 1
                else
                  np = np + 1
                  sys(np)%mass = mass
                  sys(np)%pos  = pos
                  sys(np)%vel  = [0.,0.,0.]
                endif
              endif
            enddo
          enddo
        enddo
        if(l.eq.0) then
          allocate(sys(np))
          write(0,'(i0," particles total")') np
!         mass = 4.*pi/3.*rad**3*rho/real(np)
          mass = 8.*rho*dx**3
          write(0,'("particle mass = ",es12.4," Msun")') mass * m_unit/M_sun
        endif
      enddo

! say how many threads we're going to use
      !$omp parallel
      !$omp master
      write(0,'("Using ",i0," OpenMP threads.")') omp_get_num_threads()
      !$omp end master
      !$omp end parallel

! format for output file names
    2 format("T9-",i5.5,".out")

! output format is one line for time and for each particle one line of 6 values
!   3 format("# ",f12.6/(9es16.8))
    4 format("# ",f12.6/(6es16.8))

      t    = 0.
      it   = 0
      io   = 1
      tout = io*dtout
!     tout = t+dtout
      dt   = dt0
      call cpu_time(tcpu0)

      call accel(sys,np,G,rs2)
      write(0,'(es12.4,2i10)') t,it,0

! output initial conditions
      write(filename,2) 0
      open(1,file=filename)
!     write(1,3) t,(sys(i)%pos,sys(i)%vel,sys(i)%mass,sys(i)%Epot,sys(i)%Ekin,i=1,np)
      write(1,4) t,(sys(i)%pos,sys(i)%vel,i=1,np)
      close(1)

! main loop
      do while(t.lt.tend)

! calculate next timestep size (we already have acc from previous step)
        it = it + 1
        if(var_dt) then
          accmax = 0.
          do i=1,np
            accmax = max(accmax,veclen2(sys(i)%acc))
          enddo
          dt = dt0/sqrt(accmax)
!         if(dtout.gt.0.) dt = min(dt,dtout)
          t  = t + dt
        else
          t = it*dt
        endif
        if(dt.lt.dtmin) then
          write(0,'("step size underflow, dt = ",es12.4," at t = ",es12.4)') dt,t
          stop 1
        endif

! kick-drift-kick leapfrog
        do i=1,np
          sys(i)%vel  = sys(i)%vel + sys(i)%acc*dt*.5 ! kick
          sys(i)%pos  = sys(i)%pos + sys(i)%vel*dt    ! drift
        enddo
        call accel(sys,np,G,rs2)
        do i=1,np
          sys(i)%vel  = sys(i)%vel + sys(i)%acc*dt*.5 ! kick
!         sys(i)%Ekin = .5*sys(i)%mass*veclen2(sys(i)%vel)
        enddo

! output only every dtout time-units
        if(t.ge.tout) then
          write(0,'(es12.4,2i10)') t,it,io
          write(filename,2) io
          open(1,file=filename)
!         write(1,3) t,(sys(i)%pos,sys(i)%vel,sys(i)%mass,sys(i)%Epot,sys(i)%Ekin,i=1,np)
          write(1,4) t,(sys(i)%pos,sys(i)%vel,i=1,np)
          close(1)
          io   = io + 1
          tout = io*dtout
!         tout = t+dtout
        endif

! end of main loop
      enddo
      call cpu_time(tcpu1)
      write(0,1) tcpu1-tcpu0,it,io,t/it
    1 format("# cpu-time ",f0.2," seconds"/ &
        "# ",i0," steps (",i0," output)"/ &
        "# average step size ",es11.4)

      end


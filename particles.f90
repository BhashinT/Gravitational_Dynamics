
      module particles
      use vectors
      implicit none

      type particle
        real, dimension(3) :: pos, vel, acc
        real :: mass, Epot, Ekin
      end type

      contains

      subroutine accel(sys,np,G,rs2)
      implicit none

      integer :: np,i,j
      type(particle), dimension(np) :: sys
      real, dimension(3) :: dr,acc
      real :: G,rs2,dr1,dr3,Epot

!$omp parallel do private(i,j,dr,dr1,dr3,acc,Epot)
      do i=1,np
        acc  = 0.
!       Epot = 0.
        do j=1,np
          if(i.ne.j) then
            dr = sys(j)%pos-sys(i)%pos
            dr1 = sqrt(veclen2(dr)+rs2)
            dr3 = dr1**3
            acc  = acc  + sys(j)%mass*dr/dr3
!           Epot = Epot + sys(j)%mass/dr1     *sys(j)%mass/(sys(i)%mass+sys(j)%mass)
          endif
        enddo
        sys(i)%acc  =  G * acc
!       sys(i)%Epot = -G * Epot * sys(i)%mass
      enddo
!$omp end parallel do

      end subroutine

      end module


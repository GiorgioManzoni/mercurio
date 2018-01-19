      program kepler
        implicit none
        real :: x,v_x,y,v_y,dt      
        print*, "KEPLER"
        call input(x,v_x,y,v_y,dt)
        call calculate(x,v_x,y,v_y,dt)
      end program kepler
!---------------------------------------------------
      subroutine input(x,v_x,y,v_y,dt)
        real,intent(inout) :: x,v_x,y,v_y,dt
        write(*,*) "Insert the initial x position: "
        read(*,*) x
        write(*,*) "Insert the initial y position: "
        read(*,*) y
        write(*,*) "Insert the initial x velocity: "
        read(*,*) v_x
        write(*,*) "Insert the initial y velocity: "
        read(*,*) v_y
        write(*,*) "Insert the time step: "
        read(*,*) dt
        write(*,*) "-------------------------------"
        write(*,*) "INITIAL CONDITIONS: "
        write(*,*) "x   = ",x 
        write(*,*) "y   = ",y
        write(*,*) "v_x = ",v_x
        write(*,*) "v_y = ",v_y
        write(*,*) "dt  = ",dt
        write(*,*) "-------------------------------"
      end subroutine input
!-----------------------------------------------------
      subroutine calculate(x,v_x,y,v_y,dt)
        real,intent(inout) :: x,v_x,y,v_y,dt
        integer :: i=0
        real :: r 
        real, parameter :: Pi = 3.1415926535897932
        open(unit=1,file="kepler_positions.dat")
        write(1,*) "# x    y    v_x    v_y"
        do while (i<100000)
            !radius
            r = sqrt(x**2.+y**2.)
            !first update velocities
            v_x = v_x - (4.*Pi**2.*x*dt)/r**3.
            v_y = v_y - (4.*Pi**2.*y*dt)/r**3.
            !and then update the positions
            x = x + v_x*dt
            y = y + v_y*dt                        
            print*, x, y, v_x, v_y
            write(1,*) x, y, v_x, v_y
            i = i + 1
        end do
        close(1)
      end subroutine calculate
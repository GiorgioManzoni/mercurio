      program mercury
        !VARIABLES      
        implicit none
        double precision :: x,v_x,y,v_y,dt
        double precision :: x0,vx0,y0,vy0,dt0
        double precision,dimension(100) :: alpha
        double precision :: m,a
        double precision :: b,c
        double precision :: precession
        integer :: i
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
        do i=1,size(alpha)
            alpha(i:i) = 0.00005 + i*0.00001
            end do
        
        !alpha(100:100) = 0.005
        
        !alpha(1:1) =   0.0008
        !alpha(2:2) =   0.0007
        !alpha(3:3) =   0.0006    
        !alpha(4:4) =   0.0005
        !alpha(5:5) =   0.0004
        !alpha(6:6) =   0.00009  
        !alpha(7:7) =   0.001  
        !alpha(8:8) =   0.00005
        !alpha(9:9) =   0.0015
        !alpha(10:10) = 0.004
        !alpha(11:11) = 0.003
        !alpha(12:12) = 0.002
        !alpha(13:13) = 0.005
        
        
        print*, alpha(1)
        print*, alpha(2)
        print*, alpha(3)       
        print*, "MERCURY  PRECESSION"
        
        call input_fast_mercury(x,v_x,y,v_y,dt)
        x0 = x
        vx0 = v_x
        y0 = y
        vy0 = v_y
        dt0=dt
        open(unit=4,file="alpha_vs_precession_rate.dat")
        write(4,*) "# results for the linea fit of theta vs time for different alphas"
        write(4,*) "# alpha    m    b"
        write(4,*) "#################################################################"
        !write(4,*) 0.0,0.0, 0.0
        
        do i=1, size(alpha)
            call rehash_initial_condition(x0,vx0,y0,vy0,dt0) 
            !loop of 3500 is ideal to stay in the right angle range
            call calculate(x0,vx0,y0,vy0,dt0,alpha(i),3500)!The exact value for Mercury is 1.1e-8
            call least_square(m,b)
            !   mb(2:2) = least_square()
            print*,"alpha = ",alpha(i)," m = ",m, "b = ",b
            write(4,*) alpha(i),m,b
        end do
        close(4)
        print*,"------------------------------------------------------"
        print*,"------------------------------------------------------"
        call least_square_precession_mercury(a,c,precession)
        print*,"#######################################################"
        print*,"#######################################################"
        print*,"#######################################################"
        print*, "Precession of mercury [degree/years]   = ", a*(1.1e-8), " a = ",a !, "c =  ",c        
        print*, "Precession of mercury [arcsec/century] = ", a*(1.1e-8)*60*60*100
        !be careful you don't want to use the intercept c, because you are sure that
        !if alpha = 0 there is no precession so the major axes keep the same position 
        !forever so that d_theta/d_time = 0 
      end program mercury
!---------------------------------------------------
      subroutine input(x,v_x,y,v_y,dt)
        double precision,intent(inout) :: x,v_x,y,v_y,dt
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
      
      
      !---------------------------------------------------
      subroutine input_fast_mercury(x,v_x,y,v_y,dt)
        double precision,intent(out) :: x,v_x,y,v_y,dt
        x = 0.47
        y = 0.
        v_x = 0.
        v_y = 8.2
        dt = 0.0008
        write(*,*) "-------------------------------"
        write(*,*) "INITIAL CONDITIONS: "
        write(*,*) "x   = ",x 
        write(*,*) "y   = ",y
        write(*,*) "v_x = ",v_x
        write(*,*) "v_y = ",v_y
        write(*,*) "dt  = ",dt
        write(*,*) "-------------------------------"
      end subroutine input_fast_mercury
      
      subroutine rehash_initial_condition(x,v_x,y,v_y,dt)
        double precision,intent(inout) :: x,v_x,y,v_y,dt
        x = 0.47
        y = 0.
        v_x = 0.
        v_y = 8.2
        dt = 0.0008
      end subroutine rehash_initial_condition
!-----------------------------------------------------
      subroutine calculate(x,v_x,y,v_y,dt,alpha,iterations)
        double precision,intent(inout) :: x,v_x,y,v_y,dt
        double precision,intent(in) :: alpha !alpha== 1.1e-8 AU^2 (for Mercury)
        integer,intent(in) :: iterations
        integer :: i
        double precision :: r,t,dr,r_old,theta_prec_rad,theta_prec_deg !,x_temp
        double precision, parameter :: Pi = 3.1415926535897932
        open(unit=1,file="mercury_positions.dat")
        write(1,*) "# Euler-Cromer algorithm to compute the motion of a planet"
        write(1,*) "# considering the relativistic correction taken into account by alpha"
        write(1,*) "# File generated using the following initial conditions: "
        write(1,*) "######################################################################"
        write(1,*) "# x      = ",x
        write(1,*) "# y      = ",y
        write(1,*) "# v_x    = ",v_x
        write(1,*) "# v_y    = ",v_y
        write(1,*) "# dt     = ",dt
        write(1,*) "# alpha  = ",alpha
        write(1,*) "######################################################################"
        write(1,*) "# x    y    v_x    v_y    r     t"    
        
        open(unit=2,file="values_precession.dat")
        write(2,*) "# position and angles of the perihelium in each orbit "
        write(2,*) "# x   y   r   theta[deg]    t"        
        write(2,*) "######################################################"
        
        i=0
        t = 0
        r_old = 0 
        dr = 0
        do while (i<iterations)
            !old radius
            if (i>0) then
                r_old = r
                end if
            !radius
            r = sqrt(x**2.+y**2.)
            if (i>1) then
                dr_old = dr
                end if            
            dr = r - r_old
            
            if ((i>1) .and. ((dr*dr_old)<0.) .and. (dr<0.)) then 
                !if (y>=0.) then
                    theta_prec_rad = ACOS(x/r)
                !else
                !    x_temp = ABS(x)
                !    theta_prec_rad = 180. + ACOS(x_temp/r)
                !    end if
                theta_prec_deg = (theta_prec_rad*180.)/Pi
                write(2,*) x, y, r_old, theta_prec_deg, t
                end if
                
                
            !time from the beginning of the orbit
            t = t + dt
            !first update velocities
            v_x = v_x - ((4.*Pi**2.*x*dt)/r**3.)*(1.+alpha/r**2.)
            v_y = v_y - ((4.*Pi**2.*y*dt)/r**3.)*(1.+alpha/r**2.)
            !and then update the positions
            x = x + v_x*dt
            y = y + v_y*dt                        
            !print*, x, y, v_x, v_y
            write(1,*) x, y, v_x, v_y, r, t
            i = i + 1
            end do
        close(1)
        close(2)
      end subroutine calculate
 
 
!--------------------------------------------------------------------------------------     
      subroutine least_square(m,b)
        character(len=100) :: line
        integer :: IOstatus
        integer :: N
        double precision :: var1,var2,var3,var4,var5
        double precision :: sx,sxy,sx2,sy
        !if I do not use double precision, every time I find a different result for the 
        !least square coefficients
        double precision,intent(out) :: m
        double precision,intent(out) :: b
        double precision :: x_,y_
        open(unit=3,file="values_precession.dat",action="read")
        
        
        m=0
        b=0
        N = 0
        sx = 0
        sxy = 0
        sx2 = 0
        sy = 0
        
        read(3,'(A)') line
        print*,line
        read(3,'(A)') line
        print*,line
        read(3,'(A)') line
        print*,line
        do
            read(3,*,IOSTAT=IOstatus) var1, var2, var3, y_, x_ 
            print*,x_,y_
            
            N = N + 1
            sx = sx + x_
            sxy = sxy + x_*y_
            sx2 = sx2 + x_**2.
            sy = sy +y_
            
            if (IOstatus /= 0) exit
            end do

        print*," "
        print*, N, sx, sxy,sx2,sy

        b = ((sx*sxy-sx2*sy)/(sx**2.-N*sx2))
        m = ((sx*sy-N*sxy)/(sx**2.-N*sx2))

        print*,"m = ", m,"  b = ", b
        
        close(3)
      end subroutine least_square
      
      
      
      
!--------------------------------------------------------------------------------------     
      subroutine least_square_precession_mercury(a,c,precession)
        character(len=100) :: line
        integer :: IOstatus
        integer :: N
        double precision :: var1,var2,var3,var4,var5
        double precision :: sx,sxy,sx2,sy
        !if I do not use double precision, every time I find a different result for the 
        !least square coefficients
        double precision,intent(out) :: a
        double precision,intent(out) :: c
        double precision,intent(out) :: precession
        double precision :: xx,yy
        open(unit=5,file="alpha_vs_precession_rate.dat",action="read")
        
        
        m=0
        b=0
        N = 0
        sx = 0
        sxy = 0
        sx2 = 0
        sy = 0
        
        read(5,'(A)') line
        print*,line
        read(5,'(A)') line
        print*,line
        read(5,'(A)') line
        print*,line
        do
            read(5,*,IOSTAT=IOstatus) xx,yy,var3
            print*,xx,yy
            
            N = N + 1
            sx = sx + xx 
            sxy = sxy + xx*yy
            sx2 = sx2 + xx**2.
            sy = sy +yy
            
            if (IOstatus /= 0) exit
            end do

        print*," "
        print*, N, sx, sxy,sx2,sy
        print*," "

        a = ((sx*sy-N*sxy)/(sx**2.-N*sx2))
        c = ((sx*sxy-sx2*sy)/(sx**2.-N*sx2))

        print*,"a = ", a,"  c = ", c
        
        precession = a*(1.1e-8) !+c
        !be careful not to add the intercept that you are sure it is zero!!!

        close(5)
      end subroutine least_square_precession_mercury
            
      
module duffing
    implicit none

    integer, parameter :: N=100000

    real(8) :: a, b, h, yi, ui

    real(8), dimension(N) :: x,yy,u
    real(8), dimension(N,2) :: y ! Here y(1) = y e y(2) = dy/dx

    contains

    subroutine rk4()
        integer :: i
        real(8), dimension(2) :: k1,k2,k3,k4
        real(8) :: h2,h6

        x(1) = 0
        y(1,1) = yi
        y(1,2) = ui

        h2 = 0.5d0*h
        h6 = h/6.0

        do i=2,N
            k1 = f(y(i-1,:))
            k2 = f(y(i-1,:)+h2*k1)
            k3 = f(y(i-1,:)+h2*k2)
            k4 = f(y(i-1,:) + h*k3)
            
            x(i) = x(i-1) + h
            y(i,:) = y(i-1,:) + h6*(k1 + 2*(k2+k3)+k4)
        enddo

        yy = y(:,1)
        u = y(:,2)
    
    end subroutine

    subroutine reset()
        a = 0.
        b = 0.
        h = 0.
        yi = 0.
        ui = 0.
    end subroutine reset

    function f(x)
        real(8), dimension(2) :: f, x
        f(1) = x(2)
        f(2) = -a*x(1) -b*x(1)**3
    end function

end module

module duffing_t

    implicit none

    integer, parameter :: N=1000000

    real(8) :: yi, w, ti,tf

    real(8), dimension(N) :: yy,u,t
    real(8), dimension(N,2) :: y 

    contains

    subroutine rk4df_t()
        integer :: i
        real(8), dimension(2) :: k1,k2,k3,k4
        real(8) :: h,h2,h6

        h = (tf-ti)/N

        t(1) = ti
        y(1,1) = yi
        y(1,2) = 0

        h2 = 0.5d0*h
        h6 = h/6.0

        do i=2,N
            k1 = f2(y(i-1,:), t(i-1))
            k2 = f2(y(i-1,:)+h2*k1, t(i-1)+h2)
            k3 = f2(y(i-1,:)+h2*k2, t(i-1)+h2)
            k4 = f2(y(i-1,:) + h*k3, t(i-1)+h)
            
            t(i) = t(i-1) + h
            y(i,:) = y(i-1,:) + h6*(k1 + 2*(k2+k3)+k4)
        enddo

        yy = y(:,1)
        u = y(:,2)
    
    end subroutine

    FUNCTION f2(x,t)
        REAL(8), DIMENSION(2) :: f2, x
        REAL(8)               :: t
        f2(1) = x(2)
        f2(2) = -x(1)-0.1*x(1)**3-0.05*x(2)+0.2*cos(w*t)
        RETURN
    END FUNCTION

end module

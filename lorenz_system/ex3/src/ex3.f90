module lorenz

    implicit none

    real(8), parameter :: P = 3.0d0, b=1.0d0
    
    integer, parameter :: N = 10000

    real(8), dimension(N,3) :: x

    real(8) :: Tf, r

    contains  


    subroutine rk4()

        real(8),dimension(3) :: k1, k2, k3, k4 ! The k's now are 3d vectors
        real(8)              :: dt ! The step size in time
        integer :: i

        dt = Tf/N

        x(1,1) = 0.
        x(1,2) = 1.
        x(1,3) = 0.

        do i=1,n
            k1 = dt*f(x(i,:))
            k2 = dt*f(x(i,:) + (0.5*k1))
            k3 = dt*f(x(i,:) + (0.5*k2))
            k4 = dt*f(x(i,:) + k3)

            x(i+1,:) = x(i,:) + (k1 + 2.0*k2 + 2.0*k3 + k4)/6
        enddo    

    end subroutine rk4

    function f(v)
        real(8), dimension(3) :: f, v
        f(1) = -P*(v(1)-v(2))
        f(2) = -v(2) + r*v(1) - v(1)*v(3)
        f(3) = -b*v(3) + v(1)*v(2)
    end function f


end module

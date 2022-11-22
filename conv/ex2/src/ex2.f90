module euler

    implicit none
    
    integer :: N

    real(8),dimension(:),allocatable :: x(:),y(:)

    real(8) :: xi,xf,yi,h

    contains
    
    subroutine it()
        integer                 :: i
        
        N = int((xf-xi)/h)

        allocate(x(N),y(N))

        x(1) = xi
        y(1) = yi

        ! h = (xf-xi)/N

        do i=2,N
            x(i) = x(i-1) + h
            y(i) = y(i-1) + h*g(x(i-1),y(i-1))
        enddo
    end subroutine it
    
        
    function g(xv,yv)
        real(8) :: g, xv, yv
        g = (yv/(2*xv)) + (xv**2/(2*yv))
    end function g

end module euler

module midpoint

    implicit none
    
    integer :: N

    real(8),dimension(:),allocatable :: x(:),y(:)

    real(8) :: xi,xf,yi,h

    contains
    
    subroutine it()
        integer                 :: i

        real(8)                 :: h2, k1, k2

        N = int((xf-xi)/h)

        allocate(x(N),y(N))

        x(1) = xi
        y(1) = yi

        h2 = 0.5*h

        do i=2,N
            x(i) = x(i-1) + h
            k1 = h*g(x(i-1),y(i-1))
            k2 = h*g(x(i-1)+h2,y(i-1)+0.5*k1)
            y(i) = y(i-1) + k2
        enddo

    end subroutine it
    
        
    function g(xv,yv)
        real(8) :: g, xv, yv
        g = (yv/(2*xv)) + (xv**2/(2*yv))
    end function g

end module

module rk4

    implicit none
    
    integer :: N

    real(8),dimension(:),allocatable :: x(:),y(:)

    real(8) :: xi,xf,yi,h

    contains
    
    subroutine it()  
        integer                 :: i

        real(8)                 :: h2, h6, k1, k2, k3, k4

        N = int((xf-xi)/h)

        allocate(x(N),y(N))

        x(1) = xi
        y(1) = yi

        h2 = 0.5*h
        h6 = h/6.

        do i=2,N
            x(i) = x(i-1) + h
            k1 = g(x(i-1),y(i-1))
            k2 = g(x(i-1)+h2,y(i-1)+h2*k1)                
            k3 = g(x(i-1)+h2,y(i-1)+h2*k2)
            k4 = g(x(i-1)+h,y(i-1)+k3*h)

            y(i) = y(i-1) + h6*(k1+2.*k2+2.*k3+k4)

        enddo


    end subroutine it
    
        
    function g(xv,yv)
        real(8) :: g, xv, yv
        g = (yv/(2*xv)) + (xv**2/(2*yv))
    end function g

end module
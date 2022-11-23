# Comparison Between the Convergence of Different Methods

As discussed in the main [README](../../main/README.md) in this example we want to compare the three methods presented before: Euler; Mid-point and Runge-Kutta methods. For that, we will write a program in Fortran with three modules, each one for different methods, we then will run this module in python using f2py. In this README, we will discuss in details how the program was built, but it's important to know more about f2py before.

## f2py

Fortran to Python interface generator (f2py) is a way to execute modules built in Fortran on python. Modules are a sequence of instructions that can change the statement declared but not reproduce an output as in a program. You can think of this like a function in python that returns nothing. The structure of a module is:

```fortran
module name     
   [statement declarations]  
   contains
   [subroutine and function definitions]
end module
```

The statements in a module are accessible (modifiable) through python, but the statements declared on subroutines and functions are not. For example:

```fortran
module sum_x     
   implicit none ! this is important to deactivate the default statement declarations
   integer :: n ! declares n as an integer
   real(8) :: x ! declares x as an double-width real number 
   contains
   subroutine sum()
      integer :: i, j
      j = 1
      do i=1,n
         x = x + j
      enddo
   end subroutine sum
end module
```

Then n and x statements are modifiable, but we can not modify j for example, this is why we need to declare the j value in the subroutine.

The quick way to wrap this module using f2py is:

```terminal
    $ f2py -c sum_x.f90 -m sum_x
```
This command compiles and wraps sum_x.f90 (-c) to create the extension module sum_x.so (-m) in the current directory. We then can run this extension module in python as:

```python3
    from sum_x import sum_x # Importing the sum_x module in the extension module sum_x
    import rich
    
    rich.inspect(sum_x) # Inspecting the sum_x module
```
The inspection gives us the information about the module

```python3
╭──── <fortran object> ─────╮
│ def (...)                 │
│                           │
│ n : 'i'-scalar            │
│ x : 'd'-scalar            │
│ sum()                     │
│                           │
│ n = array(0, dtype=int32) │
│ x = array(0.)             │
╰───────────────────────────╯
```

So the inspection shows us that this module contains an integer named n (n: 'I'-scalar), a double-width real number named x (x: 'd'-scalar) and a subroutine named sum(). As we don't declare the value of n and x they are defined as 0.

Now we can declare the values of x and n as we want. For example

```python3
    sum_x.n = 10
    
    sum_x.x = 1
    
    rich.inspect(sum_x)
```

And the inspection shows us the new values of x and n declared.

```python3
╭───── <fortran object> ─────╮
│ def (...)                  │
│                            │
│ n : 'i'-scalar             │
│ x : 'd'-scalar             │
│ sum()                      │
│                            │
│ n = array(10, dtype=int32) │
│ x = array(1.)              │
╰────────────────────────────╯
```

Lets execute the sum() subroutine and see the results


```python3
    sum_x.sum()
    
    rich.inspect(sum_x)
```
Which  returns

```python3
╭───── <fortran object> ─────╮
│ def (...)                  │
│                            │
│ n : 'i'-scalar             │
│ x : 'd'-scalar             │
│ sum()                      │
│                            │
│ n = array(10, dtype=int32) │
│ x = array(11.)             │
╰────────────────────────────╯
```

As as we expected, that is, the subroutine added j (j=1) n times to x.

We just scratched the surface of what f2py can do, to know more access the f2py [user manual](https://numpy.org/doc/stable/f2py/).

## Building Modules in Fortran

Now that we know more about f2py we can build the modules in Fortran. The three modules have something in common: An x and y vector which will contain the solution, the initial and final of x (xi and xf respectively), the initial value of y (yi), the step size (h) and the function we want to compute $g (x, y) = \frac{y}{2x} + \frac{x^2}{2y}$. So a generic module looks like this:

```fortran
module generic_module
    implicit none
    
    integer :: N

    real(8),dimension(:),allocatable :: x(:),y(:)

    real(8) :: xi,xf,yi,h

    contains
    [A subroutine that compute the solution]
    
    function g(xv,yv)
        real(8) :: g, xv, yv
        g = (yv/(2*xv)) + (xv**2/(2*yv))
    end function g
end module
```

Note that at first we don't know the dimension of the vectors x and y, so they were declared as allocatable (this tool are not available in Fortran 77), this means that we will allocate them later, once we had declare the values of xi, xf and h. 

The subroutine that computes the solution are different, as we discussed in the main [README](../../main/README.md). For the Euler method we have:

```fortran
subroutine it()
    integer                 :: i

    N = int((xf-xi)/h)

    allocate(x(N),y(N))

    x(1) = xi
    y(1) = yi

    do i=2,N
        x(i) = x(i-1) + h
        y(i) = y(i-1) + h*g(x(i-1),y(i-1))
    enddo
end subroutine it
```

For the mid-point method:

```fortran
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
```

And for the Runge-Kutta method:

```fortran
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
```

## Analysing The Results on Python

As discussed above, we need to run the program using f2py to create an extension module:

```
$ f2py -c src/ex2.f90 -m ex2
```

We then import the modules from the extension module created along with some useful libraries:

```python3
from ex2 import euler
from ex2 import midpoint as mp
from ex2 import rk4
import numpy as np
import rich
```

As we want to compare the analytical solution to the numerical one, we need to define a function that returns the analytical solution:

```python3
def sol(x):
    y = np.sqrt((x**3 + 7*x)/2)
    return y
```

Then for each method we will define a function that set xi, xf and yi as parameters but let h as a argument.

```python3
def comp_euler(arg1):

    euler.xi = 1.
    euler.xf = 20.
    euler.yi = 2.
    euler.h = arg1
    euler.it()
    xeuler = np.zeros(euler.n)
    yeuler = np.zeros(euler.n)
    xeuler[:] = euler.x[:]
    yeuler[:] = euler.y[:]

    yeuler_sol = sol(xeuler)

    delta_euler = abs(yeuler - yeuler_sol) 

    mean_delta_euler = np.mean(delta_euler)

    return mean_delta_euler

def comp_mp(arg1):

    mp.xi = 1.
    mp.xf = 20.
    mp.yi = 2.
    mp.h = arg1

    mp.it()
    xmp = np.zeros(mp.n)
    ymp = np.zeros(mp.n)
    xmp[:] = mp.x[:]
    ymp[:] = mp.y[:]

    ymp_sol = sol(xmp)

    delta_mp = abs(ymp - ymp_sol)

    mean_delta_mp = np.mean(delta_mp)

    return mean_delta_mp

def comp_rk4(arg1):

    rk4.xi = 1.
    rk4.xf = 20.
    rk4.yi = 2.
    rk4.h = arg1

    rk4.it()
    xrk4 = np.zeros(rk4.n)
    yrk4 = np.zeros(rk4.n)
    xrk4[:] = rk4.x[:]
    yrk4[:] = rk4.y[:]

    yrk4_sol = sol(xrk4)

    delta_rk4 = abs(yrk4 - yrk4_sol) 

    mean_delta_rk4 = np.mean(delta_rk4)

    return mean_delta_rk4
```

This function returns the mean value of the difference between the numerical solution and the analytical one, so we can run it for different values of h and see their results.

Notice that we could build this program in a much more efficient way, defining xi, xf and yi as parameter and letting h as a statement, we could also build just one module with different subroutines, one for each method. But we build this way to show that we can run several modules in just one extension module.

# The Duffing Oscillator

As discussed before, to solve the Duffing Oscillator, where their equation of motion is given by a second order nonlinear differential equation $\ddot y + \alpha y + \beta y^3$, first we will construct a system of two first order differential equations:

<p align="center">
$u = \dot y$,
</p>
  
<p align="center">
$\dot u = -\alpha y - \beta y^3$
</p>

And for a Damped Duffing Oscillator, we have $\ddot y + \kappa \dot y + \alpha y + \beta y^3 = f \cos{\omega t}$, which gives the following system:

<p align="center">
$u = \dot y$
</p>

<p align="center">
$\dot u = - \alpha y - \beta y^3 - \kappa u + f \cos{\omega t} $,
</p>

As we did in the Lorenz System (see the [README](../lorenz_system/README.md)) we will implement the Runge-Kutta method in Fortran and analyse the results in Python.

## Build the Fortran Modules

As we have two different problems (a non-damped and an damped one) we will construct two modules, named duffing and duffing_t respectively

### The non-damped problem

First, we will define the variables

```fortran
implicit none

integer, parameter :: N=100000

real(8) :: a, b, h, yi, ui

real(8), dimension(N) :: x,yy,u
real(8), dimension(N,2) :: y ! Here y(1) = y e y(2) = dy/dx
```

Here, N is the number of steps, $a = \alpha$, $b = \beta$, yi is the initial condition of $y$, ui is the initial condition of $u=\dot y$  and h is the step size. y is a matrix of $N \times 2$ which contains the solution of the system in each step, yy = y(:,1), i.e. the first column of the system, u = y(:,2) i.e. the second column of the system and x is the variable of integration, remembering that $\dot y =\frac{dy}{dx}$, which is computed at each step, this is why it is a vector with size N.

So we compute the solutions through the rk4() subroutine:

```fortran
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
```
Where:

```fortran
function f(x)
    real(8), dimension(2) :: f, x
    f(1) = x(2)
    f(2) = -a*x(1) -b*x(1)**3
end function
```

Is also interesting to define a subfuction to reset the variables:

```fortran
subroutine reset()
    a = 0.
    b = 0.
    h = 0.
    yi = 0.
    ui = 0.
end subroutine reset
```


### The damped problem

Here we have some differences, we set $\alpha = 1$, $\beta = 0.1$, $\kappa = 0.05$, $f = 0.2$, ui = 0 and let $\omega$ as a variable, so we have this:


```fortran
implicit none

integer, parameter :: N=1000000

real(8) :: yi, w, ti,tf

real(8), dimension(N) :: yy,u,t
real(8), dimension(N,2) :: y 
```
Notice that instead of declaring $h$, the step size, we declare ti and tf, wich is the initial and final time (we did the change of variables $x \rightarrow t$).

The rk4() subroutine is very similar that we had before, but we compute h now:

```fortran
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
```

And now the system of functions to be solve is:

```fortran
FUNCTION f2(x,t)
    REAL(8), DIMENSION(2) :: f2, x
    REAL(8)               :: t
    f2(1) = x(2)
    f2(2) = -x(1)-0.1*x(1)**3-0.05*x(2)+0.2*cos(w*t)
    RETURN
END FUNCTION
```

## Analysing the Results With Python

First, we imported the modules and some useful libraries

```python3
from ex4 import duffing as df
from ex4 import duffing_t as dft

from matplotlib import pyplot as plt
import numpy as np
```

To analyse the orbits, i.e. $\dot y \times y$ for different values of $\beta$ we will build a function that sets the others variables and let $\beta$ as an argument:

```python3
def ex4a(arg1):
    df.a = 1
    df.b = arg1
    df.h = 0.001
    df.yi = 2
    df.ui = 0
    df.rk4()
    y = np.zeros(df.n)
    u = np.zeros(df.n)
    y[:] = df.yy[:]
    u[:] = df.u[:]
    df.reset()
    return y,u
```

To investigate the relation between the amplitude (A) and period (P) for different values of $\beta$, we will build a function that defines the variables (letting $\beta$ as an argument) as we did before, however it varies yi from 0 to 100 in 100 points and compute A and P for each one:


```python3
def ex4b(arg1):
    df.a = (2*np.pi)**2
    df.b = arg2
    df.h = 0.001
    df.ui = 0

    yi = np.linspace(0,100,100)

    amp = np.zeros(len(yi))
    T = np.zeros(len(yi))
    for i in range(len(yi)):
        df.yi = yi[i]
        df.rk4()
        y,x = df.yy,df.x
        k = 0

        for w in range(1,len(y)-1):
            if y[w]>y[w+1] and y[w] > y[w-1]:
                k += 1

        A = np.zeros(k,int)

        j=0

        for w in range(1,len(y)-1):
            if y[w]>y[w+1] and y[w] > y[w-1]:
                A[j] = int(w)
                j += 1

        amp[i] = np.mean(y[A])
        T_temp = np.zeros(len(A))
        for z in range(len(A)-1):
            T_temp[z] = x[A[z+1]] - x[A[z]]
        T[i] = np.mean(T_temp)
    return T, amp
```

Now, to investigate the relations between $\omega$ and A for different values of yi for a Damped Duffing Oscillator we will build a function that sets the variables ti and tf and let yi as an argument, then it will compute A for different values of $\omega$ varying from 1/3 to 3 in 100 points. So it returns $\omega$ and A.

```python3
def ex4c(arg1):
    dft.yi = arg1
    dft.ti = 0
    dft.tf = 10000

    W = np.linspace(1/3,3,100)
    amp = np.zeros(100)

    for j in range(len(W)):
        dft.w = W[j]
        dft.rk4df_t()
        t,y = dft.t,dft.yy
        i = t >= 300 # Setting the relaxation time
        t = t[i]
        y = y[i]
        amp[j] = max(y)

    return W, amp
```

Now we just need to plot the results

```python3
y_0,u_0 = ex4a(0)
y_1,u_1 = ex4a(1/10)
y_2,u_2 = ex4a(1/5)
y_3,u_3 = ex4a(1/2)

plt.figure()
plt.title(r'$y \times \dot{y}$ ' 'for ' r'$\alpha = 1$' ' and different values of 'r'$\beta$')
plt.plot(y_0,u_0,label=r'$\beta = 0.$')
plt.plot(y_1,u_1,label=r'$\beta = 1/10$')
plt.plot(y_2,u_2,label=r'$\beta = 1/5$')
plt.plot(y_3,u_3,label=r'$\beta = 1/2$')
plt.xlabel('y')
plt.ylabel(r'$\dot{y}$')
plt.legend()
plt.savefig('do_1.png',dpi=300)

T_0,amp_0 = ex4b((2*np.pi)**2,0)
T_1,amp_1 = ex4b((2*np.pi)**2,0.1)
T_2,amp_2 = ex4b((2*np.pi)**2,0.2)
T_3,amp_3 = ex4b((2*np.pi)**2,0.3)
T_4,amp_4 = ex4b((2*np.pi)**2,0.4)
T_5,amp_5 = ex4b((2*np.pi)**2,0.5)

plt.figure()
plt.title('Amplitude ' r'$\times$' ' Period for ' r'$\alpha=(2\pi)^2$' ' and different values of 'r'$\beta$')
plt.plot(amp_0,T_0,label=r'$\beta = 0.$')
plt.plot(amp_1,T_1,label=r'$\beta = 0.1$')
plt.plot(amp_2,T_2,label=r'$\beta = 0.2$')
plt.plot(amp_3,T_3,label=r'$\beta = 0.3$')
plt.plot(amp_4,T_4,label=r'$\beta = 0.4$')
plt.plot(amp_5,T_5,label=r'$\beta = 0.5$')
plt.ylabel('Period')
plt.xlabel('Amplitude')
plt.legend()
plt.savefig('do_2.png',dpi=300)

W5,amp5 = ex4c(5.)
W0,amp0 = ex4c(0.)

plt.figure()
plt.title(r'$\omega$ ' r'$\times$ ' 'A '   'for ' r'$\dot{y}(0) = 0$ ' 'and different values of y(0)')
plt.plot(W0,amp0,label='y(0) = 0')
plt.plot(W5,amp5,label='y(0) = 5')
plt.xlabel('Oscillation ' r'$\omega$')
plt.ylabel('Amplitude (A)')
plt.legend()
plt.savefig('do_3.png',dpi=300)
```

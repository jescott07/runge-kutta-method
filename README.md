# Computacional Methods to Solve Differential Equations

In this repository we will show how to use different methods to solve differential equations. In the introduction, we will discuss about differential equations and how we can use different methods to solve them. In a second section we show the convergence of different methods and discuss which is more efficient. In the third section we will study the Lorenz System and in the last sections we will discuss about the duffing oscillator.

## Introduction

In a simple equation we have an equality that relates one variable to others. For example, the equation

<p align="center">
$z=x^2+y$,
</p>

that can be represented as a function of $x$ and $y$

<p align="center">
$f(x,y) = x^2 + y$.
</p>

A differential equation for other hand tells how one or more functions, relate to their derivatives, i.e., the rate of change of one variable in relation to another. For example:

<p align="center">
$f(x, y) = \frac{dx}{dy}$,
</p>

is a first order differential equation. An ordinary differential equation (ODE) is an equation ho involves an ordinary, lets say n, order differential equation. For example, a linear ODE with order N is:

<p align="center">
$f_{N}(x)y^{{(N)}}(x)+f_{{N-1}}(x)y^{{(N-1)}}(x)+\ldots +f_{1}(x)y'(x)+f_{0}(x)y(x)=q(x)$,
</p>

where $y^{(N)}$ represents the N's derivative of $y$.

We can rewrite an ODE of order N as a set of N coupled first-order differential equations, i.e.

<p align="center">
$\frac{dy_i(x)}{dx}=f_i(x,y_1,\ldots,y_N), \qquad i=1,\ldots,N$,
</p>


so the problem of solving a N order ODE reduces to solving a system of N first order differential equation. In math we study a set of methods to solve analytically a first order differential equations, i.e., find a solution $f(x)$ in the example above. But if we have a more complicate function, or a set of functions they will not have an analytical solution, but we can build a computational one. In this section we will show three methods to computational integrate a first order differential equation, the Euler method, the mid-point method and the fourth order Runge-Kutta method. Once we had learned this methods we will show some examples and how to apply the different methods to them.

### The Euler Method

Before discussing the Euler Method, we have to discretize the solution space, i.e., if we want to find a solution $f(x)$ given some first order differential equation, we need to discretize x, for example:

<p align="center">
$x_{n+1} = x_n + h$,
</p>

where, h is the step size of the solution and $x_n$ goes from the initial value $x_i$ to the final value $x_f = x_i + N \cdot h$, where $N$ is the number of iterations.

Thus, given the following equation:

<p align="center">
$f (x, y) = \frac{dy}{dx}$,
</p>

We can discretize it as follows:

<p align="center">
$f(x_n,y_n) = \frac{\Delta y}{\Delta x} = \frac{y_{n+1} - y_n}{x_{n+1} - x_n} = \frac{y_{n+1} - y_n}{(x_n + h) - x_n}$,
</p>

<p align="center">
$\Rightarrow y_{n+1} = y_n + h \cdot f(x_n,y_n)$.
</p>

Thus, knowing $f(x, y)$ and $y_i$, i.e., the initial condition, we can advance the solution to obtain the complete solution of the differential equation. Notice that the expansion that we did before is a first order expansion in power series, thus the error of this method is $\mathcal{O}(h^2)$.

Even though this method is not used in practice, didactically it is important to understand other more complex methods.

One way to implement this method is through the pseudo code:

**Algorithm** *Euler*

**Input** $f(x, y) = \frac{dy}{dx}$: First order differential equation to be integrated (function), $x(0)$, $y(0)$: initial conditions (float), $h$: Integration step (float), $N$: Number of integration steps (positive integer).

**Output** $x(N)$, $y(N)$: Vectors with size N containing the solutions of the first order differential equation.

**Compute y(i+1), x(i+1) for each x(i) and y(i) from $f(x, y)$ and $h$ as discussed above.**

1. Define i = 0.
2. Define $x(0) = x_i$ and $y(0) = y_i$ (the initial condition inputs)
3. **do**

      y(i+1) = y(i) + f(x(i),y(i))
      
      x(i+1) = x(i) + h
      
      i = i + 1
      
    **while** i $\neq$ N + 1
    
 4. Return x(N) and y(N).

### The Mid Point Method

As we saw in the section above, the Euler method is not used in practice due to the $\mathcal{O}(h^2)$ error. So we can try a "second order" Euler method, i.e. instead of take a big step in $y_{n+1}$ we will take a mid point as a trial step, in practice what we need to do is to expand the solution in second order power series. Thus:

<p align="center">
$y_{n+1} = y_n + h\cdot f(x_n + \frac{1}{2}h,y_n + \frac{1}{2}k_1) + \mathcal{O}(h^3) $,
</p>


where

<p align="center">
$k_1 = h \cdot f(x_n,y_n)$.
</p>

This method is also known as second-order Runge-Kutta. Notice that now the error is $\mathcal{O}(h^3)$ which is way better them we had before. 

One way to implement this method is through the pseudo code:

**Algorithm** *MidPoint*

**Input** $f(x, y) = \frac{dy}{dx}$: First order differential equation to be integrated (function), $x(0)$, $y(0)$: initial conditions (float), $h$: Integration step (float), $N$: Number of integration steps (positive integer).

**Output** $x(N)$, $y(N)$: Vectors with size N containing the solutions of the first order differential equation.

**Compute y(i+1), x(i+1) for each x(i) and y(i) from $f(x, y)$ and $h$ as discussed above.**

1. Define i = 0.
2. Define h2 = 0.5 $\cdot$ h
3. Define $x(0) = x_i$ and $y(0) = y_i$ (the initial condition inputs)
4. **do**

      $k_1$ = h $\cdot$ f (x(i), y(i))
      $k_2$ = h $\cdot$ f (x(i) + h2, y(i) + 0.5 $\cdot k_1$)

      y(i+1) = y(i) + $k_2$
      
      x(i+1) = x(i) + h
      
      i = i + 1
      
    **while** I $\neq$ N + 1
    
 4. Return x(N) and y(N).

### The Runge-Kutta Method

As the content presented above suggests, higher expansion will lead to minor errors, in fact an n order power series expansion leads to $\mathcal{O}(h^{n+1})$ error. There is many ways that we can evaluate $f(x, y)$, the fourth-order Runge-Kutta Method (or simply Runge-Kutta method) evaluate the derivative four times in each step: One in $y_n$, one in $y_ {n+1}$ and twice in the midpoint. The solution is computed as follows:

<p align="center">
$y_{n+1} = y_n + \frac{k_1}{6} + \frac{k_2}{3} + \frac{k_3}{3} + \frac{k_4}{6} + \mathcal{O}(h^5)$
</p>

Where

<p align="center">
$k_4 = h \cdot f(x_n + h,y_n + k_3)$,
</p>

<p align="center">
$k_3 = h \cdot f(x_n + \frac{h}{2},y_n + \frac{k_2}{2})$,
</p>

<p align="center">
$k_2 = h \cdot f(x_n + \frac{h}{2},y_n + \frac{k_1}{2})$
</p>

and

<p align="center">
$k_1 = h \cdot f(x_n,y_n)$.
</p>

One way to implement this method is through the pseudo code:

**Algorithm** *rk4*

**Input** $f(x, y) = \frac{dy}{dx}$: First order differential equation to be integrated (function), $x(0)$, $y(0)$: initial conditions (float), $h$: Integration step (float), $N$: Number of integration steps (positive integer).

**Output** $x(N)$, $y(N)$: Vectors with size N containing the solutions of the first order differential equation.

**Compute y(i+1), x(i+1) for each x(i) and y(i) from $f(x,y) $ and $h$ as discussed above.**

1. Define i = 0.
2. Define h2 = 0.5 $\cdot$ h
3. Define h6 = h/6
4. Define $x(0) = x_i$ and $y(0) = y_i$ (the initial condition inputs)
5. **do**

      $k_1$ = h $\cdot$ f(x(i), y(i))
      
      $k_2$ = h $\cdot$ f(x(i) + h2, y(i) + 0.5 $\cdot k_1$)
      
      $k_3$ = h $\cdot$ f(x(i) + h2, y(i) + 0.5 $\cdot k_2$)
      
      $k_4$ = h $\cdot$ f(x(i) + h, y(i) + $k_3$)

      y(i+1) = y(i) + h6 $\cdot$ ($k_1 + 2 \cdot (k_2 + k_3) + k_4$)
      
      x(i+1) = x(i) + h
      
      i = i + 1
      
    **while** i $\neq$ N + 1
    
 6. Return x(N) and y(N).

There is much more to discuss about this methods, for further reading I suggest the chapter 16 of the book Numerical Recipier in F77 from WILLIAM H. PRESS, SAUL A. TEUKOLSKY, WILLIAM T. VETTERLING and BRIAN P. FLANNERY.

## Comparison Between the Convergence of Different Methods

In this section we will analyse the convergence of different methods, Euler, mid point and Runge-Kutta. For this, is useful to integrate a differential equation with an analytical solution.

Consider the following differential equation:

<p align="center">
$\frac{dy}{dx} = \frac{y}{2x} + \frac{x^2}{y}$,
</p>

which has the following analytical solution:

<p align="center">
$y(x) = \sqrt{\frac{x^2 + 7x} {2}} $.
</p>

To analyze the convergence of the series, we write a program in $\texttt{Fortran 90}$ with three modules, which compute the solution of the differential equation above using the three different methods discussed, one for each module. We then use $\texttt{f2py}$ from $\texttt{numpy}$ to execute these modules in $\texttt{Pyhon}$. Then we write a code in $\texttt{Python}$ that compares the numerical solution to the analytical one for different step sizes $h$ and give the error $\epsilon$ for each method for a given $h$.

We see that for the Runge-Kutta method, a $h = 10^{-2}$ generates an error of order $\epsilon \sim 10^{-11}$ to achieve the same error we need a $h \sim 10^ {-5} $ for the midpoint method and a $h \sim 10^{-8}$ for the Euler method. Notice that if we fix $x_i$ and $x_f$ (which is the case here) we need more integration steps given a smaller $h$. This analysis shows the robustness of the Runge-Kutta method compared to others methods.

To see in details the codes, please see the README in in the directory called $\texttt{conv}$ in this repository.

## The Lorenz System 

The Lorenz System is a set of three non-lean differential equations, as shown in the system below:

<p align="center">
$\frac{dy_1}{dx} = -P(y_1 - y_2)$
</p>
      
<p align="center">
$\frac{dy_2}{dx} = -y_2 + r y_1 - y_1 y_3$
</p>
      
<p align="center">
$\frac{dy_3}{dx} = -b y_3 + y_1 y_2$
</p>

Where, $P$, $b$ and $r$ are the parameters of this system. This model was first introduced in the study of convection in the atmosphere and it was a breaking point in the chaos studies in deterministic dynamical systems. If we fix $P$ and $b$ the system will determine the character of the solution.

For $y_1 = y_2 = \pm \sqrt{r-1}$, $y_3 = r-1$ and setting $P=3$ and $b=1$ we have:

<p align="center">
$\frac{dy_1}{dx} = - 3(\sqrt{r-1} - \sqrt{r-1})$,
</p>
      
<p align="center">
$\frac{dy_2}{dx} = -\sqrt{r-1} + r\sqrt{r-1} - \sqrt{r-1} (r-1) $,
</p>
      
<p align="center">
$\frac{dy_3}{dx} = -(r-1) + \sqrt{r-1} \sqrt{r-1}$.
</p>

Thus, for this condition we have the stationary solution, i.e.

<p align="center">
$\frac{dy_1}{dx} = 0$,
</p>
      
<p align="center">
$\frac{dy_2}{dx} = 0$,
</p>
      
<p align="center">
$\frac{dy_3}{dx} = 0$.
</p>

Using the theory of linear stability, this solution can be shown to be stable for:

<p align="center">
$r < r_c = \frac{P(P+b+3)}{P-b-1}$,
</p>

for $P=3$ and $b=-1$, $r_c = 21$ what gives the limit of stability for this parameter. 

To solve the system of differential equations given by the Lorenz system, we write a module in $\texttt{Fortran 90} $, using the Runge-Kutta method and execute it in $\texttt{Python}$ as explained in the section above. In order to analyze the convergence of the solutions, we fixed $P = 3$, $b = 1$ and plotted the trajectory in the $y_1 \times y_2$ plane for two different values of $r$, $r = 10$ and $r = 30$ (as show in the figures below). The initial conditions were set as:

<p align="center">
$y_1(0) = 0$,
</p>

<p align="center">
$y_2(0) = 1$,
</p>

<p align="center">
$y_3(0) = 0$,
</p>

The conditions for the integration was $x_i = 0$, $x_f = 100$ and $h = 10^ {-4} $.

To see in details how the programs were constructed, please see the README in the lorenz_system directory in this repository.

![](/../main/lorenz_system/ls_1.png)


![](/../main/lorenz_system/ls_2.png)

## The Duffing Oscillator

The Duffing Oscillator represents a non-harmonic oscillator, their equation of motion is given by:

<p align="center">
$\ddot y + \alpha y + \beta y^3 = 0$,
</p>

where $\ddot y$ represents the second derivative of y and $\alpha$ and $\beta$ are parameters to be analysed, where $\alpha > 0$.

To solve this differential equation, we write a program in $\texttt{Fortran 90} $ and execute their modules (which use the Runge-Kutta method) in $\texttt{Python}$, to more details see the README in duffing_oscilator folder in this repository.

### Analysing the orbits with different parameters

To analyse the orbits in the plane $y \times \dot y$ we define the initial conditions $\dot y (0) = 0$, $y(0) = 2$ and fixed $\alpha = 1$ varying $\beta$ = 1/10, 1/5, 1/2 and 0. Thus we executed the module for this condition, and plotted the results, as show in the figure below.

![](/.. /main/duffing_oscillator/do_1.png)

We can see that for every $\beta$ the system is conservative, i.e. form a closed orbit, as we expected. Their shape goes to an elliptical one as $\beta \rightarrow 0$, that makes sense since in a harmonic oscillator we expect an elliptical orbit, so we can think the $\beta$ as a "non-harmonicity" parameter.

### Period $\times$ Amplitude

To investigate how the period (P) changes with the amplitude (A) we define the initial parameters as before, but with $\alpha = (2\pi)^2$ and $\beta$ varying from 0 to 0.5 in 5 steps. We then execute the module for 100 periods and calculate the amplitude for each one, we did this for each $\beta$ and plotted the results, as show in the figure below.

![](/.. /main/duffing_oscillator/do_2.png)

We can see that, for every $\beta$, P goes to $(2\pi)^2$ for $A \rightarrow 0$, now, for $\beta = 0$ P is constant as we expected since in this case we have a harmonic oscillator. But, for $\beta > 0$ A decays with P and this decay is greater for higher values of $\beta$.

### Damped Duffing Oscillator

In a damped duffing oscillator, the system is subjected to a periodic external force, that is, the equation is given by:

<p align="center">
$\ddot y + \kappa \dot y + \alpha y + \beta y^3 = f \cos{\omega t}$,
</p>

Solving this differential equation is equivalent to solve the following system of equations:

<p align="center">
$u = \dot y$
</p>

<p align="center">
$\dot u = - \alpha y - \beta y^3 - \kappa u - f \cos{\omega t} $,
</p>

To investigate how A changes with $\omega$, as did before, we execute the module (different from the used above, see the README in duffing_oscillator folder) for 1000 values of $\omega$ varying from 1/3 to 3 and fixing the other parameters, $\alpha = 1$, $\beta = 0.1$, $\kappa = 0.05$ and $f=0.02$. We did this for two different initial conditions a) $y(0) = 0$ and $\dot y (0) =0$ b) $y(0) = 5$ and $\dot y (0) = 0$. The figure below shows the comparison between $A \times \omega$

![](/../main/duffing_oscillator/do_3.png)

From this figure we can see that in principle A changes with $\omega$ regardless of initial conditions, but there is a point where A its maximum in some $\omega$,that is, at the resonant frequency of the system ( $\omega_r$ ). We can see also that $\omega_f$ depends on the initial conditions and in both cases A quickly decays for $\omega > \omega_r$. To read more about the Duffing equation [see](https://en.wikipedia.org/wiki/Duffing_equation) mainly about hysteresis that was not fully investigated in this work.

<!--To analyze the convergence of the series we write a program in $\texttt{Fortran 90}$ with three modules, which compute the solution of the diferencial equation above using the three different methods discussed, one for each module. We then use $\texttt{f2py}$ from $\texttt{numpy}$ to execute this modules in $\texttt{Pyhon}$, to do this we run the folowing comand on the terminal $\texttt{f2py -c directory/FileName.f90 -m ProgramName}$ we them import the modules of the $\texttt{Fortran}$ program on $\texttt{Pyhon}$. To know more about $\texttt{f2py}$ acess the [user guide](https://numpy.org/doc/stable/f2py/).-->

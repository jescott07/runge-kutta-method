# Runge-Kutta method

In this repository I will show how to use Runge-Kutta and other methods to solve differential equations.

## Introduction

In a simple equation we have a equality that relates one variable to others. For example, the equation

<p align="center">
$z=x^2+y$
</p>

can be represented as a fuction of $x$ and $y$

<p align="center">
$f(x,y) = x^2 + y$.
</p>

A diferencial equation for other hand tells how one or more fuctions, relate to their derivatives, i.e, the rate of chage of one variable in relation to other. For example:

<p align="center">
$f(x,y) = \frac{dx}{dy}$,
</p>



is a first order diferencial equation. A ordinary differencial equation (ODE) is a equation ho involves a ordinary, lets sey n, order diferencial equation. For exemple, a linear ODE with order N is:

<p align="center">
$f_{N}(x)y^{{(N)}}(x)+f_{{N-1}}(x)y^{{(N-1)}}(x)+\ldots +f_{1}(x)y'(x)+f_{0}(x)y(x)=q(x)$,
</p>

where $y^{(N)}$ representes the n's derivative of $y$.

We can rewritte a ODE of order N as a set of N coupled first-order differential equations, i.e

<p align="center">
$\frac{dy_i(x)}{dx}=f_i(x,y_1,\ldots,y_N), \qquad i=1,\ldots,N$,
</p>


so the problem of solving a N order ODE reduce to solve a sistem of N first order diferencial equation. In math we study a set of methods to solve analyticaly a first order diferencial equations, i.e, find a solution $f(x)$ in the example above. But if we have a more complicate function, or a set of fuctions they will not have a analytical solution, but we can bulid a computacional one. In this section I will show three methods to integrate computacionaly a first order differencial equation, the Euler method, the mid-point method and the fourth order Runge-Kutta method. Once we had learn this methods I will show some examples and how to aply the diferent methods to them.

### The Euler Method

Before discussing the Euler Method, we have to discretize the solution space, i.e, if we want to find a solution $f(x)$ given some first order diferencial equation, we need to discretize x, for example:

<p align="center">
$x_{n+1} = x_n + h$,
</p>

were, h is the step size of the solution and $x_n$ goes from the inicial value $x_i$ to the final value $x_f = x_i + N \cdot h$, where $N$ is the number of iterations.

Thus, given the folowing equation:

<p align="center">
$f(x,y) = \frac{dy}{dx}$,
</p>

we can discretize it as follows:

<p align="center">
$f(x_n,y_n) = \frac{\Delta y}{\Delta x} = \frac{y_{n+1} - y_n}{x_{n+1} - x_n} = \frac{y_{n+1} - y_n}{(x_n + h) - x_n}$,
</p>

Thus

<p align="center">
$y_{n+1} = y_n + h \cdot f(x_n,y_n)$.
</p>

Thus, knowing $f(x,y)$ and $y_i$, ie, the inicial condition, we can advance the solution to obtein the complete solution of the diferencial equation. Notice that the expansion that we did before is a first order expansion in power series, thus the error of this method is $\mathcal{O}(h^2)$.

Even though this method is not used in practice, didactically it is important to understand other more complex methods.

One way to implement this method is through the pseudo code:

**Algorithm** *Euler*

**Input** $f(x,y) = \frac{dy}{dx}$: First order differencial equation to me integrated (function), $x(0)$, $y(0)$: Inicial conditions (float), $h$: Integration step (float), $N$: Number of integrations steps (positive integer).

**Output** $x(N)$, $y(N)$: Vectors with size N contening the solutions of the first order differencial equation.

**Compute y(i+1), x(i+1) for each x(i) and y(i) from $f(x,y)$ and $h$ as discusted above.**

1. Define i = 0.
2. Define $x(0) = x_i$ and $y(0) = y_i$ (the inicial condicions inputs)
3. **do**

      y(i+1) = y(i) + f(x(i),y(i))
      
      x(i+1) = x(i) + h
      
      i = i + 1
      
    **while** i $\neq$ N + 1
    
 4. Return x(N) and y(N).

### The Mid Point Method

As we saw in the section above, the Euler method is not used in practice due to the $\mathcal{O}(h^2)$ error. So we can try a "second order" Euler method, i.e. insted of take a big step in $y_{n+1}$ we will take a mid point as a trial step, in practice what we need to do is to expand the solution in second order power series. Thus:

<p align="center">
$y_{n+1} = y_n + h\cdot f(x_n + \frac{1}{2}h,y_n + \frac{1}{2}k_1) + \mathcal{O}(h^3) $,
</p>


where

<p align="center">
$k_1 = h \cdot f(x_n,y_n)$.
</p>

This method is also know as second-order Runge-Kutta. Notice that now the error is $\mathcal{O}(h^3)$ whic is way beter them we had before. 

One way to implement this method is through the pseudo code:

**Algorithm** *MidPoint*

**Input** $f(x,y) = \frac{dy}{dx}$: First order differencial equation to me integrated (function), $x(0)$, $y(0)$: Inicial conditions (float), $h$: Integration step (float), $N$: Number of integrations steps (positive integer).

**Output** $x(N)$, $y(N)$: Vectors with size N contening the solutions of the first order differencial equation.

**Compute y(i+1), x(i+1) for each x(i) and y(i) from $f(x,y)$ and $h$ as discusted above.**

1. Define i = 0.
2. Define h2 = 0.5 $\cdot$ h
3. Define $x(0) = x_i$ and $y(0) = y_i$ (the inicial condicions inputs)
4. **do**

      $k_1$ = h $\cdot$ f (x(i), y(i))
      $k_2$ = h $\cdot$ f (x(i) + h2, y(i) + 0.5 $\cdot k_1$)

      y(i+1) = y(i) + $k_2$
      
      x(i+1) = x(i) + h
      
      i = i + 1
      
    **while** i $\neq$ N + 1
    
 4. Return x(N) and y(N).

### The Runge-Kutta Method

As the content presented above suggests, higher expanctions will leed to minor errors, in fact a n order power series expaction leads to $\mathcal{O}(h^{n+1})$ error. There is many ways that we can evaluate $f(x,y)$, the fouth-order Runge-Kutta Method (or simply Runge-Kutta method) evaluate the derivative four times in each step: One in $y_n$, one in $y_{n+1} and twice in the mid point. The solution is computed as folows:

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

**Algorithm** *MidPoint*

**Input** $f(x,y) = \frac{dy}{dx}$: First order differencial equation to me integrated (function), $x(0)$, $y(0)$: Inicial conditions (float), $h$: Integration step (float), $N$: Number of integrations steps (positive integer).

**Output** $x(N)$, $y(N)$: Vectors with size N contening the solutions of the first order differencial equation.

**Compute y(i+1), x(i+1) for each x(i) and y(i) from $f(x,y)$ and $h$ as discusted above.**

1. Define i = 0.
2. Define h2 = 0.5 $\cdot$ h
3. Define h6 = h/6
4. Define $x(0) = x_i$ and $y(0) = y_i$ (the inicial condicions inputs)
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

There is much more to dicuss about this methods, for forther reading I sugest the chapter 16 of the book Numerical Recipier in F77 from WILLIAM H. PRESS, SAUL A. TEUKOLSKY, WILLIAM T. VETTERLING and BRIAN P. FLANNERY.




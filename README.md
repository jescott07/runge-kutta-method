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

Before w

Lets take the foloing first order differencial equation:

<p align="center">
$f(x,y) = \frac{dx}{dy}$.
</p>



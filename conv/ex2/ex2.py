from mailbox import linesep
import numpy as np
import matplotlib.pyplot as plt
from ex2 import euler
from ex2 import midpoint as mp
from ex2 import rk4
import rich

def sol(x):
    y = np.sqrt((x**3 + 7*x)/2)
    return y

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

    plt.plot(xeuler,yeuler)
    plt.plot(xeuler,yeuler_sol)
    plt.show()

    return mean_delta_euler

teste1 = comp_euler(0.01)

teste2 = comp_euler(0.1)

print(teste1,teste2)

# rich.inspect(euler)

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

    # plt.plot(xrk4,yrk4)
    # plt.plot(xrk4,yrk4_sol)
    # plt.show()    

    return mean_delta_rk4

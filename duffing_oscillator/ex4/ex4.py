from cProfile import label
from ex4 import duffing as df
from ex4 import duffing_t as dft

from matplotlib import pyplot as plt
import numpy as np

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

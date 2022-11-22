from ex3 import lorenz as lr

from matplotlib import pyplot as plt



lr.r = 10
lr.tf = 100

lr.rk4()


x1a,x2a,x3a = lr.x[:,0],lr.x[:,1],lr.x[:,2]

plt.figure()
plt.title('Lorenz System: ' r'$y_1 \times y_2$ ' 'for r = 10')
plt.xlabel(r'$y_1$')
plt.ylabel(r'$y_2$')
plt.plot(x1a,x2a)
plt.savefig('ls_1.png',dpi=300)
plt.show()
plt.close()

lr.r = 30
lr.tf = 100

lr.rk4()

x1b,x2b,x3b = lr.x[:,0],lr.x[:,1],lr.x[:,2]

plt.figure()
plt.title('Lorenz System: ' r'$y_1 \times y_2$ ' 'for r = 30')
plt.xlabel(r'$y_1$')
plt.ylabel(r'$y_2$')
plt.plot(x1b,x2b)
plt.savefig('ls_2.png',dpi=300)
plt.show()
plt.close()

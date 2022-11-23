from ex3 import lorenz as lr

from matplotlib import pyplot as plt

lr.r = 10
lr.tf = 100

lr.rk4()

x1a,x2a,x3a = lr.x[:,0],lr.x[:,1],lr.x[:,2]

fig = plt.figure()
ax = fig.gca(projection='3d')
ax.set_title('Lorenz System: for r = 30')
ax.set_xlabel('x1')
ax.set_ylabel('x2')
ax.set_zlabel('x3')
ax.plot(x1a, x2a, x3a)
plt.savefig('ls_1.png',dpi=300)
plt.close()

lr.r = 30
lr.tf = 100

lr.rk4()

x1b,x2b,x3b = lr.x[:,0],lr.x[:,1],lr.x[:,2]

fig = plt.figure()
ax = fig.gca(projection='3d')
ax.set_title('Lorenz System: for r = 30')
ax.set_xlabel('x1')
ax.set_ylabel('x2')
ax.set_zlabel('x3')
ax.plot(x1b, x2b, x3b)
plt.savefig('ls_2.png',dpi=300)
plt.close()
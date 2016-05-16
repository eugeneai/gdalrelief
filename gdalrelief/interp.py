import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate
x = np.arange(-5.01, 5.01, 0.25)
y = np.arange(-5.01, 5.01, 0.25)
xx, yy = np.meshgrid(x, y)
z = np.sin(xx**2+yy**3)
f = interpolate.interp2d(x, y, z, kind='cubic')

prec=0.03
xnew = np.arange(-2.01, 2.21, prec)
ynew = np.arange(-2.01, 2.21, prec)
znew = f(xnew, ynew)
# zz = np.sin(xnew**2+ynew**2)
#plt.plot(x, z[0, :], 'ro-', xnew, znew[0, :], 'b-')
#plt.pcolormesh(x,y,f(x,y))
#plt.plot_surface(x,y,f(x,y))
#plt.show()


from mpl_toolkits.mplot3d.axes3d import Axes3D
import matplotlib.pyplot as plt


# imports specific to the plots in this example
import numpy as np
from matplotlib import cm
from mpl_toolkits.mplot3d.axes3d import get_test_data

# Twice as wide as it is tall.
fig = plt.figure(figsize=plt.figaspect(1))

#---- First subplot
ax = fig.add_subplot(1, 1, 1, projection='3d')
# X = np.arange(-5, 5, 0.25)
# Y = np.arange(-5, 5, 0.25)
xnew, ynew = np.meshgrid(xnew, ynew)
# R = np.sqrt(X**2 + Y**2)
# Z = np.sin(R)

print (znew)

surf = ax.plot_surface(xnew, ynew, znew, rstride=1, cstride=1, cmap=cm.coolwarm,
                       linewidth=0, antialiased=True)
ax.set_zlim3d(-1.01, 1.01)

fig.colorbar(surf, shrink=0.5, aspect=10)

#---- Second subplot
#ax = fig.add_subplot(1, 2, 2, projection='3d')
#X, Y, Z = get_test_data(0.05)
#ax.plot_wireframe(X, Y, Z, rstride=10, cstride=10)

plt.show()

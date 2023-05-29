import matplotlib.pyplot as plt
from matplotlib import cm, colors
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
from scipy.special import sph_harm
from matplotlib.patches import Circle, PathPatch
import mpl_toolkits.mplot3d.art3d as art3d

def f(phi1, phi2):
    mu = 0.75
    landa = 0.4
    return -0.5*mu**2*(phi1**2+phi2**2) + 1.0/4.0*landa**2*(phi1**2+phi2**2)**2

x = np.linspace(-2.0, 2.0, 200)
y = np.linspace(-2.0, 2.0, 200)

mu = 0.75
landa = 0.4
radius = mu/landa
X, Y = np.meshgrid(x, y)
Z = f(X, Y)

p = Circle((0, radius), 2.2)


fig = plt.figure()
ax = plt.axes(projection='3d')
ax.contour3D(X, Y, Z, 200, cmap='jet')
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z');

ax.add_patch(p)
art3d.pathpatch_2d_to_3d(p, z=-0.4, zdir="z")
# ax.set_yscale('log')
plt.show()

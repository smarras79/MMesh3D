import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

plt.close('all')

print_lables=False

#
# Load coords
#
# 1) Low order:
coords_lo = np.loadtxt('COORDS_LO.dat', usecols=range(4))

x=coords_lo[:,0];
y=coords_lo[:,1];
z=coords_lo[:,2];
ip=coords_lo[:,3]+1;

fig  = plt.figure()
ax3d = fig.add_subplot(projection='3d')

ax3d.scatter(x, y, z, marker='o')

if (print_lables == True):
    for xcoords, ycoords, zcoords, label in zip(x, y, z, ip):
        ax3d.text(xcoords, ycoords, zcoords, int(label))


# 2) high order edges
del x, y, z, ip
coords_ho = np.loadtxt('COORDS_HO_edges.dat', usecols=range(4))

x=coords_ho[:,0];
y=coords_ho[:,1];
z=coords_ho[:,2];
ip=coords_ho[:,3]+1;

ax3d.scatter(x, y, z, marker='x')

if (print_lables == True):
    for xcoords, ycoords, zcoords, label in zip(x, y, z, ip):
        ax3d.text(xcoords, ycoords, zcoords, int(label))


# 3) high order faces
del coords_ho
del x, y, z, ip
coords_ho = np.loadtxt('COORDS_HO_faces.dat', usecols=range(4))

x=coords_ho[:,0];
y=coords_ho[:,1];
z=coords_ho[:,2];
ip=coords_ho[:,3]+1;

ax3d.scatter(x, y, z, marker='s')

if (print_lables == True):
    for xcoords, ycoords, zcoords, label in zip(x, y, z, ip):
        ax3d.text(xcoords, ycoords, zcoords, int(label))


# 4 high order internal
del coords_ho
del x, y, z, ip
coords_ho = np.loadtxt('COORDS_HO_vol.dat', usecols=range(4))

x=coords_ho[:,0];
y=coords_ho[:,1];
z=coords_ho[:,2];
ip=coords_ho[:,3]+1;

ax3d.scatter(x, y, z, marker='s')

if (print_lables == True):
    for xcoords, ycoords, zcoords, label in zip(x, y, z, ip):
        ax3d.text(xcoords, ycoords, zcoords, int(label))
        
        
plt.show()

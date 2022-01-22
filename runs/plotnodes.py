import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

coords_lo = np.loadtxt('COORDS_LO.dat', usecols=range(4))

x=coords_lo[:,0];
y=coords_lo[:,1];
z=coords_lo[:,2];
ip=coords_lo[:,3]+1;

fig  = plt.figure()
ax3d = fig.add_subplot(projection='3d')

ax3d.scatter(x, y, z, marker='o')

for xcoords, ycoords, zcoords, label in zip(x, y, z, ip):
    ax3d.text(xcoords, ycoords, zcoords, int(label))

plt.show()
#for i in x: 
#    scatter(x(i), y(i), z(i), 140, 'Filled'); hold on
#    print(i)
   # ipstr   = ip(i);
    
   # dx    = 0.015*max(x); dy = 0.015*max(y); 0.001*max(z); % displacement so the text
   # tx = text(x(i)+dx, y(i)+dy, z(i)+dz, int2str(ipstr), 'FontSize', 12);
#end

from sphericalRay import rayTable
from sphericalRay import geographicalRay
from sphericalRay import cartesianRay
import matplotlib.pyplot as plt
import numpy as np
"""
   read ray table file
"""
rtab = rayTable('ray-table-500-6000.hdf','ray.log')
"""
   Basic rays for various distances
"""
bray = []
nray = 10
for i in range(nray):
    delta = np.pi/5+(i-1)*np.pi/18
    bray.insert(i,rtab.getBasicRay(delta,6220000.0,6370998.0))   
"""
   Put source at lat = 63 and lon = 10, azimuth to 20 degrees
"""
lats = np.deg2rad(63.0)
lons = np.deg2rad(10.0)
xi = np.deg2rad(20.0)
"""
   Geographical rays for various distances
"""
ggray = []
for i in range(10):
    ggray.insert(i,geographicalRay(bray[i],xi,lats,lons))   
"""
   Cartesian rays for various distances
"""
ctray = []
for i in range(10):
    ctray.insert(i,cartesianRay(bray[i],xi,lats,lons))   
"""
   open figure
"""
fig = plt.figure(figsize=(15,12))
ax = fig.add_subplot(3,2,1)
ax.grid(b=True,which='major',axis='both',ls=':')
ax.set_ylim([3000000,6371000])
ax.set_ylabel("Radius [m]")
ax.set_xlabel("Epicentral distance [rad]")
for i in range(nray):
    ax.plot(bray[i].delta,bray[i].rad,'-b')

bx = fig.add_subplot(3,2,2)
bx.set_ylim([3000000,6371000])
bx.grid(b=True,which='major',axis='both',ls=':')
bx.set_ylabel("Radius [m]")
bx.set_xlabel("Latitude [rad]")
for i in range(nray):
    bx.plot(ggray[i].lat,ggray[i].rad,'-b')

cx = fig.add_subplot(3,2,3)
cx.set_ylim([3000000,6371000])
cx.grid(b=True,which='major',axis='both',ls=':')
cx.set_ylabel("Radius [m]")
cx.set_xlabel("Longitude [rad]")
for i in range(10):
    cx.plot(ggray[i].lon,ggray[i].rad,'-b')

dx = fig.add_subplot(3,2,4)
dx.grid(b=True,which='major',axis='both',ls=':')
dx.set_ylabel("Longitude [rad]")
dx.set_xlabel("Latitude [rad]")
for i in range(10):
    dx.plot(ggray[i].lat,ggray[i].lon,'-b')

ex = fig.add_subplot(3,2,5)
ex.grid(b=True,which='major',axis='both',ls=':')
ex.set_ylabel("Y [m]")
ex.set_xlabel("X [m]")
for i in range(10):
    ex.plot(ctray[i].x,ctray[i].y,'-b')

fx = fig.add_subplot(3,2,6)
fx.grid(b=True,which='major',axis='both',ls=':')
fx.set_ylabel("Z [m]")
fx.set_xlabel("X [m]")
for i in range(10):
    fx.plot(ctray[i].x,ctray[i].z,'-b')

plt.savefig("basic-geographical-cartesian-ray.png", dpi=100,format='png')
plt.show()


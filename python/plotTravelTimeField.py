#!/usr/bin/env python3
#
#  Plot a 2D travel time field defined on a regular delta-radius grid
# -------------------------------------------------------------------
import h5py
import numpy as np
import matplotlib.pyplot as plt
from argparse import ArgumentParser
from argparse import ArgumentDefaultsHelpFormatter


############################################################
# MAIN
############################################################
ap = ArgumentParser(description="Plot 2D travel time field or time difference", formatter_class=ArgumentDefaultsHelpFormatter)
ap.add_argument("ttfile", help="HDF file with travel time grid")
ap.add_argument("--niso", help="Number of isolines", default='50')
ap.add_argument("--second", help="Second compatible file with travel times of other phase", default='None')
ap.add_argument("--tmax", help="Plot travel time grid up to tmax", default='1.e6')
args = ap.parse_args()
ttfile = args.ttfile
niso = int(args.niso)
ttfile2 = args.second
tmax = float(args.tmax)

#
#  read travel time and slowness on grid
#  ttime[jr,jd] = travel time at grid point rmin+(jr-1)*dr and delmin+(jd-1)*dstep
#
fid = h5py.File(ttfile,'r')
ndel, nr = fid.attrs['dimensions']
dstep, dr = fid.attrs['grid_spacings_rad_m']
delmin, rmin = fid.attrs['minimum_values_rad_m']
ttime1 = np.copy(fid['ttime'])
slowness1 = np.copy(fid['slowness'])
delmax = delmin+(ndel-1)*dstep

#
#  read second file if present
#
if ttfile2 != 'None':
    fid2 = h5py.File(ttfile2,'r')
    ttime2 = np.copy(fid2['ttime'])
    slowness2 = np.copy(fid2['slowness'])
    ttime = np.where(ttime2-ttime1 < tmax,ttime2-ttime1,0.0)
    slowness = slowness2-slowness1
else:
    ttime = np.where(ttime1 < tmax,ttime1,0.0)
    slowness = slowness1

#
#  set up Quadmesh defining x and y coordinates of grid points
#  see matplotlib docu on Quadmesh for a description
#
delta = np.linspace(delmin, delmax, num=ndel)
xg = np.zeros((nr,ndel))
yg = np.zeros((nr,ndel))
#
#  grid points as cell centers
#
for k in range(0,nr):
    r = rmin+k*dr
    xg[k,:] = r*np.sin(delta-0.5*(delmin+delmax))
    yg[k,:] = r*np.cos(delta-0.5*(delmin+delmax))
#
#  cell edges
#
xedge = np.zeros((nr+1,ndel+1))
yedge = np.zeros((nr+1,ndel+1))
for k in range(0,nr+1):
    r = rmin-0.5*dr+k*dr
    xedge[k,0:ndel] = r*np.sin(delta-0.5*dstep-0.5*(delmin+delmax))
    yedge[k,0:ndel] = r*np.cos(delta-0.5*dstep-0.5*(delmin+delmax))
    xedge[k,ndel] = r*np.sin(delta[ndel-1]+0.5*dstep-0.5*(delmin+delmax))
    yedge[k,ndel] = r*np.cos(delta[ndel-1]+0.5*dstep-0.5*(delmin+delmax))
#
#  plotting
#
fig = plt.figure(figsize=(20,10))
ax1 = fig.add_subplot(2, 1, 1)
ax1.set_xlabel('x-coordinate')
ax1.set_ylabel('y-coordinate')
ax1.axis('equal')
ax1.set_title('Slowness field')
vmin = np.min(slowness)
vmax = np.max(slowness)
cf = ax1.pcolormesh(xedge, yedge, slowness, cmap='coolwarm_r', vmin=vmin, vmax=vmax, shading='flat')
ax1.contour(xg,yg,slowness,niso)
fig.colorbar(cf, shrink=1.0, pad=0.03)
print('Extrema of slowness: ',vmin,vmax)
#
#
ax2 = fig.add_subplot(2, 1, 2)
ax2.set_xlabel('x-coordinate')
ax2.set_ylabel('y-coordinate')
ax2.axis('equal')
ax2.set_title('Travel time field')
vmin = np.min(ttime)
vmax = np.max(ttime)
cf = ax2.pcolormesh(xedge, yedge, ttime, cmap='coolwarm_r', vmin=vmin, vmax=vmax, shading='flat')
ax2.contour(xg,yg,ttime,niso)
fig.colorbar(cf, shrink=1.0, pad=0.03)
print('Extrema of travel time: ',vmin,vmax)

fig.savefig(ttfile.split('.')[0]+'.png', dpi=300, bbox_inches='tight')

plt.show()




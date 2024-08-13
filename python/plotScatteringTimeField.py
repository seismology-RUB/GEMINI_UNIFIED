#!/usr/bin/env python3
#
#  Plot the difference between the arrival time of a scattered wave
#  and the arrival time of the teleseismic wave at a given reciever
# -------------------------------------------------------------------
import h5py
import numpy as np
import matplotlib.pyplot as plt
from argparse import ArgumentParser
from argparse import ArgumentDefaultsHelpFormatter


############################################################
# MAIN
############################################################
ap = ArgumentParser(description="Plot time difference between a scattered wave and the direct wave at a receiver",
                    formatter_class=ArgumentDefaultsHelpFormatter)
ap.add_argument("tele", help="HDF file with grid of teleseismic travel times")
ap.add_argument("localp", help="HDF file with grid of local P travel times")
ap.add_argument("locals", help="HDF file with grid of local S travel times")
ap.add_argument("--niso", help="Number of isolines", default='50')
ap.add_argument("--delta", help="Epicentral distance of receiver in degrees", default='90')
ap.add_argument("--tsmax", help="Max effective scattering time plotted", default='1.e6')
args = ap.parse_args()
telefile = args.tele
localfile_p = args.localp
localfile_s = args.locals
niso = int(args.niso)
tsmax = float(args.tsmax)
deltarec = np.deg2rad(float(args.delta))

#
#  read teleseismic travel time grid
#  ttime[jr,jd] = travel time at grid point rmin+(jr-1)*dr and delmin+(jd-1)*dstep
#
fid = h5py.File(telefile,'r')
ndel, nr = fid.attrs['dimensions']
dstep, dr = fid.attrs['grid_spacings_rad_m']
delmin, rmin = fid.attrs['minimum_values_rad_m']
ttime_tele = np.copy(fid['ttime'])
delmax = delmin+(ndel-1)*dstep

#
#  read local travel time grid for P
#  ttime[jr,jd] = travel time at grid point rmin+(jr-1)*dr and delmin+(jd-1)*dstep
#
fid2 = h5py.File(localfile_p,'r')
ttime_local_p = np.copy(fid2['ttime'])
ndel2, nr2 = fid2.attrs['dimensions']
dstep2, dr2 = fid2.attrs['grid_spacings_rad_m']
delmin2, rmin2 = fid2.attrs['minimum_values_rad_m']
if ndel != ndel2 or nr != nr2:
    print(ndel,ndel2,nr,nr2)
    raise Exception('Incompatible grid dimensions)')
if abs(dstep-dstep2) > 1.e-5 or abs(dr-dr2) > 1.e-5:
    print(dstep,dstep2,dr,dr2,rmin,rmin2)
    raise Exception('Incompatible grid spacings')
if abs(rmin-rmin2) > 1.e-5:
    print(rmin,rmin2)
    raise Exception('Incompatible grid minimum radii')


#
#  read local travel time grid for S
#  ttime[jr,jd] = travel time at grid point rmin+(jr-1)*dr and delmin+(jd-1)*dstep
#
fid3 = h5py.File(localfile_s,'r')
ttime_local_s = np.copy(fid3['ttime'])
ndel3, nr3 = fid2.attrs['dimensions']
dstep3, dr3 = fid2.attrs['grid_spacings_rad_m']
delmin3, rmin3 = fid2.attrs['minimum_values_rad_m']
if ndel != ndel3 or nr != nr3:
    print(ndel,ndel3,nr,nr3)
    raise Exception('Incompatible grid dimensions)')
if abs(dstep-dstep3) > 1.e-5 or abs(dr-dr3) > 1.e-5:
    print(dstep,dstep3,dr,dr3,rmin,rmin3)
    raise Exception('Incompatible grid spacings')
if abs(rmin-rmin2) > 1.e-5:
    print(rmin,rmin3)
    raise Exception('Incompatible grid minimum radii')

#
#  find travel time of teleseismic wave at receiver (r=rmax and delta = deltarec)
#  delmin+jd*dstep <= deltarec  --> jd <= (deltarec-delmin)/dstep = floor((deltarec-delmin)/dstep)
#
jrec_tele = int(np.floor((deltarec-delmin)/dstep))
ttrec = ttime_tele[nr-1,jrec_tele]
print('Travel time of direct teleseismic wave = ',ttrec)

#
#  travel time of a wave scattered from a grid point =
#  teleseismic travel time to scatterer plus local travel time from scatterer to station 
#  Walk through scattering points:
#
ttstreu_p = np.zeros((nr,ndel))
ttstreu_s = np.zeros((nr,ndel))
for jd in range(ndel):
    dels_tele = delmin+jd*dstep
    dels_local = abs(deltarec-dels_tele)
    js_local = int(np.floor((dels_local-delmin2)/dstep2))
    for k in range(nr):
        tspp = ttime_tele[k,jd]+ttime_local_p[k,js_local]-ttrec
        tsps = ttime_tele[k,jd]+ttime_local_s[k,js_local]-ttrec
        if tsps <= tsmax:
            ttstreu_p[k,jd] = tspp
            ttstreu_s[k,jd] = tsps
        elif tsps > tsmax and tspp <= tsmax:
            ttstreu_p[k,jd] = tspp
            ttstreu_s[k,jd] = tspp
        else:
            ttstreu_p[k,jd] = 0.0
            ttstreu_s[k,jd] = 0.0
#
#-----------------   PLOTTING  ---------------------------------
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
ax1.set_title('Arrival time difference between scattered P and direct P wave')
vmin = np.min(ttstreu_p)
vmax = np.max(ttstreu_p)
cf = ax1.pcolormesh(xedge, yedge, ttstreu_p, cmap='coolwarm_r', vmin=vmin, vmax=vmax, shading='flat')
ax1.contour(xg,yg,ttstreu_p,niso)
fig.colorbar(cf, shrink=1.0, pad=0.03)
print('Extrema of PP-scattering time: ',vmin,vmax)

ax2 = fig.add_subplot(2, 1, 2)
ax2.set_xlabel('x-coordinate')
ax2.set_ylabel('y-coordinate')
ax2.axis('equal')
ax2.set_title('Arrival time difference between scattered S and direct P wave')
vmin = np.min(ttstreu_s)
vmax = np.max(ttstreu_s)
cf = ax2.pcolormesh(xedge, yedge, ttstreu_s, cmap='coolwarm_r', vmin=vmin, vmax=vmax, shading='flat')
ax2.contour(xg,yg,ttstreu_s,niso)
fig.colorbar(cf, shrink=1.0, pad=0.03)
print('Extrema of PS-scattering time: ',vmin,vmax)

fig.savefig('ttstreu.png', dpi=300, bbox_inches='tight')

plt.show()


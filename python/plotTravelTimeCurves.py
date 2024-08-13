#!/usr/bin/env python3
#
#  Plot delta and ttime and tau versus slowness or turning point radius
# -------------------------------------------------------------------
import h5py
import numpy as np
import matplotlib.pyplot as plt
from argparse import ArgumentParser
from argparse import ArgumentDefaultsHelpFormatter


############################################################
# MAIN
############################################################
ap = ArgumentParser(description="Plot travel time curves",
                    formatter_class=ArgumentDefaultsHelpFormatter)
ap.add_argument("hdffile", help="HDF file with ray table")
ap.add_argument("--re", help="receiver radius in km", default='6345')
ap.add_argument("--slow", help="plot delta and tt versus slowness", action='store_true')
args = ap.parse_args()
hdffile = args.hdffile
plotslow = args.slow
re = float(args.re)*1.e3

#
#  read out ray table file
#  raytable[jtp,jnr,0] = delta(jtp,jr)
#  raytable[jtp,jnr,1] = ttime(jtp,jr)
#
fid = h5py.File(hdffile,'r')
slowness = np.copy(fid['rayParameters'])
rtp = np.copy(fid['turningPointRadii'])
rnod = np.copy(fid['receiverRadii'])
raytable = np.copy(fid['rayTable'])
nnod = np.size(rnod)

je = np.nonzero(rnod > re)[0][0]-1               # index of receiver in rnod array just below re
print('Closest receiver below re: ',rnod[je])
re = rnod[je]

jtpe = np.nonzero(rtp < re)[0]                   # array of indices of turning points below receiver

#
#  delta, ttime and tau from surf to tp to receiver
#
delta_surf_tp_rec = np.rad2deg(raytable[jtpe,je,0]+raytable[jtpe,nnod-1,0])
tt_surf_tp_rec = raytable[jtpe,je,1]+raytable[jtpe,nnod-1,1]
tau_surf_tp_rec = tt_surf_tp_rec-slowness[jtpe]*np.deg2rad(delta_surf_tp_rec)

#
#  delta, ttime and tau from surf to receiver
#
delta_surf_rec = np.rad2deg(raytable[jtpe,nnod-1,0]-raytable[jtpe,je,0])
tt_surf_rec = raytable[jtpe,nnod-1,1]-raytable[jtpe,je,1]
tau_surf_rec = tt_surf_rec-slowness[jtpe]*np.deg2rad(delta_surf_rec)

fig = plt.figure(figsize=(20,10))
ax1 = fig.add_subplot(2, 3, 1)
ax1.set_ylabel('Epicentral distance from surface via turning point (deg)')
if plotslow:
    ax1.set_xlabel('Slowness (s/rad)')
    ax1.plot(slowness[jtpe], delta_surf_tp_rec, color='b')
else:
    ax1.set_xlabel('Turning point radius (m)')
    ax1.plot(rtp[jtpe], delta_surf_tp_rec, color='b')

ax2 = fig.add_subplot(2, 3, 2)
ax2.set_ylabel('Travel time from surface via turning point (s)')
if plotslow:
    ax2.set_xlabel('Slowness (s/rad)')
    ax2.plot(slowness[jtpe], tt_surf_tp_rec, color='b')
else:
    ax2.set_xlabel('Turning point radius (m)')
    ax2.plot(rtp[jtpe], tt_surf_tp_rec, color='b')

ax3 = fig.add_subplot(2, 3, 4)
ax3.set_ylabel('Epicentral distance from surface (s)')
if plotslow:
    ax3.set_xlabel('Slowness (s/rad)')
    ax3.plot(slowness[jtpe], delta_surf_rec, color='b')
else:
    ax3.set_xlabel('Turning point radius (m)')
    ax3.plot(rtp[jtpe], delta_surf_rec, color='b')

ax4 = fig.add_subplot(2, 3, 5)
ax4.set_ylabel('Travel time from surface (s)')
if plotslow:
    ax4.set_xlabel('Slowness (s/rad)')
    ax4.plot(slowness[jtpe], tt_surf_rec, color='b')
else:
    ax4.set_xlabel('Turning point radius (m)')
    ax4.plot(rtp[jtpe], tt_surf_rec, color='b')

ax5 = fig.add_subplot(2, 3, 3)
ax5.set_ylabel('Tau from surface via turning point (s)')
if plotslow:
    ax5.set_xlabel('Slowness (s/rad)')
    ax5.plot(slowness[jtpe], tau_surf_tp_rec, color='b')
else:
    ax5.set_xlabel('Turning point radius (m)')
    ax5.plot(rtp[jtpe], tau_surf_tp_rec, color='b')

ax6 = fig.add_subplot(2, 3, 6)
ax6.set_ylabel('Tau from surface (s)')
if plotslow:
    ax6.set_xlabel('Slowness (s/rad)')
    ax6.plot(slowness[jtpe], tau_surf_rec, color='b')
else:
    ax6.set_xlabel('Turning point radius (m)')
    ax6.plot(rtp[jtpe], tau_surf_rec, color='b')

plt.show()
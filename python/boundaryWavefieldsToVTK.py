#!/usr/bin/env python3
#
#  Write snapshots of boundary wavefields to VTK files
#------------------------------------------------------
#For Anaconda users import the following command:
from pyevtk.hl import pointsToVTK
#users who have installed evtk in a different way may use
#from evtk.hl import pointsToVTK
from argparse import ArgumentParser
from argparse import ArgumentDefaultsHelpFormatter
from boundaryWavefield import BoundaryWavefield
#----------------------------------------------------------
#  read command line arguments
#
ap = ArgumentParser(description = "Write Gemini boundary wavefields to VTK",
                    formatter_class=ArgumentDefaultsHelpFormatter)
ap.add_argument("hdffile", help = "Name of HDF file with boundary data (synthetic seismograms)")
ap.add_argument("eventid", help = "Event ID")
ap.add_argument("vtkbasename", help = "VTK file basename")
ap.add_argument("--timestep", help = "start,stop,step for timestepping", default = "0,20,1")
args = ap.parse_args()
#
# variable arguments
#
hdffile = args.hdffile
eventid = args.eventid
vtkbasename = args.vtkbasename
timestep = args.timestep.split(",")
ita = int(timestep[0])
itb = int(timestep[1])
its = int(timestep[2])

bw = BoundaryWavefield(hdffile,eventid)
for it in range(ita,itb,its):
    print("Time step:",it)
    (coords,veltrac) = bw.readCoordsVelocityTraction(eventid,it,vtkbasename)
    vtkname = vtkbasename+"_{:08d}".format(it)
    pointsToVTK(vtkname,coords[0,:],coords[1,:],coords[2,:],
                data = {"1_vx":veltrac[0,:],"2_vy":veltrac[1,:],"3_vz":veltrac[2,:],
                        "4_tx":veltrac[3,:],"5_ty":veltrac[4,:],"6_tz":veltrac[5,:]})



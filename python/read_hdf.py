#!/usr/bin/env python3
#
#  Read HDF using python
#-------------------------------------------------------------------
import h5py
from evtk.hl import pointsToVTK
import numpy as np
    
fid = h5py.File('prep_specfem_coupling.hdf','r')
[dimx,dimy,dimz] = fid.attrs["boxElementDimensions"]
print(dimx,dimy,dimz)
[clat,clon] = fid.attrs["centerLatLon"]
print(clat,clon)
print("Attributes:")
for att in fid.attrs:
    print(att)
print("Attribute items:")
for att in fid.attrs.items():
    print(att)
print("Attribute values:")
for att in fid.attrs.values():
    print(att)
print("Attribute keys:")
for att in fid.attrs.keys():
    print(att)
print("File group:")
for obj in fid.values():
    print("Object: ",obj)
    print("Name: ",obj.name)
    print("Shape: ",obj.shape)
    print("Size: ",obj.size)
    print("Dtype: ",obj.dtype)
print("Number of group members: ",len(fid.keys()))
print("Some values of first dataset: ")
print(fid["gllLatLonBottom"][0:20,1])
print("gllRadii: ")
print(fid["gllRadii"][:])

x = fid["gllLatLonBottom"][:,0]
y = fid["gllLatLonBottom"][:,1]
z = np.zeros((np.size(x)),dtype = "float32")
temp = np.sqrt(x*x+y*y)
basename = "specfem"
vtkfile = basename+"_{:08d}".format(25)
print(vtkfile)
pointsToVTK(vtkfile, x,y,z,data = {"1_temp" : temp})

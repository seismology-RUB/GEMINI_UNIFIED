import h5py
import numpy as np

class SpecfemBox:
    """
    Define an object to store the data relevant to the specfem box
    """
    def __init__(self, file = "specfem_boundary_data.hdf", readLatLonOnly=False):
        self.boundaryPoints = []
        box = h5py.File(file, "r")
        if readLatLonOnly:
            self.centerLat = box.attrs["centerLatLon"][0]
            self.centerLon = box.attrs["centerLatLon"][1]
            self.dimension = box.attrs["boxElementDimensions"]
        else:
            self.centerLat = box.attrs["centerLatLon"][0]
            self.centerLon = box.attrs["centerLatLon"][1]
            self.dimension = box.attrs["boxElementDimensions"]
            self.boundaryPoints.extend(box["gllRadLatLonNormalsBottom"])
            self.boundaryPoints.extend(box["gllRadLatLonNormalsSides"])
            self.SortingIndices = box["gllSortingIndicesSides"]
            # Save only the  radius, lattitude and longitude of the boundary points
            self.boundaryPoints = np.asarray(self.boundaryPoints).transpose()[:3]
            self.uniqueRadii = box["gllUniqueRadii"]


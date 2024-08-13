from seismicArray import Profile
from specfemBox import SpecfemBox
from boundaryWavefield import BoundaryWavefield

"""
Extract seismogram data from the hdf file created by GEMINI.
"""

folder = "/Users/tommy/git/Gemini/shyam_test/comparisonGeminiHybrid/"
geminiFolder = "gemini/"  # Der muss noch dynamisch gemacht werden. Die Ordner "2019.12" und "e0001.001.12" können aus der EventID erzeugt werden
# For testing purposes
hdffile = folder + geminiFolder + "specfem_coupling_seismograms.hdf"
eventid = "171117_223425"
side = "right"                    # left, right, front, back or bottom
direction = "LtoR"                  # left to right (LtoR), front to back (FtoB)
profileType = "horizontal"         # horizontal or vertical
target = 0.                    # x (front or back) or y (left or right) coordinate for depth profiles
buffer = 20.                     # Buffer for the start-time of the seismogram in seconds (will be an input later on)
rotate = True


box = SpecfemBox(file = folder + "specfem_boundary_data.hdf", readLatLonOnly=False)
bw = BoundaryWavefield(hdffile, eventid, si=True)

array = Profile(bw, profileType, side, direction=direction)
if profileType == "vertical":
    array.locateTarget(target)
    array.extractProfile(box.centerLat, box.centerLon, rotate = rotate)
else:
    if side == "bottom":
        rad = box.uniqueRadii[0]
        array.locateTarget(target)
    else:
        rad = box.uniqueRadii[-1]
    array.extractProfile(box.centerLat, box.centerLon, radius=rad, rotate = rotate)

array.writeStreamData(folder, geminiFolder)
array.writeArrayToVTK(folder)
array.writeSpecfemStationFile(folder)
array.writeGeminiStationFile(folder)


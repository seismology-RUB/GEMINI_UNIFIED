"""
Module to create station files for an array of seismic Stations for Gemini and Specfem

"""
import utils
from obspy import read, Trace
from seismicStation import SeismicStation
from pyevtk.hl import pointsToVTK
import numpy as np

class SeismicArray:
    """
    Class to handle operations for an array of seismic stations
    """

    def __init__(self, name = "", si = True, eventid = ""):
        self.name = name
        self.stations = []
        self.test = {}

        if si:
            self.multi = 1000.
        else:
            self.multi = 1.
        self.rEarth = 6371.*self.multi

        self.outfile = self.name

        if eventid:
            self.outfile += "_" + eventid

    def writeArrayToVTK(self, folder=""):
        """

        :param folder:
        :return:
        """
        filename = folder + self.name
        print("writing file: " + filename + ".vtu")
        xs = np.asarray([station.x for station in self.stations])
        ys = np.asarray([station.y for station in self.stations])
        zs = np.asarray([station.z for station in self.stations])
        print(len(xs), len(ys), len(zs))
        print(xs.dtype.itemsize , ys.dtype.itemsize , zs.dtype.itemsize)
        pointsToVTK(filename, xs, ys, zs)

    def writeSpecfemStationFile(self, folder=""):
        """

        :param folder:
        :return:
        """
        print("Writing STATIONS file for specfem use")
        with open(folder + "STATIONS_" + self.name, "w") as specfemStations:
            for station in self.stations:
                values = [station.networkCode, station.name, station.y, station.x, station.z, self.rEarth + station.z]
                specfemStations.writelines(' '.join(str(j) for j in values) + "\n")


    def writeGeminiStationFile(self, folder=""):
        """

        :param folder:
        :return:
        """

        print("Writing GEMINI stations file")
        with open(folder + self.name + ".sta", "w") as geminiStations:
            geminiStations.writelines("S \n")
            for station in self.stations:
                values = [station.name, station.networkCode, station.lat, station.lon, station.elevation]
                geminiStations.writelines(' '.join(str(j) for j in values) + "\n")

    def writeStreamData(self, folder, pylotfolder):
        """
        Write the obspy stream objects into binary files readable by PyLot.
        Write point data as vtk to view with Paraview.
        :param stream_access: data dictionary
        :param side: which side of the box
        :param profileType: vertical or horzontal profile
        :return:
        """

        self.outfile += ".dat"
        print("writing file: " + self.outfile)
        self.stream.write(folder + pylotfolder + self.outfile, format='PICKLE')

    def addStation(self, station):
        """
        Add a station to the list of stations in the array.

        :param station: object of type SeismicStation
        :return:
        """

        self.stations.append(station)



class Profile(SeismicArray):
    def __init__(self, bw, profile = "", side = "", direction=None, eventid="" ):
        """
        :param profile: Type of the desired profile. Options are: horizontal and vertical
        :param side: Specify the side of the box on which the profile ist created.
                    Options are: left, right, front, back.

        :param group: A paramter that should be inherited from the boundary wavefield class
        :param dt: A paramter that should be inherited from the boundary wavefield class
        """

        if profile == "horizontal" or profile == "vertical":
            self.profile = profile
        else:
            raise ValueError("Input 'horizontal' or 'vertical' as type of the profile")

        if side is None:
            raise ValueError(
                "A side of the box must be specified for the profile. Options are 'left', 'right', 'front', 'back' or 'bottom'.")
        else:
            self.side = side

        if self.side == "bottom":
            if direction is None:
                raise ValueError("Set a direction of the profile: Left to right (LtoR) or front to back (FtoB).")
            else:
                self.direction = direction

        name = self.profile + "-profile_side_" + self.side + "_"
        if self.side == "bottom":
            name += self.direction + "_"
        super().__init__(name=name, eventid=eventid)
        self.bw = bw
        self.target = 0.
        self.stream = read()
        self.stream.clear()
        #todo radius need to be passed to here from the calling routine
        self.radius = 0.

        if self.profile == "vertical":
            self.outfile += "target_" + str(int(self.target)) + "_velocity"
        elif self.profile == "horizontal":
            self.outfile +=  "radius_" + str(int(self.radius)) + "_velocity"
            if self.side == "bottom":
                self.outfile += "_direction_" + self.direction + "_target_" + str(int(self.target)) + "_velocity"
        else:
            raise ValueError("The profile is wronly specified: ", self.profile)
        if self.bw.si_b:
            self.outfile += "_si"

    def locateTarget(self, target):
        """
        Based on the user input: Find the Grid point closest to the location specified as input and add the new
            target location to the objet.

        :param self: BoundaryWaveField object.
        :param target: Targets other coordinate (if front/back that means x, if left/right that means y)
        """
        distanceToTarget = 2500. * self.multi  # If the point is in the box the coordinate can be at maximum 2400km apart from the target
        newTarget = 0.
        for proc in self.bw.group.values():
            coordinates = proc["coordinates"]
            for coords in coordinates:
                xs, ys, zs = coords[0], coords[1], coords[2]
                xc, yc, zc = utils.mapSphericalChunkToCartesianBox(xs, ys, zs, si=self.bw.si_b)
                if self.side == "left" or self.side == "right" or (self.side == "bottom" and self.direction == "LtoR"):
                    if target < 900 * self.multi and target > -900 * self.multi:
                        delta_y = abs(target - yc)
                        if delta_y < distanceToTarget:
                            distanceToTarget = delta_y
                            newTarget = yc
                    else:
                        raise ValueError("Target coordinate not in Box!")
                elif self.side == "front" or self.side == "back" or (
                        self.side == "bottom" and self.direction == "FtoB"):
                    if target < 1200 * self.multi and target > -1200 * self.multi:
                        delta_x = abs(target - xc)
                        if delta_x < distanceToTarget:
                            distanceToTarget = delta_x
                            newTarget = xc
                    else:
                        raise ValueError("Target coordinate not in Box!")
                else:
                    raise ValueError("Wrong side: ", self.side)
        self.target = newTarget

    def extractProfile(self, boxCenterLat, boxCenterLon, radius=None, rotate=True):
        """
        Extract seismograms for a specific profile from the boundary wavefield data

        :param self: boundaryWavefield object
        :param start: Tuple containing the start Time for the seismogram and the corresponding index for the
                      data array.
        :param radius: A Radius must be specified for a horizontal profile (optional)
        :return: unsorted data along the desired profile
        """
        tmp = []
        eps = 0.2 * self.multi
        switch = {"front": -675 * self.multi, "back": 675 * self.multi, "left": -900 * self.multi,
                  "right": 900 * self.multi}

        if self.profile == "horizontal" and radius == None:
            raise ValueError("A radius must be specified for a horizontal profile")
        else:
            self.radius = radius
        cnt = 0
        for proc in self.bw.group.values():
            coordinates = proc["coordinates"]
            data = proc["velocityTraction"]
            for coords, index in zip(coordinates, range(len(coordinates))):
                if self.profile == "horizontal":
                    r = coords[3]
                    if abs(r - self.radius) < eps:
                        xs, ys, zs = coords[0], coords[1], coords[2]
                        xc, yc, zc = utils.mapSphericalChunkToCartesianBox(xs, ys, zs, si=self.bw.si_b)
                        if self.side == "left" or self.side == "right":
                            if abs(xc - switch[self.side]) < eps:
                                if ys not in tmp:
                                    cnt += 1
                                    self.addStation(data, index, coords, rotate, boxCenterLat, boxCenterLon, cnt)
                                    tmp.append(ys)
                        elif self.side == "front" or self.side == "back":
                            if abs(yc - switch[self.side]) < eps:
                                if xs not in tmp:
                                    cnt += 1
                                    self.addStation(data, index, coords, rotate, boxCenterLat, boxCenterLon, cnt)
                                    tmp.append(xs)
                        elif self.side == "bottom":
                            if self.direction == "LtoR":
                                if abs(yc - self.target) < eps:
                                    if xs not in tmp:
                                        cnt += 1
                                        self.addStation(data, index, coords, rotate, boxCenterLat, boxCenterLon, cnt)
                                        tmp.append(xs)
                            elif self.direction == "FtoB":
                                if abs(xc - self.target) < eps:
                                    if ys not in tmp:
                                        cnt += 1
                                        self.addStation(data, index, coords, rotate, boxCenterLat, boxCenterLon, cnt)
                                        tmp.append(ys)
                        else:
                            raise ValueError(
                                "The side of the box must be correctly specified: " + self.side + " Options are 'left', 'right', 'front', 'back' or 'bottom'.")
                elif self.profile == "vertical":
                    xs, ys, zs = coords[0], coords[1], coords[2]
                    xc, yc, zc = utils.mapSphericalChunkToCartesianBox(xs, ys, zs, si=self.bw.si_b)
                    if self.side == "left" or self.side == "right":
                        if abs(xc - switch[self.side]) < eps and abs(yc - self.target) < eps:
                            if zs not in tmp:
                                self.addStation(data, index, coords, rotate, boxCenterLat, boxCenterLon, cnt)
                                tmp.append(zs)
                    elif self.side == "front" or self.side == "back":
                        if abs(yc - switch[self.side]) < eps and abs(xc - self.target) < eps:
                            if zs not in tmp:
                                self.addStation(data, index, coords, rotate, boxCenterLat, boxCenterLon, cnt)
                                tmp.append(zs)
                    else:
                        raise ValueError(
                            "The side of the box must be correctly specified: " + self.side + " Options are 'left', 'right', 'front', 'back' or 'bottom'.")
                else:
                    raise ValueError("The profile is wronly specified: ", self.profile)

    def addStation(self, data, index, coords, rotate, boxCenterLat, boxCenterLon, cnt):
        """
        This function adds a station and its corresponding data to a given seismic array. It assumes that the station is
        created during the creation of the array and no data on it exists beforehand. In addition it is assumed the the
        input data is provided by SPECFEM

        :param data: list - Seismogram data to add to the station
        :param index: int -
        :param coords: list - Cartesian coordinates of the station inside the SPECFEM box
        :param rotate: bool - Defines if the data is rotated to the ZNE coordinate system
        :param boxCenterLon, boxCenterLat: Center coordinates of the SPECFEM box.
        :param cnt: int - Station counter
        :return:
        """

        station = SeismicStation(name = "ST"+ str(cnt).zfill(3) , lat = coords[4], lon=coords[5], networkCode = "PR")
        station.r = coords[3]
        station.specfem_data = dict(zip(station.components, data[self.bw.start[0]:, index, :3].transpose()))
        station.x, station.y, station.z = coords[0], coords[1], coords[2]
        station.dt = self.bw.dt
        if rotate:
            station.calculatePropagationDirection(self.bw.event)
            station.calculateEpicentralCoordinates(self.bw.event)
            station.rotateToZNE(self.bw.event, boxCenterLat, boxCenterLon)

        station.addTraces(self.bw.event, self.stream, type="SPECFEM")
        self.stations.append(station)

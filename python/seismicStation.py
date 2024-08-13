"""
Class to handle the information on seismic stations

"""
import numpy as np
from obspy import read, Trace
from obspy.clients.iris import Client
import utils
import axesRotation as ar
import h5py
from matplotlib.lines import Line2D
import matplotlib.pyplot as plt

class SeismicStation:
    """
    Class to handle the seismogram files created by specfem
    """

    rEarth = 6371000.

    def __init__(self, name="", lat=0., lon=0., elevation=0., institute_name="undef", networkCode="XX"):

        self.name = name
        self.lat = lat
        self.lon = lon
        self.elevation = elevation
        self.institute_name = institute_name
        self.networkCode = networkCode

        self.specfem_data = {}
        self.gemini_data = {}
        self.axisem_data = {}
        self.r = SeismicStation.rEarth + self.elevation
        self.delta = 0.
        self.xi = 0.
        self.coordinateSystem = ""
        self.x = 0
        self.y = 0
        self.z = 0
        self.startTime = 0.
        self.dt = 0.
        self.components = ["Z", "N", "E"]#["BXX", "BXY", "BXZ"]

    def readStationData(self, stationData):
        """
        Read station parameter from input file

        :param stationData: Line from the input file containing the event data
        """
        # This routine currently does not support the added buffer time. this has to be changed!
        self.name = stationData[0]
        self.lat = float(stationData[2])
        self.lon = float(stationData[3])
        self.elevation = float(stationData[4])
        self.networkCode = stationData[1]

    def readCoordSystem(self, stationData):
        """
        :param stationData: This is the first line in the event List
        :return:
        """
        self.coordinateSystem = stationData[0]

    def readSpecfemData(self, event, folder, type="velocity", stdBasename=""):
        """

        :param folder: Folder where the data is placed
        :param type: Type of seismogram: velocity, displacement, acceleration
        :param standartBasename: standart convention for filenaming
        :return:
        """
        component = ["BXX", "BXY", "BXZ"]
        seismoType = {"velocity":"v", "acceleration":"a", "displacement":"d"}
        for comp in component:
            if stdBasename == "":
                self.basename = event.eventId + "." + self.networkCode + "." + self.name
            else:
                self.basename = stdBasename
            filename = self.basename + "." + comp + ".sem" + seismoType[type]
            tmp = np.loadtxt(folder + filename).transpose()
            self.dt = (tmp[0][-1] - tmp[0][0])/(len(tmp[0])-1)
            self.specfem_data[comp] = tmp[1]

    def readAxisemData(self, folder, type):
     """

     :param folder: Folder where the data is placed
     :param type: Type of seismogram: velocity, displacement, acceleration
     :return:
     """

     standardName = "_post_mij_conv0000_"
     component = ["Z","N","E"]
     seismoType = {"velocity":"velo", "acceleration":"accel", "displacement":"disp"}

     for comp in component:
         self.basename = self.name + "_" + self.networkCode
         filename = self.basename + seismoType[type] + standardName + comp + ".dat"
         tmp = np.loadtxt(folder + filename).transpose()
         self.dt = (tmp[0][-1] - tmp[0][0])/(len(tmp[0])-1)
         self.axisem_data[comp] = tmp[1]


    def readGeminiData(self, event, file):
        """
        Read the data from the Gemini simulation at the station

        :param event: object of type seismicEvent
        :param file: full path to the file containing the synthetic seismograms
        :return:
        """
        gemini_syn_file = h5py.File(file, "r")
        group = gemini_syn_file[event.eventId]
        self.dataset = group[self.networkCode+ "."+self.name]
        self.dt = self.dataset.attrs["Time"][1]
        ind = int(event.cutoffTime / self.dt)
        self.gemini_data = dict(zip(self.components, self.dataset[:, ind:]))
        self.recordLength = len(self.dataset[0, :]) * self.dt

    def addTraces(self, event, stream, type = "SPECFEM", singleStation = False):
        """
        Methode to store the profile data as obspy stream object and prepare the coordinates of the
        points in the profile such that they can be written to a vtk file for plotting

        :param self: BoundaryWaveField object
        """

        #Stream objects to be used for a single station
        if singleStation:
            self.stream = read()
            self.stream.clear()

        # Add data to stream object
        for comp in self.components:
            trace = Trace()
            if type == "SPECFEM":
                trace.data = self.specfem_data[comp]
            elif type == "GEMINI":
                trace.data = self.gemini_data[comp]
            else:
                raise ValueError("Wrong type of data selected. Please select either SPECFEM or GEMINI")

            trace.stats.network = self.networkCode
            trace.stats.station = self.name
            trace.stats.channel = comp.upper()
            trace.stats.starttime = event.recordStart
            trace.stats.delta = self.dt
            trace.stats.sampling_rate = 1 / self.dt
            if singleStation:
                self.stream.append(trace)
            else:
                stream.append(trace)

    def writeStationStreamData(self, outfile):
        """
        Methode to store the profile data as obspy stream object and prepare the coordinates of the
        points in the profile such that they can be written to a vtk file for plotting

        :param event: object of type seismicEvent
        :param outfile: Full path to where the data is to be saved
        """
        self.stream.write(outfile, format='PICKLE')

    def calculateDistaz(self, event):
        """
        Calculate the distance, azimuth and backazimuth for the station for the observed event.

        :param event: Object of type seismicEvent
        :return: Adds the distance to the event, the stations azimuth and backazimuth to the object.
        """
        client = Client()
        result = client.distaz(stalat=self.lat, stalon=self.lon, evtlat=event.lat, evtlon=event.lon)
        self.distance = result["distance"]
        self.azimuth = result["azimuth"]
        self.backazimuth = result["backazimuth"]

    def calculateArrivalTimes(self, event, model, phaseList=[]):
        """
        Calculate the arrival time of various Phases and sore them for later use.

        :param event: seismicEvent Object
        :param model: TauPy velocity model.
        :param phaseList: A list of phases for which the arrival times are to be calculated.
        :return:
        """

        self.calculateDistaz(event)
        recDepth = 0
        if not phaseList:
            self.arrivals = model.get_travel_times(source_depth_in_km=event.depth, distance_in_degree=self.distance,
                                                   receiver_depth_in_km=recDepth)
        else:
            self.arrivals = model.get_travel_times(source_depth_in_km=event.depth, distance_in_degree=self.distance,
                                                   phase_list=phaseList, receiver_depth_in_km=recDepth)

    def calculatePropagationDirection(self, event):
        """
        Calculate propagation direction at receiver

        :param: event: seismic_event object
        :return: self.propdir: propagation direction of wave at station (from south over east)
        """

        # Convert Global coordinates given in degree to and rotate latitude by 180 degree to make rotation into
        # different coordinate systems easier.
        sta_colat = 0.5 * np.pi - np.deg2rad(self.lat)
        sta_lon = np.deg2rad(self.lon)
        ev_colat = 0.5 * np.pi - (np.deg2rad(event.lat))
        ev_lon =np.deg2rad(event.lon)
        # Calculate propagation direction. The epicentral distance is calculated as a byproduct.
        delta, propdir = utils.geo2epi(ev_colat, ev_lon, sta_colat, sta_lon)      # seen from station
        self.propdir = propdir - np.pi # propagation direction

    def calculateEpicentralCoordinates(self, event):
        """
        Calculate the epicentral coordinates from station and event information.

        :param event: object of type seismicEvent
        :return:
        """
        self.r, self.delta, self.xi = utils.calculateEpicentralCoordinates(self.r, self.lat, self.lon, event)
        # Convert to rad for later use.
        self.delta = np.deg2rad(self.delta)
        self.xi = np.deg2rad(self.xi)


    def rotateToBox(self, event, boxCenterLat, boxCenterLon):
        """
        Rotate the Gemini ZNE data to Specfem XYZ data

        :param event: object of type seismicEvent
        :param boxCenterLat, boxCenterLon: Lat, Lon Coordinates of the center of the Specfem box
        :return:
        """
        # Convert data dictionary to numpy array to rotate the data to a new coordinate system
        tmp = np.asarray([item for item in self.gemini_data.values()])

        uz, un, ue = tmp[0], tmp[1], tmp[2]

        # Convert Global coordinates given in degree to and rotate latitude by 180 degree to make rotation into
        # different coordinate systems easier.
        thetas = np.pi / 2 - np.deg2rad(event.lat)
        phis = np.deg2rad(event.lon)
        thetac = np.pi / 2 - np.deg2rad(boxCenterLat)
        phic = np.deg2rad(boxCenterLon)

        tmp[0], tmp[1], tmp[2] = ar.vectorRLTfromZNE(self.propdir, uz, un, ue)
        uc = ar.vectorLCfromLS(self.delta, self.xi, tmp)
        ug = ar.vectorGCfromLC(thetas, phis, uc)
        ub = ar.vectorLCfromGC(thetac, phic, ug)
        tmp = ar.vectorRCfromLC(0.5 * np.pi, ub)

        # Fill the data dictionary with the rotated data
        self.gemini_data = dict(zip(["BXX", "BXY", "BXZ"], tmp))

    def rotateToZNE(self, event, boxCenterLat, boxCenterLon):
        """
        Rotate the Specfem-Data from the Box cartesian system to the global ZNE system.

        :param event: object of type seismicEvent
        :param boxCenterLat, boxCenterLon: Lat, Lon Coordinates of the center of the Specfem box
        :return:
        """

        # Convert Global coordinates given in degree to and rotate latitude by 180 degree to make rotation into
        # different coordinate systems easier.
        thetas = np.pi / 2 - np.deg2rad(event.lat)
        phis = np.deg2rad(event.lon)
        thetac = np.pi / 2 - np.deg2rad(boxCenterLat)
        phic = np.deg2rad(boxCenterLon)

        # Convert data dictionary to numpy array to rotate the data to a new coordinate system
        tmp = np.asarray([item for item in self.specfem_data.values()])

        ul = ar.vectorLCfromRC(np.pi / 2., tmp)
        ug = ar.vectorGCfromLC(thetac, phic, ul)
        uc = ar.vectorLCfromGC(thetas, phis, ug)
        tmp = ar.vectorLSfromLC(self.delta, self.xi, uc)

        ur, ul, ut = tmp[0], tmp[1], tmp[2]
        tmp[0], tmp[1], tmp[2] = ar.vectorZNEfromRLT(self.propdir, ur, ul, ut)

        # Fill the data dictionary with the rotated data
        self.specfem_data = dict(zip(self.components, tmp))

    def plotIndividualStationcomparison(self, event, folder, rotate, rotateToBox):
        """

        :param event: object of type seismicEvent
        :param folder: Where to save the plot to
        :param rotate, rotateToBox: Bools to select a specific type of rotation
        :return:
        """
        rotDict = {True: "rot", False: "unrot"}
        boxDict = {True: "toBox", False: "fromBox"}
        offset = max(map(max, self.dataset)) * 2
        minimum = min(map(min, self.dataset))
        fig = plt.figure()
        ax = fig.add_axes([0.1, 0.1, 0.85, 0.85])
        ax.set_ylim(1.5 * minimum, offset * 2.5)
        ax.set_yticks([1, offset, offset * 2])
        if rotateToBox:
            ax.set_yticklabels(["BXX", "BXY", "BXZ"])
        else:
            ax.set_yticklabels(["BXZ", "BXN", "BXE"])
        ax.set_title("Seismograms recorded at " + self.name)
        ax.set_xlabel("Seconds after " + event.originTimeUTC.isoformat())
        ax.set_ylim(minimum * 1.35, offset * 2.6)
        for phase in self.arrivals:
            if phase.time <= self.recordLength:
                wave = phase.time
                ax.add_line(Line2D([wave, wave], [minimum * 1.05, offset * 2.5], color='b', linewidth=1, ls='--'))
                ax.annotate(phase.name, xy=(wave - 5, 1.2 * minimum), xytext=(wave - 5, 1.2 * minimum), color='b',
                            rotation=90, fontsize=5)
        gemini_time = [t * self.dt + event.cutoffTime for t in range(len(self.gemini_data[0]))]
        specfem_time = [t * 0.03 + event.cutoffTime for t in range(len(self.specfem_data[0]))]
        for i in range(3):
            if i == 0:
                ax.plot(gemini_time, self.gemini_data[i, :] + (offset * i), c="red", lw=1, label="Gemini")
                ax.plot(specfem_time, self.specfem_data[i] + (offset * i), c="k", lw=1, label="Specfem")
            else:
                ax.plot(gemini_time, self.gemini_data[i, :] + (offset * i), c="red", lw=1)
                ax.plot(specfem_time, self.specfem_data[i] + (offset * i), c="k", lw=1)
        leg = ax.legend(loc="lower left")
        filename = self.basename + "_plot_" + rotDict[rotate] + "_" + event.eventId
        if rotate:
            filename = filename + "_" + boxDict[rotateToBox]
        filename += ".png"
        plt.savefig(folder + "comparison/" + filename, dpi=300)
        #plt.show()
        plt.close(fig)

    def calculateBoxCoordinates(self, boxCenterLat, boxCenterLon):
        """
        Calculate the cartesian coordinates for the station in the specfem box.

        :param boxCenterLat, boxCenterLon: Lat, Lon Coordinates of the center of the Specfem box
        :return:
        """

        r = 6371000. + self.elevation
        sta_colat = 0.5 * np.pi - np.deg2rad(self.lat)
        sta_lon = np.deg2rad(self.lon)
        xg, yg, zg = ar.coordinatesLCfromLS(r, sta_colat, sta_lon)
        xc, yc, zc = ar.coordinatesLCfromGC(0.5 * np.pi - np.deg2rad(boxCenterLat), np.deg2rad(boxCenterLon), xg, yg, zg)
        self.x, self.y, self.z = ar.coordinatesRCfromLC(np.pi / 2, xc, yc, zc)



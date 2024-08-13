from obspy import UTCDateTime
import numpy as np
import utils
from obspy import read, Trace
from seismicStation import SeismicStation


class SeismicEvent(object):
    """
    Container for seismic event information
    """

    component = ["z", "n", "e"]

    def __init__(self,eventId = "",csys = "",sourceType = 1,location= [],magnitude = 5,
                 originTime = "20191207_000000_000000",excitation = [], cutoffTimeFront = 0.):
        self.eventId = eventId
        self.csys = csys
        self.sourceType = sourceType
        self.location = location
        if not self.location == []:
            self.lat = self.location[0]
            self.lon = self.location[1]
            self.depth = self.location[2]
        self.magnitude = magnitude
        self.originTime = originTime
        self.excitation = excitation
        self.originTimeUTC = self.__convertToUTC()
        self.recordStart = self.__convertToUTC() + cutoffTimeFront
        self.cutoffTime = cutoffTimeFront

    def __convertToUTC(self):
        """
        Method to convert to originTime string to a UTCDateTime object
        :return: UTCDateTime version originTime of the event
        """
        return UTCDateTime(int(self.originTime[:4]), int(self.originTime[4:6]), int(self.originTime[6:8]),
                           int(self.originTime[9:11]), int(self.originTime[11:13]), int(self.originTime[13:15]),
                           int(self.originTime[16:]))

    def readFromFile(self, eventData, cutSeismo = True):
        """
        Read Eventparameter from input file

        :param eventData: Line from the input file containing the event data
        """

        if len(eventData) >= 6:
            self.originTime = eventData[1]
            self.eventId = eventData[0]
            self.magnitude = eventData[5]
            self.lat = float(eventData[2])
            self.lon = float(eventData[3])
            self.depth = float(eventData[4])
            self.cutoffTime = float(eventData[6])
            self.recordStart += self.cutoffTime
        if len(eventData) > 8:
            self.sourceType = eventData[7]
            self.excitation = eventData[8:-1]
        if len(eventData) == 7 or len(eventData) == 14 or len(eventData) == 11:
            self.recordStart = self.__convertToUTC()


    def readFromCMT(self, filename):
     """
     Reading event parameters from a CMT type file

     :param filename: filename containing CMT data
     """

     with open(filename,'r') as f:
         for line in f:
             a = line.split()
             if(a[0] == 'event'):
                 self.eventId = a[2]
             elif(a[0] == 'latitude:'):
                 self.lat = float(a[1])
             elif(a[0] == 'longitude:'):
                 self.lon = float(a[1])
             elif(a[0] == 'depth:'):
                 self.depth = float(a[1])

    def readCoordSystem(self, eventData):
        """
        :param eventData: In this special case this is the first line in the event List
        :return:
        """
        self.csys = eventData[0]

    def calculateCutoffTime(self, box, model, buffer):
        """
        Method to define a start time at which the synthetic seismogram for the current event is written to the hdf file
        :param box : object containing parameter for the specfem box
        :param model: TauPy velocity model.
        :param buffer: Since the seismogram should not start directly with the waveform a buffer added.
        :param phaseList: The arrival times of those phases will be calculated.
        :return:
        """
        self.__rDeltaXi = []
        for r, lat, lon in zip(box.boundaryPoints[0], box.boundaryPoints[1], box.boundaryPoints[2]):
            self.__rDeltaXi.append(utils.calculateEpicentralCoordinates(r, lat, lon, self))

        minDist = 180.  # The allowed max angle is 180 degrees, as that is directly opposite on the Earth
        minR = 6371000.  # Since the closest receiver is always subterran, we may start the search at 0.

        for r, delta, xi in self.__rDeltaXi:
            if delta <= minDist and r <= minR:
                minDist = delta
                minR = r
        maxDepth = 6371000. - minR

        self.arrivals = model.get_travel_times(source_depth_in_km=self.depth/1000, distance_in_degree=minDist,
                                               receiver_depth_in_km=maxDepth/1000)
        pwave = self.arrivals[0]
        self.cutoffTime = pwave.time - buffer

    def createStreamData(self, boxCenterLat, boxCenterLon, stationFile, geminiFile = "", dataFolder = "", rotate = True, datatype = "displacement", typ = "SPECFEM"):
        """
        :param boxCenterLon, boxCenterLat: Center coordinates of the SPECFEM box.
        :param box: Object of Type SpecfemBox
        :param stationFile: File containing a list of stations and corresponding infos. Currently based on the Gemini version
        :param geminiFile: File containing the synthetic Seismograms calculated with GEMINI. Only needed if type == GEMINI
        :param dataFolder: Path to the Seismograms generated by Specfem. Only needed if type == SPECFEM
        :param rotate: Depends on the type of data: True for SPECFEM means rotate to ZNE, for GEMINI to XYZ.
        :param datatype: Specfiy wich datatype is selected: inputs are "velocity" or "displacement"
        :param typ: Specify which type of data is read (GEMINI or SPECFEM)
        :return:
        """

        stream = read()
        stream.clear()

        #Select filename for the ouput file based on data type
        if typ == "SPECFEM":
            outfile = "specfemData" + self.eventId + ".dat"
        elif typ == "GEMINI":
            outfile = "geminiData" + self.eventId + ".dat"

        # Go over all Stations in the stations-file
        with open(stationFile, "r") as infile:
            for line in infile:
                stationData = line.split()
                if stationData[0] == "S" or stationData[0] == "C":
                    continue
                else:
                    # Create a Station object and read basic information from station file
                    station = SeismicStation()
                    station.readStationData(stationData)
                    #print("Processing Station: " + station.name + " at location: ", station.lat, station.lon,
                    #      station.elevation)
                    # Calculate Propagation direction and Epicentral coordinates in case rotation is needed.
                    station.calculatePropagationDirection(self)
                    station.calculateEpicentralCoordinates(self)

                    # Read and rotate Data according to input type
                    if typ == "SPECFEM":
                        station.readSpecfemData(self, dataFolder, datatype)
                        if rotate:
                            station.rotateToZNE(self, boxCenterLat, boxCenterLon)
                    elif typ == "GEMINI":
                        station.readGeminiData(self, geminiFile)
                        if rotate:
                            station.rotateToBox(self, boxCenterLat, boxCenterLon)
                    else:
                        raise ValueError("Wrong datatype entered. Enter GEMINI or SPECFEM.")
                    # Add station data to event stream
                    station.addTraces(self, stream, typ)

        print("writing file: " + dataFolder + outfile)
        # Write event stream
        stream.write(dataFolder + outfile, format="PICKLE")

    def createAxisemStreamData(self, box, folder, stationFile, axisemlocation, pylot, rotate = False):
        stream = read()
        stream.clear()
        outfile = "axisemStreamData.dat"
        with open(folder + stationFile) as stations:
            for line in stations:
                stationInfo = line.split(" ")
                for i in range(3):
                    comp = SeismicEvent.component[i]
                    filename = axisemlocation + stationInfo[0] + "_" + stationInfo[1] + "_disp_post_mij_conv0000_" + comp.upper() + ".dat"
                    tmp = np.loadtxt(folder + filename).transpose()
                    tr = Trace()
                    tr.data = tmp[1]
                    tr.stats.network = stationInfo[1]
                    tr.stats.station = stationInfo[0]
                    tr.stats.channel = comp.upper()
                    tr.stats.starttime = self.recordStart
                    tr.stats.delta = abs(tmp[0][0] - tmp[0][1])
                    tr.stats.sampling_rate = 1 / abs(tmp[0][0] - tmp[0][1])
                    stream.append(tr)
        print("writing file: " + folder + axisemlocation + pylot + outfile)
        stream.write(folder + axisemlocation + pylot + outfile, format = "PICKLE")

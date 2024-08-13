import h5py
import numpy as np
from seismicEvent import SeismicEvent
#For Anaconda users import the following command:
from pyevtk.hl import pointsToVTK



class BoundaryWavefield(object):
    """
    Container for boundary data involved in coupling SPECFEM to GEMINI
    for a selected event.
    """
    def __init__(self,hdffile,evkey,si=False):
        """ 
        Class instance initialized with valid HDF file.
        Opens HDF file. Keeps file ID. Reads attributes
        and other essential information about the dataset

        :param hdffile: Name of the hdffile containing tha data
        :param evkey: EventID
        :param profile (optinal): Specify what kind of profile is to be extracted from the dataset
                                  Options are: horizontal or vertical
        :param side (optinal): Specifies the side of the boundary box on which the profile is selected.
                               Options are: left, right, front, back
        """
        self.hdffile = hdffile
        self.fid =  h5py.File(hdffile,'r')
        self.group = self.fid[evkey]
        [self.dt,self.tlen] = self.group.attrs["samplingIntervalLength"]
        #self.dt = self.fid.attrs["samplingInterval"]
        self.event = self.readEvent(evkey)
        self.nbptot = self.getNumBoundaryPoints(evkey)
        self.start = (0, self.event.originTimeUTC)
        self.si_b = si
        if si:
            self.si = 1000.
        else:
            self.si = 1.

    def readEvent(self,evkey):
        """
        Return SeismicEvent object associated with given key in file
        """
        csys = self.group.attrs["Coordinate_system"]
        location = self.group.attrs["Coordinates"]
        originTime = self.group.attrs["Date_and_time"]
        magnitude = self.group.attrs["Magnitude"]
        styp = self.group.attrs["Source_type"]
        cutoffTimeFront = self.group.attrs.get("timing")[-1]
        if styp == 1:
            excitation = self.group.attrs["Moment"]
        else:
            exitation = self.group.attrs["Force"]
        return SeismicEvent(evkey,csys,styp,location,magnitude,
                            originTime,excitation,cutoffTimeFront)

    def getNumBoundaryPoints(self,evkey):
        """
        Get number of boundary points involved for event evkey
        """
        nbptot = 0
        for proc in self.group.values():                                     # get total number of boundary points
            nb = proc.attrs["numBoundaryPointsPerProc"]
            nbptot = nbptot+np.sum(nb)
        return nbptot        

    def readCoordsVelocityTraction(self,evkey,timeStep,vtkbasename):
        """
        Read coordinates, velocity and traction data 
        for a given event and time step and write to VTK.
        Collect data of all processes. Append time step index
        to VTK base name.
        """
        veltrac = np.zeros((6,self.nbptot))                                # allocate a big np.array for values of all procs
        coords = np.zeros((3,self.nbptot))
        ja = 0
        jb = 0
        for proc in self.group.values():                                   # collect data from processes
            [nbpvert,nbpbot] = proc.attrs["numBoundaryPointsPerProc"]      # number of points per process
            nbp = nbpvert+nbpbot
            print("read data of process: ",proc.name,nbp)
            jb = jb+nbp
            p = proc["coordinates"]
            coords[:,ja:jb] = np.transpose(p[0:nbp,0:3])
            p = proc["velocityTraction"]
            veltrac[:,ja:jb] = np.transpose(p[timeStep,0:nbp,:])
            ja = ja+nbp
        return (coords,veltrac)

    def setStartTime(self, buffer, model):
        """
        Sets a starttime for the seismogram to remove unwanted data for forward simulation with specfem

        :param buffer: padding time for the start time of the seismogram
        :param model: TauPy velocity model
        :param phaseList: A list of Phase Names. Arrival times for these phases are extracted from the TauPy model
        sets: (startIndex, startTime): startTime: UTCDate time object of the start time of the seismogram, startIndex: Corresponding index for the data array stored in the bw object.
        """
        #Testabschnitt
        # Der abschnitt wird noch darfurch ersetzt, dass diese informationen in Zukunft mit dem hdf5 file mitgeliefter werden.
        minDist = 180.  # The allowed max angle is 180 degrees, as that is directly opposite on the Earth
        maxDepth = 0.  # Since the closest receiver is always subterran, we may start the search at 0.
        for proc in self.group.values():
            coordinates = proc["coordinates"]
            for coords in coordinates:
                recDepth = 6371.*self.si - coords[3]
                if coords[-2] <= minDist and recDepth > maxDepth:
                    minDist = coords[-2]
                    maxDepth = recDepth

        print(minDist, maxDepth)
        # ende Testabschnitt

        arrivals = model.get_travel_times(source_depth_in_km=self.event.depth, distance_in_degree=minDist,
                                          receiver_depth_in_km=maxDepth)
        print(arrivals)
        pwave = arrivals[0]
        # To get the right portion from the data-array stored in the hd5 file an index at the desired cutoff time is needed.
        # pwave.time = P-wave arrival time; 50s buffer, divide by dt to get the approximate position in the original array
        cutoff = pwave.time - buffer
        startIndex = int(cutoff / self.dt)
        startTime = self.event.originTimeUTC + cutoff
        self.start = (startIndex, startTime)

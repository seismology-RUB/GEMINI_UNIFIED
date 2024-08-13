#------------------------------------------------------------------------------------
#  Python class for working with rays and travel times in spherically symmetric media
#  Ray table is organized as follows:
#     raytable[slowness index,radial node index,delta_or_tt index]
#  For example:
#  travel time for slowness 7 and radial node 5 is tt = raytable[7,5,0]
#  epicentral distance for slowness 7 and radial node 5 is delta = raytable[7,5,1]
#------------------------------------------------------------------------------------
import h5py
import numpy as np
import axesRotation as ar

class basicRay:
    '''
       Basic ray object with:
       p: Slowness
       tt, delta, rad: numpy arrays with travel time, epicentral distance and radii of ray points
    '''
    def __init__(self,p,tt,rad,delta):
        self.p = p
        self.tt = tt
        self.delta = delta
        self.rad = rad

class epicentralRay(basicRay):
    '''
       Ray object in epicentral coordinates:
       Input: basic ray: br
       xi: azimuth of ray at source (counted ccw from south over east)
       lons: geographical longitude of source
       lats: geographical latitude of source
    '''
    def __init__(self,br,xi,lats,lons):
        self.lats = lats
        self.lons = lons
        self.xi = xi
        super().__init__(br.p,br.tt,br.rad,br.delta)

class geographicalRay(basicRay):
    '''
       Ray object in geographical coordinates:
       Input basic ray: br
       xi: azimuth of ray at source (counted ccw from south over east)
       lons: geographical longitude of source
       lats: geographical latitude of source
    '''
    def __init__(self,br,xi,lats,lons):
        super().__init__(br.p,br.tt,br.rad,br.delta)
        xxi = xi*np.ones(np.size(br.delta))
        xl,yl,zl = ar.coordinatesLCfromLS(br.rad, br.delta, xxi)
        xg,yg,zg = ar.coordinatesGCfromLC(0.5*np.pi-lats, lons, xl, yl, zl)
        r,theta,self.lon = ar.coordinatesLSfromLC(xg, yg, zg)
        self.lat = 0.5*np.pi-theta

class cartesianRay(basicRay):
    '''
       Ray object in global Cartesian coordinates:
       Input: basic ray: br
       xi: azimuth of ray at source (counted ccw from south over east)
       lons: geographical longitude of source
       lats: geographical latitude of source
    '''
    def __init__(self,br,xi,lats,lons):
        super().__init__(br.p,br.tt,br.rad,br.delta)
        xxi = xi*np.ones(np.size(br.delta))
        xl,yl,zl = ar.coordinatesLCfromLS(br.rad, br.delta, xxi)
        self.x,self.y,self.z = ar.coordinatesGCfromLC(0.5*np.pi-lats, lons, xl, yl, zl)

class rayTable:

    def __init__(self,ray_table_file,logfile = "ray.log"):
        '''
            Read ray table and fill array
        '''        
        self.fid = h5py.File(ray_table_file,'r')
        self.slowness = np.array(self.fid['rayParameters'])
        self.rtp = np.array(self.fid['turningPointRadii'])
        self.rnod = np.array(self.fid['receiverRadii'])
        self.raytable = np.array(self.fid['rayTable'])
        self.ntp = np.size(self.rtp)
        self.nnod = np.size(self.rnod)
        self.drt = self.rtp[1]-self.rtp[0]
        self.log = open(logfile,"w")

    def getBasicRay(self,delta,rs,re):
        '''
           Get a basic ray, defined by its travel time and slowness for given delta, source and receiver radius.
           Return 2D array ray[0:1,0:npray] containing epicentral distance (rad) and radii (m).
           Also returns ray parameters [slowness (s/rad), delta (rad), tt (s).
        '''
        """
           Find indices in rnod-array closest to rs and re
        """
        js = np.abs(self.rnod-rs).argmin()
        je = np.abs(self.rnod-re).argmin()
        self.log.write("Location of source:   rs = {0:12f} js = {1:5d} rnod[js] = {2:12f}\n".format(rs,js,self.rnod[js]))
        self.log.write("Location of receiver: re = {0:12f} je = {1:5d} rnod[je] = {2:12f}\n".format(re,je,self.rnod[je]))
        """
           Calculate delta versus slowness and find slowness index closest to delta
        """
        delta_vs_p = self.raytable[:,js,0]+self.raytable[:,je,0]
        if delta < np.min(delta_vs_p):
            raise Exception("Epicentral distance too small: you are in the range of triplications")
        jsl = np.abs(delta_vs_p-delta).argmin()
        rturn = self.rtp[jsl]
        p = self.slowness[jsl]
        deltos = self.raytable[jsl,js,0]
        tttos = self.raytable[jsl,js,1]
        self.log.write("Slowness and turning radius: jsl = {0:5d} p = {1:12f} rturn = {2:12f}\n".format(jsl,p,rturn))
        self.log.write("Distance in and out: delta_in = {0:12f} delta_out = {1:12f}\n".format(delta,delta_vs_p[jsl]))
        """
           Treat diffraction case when delta is larger than for CMB as turning radius
        """
        if delta > delta_vs_p[0]:
            deldiff = delta-delta_vs_p[0]
            self.log.write("Diffraction: distance overshoot: {0:12f} \n".format(deldiff))
            diffraction = True
        else:
            diffraction = False
            deldiff = 0.0
        """
           find index in rnod-array just above turning radius
        """
        jt = np.min(np.where(self.rnod > rturn))
        self.log.write("Location of turning point: rturn = {0:12f} jt = {1:5d} rnod[jt] = {2:12f}\n".format(rturn,jt,self.rnod[jt]))
        """
           Compose ray from points rnod[js:jt-1:-1] plus rtp plus rnod[jt:je+1]
           Total number of ray points: js-jt+1 + 1 + je+1-jt = js+je-2*jt+3
           Start with delta=0 at source
           Consider additional points in case of diffraction
           We add nadd intervals but do not add the end point 
           because it is the first point of the upgoing ray.
        """
        npray = js+je-2*jt+3
        if diffraction:
            nadd = int(np.ceil(deldiff*rturn/self.drt))
            delstep = deldiff/nadd
            npray = npray+nadd-1 
        """
           allocate space for ray(delta,rad,tt)
        """
        delr = np.zeros(npray)
        rad = np.zeros(npray)
        tt = np.zeros(npray)
        self.log.write("Generate: {0:5d} points on ray\n\n".format(npray))
        """ 
           ray on source side
        """
        jend = js-jt+1
        delr[0:jend] = deltos-self.raytable[jsl,js:jt-1:-1,0]
        rad[0:jend] = self.rnod[js:jt-1:-1]
        tt[0:jend] = tttos-self.raytable[jsl,js:jt-1:-1,1]
        """ 
           ray at turning point 
        """
        delr[jend] = deltos
        rad[jend] = rturn
        tt[jend] = tttos
        """
           in case of a diffraction, add overshoot distance to ray
        """
        if diffraction:
           delr[jend+1:jend+nadd] = [deltos+j*delstep for j in range(1,nadd)]
           rad[jend+1:jend+nadd] = rturn
           tt[jend+1:jend+nadd] = [tttos+j*delstep*p for j in range(1,nadd)]
           jend = jend+nadd-1
        """ 
           ray on receiver side 
        """
        delr[jend+1:npray] = deltos+self.raytable[jsl,jt:je+1,0]+deldiff
        rad[jend+1:npray] = self.rnod[jt:je+1]
        tt[jend+1:npray] = tttos+self.raytable[jsl,jt:je+1,1]+deldiff*p

        return basicRay(p,tt,rad,delr)


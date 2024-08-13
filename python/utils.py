#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import axesRotation as ar

def mapCartesianBoxToSphericalChunk(x,y,z, si = False):
    """
        :param x,y,z: cartesian coordinates of point in cartesian box on input
        :return xs,ys,zs: cartesian coordinates of point in spherial chunk on output
    """
    if si:
        R_EARTH = 6371000.
    else:
        R_EARTH = 6371.
    r = R_EARTH + z                             # true distance of point from center of sphere
    lat = y/R_EARTH                             # y interpreted as arclength along 0-degree meridian
    ys = r*np.sin(lat)                          # new y-value in spherical chunk
    rp = r*np.cos(lat)                          # radius of parallel at latiude = lat

    lon = x/R_EARTH                             # x interpreted as arclength along equator
    xs = rp*np.sin(lon)                         # use radius of parallel at latitude lat
    zs = np.sqrt(r**2 - xs**2 - ys**2)-R_EARTH

    return xs, ys, zs

def mapSphericalChunkToCartesianBox(xs,ys,zs,si=False):
    """
        map spherical chunk to cartesian box (MD, WF)
        assumes that origin is at the equator at lon=0 and lat=0 and at surface

        :param  xs,ys,zs: cartesian coordinates of point in spherial chunk on input
        :return x,y,z: cartesian coordinates of point in cartesian box on output
    """

    if si:
        R_EARTH = 6371000.
    else:
        R_EARTH = 6371.
    r = np.sqrt((zs+R_EARTH)**2+xs**2+ys**2)
    lat = np.arcsin(ys/r)
    rp = r*np.cos(lat)
    lon = np.arcsin(xs/rp)
    x = R_EARTH*lon
    y = R_EARTH*lat
    z = r-R_EARTH

    return x,y,z

def geo2epi(thgeo,phigeo,thpol,phipol):
    """
    computes epicentral coordinates of a point on the unit
    sphere with respect to pole thpol,phipol that has
    geographical coordinates thgeo, phigeo
    :param thgeo, phigeo: geographical coordinates [rad]
    :param :
    :param thpol:
    :param phipol:
    :param theta:
    :param phi:
    :return: theta between 0 and pi and phi between 0 and 2*pi
    """

    cthpol=np.cos(thpol)
    sthpol=np.sin(thpol)
    cthgeo=np.cos(thgeo)
    sthgeo=np.sin(thgeo)
    ctheta=cthpol*cthgeo+sthpol*sthgeo*np.cos(phigeo-phipol)
    if(abs(ctheta) > 1.):
        ctheta = sign(1.,ctheta)
    theta=np.arccos(ctheta)
    if(theta < 1.-4*np.pi or theta > 0.999*np.pi):
        phi=phipol
        return theta, phi
    stheta=np.sin(theta)
    sa=sthgeo*np.sin(phigeo-phipol)/stheta
    ca=(cthgeo-ctheta*cthpol)/(stheta*sthpol)
    if(abs(sa) > 1.):
        sa=sign(1.,sa)

    alfa=np.arcsin(sa)
    if(ca >= 0.):
        phi=alfa
    if(ca < 0.):
        if(sa>=0.):
            phi=np.pi-alfa
        if(sa < 0.):
            phi=-np.pi-alfa
    phi=np.pi-phi
    return theta, phi

def calculateEpicentralCoordinates(r, lat, lon , event):
    """
    Calculate Epicentral coordinates for a point on a sphere with regard too a given seismic event.


    :param r: Radius of a point or station
    :param lat: Latitude of a point or station (in degree)
    :param lon: Longitude of a point or station (in degree)
    :param event: object of type seismic event
    :return: r, delta, xi (in degree)
    """

    #Convert the geographical coordinates from degree top rad and rotate latitude.
    sta_colat = np.pi / 2. - np.deg2rad(lat)
    sta_lon = np.deg2rad(lon)

    ev_colat = np.pi / 2 - np.deg2rad(event.lat)
    ev_lon = np.deg2rad(event.lon)

    xg, yg, zg = ar.coordinatesLCfromLS(r, sta_colat, sta_lon)
    x,y,z = ar.coordinatesLCfromGC(ev_colat, ev_lon, xg, yg, zg)
    r, delta, xi = ar.coordinatesLSfromLC(x,y,z)
    return r, delta, xi


def sign(a,b):
    """
    Implementation of the fortran function SIGN
    :param a: Int or float
    :param b: Int or float
    :return: The value of A with the sign of B.
    """
    if b < 0:
        return a*(-1)
    else:
        return a


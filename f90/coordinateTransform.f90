!----------------------------------------------------------------------------
!   Copyright 2016 Wolfgang Friederich (Ruhr-Universitaet Bochum, Germany)
!
!   This file is part of GEMINI_UNIFIED version 1.0.
!
!   GEMINI_UNIFIED version 1.0 is free software: you can redistribute it and/or modify
!   it under the terms of the GNU General Public License as published by
!   the Free Software Foundation, either version 2 of the License, or
!   (at your option) any later version.
!
!   GEMINI_UNIFIED version 1.0 is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU General Public License for more details.
!
!   You should have received a copy of the GNU General Public License
!   along with GEMINI_UNIFIED version 1.0.  If not, see <http://www.gnu.org/licenses/>.
!-------------------------------------------------------------------------------
!  Module for coordinate transformations
!-------------------------------------------------------------------------------
 module coordinateTransform
    use mathConstants
    implicit none
!
 contains
!--------------------------------------------------------------------------------------
!> \brief Transformation matrix from local to global cartesian basis vectors on a sphere
!! Local cartesian CS has 3-axis through point on sphere and 1-axis pointing through
!! meridian of point at 90 degrees to 3-axis.
!! Global system has 3-axis through north pole and 1-axis through Greenwich meridian
!! Global system basis vectors: e_i
!! Local system basis vectors:  k_n
!! Transformation: e_i = T_in k_n
!> \param lat Latitude in rad of point on sphere
!> \param lon Longitude in rad of point on sphere
!
    function localToGlobalCartesianOnSphereCoordinateTransform(lat,lon) result(tm)
    real, dimension(:,:), pointer :: tm
    real :: lon,lat
    double precision :: colat,ct,st,cf,sf
    colat = 0.5*mc_pid-dble(lat) 
    ct = dcos(colat)
    st = dsin(colat)
    cf = dcos(dble(lon))
    sf = dsin(dble(lon))
    allocate(tm(3,3))
    tm(1,1)=ct*cf
    tm(1,2)=-sf
    tm(1,3)=st*cf
    tm(2,1)=ct*sf
    tm(2,2)=cf
    tm(2,3)=st*sf
    tm(3,1)=-st
    tm(3,2)=0.
    tm(3,3)=ct
    end function localToGlobalCartesianOnSphereCoordinateTransform
!-----------------------------------------------------------------------------------------
!> \brief Transformation matrix from local spherical basis vectors to global cartesian ones
!! Local basis vectors at some point on sphere with coordinates (lat,lon) in rad
!! Global cartesian basis vectors with 3-axis through north pole and 1-axis through Greenwich meridian
!! Local spherical basis vectors with 1-axis pointing vertical at point, 2-axis tangential south
!! Global system basis vectors: e_k
!! Local system basis vectors:  f_s
!! Transformation: e_k = T_ks f_s
!! This is essentially the same transform as localToGlobalCartesianCoordinateTransform but with
!! the order of basis vectors changed
!> \param lat Latitude in rad of point on sphere
!> \param lon Longitude in rad of point on sphere
!
    function localSphericalToGlobalCartesianOnSphereCoordinateTransform(lat,lon) result(tm)
    real, dimension(:,:), pointer :: tm
    real :: lon,lat
    double precision :: colat,ct,st,cf,sf
    colat = 0.5*mc_pid-dble(lat) 
    ct = dcos(colat)
    st = dsin(colat)
    cf = dcos(dble(lon))
    sf = dsin(dble(lon))
    allocate(tm(3,3))
    tm(1,1) = st*cf
    tm(2,1) = st*sf
    tm(3,1) = ct
    tm(1,2) = ct*cf
    tm(2,2) = ct*sf
    tm(3,2) = -st
    tm(1,3) = -sf
    tm(2,3) = cf
    tm(3,3) = 0.
    end function localSphericalToGlobalCartesianOnSphereCoordinateTransform
!-----------------------------------------------------------------------------------------
!> \brief Transformation matrix from local rotated basis vectors to global unrotated ones
!! Local cartesian basis vectors rotated mathematically positive by angle phi in rad around z-axis
!! Local basis vectors with 1-axis = z = up, 2-axis = x = south and 3-axis = y = east
!! Global system basis vectors: e_k
!! Local system basis vectors:  f_s
!! Transformation: e_k = T_ks f_s
!> \param phi Rotation of local basis vectors around 1-axis in rad
!
    function localToGlobalCartesianInPlaneCoordinateTransform(phi) result(tm)
    real, dimension(:,:), pointer :: tm
    real :: phi,cf,sf
    cf = cos(phi)
    sf = sin(phi)
    allocate(tm(3,3))
    tm(1,1) = 0.
    tm(2,1) = 0.
    tm(3,1) = 1.
    tm(1,2) = cf
    tm(2,2) = sf
    tm(3,2) = 0.
    tm(1,3) = -sf
    tm(2,3) = cf
    tm(3,3) = 0.
    end function localToGlobalCartesianInPlaneCoordinateTransform
!
 end module

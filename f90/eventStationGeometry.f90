! ===============================================================================
!  Module handling event-station geometry
! ===============================================================================
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
!---------------------------------------------------------------------
!   Compute epicentral distance, azimuth and propagation direction
!-----------------------------------------------------------------------------
module eventStationGeometry
    use seismicStation
    use seismicEvent
    use mathConstants
    implicit none
contains
!----------------------------------------------------------------------------------------------
!  calculate epicentral distance in km and azimuth in rad on sphere for given event and station
!  event:           seismic_event object
!  station:         seismic_station object
!  rearth:          earth radius in m
!  delta:           epicenral distance in m
!  phi:             azimuth of station as seen from event in rad
!
   subroutine deltaMeterPhiRadSphericalEventStationGeometry(event,station,rearth,delta,phi)
   type (seismic_event) :: event
   type (seismic_station) :: station
   double precision :: rearth,delta,phi
   double precision :: sta_colat,sta_lon,ev_colat,ev_lon
!
   sta_colat = 0.5*mc_pid-(.latrad.station)
   sta_lon = .lonrad.station
   ev_colat  = 0.5*mc_pid-(.slatrad.event)
   ev_lon  = .slonrad.event
   call geo2epi(sta_colat,sta_lon,ev_colat,ev_lon,delta,phi)            ! seen from event
   delta = delta*rearth                                                 ! convert to m
   end subroutine deltaMeterPhiRadSphericalEventStationGeometry
!---------------------------------------------------------------------------
!  calculate propagation direction at receiver
!  event:           seismic_event object
!  station:         seismic_station object
!  propdir:         propagation direction of wave at station (from south over east)
!
   subroutine propdirRadSphericalEventStationGeometry(event,station,propdir)
   type (seismic_event) :: event
   type (seismic_station) :: station
   double precision :: propdir
   double precision :: delta,sta_colat,sta_lon,ev_colat,ev_lon
!   
   sta_colat = 0.5*mc_pid-(.latrad.station)
   sta_lon = .lonrad.station
   ev_colat  = 0.5*mc_pid-(.slatrad.event)
   ev_lon  = .slonrad.event
   call geo2epi(ev_colat,ev_lon,sta_colat,sta_lon,delta,propdir)      ! seen from station
   propdir = propdir-mc_pid                                            ! propagation direction
   end subroutine propdirRadSphericalEventStationGeometry
!------------------------------------------------------------------------------------------------
!  calculate epicentral distance in meter and azimuth in rad on halfspace
!  for given event and station
!  event:           seismic_event object
!  station:         seismic_station object
!  delta:           epicenral distance in meter
!  phi:             azimuth of station as seen from event in rad
!
   subroutine deltaMeterPhiRadHalfspaceEventStationGeometry(event,station,delta,phi)
   type (seismic_event) :: event
   type (seismic_station) :: station
   double precision :: delta,phi
   double precision :: pythag
!
   delta = pythag(.lat.station-.slat.event,.lon.station-.slon.event)
   phi = atan2(.lon.station-.slon.event,.lat.station-.slat.event)
   end subroutine deltaMeterPhiRadHalfspaceEventStationGeometry
!
end module eventStationGeometry

!----------------------------------------------------------------------------
!   Copyright 2016 Wolfgang Friederich (Ruhr-Universitaet Bochum, Germany)
!
!   This file is part of GEMINI_UNIFIED version 1.0.
!
!   GEMINI_UNIFIED version 1.0 are free software:
!   you can redistribute it and/or modify
!   it under the terms of the GNU General Public License as published by
!   the Free Software Foundation, either version 2 of the License, or
!   (at your option) any later version.
!
!   GEMINI_UNIFIED version 1.0 are is distributed
!   in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU General Public License for more details.
!
!   You should have received a copy of the GNU General Public License
!   along with GEMINI_UNIFIED version 1.0.
!   If not, see <http://www.gnu.org/licenses/>.
!--------------------------------------------------------------
!  computes epicentral coordinates of a point on the unit
!  sphere with respect to pole thpol,phipol that has 
!  geographical coordinates thgeo, phigeo
!  returns theta between 0 and pi
!  and phi between 0 and 2*pi
!--------------------------------------------------------------
    subroutine geo2epi(thgeo,phigeo,thpol,phipol,theta,phi)
    double precision :: theta,phi,thpol,phipol,thgeo,phigeo
    double precision :: cthpol,sthpol,ctheta,stheta,cthgeo,sthgeo
    double precision :: pi,sa,ca,alfa,one
!
    pi=4.d0*datan(1.d0)
    one=1.d0
    cthpol=dcos(thpol)
    sthpol=dsin(thpol)
    cthgeo=dcos(thgeo)
    sthgeo=dsin(thgeo)
    ctheta=cthpol*cthgeo+sthpol*sthgeo*dcos(phigeo-phipol)
    if(abs(ctheta).gt.1.d0) ctheta=sign(1.d0,ctheta)
    theta=dacos(ctheta)
    if(theta.lt.1.d-4*pi.or.theta.gt.0.999d0*pi) then
        phi=phipol
        return
    endif
    stheta=dsin(theta)
    sa=sthgeo*dsin(phigeo-phipol)/stheta
    ca=(cthgeo-ctheta*cthpol)/(stheta*sthpol)
    if(abs(sa).gt.1.d0) then
        sa=sign(one,sa)
    endif
    alfa=dasin(sa)
    if(ca.ge.0.d0) phi=alfa
    if(ca.lt.0.d0) then
        if(sa.ge.0.d0) phi=pi-alfa
        if(sa.lt.0.d0) phi=-pi-alfa
    endif
    phi=pi-phi
!
    end subroutine geo2epi

! ====================================================================================
!  Calculate various wave number integrals using Green FK spectra
! ====================================================================================
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
!----------------------------------------------------------------------
!  This module provides routines for computing various wavenumber integrals
!  needed for the calculation of displacement spectra.
!----------------------------------------------------------------------
module wavenumberIntegrals
   use greenFKSpectra
   use besselFunctions
   use mathConstants
   implicit none
!                                 FORCE EXCITATION
!
!  Mapping from the 23 wavenumber integrals to the two frequency-wavenumber indices:
!
    integer, dimension(23,2) :: force_component_jump_from_wavenumber_integral = reshape( &
        & (/ 1,1,3,3,7,7,3, 5,5,6,6,9,9,6, 1,1,3,3,7,7,3, 2,2, &
        &    1,2,1,2,1,1,2, 1,2,1,2,1,1,2, 1,2,1,2,1,1,2, 1,2 /) &
        & ,(/ 23,2 /))
!                                 MOMENT EXCITATION
!
!  Mapping from the 46 wavenumber integrals to the two frequency-wavenumber indices:
!
    integer, dimension(46,2) :: moment_component_jump_from_wavenumber_integral = reshape( &
         & (/ 1,1,1,1,3,3,3,3,7,7,7,7,3,3, 5,5,5,5,6,6,6,6,9,9,9,9,6,6, 1,1,1,1,3,3,3,3,7,7,7,7,3,3, 2,2,2,2,  &
         &    1,2,3,4,1,2,3,4,1,2,1,2,3,4, 1,2,3,4,1,2,3,4,1,2,1,2,3,4, 1,2,3,4,1,2,3,4,1,2,1,2,3,4, 1,2,3,4 /) &
         & ,(/ 46,2 /))
!
!  Maximum wavenumber index depending on source type
!
    integer, dimension(0:1) :: maxIndexWavenumberIntegrals = (/23,46/)
!    
contains
!--------------------------------------------------------------------------
!  Compute values of Bessel functions at wavenumber points
!
!  nwn:         number of wabenumbers
!  dwn:         wavenumber spacing
!  x:           epicentral distance in km
!  besselj:     array with Bessel functions
!
    subroutine besselWavenumberIntegrals(nwn,dwn,x,besselj)
    integer :: nwn
    double precision :: dwn,x
    double precision, dimension(:,0:) :: besselj
    double precision :: wn
    integer :: j
!
    wn = dwn
    do j = 2,nwn
        besselj(j,0)=dbessj0(wn*x)
        besselj(j,1)=dbessj1(wn*x)
        besselj(j,2)=-besselj(j,0)+2./(wn*x)*besselj(j,1)
        wn = wn+dwn
    enddo
    end subroutine besselWavenumberIntegrals
!-----------------------------------------------------------
!  Compute wavenumber integrals for single force source
!
!  kint= 1:       R^2/(2*pi)*int ( dk k U^1 J_0 )  
!  kint= 2:       R^2/(2*pi)*int ( dk k (-2/kR) U^2 J_1 )  
!  kint= 3:       R^2/(2*pi)*int ( dk k (-kR) V^1 J_1 )  
!  kint= 4:       R^2/(2*pi)*int ( dk k (-2 V^2)(J_0-J_1/(kx)) )  
!  kint= 5:       R^2/(2*pi)*int ( dk k (-2 W^1) J_1/(kx) )  
!  kint= 6:       R^2/(2*pi)*int ( dk k (+2 W^1)(J_0-J_1/(kx)) )  
!  kint= 7:       R^2/(2*pi)*int ( dk k (+2 V^2) J_1/(kx) )
!
!  kint= 8:       R^2/(2*pi)*int ( dk k U^1_r J_0 )  
!  kint= 9:       R^2/(2*pi)*int ( dk k (-2/kR) U^2_r J_1 )  
!  kint=10:       R^2/(2*pi)*int ( dk k (-kR) V^1_r J_1 )  
!  kint=11:       R^2/(2*pi)*int ( dk k (-2 V^2_r)(J_0-J_1/(kx)) )  
!  kint=12:       R^2/(2*pi)*int ( dk k (-2 W^1_r) J_1/(kx) )  
!  kint=13:       R^2/(2*pi)*int ( dk k (+2 W^1_r)(J_0-J_1/(kx)) )  
!  kint=14:       R^2/(2*pi)*int ( dk k (+2 V^2_r) J_1/(kx) )
!
!  kint=15:       R^2/(2*pi)*int ( dk k U^1 (-kR) J_1 )  
!  kint=16:       R^2/(2*pi)*int ( dk k -2 U^2 (J_0-J_1/(kx)) )  
!  kint=17:       R^2/(2*pi)*int ( dk k -(kR)^2 V^1 (J_0-J_1/(kx)) )  
!  kint=18:       R^2/(2*pi)*int ( dk k (-2 V^2)(kR)(-J_1+J_2/(kx)) )  
!  kint=19:       R^2/(2*pi)*int ( dk k (+2 W^1)(kR) J_2/(kx) )  
!  kint=20:       R^2/(2*pi)*int ( dk k (+2 W^1)(kR)(-J_1+J_2/(kx)) )  
!  kint=21:       R^2/(2*pi)*int ( dk k (-2 V^2)(kR) J_2(kx)/(kx) )
!
!  kint=22:       R^2/(2*pi)*int ( dk k R^1 J_0 )  
!  kint=23:       R^2/(2*pi)*int ( dk k (-2/kR) R^2 J_1 )  
!
!  for one receiver at xx and one frequency using precomputed
!  values for the Bessel functions. This is possible because
!  the Bessel functions do not depend on frequency
!  and the wavenumber sampling is independent of frequency.
!  A complex value zdis is returned containing the
!  value of the desired integral.
!
!  Note: Unit of zdis is returned in meters/N
!
!  gfk:               greenFKSpectra object containing all meta information on Green FK spectra
!  gwnspr,gwnspi:     extracted wavenumber spectrum for selected frequency and wavenumber integral (double precision and imag)
!  kint:              index of requested wavenumber integral
!  xx:                epicentral distance in meters
!  besselj:           precomputed Bessel functions for all wavenumbers and desired distance
!  tapfrac:           fraction of wavenumber spectrum which is tapered
!  zdis:              complex value of wavenumber integral
!
    subroutine forceWavenumberIntegrals(gfk,gwnspr,gwnspi,kint,xx,besselj,tapfrac,zdis)
    type (green_fk_spectra) :: gfk
    real, dimension(:) :: gwnspr,gwnspi
    integer :: kint
    double precision :: xx,tapfrac
    double precision, dimension(:,0:) :: besselj
    double complex :: zdis
    double precision, dimension(:), allocatable :: wnbesj,greentap,wn
    double precision :: re,x,wn3,wne,wn3besy,xdk,wnbesy,by0,by1,by2
    integer :: j3,j33,j,nwn
    double complex :: zsum
!
!  number of wavenumbers
!
    nwn = size(gwnspr)
    re = gfk%rearth
!
!  guard against too small distances
!  if x is too small j3 may become larger than representable by integer*4
!  therefore choose x such that j3 = 1e9
!
    if(xx < 3.d0/(gfk%dwn*1.d9)) then
       x = 3.d0/(gfk%dwn*1.d9)
    else
        x = xx
    endif
!
!  find index of first wn-point greater than 3
!
    wn3=3.d0/x
    j3=int(wn3/gfk%dwn)+2
    wn3=(j3-1)*gfk%dwn
!
!  Bessel terms over whole interval
!  division by zero is unproblematic because we avoid wn=0
!  Bessel terms at k=0 vanish anyway
!  and x shouldn't be zero either
!
    allocate(wnbesj(nwn),wn(nwn))
    forall (j = 1:nwn) wn(j) = (j-1)*gfk%dwn
    select case (kint)
       case (1,8,22);  forall(j = 2:nwn) wnbesj(j)=wn(j)*besselj(j,0)
       case (2,9,23);  forall(j = 2:nwn) wnbesj(j)=-2*besselj(j,1)/re
       case (3,10,15); forall(j = 2:nwn) wnbesj(j)=-wn(j)*wn(j)*besselj(j,1)*re
       case (4,11,16); forall(j = 2:nwn) wnbesj(j)=-2.*wn(j)*( besselj(j,0)-besselj(j,1)/(wn(j)*x) )
       case (5,12);    forall(j = 2:nwn) wnbesj(j)=-2.*wn(j)*besselj(j,1)/(wn(j)*x)
       case (6,13);    forall(j = 2:nwn) wnbesj(j)=+2.*wn(j)*( besselj(j,0)-besselj(j,1)/(wn(j)*x) )
       case (7,14);    forall(j = 2:nwn) wnbesj(j)=+2.*wn(j)*besselj(j,1)/(wn(j)*x)
       case (17);      forall(j = 2:nwn) wnbesj(j)=-wn(j)*(wn(j)*re)**2*( besselj(j,0)-besselj(j,1)/(wn(j)*x) )
       case (18);      forall(j = 2:nwn) wnbesj(j)=-2*wn(j)*wn(j)*re*( -besselj(j,1)+besselj(j,2)/(wn(j)*x) )
       case (19);      forall(j = 2:nwn) wnbesj(j)=+2*wn(j)*wn(j)*re*besselj(j,2)/(wn(j)*x)
       case (20);      forall(j = 2:nwn) wnbesj(j)=+2*wn(j)*wn(j)*re*( -besselj(j,1)+besselj(j,2)/(wn(j)*x) )
       case (21);      forall(j = 2:nwn) wnbesj(j)=-2*wn(j)*wn(j)*re*besselj(j,2)/(wn(j)*x)
    end select
!
!  Neumann terms at wn3
!
    by0 = dbessy0(wn3*x)
    by1 = dbessy1(wn3*x)
    by2 = -by0+2./(wn3*x)*by1
    wn3besy = 0.d0
    select case (kint)
       case (1,8,22);  wn3besy=wn3*by0
       case (2,9,23);  wn3besy=-2*by1/re
       case (3,10,15); wn3besy=-wn3*wn3*by1*re
       case (4,11,16); wn3besy=-2.*wn3*(by0-by1/(wn3*x))
       case (5,12);    wn3besy=-2.*wn3*by1/(wn3*x)
       case (6,13);    wn3besy=+2.*wn3*(by0-by1/(wn3*x))
       case (7,14);    wn3besy=+2.*wn3*by1/(wn3*x)
       case (17);      wn3besy=-wn3*(wn3*re)**2*(by0-by1/(wn3*x))
       case (18);      wn3besy=-2*wn3*wn3*re*(-by1+by2/(wn3*x))
       case (19);      wn3besy=+2*wn3*wn3*re*by2/(wn3*x)
       case (20);      wn3besy=+2*wn3*wn3*re*(-by1+by2/(wn3*x))
       case (21);      wn3besy=-2*wn3*wn3*re*by2/(wn3*x)
    end select
!
!  perform integration for each frequency
!  Use trapezoidal rule if wn*x < 3 and Filon else
!
    allocate(greentap(nwn))
    call taperWavenumberIntegrals(nwn,dble(gfk%dwn),tapfrac,greentap)
    j33 = min(j3,nwn)
    zdis = 0.d0
    do j = 2,j33-1
       zdis = zdis+dcmplx(gwnspr(j),gwnspi(j))*wnbesj(j)*greentap(j)
    enddo
    zdis = zdis+0.5*dcmplx(gwnspr(j33),gwnspi(j33))*wnbesj(j33)*greentap(j33)
    zdis = zdis*gfk%dwn
    if (j33 .eq. nwn) goto 11
!
!  use the Filon integration for kx > 3
!  Neumann terms at the end of integration interval (k_N)
!
    wne=(nwn-1)*gfk%dwn
    by0 = dbessy0(wne*x)
    by1 = dbessy1(wne*x)
    by2 = -by0+2./(wne*x)*by1
    wnbesy = 0.d0
    select case (kint)
       case (1,8,22);  wnbesy=wne*by0
       case (2,9,23);  wnbesy=-2*by1/re
       case (3,10,15); wnbesy=-wne*wne*by1*re
       case (4,11,16); wnbesy=-2.*wne*(by0-by1/(wne*x))
       case (5,12);    wnbesy=-2.*wne*by1/(wne*x)
       case (6,13);    wnbesy=+2.*wne*(by0-by1/(wne*x))
       case (7,14);    wnbesy=+2.*wne*by1/(wne*x)
       case (17);      wnbesy=-wne*(wne*re)**2*(by0-by1/(wne*x))
       case (18);      wnbesy=-2*wne*wne*re*(-by1+by2/(wne*x))
       case (19);      wnbesy=+2*wne*wne*re*by2/(wne*x)
       case (20);      wnbesy=+2*wne*wne*re*(-by1+by2/(wne*x))
       case (21);      wnbesy=-2*wne*wne*re*by2/(wne*x)
    end select
!
    xdk = x*gfk%dwn
    zsum = dcmplx(gwnspr(j3),gwnspi(j3))*wnbesj(j3)*greentap(j3)
    do j = j3+1,nwn-1
       zsum = zsum+2.*dcmplx(gwnspr(j),gwnspi(j))*wnbesj(j)*greentap(j)
    enddo
    zsum = zsum+dcmplx(gwnspr(nwn),gwnspi(nwn))*wnbesj(nwn)*greentap(nwn)
    zdis = zdis+zsum/x*(1-cos(xdk))/xdk
    zdis = zdis+(dcmplx(gwnspr(nwn),gwnspi(nwn))*wnbesy*greentap(nwn) &
              & -dcmplx(gwnspr(j3),gwnspi(j3))*wn3besy*greentap(j3))/x*(1.-sin(xdk)/xdk)
!
!  multiply by R^2/(2*pi)
!
 11 zdis = zdis*re**2/(2.*mc_pid)
!
    deallocate(greentap,wnbesj)
    end subroutine forceWavenumberIntegrals
!-----------------------------------------------------------
!  Compute wavenumber integrals for moment tensor source
!
!  kint=1:       R^2/(2*pi)*int ( dk k U^1 J_0 )  
!  kint=2:       R^2/(2*pi)*int ( dk k U^2 J_0 )
!  kint=3:       R^2/(2*pi)*int ( dk k U^3 2/(kR) J_1 )
!  kint=4:       R^2/(2*pi)*int ( dk k U^4 2J_2 )
!  kint=5:       R^2/(2*pi)*int ( dk k V^1 -kR J_1 )
!  kint=6:       R^2/(2*pi)*int ( dk k V^2 -kR J_1 )
!  kint=7:       R^2/(2*pi)*int ( dk k V^3 2 (J_0-J_1/(kx)) )
!  kint=8:       R^2/(2*pi)*int ( dk k V^4 2kR (J_1-2J_2/(kx)) )
!  kint=9:       R^2/(2*pi)*int ( dk k W^1 2J_1/(kx) )
!  kint=10:      R^2/(2*pi)*int ( dk k W^2 4kRJ_2/(kx) )
!  kint=11:      R^2/(2*pi)*int ( dk k W^1 2(J_0-J_1/(kx)) )
!  kint=12:      R^2/(2*pi)*int ( dk k W^2 -2kR(J_1-2J_2/(kx)) )
!  kint=13:      R^2/(2*pi)*int ( dk k V^3 2J_1/(kx) )
!  kint=14:      R^2/(2*pi)*int ( dk k V^4 -4kRJ_2/(kx) )
!
!  kint=15:      R^2/(2*pi)*int ( dk k U^1_r J_0 )
!  kint=16:      R^2/(2*pi)*int ( dk k U^2_r J_0 )
!  kint=17:      R^2/(2*pi)*int ( dk k U^3_r 2/(kR) J_1 )
!  kint=18:      R^2/(2*pi)*int ( dk k U^4_r 2J_2 )
!  kint=19:      R^2/(2*pi)*int ( dk k V^1_r -kR J_1 )
!  kint=20:      R^2/(2*pi)*int ( dk k V^2_r -kR J_1 )
!  kint=21:      R^2/(2*pi)*int ( dk k V^3_r 2(J_0-J_1/(kx)) )
!  kint=22:      R^2/(2*pi)*int ( dk k V^4_r 2kR(J_1-2J_2/(kx)) )
!  kint=23:      R^2/(2*pi)*int ( dk k W^1_r 2J_1/(kx) )
!  kint=24:      R^2/(2*pi)*int ( dk k W^2_r 4kRJ_2/(kx) )
!  kint=25:      R^2/(2*pi)*int ( dk k W^1_r 2(J_0-J_1/(kx)) )
!  kint=26:      R^2/(2*pi)*int ( dk k W^2_r -2kR(J_1-2J_2/(kx)) )
!  kint=27:      R^2/(2*pi)*int ( dk k V^3_r 2J_1/(kx) )
!  kint=28:      R^2/(2*pi)*int ( dk k V^4_r -4kRJ_2/(kx) )
!
!  kint=29:      R^2/(2*pi)*int ( dk k U^1 -kR J_1 )  
!  kint=30:      R^2/(2*pi)*int ( dk k U^2 -kR J_1 )
!  kint=31:      R^2/(2*pi)*int ( dk k U^3 2 (J_0-J_1/(kx)) )
!  kint=32:      R^2/(2*pi)*int ( dk k U^4 2kR (J_1-2J_2/(kx)) )
!  kint=33:      R^2/(2*pi)*int ( dk k V^1 -(kR)^2 (J_0-J_1/(kx)) )
!  kint=34:      R^2/(2*pi)*int ( dk k V^2 -(kR)^2 (J_0-J_1/(kx)) )
!  kint=35:      R^2/(2*pi)*int ( dk k V^3 -2kR (J_1-J_2/(kx)) )
!  kint=36:      R^2/(2*pi)*int ( dk k V^4 2(kR)^2 ((6/(kx)^2 -1)J_2-J_1/(kx)) )
!  kint=37:      R^2/(2*pi)*int ( dk k W^1 -2kR J_2/(kx) )
!  kint=38:      R^2/(2*pi)*int ( dk k W^2 4(kR)^2 (-3J_2/(kx)^2 + J_1/(kx)) )
!  kint=39:      R^2/(2*pi)*int ( dk k W^1 -2kR (J_1-J_2/(kx)) )
!  kint=40:      R^2/(2*pi)*int ( dk k W^2 -2(kR)^2 ((6/(kx)^2 -1)J_2-J_1/(kx)) )
!  kint=41:      R^2/(2*pi)*int ( dk k V^3 -2kR J_2/(kx) )
!  kint=42:      R^2/(2*pi)*int ( dk k V^4 -4(kR)^2 (-3J_2/(kx)^2 + J_1/(kx)) )
!
!  kint=43:       R^2/(2*pi)*int ( dk k R^1 J_0 )  
!  kint=44:       R^2/(2*pi)*int ( dk k R^2 J_0 )
!  kint=45:       R^2/(2*pi)*int ( dk k R^3 2/(kR) J_1 )
!  kint=46:       R^2/(2*pi)*int ( dk k R^4 2J_2 )
!
!  for one receiver at x and one frequency using precomputed
!  values for the Bessel functions. This is possible because
!  the Bessel functions do not depend on frequency
!  and the wavenumber sampling is independent of frequency.
!  A complex value zdis is returned containing the
!  value of the desired integral.
!
!  Note: Unit of zdis is returned in meters/(Nm)
!
!  gfk:               greenFKSpectra object containing all meta information on Green FK spectra
!  gwnspr,gwnspi:     extracted wavenumber spectrum for selected frequency and wavenumber integral (double precision and imag)
!  kint:              index of requested wavenumber integral
!  xx:                epicentral distance in meters
!  besselj:           precomputed Bessel functions for all wavenumbers and desired distance
!  tapfrac:           fraction of wavenumber spectrum which is tapered
!  zdis:              complex value of wavenumber integral
!
    subroutine momentWavenumberIntegrals(gfk,gwnspr,gwnspi,kint,xx,besselj,tapfrac,zdis)
    type (green_fk_spectra) :: gfk
    real, dimension(:) :: gwnspr,gwnspi
    integer :: kint
    double precision :: xx,tapfrac
    double precision, dimension(:,0:) :: besselj
    double complex :: zdis
    double precision, dimension(:), allocatable :: wnbesj,greentap,wn
    double precision :: re,x,wn3,wne,wn3besy,xdk,wnbesy,by0,by1,by2
    integer :: j3,j33,j,nwn
    double complex :: zsum
!
!  number of wavenumbers and rearth
!
    nwn = size(gwnspr)
    re = gfk%rearth
!
!  guard against too small distances
!  if x is too small j3 may become larger than representable by integer*4
!  therefore choose x such that j3 < 1e9
!
    if(xx < 3.d0/(gfk%dwn*1.d9)) then
        x = 3.d0/(gfk%dwn*1.d9)
        print *,'<wnintMomentGreenFrequencyWavenumber>: WARNING: distance set to ',x,' instead less'
    else
        x = xx
    endif
    wn3=3.d0/x
!
!  find index of first wn-point greater than 3
!
    j3=int(wn3/gfk%dwn)+2
    wn3=(j3-1)*gfk%dwn
!
!  Bessel terms over whole interval
!  division by zero is unproblematic because we avoid wn=0
!  Bessel terms at k=0 vanish anyway
!  and x shouldn't be zero either
!
    allocate(wnbesj(nwn),wn(nwn))
    forall (j = 1:nwn) wn(j) = (j-1)*gfk%dwn
    select case (kint)
    case (1,15,43); forall(j = 1:nwn) wnbesj(j) = wn(j)*besselj(j,0)
    case (2,16,44); forall(j = 1:nwn) wnbesj(j) = wn(j)*besselj(j,0)
    case (3,17,45); forall(j = 1:nwn) wnbesj(j) = 2./re*besselj(j,1)
    case (4,18,46); forall(j = 1:nwn) wnbesj(j) = wn(j)*2.*besselj(j,2)
    case (5,19,29); forall(j = 1:nwn) wnbesj(j) = wn(j)*(-wn(j)*re)*besselj(j,1)
    case (6,20,30); forall(j = 1:nwn) wnbesj(j) = wn(j)*(-wn(j)*re)*besselj(j,1)
    case (7,21,31); forall(j = 1:nwn) wnbesj(j) = wn(j)*2.*(besselj(j,0)-besselj(j,1)/(wn(j)*x))
    case (8,22,32); forall(j = 1:nwn) wnbesj(j) = wn(j)*2.*wn(j)*re*(besselj(j,1)-2.*besselj(j,2)/(wn(j)*x))
    case (9,23);    forall(j = 1:nwn) wnbesj(j) = wn(j)*2.*besselj(j,1)/(wn(j)*x)
    case (10,24);   forall(j = 1:nwn) wnbesj(j) = wn(j)*4.*wn(j)*re*besselj(j,2)/(wn(j)*x)
    case (11,25);   forall(j = 1:nwn) wnbesj(j) = wn(j)*2.*(besselj(j,0)-besselj(j,1)/(wn(j)*x))
    case (12,26);   forall(j = 1:nwn) wnbesj(j) = -wn(j)*2.*wn(j)*re*(besselj(j,1)-2.*besselj(j,2)/(wn(j)*x))
    case (13,27);   forall(j = 1:nwn) wnbesj(j) = wn(j)*2.*besselj(j,1)/(wn(j)*x)
    case (14,28);   forall(j = 1:nwn) wnbesj(j) = -wn(j)*4.*wn(j)*re*besselj(j,2)/(wn(j)*x)
    case (33,34);   forall(j = 1:nwn) wnbesj(j) = -wn(j)*(wn(j)*re)**2*(besselj(j,0)-besselj(j,1)/(wn(j)*x))
    case (35,39);   forall(j = 1:nwn) wnbesj(j) = -wn(j)*2*wn(j)*re*(besselj(j,1)-besselj(j,2)/(wn(j)*x))
    case (36);      forall(j = 1:nwn) wnbesj(j) = +wn(j)*2*(wn(j)*re)**2*((6./(wn(j)*x)**2-1.)*besselj(j,2)-besselj(j,1)/(wn(j)*x))
    case (37,41);   forall(j = 1:nwn) wnbesj(j) = -wn(j)*2*wn(j)*re*besselj(j,2)/(wn(j)*x)
    case (38);      forall(j = 1:nwn) wnbesj(j) = +wn(j)*4.*(wn(j)*re)**2*(-3.*besselj(j,2)/(wn(j)*x)**2+besselj(j,1)/(wn(j)*x))
    case (40);      forall(j = 1:nwn) wnbesj(j) = -wn(j)*2*(wn(j)*re)**2*((6./(wn(j)*x)**2-1.)*besselj(j,2)-besselj(j,1)/(wn(j)*x))
    case (42);      forall(j = 1:nwn) wnbesj(j) = -wn(j)*4.*(wn(j)*re)**2*(-3.*besselj(j,2)/(wn(j)*x)**2+besselj(j,1)/(wn(j)*x))
    end select
!
!  Neumann terms at wn3
!
    by0 = dbessy0(wn3*x)
    by1 = dbessy1(wn3*x)
    by2 = -by0+2./(wn3*x)*by1
    wn3besy = 0.d0
    select case (kint)
    case (1,15,43); wn3besy = wn3*by0
    case (2,16,44); wn3besy = wn3*by0
    case (3,17,45); wn3besy = 2./re*by1
    case (4,18,46); wn3besy = wn3*2.*by2
    case (5,19,29); wn3besy = wn3*(-wn3*re)*by1
    case (6,20,30); wn3besy = wn3*(-wn3*re)*by1
    case (7,21,31); wn3besy = wn3*2.*(by0-by1/(wn3*x))
    case (8,22,32); wn3besy = wn3*2.*wn3*re*(by1-2.*by2/(wn3*x))
    case (9,23);    wn3besy = wn3*2.*by1/(wn3*x)
    case (10,24);   wn3besy = wn3*4.*wn3*re*by2/(wn3*x)
    case (11,25);   wn3besy = wn3*2.*(by0-by1/(wn3*x))
    case (12,26);   wn3besy = -wn3*2.*wn3*re*(by1-2.*by2/(wn3*x))
    case (13,27);   wn3besy = wn3*2.*by1/(wn3*x)
    case (14,28);   wn3besy = -wn3*4.*wn3*re*by2/(wn3*x)
    case (33,34);   wn3besy = -wn3*(wn3*re)**2*(by0-by1/(wn3*x))
    case (35,39);   wn3besy = -wn3*2*wn3*re*(by1-by2/(wn3*x))
    case (36);      wn3besy = +wn3*2*(wn3*re)**2*((6./(wn3*x)**2-1.)*by2-by1/(wn3*x))
    case (37,41);   wn3besy = -wn3*2*wn3*re*by2/(wn3*x)
    case (38);      wn3besy = +wn3*4.*(wn3*re)**2*(-3.*by2/(wn3*x)**2+by1/(wn3*x))
    case (40);      wn3besy = -wn3*2*(wn3*re)**2*((6./(wn3*x)**2-1.)*by2-by1/(wn3*x))
    case (42);      wn3besy = -wn3*4.*(wn3*re)**2*(-3.*by2/(wn3*x)**2+by1/(wn3*x))
    end select
!
!  perform integration for each frequency
!  Use trapezoidal rule if wn*x < 3 and Filon else
!
    allocate(greentap(nwn))
    call taperWavenumberIntegrals(nwn,dble(gfk%dwn),tapfrac,greentap)
    j33 = min(j3,nwn)
    zdis = 0.d0
    do j = 2,j33-1
       zdis = zdis+dcmplx(gwnspr(j),gwnspi(j))*wnbesj(j)*greentap(j)
    enddo
    zdis = zdis+0.5*dcmplx(gwnspr(j33),gwnspi(j33))*wnbesj(j33)*greentap(j33)
    zdis = zdis*gfk%dwn
    if (j33 .eq. nwn) goto 11
!
!  use the Filon integration for kx > 3
!  Neumann terms at the end of integration interval (k_N)
!
    wne=(nwn-1)*gfk%dwn
    by0=dbessy0(wne*x)
    by1=dbessy1(wne*x)
    by2=-by0+2./(wne*x)*by1
    wnbesy = 0.d0
    select case (kint)
    case (1,15,43); wnbesy = wne*by0
    case (2,16,44); wnbesy = wne*by0
    case (3,17,45); wnbesy = 2./re*by1
    case (4,18,46); wnbesy = wne*2.*by2
    case (5,19,29); wnbesy = wne*(-wne*re)*by1
    case (6,20,30); wnbesy = wne*(-wne*re)*by1
    case (7,21,31); wnbesy = wne*2.*(by0-by1/(wne*x))
    case (8,22,32); wnbesy = wne*2.*wne*re*(by1-2.*by2/(wne*x))
    case (9,23);    wnbesy = wne*2.*by1/(wne*x)
    case (10,24);   wnbesy = wne*4.*wne*re*by2/(wne*x)
    case (11,25);   wnbesy = wne*2.*(by0-by1/(wne*x))
    case (12,26);   wnbesy = -wne*2.*wne*re*(by1-2.*by2/(wne*x))
    case (13,27);   wnbesy = wne*2.*by1/(wne*x)
    case (14,28);   wnbesy = -wne*4.*wne*re*by2/(wne*x)
    case (33,34);   wnbesy = -wne*(wne*re)**2*(by0-by1/(wne*x))
    case (35,39);   wnbesy = -wne*2*wne*re*(by1-by2/(wne*x))
    case (36);      wnbesy = +wne*2*(wne*re)**2*((6./(wne*x)**2-1.)*by2-by1/(wne*x))
    case (37,41);   wnbesy = -wne*2*wne*re*by2/(wne*x)
    case (38);      wnbesy = +wne*4.*(wne*re)**2*(-3.*by2/(wne*x)**2+by1/(wne*x))
    case (40);      wnbesy = -wne*2*(wne*re)**2*((6./(wne*x)**2-1.)*by2-by1/(wne*x))
    case (42);      wnbesy = -wne*4.*(wne*re)**2*(-3.*by2/(wne*x)**2+by1/(wne*x))
    end select
!
    xdk = x*gfk%dwn
    zsum = dcmplx(gwnspr(j3),gwnspi(j3))*wnbesj(j3)*greentap(j3)
    do j = j3+1,nwn-1
       zsum = zsum+2.*dcmplx(gwnspr(j),gwnspi(j))*wnbesj(j)*greentap(j)
    enddo
    zsum = zsum+dcmplx(gwnspr(nwn),gwnspi(nwn))*wnbesj(nwn)*greentap(nwn)
    zdis = zdis+zsum/x*(1-cos(xdk))/xdk
    zdis = zdis+(dcmplx(gwnspr(nwn),gwnspi(nwn))*wnbesy*greentap(nwn) &
              & -dcmplx(gwnspr(j3),gwnspi(j3))*wn3besy*greentap(j3))/x*(1.-sin(xdk)/xdk)
!
!  multiply by R^2/(2*pi)
!
 11 zdis = zdis*re**2/(2.d0*mc_pid)
!
    deallocate(greentap,wnbesj)
    end subroutine momentWavenumberIntegrals
!-------------------------------------------------------------
!  Compute a cos**2 taper for given frequency
!  internal use
!
    subroutine taperWavenumberIntegrals(nwn,dwn,tapfrac,greentap)
    integer :: nwn
    double precision :: dwn,tapfrac
    double precision, dimension(:) :: greentap
    double precision :: c1,c2
    integer :: jt,j
!
    jt = min(nint(nwn*(1.-tapfrac)),nwn)
    jt = max(1,jt)
    c1 = (nwn-1)*dwn*tapfrac*2./mc_pi
    c2 = (1.-tapfrac)/tapfrac*0.5*mc_pi
    greentap(1:jt-1) = 1.
    forall (j = jt:nwn) greentap(j)=cos( (j-1)*dwn/c1-c2 )**2
    end subroutine taperWavenumberIntegrals
!-------------------------------------------------------------------------
!  get wavenumber intergral index array
!  depending on dsvmask and source type    
!  nwint:    total number of required wavenumber integrals (out)
!  kint:     index array of length nwint with integral indices (out)
!
    subroutine getIndexArrayWavenumberIntegrals(gfk,nwint,kint)
    type (green_fk_spectra) :: gfk
    integer :: nwint,i
    integer, dimension(:), allocatable :: kint
!
!  Moment tensor source
!
    if (gfk%istyp == 1) then
        if (gfk%dsvmask(2) == 1 .and. gfk%dsvmask(5) == 1) then    ! stresses and derivatives
            nwint = 46
            allocate(kint(nwint))
            kint = (/ (i,i=1,46) /)
        else if (gfk%dsvmask(2) == 1 .and. gfk%dsvmask(5) == 0) then   ! stresses without derivatives
            nwint = 18
            allocate(kint(nwint))
            kint = (/ (i,i=1,14),(i,i=43,46) /)
        else if (gfk%dsvmask(2) == 0 .and. gfk%dsvmask(5) == 1) then   ! derivatives without stresses
            nwint = 42
            allocate(kint(nwint))
            kint = (/ (i,i=1,42) /)
        else if (gfk%dsvmask(2) == 0 .and. gfk%dsvmask(5) == 0) then     ! neither stresses nor derivatives
            nwint = 14
            allocate(kint(nwint))
            kint = (/ (i,i=1,14) /)
        endif
!
!  Force source
!
    else if (gfk%istyp == 0) then
        if (gfk%dsvmask(2) == 1 .and. gfk%dsvmask(5) == 1) then    ! stresses and derivatives
            nwint = 23
            allocate(kint(nwint))
            kint = (/ (i,i=1,23) /)
        else if (gfk%dsvmask(2) == 1 .and. gfk%dsvmask(5) == 0) then   ! stresses without derivatives
            nwint = 9
            allocate(kint(nwint))
            kint = (/ (i,i=1,7),(i,i=22,23) /)
        else if (gfk%dsvmask(2) == 0 .and. gfk%dsvmask(5) == 1) then   ! derivatives without stresses
            nwint = 21
            allocate(kint(nwint))
            kint = (/ (i,i=1,21) /)
        else if (gfk%dsvmask(2) == 0 .and. gfk%dsvmask(5) == 0) then     ! neither stresses nor derivatives
            nwint = 7
            allocate(kint(nwint))
            kint = (/ (i,i=1,7) /)
        endif
    endif
    end subroutine getIndexArrayWavenumberIntegrals
!
end module wavenumberIntegrals

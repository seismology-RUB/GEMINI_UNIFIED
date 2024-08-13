! ====================================================================================
!  Calculate various sums over harmonci degree using Green Om-El spectra
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
module harmonicDegreeSums
   use greenFKSpectra
   use legendreFunctions
   use mathConstants
   implicit none
!                                 FORCE EXCITATION
!
!  Mapping from the 23 harmonic degree sums to the two Green spectra indices:
!
   integer, dimension(23,2) :: force_component_jump_from_harmonic_degree_sum = reshape( &
        & (/ 1,1,3,3,7,7,3, 5,5,6,6,9,9,6, 1,1,3,3,7,7,3, 2,2, &
        &    1,2,1,2,1,1,2, 1,2,1,2,1,1,2, 1,2,1,2,1,1,2, 1,2 /) &
        & ,(/ 23,2 /))
!                                 MOMENT EXCITATION
!
!  Mapping from the 46 harmonic degree sums to the two Green spectra indices:
!
    integer, dimension(46,2) :: moment_component_jump_from_harmonic_degree_sum = reshape( &
         & (/ 1,1,1,1,3,3,3,3,7,7,7,7,3,3, 5,5,5,5,6,6,6,6,9,9,9,9,6,6, 1,1,1,1,3,3,3,3,7,7,7,7,3,3, 2,2,2,2,  &
         &    1,2,3,4,1,2,3,4,1,2,1,2,3,4, 1,2,3,4,1,2,3,4,1,2,1,2,3,4, 1,2,3,4,1,2,3,4,1,2,1,2,3,4, 1,2,3,4 /) &
         & ,(/ 46,2 /))
!
!  Maximum wavenumber index depending on source type
!
    integer, dimension(0:1) :: maxIndexHarmonicDegreeSums = (/23,46/)
!    
contains
!-----------------------------------------------------------
!  Compute harmonic degree sums for single force source
!
!  for one receiver and one frequency using precomputed
!  values for the Legendre functions. This is possible because
!  the Legendre functions do not depend on frequency
!  and the harmonic degree spacing is independent of frequency.
!  A complex value zdis is returned containing the
!  value of the desired sum.
!
!  Note: Unit of zdis is returned in meters/N
!
!  gfk:               greenFKSpectra object containing all meta information on Green FK spectra
!  gwnspr,gwnspi:     extracted wavenumber spectrum for selected frequency and wavenumber integral (real and imag)
!  kint:              index of requested harmonic degree sum
!  alf:               legendre functions object with precomputed Legendre functions for all l,m and desired theta
!  tapfrac:           fraction of wavenumber spectrum which is tapered
!  zdis:              complex value of harmonic degree sum
!
!  kint= 1:       sum_l (gamma_l^2 P_l0 U_l^1)                Fr
!  kint= 2:       sum_l (gamma_l^2 P_l1 U_l^2 (-2/elp1))      (Ft*cf + Ff*sf)
!  kint= 3:       sum_l (gamma_l^2 dpdt_l0 V_l^1)             Fr
!  kint= 4:       sum_l (gamma_l^2 dpdt_l1 V_l^2 (-2/elp1))   (Ft*cf + Ff*sf)
!  kint= 5:       sum_l (gamma_l^2 pdst_l1 W_l^1 (-2/elp1))   (Ft*cf + Ff*sf)
!  kint= 6:       sum_l (gamma_l^2 dpdt_l1 W_l^1 (+2/elp1))   (Ft*sf - Ff*cf)
!  kint= 7:       sum_l (gamma_l^2 pdst_l1 V_l^2 (+2/elp1))   (Ft*sf - Ff*cf)
!
!  kint= 8:       sum_l (gamma_l^2 P_l0 UR_l^1)                Fr 
!  kint= 9:       sum_l (gamma_l^2 P_l1 UR_l^2 (-2/elp1))      (Ft*cf + Ff*sf)
!  kint=10:       sum_l (gamma_l^2 dpdt_l0 VR_l^1)             Fr
!  kint=11:       sum_l (gamma_l^2 dpdt_l1 VR_l^2 (-2/elp1))   (Ft*cf + Ff*sf)
!  kint=12:       sum_l (gamma_l^2 pdst_l1 WR_l^1 (-2/elp1))   (Ft*cf + Ff*sf)
!  kint=13:       sum_l (gamma_l^2 dpdt_l1 WR_l^1 (+2/elp1))   (Ft*sf - Ff*cf)
!  kint=14:       sum_l (gamma_l^2 pdst_l1 VR_l^2 (+2/elp1))   (Ft*sf - Ff*cf)
!
!  kint=15:       sum_l (gamma_l^2 dPdt_l0 U_l^1)               Fr 
!  kint=16:       sum_l (gamma_l^2 dPdt_l1 U_l^2 (-2/elp1))     (Ft*cf + Ff*sf)
!  kint=17:       sum_l (gamma_l^2 d2pdt2_l0 V_l^1)             Fr
!  kint=18:       sum_l (gamma_l^2 d2pdt2_l1 V_l^2 (-2/elp1))   (Ft*cf + Ff*sf)
!  kint=19:       sum_l (gamma_l^2 dtpdst_l1 W_l^1 (-2/elp1))   (Ft*cf + Ff*sf)
!  kint=20:       sum_l (gamma_l^2 d2pdt2_l1 W_l^1 (+2/elp1))   (Ft*sf - Ff*cf)
!  kint=21:       sum_l (gamma_l^2 dtpdst_l1 V_l^2 (+2/elp1))   (Ft*sf - Ff*cf)
!
    subroutine forceHarmonicDegreeSums(gwnspr,gwnspi,kint,alf,tapfrac,zdis)
    real, dimension(:) :: gwnspr,gwnspi
    integer :: kint
    double precision :: tapfrac
    double precision :: fourpi
    type (legendre_functions) :: alf
    double complex :: zdis
    double precision, dimension(:), allocatable :: greentap
    integer :: l,lmax
!
!  number of harmonic degrees (lmax)
!
    lmax = size(gwnspr)-1
    fourpi = 4.d0*mc_pid
!
!  perform summation for each frequency
!
    allocate(greentap(0:lmax))
    call taperHarmonicDegreeSums(lmax,tapfrac,greentap)
    zdis = 0.d0
    select case (kint)
    case (1,8)                  ! Fr, ur, d(ur)/dr
       do l = 0,lmax
          zdis = zdis+alf%plm(l,0)*dcmplx(gwnspr(l+1),gwnspi(l+1))*(2.d0*l+1.d0)/fourpi*greentap(l)
       enddo
    case (2,9)                  ! (Ft*cf+Ff*sf), ur, d(ur)/dr
       do l = 1,lmax
          zdis = zdis+alf%plm(l,1)*(-2.d0/float(l*(l+1)))*dcmplx(gwnspr(l+1),gwnspi(l+1))*(2.d0*l+1.d0)/fourpi*greentap(l)
       enddo
    case (3,10,15)              ! Fr, ut, d(ut)/dr, d(ur)/dt
       do l = 1,lmax
          zdis = zdis+alf%dpdt(l,0)*dcmplx(gwnspr(l+1),gwnspi(l+1))*(2.d0*l+1.d0)/fourpi*greentap(l)
       enddo
    case (4,11,16)              ! (Ft*cf+Ff*sf), ut, d(ut)/dr, d(ur)/dt 
       do l = 1,lmax
          zdis = zdis+alf%dpdt(l,1)*(-2.d0/float(l*(l+1)))*dcmplx(gwnspr(l+1),gwnspi(l+1))*(2.d0*l+1.d0)/fourpi*greentap(l)
       enddo
    case (5,12)                 ! (Ft*cf+Ff*sf), ut, d(ut)/dr
       do l = 1,lmax
          zdis = zdis+alf%pdst(l,1)*(-2.d0/float(l*(l+1)))*dcmplx(gwnspr(l+1),gwnspi(l+1))*(2.d0*l+1.d0)/fourpi*greentap(l)
       enddo
    case (6,13)                 ! (Ft*sf-Ff*cf), uf, d(uf)/dr
       do l = 1,lmax
          zdis = zdis+alf%dpdt(l,1)*(+2.d0/float(l*(l+1)))*dcmplx(gwnspr(l+1),gwnspi(l+1))*(2.d0*l+1.d0)/fourpi*greentap(l)
       enddo
    case (7,14)                 ! (Ft*sf-Ff*cf), uf, d(uf)/dr
       do l = 1,lmax
          zdis = zdis+alf%pdst(l,1)*(+2.d0/float(l*(l+1)))*dcmplx(gwnspr(l+1),gwnspi(l+1))*(2.d0*l+1.d0)/fourpi*greentap(l)
       enddo
    case (17)                   ! Fr, d(ut)/dt
       do l = 1,lmax
          zdis = zdis+alf%d2pdt2(l,0)*dcmplx(gwnspr(l+1),gwnspi(l+1))*(2.d0*l+1.d0)/fourpi*greentap(l)
       enddo
    case (18)                   ! (Ft*cf+Ff*sf), d(ut)/dt 
       do l = 1,lmax
          zdis = zdis+alf%d2pdt2(l,1)*(-2.d0/float(l*(l+1)))*dcmplx(gwnspr(l+1),gwnspi(l+1))*(2.d0*l+1.d0)/fourpi*greentap(l)
       enddo
    case (19)                   ! (Ft*cf+Ff*sf), d(ut)/dt
       do l = 1,lmax
          zdis = zdis+alf%dtpdst(l,1)*(-2.d0/float(l*(l+1)))*dcmplx(gwnspr(l+1),gwnspi(l+1))*(2.d0*l+1.d0)/fourpi*greentap(l)
       enddo
    case (20)                   ! (Ft*sf-Ff*cf), d(uf)/dt
       do l = 1,lmax
          zdis = zdis+alf%d2pdt2(l,1)*(+2.d0/float(l*(l+1)))*dcmplx(gwnspr(l+1),gwnspi(l+1))*(2.d0*l+1.d0)/fourpi*greentap(l)
       enddo
    case (21)                   ! (Ft*sf-Ff*cf), d(uf)/dt
       do l = 1,lmax
          zdis = zdis+alf%dtpdst(l,1)*(+2.d0/float(l*(l+1)))*dcmplx(gwnspr(l+1),gwnspi(l+1))*(2.d0*l+1.d0)/fourpi*greentap(l)
       enddo
    end select
!
    deallocate(greentap)
    end subroutine forceHarmonicDegreeSums
!-----------------------------------------------------------
!  Compute harmonic degree sums for moment tensor source
!
!  for one receiver and one frequency using precomputed
!  values for the Legendre functions. This is possible because
!  the Legendre functions do not depend on frequency
!  and the harmonic degree sampling is independent of frequency.
!  A complex value zdis is returned containing the
!  value of the desired sum.
!
!  Note: Unit of zdis is returned in meters/(Nm)
!
!  gwnspr,gwnspi:     extracted wavenumber spectrum for selected frequency and wavenumber integral (real and imag)
!  kint:              index of requested wavenumber integral
!  alf:               precomputed Legendre functions for all l,m and desired theta
!  tapfrac:           fraction of wavenumber spectrum which is tapered
!  zdis:              complex value of wavenumber integral
!
!  kint = 1:  sum_l (gamma_l^2 P(l,0) U_l^1)
!  kint = 2:  sum_l (gamma_l^2 P(l,0) U_l^2)
!  kint = 3:  sum_l (gamma_l^2 P(l,1) 2/elp1 U_l^3)
!  kint = 4:  sum_l (gamma_l^2 P(l,2) 2/elp1 U_l^4)
!  kint = 5:  sum_l (gamma_l^2 dpdt(l,0) V_l^1    
!  kint = 6:  sum_l (gamma_l^2 dpdt(l,0) V_l^2
!  kint = 7:  sum_l (gamma_l^2 dpdt(l,1) 2/elp1 V_l^3
!  kint = 8:  sum_l (gamma_l^2 dpdt(l,2) 2/elp1 V_l^4
!  kint = 9:  sum_l (gamma_l^2 pdst(l,1) 2/elp1 W_l^1
!  kint = 10: sum_l (gamma_l^2 pdst(l,2) 4/elp1 W_l^2
!  kint = 11: sum_l (gamma_l^2 dpdt(l,1) 2/elp1 W_l^1 
!  kint = 12: sum_l (gamma_l^2 dpdt(l,2) -2/elp1 W_l^2
!  kint = 13: sum_l (gamma_l^2 pdst(l,1) 2/elp1 V_l^3
!  kint = 14: sum_l (gamma_l^2 pdst(l,2) -4/elp1 V_l^4
!
!  kint = 15: sum_l (gamma_l^2 P(l,0) UR_l^1)
!  kint = 16: sum_l (gamma_l^2 P(l,0) UR_l^2)
!  kint = 17: sum_l (gamma_l^2 P(l,1) 2/elp1 UR_l^3)
!  kint = 18: sum_l (gamma_l^2 P(l,2) 2/elp1 UR_l^4)
!  kint = 19: sum_l (gamma_l^2 dpdt(l,0) VR_l^1    
!  kint = 20: sum_l (gamma_l^2 dpdt(l,0) VR_l^2
!  kint = 21: sum_l (gamma_l^2 dpdt(l,1) 2/elp1 VR_l^3
!  kint = 22: sum_l (gamma_l^2 dpdt(l,2) 2/elp1 VR_l^4
!  kint = 23: sum_l (gamma_l^2 pdst(l,1) 2/elp1 WR_l^1
!  kint = 24: um_l (gamma_l^2 pdst(l,2) 4/elp1 WR_l^2
!  kint = 25: sum_l (gamma_l^2 dpdt(l,1) 2/elp1 WR_l^1 
!  kint = 26: sum_l (gamma_l^2 dpdt(l,2) -2/elp1 WR_l^2
!  kint = 27: sum_l (gamma_l^2 pdst(l,1) 2/elp1 VR_l^3
!  kint = 28: sum_l (gamma_l^2 pdst(l,2) -4/elp1 VR_l^4
!
!  kint = 29:  sum_l (gamma_l^2 dPdt(l,0) U_l^1)
!  kint = 30:  sum_l (gamma_l^2 dPdt(l,0) U_l^2)
!  kint = 31:  sum_l (gamma_l^2 dPdt(l,1) 2/elp1 U_l^3)
!  kint = 32:  sum_l (gamma_l^2 dPdt(l,2) 2/elp1 U_l^4)
!  kint = 33:  sum_l (gamma_l^2 d2pdt2(l,0) V_l^1    
!  kint = 34:  sum_l (gamma_l^2 d2pdt2(l,0) V_l^2
!  kint = 35:  sum_l (gamma_l^2 d2pdt2(l,1) 2/elp1 V_l^3
!  kint = 36:  sum_l (gamma_l^2 d2pdt2(l,2) 2/elp1 V_l^4
!  kint = 37:  sum_l (gamma_l^2 dtpdst(l,1) 2/elp1 W_l^1
!  kint = 38: sum_l (gamma_l^2 dtpdst(l,2) 4/elp1 W_l^2
!  kint = 39: sum_l (gamma_l^2 d2pdt2(l,1) 2/elp1 W_l^1 
!  kint = 40: sum_l (gamma_l^2 d2pdt2(l,2) -2/elp1 W_l^2
!  kint = 41: sum_l (gamma_l^2 dtpdst(l,1) 2/elp1 V_l^3
!  kint = 42: sum_l (gamma_l^2 dtpdst(l,2) -4/elp1 V_l^4
!
    subroutine momentHarmonicDegreeSums(gwnspr,gwnspi,kint,alf,tapfrac,zdis)
    real, dimension(:) :: gwnspr,gwnspi
    integer :: kint
    double precision :: tapfrac
    double precision :: fourpi
    type (legendre_functions) :: alf
    double complex :: zdis
    double precision, dimension(:), allocatable :: greentap
    integer :: l,lmax
!
!  number of wavenumbers and rearth
!
    fourpi = 4.d0*mc_pid
    lmax = size(gwnspr)-1
!
!  perform summation over Legendre functions
!
    allocate(greentap(0:lmax))
    call taperHarmonicDegreeSums(lmax,tapfrac,greentap)
    zdis = 0.d0
    select case (kint)
    case (1,15)             ! Mrr, ur, d(ur)/dr
       do l = 0,lmax
          zdis = zdis+alf%plm(l,0)*dcmplx(gwnspr(l+1),gwnspi(l+1))*(2.d0*l+1.d0)/fourpi*greentap(l)
       enddo
    case (2,16)             ! Mtt+Mff, ur, d(ur)/dr
       do l = 0,lmax
          zdis = zdis+alf%plm(l,0)*dcmplx(gwnspr(l+1),gwnspi(l+1))*(2.d0*l+1.d0)/fourpi*greentap(l)
       enddo
    case (3,17)             ! Mrt*cf+Mrf*sf, ur, d(ur)/dr
       do l = 1,lmax
          zdis = zdis+alf%plm(l,1)*2.d0/float(l*(l+1))*dcmplx(gwnspr(l+1),gwnspi(l+1))*(2.d0*l+1.d0)/fourpi*greentap(l)
       enddo
    case (4,18)             ! (Mff-Mtt)*c2f-2Mtf*sin2f, ur, d(ur)/dr
       do l = 2,lmax
          zdis = zdis+alf%plm(l,2)*2.d0/float(l*(l+1))*dcmplx(gwnspr(l+1),gwnspi(l+1))*(2.d0*l+1.d0)/fourpi*greentap(l) 
       enddo
    case (5,19,29)          ! Mrr, ut(V), d(ut)/dr, d(ur)/dt
       do l = 1,lmax
          zdis = zdis+alf%dpdt(l,0)*dcmplx(gwnspr(l+1),gwnspi(l+1))*(2.d0*l+1.d0)/fourpi*greentap(l)
       enddo
    case (6,20,30)          ! Mtt+Mff, ut(V), d(ut)/dr, d(ur)/dt
       do l = 1,lmax
          zdis = zdis+alf%dpdt(l,0)*dcmplx(gwnspr(l+1),gwnspi(l+1))*(2.d0*l+1.d0)/fourpi*greentap(l)
       enddo
    case (7,21,31)          ! Mrt*cf+Mrf*sf, ut(V), d(ut)/dr, d(ur)/dt
       do l = 1,lmax
          zdis = zdis+alf%dpdt(l,1)*2.d0/float(l*(l+1))*dcmplx(gwnspr(l+1),gwnspi(l+1))*(2.d0*l+1.d0)/fourpi*greentap(l)
       enddo
    case (8,22,32)          ! (Mff-Mtt)*c2f-2Mtf*sin2f, ut(V), d(ut)/dr, d(ur)/dt
       do l = 2,lmax
          zdis = zdis+alf%dpdt(l,2)*2.d0/float(l*(l+1))*dcmplx(gwnspr(l+1),gwnspi(l+1))*(2.d0*l+1.d0)/fourpi*greentap(l) 
       enddo
    case (9,23)             ! Mrf*sf+Mrt*cf, ut(W), d(ut)/dr
       do l = 1,lmax
          zdis = zdis+alf%pdst(l,1)*2.d0/float(l*(l+1))*dcmplx(gwnspr(l+1),gwnspi(l+1))*(2.d0*l+1.d0)/fourpi*greentap(l)
       enddo
    case (10,24)            ! (Mff-Mtt)*c2f-2Mtf*sin2f, ut(W), d(ut)/dr
       do l = 2,lmax
          zdis = zdis+alf%pdst(l,2)*4.d0/float(l*(l+1))*dcmplx(gwnspr(l+1),gwnspi(l+1))*(2.d0*l+1.d0)/fourpi*greentap(l) 
       enddo
    case (11,25)             ! Mrf*cf-Mrt*sf, uf(W), d(uf)/dr
       do l = 1,lmax
          zdis = zdis+alf%dpdt(l,1)*2.d0/float(l*(l+1))*dcmplx(gwnspr(l+1),gwnspi(l+1))*(2.d0*l+1.d0)/fourpi*greentap(l)
       enddo
    case (12,26)             ! (Mff-Mtt)*s2f+2Mtf*c2f, uf(W), d(uf)/dr
       do l = 2,lmax
          zdis = zdis+alf%dpdt(l,2)*(-2.d0)/float(l*(l+1))*dcmplx(gwnspr(l+1),gwnspi(l+1))*(2.d0*l+1.d0)/fourpi*greentap(l) 
       enddo
    case (13,27)             ! Mrf*cf-Mrt*sf, uf(V), d(ut)/dr
       do l = 1,lmax
          zdis = zdis+alf%pdst(l,1)*2.d0/float(l*(l+1))*dcmplx(gwnspr(l+1),gwnspi(l+1))*(2.d0*l+1.d0)/fourpi*greentap(l)
       enddo
    case (14,28)             ! (Mff-Mtt)*s2f+2Mtf*c2f, uf(V), d(ut)/dr
       do l = 2,lmax
          zdis = zdis+alf%pdst(l,2)*(-4.d0)/float(l*(l+1))*dcmplx(gwnspr(l+1),gwnspi(l+1))*(2.d0*l+1.d0)/fourpi*greentap(l) 
       enddo
    case (33)                ! Mrr, d(ut)/dt
       do l = 1,lmax
          zdis = zdis+alf%d2pdt2(l,0)*dcmplx(gwnspr(l+1),gwnspi(l+1))*(2.d0*l+1.d0)/fourpi*greentap(l)
       enddo
    case (34)                ! Mtt+Mff, d(ut)/dt
       do l = 1,lmax
          zdis = zdis+alf%d2pdt2(l,0)*dcmplx(gwnspr(l+1),gwnspi(l+1))*(2.d0*l+1.d0)/fourpi*greentap(l)
       enddo
    case (35)                ! Mrt*cf+Mrf*sf, d(ut)/dt
       do l = 1,lmax
          zdis = zdis+alf%d2pdt2(l,1)*2.d0/float(l*(l+1))*dcmplx(gwnspr(l+1),gwnspi(l+1))*(2.d0*l+1.d0)/fourpi*greentap(l)
       enddo
    case (36)                ! (Mff-Mtt)*c2f-2Mtf*sin2f, d(ut)/dt
       do l = 2,lmax
          zdis = zdis+alf%d2pdt2(l,2)*2.d0/float(l*(l+1))*dcmplx(gwnspr(l+1),gwnspi(l+1))*(2.d0*l+1.d0)/fourpi*greentap(l) 
       enddo
    case (37)                ! Mrf*sf+Mrt*cf, d(ut)/dt
       do l = 1,lmax
          zdis = zdis+alf%dtpdst(l,1)*2.d0/float(l*(l+1))*dcmplx(gwnspr(l+1),gwnspi(l+1))*(2.d0*l+1.d0)/fourpi*greentap(l)
       enddo
    case (38)                ! (Mff-Mtt)*c2f-2Mtf*sin2f, d(ut)/dt
       do l = 2,lmax
          zdis = zdis+alf%dtpdst(l,2)*4.d0/float(l*(l+1))*dcmplx(gwnspr(l+1),gwnspi(l+1))*(2.d0*l+1.d0)/fourpi*greentap(l) 
       enddo
    case (39)                ! Mrf*cf-Mrt*sf, d(uf)/dt
       do l = 1,lmax
          zdis = zdis+alf%d2pdt2(l,1)*2.d0/float(l*(l+1))*dcmplx(gwnspr(l+1),gwnspi(l+1))*(2.d0*l+1.d0)/fourpi*greentap(l)
       enddo
    case (40)                ! (Mff-Mtt)*s2f+2Mtf*c2f, d(uf)/dt
       do l = 2,lmax
          zdis = zdis+alf%d2pdt2(l,2)*(-2.d0)/float(l*(l+1))*dcmplx(gwnspr(l+1),gwnspi(l+1))*(2.d0*l+1.d0)/fourpi*greentap(l) 
       enddo
    case (41)                ! Mrf*cf-Mrt*sf, d(uf)/dt
       do l = 1,lmax
          zdis = zdis+alf%dtpdst(l,1)*2.d0/float(l*(l+1))*dcmplx(gwnspr(l+1),gwnspi(l+1))*(2.d0*l+1.d0)/fourpi*greentap(l)
       enddo
    case (42)                ! (Mff-Mtt)*s2f+2Mtf*c2f, d(uf)/dt
       do l = 2,lmax
          zdis = zdis+alf%dtpdst(l,2)*(-4.d0)/float(l*(l+1))*dcmplx(gwnspr(l+1),gwnspi(l+1))*(2.d0*l+1.d0)/fourpi*greentap(l) 
       enddo
    end select
!
    deallocate(greentap)
    end subroutine momentHarmonicDegreeSums
!-------------------------------------------------------------
!  Compute a cos**2 taper for given frequency
!  internal use
!
    subroutine taperHarmonicDegreeSums(lmax,tapfrac,greentap)
    integer :: lmax
    double precision :: tapfrac
    double precision, dimension(0:) :: greentap
    double precision :: c1,c2
    integer :: jt,j
!
    jt = min(nint(lmax*(1.d0-tapfrac)),lmax)
    jt = max(1,jt)
    c1 = lmax*tapfrac*2.d0/mc_pid
    c2 = (1.-tapfrac)/tapfrac*0.5*mc_pid
    greentap(0:jt-1) = 1.d0
    forall (j = jt:lmax) greentap(j)=cos( float(j)/c1-c2 )**2
    end subroutine taperHarmonicDegreeSums
!-------------------------------------------------------------------------
!  get harmonic degree sum index array
!  depending on dsvmask and source type    
!  nwint:    total number of harmonic degree sums (out)
!  kint:     index array of length nwint with harmonic degree sum indices (out)
!
    subroutine getIndexArrayHarmonicDegreeSums(gfk,nwint,kint)
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
    end subroutine getIndexArrayHarmonicDegreeSums
!
end module harmonicDegreeSums

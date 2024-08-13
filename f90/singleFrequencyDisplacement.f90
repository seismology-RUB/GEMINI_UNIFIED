! ====================================================================================
!  Compute single frequency displacement using wavenumber integrals
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
!  Module with routines to calculate single frequency displacements
!  and displacement derivatives using precalculated wavenumber integrals
!  for given force or moment tensor excitation
!----------------------------------------------------------------------
 module singleFrequencyDisplacement
    implicit none
!
 contains
!-------------------------------------------------------------------
!  single frequency displacement for single force and ZRT components
!  from precomputed wavenumber integrals
!  phi:  receiver azimuth in rad
!
    subroutine forceSingleFrequencyDisplacement(phi,force,wnint,zsp)
    double precision :: phi
    double precision, dimension(:) :: force       !  force components
    double complex, dimension(:) :: wnint         !  wavenumber integrals
    double complex, dimension(:) :: zsp           !  displacement
    double precision :: sphi,cphi
!
    sphi = sin(phi); cphi = cos(phi)
    zsp(1) = wnint(1)*force(1)+wnint(2)*(force(2)*cphi+force(3)*sphi)
    zsp(2) = wnint(3)*force(1)+(wnint(4)+wnint(5))*(force(2)*cphi+force(3)*sphi)
    zsp(3) = (wnint(6)+wnint(7))*(force(2)*sphi-force(3)*cphi)
    end subroutine forceSingleFrequencyDisplacement
!---------------------------------------------------------------
!  radial derivative of displacement spectrum for single force and one component
!  comp: component (Z,R,T)
!  phi:  receiver azimuth in rad
!
    subroutine forceRadDerivSingleFrequencyDisplacement(phi,force,wnint,zsp)
    double precision :: phi
    double precision, dimension(:) :: force       !  force components
    double complex, dimension(:) :: wnint         !  wavenumber integrals
    double complex, dimension(:) :: zsp           !  displacement
    double precision :: sphi,cphi
!
    sphi = sin(phi); cphi = cos(phi)
    zsp(1) = wnint(8)*force(1)+wnint(9)*(force(2)*cphi+force(3)*sphi)
    zsp(2) = wnint(10)*force(1)+(wnint(11)+wnint(12))*(force(2)*cphi+force(3)*sphi)
    zsp(3) = (wnint(13)+wnint(14))*(force(2)*sphi-force(3)*cphi)
    end subroutine forceRadDerivSingleFrequencyDisplacement
!-------------------------------------------------------------------------
!  1/r*theta derivative of displacement for single force and one component
!  comp: component (Z,R,T)
!  r:    receiver radius in m
!  phi:  receiver azimuth in rad
!
    subroutine forceThetaDerivSingleFrequencyDisplacement(r,phi,force,wnint,zsp)
    double precision :: r,phi
    double precision, dimension(:) :: force       !  force components
    double complex, dimension(:) :: wnint         !  wavenumber integrals vs frequency
    double complex, dimension(:) :: zsp           !  displacement spectrum
    double precision :: sphi,cphi
!
    sphi = sin(phi); cphi = cos(phi)
    zsp(1) = wnint(15)*force(1)+wnint(16)*(force(2)*cphi+force(3)*sphi)
    zsp(2) = wnint(17)*force(1)+(wnint(18)+wnint(19))*(force(2)*cphi+force(3)*sphi)
    zsp(3) = (wnint(20)+wnint(21))*(force(2)*sphi-force(3)*cphi)
    zsp = zsp/r
    end subroutine forceThetaDerivSingleFrequencyDisplacement
!---------------------------------------------------------------
!  1/(r*sin(theta))*phi-derivative of displacement for single force and one component
!  comp: component (Z,R,T)
!  r:      receiver radius in m
!  delta:  epicentral distance in rad
!  phi:    receiver azimuth in rad
!
    subroutine forcePhiDerivSingleFrequencyDisplacement(r,delta,phi,force,wnint,zsp)
    double precision :: r,delta,phi
    double precision, dimension(:) :: force            !  force components
    double complex, dimension(:) :: wnint              !  wavenumber integrals vs frequency
    double complex, dimension(:) :: zsp                !  displacement spectrum
    double precision :: sphi,cphi,rst
!
    sphi = sin(phi); cphi = cos(phi)
    rst = r*sin(delta)
    zsp(1) = wnint(2)*(-force(2)*sphi+force(3)*cphi)
    zsp(2) = (wnint(4)+wnint(5))*(-force(2)*sphi+force(3)*cphi)
    zsp(3) = (wnint(6)+wnint(7))*(+force(2)*cphi+force(3)*sphi)
    zsp = zsp/rst
    end subroutine forcePhiDerivSingleFrequencyDisplacement
!--------------------------------------------------------------------
!  single frequency displacement for moment tensor
!  phi:  receiver azimuth in rad
!
    subroutine momentSingleFrequencyDisplacement(phi,mt,wnint,zsp)
    double precision :: phi
    double precision, dimension(:) :: mt               !  moment tensor components
    double complex, dimension(:) :: zsp                !  displacement spectrum
    double complex, dimension(:) :: wnint              !  wavenumber integrals vs frequency
    double precision :: sphi,cphi,c2phi,s2phi
!
    sphi = sin(phi); cphi = cos(phi)
    s2phi = sin(2.*phi); c2phi = cos(2.*phi)
    zsp(1) = wnint(1)* mt(1) + wnint(2)*(mt(2)+mt(3)) &
           & +wnint(3)*(mt(4)*cphi+mt(5)*sphi) &
           & +wnint(4)*((mt(3)-mt(2))*c2phi-2.*mt(6)*s2phi)
    zsp(2) = wnint(5)*mt(1) + wnint(6)*(mt(2)+mt(3)) &
           & +(wnint(7)+wnint(9))*(mt(4)*cphi+mt(5)*sphi) &
           & +(wnint(8)+wnint(10))*((mt(3)-mt(2))*c2phi-2.*mt(6)*s2phi)
    zsp(3) = (wnint(11)+wnint(13))*(mt(5)*cphi-mt(4)*sphi) &
           & +(wnint(12)+wnint(14))*((mt(3)-mt(2))*s2phi+2.*mt(6)*c2phi)
    end subroutine momentSingleFrequencyDisplacement
!--------------------------------------------------------------------
!  radial derivatives of displacement for moment tensor
!  phi:  receiver azimuth from south in rad
!
    subroutine momentRadDerivSingleFrequencyDisplacement(phi,mt,wnint,zsp)
    double precision :: phi
    double precision, dimension(:) :: mt          !  moment tensor components
    double complex, dimension(:) :: wnint         !  wavenumber integrals vs frequency
    double complex, dimension(:) :: zsp           !  displacement spectrum
    double precision :: sphi,cphi,c2phi,s2phi
!
    sphi = sin(phi); cphi = cos(phi)
    s2phi = sin(2.*phi); c2phi = cos(2.*phi)
    zsp(1) = wnint(15)* mt(1) + wnint(16)*(mt(2)+mt(3)) &
          & +wnint(17)*(mt(4)*cphi+mt(5)*sphi) &
          & +wnint(18)*((mt(3)-mt(2))*c2phi-2.*mt(6)*s2phi)
    zsp(2) = wnint(19)*mt(1) &
           & +wnint(20)*(mt(2)+mt(3)) &
           & +(wnint(21)+wnint(23))*(mt(4)*cphi+mt(5)*sphi) &
           & +(wnint(22)+wnint(24))*((mt(3)-mt(2))*c2phi-2.*mt(6)*s2phi)
    zsp(3) = (wnint(25)+wnint(27))*(mt(5)*cphi-mt(4)*sphi) &
           & +(wnint(26)+wnint(28))*((mt(3)-mt(2))*s2phi+2.*mt(6)*c2phi)
    end subroutine momentRadDerivSingleFrequencyDisplacement
!--------------------------------------------------------------------
!  1/r*theta derivatives of displacement for moment tensor
!  comp: component (Z,R,T)
!  r:    radius of receiver in m
!  phi:  receiver azimuth from south in rad
!
    subroutine momentThetaDerivSingleFrequencyDisplacement(r,phi,mt,wnint,zsp)
    double precision :: r,phi
    double precision, dimension(:) :: mt          !  moment tensor components
    double complex, dimension(:) :: wnint         !  wavenumber integrals vs frequency
    double complex, dimension(:) :: zsp           !  displacement spectrum
    double precision :: sphi,cphi,c2phi,s2phi
!
    sphi = sin(phi); cphi = cos(phi)
    s2phi = sin(2.d0*phi); c2phi = cos(2.d0*phi)
    zsp(1) = wnint(29)* mt(1) + wnint(30)*(mt(2)+mt(3)) &
          & +wnint(31)*(mt(4)*cphi+mt(5)*sphi) &
          & +wnint(32)*((mt(3)-mt(2))*c2phi-2.*mt(6)*s2phi)
    zsp(2) = wnint(33)*mt(1) + wnint(34)*(mt(2)+mt(3)) &
           & +(wnint(35)+wnint(37))*(mt(4)*cphi+mt(5)*sphi) &
           & +(wnint(36)+wnint(38))*((mt(3)-mt(2))*c2phi-2.*mt(6)*s2phi)
    zsp(3) = (wnint(39)+wnint(41))*(mt(5)*cphi-mt(4)*sphi) &
           & +(wnint(40)+wnint(42))*((mt(3)-mt(2))*s2phi+2.*mt(6)*c2phi)
    zsp = zsp/r
    end subroutine momentThetaDerivSingleFrequencyDisplacement
!--------------------------------------------------------------------
!  1/(r*sin(theta))*phi derivatives of displacement for moment tensor
!  comp:  component (Z,R,T)
!  r:     receiver radius in m
!  delta: epicentral distance in rad
!  phi:   receiver azimuth from south
!
    subroutine momentPhiDerivSingleFrequencyDisplacement(r,delta,phi,mt,wnint,zsp)
    double precision :: r,delta,phi
    double precision, dimension(:) :: mt            !  moment tensor components
    double complex, dimension(:) :: wnint           !  wavenumber integrals vs frequency
    double complex, dimension(:) :: zsp             !  displacement spectrum
    double precision :: sphi,cphi,c2phi,s2phi,rst
!
    sphi = sin(phi); cphi = cos(phi)
    s2phi = sin(2.*phi); c2phi = cos(2.*phi)
    rst = r*sin(delta)
!
    zsp(1) = wnint(3)*(-mt(4)*sphi+mt(5)*cphi) &
          & +wnint(4)*(-2.*(mt(3)-mt(2))*s2phi-4.*mt(6)*c2phi)
    zsp(2) = (wnint(7)+wnint(9))*(-mt(4)*sphi+mt(5)*cphi) &
          & +(wnint(8)+wnint(10))*(-2.*(mt(3)-mt(2))*s2phi-4.*mt(6)*c2phi)
    zsp(3) = (wnint(11)+wnint(13))*(-mt(5)*sphi-mt(4)*cphi) &
          & +(wnint(12)+wnint(14))*(+2.*(mt(3)-mt(2))*c2phi-4.*mt(6)*s2phi)
    zsp = zsp/rst
    end subroutine momentPhiDerivSingleFrequencyDisplacement
!-------------------------------------------------------------------------
!  Spherical gradient tensor of displacement field as described
!  in born-scattering.tex
!  zspdr:    radial derivative of displacement vector
!  zspdt:    1/r * delta-derivative of displacement vector
!  zspdf:    1/(r*sin(delta)) * phi-derivative of displacement vector
!  r:        receiver radius in m
!  delta:    epicentral distance in rad
!  gradus:   spherical gradient of displacement (3x3-matrix)
!    
    subroutine sphericalGradientSingleFrequencyDisplacement(zsp,zspdr,zspdt,zspdf,r,delta,gradus)
    double complex, dimension(:) :: zsp,zspdr,zspdt,zspdf
    double complex, dimension(:,:) :: gradus
    double precision :: r,delta,rtandel
!
    rtandel = r*tan(delta)
    gradus(1,1) = zspdr(1)                           ! d_r u_r
    gradus(1,2) = zspdr(2)                           ! d_r u_t
    gradus(1,3) = zspdr(3)                           ! d_r u_f
    gradus(2,1) = zspdt(1)-zsp(2)/r                  ! 1/r d_t u_r - 1/r u_t
    gradus(2,2) = zspdt(2)+zsp(1)/r                  ! 1/r d_t u_t + 1/r u_r
    gradus(2,3) = zspdt(3)                           ! 1/r d_t u_f
    gradus(3,1) = zspdf(1)-zsp(3)/r                  ! 1/(rst) d_f u_r-1/r u_f
    gradus(3,2) = zspdf(2)-zsp(3)/rtandel            ! 1/(rst) d_f u_t-cot(delta)/r u_f
    gradus(3,3) = zspdf(3)+zsp(1)/r+zsp(2)/rtandel   ! 1/rst d_f u_f+1/r u_r +cot(delta)/r u_t
    end subroutine sphericalGradientSingleFrequencyDisplacement
!-------------------------------------------------------------------------
!  Spherical strain components
!  zspdr:    radial derivative of displacement vector
!  zspdt:    1/r * delta-derivative of displacement vector
!  zspdf:    1/(r*sin(delta)) * phi-derivative of displacement vector
!  r:        receiver radius in m
!  delta:    epicentral distance in rad
!  strain:   strain components (6-vector): err,ett,eff,ert,erf,etf
!    
    subroutine sphericalStrainSingleFrequencyDisplacement(zsp,zspdr,zspdt,zspdf,r,delta,strain)
    double complex, dimension(:) :: zsp,zspdr,zspdt,zspdf
    double complex, dimension(:,:) :: strain
    double precision :: r,delta,rtandel
!
    rtandel = r*tan(delta)
    strain(1,1) = zspdr(1)                                ! err = d_r u_r
    strain(2,2) = zspdt(2)+zsp(1)/r                       ! ett = 1/r d_t u_t + 1/r u_r
    strain(3,3) = zspdf(3)+zsp(1)/r+zsp(2)/rtandel        ! eff = 1/rst d_f u_f+1/r u_r+cot(delta)/r u_t
    strain(1,2) = 0.5d0*(zspdt(1)+zspdr(2)-zsp(2)/r)      ! ert = 1/2*(1/r d_t u_r+d_r u_t-u_t/r)
    strain(1,3) = 0.5d0*(zspdf(1)+zspdr(3)-zsp(3)/r)      ! erf = 1/2*(1/rst d_f u_r+d_r u_f-u_f/r)
    strain(2,3) = 0.5d0*(zspdf(2)+zspdt(3)-zsp(3)/rtandel)! etf = 1/2*(1/rst d_f u_t+d_t u_f-cot(del)/r u_f)
    strain(3,2) = strain(2,3)
    strain(3,1) = strain(1,3)
    strain(2,1) = strain(1,2)
    end subroutine sphericalStrainSingleFrequencyDisplacement
!-------------------------------------------------------------------------
!  Spherical stress components from strains
!  strain:   strain components (3x3 matrix)
!  elcon:    complex elastic constants (A,C,F,L,N,kap,mu)
!  stress:   stress tensor (3x3-matrix)
!  Trr = C Err + F (Ett+Eff)
!  Ttt = F Err + A (Ett+Eff) - 2N Eff   
!  Tff = F Err + A (Ett+Eff) - 2N Ett
!  Trt = 2L Ert
!  Trf = 2L Erf
!  Ttf = 2N Etf    
!    
    subroutine sphericalStressSingleFrequencyDisplacement(strain,elcon,stress)
    double complex, dimension(:,:) :: strain,stress
    double complex, dimension(:) :: elcon
!
    stress(1,1) = elcon(2)*strain(1,1)+elcon(3)*(strain(2,2)+strain(3,3))
    stress(2,2) = elcon(3)*strain(1,1)+elcon(1)*(strain(2,2)+strain(3,3))-2.d0*elcon(5)*strain(3,3)
    stress(3,3) = elcon(3)*strain(1,1)+elcon(1)*(strain(2,2)+strain(3,3))-2.d0*elcon(5)*strain(2,2)
    stress(1,2) = 2.d0*elcon(4)*strain(1,2)
    stress(1,3) = 2.d0*elcon(4)*strain(1,3)
    stress(2,3) = 2.d0*elcon(5)*strain(2,3)
    stress(2,1) = stress(1,2)
    stress(3,1) = stress(1,3)
    stress(3,2) = stress(2,3)
    end subroutine sphericalStressSingleFrequencyDisplacement
!
end module singleFrequencyDisplacement

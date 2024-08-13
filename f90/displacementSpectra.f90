! ====================================================================================
!  Compute ZRT displacement and strain spectra from GreenFKSpectra
! ====================================================================================
!----------------------------------------------------------------------------
!   Copyright 2020 Wolfgang Friederich (Ruhr-Universitaet Bochum, Germany)
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
!-----------------------------------------------------------------------------
!  Module with routines to calculate frequency displacement and strain spectra
!  from Green FK spectra either using Bessel of Legendre functions.
!----------------------------------------------------------------------------
 module displacementSpectra
    use singleFrequencyDisplacement
    use legendreFunctions
    use wavenumberIntegrals
    use harmonicDegreeSums
    use greenFKSpectra
    implicit none
    integer, parameter :: NCOMP = 3
    double precision, parameter :: TAPFRAC = 0.2d0
!
contains
!------------------------------------------------------------------------------------------------------------    
!  Compute spherical displacement components using expansion into spherical harmonics
!  gfk: GreenFkSpectra object with meta information
!  gfkspr: real part of GFK-spectra for fixed radial node
!  gfkspi: imag part of GFK-spectra for fixed radial node
!  fomt: array with either force or moment tensor components
!  alf: object with Legendre functions versus harmonic degree evaluated at receiver epicentral coordinates
!  re, epidis, phi: epicentral receiver coordinates (units: m, rad, rad)
!  propdir: propagation direction measured at receiver position (different from phi)
!  zsp: array with displacement spectra versus frequency
!  errmsg: error message object
!  strainsp: array with strain spectra versus frequency (optional)
!
    subroutine sphericalHarmonicsDisplacementSpectra(gfk,gfkspr,gfkspi,fomt,alf,re,delta,phi,zsp,errmsg,strainsp)
    type (green_fk_spectra) :: gfk
    real, dimension(:,:) :: gfkspr,gfkspi
    double precision, dimension(:) :: fomt
    type (legendre_functions) :: alf
    double precision :: re,delta,phi
    double complex, dimension(:,:), allocatable :: zsp
    double complex, dimension(:,:,:), optional, allocatable :: strainsp
    type (error_message) :: errmsg
    integer :: nwe,nwa,jf,j,nwint,isp,jsp,offset
    integer, dimension(:), allocatable :: kint
    double complex, dimension(:), allocatable :: hdsum
!
    allocate(zsp(gfk%nf2,NCOMP))
    zsp = 0.d0
    if (present(strainsp)) then
       allocate(strainsp(gfk%nf2,NCOMP,NCOMP))
       strainsp = 0.d0
    endif
!
!  get array of harmonic degree sum indices as required by dsvmask and istyp
!   
    call getIndexArrayHarmonicDegreeSums(gfk,nwint,kint)
    allocate(hdsum(maxIndexHarmonicDegreeSums(gfk%istyp)))
!
!  start frequency loop to compute displacement vs frequency
!  must be within event loop because source depth may change   
!
    nwe = 0
    do jf = gfk%nf1,gfk%nf2
       nwa = nwe+1                        ! select wavenumber range for frequency
       nwe = nwe+gfk%nwn(jf)
    !
    !  compute harmonic degree sums, skip if dsvmask(isp) = 0
    !
       do j = 1,nwint
          hdsum(kint(j)) = 0.d0           ! zero value in case dsvmask(isp) = 0
          if (gfk%istyp == 1) then
             isp = moment_component_jump_from_harmonic_degree_sum(kint(j),1)
             if (gfk%dsvmask(isp) == 0) cycle
             jsp = moment_component_jump_from_harmonic_degree_sum(kint(j),2)
             offset = positionDataGreenFKSpectra(gfk,isp,jsp)
             call momentHarmonicDegreeSums(gfkspr(nwa:nwe,offset),gfkspi(nwa:nwe,offset),&
                  & kint(j),alf,TAPFRAC,hdsum(kint(j)))
          else
             isp = force_component_jump_from_harmonic_degree_sum(kint(j),1)
             if (gfk%dsvmask(isp) == 0) cycle
             jsp = force_component_jump_from_harmonic_degree_sum(kint(j),2)
             offset = positionDataGreenFKSpectra(gfk,isp,jsp)
             call forceHarmonicDegreeSums(gfkspr(nwa:nwe,offset),gfkspi(nwa:nwe,offset),&
                  & kint(j),alf,TAPFRAC,hdsum(kint(j)))
          endif
       enddo                              ! end harmonic degree sums loop
    !
    !  evaluate single frequency displacements
    !
       if (present(strainsp)) then
          call evalSingleFrequencyDisplacementSpectra(gfk%istyp,re,delta,phi,hdsum,fomt,zsp(jf,:),strainsp(jf,:,:))
       else
          call evalSingleFrequencyDisplacementSpectra(gfk%istyp,re,delta,phi,hdsum,fomt,zsp(jf,:))
       endif
    enddo
!
    deallocate(hdsum,kint)
    end subroutine sphericalHarmonicsDisplacementSpectra
!------------------------------------------------------------------------------------------------------------    
!  Compute spherical displacement components using expansion into Bessel functions
!  gfk: GreenFkSpectra object with meta information
!  gfkspr: real part of GFK-spectra for fixed radial node
!  gfkspi: imag part of GFK-spectra for fixed radial node
!  bessel: array with Bessel functions versus wavenumber evaluated at receiver epicentral distance
!  re, delta, phi: epicentral receiver coordinates (units: m, rad, rad)
!  propdir: propagation direction measured at receiver position (different from phi)
!  zsp: array with displacement spectra versus frequency
!  errmsg: error message object
!  strainsp: array with strain spectra versus frequency (optional)
!
    subroutine besselDisplacementSpectra(gfk,gfkspr,gfkspi,fomt,besselj,re,delta,phi,zsp,errmsg,strainsp)
    type (green_fk_spectra) :: gfk
    real, dimension(:,:) :: gfkspr,gfkspi
    double precision, dimension(:) :: fomt
    double precision, dimension(:,0:) :: besselj    
    double precision :: re,delta,phi
    double complex, dimension(:,:), allocatable :: zsp
    double complex, dimension(:,:,:), optional, allocatable :: strainsp
    type (error_message) :: errmsg
    integer :: nwe,nwa,jf,j,nwint,isp,jsp,offset
    integer, dimension(:), allocatable :: kint
    double complex, dimension(:), allocatable :: wnint
!
    allocate(zsp(gfk%nf2,NCOMP))
    zsp = 0.d0
    if (present(strainsp)) then
       allocate(strainsp(gfk%nf2,NCOMP,NCOMP))
       strainsp = 0.d0
    endif
!
!  get array of harmonic degree sum indices as required by dsvmask and istyp
!   
    call getIndexArrayWavenumberIntegrals(gfk,nwint,kint)
    allocate(wnint(maxIndexWavenumberIntegrals(gfk%istyp)))
!
!  start frequency loop to compute displacement vs frequency
!  must be within event loop because source depth may change   
!
    nwe = 0
    do jf = gfk%nf1,gfk%nf2
       nwa = nwe+1                        ! select wavenumber range for frequency
       nwe = nwe+gfk%nwn(jf)
    !
    !  compute harmonic degree sums, skip if dsvmask(isp) = 0
    !
       do j = 1,nwint
          wnint(kint(j)) = 0.d0          ! zero value in case dsvmask(isp) = 0
          if (gfk%istyp == 1) then
             isp = moment_component_jump_from_wavenumber_integral(kint(j),1)
             if (gfk%dsvmask(isp) == 0) cycle
             jsp = moment_component_jump_from_wavenumber_integral(kint(j),2)
             offset = positionDataGreenFKSpectra(gfk,isp,jsp)
             call momentWavenumberIntegrals(gfk,gfkspr(nwa:nwe,offset),gfkspi(nwa:nwe,offset),kint(j),&
                  delta*gfk%rearth,besselj(:,0:),TAPFRAC,wnint(kint(j)))
          else
             isp = force_component_jump_from_wavenumber_integral(kint(j),1)
             if (gfk%dsvmask(isp) == 0) cycle
             jsp = force_component_jump_from_wavenumber_integral(kint(j),2)
             offset = positionDataGreenFKSpectra(gfk,isp,jsp)
             call forceWavenumberIntegrals(gfk,gfkspr(nwa:nwe,offset),gfkspi(nwa:nwe,offset),kint(j),&
                  delta*gfk%rearth,besselj(:,0:),TAPFRAC,wnint(kint(j)))
          endif
       enddo                      ! end harmonic degree sums loop
    !
    !  evaluate single frequency displacements
    !
       if (present(strainsp)) then
          call evalSingleFrequencyDisplacementSpectra(gfk%istyp,re,delta,phi,wnint,fomt,zsp(jf,:),strainsp(jf,:,:))
       else
          call evalSingleFrequencyDisplacementSpectra(gfk%istyp,re,delta,phi,wnint,fomt,zsp(jf,:))
       endif
    enddo
!
    deallocate(wnint,kint)
    end subroutine besselDisplacementSpectra
!-------------------------------------------------------------------------------------------------------
!  Evaluate single frequency displacement
!  istyp :: source type
!  re,delta,phi: epicentral coordinates
!  hdsum: harmonic degree sums or wavenumber integrals
!  fomt: force or moment tensor elements
!  usp: displacement spectra components
!  strainsp: strain spectra (optional)    
!
    subroutine evalSingleFrequencyDisplacementSpectra(istyp,re,delta,phi,hdsum,fomt,usp,strainsp)
    integer :: istyp
    double precision :: re,delta,phi
    double complex, dimension(:) :: hdsum
    double precision, dimension(:) :: fomt
    double complex, dimension(:) :: usp
    double complex, dimension(:,:), optional :: strainsp
    double complex, dimension(NCOMP) :: uspdr,uspdt,uspdf
!     
    select case (istyp)
    case (0)
       call forceSingleFrequencyDisplacement(phi,fomt,hdsum,usp)
       if (present(strainsp)) then
          call forceRadDerivSingleFrequencyDisplacement(phi,fomt,hdsum,uspdr)
          call forceThetaDerivSingleFrequencyDisplacement(re,phi,fomt,hdsum,uspdt)
          call forcePhiDerivSingleFrequencyDisplacement(re,delta,phi,fomt,hdsum,uspdf)
       endif
    case (1)
       call momentSingleFrequencyDisplacement(phi,fomt,hdsum,usp)
       if (present(strainsp)) then
          call momentRadDerivSingleFrequencyDisplacement(phi,fomt,hdsum,uspdr)
          call momentThetaDerivSingleFrequencyDisplacement(re,phi,fomt,hdsum,uspdt)
          call momentPhiDerivSingleFrequencyDisplacement(re,delta,phi,fomt,hdsum,uspdf)
       endif
    end select
    if (present(strainsp)) then
       call sphericalStrainSingleFrequencyDisplacement(usp,uspdr,uspdt,uspdf,re,delta,strainsp)
    endif
    end subroutine evalSingleFrequencyDisplacementSpectra
!
end module displacementSpectra

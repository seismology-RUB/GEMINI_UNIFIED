! ====================================================================
! Computes source terms for single force and moment tensor components
! ====================================================================
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
!--------------------------------------------------------------------------
!  Evaluates expressions for source terms as described in gemini-displ.tex.
!  Nodes do not fall on layer boundaries. See externalRadialNodes.f90.
!--------------------------------------------------------------------------
module sourceTerms
    use string
    use geminiIntegrationEnvironment
    use nodeEarthmodel
    use splineEarthmodelCoefficients
    use externalRadialNodes
    implicit none
!
contains
!--------------------------------------------------------------------------
!  Return number of source terms depending on source type
!!  sourcetype:     source type (in)
!!  njs:            number of spheroidal source vectors (out)
!!  njt:            number of toroidal source vectors (out)
!
  subroutine getNumberSourceTerms(sourcetype,state_src,njs,njt)
     character (len=*) :: sourcetype
     integer :: njs,njt,state_src
     if (sourcetype.equal.'FORCE') then
        njs = 2; njt = 1
     else
        njs = 4; njt = 2
     endif
  end subroutine getNumberSourceTerms
!-----------------------------------------------------------------------------------
!  Calculate source terms for either force or moment source at one
!  external node
!  sourcetype:          type of source (in)
!  rs:                  radius of source node (in)
!  state:               fluid (1) or solid (0) at source
!  nl:                  layer of source
!  zsph:                spheroidal source terms (out)
!  ztor:                toroidal source terms (out)
!
  subroutine oneNodeSourceTerms(sourcetype,elp1,om,rs,state,nl,zsph,ztor)
     character(len=*) :: sourcetype
     double precision :: rs,elp1,om
     integer :: state,nl
     double complex, dimension(:,:), allocatable :: zsph    ! (comp,jump)
     double complex, dimension(:,:), allocatable :: ztor    ! (comp,jump)
     double complex, dimension(7) :: zelcon
     double precision :: ro,rr2,rr3
   !
     rr2 = 1./(rs*rs)
     rr3 = rr2/rs
   !
   !  jump vectors for moment tensor source
   !
     if (sourcetype.equal.'MOMENT') then
        if (state == 1) then                                                   ! moment in ocean
           allocate(zsph(2,2))
           call evalComplexElasticConstantsSplineEarthmodel(rs,nl,ro,zelcon)
           zsph(1,1) = rr2*(1.d0/zelcon(2)-elp1*rr2/(om**2*ro))                ! m=0, Mrr, gamma_l not included
           zsph(2,1) = rr3*2.d0
           zsph(1,2) = 0.5*elp1*rr2*rr2/(om**2*ro)                             ! m=0, basis solution for (Mtt+Mff), gamm_l not included
           zsph(2,2) = -rr3                           
        else                                                                   ! moment in solid
           allocate(zsph(4,4))
           allocate(ztor(2,2))
           ztor = dcmplx(0.d0,0.d0)
           zsph = dcmplx(0.d0,0.d0)
           call evalComplexElasticConstantsSplineEarthmodel(rs,nl,ro,zelcon)
           zsph(1,1) = rr2/zelcon(2)                  ! m=0, basis solution for Mrr
           zsph(2,1) = rr3*2.d0*zelcon(3)/zelcon(2)
           zsph(4,1) = -rr3*zelcon(3)/zelcon(2)
           zsph(2,2) = -rr3                           ! m=0, basis solution for (Mtt+Mff)
           zsph(4,2) = 0.5*rr3
           zsph(3,3) = 0.5*rr2/zelcon(4)       ! m=+-1, basis solution for +-Mrt + iMrf, sqrt(1/[l(l+1)]) not included!
           zsph(4,4) = 0.25*rr3                ! m=+-2, basis solution for Mff-Mtt +- 2iMtf, sqrt[(l+2)(l-1)]/[l(l+1)]] is not included
           ztor(1,1) = 0.5*rr2/zelcon(4)       ! m=+-1, basis solution for +-Mrf + iMrt, sqrt(1/[l(l+1)]) not included!
           ztor(2,2) = 0.25*rr3                ! m=+-2, basis solution for 2Mtf -+ i(Mff-Mtt), sqrt[(l+2)(l-1)]/[l(l+1)]] is not included
        endif
   !
   !  Jump vectors for single force source.
   !
     else
        if (state == 1) then                       ! force in ocean
           allocate(zsph(2,2))
           call evalComplexElasticConstantsSplineEarthmodel(rs,nl,ro,zelcon)
           zsph(1,1) = 0.d0                        ! m=0, Fr, gamma_l extracted
           zsph(2,1) = -rr2
           zsph(1,2) = elp1*rr3/(om**2*ro)         ! m=1, +-Ft -i*Ff, gamma_l*sqrt(1/[l(l+1)] extracted
           zsph(2,2) = 0.d0
        else                                       ! force in solid
           allocate(zsph(4,2))
           allocate(ztor(2,1))
           zsph = dcmplx(0.d0,0.d0)
           ztor = dcmplx(0.d0,0.d0)
           zsph(2,1) = -rr2                        ! m=0, Fr
           zsph(4,2) = 0.5*rr2                     ! m=+-1,  gamma_l*sqrt(1/[l(l+1)] extracted
           ztor(2,1) = 0.5*rr2                     ! m=+-1,  -+Ff-i*Ft, gamma_l*sqrt(1/[l(l+1)] extracted
        endif
     endif
  end subroutine oneNodeSourceTerms
!
end module sourceTerms

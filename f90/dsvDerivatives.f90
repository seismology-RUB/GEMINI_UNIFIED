! ============================================================================
!  Compute radial derivatives of FK Green functions
! ============================================================================
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
!----------------------------------------------------------------------------
module dsvDerivatives
    use geminiIntegrationEnvironment
    use baseOdeSystem
    use nodeEarthmodel
    implicit none
!
contains
!-------------------------------------------------------------------------------------------
!  compute derivatives of Green functions
!  gf:      array of rank 2 with Green function components U,R,V,S or W,T or U,R
!           for different jump vectors
!  r:       radius where derivative is evaluated
!  intenv:  integration environment (same as linked with sode)
!  sode:    differential equation system used to calculate derivatives
!  gfd:     derivatives of U,V or U or W for different jump vectors
!
     subroutine computeDsvDerivatives(gf,r,intenv,sode,gfd,errmsg)
     double complex, dimension(:,:) :: gf
     double precision :: r
     type (gemini_integration_environment) :: intenv
     class (base_ode_system) :: sode
     double complex, dimension(:,:) :: gfd
     type (error_message) :: errmsg
     integer :: top,nl,ios,nc,nj,nvar
     double complex, dimension(:,:), allocatable :: a
     character(len=21) :: myname = 'computeDsvDerivatives'
!
     nvar = sode%getNvar()
     nc = size(gf,1)
     nj = size(gf,2)
     if (nc .ne. nvar) then
        call add(errmsg,2,'number of DSV-components do not fit to ODE',myname)
        return
     endif
     if (size(gfd,1) .ne. nc) then
        call add(errmsg,2,'incorrect number of components for derivatives',myname)
        return
     endif
     if (size(gfd,2) .ne. nj) then
        call add(errmsg,2,'incorrect number of jump vectors for derivatives',myname)
        return
     endif
     allocate(a(nvar,nvar))
     call getLayerIndexNodeEarthmodel(intenv%nem,r,nl,top)
     call intenv%setInterval(nl)
     call sode%getSysmatComplex(r,a,ios)         ! returns sytem matrix
     gfd(1:nvar,:) = matmul(a,gf)
     deallocate(a)
     end subroutine computeDsvDerivatives
!
end module dsvDerivatives

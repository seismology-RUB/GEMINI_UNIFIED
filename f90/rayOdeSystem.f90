! ===============================================================================
!  Implementation of the ODE system for ray tracing in spherically symmetric media
! ===============================================================================
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
!-------------------------------------------------------------
!  ODE system for integration of ray equations.
!  Extends the "baseOdeSystem" class. Also has a pointer to the
!  integration environment. A copy of the integration environment
!  is not sufficient here because the environment (layer index)
!  changes dynamically elsewhere.
!------------------------------------------------------------
!-------------------------------------------------------------
!> \brief System for integration of minors
!------------------------------------------------------------
module rayOdeSystem
    use baseOdeSystem
    use rayIntegrationEnvironment
    use splineEarthmodelCoefficients                 ! global module with precalculated spline coefficients
    use errorMessage
    implicit none
    type, extends(base_ode_system) :: ray_ode_system
       private
       type (ray_integration_environment), pointer :: intenv => null()    ! pointer to integration environment
    contains
        procedure :: getRhsComplex => getRhsComplexRayOdeSystem
        procedure :: getRhsReal => getRhsRealRayOdeSystem
        procedure :: getSysmatComplex => getSysmatComplexRayOdeSystem
        procedure :: getSysmatReal => getSysmatRealRayOdeSystem
        procedure :: derivativeReal => derivativeRealRayOdeSystem
        procedure :: derivativeComplex => derivativeComplexRayOdeSystem
        procedure :: createRayOdeSystem
        procedure :: deallocRayOdeSystem
    end type ray_ode_system
contains
!-----------------------------------------------------------------
!  Create ray ODE system
!
    subroutine createRayOdeSystem(this,intenv,errmsg)
    class (ray_ode_system) :: this
    type (ray_integration_environment), target :: intenv
    type (error_message) :: errmsg
!
    call addTrace(errmsg,'createRayOdeSystem')
    this%nvar = 2
    this%vartyp = 'R'
    this%linear = .false.
    this%intenv => intenv
    end subroutine createRayOdeSystem
!----------------------------------------------------------------
!> \brief Deallocate
!
    subroutine deallocRayOdeSystem(this)
    class (ray_ode_system) :: this
    if (associated(this%intenv)) nullify(this%intenv)
    end subroutine deallocRayOdeSystem
!-------------------------------------------------------------------
!> \brief Get system matrix (not used)
!! x:      location where rhs of ODE system is evaluated (in)
!! sysmat: matrix of ODE system (out)
!! ios:    error indicator (out)
!
    subroutine getSysmatComplexRayOdeSystem(this,x,sysmat,ios)
    class (ray_ode_system) :: this
    double precision :: x
    double complex, dimension(:,:) :: sysmat
    integer :: ios
    sysmat = 0.d0; ios = 0
    end subroutine getSysmatComplexRayOdeSystem
!-------------------------------------------------------------------
!> \brief Get right hand side (complex, not used)
!! x:   location where rhs of ODE system is evaluated (in)
!! rhs: right hand side of ODE system (out)
!! ios: error indicator (out)
!! This routine does nothing here.
!
    subroutine getRhsComplexRayOdeSystem(this,x,rhs,ios)
    class (ray_ode_system) :: this
    double precision :: x
    double complex, dimension(:) :: rhs
    integer :: ios
    rhs = dcmplx(0.d0); ios = 0
    end subroutine getRhsComplexRayOdeSystem
!-------------------------------------------------------------------
!> \brief Get right hand side (real, not used)
!! x:   location where rhs of ODE system is evaluated (in)
!! rhs: right hand side of ODE system (out)
!! ios: error indicator (out)
!! This routine does nothing here.
!
    subroutine getRhsRealRayOdeSystem(this,x,rhs,ios)
    class (ray_ode_system) :: this
    double precision :: x
    double precision, dimension(:) :: rhs
    integer :: ios
    rhs = 0.d0; ios = 0
    end subroutine getRhsRealRayOdeSystem
!-------------------------------------------------------------------
!> \brief Get system matrix (real, not used)
!! x:      location where rhs of ODE system is evaluated (in)
!! sysmat: matrix of ODE system (out)
!! ios:    error indicator (out)
!
    subroutine getSysmatRealRayOdeSystem(this,x,sysmat,ios)
    class (ray_ode_system) :: this
    double precision :: x
    double precision, dimension(:,:) :: sysmat
    integer :: ios
    sysmat = 0.d0; ios = 0
    end subroutine getSysmatRealRayOdeSystem
!--------------------------------------------------------------------
!> \brief Get derivatives (complex, not used)
!! x:      location where rhs of ODE system is evaluated (in)
!! y:      current function vector (in)
!! dy:     derivative vector (out)
!! ios:    error indicator (out)
!
    subroutine derivativeComplexRayOdeSystem(this,x,y,dy,ios)
    class (ray_ode_system) :: this
    double precision :: x
    double complex, dimension(:) :: y,dy
    integer :: ios
    dy = 0.d0; ios = 0
    end subroutine derivativeComplexRayOdeSystem
!--------------------------------------------------------------------
!> \brief Get derivatives (real, used)
!! r:      current radius of ray (in)
!! y:      current function vector (in), (delta(rad) and T(s))
!! dy:     derivative vector (out)
!! ios:    error indicator (out)
!
    subroutine derivativeRealRayOdeSystem(this,x,y,dy,ios)
    class (ray_ode_system) :: this
    double precision :: x
    double precision, dimension(:) :: y,dy
    integer :: ios
    double precision :: p,q,ur2,v
!
    p = this%intenv%slowness                ! slowness of ray in s/rad (ray parameter)
    if (trim(this%intenv%raytype) == 'P') then
       v = getPVelocitySplineEarthmodel(this%intenv%nl,x)
    else if (trim(this%intenv%raytype) == 'S') then
       v = getSVelocitySplineEarthmodel(this%intenv%nl,x)
    else
       ios = 1
       return
    endif
!
    ur2 = (x/v)**2
    q = sqrt(ur2-p**2)
    dy(1) = p/(q*x)               ! dDelta/dr
    dy(2) = ur2/(q*x)             ! dT/dr
    ios = 0
!
    end subroutine derivativeRealRayOdeSystem
!
end module rayOdeSystem

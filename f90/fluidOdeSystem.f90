! ======================================================================
!  Implementation of the ODE system for a fluid medium
! ======================================================================
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
!-------------------------------------------------------------
!  ODE system for integration DSV-equations in fluid medium.
!  Extends the "baseOdeSystem" class. Also has a pointer to the
!  integration environment. A copy of the integration environment
!  is not sufficient here because the environment (layer index)
!  changes dynamically elsewhere.
!------------------------------------------------------------
module fluidOdeSystem
    use geminiIntegrationEnvironment
    use splineEarthmodelCoefficients                     ! global module with precalculated spline coefficients
    use errorMessage
    use baseOdeSystem
    implicit none
    type, extends(base_ode_system) :: fluid_ode_system
       private
       type (gemini_integration_environment), pointer :: intenv => null()  ! pointer important here because
    contains                                                               ! members will be updated elsewhere
        procedure :: getRhsComplex => getRhsComplexFluid
        procedure :: getRhsReal => getRhsRealFluid
        procedure :: getSysmatComplex => getSysmatComplexFluid
        procedure :: getSysmatReal => getSysmatRealFluid
        procedure :: derivativeReal => derivativeRealFluid
        procedure :: derivativeComplex => derivativeComplexFluid
        procedure :: createFluid
        procedure :: deallocFluid
        procedure :: getHorcompComplexFluid
    end type fluid_ode_system
contains
!-----------------------------------------------------------------
!> \brief Create fluid ODE system
!! intenv:   object representing integration environment
!! errmsg:   error messsage object
!
    subroutine createFluid(this,intenv,errmsg)
    class (fluid_ode_system) :: this
    type (gemini_integration_environment), target :: intenv
    type (error_message) :: errmsg
!
    call addTrace(errmsg,'createFluid')
    this%nvar = 2
    this%vartyp = 'C'
    this%linear = .true.
    this%intenv => intenv            ! pointer important here
    end subroutine createFluid
!----------------------------------------------------------------
!> \brief Deallocate
!
    subroutine deallocFluid(this)
    class (fluid_ode_system) :: this
    if (associated(this%intenv)) nullify(this%intenv)  ! only nullify, no dellaoc
    end subroutine deallocFluid
!-------------------------------------------------------------------
!> \brief Get system matrix (complex)
!! x:      location where rhs of ODE system is evaluated (in)
!! sysmat: matrix of ODE system (out)
!! ios:    error indicator (out)
!
    subroutine getSysmatComplexFluid(this,x,sysmat,ios)
    class (fluid_ode_system) :: this
    double precision :: x
    double complex, dimension(:,:) :: sysmat
    integer :: ios
    double complex, dimension(7) :: zelcon
    double complex:: zom,zomro
    double precision :: dll1,ro
!
    call evalComplexElasticConstantsSplineEarthmodel(x,this%intenv%nl,ro,zelcon)
    zom = dcmplx(this%intenv%omre,this%intenv%omim)
    dll1 = this%intenv%dll1
    zomro = zom*zom*ro
    sysmat(1,1) = -2.d0/x
    sysmat(1,2) = 1.d0/zelcon(2) - dll1/(zomro*x*x)
    sysmat(2,1) = -zomro
    sysmat(2,2) = 0.d0
    end subroutine getSysmatComplexFluid
!-------------------------------------------------------------------
!> \brief Get right hand side (complex)
!! x:   location where rhs of ODE system is evaluated (in)
!! rhs: right hand side of ODE system (out)
!! ios: error indicator (out)
!! This routine does nothing here.
!
    subroutine getRhsComplexFluid(this,x,rhs,ios)
    class (fluid_ode_system) :: this
    double precision :: x
    double complex, dimension(:) :: rhs
    integer :: ios
    rhs = dcmplx(0.d0)
    ios = 0
    end subroutine getRhsComplexFluid
!-------------------------------------------------------------------
!> \brief Get right hand side (real, not used)
!! x:   location where rhs of ODE system is evaluated (in)
!! rhs: right hand side of ODE system (out)
!! ios: error indicator (out)
!! This routine does nothing here.
!
    subroutine getRhsRealFluid(this,x,rhs,ios)
    class (fluid_ode_system) :: this
    double precision :: x
    double precision, dimension(:) :: rhs
    integer :: ios
    rhs = 0.d0; ios = 0
    end subroutine getRhsRealFluid
!-------------------------------------------------------------------
!> \brief Get system matrix (real, not used)
!! x:      location where rhs of ODE system is evaluated (in)
!! sysmat: matrix of ODE system (out)
!! ios:    error indicator (out)
!
    subroutine getSysmatRealFluid(this,x,sysmat,ios)
    class (fluid_ode_system) :: this
    double precision :: x
    double precision, dimension(:,:) :: sysmat
    integer :: ios
    sysmat = 0.d0; ios = 0
    end subroutine getSysmatRealFluid
!--------------------------------------------------------------------
!> \brief Get derivatives (complex, not used)
!! x:      location where rhs of ODE system is evaluated (in)
!! y:      current function vector (in)
!! dy:     derivative vector (out)
!! ios:    error indicator (out)
!
    subroutine derivativeComplexFluid(this,x,y,dy,ios)
    class (fluid_ode_system) :: this
    double precision :: x
    double complex, dimension(:) :: y,dy
    integer :: ios
    dy = 0.d0; ios = 0
    end subroutine derivativeComplexFluid
!--------------------------------------------------------------------
!> \brief Get derivatives (real, not used)
!! x:      location where rhs of ODE system is evaluated (in)
!! y:      current function vector (in)
!! dy:     derivative vector (out)
!! ios:    error indicator (out)
!
    subroutine derivativeRealFluid(this,x,y,dy,ios)
    class (fluid_ode_system) :: this
    double precision :: x
    double precision, dimension(:) :: y,dy
    integer :: ios
    dy = 0.d0; ios = 0
    end subroutine derivativeRealFluid
!----------------------------------------------------------
!  compute horizontal component and derivative
!  V = -R/(om**2*r)
!  dV/dr = -1/r*U+1/r*V = -1/r*(U-V)
!
    subroutine getHorcompComplexFluid(this,x,y,v,vr)
    class (fluid_ode_system) :: this
    double precision :: x
    double complex, dimension(:,:) :: y
    double complex, dimension(:) :: v,vr
    double complex, dimension(7) :: zelcon
    double complex:: zom,zomro
    double precision :: ro
    integer :: nl,top
!
    call getLayerIndexNodeEarthmodel(this%intenv%nem,x,nl,top)
    call evalComplexElasticConstantsSplineEarthmodel(x,nl,ro,zelcon)
    zom = dcmplx(this%intenv%omre,this%intenv%omim)
    zomro = zom*zom*ro
!
    v = -y(2,:)/(x*zomro)
    vr = -1.d0/x*(y(1,:)-v(:))
    end subroutine getHorcompComplexFluid
!
end module fluidOdeSystem

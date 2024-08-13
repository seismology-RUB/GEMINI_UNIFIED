! ======================================================================
!  Implementation of the ODE system for toroidal motions
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
!  ODE system for integration of toroidal motions
!  Extends the "baseOdeSystem" class. Also has a pointer to the
!  integration environment. A copy of the integration environment
!  is not sufficient here because the environment (layer index)
!  changes dynamically elsewhere.
!------------------------------------------------------------
!-------------------------------------------------------------
!> \brief System for integration of minors
!------------------------------------------------------------
module toroidalOdeSystem
    use baseOdeSystem
    use geminiIntegrationEnvironment
    use splineEarthmodelCoefficients                     ! global module with precalculated spline coefficients
    use errorMessage
    implicit none
    type, extends(base_ode_system) :: toroidal_ode_system
       private
       type (gemini_integration_environment), pointer :: intenv => null()    ! pointer to integration environment
    contains
        procedure :: getRhsComplex => getRhsComplexToroidal
        procedure :: getRhsReal => getRhsRealToroidal
        procedure :: getSysmatComplex => getSysmatComplexToroidal
        procedure :: getSysmatReal => getSysmatRealToroidal
        procedure :: derivativeReal => derivativeRealToroidal
        procedure :: derivativeComplex => derivativeComplexToroidal
        procedure :: createToroidal
        procedure :: deallocToroidal
    end type toroidal_ode_system
contains
!-----------------------------------------------------------------
!> \brief Create toroidal ODE system
!! x:      location where rhs of ODE system is evaluated (in)
!! sysmat: matrix of ODE system (out)
!! ios:    error indicator (out)
!
    subroutine createToroidal(this,intenv,errmsg)
    class (toroidal_ode_system) :: this
    type (gemini_integration_environment), target :: intenv
    type (error_message) :: errmsg
!
    call addTrace(errmsg,'createToroidal')
    this%nvar = 2
    this%vartyp = 'C'
    this%linear = .true.
    this%intenv => intenv
    end subroutine createToroidal
!----------------------------------------------------------------
!> \brief Deallocate
!
    subroutine deallocToroidal(this)
    class (toroidal_ode_system) :: this
    if (associated(this%intenv)) nullify(this%intenv)
    end subroutine deallocToroidal
!-------------------------------------------------------------------
!> \brief Get system matrix
!! x:      location where rhs of ODE system is evaluated (in)
!! sysmat: matrix of ODE system (out)
!! ios:    error indicator (out)
!
    subroutine getSysmatComplexToroidal(this,x,sysmat,ios)
    class (toroidal_ode_system) :: this
    double precision :: x
    double complex, dimension(:,:) :: sysmat
    integer :: ios
    double complex, dimension(7) :: zelcon
    double complex :: zom,zomro
    double precision :: dll1,rr,rr2,ro
!
    call evalComplexElasticConstantsSplineEarthmodel(x,this%intenv%nl,ro,zelcon)
!
    zom = dcmplx(this%intenv%omre,this%intenv%omim)
    dll1 = this%intenv%dll1
    zomro = zom*zom*ro
    rr = 1.d0/x
    rr2 = rr*rr
    sysmat(1,1) = rr
    sysmat(1,2) = 1.d0/zelcon(4)
    sysmat(2,1) = -zomro-zelcon(5)*rr2*(2.d0-dll1)
    sysmat(2,2) = -3.*rr
    ios = 0
    end subroutine getSysmatComplexToroidal
!-------------------------------------------------------------------
!> \brief Get right hand side
!! x:   location where rhs of ODE system is evaluated (in)
!! rhs: right hand side of ODE system (out)
!! ios: error indicator (out)
!! This routine does nothing here.
!
    subroutine getRhsComplexToroidal(this,x,rhs,ios)
    class (toroidal_ode_system) :: this
    double precision :: x
    double complex, dimension(:) :: rhs
    integer :: ios
    rhs = dcmplx(0.d0); ios = 0
    end subroutine getRhsComplexToroidal
!-------------------------------------------------------------------
!> \brief Get right hand side (real, not used)
!! x:   location where rhs of ODE system is evaluated (in)
!! rhs: right hand side of ODE system (out)
!! ios: error indicator (out)
!! This routine does nothing here.
!
    subroutine getRhsRealToroidal(this,x,rhs,ios)
    class (toroidal_ode_system) :: this
    double precision :: x
    double precision, dimension(:) :: rhs
    integer :: ios
    rhs = 0.d0; ios = 0
    end subroutine getRhsRealToroidal
!-------------------------------------------------------------------
!> \brief Get system matrix (real, not used)
!! x:      location where rhs of ODE system is evaluated (in)
!! sysmat: matrix of ODE system (out)
!! ios:    error indicator (out)
!
    subroutine getSysmatRealToroidal(this,x,sysmat,ios)
    class (toroidal_ode_system) :: this
    double precision :: x
    double precision, dimension(:,:) :: sysmat
    integer :: ios
    sysmat = 0.d0; ios = 0
    end subroutine getSysmatRealToroidal
!--------------------------------------------------------------------
!> \brief Get derivatives (complex, not used)
!! x:      location where rhs of ODE system is evaluated (in)
!! y:      current function vector (in)
!! dy:     derivative vector (out)
!! ios:    error indicator (out)
!
    subroutine derivativeComplexToroidal(this,x,y,dy,ios)
    class (toroidal_ode_system) :: this
    double precision :: x
    double complex, dimension(:) :: y,dy
    integer :: ios
    dy = 0.d0; ios = 0
    end subroutine derivativeComplexToroidal
!--------------------------------------------------------------------
!> \brief Get derivatives (real, not used)
!! x:      location where rhs of ODE system is evaluated (in)
!! y:      current function vector (in)
!! dy:     derivative vector (out)
!! ios:    error indicator (out)
!
    subroutine derivativeRealToroidal(this,x,y,dy,ios)
    class (toroidal_ode_system) :: this
    double precision :: x
    double precision, dimension(:) :: y,dy
    integer :: ios
    dy = 0.d0; ios = 0
    end subroutine derivativeRealToroidal
!
end module toroidalOdeSystem

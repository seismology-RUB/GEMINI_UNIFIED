! ==================================================================================
!  Implementation of negative transposed (adjoint) ODE system for spheroidal motions
! ==================================================================================
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
!  ODE system for integration of negative transposed system.
!  Extends the "baseOdeSystem" class. Also has a pointer to the
!  integration environment. A copy of the integration environment
!  is not sufficient here because the environment (layer index)
!  changes dynamically elsewhere.
!------------------------------------------------------------
module solidtOdeSystem
    use baseOdeSystem
    use geminiIntegrationEnvironment
    use splineEarthmodelCoefficients                     ! global module with precalculated spline coefficients
    use errorMessage
    implicit none
    type, extends(base_ode_system) :: solidt_ode_system
       private
       type (gemini_integration_environment), pointer :: intenv => null()    ! pointer to integration environment
    contains
        procedure :: getRhsComplex => getRhsComplexSolidt
        procedure :: getRhsReal => getRhsRealSolidt
        procedure :: getSysmatComplex => getSysmatComplexSolidt
        procedure :: getSysmatReal => getSysmatRealSolidt
        procedure :: derivativeReal => derivativeRealSolidt
        procedure :: derivativeComplex => derivativeComplexSolidt
        procedure :: createSolidt
        procedure :: deallocSolidt
    end type solidt_ode_system
contains
!-----------------------------------------------------------------
!> \brief Create solidt ODE system
!! intenv:   object representing integration environment
!! errmsg:   error messsage object
!
  subroutine createSolidt(this,intenv,errmsg)
     class (solidt_ode_system) :: this
     type (gemini_integration_environment), target :: intenv
     type (error_message) :: errmsg
   !
     call addTrace(errmsg,'createSolidt')
     this%nvar = 4
     this%vartyp = 'C'
     this%linear = .true.
     this%intenv => intenv
  end subroutine createSolidt
!----------------------------------------------------------------
!> \brief Deallocate
!
  subroutine deallocSolidt(this)
     class (solidt_ode_system) :: this
     if (associated(this%intenv)) nullify(this%intenv)  ! only nullify, no dellaoc
  end subroutine deallocSolidt
!-------------------------------------------------------------------
!> \brief Get system matrix
!! x:      location where rhs of ODE system is evaluated (in)
!! sysmat: matrix of ODE system (out)
!! ios:    error indicator (out)
!
  subroutine getSysmatComplexSolidt(this,x,sysmat,ios)
     class (solidt_ode_system) :: this
     double precision :: x
     double complex, dimension(:,:) :: sysmat
     integer :: ios
     double complex :: zrpc,z1,z2,zdlz1,zom,zomro
     double complex, dimension(7) :: zelcon
     double precision :: dll1,rr,rr2,ro
   !
     call evalComplexElasticConstantsSplineEarthmodel(x,this%intenv%nl,ro,zelcon)
   !
     zom = dcmplx(this%intenv%omre,this%intenv%omim)
     dll1 = this%intenv%dll1
     rr = 1.d0/x
     rr2 = rr*rr
     zrpc = 1.d0/zelcon(2)
     z1 = zelcon(1)-zelcon(3)**2*zrpc-zelcon(5)
     z2 = zelcon(3)*rr*zrpc
     zdlz1 = dll1*z1
     zomro = zom*zom*ro
     sysmat(1,1) = 2.d0*z2
     sysmat(2,1) = -zrpc
     sysmat(3,1) = -dll1*z2
     sysmat(4,1) = 0.d0
     sysmat(1,2) = zomro-4.d0*z1*rr2
     sysmat(2,2) = -sysmat(1,1)+2.d0*rr
     sysmat(3,2) = 2.d0*rr2*zdlz1
     sysmat(4,2) = -dll1*rr
     sysmat(1,3) = rr
     sysmat(2,3) = 0.d0
     sysmat(3,3) = -sysmat(1,3)
     sysmat(4,3) = -1.d0/zelcon(4)
     sysmat(1,4) = 2.d0*z1*rr2
     sysmat(2,4) = z2
     sysmat(3,4) = zomro - ( zdlz1 + (dll1-2.d0)*zelcon(5) )*rr2
     sysmat(4,4) = 3.d0*rr
     ios = 0
  end subroutine getSysmatComplexSolidt
!-------------------------------------------------------------------
!> \brief Get right hand side
!! x:   location where rhs of ODE system is evaluated (in)
!! rhs: right hand side of ODE system (out)
!! ios: error indicator (out)
!! This routine does nothing here.
!
  subroutine getRhsComplexSolidt(this,x,rhs,ios)
     class (solidt_ode_system) :: this
     double precision :: x
     double complex, dimension(:) :: rhs
     integer :: ios
     rhs = dcmplx(0.d0); ios = 0
  end subroutine getRhsComplexSolidt
!-------------------------------------------------------------------
!> \brief Get right hand side (real, not used)
!! x:   location where rhs of ODE system is evaluated (in)
!! rhs: right hand side of ODE system (out)
!! ios: error indicator (out)
!! This routine does nothing here.
!
  subroutine getRhsRealSolidt(this,x,rhs,ios)
     class (solidt_ode_system) :: this
     double precision :: x
     double precision, dimension(:) :: rhs
     integer :: ios
     rhs = 0.d0; ios = 0
  end subroutine getRhsRealSolidt
!-------------------------------------------------------------------
!> \brief Get system matrix (real, not used)
!! x:      location where rhs of ODE system is evaluated (in)
!! sysmat: matrix of ODE system (out)
!! ios:    error indicator (out)
!
  subroutine getSysmatRealSolidt(this,x,sysmat,ios)
     class (solidt_ode_system) :: this
     double precision :: x
     double precision, dimension(:,:) :: sysmat
     integer :: ios
     sysmat = 0.d0; ios = 0
  end subroutine getSysmatRealSolidt
!--------------------------------------------------------------------
!> \brief Get derivatives (complex, not used)
!! x:      location where rhs of ODE system is evaluated (in)
!! y:      current function vector (in)
!! dy:     derivative vector (out)
!! ios:    error indicator (out)
!
  subroutine derivativeComplexSolidt(this,x,y,dy,ios)
     class (solidt_ode_system) :: this
     double precision :: x
     double complex, dimension(:) :: y,dy
     integer :: ios
     dy = 0.d0; ios = 0
  end subroutine derivativeComplexSolidt
!--------------------------------------------------------------------
!> \brief Get derivatives (real, not used)
!! x:      location where rhs of ODE system is evaluated (in)
!! y:      current function vector (in)
!! dy:     derivative vector (out)
!! ios:    error indicator (out)
!
  subroutine derivativeRealSolidt(this,x,y,dy,ios)
     class (solidt_ode_system) :: this
     double precision :: x
     double precision, dimension(:) :: y,dy
     integer :: ios
     dy = 0.d0; ios = 0
  end subroutine derivativeRealSolidt
!
end module solidtOdeSystem

! ===========================================================================
!  Base class for handling systems of ordinary differential equations
! ===========================================================================
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
!  Defines an abstract type "baseOdeSystem" which specifies abstract interfaces
!  for procedures required to evaluate derivatives. The interfaces reflect linear
!  ODE systems where a sytem matrix can be defined and non-linear ones, where
!  this is not the case. The interfaces offer support for either real-valued 
!  or complex-valued solution vectors. To integrate a specific ODE, a concrete
!  type must be written which extends the base class and provides realizations 
!  of the abstract interfaces listed below. 
!-----------------------------------------------------------------------
module baseOdeSystem
    implicit none
    type, abstract :: base_ode_system
       integer :: nvar                         ! number of function variables
       character (len=1) :: vartyp             ! type of function variables
       logical :: linear                       ! ODE system is linear
    contains
        procedure (getSysmatComplexBaseOdeSystem), deferred :: getSysmatComplex
        procedure (getSysmatRealBaseOdeSystem), deferred :: getSysmatReal
        procedure (getRhsComplexBaseOdeSystem), deferred :: getRhsComplex
        procedure (getRhsRealBaseOdeSystem), deferred :: getRhsReal
        procedure (derivativeComplexBaseOdeSystem), deferred :: derivativeComplex
        procedure (derivativeRealBaseOdeSystem), deferred :: derivativeReal
        procedure :: isLinear
        procedure :: getNvar
        procedure :: getVartyp
    end type base_ode_system
!------------------------------------------------------------------
! Interface for obtaining complex system matrix
! x:      x-value where sysmat is computed
! sysmat: system matrix returned by procedure
! ios:    error indicator
!
    abstract interface
       subroutine getSysmatComplexBaseOdeSystem(this,x,sysmat,ios)
       import base_ode_system
       class (base_ode_system) :: this
       double precision :: x
       double complex, dimension(:,:) :: sysmat
       integer :: ios
       end subroutine getSysmatComplexBaseOdeSystem
    end interface
!------------------------------------------------------------------
! Interface for obtaining real system matrix
! x:      x-value where sysmat is computed
! sysmat: system matrix returned by procedure
! ios:    error indicator
!
    abstract interface
       subroutine getSysmatRealBaseOdeSystem(this,x,sysmat,ios)
       import base_ode_system
       class (base_ode_system) :: this
       double precision :: x
       double precision, dimension(:,:) :: sysmat
       integer :: ios
       end subroutine getSysmatRealBaseOdeSystem
    end interface
!------------------------------------------------------------------
! Interface for obtaining complex right hand side source term
! x:      x-value where right-hand side is computed
! rhs:    right-hand side vector returned by procedure
! ios:    error indicator
!
    abstract interface
       subroutine getRhsComplexBaseOdeSystem(this,x,rhs,ios)
       import base_ode_system
       class (base_ode_system) :: this
       double precision :: x
       double complex, dimension(:) :: rhs
       integer :: ios
       end subroutine getRhsComplexBaseOdeSystem
    end interface
!------------------------------------------------------------------
! Interface for obtaining real right hand side source term
! x:      x-value where right-hand side is computed
! rhs:    right-hand side vector returned by procedure
! ios:    error indicator
!
    abstract interface
       subroutine getRhsRealBaseOdeSystem(this,x,rhs,ios)
       import base_ode_system
       class (base_ode_system) :: this
       double precision :: x
       double precision, dimension(:) :: rhs
       integer :: ios
       end subroutine getRhsRealBaseOdeSystem
    end interface
!------------------------------------------------------------------
! Interface for obtaining complex derivative in non-linear case or linear
! case without use of system matrix
! x:      x-value where derivative is computed (input)
! y:      current function values needed for derivative (input) 
! dy:     derivative vector returned by procedure
! ios:    error indicator
!
    abstract interface
       subroutine derivativeComplexBaseOdeSystem(this,x,y,dy,ios)
       import base_ode_system
       class (base_ode_system) :: this
       double precision :: x
       double complex, dimension(:) :: y,dy
       integer :: ios
       end subroutine derivativeComplexBaseOdeSystem
    end interface
!------------------------------------------------------------------
! Interface for obtaining real derivative in non-linear case or linear
! case without use of system matrix
! x:      x-value where derivative is computed (input)
! y:      current function values needed for derivative (input) 
! dy:     derivative vector returned by procedure
! ios:    error indicator
!
    abstract interface
       subroutine derivativeRealBaseOdeSystem(this,x,y,dy,ios)
       import base_ode_system
       class (base_ode_system) :: this
       double precision :: x
       double precision, dimension(:) :: y,dy
       integer :: ios
       end subroutine derivativeRealBaseOdeSystem
    end interface
!
contains
!--------------------------------------------------------------------
!> \brief get number of function variables
!
    function getNvar(this) result(res)
    class (base_ode_system) :: this
    integer :: res
    res = this%nvar
    end function getNvar
!--------------------------------------------------------------------
!> \brief is ode system linear ?
!
    function isLinear(this) result(res)
    class (base_ode_system) :: this
    logical :: res
    res = this%linear
    end function isLinear
!--------------------------------------------------------------------
!> \brief get variable type
!
    function getVartyp(this) result(res)
    class (base_ode_system) :: this
    character(len=1) :: res
    res = this%vartyp
    end function getVartyp
end module baseOdeSystem 

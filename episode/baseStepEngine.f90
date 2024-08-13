! ===========================================================================
!  Base class for engine doing one single (accuracy-controlled) integration step
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
!-----------------------------------------------------------------------------
!--------------------------------------------------------------------------
!  Defines an abstract base class for an integration method (such as Runge-Kutta
!  or Bulirsch-Stoer) which does one single, potentially accuracy-controlled,
!  integration step. Several of such steps may be required to straddle one sub-interval
!  (e.g. from node to node). User-provided methods should extend this base class and
!  implement procedures according to the abstract interfaces defined below.
!  Since the integration method needs to compute derivatives, this base class
!  "has" a polymorphic member specifying the ODE to be solved.
!-------------------------------------------------------------
module baseStepEngine
    use errorMessage
    use baseOdeSystem
    implicit none
!-------------------------------------------------------------------------------------------------------------
    type, abstract :: base_step_engine
       class (base_ode_system), pointer :: sode                                     ! polymorphic object representing SODE to be solved
       double precision :: eps                                                      ! desired accuracy
       integer :: secrecy                                                           ! secrecy level below which debugging output is produced
       integer :: nvar                                                              ! number of solution variables
       type (error_message) :: errmsg                                               ! object owned error message
       logical :: success                                                           ! step was successful
    contains
       procedure (createBaseStepEngine), deferred :: create                         ! create integration engine
       procedure (deallocBaseStepEngine), deferred :: dealloc                       ! dealloc integration engine
       procedure (doFullStepBaseStepEngine), deferred :: doFullStep                 ! carry out one, accuracy-controlled step
       procedure (getSolutionBaseStepEngine), deferred :: getSolution               ! get solution value at end of step
       procedure (adjustStepsizeBaseStepEngine), deferred :: adjustStepsize         ! optionally modify step size after a step
       procedure :: printErrmsg
       procedure :: getErrlevel
       procedure :: successful
    end type base_step_engine
!----------------------------------------------------------------
!  Interface for construction of integration engine object
!  sode:    Polymorphic object representing sytem of ODEs (in)
!  eps:     Desired accuracy of integration (in)
!  secrecy: Level below which debugging output is produced (in)
!
    abstract interface
       subroutine createBaseStepEngine(this,sode,eps,secrecy)
       import base_step_engine
       import base_ode_system
       class (base_step_engine) :: this
       class (base_ode_system), target :: sode
       double precision :: eps
       integer :: secrecy
       end subroutine createBaseStepEngine
    end interface
!--------------------------------------------------------------
!  Interface for doing one accuracy-controlled integration step
!  h:       trial step size (in)
!  xstart:  starting point (in)
!  yr:      real part of initial solution value (in)
!  yi:      imaginary part of initial solution value (in)
!
    abstract interface
       subroutine doFullStepBaseStepEngine(this,h,xstart,yr,yi)
       import base_step_engine
       import base_ode_system
       class (base_step_engine) :: this
       double precision :: h,xstart
       double precision, dimension(:) :: yr,yi
       end subroutine doFullStepBaseStepEngine
    end interface
!--------------------------------------------------------------
!  Interface for deallocation
!
    abstract interface
       subroutine deallocBaseStepEngine(this)
       import base_step_engine
       class (base_step_engine) :: this
       end subroutine deallocBaseStepEngine
    end interface
!--------------------------------------------------------------
!  Interface for recovering solution vector
!  yr:      real part of solution value (out)
!  yi:      imaginary part of solution value (out)
!
    abstract interface
       subroutine getSolutionBaseStepEngine(this,yr,yi)
       import base_step_engine
       class (base_step_engine) :: this
       double precision, dimension(:) :: yr,yi
       end subroutine getSolutionBaseStepEngine
    end interface
!-------------------------------------------------------------
!  Interface for adjusting step size
!  h:       new step size (out)
!  success: tells whether step was successful or not (out)
!
    abstract interface
       subroutine adjustStepsizeBaseStepEngine(this,h,success)
       import base_step_engine
       class (base_step_engine) :: this
       double precision :: h
       logical :: success
       end subroutine adjustStepsizeBaseStepEngine
    end interface
!-------------------------------------------------------------
    interface operator (.errlevel.); module procedure getErrlevel; end interface
    interface operator (.success.); module procedure successful; end interface
!
contains
!-------------------------------------------------------------------------------
!> \brief print error message
!
    subroutine printErrmsg(this)
    class (base_step_engine) :: this
    call print(this%errmsg)
    end subroutine printErrmsg
!-------------------------------------------------------------------------------
!> \brief print error message
!
    function getErrlevel(this) result(res)
    class (base_step_engine), intent(in) :: this
    integer :: res
    res = .level.(this%errmsg)
    end function getErrlevel
!---------------------------------------------------------------
!> \brief Step successful ?
!
    function successful(this) result(res)
    class (base_step_engine), intent(in) :: this
    logical :: res
    res = this%success
    end function successful
end module

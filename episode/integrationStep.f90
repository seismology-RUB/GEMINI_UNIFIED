! ===========================================================================
!  Module organizing integration through sub-intervals
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
!-------------------------------------------------------------
!  Organizes the integration through one sub-interval which may
!  be split into several individual steps depending on accuracy
!  requirements. This module acts a kond of driver for the 
!  integration engine. Since the solution is only stored at the
!  nodes, there is no need for intermediate storage of the solution
!  here. This module makes use of a polymorphic object for the
!  ODE system and for the integration engine.
!--------------------------------------------------------------
module integrationStep
    use baseStepEngine
    use baseOdeSystem
    implicit none
    integer, parameter :: secrecy_integration_step = 1                              ! print screen output if secrecy <= this value
    double precision, parameter :: minsize_integration_step = 1.d-12                ! minimal step size
    integer, parameter :: maxnum_integration_step = 200000                          ! maximum number of steps for interval
    type integration_step
       class (base_ode_system), pointer :: sode                                     ! polymorphic object representing system of ODEs
       class (base_step_engine), pointer :: steng                                   ! polymorphic object representing stepping engine
       double precision :: eps                                                      ! desired accuracy
       integer :: secrecy                                                           ! secrecy level below which debugging info is produced
       type (error_message) :: errmsg                                               ! object owned error message
       double precision :: hdid                                                     ! real step size
       double precision :: htry                                                     ! trial step size
       integer :: nvar                                                              ! number of variables
       double precision, dimension(:), allocatable :: yr                            ! real inital function values, also end values
       double precision, dimension(:), allocatable :: yi                            ! imaginary initial function values, also end values
    end type integration_step
    interface dealloc; module procedure deallocIntegrationStep; end interface
    interface operator (.errlevel.); module procedure getErrlevelIntegrationStep; end interface
!
contains
!-----------------------------------------------------------
!> \brief create real integration step object
!  htry:     trial step size
!  sode:     object representing ODEs
!  steng:    object representing integration engine (stepping engine)
!  eps:      desired accuracy
!  secrecy:  level below which debugging info is produced
!
    subroutine createIntegrationStep(this,htry,sode,steng,eps,secrecy)
    type (integration_step) :: this
    class (base_ode_system), target :: sode
    class (base_step_engine), target :: steng
    double precision :: eps,htry
    integer :: secrecy
    character(len=21) :: myname = 'createIntegrationStep'
!
    call new(this%errmsg,myname)
    this%htry = htry
    this%hdid = htry
    this%sode => sode
    this%steng => steng
    this%eps = eps
    this%secrecy = secrecy
    this%nvar = this%sode%getNvar()
    allocate(this%yr(this%nvar),this%yi(this%nvar))
    end subroutine createIntegrationStep
!----------------------------------------------------------------
!> \brief Deallocate 
!
    subroutine deallocIntegrationStep(this)
    type (integration_step) :: this
    if (allocated(this%yr)) deallocate(this%yr)
    if (allocated(this%yi)) deallocate(this%yi)
    if (associated(this%sode)) nullify(this%sode)
    if (associated(this%steng)) nullify(this%steng)
    if (this%secrecy <= secrecy_integration_step) call print(this%errmsg)
    call dealloc(this%errmsg)
    end subroutine deallocIntegrationStep
!----------------------------------------------------------------------
!> \brief do integration sub-interval by accumulating integration steps
!  x1:    begin of sub-interval (in)
!  x2:    end of sub-interval (in)
!  yr:    real part of initial solution vector (in)
!  yi:    imaginary part of initial solution vector (in) 
!
    subroutine doIntegrationStep(this,x1,x2,yr,yi)
    type (integration_step) :: this
    character(len=17) :: myname = 'doIntegrationStep'
    double precision :: x1,x2
    double precision, dimension(:) :: yr,yi
    double precision, dimension(:), allocatable :: yrs,yis
    double precision :: x,h
    integer :: nstep
!
!  check if x2=x1, because integration sub-interval may start 
!  (or end) with a node or where we want the solution
!  we do not need to call step engine but need the solution
!
    if (abs(x2-x1) < minsize_integration_step) then                       
       this%yr = yr; this%yi = yi                      
       return
    endif
!
!  Initialize step engine object
!  use current this%hdid as engine step size
!
    call this%steng%create(this%sode,this%eps,this%secrecy)
    if (.errlevel.(this%steng) == 2) then
       call add(this%errmsg,2,'Problems creating step engine',myname)
       call this%steng%printErrmsg(); return
    endif
!
!  Initialize begin of interval, step length
!  and copy initial solution to temporary array
!
    x = x1
    h = sign(this%hdid,x2-x1)
    allocate(yrs(this%nvar),yis(this%nvar))
    yrs = yr; yis = yi
!
!  Loop over sequence of steps done by engine
!
    nstep = 0
    do while(nstep < maxnum_integration_step)
       if((x+h-x2)*(x+h-x1).gt.0.d0) then
          h = x2-x                                                ! avoid overshoot
       endif
       call this%steng%doFullStep(h,x,yrs,yis)                    ! do integration step
       if (.errlevel.(this%steng) == 2) then
          call add(this%errmsg,2,'Problems performing integration step',myname)
          call this%steng%printErrmsg(); return
       endif
       nstep = nstep+1
!
!  Engine step was successful, copy solution to ystart, update x,
!  and go on with next step if not finished
!
       if (.success.(this%steng)) then                                       ! step was successful
          x = x+h                                                            ! update x
          call this%steng%getSolution(yrs,yis)                               ! copy updated solution into temporary array
          if (this%secrecy <= secrecy_integration_step) then
             call printStengResultIntegrationStep(this,h,x,yrs,yis)
          endif
          if (abs(x-x2) < minsize_integration_step) then                     ! sub-interval end was reached
             this%yr = yrs; this%yi = yis                                    ! keep solution at end of sub-interval
             goto 1                                                          ! keep hdid used in previous steps
          else                                                               ! next engine step
             call this%steng%adjustStepsize(h,.true.)                        ! sub-interval end not yet reached
             this%hdid = h                                                   ! allow for dynamic adjustment of step size
          endif
!
!  Engine step failed, htry was too large,
!  drastically reduce step size and try again using original values
!  for x and yrs (no update, keep steng object)
!  Before retrying check whether minimum step size was reached
!
       else
          if (this%secrecy <= secrecy_integration_step) then 
             print *,'doIntegrationStep: engine step not successful, decrease h'
             print *,'doIntegrationStep: x1,x2,htry,hdid: ',x1,x2,this%htry,this%hdid
             print *,'Current h: ',h
          endif
          call this%steng%adjustStepsize(h,.false.)
          if (this%secrecy <= secrecy_integration_step) print *,'Next h: ',h
          if (abs(h) < minsize_integration_step) then
             call add(this%errmsg,2,'Steps size smaller than minimum',myname)
             goto 1
          endif
       endif
    enddo         ! end of engine step loop
!
!  clean up
!
    call add(this%errmsg,2,'Maximum number of steps reached',myname)
1   call this%steng%dealloc()
    if (allocated(yrs)) deallocate(yrs)
    if (allocated(yis)) deallocate(yis)
    end subroutine doIntegrationStep
!----------------------------------------------------------------------
!> \brief print result of successful integration step
!
    subroutine printStengResultIntegrationStep(this,h,x,yr,yi)
    type (integration_step) :: this
    double precision :: h,x
    double precision, dimension(:) :: yr,yi
!
    print *,'Result after successful Bulirsch step'
    print *,'h = ',h
    print *,'x = ',x
    if (this%sode%getVartyp() == 'C') then
       print *,'real(ysol) = ',yr
       print *,'imag(ysol) = ',yi
    else if (this%sode%getVartyp() == 'R') then
       print *,'real(ysol) = ',yr
    endif
    end subroutine printStengResultIntegrationStep
!-------------------------------------------------------------------------------
!> \brief print error message
!
    subroutine printErrmsgIntegrationStep(this)
    type (integration_step) :: this
    call print(this%errmsg)
    end subroutine printErrmsgIntegrationStep
!-------------------------------------------------------------------------------
!> \brief print error message
!
    function getErrlevelIntegrationStep(this) result(res)
    type (integration_step), intent(in) :: this
    integer :: res
    res = .level.(this%errmsg)
    end function getErrlevelIntegrationStep
!---------------------------------------------------------------------
!> \brief Get solution at end of step (complex)
!
    subroutine getSolutionIntegrationStep(this,yr,yi)
    type (integration_step) :: this
    double precision, dimension(:) :: yr,yi
    yr = this%yr; yi = this%yi
    end subroutine getSolutionIntegrationStep
!
end module integrationStep

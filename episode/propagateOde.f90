! ====================================================================
!  General propagator for systems of ordinary differential equations
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
!--------------------------------------------------------------------------------
!  To use this propagator, you need to provide concrete classes for
!  - the ODE-system you want to solve according to specifications in baseOdeSystem.f90,
!  - the integration environment according to specifications in baseIntegrationEnvironment.f90,
!  - the stepping engine for solving the ODE according to specs in baseStepEngine.f90.
!
!  The propagator drives the integration over intervals (so-called "layers") within 
!  which the derivatives of the solution are continuous. It handles stopping and 
!  continuation at interval boundaries as well as at intermediate nodes. The solution 
!  at intermediate nodes can be stored. For doing the integration, the propagator uses 
!  the integrationStep-module which, for each subinterval, organizes a sequence of 
!  smaller steps that are in turn carried out by the stepping engine. The propagator 
!  handles complex and real valued ODEs and can work with any stepping engine that 
!  provides interfaces defined in baseStepEngine.f90.
!-------------------------------------------------------------------------------------
module propagateOde
    use integrationStep
    use baseOdeSystem
    use baseIntegrationEnvironment
    use baseStepEngine
    use errorMessage
    implicit none
    interface operator (.errlevel.); module procedure getErrlevelPropagateOde; end interface
    interface getSolutionPropagateOde
       module procedure getComplexSolutionPropagateOde
       module procedure getRealSolutionPropagateOde
    end interface
    interface getSolutionAtNodesPropagateOde
       module procedure getComplexSolutionAtNodesPropagateOde
       module procedure getRealSolutionAtNodesPropagateOde
    end interface
    interface dealloc; module procedure deallocPropagateOde; end interface
    integer, parameter :: secrecy_propagate_ode = 2                   ! print screen output if secrecy <= this value
    double precision, parameter :: maxy_propagate_ode = 1.d10         ! limit for norm of y triggering normalization    
    double precision, parameter :: lg_maxy_propagate_ode = 10.d0      ! log10-limit for norm of y triggering normalization    
    type :: propagate_ode
       private
       class (base_ode_system), pointer :: sode                       ! Polymorphic ODE system object
       class (base_integration_environment), pointer :: intenv        ! Polymorphic integration environment
       class (base_step_engine), pointer :: steng                     ! Polymorphic single step integration engine
       double precision :: eps                                        ! desired error bound
       double precision :: lgnorm                                     ! current log-normalization factor
       integer :: secrecy                                             ! secrecy level
       integer :: nvar                                                ! number of function values
       integer :: jl,jr                                               ! global indices of first and last storage node
       integer :: jld,jrd                                             ! global indices of first and last normalization disco
       logical :: do_norm                                             ! do normalization of solution
       double precision, dimension(:), allocatable :: yr,yi           ! holds intermediate function values
       double precision, dimension(:,:), allocatable :: ynr,yni       ! holds function values at nodes (ynr(i,node))
       double precision, dimension(:), allocatable :: lgynorm         ! lg-normalization factor at nodes
       type (error_message) :: errmsg                                 ! object owned error message
    end type
!
contains
!------------------------------------------------------------
!> \brief Create a propagator object
!!  sode:    concrete type extending baseOdeSystem
!!  intenv:  concrete type extending baseIntegrationEnvironment
!!  steng:   concrete type extending baseStepEngine
!!  eps:     desired relative solution accuracy
!!  secrecy: level up to which debugging output is produced
!!  do_norm: do normalization flag
!!  lgnorm:  initial value of log normalization factor
!
    subroutine createPropagateOde(this,sode,intenv,steng,eps,secrecy,do_norm,lgnorm)
    type (propagate_ode) :: this
    class (base_ode_system), target :: sode
    class (base_integration_environment), target :: intenv
    class (base_step_engine), target :: steng
    double precision :: eps,lgnorm
    integer :: secrecy
    logical :: do_norm
    character(len=18) :: myname = 'createPropagateOde'
!
    call new(this%errmsg,myname)
    this%sode => sode
    this%intenv => intenv
    this%steng => steng
    this%eps = eps
    this%lgnorm = lgnorm
    this%do_norm = do_norm
    this%secrecy = secrecy
    this%nvar = sode%getNvar()
    allocate(this%yr(this%nvar),this%yi(this%nvar))
    end subroutine createPropagateOde
!----------------------------------------------------------------
!> \brief Deallocate 
!
    subroutine deallocPropagateOde(this)
    type (propagate_ode) :: this
    if (allocated(this%yr)) deallocate(this%yr)
    if (allocated(this%yi)) deallocate(this%yi)
    if (allocated(this%ynr)) deallocate(this%ynr)
    if (allocated(this%yni)) deallocate(this%yni)
    if (allocated(this%lgynorm)) deallocate(this%lgynorm)
    if (associated(this%sode)) nullify(this%sode)
    if (associated(this%intenv)) nullify(this%intenv)
    if (associated(this%steng)) nullify(this%steng)
    if (this%secrecy <= secrecy_propagate_ode) call print(this%errmsg)
    call dealloc(this%errmsg)
    end subroutine deallocPropagateOde
!----------------------------------------------------------------
!> \brief Propagate ODE system over specific large interval
!!  Use ODE system, integration environment and stepping engine
!!  as specified on construction of propagateOde object.
!!  x1:      start of interval (input)
!!  x2:      end of interval (input)
!!  htry:    inital step size for integration (input)
!!  yrstart: real part of initial values of function vector (input)
!!  yistart: imaginary part of initial values of function vector (input)
!
    subroutine doPropagateOde(this,x1,x2,htry,yrstart,yistart)
    type (propagate_ode) :: this
    double precision, intent(in) :: x1,x2,htry
    double precision, dimension(:), intent(in) :: yrstart,yistart
    character(len=14) :: myname = 'doPropagateOde'
    double precision, dimension(:), allocatable :: xd,xn       ! location of discontinuities and nodes
    integer :: jd,jnod,nd,nnod,idir,jl,jr
    double precision, dimension(3) :: steps                    ! 3 possible lengths of next subinterval
    double precision :: xa,xe                                  ! intermediate integration intervals
    double precision :: ynorm                                  ! solution norm
    type (integration_step) :: is                              ! stepping object for intermediate intervals
!
    call addTrace(this%errmsg,myname)
    this%yr = yrstart; this%yi = yistart                       ! copy initial value to final solution value
!
!  deallocate solution at nodes in case previous run on
!  different integration range was done
!
    if (allocated(this%ynr)) deallocate(this%ynr)
    if (allocated(this%yni)) deallocate(this%yni)
    if (allocated(this%lgynorm)) deallocate(this%lgynorm)
!
!  nodes and discontinuities within integration range, routines return
!  descending order if x2 < x1, else ascending order
!  ascending  case: ---x1---xd(1)---------------xd(n)----x2------>
!  descending case: ---x2---xd(n)---------------xd(1)----x1 ----->
!
    nnod = 0; nd = 0
    call this%intenv%getNodes(x1,x2,jl,jr,xn,this%errmsg)
    this%jl = jl; this%jr = jr                                       ! assign values in any case, even if no nodes in interval, getNodes takes care
    if (allocated(xn)) then
       nnod = size(xn)
       allocate(this%ynr(this%nvar,nnod),this%yni(this%nvar,nnod))
       allocate(this%lgynorm(nnod))
    endif
    call this%intenv%getDiscontinuities(x1,x2,jl,jr,xd,this%errmsg)  ! Solution not stored at discos
    this%jld = jl; this%jrd = jr                                     ! assign values in any case, even in no disco in interval
    if (allocated(xd)) then
       nd = size(xd)
    endif
    idir = int(sign(1.d0,x2-x1))                                     ! Integration direction
!
!  Check whether discontinuities respectively nodes lie within integration range
!
    if (nd > 0) then
       if ((x1-xd(1))*idir > 0 .or. (xd(nd)-x2)*idir > 0) then
          call add(this%errmsg,2,'Discontinuities not entirely within integration range',myname)
          return
       endif
    endif
    if (nnod > 0) then
       if ((x1-xn(1))*idir > 0 .or. (xn(nnod)-x2)*idir > 0) then
          call add(this%errmsg,2,'Nodes not entirely within integration range',myname)
          return
       endif
    endif
!
!  Check if integration step size is zero. IF yes, and we start on a node, we need
!  to store the starting values as result.
!
    if (dabs(x2-x1) < minsize_integration_step .and. nnod == 0) then
       call add(this%errmsg,0,'Integration interval is zero and there is no node, do nothing',myname)
       return
    else if (dabs(x2-x1) < minsize_integration_step .and. nnod > 1) then
       call add(this%errmsg,2,'Integration interval is zero, but there are more than two nodes???',myname)
       return
    else if (dabs(x2-x1) < minsize_integration_step .and. nnod == 1) then
       this%ynr(:,1) = this%yr(:); this%yni(:,1) = this%yi(:)
       jnod = 2; jd = 1                                           ! Indices of next node and disco to be reached
       goto 1                                                     ! leave routine after index check
    endif       
!
!  Initialize integration step object, used for the all following subintervals
!  between discos or nodes or discos and nodes
!  -----N----D---N---N---N---N---D---N---N---ND----
!
    call createIntegrationStep(is,htry,this%sode,this%steng,this%eps,this%secrecy)
!
!  Initialize indices of nodes and discos to be reached next
!  Set start of subinterval
!
    jd = 1; jnod = 1
    xa = x1
!
!  Start the loop over integration sub-intervals
!
    do while (dabs(x2-xa) > minsize_integration_step)
!
!  Find out whether the next stop is a discontinuity or a node or 
!  the end of the interval. It could also happen that a disco coincides
!  with a node. We do not take care of a disco coinciding 
!  with either x1 or x2 because discos are assumed to lie
!  within the open integration interval. This must be ensured by the
!  integration environment. Nodes may lie on the interval ends.
!  There are three possible intervall lengths: 
!  (1) |xn(jnod)-xa|; (2) |xd(jd)-xa|; (3) |x2-xa|.
!  We must find the shortest one. 
!
       steps(3) = dabs(x2-xa)                 ! biggest possible step
       steps(1) = 2.*dabs(x2-x1)              ! initialize with value greater than steps(3)
       steps(2) = 3.*dabs(x2-x1)              ! initialize with value greater than steps(3) and steps(1)
       if (jnod <= nnod) then                 ! not all nodes reached yet
          steps(1) = dabs(xn(jnod)-xa)
       endif
       if (jd <= nd) then                     ! not all disco reached yet
          steps(2) = dabs(xd(jd)-xa)
       endif
       xe = xa+idir*minval(steps)             ! take the smallest of the 3
!
!  integrate over sub-interval using current solution values
!
       call doIntegrationStep(is,xa,xe,this%yr,this%yi)
       if (.errlevel.is == 2) then
          call printErrmsgIntegrationStep(is)
          call add(this%errmsg,2,'Problems doing integration step',myname)
          return
       endif
!
!  update current solution value
!
       call getSolutionIntegrationStep(is,this%yr,this%yi)
       if (this%secrecy <= secrecy_propagate_ode) call printResultPropagateOde(this,xa,xe)
!
!  check for normalization
!
       if (this%do_norm) then                                        ! check for normalization
          ynorm = sum(abs(this%yr)+abs(this%yi))/this%nvar           ! avoid squaring
          if (ynorm > maxy_propagate_ode) then
             this%yr = this%yr/ynorm                                 ! renormalize solution
             this%yi = this%yi/ynorm
             this%lgnorm = this%lgnorm+log10(ynorm)                  ! update current normalization factor 
             if (this%secrecy <= secrecy_propagate_ode) then
                print *,myname,': after normalization'
                call printResultPropagateOde(this,xa,xe)
             endif
          endif
       endif
!
!  What was at the end of current interval?
!  if it was a node, store result and increase jnod
!  if it was a disco, increase jd
!  do both if it was both
!
       select case (minloc(steps,1))
       case (1)                                                         ! it was a node
          this%ynr(:,jnod) = this%yr(:); this%yni(:,jnod) = this%yi(:)  ! store current solution
          this%lgynorm(jnod) = this%lgnorm                              ! store normalization factor
          jnod = jnod+1                                                 ! next node to be reached
          if (dabs(steps(1)-steps(2)) < minsize_integration_step) then  ! and also a disco
             call this%intenv%nextInterval(idir)                        ! we go to the next layer
             jd = jd+1                                                  ! next disco to be reached
             if (this%secrecy <= secrecy_propagate_ode) then
                print *,'Interval end was a node and a disco: ',jnod-1,' Disco ID: ',jd-1
                print *,'Now move to layer: ',this%intenv%nl
             endif
          else
             if (this%secrecy <= secrecy_propagate_ode) then
                print *,'Interval end was a node: ',jnod-1,' Disco ID: ',jd-1
                print *,'Stay in layer: ',this%intenv%nl
             endif
          endif
       case (2)                                                            ! it was a discontinuity
          jd = jd+1                                                        ! next disco to be reached
          call this%intenv%nextInterval(idir)                              ! we go to next layer
          if (dabs(steps(1)-steps(2)) < minsize_integration_step) then     ! and also a node
             this%ynr(:,jnod) = this%yr(:); this%yni(:,jnod) = this%yi(:)  ! store current solution
             this%lgynorm(jnod) = this%lgnorm                              ! store normalization factor
             jnod = jnod+1                                                 ! next mode to be reached
             if (this%secrecy <= secrecy_propagate_ode) then
                print *,'Interval end was a disco and a node: ',jnod-1,' Disco ID: ',jd-1
                print *,'Now move to layer: ',this%intenv%nl
             endif
          else
             if (this%secrecy <= secrecy_propagate_ode) then
                print *,'Interval end was a disco: ',jd-1,' Node ID: ',jnod-1
                print *,'Now move to layer: ',this%intenv%nl
             endif
          endif
       end select
       xa = xe                                                          ! end of sub_interval becomes new start
    enddo                                                               ! next sub-interval
!
!  propagator at end of integration range
!  deallocate integration step
!  check if all nodes and discos have been passed
!
    call dealloc(is)
 1  if (jnod-1 /= nnod) then
       call add(this%errmsg,2,'Not all nodes have been passed',myname)
    endif
    if (jd-1 /= nd) then
       call add(this%errmsg,2,'Not all discos have been passed',myname)
    endif
!
!  print out values at nodes
!
    if (this%secrecy <= secrecy_propagate_ode) then
       if (nnod > 0) call printResultsNodesPropagateOde(this,xn)
    endif
    end subroutine doPropagateOde
!-------------------------------------------------------------------------------
!> \brief Print results for current interval
!
    subroutine printResultPropagateOde(this,x1,x2)
    type (propagate_ode) :: this
    double precision :: x1,x2
!
    print *,'Results of propagate step'
    print *,'x1 = ',x1
    print *,'x2 = ',x2
    if (this%sode%getVartyp() == 'C') then
       print *,'real(ysol) = ',this%yr
       print *,'imag(ysol) = ',this%yi
    else if (this%sode%getVartyp() == 'R') then
       print *,'ysol = ',this%yr
    endif
    end subroutine printResultPropagateOde
!-------------------------------------------------------------------------------
!> \brief Print results at nodes for entire integration range
!
    subroutine printResultsNodesPropagateOde(this,xn)
    type (propagate_ode) :: this
    double precision, dimension(:) :: xn
    integer :: j
!
    print *,'Results at nodes'
    do j = 1,size(xn)
       print *,'Index: ',j
       print *,'xn = : ',xn(j)
       if (this%sode%getVartyp() == 'C') then
          print *,'real(ysol) = ',this%ynr(:,j)
          print *,'imag(ysol) = ',this%yni(:,j)
       else if (this%sode%getVartyp() == 'R') then
          print *,'ysol = ',this%ynr(:,j)
       endif
    enddo
    end subroutine printResultsNodesPropagateOde
!-------------------------------------------------------------------------------
!> \brief print error message
!
    subroutine printErrmsgPropagateOde(this)
    type (propagate_ode) :: this
    call print(this%errmsg)
    end subroutine printErrmsgPropagateOde
!-------------------------------------------------------------------------------
!> \brief Get error level of propagateOde error message object
!
    function getErrlevelPropagateOde(this) result(res)
    type (propagate_ode), intent(in) :: this
    integer :: res
    res = .level.(this%errmsg)
    end function getErrlevelPropagateOde
!---------------------------------------------------------------------
!> \brief Get solution at end of integration range (complex)
!
    subroutine getComplexSolutionPropagateOde(this,y)
    type (propagate_ode) :: this
    double complex, dimension(:) :: y
    y = dcmplx(this%yr,this%yi)
    end subroutine getComplexSolutionPropagateOde
!---------------------------------------------------------------------
!> \brief Get solution at end of integration range (real)
!
    subroutine getRealSolutionPropagateOde(this,y)
    type (propagate_ode) :: this
    double precision, dimension(:) :: y
    y = this%yr
    end subroutine getRealSolutionPropagateOde
!---------------------------------------------------------------------
!> \brief Get solution at nodes (double)
!
    subroutine getRealSolutionAtNodesPropagateOde(this,y)
    type (propagate_ode) :: this
    double precision, dimension(:,:) :: y
    y = this%ynr
    end subroutine getRealSolutionAtNodesPropagateOde
!---------------------------------------------------------------------
!> \brief Get solution at nodes (double)
!   Do not correct solution for normalizations done during integration    
!
    subroutine getComplexSolutionAtNodesPropagateOde(this,y)
    type (propagate_ode) :: this
    double complex, dimension(:,:) :: y
    y = dcmplx(this%ynr,this%yni)
    end subroutine getComplexSolutionAtNodesPropagateOde
!---------------------------------------------------------------------
!> \brief Get original indices of leftmost and rightmost node
!
    subroutine getIndexLimitsNodesPropagateOde(this,jl,jr)
    type (propagate_ode) :: this
    integer :: jl,jr
    jl = this%jl; jr = this%jr
    end subroutine getIndexLimitsNodesPropagateOde
!---------------------------------------------------------------------
!> \brief Get log normalization factors at nodes (real)
!
    subroutine getLogNormalizationAtNodesPropagateOde(this,lgy)
    type (propagate_ode) :: this
    double precision, dimension(:) :: lgy
    lgy = this%lgynorm
    end subroutine getLogNormalizationAtNodesPropagateOde
!---------------------------------------------------------------------
!  Get current value of normalization factor
!
    subroutine getLogNormalizationPropagateOde(this,lgend)
    type (propagate_ode) :: this
    double precision :: lgend
    lgend = this%lgnorm
    end subroutine getLogNormalizationPropagateOde
!
end module propagateOde

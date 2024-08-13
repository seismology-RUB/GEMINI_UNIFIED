! ====================================================================================
!  Computes solution of negative transposed spheroidal system using EPISODE propagator
! ====================================================================================
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
!--------------------------------------------------------------
!! Does an integration of the negative transposed solid ODE to later obtain
!! Green functions by combining this solution with the minor matrix
!! Integration is either done up- or downwards from a potential source node
!! to the nodes. Depending on the initial conditions at the source
!! node, the Green function satisfies a jump condition there.
!------------------------------------------------------------------
module spheroidalAdjoint
    use propagateOde
    use solidtOdeSystem
    use geminiIntegrationEnvironment
    use complexBulirschStep
    use splineEarthmodelCoefficients    ! needs to evaluate material parameters
    use errorMessage
    implicit none
    interface dealloc; module procedure deallocSpheroidalAdjoint; end interface
    integer, parameter :: secrecy_spheroidal_adjoint = 3             ! print info if secrecy <= this value
    type spheroidal_adjoint
       type (gemini_integration_environment) :: intenv                   ! gemini as integration environment
       double complex, dimension(:,:), allocatable :: adj                ! adjoint solution at nodes
    end type
!
contains
!-----------------------------------------------------------------------------------------
!>  \brief Compute spheroidal adjoint solutions
!!  ystart:     initial values of solution vector (overwritten with values at end on retun)
!!  ra:         start of integration interval
!!  re:         end of integration interval
!!  gem_intenv: gemini integration environment object (in)
!!  eps:        desired accuracy (in)
!!  secrecy:    anti-debuggung level (in)
!
  subroutine computeSpheroidalAdjoint(this,ystart,ra,re,gem_intenv,eps,secrecy,errmsg)
     type (spheroidal_adjoint) :: this
     type (gemini_integration_environment) :: gem_intenv
     double precision :: eps
     integer :: secrecy
     type (error_message) :: errmsg
     type (solidt_ode_system) :: sode
     type (complex_bulirsch_step) :: cbs                              ! complex bs-step as step engine
     type (propagate_ode) :: prop
     double complex, dimension(:) :: ystart
     double complex, dimension(:,:) , allocatable :: zyn
     double precision :: ra,re,v,htry
     character(len=24) :: myname = 'computeSpheroidalAdjoint'
     integer :: nl,top,j,jl,jr,nnod
   !
     call addTrace(errmsg,myname)
   !
   !  get number of nodes and do some initial checks
   !
     nnod = .nnod.(gem_intenv%exnod)
     this%intenv = gem_intenv            ! copy sufficient here, this%intenv will be updated
   !
   !  set up ODE system
   !  this%sode must point to this%intenv because the latter  
   !  (member nl) will be updated during integration (and not gem_intenv)
   !
     call sode%createSolidt(this%intenv,errmsg)
     if (.level.(errmsg) == 2) return
   !
   !  allocate space and zero
   !
     allocate(this%adj(4,nnod))
     this%adj = dcmplx(0.d0,0.d0)
   !
   !  Establish propagator object for integration
   !  The already created solidtOdeSystem object is passed
   !  to the polymorphic base ODE system of propagateOde.
   !  Also, the geminiIntegrationEnvironment object is passed to the
   !  polymorphic integration environment of propagateOde.
   !  The up to now only declared Bulirsch-Store integration object is passed
   !  to the polymorphic base integration engine of propagateOde. It will
   !  be created when running doIntegrationStep.
   !
     call createPropagateOde(prop,sode,this%intenv,cbs,eps,secrecy,.false.,0.d0)
     if (.errlevel.prop == 2) then
        call printErrmsgPropagateOde(prop)
        call dealloc(prop)
        call add(errmsg,2,'Problems with spheroidal adjoint creating propagator',myname)
        return
     endif
   !
   !  integration from ra to re
   !
     if (secrecy <= secrecy_spheroidal_adjoint) then
        print *,'spheroidalAdjoint: integration from ',ra,' to ',re
     endif
     call getLayerIndexNodeEarthmodel(this%intenv%nem,ra,nl,top)
     if (top == 1 .and. ra < re) nl = min(nl+1,.nlay.(this%intenv%nem))
     call this%intenv%setInterval(nl)
   !
     v = getSVelocitySplineEarthmodel(nl,ra)   ! estimate for step size
     if (v < epsilon(1.d0)) v = getPVelocitySplineEarthmodel(nl,ra)
     htry = 0.25*v*mc_two_pid/this%intenv%omre
   !
     call doPropagateOde(prop,ra,re,htry,real(ystart),imag(ystart))
     if (.errlevel.prop == 2) then
        call printErrmsgPropagateOde(prop)
        call add(errmsg,2,'Problems with doing spheroidal adjoint propagation',myname)
        call dealloc(prop); return
     endif
   !
   !  extract adjoint solution from propagator
   !
     call getIndexLimitsNodesPropagateOde(prop,jl,jr)
     if (jr-jl < 0) then
        call add(errmsg,2,'Either problems with node indices or no nodes in integration range',myname)
        return
     endif
     allocate(zyn(4,jr-jl+1))
     call getComplexSolutionAtNodesPropagateOde(prop,zyn)
     if (re > ra) then
        do j = jl,jr
           this%adj(:,j) = zyn(:,j-jl+1)
        enddo
     else
        do j = jl,jr
           this%adj(1:4,j) = zyn(:,jr-j+1)                  ! inverse order in zyn
        enddo
     endif
!
!  get solution at end of interval
!
     call getComplexSolutionPropagateOde(prop,ystart)
!
!  Debug output
!
     if (secrecy <= secrecy_spheroidal_adjoint) then
        print *,'Adjoint: j,this%adj:'
        do j = jl,jr
           write(6,'(i6,8e15.3)') j,this%adj(1:4,j)
        enddo
     endif
     deallocate(zyn)
   !
     call dealloc(prop)
     call sode%deallocSolidt()
  end subroutine computeSpheroidalAdjoint
!---------------------------------------------------------------
!> \brief Deallocate spheroidalIntegration object
!
  subroutine deallocSpheroidalAdjoint(this)
     type (spheroidal_adjoint) :: this
     if (allocated(this%adj)) deallocate(this%adj)
  end subroutine deallocSpheroidalAdjoint
!
end module spheroidalAdjoint

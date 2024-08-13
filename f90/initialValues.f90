! ======================================================================
!  Compute starting radius and initial values for integrations
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
!-------------------------------------------------------
! Provide starting radius and initial values for
! integration of minors and Green function
!-------------------------------------------------------
module initialValues
    use splineEarthmodelCoefficients
    use geminiIntegrationEnvironment
    use nodeEarthmodel
    use errorMessage
    implicit none
!
contains
!----------------------------------------------------------
!> \brief Calculate starting radius at bottom of model
!! for given integration environment
!! Makes use of exponentail decay of wavefield below
!! turning radius where F(r) = dll1-(om*r/v)**2 > 0
!! Procedure:
!! + Search radius from bottom up until F(r) goes
!!   from positive to negative -> turning point
!! + With Q(r) = sqrt(F(r))/r, the solution behaves as
!!   y = A*exp(-int Q(r) dr), or ln (y/A) = -int Q(r) dr. Therefore
!! + integrate Q(r) downwards from turning point,
!!   (into positive range) until integral exceeds -ln(eps)
!!   Hence, ln (y/A) < ln(eps) or y/A < eps.
!! + if the turning point is in the bottom halfspace or
!!   if we run into the halfspace during integration of Q(r),
!!      + we choose a starting radius just below its top and
!!        use analytic solutions for a homogeneous sphere as
!!        initial values;
!!      + since external nodes in the halfspace are not allowed
!!        rstart cannot be above a node in this case
!! + else if the starting radius is above the halfspace
!!      + we start from there
!!      + if a node is below, values there can be set to zero
!! + else if there is no turning point at all (F(r) always positive)
!!      + then we assume the turning point at the surface
!!        and integrate downward from there to find rstart.
!!      + Solutions at nodes below can be set to zero.   
!!  It follow that regardless where the starting radius is, we
!!  can zero the solution at all nodes below this radius.
!!
!!  For toroidal motion, we in addition make sure that we
!!  always stay above or at the core-mantle boundary
!!
!!  intenv:   integration environment
!!  eps:      desired accuracy
!!  mode:     motion type, either 'sph' or 'tor'
!!  errmsg:   error message object
!---------------------------------------------------------------------------
    function startRadiusInitialValues(intenv,epsin,mode,errmsg) result(res)
    type (gemini_integration_environment) :: intenv
    double precision :: epsin
    character (len=*) :: mode
    type (error_message) :: errmsg
    double precision :: res
    character(len=24) :: myname = 'startRadiusInitialValues'
    type (node_earthmodel) :: nem
    double complex, dimension(7) :: zelcon
    double precision :: r,rhs,roc,ro,v,dll1,om,rearth,dr,p,pold,rturn,rmax,eps
    integer :: nl,top,noc
!
    call addTrace(errmsg,myname)
    nem = intenv%nem               ! node earth model
    dll1 = intenv%dll1
    om = intenv%omre
    rearth = .rearth.nem
    rhs = .rhs.nem            ! radius of top of halfspace
    roc = .roc.nem            ! radius of outer core
    noc = .noc.nem            ! layer index of outer core
    eps = 1.d-5               ! keep eps at 1.d-5 for starting radius (works perfect)
!
    if (mode == 'tor' .and. .hascore.nem) then  ! start at CMB for toroidal motion
       r = roc
       nl = noc+1
    else
       r = rhs
       nl = 1
    endif
    dr = (rearth-r)/400.d0
!
!  check if squared vertical slowness is already negative
!
    call evalComplexElasticConstantsSplineEarthmodel(r,nl,ro,zelcon)
    v = dsqrt( real(zelcon(7))/ro )                            ! vs
    if (nem.isfluid.nl) v = dsqrt( real(zelcon(6))/ro )        ! vp
    p = verticalSlownessInitialValues(om,dll1,r,v)
    pold = p
    if (p <= 0.d0) then
       if (mode == 'tor' .and. .hascore.nem) then
          call add(errmsg,1,'Turning point is at CMB, start there',myname)
          res = roc
       else
          call add(errmsg,1,'Turning point is in half space, start there',myname)
          res = rhs-0.1*dr
       endif
       return
    endif
!
!  if not search for turning point above
!
    if (mode =='tor' .and. .hasocean.nem) then
       rmax = .ruc.nem                                ! limit rstart to below ocean bottom
    else
       rmax = rearth
    endif
    do while (p > 0.d0)
       pold = p
       if (dabs(1.d0-r/rmax) < epsilon(1.d0)) goto 2  ! reached surface or ocean bottom?
       r = min(r+dr,rmax)
       call getLayerIndexNodeEarthmodel(nem,r,nl,top)
       call evalComplexElasticConstantsSplineEarthmodel(r,nl,ro,zelcon)
       v = dsqrt( real(zelcon(7))/ro )                                    ! vs
       if (nem.isfluid.nl) v = dsqrt( real(zelcon(6))/ro )                ! vp
       p = verticalSlownessInitialValues(om,dll1,r,v)
    enddo
!
!  turning point above halfspace found: rturn = r-F(r)/F'(r)
!
    call add(errmsg,0,'Turning point above bottommost radius',myname)
    rturn = r-dr*p/(p-pold)
    call integrateSlownessInitialValues(nem,om,dll1,rturn,dr,eps,mode,res)
    return
!
2   call add(errmsg,1,'No turning point at all',myname)
    rturn = rmax
    call integrateSlownessInitialValues(nem,om,dll1,rturn,dr,eps,mode,res)
    return
    end function startRadiusInitialValues
!--------------------------------------------------------------
!> \brief Evaluate vertical slowness
!
    function verticalSlownessInitialValues(om,dll1,r,v) result(res)
    double precision :: res,om,dll1,r,v
    res = dll1/(r*r)-(om/v)**2
    end function verticalSlownessInitialValues
!-----------------------------------------------------------------
!> \brief Integrate slowness function downwards
!
    subroutine integrateSlownessInitialValues(nem,om,dll1,rturn,dr,eps,mode,rs)
    type (node_earthmodel) :: nem
    double precision :: rturn,dr,eps,rs,om,dll1
    character (len=*) :: mode
    double precision :: q,r,rhs,roc,qmax,v,rbot
    double complex, dimension(7) :: zelcon
    double precision :: ro
    integer :: nl,top
!
    q = 0.d0
    r = rturn
    rhs = .rhs.nem
    roc = .roc.nem
    if (mode == 'tor' .and. .hascore.nem) then
       rbot = roc
    else
       rbot = rhs-0.1*dr
    endif
    qmax = -log(eps)
    do while (q < qmax)
       if (r-dr < rbot) then                     ! enough, we are in halfspace (or in core)
          r = rbot                               ! put rs into halfspace (or at CMB)
          exit
       endif
       r = r-dr
       call getLayerIndexNodeEarthmodel(nem,r,nl,top)
       call evalComplexElasticConstantsSplineEarthmodel(r,nl,ro,zelcon)
       v = dsqrt( real(zelcon(7))/ro )                                    ! vs
       if (nem.isfluid.nl) v = dsqrt( real(zelcon(6))/ro )                ! vp
       q = q+dr*sqrt(verticalSlownessInitialValues(om,dll1,r,v))
    enddo
    rs = r
    end subroutine integrateSlownessInitialValues
!------------------------------------------------------------------------
!> brief Initial values at rstart at bottom of model
!!
!!  gem_intenv:    integration environment (in)
!!  rstart:        starting radius (in)
!!  ystart:        initial solution vector (out)
!!  htry:          initial step size (out)
!!  errmsg:        error message object (inout)
!
    subroutine minorsBottomInitialValues(gem_intenv,rstart,ystart,htry,errmsg)
    type (gemini_integration_environment) :: gem_intenv
    double precision :: rstart,htry
    integer :: nl
    double complex, dimension(:) :: ystart
    type (error_message) :: errmsg
    double complex :: zom,zxa2,zxb2,zeta,zetb,za,zb
    double complex, dimension(7) :: zelcon
    double precision :: ro,elp1,el,rsq,v
    logical :: liquid_flag
!
    call addTrace(errmsg,'minorsBottomInitialValues')
!
!  Get equivalent isotropic moduli at the starting radius
!
    nl = gem_intenv%nl
    call evalComplexElasticConstantsSplineEarthmodel(rstart,nl,ro,zelcon)
    liquid_flag = gem_intenv%nem.isfluid.nl
    el=-0.5d0+dsqrt(0.25d0+gem_intenv%dll1)
    elp1 = gem_intenv%dll1
    zom = dcmplx(gem_intenv%omre,gem_intenv%omim)
    rsq = rstart*rstart
    v = dsqrt( real(zelcon(7))/ro )
    if (liquid_flag) v = dsqrt( real(zelcon(6))/ro )
    htry = 0.25*v*mc_two_pid/gem_intenv%omre
!
!  Ratio of spherical Bessel/Hankel(first kind) functions
!
    zxa2 = ro*zom*zom*rsq/(zelcon(6)+4.d0*zelcon(7)/3.d0)
    if (el < epsilon(1.d0)) then
       call zetlInitialValues(zxa2,el,zeta,errmsg)
       if (.level.errmsg == 2) return
    else
       if (abs(zxa2/elp1) <= 1.d0) then
          call zetlInitialValues(zxa2,el,zeta,errmsg)
          if (.level.errmsg == 2) return
       else
          call sphancfInitialValues(sqrt(zxa2),el,zeta,1,errmsg)
          if (.level.errmsg == 2) return
       endif
       if (.not. liquid_flag) then
          zxb2 = ro*zom*zom*rsq/zelcon(7)
          if (abs(zxb2/elp1) <= 1.d0) then
             call zetlInitialValues(zxb2,el,zetb,errmsg)
             if (.level.errmsg == 2) return
          else
             call sphancfInitialValues(sqrt(zxb2),el,zetb,1,errmsg)
             if (.level.errmsg == 2) return
          endif
       endif
    endif
!
!  initial values
!
    if (el < epsilon(1.d0)) then       ! also covers the liquid case
       ystart(1)=-zeta
       ystart(2)=-ro*rstart*zom*zom + 4.d0*zelcon(7)/rstart*zeta
    else
       if (liquid_flag) then
          ystart(1) = -(1.d0-zeta/el)/(zom*zom*ro*rstart)
          ystart(2) = 1.d0/el
       else
          za = zeta/el
          zb = zetb/(el+1.d0)
          ystart(2) = (-za+zb*(za-1.d0))/el
          ystart(1) = zelcon(7)/rstart*(zxb2/el+2.d0*elp1*ystart(2))
          ystart(3) = zelcon(7)/rstart*(-2.d0*ystart(2)+zxb2/elp1*(za-1.d0))
          ystart(4) = -zelcon(7)/rstart*(zxb2/(el*el)*(1.d0-zb)+4.d0*ystart(2))
          ystart(5) = zelcon(7)*zelcon(7)/rsq*(-4.d0*(el-1.d0)*(el+2.d0)*ystart(2) &
               +zxb2/el*(zxb2/elp1-2.d0*(el-1.d0)*(2.d0*el+1.d0)/elp1 &
               -4.d0/(el+1.d0)*za-2.d0/el*zb))
       endif
    endif
    end subroutine minorsBottomInitialValues
!-------------------------------------------------------------------------------
!> \brief Initial values for minors at the surface
!
!  Calculate starting values for minor integration at the
!  surface (z=0). Treats both free surface and full-space models
!  with an homomgeneous spherical shell above the layer stack.
!!  gem_intenv:    integration environment (in)
!!  ystart:        initial solution vector (out)
!!  htry:          initial step size (out)
!!  errmsg:        error message object (inout)
!---------------------------------------------------------------
    subroutine minorsSurfaceInitialValues(gem_intenv,ystart,htry,errmsg)
    type (gemini_integration_environment) :: gem_intenv
    double complex, dimension(:) :: ystart
    double precision :: htry
    type (error_message) :: errmsg
    integer :: nlay
    double complex :: zom,zxa2,zeta,zxb2,zetb,za,zb
    double complex, dimension(7) :: zelcon
    double precision :: elp1,ro,el,rsq,rearth,v
    logical :: liquid_flag
!
    call addTrace(errmsg,'minorsSurfaceInitialValues')
    nlay = .nlay.(gem_intenv%nem)
    elp1 = gem_intenv%dll1
    el=-0.5d0+dsqrt(0.25d0+elp1)
    rearth = .rearth.(gem_intenv%nem)
    zom = dcmplx(gem_intenv%omre,gem_intenv%omim)
    liquid_flag = gem_intenv%nem.isfluid.nlay
    call evalComplexElasticConstantsSplineEarthmodel(rearth,nlay,ro,zelcon)
    v = dsqrt( real(zelcon(7))/ro )
    if (liquid_flag) v = dsqrt( real(zelcon(6))/ro )
    htry = 0.25*v*mc_two_pid/gem_intenv%omre
!
!  outer sphere exists
!
    if (sec_fsflag == 1) then
       rsq = rearth*rearth
!
!  Ratio of spherical Hankel of the second kind functions
!
       zxa2 = ro*zom*zom*rsq/(zelcon(6)+4.d0*zelcon(7)/3.d0)
       call sphancfInitialValues(sqrt(zxa2),el,zeta,2,errmsg)
       if (.level.errmsg == 2) return
       if (.not. liquid_flag) then
          zxb2 = ro*zom*zom*rsq/zelcon(7)
          call sphancfInitialValues(sqrt(zxb2),el,zetb,2,errmsg)
          if (.level.errmsg == 2) return
       endif
!
!  Starting values for integration DOWNWARDS.
!
       if (el < epsilon(1.d0)) then
          ystart(1)=-zeta
          ystart(2)=-ro*rearth*zom*zom + 4.d0*zelcon(7)/rearth*zeta
       else
          if (liquid_flag) then
             ystart(1) = -(1.d0-zeta/el)/(zom*zom*ro*rearth)
             ystart(2) = 1.d0/el
          else
             za = zeta/el
             zb = zetb/(el+1.d0)
             ystart(2) = (-za+zb*(za-1.d0))/el
             ystart(1) = zelcon(7)/rearth*(zxb2/el+2.d0*elp1*ystart(2))
             ystart(3) = zelcon(7)/rearth*(-2.d0*ystart(2)+zxb2/elp1*(za-1.d0))
             ystart(4) = -zelcon(7)/rearth*(zxb2/(el*el)*(1.d0-zb)+4.d0*ystart(2))
             ystart(5) = zelcon(7)*zelcon(7)/rsq*(-4.d0*(el-1.d0)*(el+2.d0)*ystart(2) &
                  +zxb2/el*(zxb2/elp1-2.d0*(el-1.d0)*(2.d0*el+1.d0)/elp1 &
                  -4.d0/(el+1.d0)*za-2.d0/el*zb))
          endif
       endif
!
!  free surface
!
    else
       if (el < epsilon(1.0)) then
          ystart(1) = 1.d0
          ystart(2) = 0.d0
       else
          if (liquid_flag) then
             ystart(1) = 1.d0
             ystart(2) = 0.d0
          else
             ystart(1) = 0.d0
             ystart(2) = 1.d0
             ystart(3) = 0.d0
             ystart(4) = 0.d0
             ystart(5) = 0.d0
          endif
       endif
    endif
    end subroutine minorsSurfaceInitialValues
!--------------------------------------------------------------------
!> \brief Initial values for toroidal integration at bottom
!!  gem_intenv:    integration environment (in)
!!  rstart:        starting radius (in)
!!  ystart:        initial solution vector (out)
!!  htry:          initial step size (out)
!!  errmsg:        error message object (inout)
!
    subroutine toroidalBottomInitialValues(gem_intenv,rstart,ystart,htry,errmsg)
    type (gemini_integration_environment) :: gem_intenv
    double precision :: rstart,htry
    double complex, dimension(:) :: ystart
    type (error_message) :: errmsg
    character(len=27) :: myname = 'toroidalBottomInitialValues'
    integer :: nl
    double complex :: zom,zxb2,zetb
    double complex, dimension(7) :: zelcon
    double precision :: ro,el,elp1,rsq,v
!
    call addTrace(errmsg,myname)
    nl = gem_intenv%nl
    call evalComplexElasticConstantsSplineEarthmodel(rstart,nl,ro,zelcon)
    if (.level.errmsg == 2) return
    el=-0.5d0+dsqrt(0.25d0+gem_intenv%dll1)
    elp1 = gem_intenv%dll1
    zom = dcmplx(gem_intenv%omre,gem_intenv%omim)
    rsq = rstart*rstart
    v = dsqrt( real(zelcon(7))/ro )
    htry = 0.25*v*mc_two_pid/gem_intenv%omre
!
    if (rstart > .roc.(gem_intenv%nem)) then
       zxb2=ro*zom*zom*rsq/zelcon(7)
       if (abs(zxb2/elp1) < 1.d0) then
          call zetlInitialValues(zxb2,el,zetb,errmsg)
          if (.level.errmsg == 2) return
       else
          call sphancfInitialValues(sqrt(zxb2),el,zetb,1,errmsg)
          if (.level.errmsg == 2) return
       endif
       ystart(1) = 1.d0/el
       ystart(2) = zelcon(7)/(rstart*el)*(el-1.d0-zetb)
    else
       ystart(1) = 1.d0/el
       ystart(2) = 0.d0
    endif
    end subroutine toroidalBottomInitialValues
!-------------------------------------------------------------------------------
!> \brief Initial values for minors at the surface
!
!  Calculate starting values for toroidal integration at the
!  surface (z=0). Treats both free surface and full-space models
!  with an homomgeneous spherical shell above the layer stack.
!!  gem_intenv:    integration environment
!!  ystart:        initial solution vector (out)
!!  htry:          initial step size (out)
!!  rsurf:         surface radius is either true surface or ocean bottom (out)
!!  nl:            layer where downward integration starts (either nlay or nlay-1 if model has ocean) (out)
!!  errmsg:        error message object (inout)
!----------------------------------------------------------------------------------
    subroutine toroidalSurfaceInitialValues(gem_intenv,ystart,htry,rsurf,nl,errmsg)
    type (gemini_integration_environment) :: gem_intenv
    double precision :: htry
    double complex, dimension(:) :: ystart
    type (error_message) :: errmsg
    type (node_earthmodel) :: nem
    character(len=27) :: myname = 'toroidalSurfaceInitialValues'
    integer :: nlay,nl
    double complex :: zom,zxb2,zetb
    double complex, dimension(7) :: zelcon
    double precision :: ro,el,rsq,rearth,v,rsurf
!
    call addTrace(errmsg,myname)
    nem = gem_intenv%nem
    nlay = .nlay.nem
    if (.hasocean.nem) then       ! start integration at ocean bottom
       nl = nlay-1
       rsurf = .ruc.nem
    else
       nl = nlay
       rsurf = .rearth.nem
    endif
    el=-0.5d0+dsqrt(0.25d0+gem_intenv%dll1)
    zom = dcmplx(gem_intenv%omre,gem_intenv%omim)
    call evalComplexElasticConstantsSplineEarthmodel(rsurf,nl,ro,zelcon)
    if (.level.errmsg == 2) return
    v = dsqrt( real(zelcon(7))/ro )
    htry = 0.25*v*mc_two_pid/gem_intenv%omre
!
!  outer sphere exists (l=0)
!
    if (sec_fsflag == 1) then
       rsq = rsurf*rsurf
       zxb2 = ro*zom*zom*rsq/zelcon(7)
       call sphancfInitialValues(sqrt(zxb2),el,zetb,2,errmsg)
       if (.level.errmsg == 2) return
       ystart(1) = 1.d0/el
       ystart(2) = zelcon(7)/(rsurf*el)*(el-1.d0-zetb)
    else
       ystart(1) = 1.d0
       ystart(2) = 0.d0
    endif
    end subroutine toroidalSurfaceInitialValues
!--------------------------------------------------------------------------------
!> brief Spherical Hankel function ratio (first or second kind)
!
!  Computes ratio of (complex) spherical Hankel-functions
!  of the first kind using a continued fraction
! 
!    z(n):= x*h(n+1)/h(n)
!
!    z(n)=n+1-zi*x-zi*CF
!    CF = a1/(b1 +) a2/(b2 + )......
!    a(j)=(j-1/2)**2-(n+1/2)**2, b(j)=2*(x+j*zi)
!
!  The continued fraction
!  is evaluated with the MODIFIED LENTZ'S METHOD. See Press, W.H. et
!  al, 1992, Numerical Recipes, 2nd Edition, Cambridge, Chapter 5.2,
!  p. 169 for further information.
!---------------------------------------------------------------------------------
    subroutine sphancfInitialValues(zx,el,zf,kindhf,errmsg)
    double precision :: tiny,el,eps,elph2
    integer :: maxit,j,kindhf
    parameter(maxit = 10000,tiny = 1.d-30,eps = 1.d-6)
    double complex :: zx,zf,zd,zc,zdelta,za,zb
    type (error_message) :: errmsg
!
!  Evaluate CF
!
    zf = tiny
    zc = zf
    zd = 0.d0
    elph2 = (el+.5d0)**2
    do j = 1,maxit
        zb = 2.d0*(zx+j*mc_cid)
        za = (dble(j)-.5d0)**2-elph2
        zd = zb+za*zd
        if(zabs(zd).lt.tiny) zd = tiny
        zc = zb+za/zc
        if(zabs(zc).lt.tiny) zc = tiny
        zd = 1.d0/zd
        zdelta = zc*zd
        zf = zf*zdelta
        if (zabs(zdelta-1.d0) .lt. eps) goto 11
    enddo    
!
!  no convergence
!
    call add(errmsg,2,'No convergence for spherical Hankel function','sphancfInitialValues')
    print *, '<sphancfInitialValues>: zx,el,eps,zf,maxit,tiny = ',zx,el,eps,zf,maxit,tiny
    return
!
11  if (kindhf == 1) then
       zf = el+1.d0-mc_cid*zx-mc_cid*zf
    else if (kindhf == 2) then
       zf = el+1.d0+mc_cid*zx+mc_cid*zf
    else
       call add(errmsg,2,'Invalid kind for spherical Hankel function','sphancfInitialValues')
    endif
    end subroutine sphancfInitialValues
!--------------------------------------------------------------------------
!>  \brief Routine zetl required for initial values at rstart
!
!  Computes ratio of (complex) spherical Bessel-functions 
!    z(n):= x*j(n+1)/j(n)
!
!  by calculation of the recurrence formula
!    z(n)=x**2/(2*n+3 -z(n+1)).
!  (Takeuchi and Saito (1972), Methods in computional Physics, p.243)
!
!  The recurrence formula is treated as a continued fraction
!  and evaluated with the MODIFIED LENTZ'S METHOD. See Press, W.H. et
!  al, 1992, Numerical Recipes, 2nd Edition, Cambridge, Chapter 5.2,
!  p. 169 for further information.
!
!  Author: J.R. Dalkolmo
!--------------------------------------------------------------------
    subroutine zetlInitialValues(zx2,el,z,errmsg)
    double precision :: tiny,el,eps
    integer :: maxit,i
    parameter(maxit=10000,tiny=1.d-30,eps=1.d-6)
    double complex :: zx2,z,zd,zc,zdelta 
    type (error_message) :: errmsg
!
! First iteration
!
    z = tiny
    zd = 2.d0*el+3.d0
    zc = zd+zx2/z
    zd = 1.d0/zd
    zdelta = zc*zd
    z = z*zdelta
!
! Remaining iterations
!
    do i = 2,maxit
       zd = 2.d0*(el+dble(i))+1.d0-zx2*zd
       if (zd.eq.0.d0) zd = tiny
       zc = 2.d0*(el+dble(i))+1.d0-zx2/zc
       if (zc.eq.0.d0) zc = tiny
       zd = 1.d0/zd
       zdelta = zc*zd
       z = z*zdelta
       if (abs(abs(zdelta)-1.d0) .lt. eps) return
    enddo
    call add(errmsg,2,'Iteration failed to converge in ZETL','zetlInitialValues')
    print *, 'zx2,el,eps,z,maxit,tiny = ',zx2,el,eps,z,maxit,tiny
    end subroutine zetlInitialValues
!
end module initialValues

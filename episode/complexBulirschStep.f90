! ======================================================================
!  Bulirsch-Stoer integration engine for complex-valued solution vectors
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
!-------------------------------------------------------------------------------------
!  Engine for the adaptive step-size Bulirsch-Stoer integration method
!  with Richardson extrapolation. Implements procedures as defined
!  in the base class "baseStepEngine".  
!  A BS-step comprises several modified midpoint steps with
!  different subdivisions of the interval and subsequent extrapolation
!  to step zero. Differences between extrapolations done for sequences
!  of subdivisions of different depths are used for an error estimate.
!  Step is successful if error estimate is below predefined threshold.
!  The BS-Step must be carried out over a smooth model region. 
!  Here version for double complex solution vector
!------------------------------------------------------------------
module complexBulirschStep
    use parametersBulirschStep
    use baseStepEngine
    implicit none
    type, extends(base_step_engine) :: complex_bulirsch_step    ! extends base class
       private
       integer :: donetry                                       ! number of done subdivision attempts
       double complex, dimension(:,:), allocatable :: ysol      ! sequence of solution vectors for all
                                                                ! attempted  subdivision levels
       double complex, dimension(:), allocatable :: yext        ! final extrapolated solution
       double complex, dimension(:), allocatable :: yerr        ! error of final extrapolated solution
       double complex, dimension(:,:), allocatable :: rhs       ! values of right hand side for all 
                                                                ! attempted subdivisions (2nd index)
       double complex, dimension(:,:,:), allocatable :: sm      ! values of system matrix for all 
                                                                ! attempted subdivisions (last index)
    contains
       procedure :: create => createComplexBulirsch             ! mapping of base class procedures 
       procedure :: dealloc => deallocComplexBulirsch           ! concrete implementations
       procedure :: doFullStep => doFullStepComplexBulirsch
       procedure :: adjustStepSize => adjustStepSizeComplexBulirsch
       procedure :: getSolution => getSolutionComplexBulirsch
       procedure :: doLinearModifiedMidpointComplexBulirsch
       procedure :: doNonLinearModifiedMidpointComplexBulirsch
       procedure :: extrapolateComplexBulirsch
       procedure :: printResultComplexBulirsch
    end type complex_bulirsch_step
!
contains
!--------------------------------------------------------------------
!> \brief Create a bulirschStep object with complex variables
!   sode:     polymorphic object representing ODS system
!   eps:      desired accuracy
!   secrecy:  anti-debugging level
!
    subroutine createComplexBulirsch(this,sode,eps,secrecy)
    class (complex_bulirsch_step) :: this
    class (base_ode_system), target :: sode
    double precision :: eps
    integer :: secrecy
    character(len=21) :: myname = 'createComplexBulirsch'
!
    call new(this%errmsg,myname,30)
    this%sode => sode
    this%eps = eps
    this%secrecy = secrecy
    this%nvar = this%sode%getNvar()
    allocate(this%ysol(this%nvar,maxtry_bulirsch_step))
    allocate(this%yext(this%nvar),this%yerr(this%nvar))
    if (this%sode%linear) then
       allocate(this%rhs(this%nvar,0:kgv_bulirsch_step))
       allocate(this%sm(this%nvar,this%nvar,0:kgv_bulirsch_step))
    endif
    end subroutine createComplexBulirsch
!-----------------------------------------------------------------------------------
!> \brief Deallocate object
!
    subroutine deallocComplexBulirsch(this)
    class (complex_bulirsch_step) :: this
!
    if (allocated(this%ysol)) deallocate(this%ysol)
    if (allocated(this%sm)) deallocate(this%sm)
    if (allocated(this%rhs)) deallocate(this%rhs)
    if (allocated(this%yext)) deallocate(this%yext)
    if (allocated(this%yerr)) deallocate(this%yerr)
    if (associated(this%sode)) nullify(this%sode)
    if (this%secrecy <= secrecy_bulirsch_step) call print(this%errmsg)
    call dealloc(this%errmsg)
    end subroutine deallocComplexBulirsch
!--------------------------------------------------------------------------
!> \brief Do a full Bulirsch step with several subdivisions
!! after which either a final solutions below the error bound is achieved
!! or the step size was too large and we need to shrink step size
!! for another attempt.
!  h:      trial step size (in)
!  xstart: begin of step (in)
!  yr:     real part of initial solution vector (in)
!  yi:     imaginary part of initial solution vector (in)
!
    subroutine doFullStepComplexBulirsch(this,h,xstart,yr,yi)
    class (complex_bulirsch_step) :: this
    double precision :: h,xstart
    double precision, dimension(:) :: yr,yi
    integer :: nfirst,ntry,j
    double precision :: errmax
    double complex, dimension(:), allocatable :: zyy,zdy
!
!  initialize success and number of used subdivisions
!
    this%success = .false.
    this%donetry = 0
!
!  initiallize first mintry steps and
!  allocate space for extrapolated solution
!
    ntry = mintry_bulirsch_step
    nfirst = 1
    allocate(zyy(this%nvar),zdy(this%nvar))        ! do this out of goto loop
!
!  first do mintry attempts without checking for convergence
!
1   do j = nfirst,ntry
       if (this%sode%linear) then
          call this%doLinearModifiedMidpointComplexBulirsch(j,h,xstart,yr,yi)
       else
          call this%doNonLinearModifiedMidpointComplexBulirsch(j,h,xstart,yr,yi)
       endif
    enddo
!
!  extrapolate solution to step size zero
!
    call this%extrapolateComplexBulirsch(h,zyy,zdy)
    errmax = sqrt(sum(abs(zdy)**2)/sum(abs(zyy)**2))         ! length ratio of error vector and solution vector
!    errmax = maxval(abs(zdy)/abs(zyy))
    if (this%secrecy <= secrecy_bulirsch_step) call this%printResultComplexBulirsch(xstart,h,errmax,zyy)
    if (errmax/this%eps < 1.d0) then                        ! Accuracy level reached
       this%yext = zyy                                      ! Store extrapolated result
       this%yerr = zdy                                      ! and error
       this%success = .true.
       deallocate(zyy,zdy)
       return
    endif
!
!  error is still above threshold, because we arrive here
!  try further subdivisions or resign
!
    if (ntry < maxtry_bulirsch_step) then                   ! try more subdivisions of same step
       ntry = ntry+1; nfirst = ntry; goto 1
    else
       this%success = .false.                               ! step size too large, go back to calling routine
       if (this%secrecy <= secrecy_bulirsch_step) then
          print *,'doFullStepComplexBulirsch: xstart,errmax = ',xstart,errmax
          print *,'Starting values: ',(yr(j),yi(j),j=1,this%nvar)
          print *,'Solution error: zdy :',zdy(1:this%nvar)
          print *,'Extrapolated solution: zyy :',zyy(1:this%nvar)
       endif
       deallocate(zyy,zdy)
    endif
    end subroutine doFullStepComplexBulirsch
!---------------------------------------------------------------------------------------------------------------
!> \brief Do a linear modified midpoint calculation for one subdivision attempt
!  jtry:     count of subdivision attempt (in)
!  h:        step size (in)
!  xstart:   begin of interval (in)
!  yr:       real part of initial solution vector (in)
!  yr:       imaginary part of initial solution vector (in)
!  No error checking here for efficiency reasons
!
    subroutine doLinearModifiedMidpointComplexBulirsch(this,jtry,h,xstart,yr,yi)
    class (complex_bulirsch_step) :: this
    integer :: jtry
    double precision :: h,xstart
    double precision, dimension(:) :: yr,yi
    double complex, dimension(:), allocatable :: z0,z1,z2
    integer :: j,ix,nstep,m,ios1,ios2
    double precision :: x
    double precision :: hdiv
!    double complex :: hdiv,dcone = dcmplx(1.d0,0.d0)            ! required to use BLAS
!
!  precalculate required values for system matrix and rhs
!
    do j = 1,nxdiv_bulirsch_step(jtry)
       ix = ixdiv_bulirsch_step(j,jtry)
       x = xstart+(ix*h)/kgv_bulirsch_step
       call this%sode%getSysmatComplex(x,this%sm(:,:,ix),ios2)
       call this%sode%getRhsComplex(x,this%rhs(:,ix),ios1)
    enddo
!
!  steps and step length for current subdivision
!
    nstep = seqdiv_bulirsch_step(jtry)
    hdiv = h/nstep
!
!  Modified-midpoint step:
!  allocate space for intermediate values
!  Initialize: z0 = y(x), z1 = z0+h*f(x,z0)
!  Steps: z_{m+1} = z_{m-1}+2h*f(x+mh,zm) for m=1,nstep-1
!  final value at x+hdiv: = 0.5*(z_n+z_{n-1}+h*f(x+H,z_n))
!  after loop z1 corresponds to z_n
!
    allocate(z0(this%nvar),z1(this%nvar),z2(this%nvar))
    z0 = dcmplx(yr,yi)
!    call zcopy(this%nvar,this%rhs(:,0),1,z1,1)                                        ! z1 = rhs
!    call zgemv('N',this%nvar,this%nvar,hdiv,this%sm(:,:,0),this%nvar,z0,1,hdiv,z1,1)  ! z1 = hdiv*(sm*z0)+hdiv*z1
!    call zaxpy(this%nvar,dcone,z0,1,z1,1)                                             ! z1 = 1.d0*z0+z1
    z1 = z0+hdiv*matmul(this%sm(:,:,0),z0)+hdiv*this%rhs(:,0)                         ! => z1 = z0+hdiv*(sm*z0)+hdiv*rhs
    do m = 1,nstep-1
       ix = m*kgv_bulirsch_step/nstep
!       call zcopy(this%nvar,this%rhs(:,ix),1,z2,1)                                                  ! z2 = rhs
!       call zgemv('N',this%nvar,this%nvar,2.d0*hdiv,this%sm(:,:,ix),this%nvar,z1,1,2.d0*hdiv,z2,1)  ! z2 = 2*hdiv*(sm*z1)+2*hdiv*z2
!       call zaxpy(this%nvar,dcone,z0,1,z2,1)                                                        ! z2 = 1.d0*z0+z2
       z2 = z0+2.d0*hdiv*matmul(this%sm(:,:,ix),z1)+2.d0*hdiv*this%rhs(:,ix)                        ! => z2 = z0+2*hdiv*(sm*z1)+2*hdiv*rhs
!       call zcopy(this%nvar,z1,1,z0,1)                                                ! z0 = z1
!       call zcopy(this%nvar,z2,1,z1,1)                                                ! z1 = z2
       z0 = z1; z1 = z2
    enddo
    ix = kgv_bulirsch_step
!    call zcopy(this%nvar,this%rhs(:,ix),1,z2,1)                                            ! z2 = rhs
!    call zgemv('N',this%nvar,this%nvar,hdiv,this%sm(:,:,ix),this%nvar,z1,1,hdiv,z2,1)      ! z2 = hdiv*(sm*z1)+hdiv*z2
!    call zaxpy(this%nvar,dcone,z0,1,z2,1)                                                  ! z2 = 1.d0*z0+z2
!    call zaxpy(this%nvar,dcone,z1,1,z2,1)                                                  ! z2 = 1.d0*z1+z2
!    call zscal(this%nvar,0.5d0*dcone,z2,1)                                                 ! z2 = 0.5*z2
!    call zcopy(this%nvar,z2,1,this%ysol(:,jtry),1)                                         ! ysol = z2
    this%ysol(:,jtry) = 0.5d0*(z1+z0+hdiv*matmul(this%sm(:,:,ix),z1)+hdiv*this%rhs(:,ix))  ! => ysol = 0.5*(z1+z0+hdiv*(sm*z1)+hdiv*rhs)
    deallocate(z0,z1,z2)
    this%donetry = jtry
    end subroutine doLinearModifiedMidpointComplexBulirsch
!---------------------------------------------------------------------------------------------------------------
!> \brief Do a non linear modified midpoint calculation for one subdivision attempt
!  jtry:     count of subdivision attempt (in)
!  h:        step size (in)
!  xstart:   begin of interval (in)
!  yr:       real part of initial solution vector (in)
!  yr:       imaginary part of initial solution vector (in)
!
    subroutine doNonLinearModifiedMidpointComplexBulirsch(this,jtry,h,xstart,yr,yi)
    class (complex_bulirsch_step) :: this
    integer :: jtry
    double precision :: h,xstart
    double precision, dimension(:) :: yr,yi
    double complex, dimension(:), allocatable :: z0,z1,z2
    double complex, dimension(:), allocatable :: zderiv
    integer :: nstep,m,ios1,ios2,ios3
    double precision :: hdiv,x
!
!  steps and step length for current subdivision
!
    nstep = seqdiv_bulirsch_step(jtry)
    hdiv = h/nstep
!
!  Modified-midpoint step:
!  space for intermediate values
!  Initialize: z0 = y(x), z1 = z0+h*f(x,z0)
!  Steps: z_{m+1} = z_{m-1}+2h*f(x+mh,zm) for m=1,nstep-1
!  final value at x+hdiv: = 0.5*(z_n+z_{n-1}+h*f(x+H,z_n))
!  after loop z1 corresponds to z_n
!
    allocate(z0(this%nvar),z1(this%nvar),z2(this%nvar),zderiv(this%nvar))
    z0 = dcmplx(yr,yi)
    x = xstart
    call this%sode%derivativeComplex(x,z0,zderiv,ios1)
    z1 = z0+hdiv*zderiv
    do m = 1,nstep-1
       x = xstart+m*hdiv
       call this%sode%derivativeComplex(x,z1,zderiv,ios2)
       z2 = z0+2.*hdiv*zderiv
       z0 = z1; z1 = z2
    enddo
    x = xstart+h
    call this%sode%derivativeComplex(x,z1,zderiv,ios3)
    this%ysol(:,jtry) = 0.5*(z1+z0+hdiv*zderiv)
    deallocate(z0,z1,z2,zderiv)
    this%donetry = jtry
    end subroutine doNonLinearModifiedMidpointComplexBulirsch
!---------------------------------------------------------------
!> \brief Extrapolate solution attempts to step size zero by
!!  polynomial interpolation of degree N-1 through N points
!!  defined by arrays x(1-N) and y(1-N) using Neville's
!!  algorithm (see Numerical Recipes 3.1)
!!  Returns value yy of interpolating polynomial at point xx = 0
!!  and an error estimate dy.
!!  hs:     step size (in)
!!  yy:     extrapolated function vector, allocated by calling proc (out)
!!  dy:     error of extrapolated vector, allocated by calling proc (out)
!------------------------------------------------------------------------
    subroutine extrapolateComplexBulirsch(this,hs,yy,dy)
    class (complex_bulirsch_step) :: this
    double precision :: hs
    double complex, dimension(:) :: yy,dy                  ! dimension = nvar
    double precision, dimension(:), allocatable :: x
    double complex, dimension(:,:), allocatable :: c,d
    double complex, dimension(:), allocatable :: u
    integer :: n,js,i,m,nvar,j
    double precision :: h
!
    n = this%donetry
    nvar = this%nvar
    allocate(x(n),c(nvar,n),d(nvar,n),u(nvar))
    x(1:n) = hs/seqdiv_bulirsch_step(1:n)
!
!  since x(i) are decreasing from big to small h
!  the extrapolation point xx=0 is always closest to the n-th node
!  calculate columns of Neville scheme
!  first columns, c and d's equal to y's
!  start with function value at point closest to x=0
!
    js = n
    forall (j = 1:n) c(:,j) = this%ysol(:,j)
    d = c
    yy(:) = this%ysol(:,js)
!
!  do recursion for next columns, number of c's and d's decreases by one
!  when moving from left to right (i.e. nc = nd = n-m)
!  d(m,i) = (x(i+m)-xx)*(c(m-1,i+1)-d(m-1,i))/(x(i)-x(i+m))
!  c(m,i) = (x(i)-xx)*(c(m-1,i+1)-d(m-1,i))/(x(i)-x(i+m))
!
    do m = 1,n-1
        do i = 1,n-m
            h = x(i)-x(i+m)
            u(:) = (c(:,i+1)-d(:,i))/h
            c(:,i) = x(i)*u(:)
            d(:,i) = x(i+m)*u(:)
        enddo
        yy(:) = yy(:)+d(:,js-1)
        dy(:) = d(:,js-1)
        js = js-1
    enddo
    deallocate(c,d,x,u)
    end subroutine extrapolateComplexBulirsch
!---------------------------------------------------------------------
!> \brief Adjust step size
!  h:       step size (inout)
!  success: flag fo success of integration step (in)
!
    subroutine adjustStepSizeComplexBulirsch(this,h,success)
    class (complex_bulirsch_step) :: this
    double precision :: h
    logical :: success
    if (success) then
       if (this%donetry > opttry_bulirsch_step) then
          h = h*shrink_bulirsch_step
       else if (this%donetry < opttry_bulirsch_step) then
          h = h*grow_bulirsch_step
       endif
    else
       h = 0.5d0*h/2**((maxtry_bulirsch_step-opttry_bulirsch_step)/2)
    endif
    end subroutine adjustStepSizeComplexBulirsch
!---------------------------------------------------------------------
!> \brief Get solution at end of step (complex)
!
    subroutine getSolutionComplexBulirsch(this,yr,yi)
    class (complex_bulirsch_step) :: this
    double precision, dimension(:) :: yr,yi
    yr = real(this%yext); yi = imag(this%yext)
    end subroutine getSolutionComplexBulirsch
!--------------------------------------------------------------
!> \brief print results from each subdivision step
!
    subroutine printResultComplexBulirsch(this,xstart,h,errmax,zy)
    class (complex_bulirsch_step) :: this
    double precision :: xstart,h,errmax
    double complex, dimension(:) :: zy
    print *,'Result of certain subdivision of complex Bulirsch Step'
    print *,'xstart: ',xstart
    print *,'h:      ',h
    print *,'ntry:   ',this%donetry
    print *,'err:    ',errmax
    print *,'real(ysol) = ',real(zy)
    print *,'imag(ysol) = ',aimag(zy)
    end subroutine printResultComplexBulirsch
!
end module complexBulirschStep

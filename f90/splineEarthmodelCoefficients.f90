! ======================================================================
!  Holds spline expansion coefficients for material parameters
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
!-------------------------------------------------------------------------
!  Module holding earthmodel spline coefficients for efficient access.
!
module splineEarthmodelCoefficients
    use nodeEarthmodel
    use cubicSpline
    implicit none
    double precision, dimension(:), allocatable :: sec_ro        ! density at radial knots, (g/cm^3)
    double precision, dimension(:), allocatable :: sec_ro2       ! spline coefficients of density
    double precision, dimension(:,:), allocatable :: sec_qinv    ! inverse QK,QM at radial knots
    double precision, dimension(:,:), allocatable :: sec_qinv2   ! spline coefficients of inverse QK, QM
    double complex, dimension(:,:), allocatable :: sec_elcon     ! 7 elastic constants at knots (elcon,knot)
                                                                 ! in the order A,C,F,L,N,KAPPA,MU
    double complex, dimension(:,:), allocatable :: sec_elcon2    ! spline coefficients for 7 elastic constants at knots (elcon,knot)
                                                                 ! in the order A,C,F,L,N,KAPPA,MU
    double precision, dimension(:), allocatable :: sec_rk        ! nodeEarthmodel radial nodes
    integer, dimension(:), allocatable :: sec_iktop              ! nodeEarthmodel radial node index of top of layer
    double precision, dimension(:), allocatable :: sec_drk       ! distances between nodes per layer if equidistant
    integer :: sec_nlay                                          ! number of layers in earth model
    integer :: sec_fsflag                                        ! full space flag
    integer :: sec_equidisnodes                                  ! equidistant nodes flag
contains
!---------------------------------------------------------------
!> \brief Carry out spline interpolation
!
    subroutine computeSplineEarthmodelCoefficients(nem,f,attmode,errmsg)
    type (node_earthmodel) :: nem
    double precision :: f
    character(len=*) :: attmode
    type (error_message) :: errmsg,errmsg2
    integer :: nl,j1,j2,nk,k
    character (len=35) :: myname = 'computeSplineEarthmodelCoefficients'
!
    call addTrace(errmsg,myname)
    nk = .nk.nem
    sec_nlay = .nlay.nem
    allocate(sec_ro(nk),sec_qinv(2,nk),sec_elcon(7,nk))
    call complexElasticConstantsAtKnotsNodeEarthmodel(nem,attmode,f,sec_ro,sec_qinv,sec_elcon,errmsg)
    if (.level.errmsg == 2) then
        deallocate(sec_ro,sec_qinv,sec_elcon)
        return
    endif
    allocate(sec_ro2(nk),sec_rk(nk),sec_iktop(sec_nlay),sec_drk(sec_nlay))
    allocate(sec_elcon2(7,nk),sec_qinv2(2,nk))
    sec_rk = getRkArrayNodeEarthmodel(nem)
    sec_fsflag = .fsflag.nem
    sec_equidisnodes = .equidisnodes.nem
    sec_iktop = getIktopArrayNodeEarthmodel(nem)
    sec_drk = getDrkArrayNodeEarthmodel(nem)
!
!  do not spline halfspace ( = layer 1)
!
    do nl = 2,.nlay.nem
        j1 = nem.ikbot.nl
        j2 = nem.iktop.nl
        errmsg2 = doubleCubicSpline(sec_rk(j1:j2),sec_ro(j1:j2),sec_ro2(j1:j2))
        if (.level.errmsg2 == 2) goto 1
        call dealloc(errmsg2)
        do k = 1,2
           errmsg2 = doubleCubicSpline(sec_rk(j1:j2),sec_qinv(k,j1:j2),sec_qinv2(k,j1:j2))
           if (.level.errmsg2 == 2) goto 1
           call dealloc(errmsg2)           
        enddo
        do k = 1,7
           errmsg2 = doubleComplexCubicSpline(sec_rk(j1:j2),sec_elcon(k,j1:j2),sec_elcon2(k,j1:j2))
           if (.level.errmsg2 == 2) goto 1
           call dealloc(errmsg2)
        enddo
    enddo
    return
1   call add(errmsg,2,'Cubic spline interpolation not successful',myname)
    call print(errmsg2)
    deallocate(sec_ro2,sec_qinv2,sec_elcon2)
    end subroutine computeSplineEarthmodelCoefficients
!----------------------------------------------------------------------
!> \brief Deallocate spline model coefficients
!
    subroutine deallocSplineEarthmodelCoefficients
    if (allocated(sec_ro)) deallocate(sec_ro) 
    if (allocated(sec_qinv)) deallocate(sec_qinv)
    if (allocated(sec_qinv2)) deallocate(sec_qinv2)
    if (allocated(sec_elcon)) deallocate(sec_elcon)
    if (allocated(sec_elcon2)) deallocate(sec_elcon2)
    if (allocated(sec_ro2)) deallocate(sec_ro2)
    if (allocated(sec_rk)) deallocate(sec_rk)
    if (allocated(sec_iktop)) deallocate(sec_iktop)
    if (allocated(sec_drk)) deallocate(sec_drk)
    end subroutine deallocSplineEarthmodelCoefficients
!-----------------------------------------------------------------------
!> \brief Evaluate complex elastic constants for given radius and layer
!!  Do not check r and nl here because of efficiency
!!  r:      radius where material parameters are evaluated (in)
!!  nl:     layer index belonging to given radius (in)
!!  ro:     density (out)
!!  zelcon: array with seven complex elastic constants, A,C,F,L,N,Kappa,Mu (out)
!  This routine is the most often called proedure.
!
    subroutine evalComplexElasticConstantsSplineEarthmodel(r,nl,ro,zelcon,qinv)
    double precision :: r
    integer :: nl
    double precision :: ro
    double precision, dimension(:), optional :: qinv
    double complex, dimension(:) :: zelcon
    double precision :: h,x1,x2
!    double complex :: a,b,c,d
    double precision :: a,b,c,d
    integer :: j1,j2,j,jr,m
!
!  constant parameters in bottom homomgeneous sphere
!
    if (nl == 1) then
        j = sec_iktop(1)
        ro = sec_ro(j)
        zelcon = sec_elcon(:,j)
        if (present(qinv)) qinv = sec_qinv(:,j)
        return
    endif
!
!  above halfspace
!
    j1 = sec_iktop(nl-1)+1
    j2 = sec_iktop(nl)
!
!  locate index j of node with rk(j) < r
!
    if (j2-j1 == 1) then                       ! two nodes in layer, linear case
       j = j1
    else
       if (sec_equidisnodes == 1) then             ! special tretament for equidistant nodes in layers
          j = int((r-sec_rk(j1))/sec_drk(nl))+j1
          if (j == j2) j = j2-1
       else                                        ! non-equidistant case: bisection method
          j = j1-1
          jr = j2+1
          do while (jr-j > 1)
             m = (j+jr)/2
             if (r > sec_rk(m)) then; j = m; else; jr = m; endif
          enddo
          if (j == j1-1) j = j1
       endif
    endif
!
    x1 = sec_rk(j)
    x2 = sec_rk(j+1)
    h = x2-x1
    a = (x2-r)/h
    b = 1.d0-a
    if (j2-j1 == 1) then
       ro = a*sec_ro(j)+b*sec_ro(j+1)
       zelcon = a*sec_elcon(:,j)+b*sec_elcon(:,j+1)
       if (present(qinv)) qinv = a*sec_qinv(:,j)+b*sec_qinv(:,j+1)
    else
       c = h*h*a*(a*a-1.d0)/6.d0
       d = h*h*b*(b*b-1.d0)/6.d0
       ro = a*sec_ro(j)+b*sec_ro(j+1)+c*sec_ro2(j)+d*sec_ro2(j+1)
       zelcon = a*sec_elcon(:,j)+b*sec_elcon(:,j+1)+c*sec_elcon2(:,j)+d*sec_elcon2(:,j+1)
       if (present(qinv)) qinv = a*sec_qinv(:,j)+b*sec_qinv(:,j+1)+c*sec_qinv2(:,j)+d*sec_qinv2(:,j+1)
    endif
!
!  calculate 
!  using BLAS
!
!    zelcon = 0.d0
!    call zaxpy(7,a,sec_elcon(:,j),1,zelcon,1)
!    call zaxpy(7,b,sec_elcon(:,j+1),1,zelcon,1)
!    call zaxpy(7,c,sec_elcon2(:,j),1,zelcon,1)
!    call zaxpy(7,d,sec_elcon2(:,j+1),1,zelcon,1)
    end subroutine evalComplexElasticConstantsSplineEarthmodel
!----------------------------------------------------------------------------
!> \brief Get P-wave velocity at given r and layer nl
!
    function getPVelocitySplineEarthmodel(nl,r) result(res)
    integer :: nl
    double precision :: r,res
    double precision :: ro
    double complex, dimension(7) :: zelcon
    call evalComplexElasticConstantsSplineEarthmodel(r,nl,ro,zelcon)
    res = dsqrt( real(zelcon(6)+4.d0/3.d0*zelcon(7))/ro )
    end function getPVelocitySplineEarthmodel
!----------------------------------------------------------------------------
!> \brief Get S-wave velocity at given r and layer nl
!
    function getSVelocitySplineEarthmodel(nl,r) result(res)
    integer :: nl
    double precision :: r,res
    double precision :: ro
    double complex, dimension(7) :: zelcon
    call evalComplexElasticConstantsSplineEarthmodel(r,nl,ro,zelcon)
    res = dsqrt( real(zelcon(7))/ro )
    end function getSVelocitySplineEarthmodel
!---------------------------------------------------------------------------
end module

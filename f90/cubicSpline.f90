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
!-----------------------------------------------------------
!  Module for cubic spline interpolation
!  Uses Lapack routines (sgtsv etc) for solving tridiagonal system
!-----------------------------------------------------------
 module cubicSpline
    use errorMessage
    implicit none
 contains
!----------------------------------------------------------
!> \brief Cubic spline of real function values
!! \param x Abscissa values array
!! \param y Function values array
!! \param y2 Values of second derivative (output) with y2(1) = y2(n) = 0
!
    function realCubicSpline(x,y,y2) result(errmsg)
    real, dimension(:) :: x,y,y2
    type (error_message) :: errmsg
    real, dimension(:,:), allocatable :: work(:,:)
    integer j,info,n
    character (len=15) :: myname = 'realCubicSpline'
!
    call new(errmsg,myname)
    n = size(x)
    allocate(work(n,3))
!
!  set up tridiagonal matrix
!
    do j = 2,n-1
        work(j-1,1) = (x(j)-x(j-1))/6.
        work(j,2) = (x(j+1)-x(j-1))/3.
        work(j,3) = (x(j+1)-x(j))/6.
        y2(j) = (y(j+1)-y(j))/(x(j+1)-x(j)) - (y(j)-y(j-1))/(x(j)-x(j-1))
    enddo
    work(1,2) = x(2)-x(1)
    work(n,2) = x(n)-x(n-1)
    work(1,3) = 0.0
    work(n-1,1) = 0.0
    y2(1) = 0.0
    y2(n) = 0.0
!
!  solve tridiagonal system for y2
!
    call sgtsv(n,1,work(1,1),work(1,2),work(1,3),y2,n,info)
    if (info .ne. 0) then
            call new(errmsg,2,'tridiagonal system could not be solved',myname)
            deallocate(work)
            return
    endif
    deallocate(work)
    end function realCubicSpline
!----------------------------------------------------------
!> \brief Cubic spline of double precision function values
!! \param x Abscissa values array
!! \param y function values array
!! \param y2 Values of second derivative (output) with y2(1) = y2(n) = 0
!
    function doubleCubicSpline(x,y,y2) result(errmsg)
    double precision, dimension(:) :: x,y,y2
    type (error_message) :: errmsg
    double precision, dimension(:,:), allocatable :: work(:,:)
    integer j,info,n
    character (len=17) :: myname = 'doubleCubicSpline'
!
    call new(errmsg,myname)
    n = size(x)
    allocate(work(n,3))
!
!  set up tridiagonal matrix
!
    do j = 2,n-1
        work(j-1,1) = (x(j)-x(j-1))/6.
        work(j,2) = (x(j+1)-x(j-1))/3.
        work(j,3) = (x(j+1)-x(j))/6.
        y2(j) = (y(j+1)-y(j))/(x(j+1)-x(j)) - (y(j)-y(j-1))/(x(j)-x(j-1))
    enddo
    work(1,2) = x(2)-x(1)
    work(n,2) = x(n)-x(n-1)
    work(1,3) = 0.0
    work(n-1,1) = 0.0
    y2(1) = 0.0
    y2(n) = 0.0
!
!  solve tridiagonal system for y2
!
    call dgtsv(n,1,work(1,1),work(1,2),work(1,3),y2,n,info)
    if (info .ne. 0) then
            call new(errmsg,2,'tridiagonal system could not be solved',myname)
            deallocate(work)
            return
    endif
    deallocate(work)
    end function doubleCubicSpline
!--------------------------------------------------------------------
!  Spline complex function values by splining real and imaginary
!  part separately.
!  x:    abscissa values
!  zy:   complex function values
!  zy2:  complex valued second derivatives (output)
!  y:    workspace of dimension n
!  y2:   workspace of dimension n
!  work: workspace of dimension (n x 3)
!-----------------------------------------------------------------
    function complexCubicSpline(x,zy,zy2) result(errmsg)
    real, dimension(:) :: x
    complex, dimension(:) :: zy,zy2
    type (error_message) :: errmsg
    integer n,i
    real, dimension(:), allocatable :: y,y2
    real, dimension(:,:), allocatable :: work
    complex zi
    character (len=24) :: myname = 'complexCubicSpline'
!
    call new(errmsg,myname) 
    zi = cmplx(0.0,1.0)
    n = size(x)
    allocate(y(n),y2(n),work(n,3))
!
!  real part
!
    do i = 1,n
        y(i) = real(zy(i))
    enddo
    errmsg = realCubicSpline(x,y,y2)
    if (.level.errmsg == 2) then
        call addTraceErrorMessage(errmsg,'doubleComplexCubicSpline')
        deallocate(y,y2,work)
        return
    endif
    do i = 1,n
        zy2(i) = cmplx(y2(i),0.0)
    enddo
!
!  imaginary part
!
    do i = 1,n
        y(i) = aimag(zy(i))
    enddo
    errmsg = realCubicSpline(x,y,y2)
    if (.level.errmsg == 2) then
        call addTraceErrorMessage(errmsg,'doubleComplexCubicSpline')
        deallocate(y,y2,work)
        return
    endif
    do i = 1,n
        zy2(i) = zy2(i)+zi*y2(i)
    enddo
    deallocate(y,y2,work)
    end function complexCubicSpline
!--------------------------------------------------------------------
!  Spline double complex function values by splining real and imaginary
!  part separately.
!  x:    abscissa values
!  zy:   complex function values
!  zy2:  complex valued second derivatives (output)
!  y:    workspace of dimension n
!  y2:   workspace of dimension n
!  work: workspace of dimension (n x 3)
!-----------------------------------------------------------------
    function doubleComplexCubicSpline(x,zy,zy2) result(errmsg)
    double precision, dimension(:) :: x
    double complex, dimension(:) :: zy,zy2
    type (error_message) :: errmsg
    integer n,i
    double precision, dimension(:), allocatable :: y,y2
    double precision, dimension(:,:), allocatable :: work
    double complex zi
    character (len=24) :: myname = 'doubleComplexCubicSpline'
!
    call new(errmsg,myname) 
    zi = dcmplx(0.d0,1.d0)
    n = size(x)
    allocate(y(n),y2(n),work(n,3))
!
!  real part
!
    do i = 1,n
        y(i) = real(zy(i))
    enddo
    errmsg = doubleCubicSpline(x,y,y2)
    if (.level.errmsg == 2) then
        call addTraceErrorMessage(errmsg,'doubleComplexCubicSpline')
        deallocate(y,y2,work)
        return
    endif
    do i = 1,n
        zy2(i) = dcmplx(y2(i),0.d0)
    enddo
!
!  imaginary part
!
    do i = 1,n
        y(i) = dimag(zy(i))
    enddo
    errmsg = doubleCubicSpline(x,y,y2)
    if (.level.errmsg == 2) then
        call addTraceErrorMessage(errmsg,'doubleComplexCubicSpline')
        deallocate(y,y2,work)
        return
    endif
    do i = 1,n
        zy2(i) = zy2(i)+zi*y2(i)
    enddo
    deallocate(y,y2,work)
    end function doubleComplexCubicSpline
 end module

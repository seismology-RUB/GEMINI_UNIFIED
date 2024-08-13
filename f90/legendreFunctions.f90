! ======================================================================
!  Legendre functions and derivatives
! ======================================================================
!----------------------------------------------------------------------------
!   Copyright 2019 Wolfgang Friederich (Ruhr-Universitaet Bochum, Germany)
!
!   This file is part of GEMINI_UNIFIED version 1.0.
!   This file is part of ASKI version 1.2.
!
!   GEMINI_UNIFIED version 1.0 and ASKI version 1.2 are free software:
!   you can redistribute it and/or modify
!   it under the terms of the GNU General Public License as published by
!   the Free Software Foundation, either version 2 of the License, or
!   (at your option) any later version.
!
!   GEMINI_UNIFIED version 1.0 and ASKI 1.2 are is distributed
!   in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU General Public License for more details.
!
!   You should have received a copy of the GNU General Public License
!   along with GEMINI_UNIFIED version 1.0 and ASKI version 1.2.
!   If not, see <http://www.gnu.org/licenses/>.
!-----------------------------------------------------------------
!   Module computing associated Legendre functions up to m=2
!
!----------------------------------------------------------------
module legendreFunctions
    use errorMessage
    implicit none
    interface dealloc; module procedure deallocLegendreFunctions; end interface
    type legendre_functions
        integer :: lmax                                                 ! max value for l
        integer :: mmax                                                 ! max value for m
        double precision :: theta                                       ! theta-value 
        double precision, dimension(:,:), allocatable :: plm            ! Value of Legendre functions for given theta
        double precision, dimension(:,:), allocatable :: dpdt           ! theta-derivatives of Plm
        double precision, dimension(:,:), allocatable :: pdst           ! Plm/sin(theta)
        double precision, dimension(:,:), allocatable :: d2pdt2         ! second theta derivatives
        double precision, dimension(:,:), allocatable :: dtpdst         ! d/dtheta(Plm/sin(theta))
    end type legendre_functions
!
contains
!-------------------------------------------------------------------------------
!  compute Legendre functions up to lmax and mmax at theta
!
    subroutine computeLegendreFunctions(this,lmax,mmax,theta,errmsg)
    type (legendre_functions) :: this
    type (error_message) :: errmsg
    character(len=17) :: myname = 'legendreFunctions'
    integer :: lmax,mmax,isct,j,l,m
    double precision :: theta
    double precision :: ct,st,fact,rst
    logical :: atpoles
  !
    call addTrace(errmsg,myname)
    if (lmax < 0) then
        call add(errmsg,2,'lmax is negative! Must be at last 0',myname)
        return
    endif
    if (mmax < 0 .or. mmax > 2) then
        call add(errmsg,2,'mmax is out of range! 0 <= mmax <= 2',myname)
        return
    endif
  !
    allocate(this%plm(0:lmax,0:mmax))
    allocate(this%dpdt(0:lmax,0:mmax))
    allocate(this%pdst(1:lmax,1:mmax))
    allocate(this%d2pdt2(0:lmax,0:mmax))
    allocate(this%dtpdst(1:lmax,1:mmax))
    this%plm = 0.d0
    this%dpdt = 0.d0
    this%pdst = 0.d0
    this%d2pdt2 = 0.d0
    this%dtpdst = 0.d0
    this%theta = theta
    ct = dcos(theta)
    st = dsin(theta)
    if (st .gt. 0.d0) then
       atpoles = .false.
    else 
       atpoles = .true.
    endif
  !
  !  calculate P(l,m)
  !
    do  m = 0,mmax
       this%plm(m,m)  = 1.d0
       fact = 1.d0
       do  j = 1, m
          this%plm(m,m) = this%plm(m,m)*fact*st         ! no minus sign here (compared to gemini-2.2)
          fact =  fact + 2.d0
       enddo
       this%plm(m+1,m) = ct*(2*m+1)*this%plm(m,m)
       do l = m+2,lmax
          this%plm(l,m) = (ct*dble(2*l-1)*this%plm(l-1,m)-(l+m-1)*this%plm(l-2,m))/dble(l-m)
       enddo
    enddo
  !
  !  Compute derivatives with respect to theta and plm/sin(theta)
  !
    if (.not. atpoles) then
       rst=1./st
       do m = 0,mmax
          this%dpdt(m,m) = m*ct*this%plm(m,m)*rst
          do l = m+1,lmax
             this%dpdt(l,m) = (l*ct*this%plm(l,m)-(l+m)*this%plm(l-1,m))*rst
          enddo
       enddo
       do m = 1,mmax                                    ! only m > 0 needed
          do l = m,lmax
             this%pdst(l,m) = this%plm(l,m)*rst
          enddo
       enddo
    else
       isct = sign(1.d0,ct)
       do l = 0,lmax
          this%dpdt(l,0) = 0.d0
          this%dpdt(l,1) = 0.5*dble(l*(l+1)*isct**l)    ! no minus sign here (compared to gemini-2.2) 
          this%dpdt(l,2) = 0.d0
       enddo
       do l=1,lmax
          this%pdst(l,1) = 0.5*dble(l*(l+1)*isct**(l+1))
          this%pdst(l,2) = 0.d0
       enddo
    endif
  !
  !  Compute second derivatives with respect to theta and d/dtheta(plm/sin(theta))
  !  d/dt(plm/st) = 1/st*dpdt+plm*(-ct/st**2) = 1/st*(dpdt-ct/st*plm)
  !
    if (.not. atpoles) then
       rst=1./st
       do m = 0,mmax
          do l = m,lmax
             this%d2pdt2(l,m) = -ct*rst*this%dpdt(l,m)+(m*m*rst**2-l*(l+1))*this%plm(l,m)
          enddo
       enddo
       do m = 1,mmax                             ! only m > 0 needed
          do l = m,lmax
             this%dtpdst(l,m) = rst*(this%dpdt(l,m)-ct*rst*this%plm(l,m))
          enddo
       enddo
    else
       isct = sign(1.d0,ct)
       do l = 0,lmax
          this%d2pdt2(l,0) = 0.5*dble(l*(l+1)*isct**l) 
          this%d2pdt2(l,1) = 0.d0
          this%d2pdt2(l,2) = 0.25*dble((l+2)*(l+1)*l*(l-1)*isct**l)
       enddo
       do l=1,lmax
          this%dtpdst(l,1) = 0.d0
          this%dtpdst(l,2) = 0.125*dble((l+2)*(l+1)*l*(l-1)*isct**(l+1))
       enddo
       
    endif
    end subroutine computeLegendreFunctions
    !-------------------------------------------------------------------
    !  deallocation
    !
    subroutine deallocLegendreFunctions(this)
    type (legendre_functions) :: this
    if (allocated(this%plm)) deallocate(this%plm)
    if (allocated(this%dpdt)) deallocate(this%dpdt)
    if (allocated(this%pdst)) deallocate(this%pdst)
    if (allocated(this%d2pdt2)) deallocate(this%d2pdt2)
    if (allocated(this%dtpdst)) deallocate(this%dtpdst)
    end subroutine deallocLegendreFunctions
end module legendreFunctions

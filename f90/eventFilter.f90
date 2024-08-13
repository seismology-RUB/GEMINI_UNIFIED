! =================================================================================
!  Compute transfer functions of event filter
! =================================================================================
!----------------------------------------------------------------------------
!   Copyright 2020 Wolfgang Friederich (Ruhr-Universitaet Bochum, Germany)
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
!--------------------------------------------------------------------------------------
!  Routines for calculating the transfer function of Butterworth high and low-pass
!  filters evaluated at either real or complex frequencies.
!-------------------------------------------------------------------------------------
module eventFilter
    use mathConstants
    use butterworthFilter
    use errorMessage
    implicit none
contains
!------------------------------------------------------------------------
!  Compute event filter
!  filtype : type of filter
!  filspecs: array of real numbers specifying filter properties
!  nf1,nf2: range of frequencies indices for which filter is computed
!  df: frequency sampling interval
!  sigma: imaginary part of frequency
!  hfil: filter transfer function versus frequency
!
    subroutine computeEventfilter(filtype,filspecs,nf1,nf2,df,sigma,hfil,errmsg)
    character (len=*) :: filtype
    double precision, dimension(:) :: filspecs
    integer :: nf1,nf2
    double precision :: df,sigma
    double complex, dimension(:), allocatable :: hfil
    type (error_message) :: errmsg
    integer :: nord,nf
    double precision :: fc
    double complex, dimension(:), allocatable :: hlow,hhigh
    character (len=18) :: myname = 'computeEventFilter'
!
    nf = nf2-nf1+1
    allocate(hfil(nf2))
    hfil = 0.d0
    select case (filtype)
    case ('BUTTERWORTH_BANDPASS')
   !
   !  calculate frequency response of HIGH-pass Butterworth filter
   !  at selected complex frequencies with fixed imaginary part
   !
       nord = nint(filspecs(1))
       fc   = filspecs(2)
       call highPassButterworthFilter(nord,nf,fc,(nf1-1)*df,df,sigma,hhigh)
   !
   !  calculate frequency response of LOW-pass Butterworth filter
   !  at selected complex frequencies with fixed imaginary part
   !
       nord = nint(filspecs(3))
       fc   = filspecs(4)
       call lowPassButterworthFilter(nord,nf,fc,(nf1-1)*df,df,sigma,hlow)
   !
   !  multiply high and low pass and write to file
   !
       hfil(nf1:nf2) = hhigh*hlow
       deallocate(hlow,hhigh)
   !
   !  default case if filter type is not implemented or invalid
   !          
    case default
       call add(errmsg,2,'Unknown event filter type specified',myname)
       return
    end select
    end subroutine computeEventfilter
!
end module eventFilter

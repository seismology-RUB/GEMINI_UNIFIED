! =================================================================================
!  Compute transfer function of Butterworth filters at real and complex frequencies
! =================================================================================
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
!--------------------------------------------------------------------------------------
!  Routines for calculating the transfer function of Butterworth high and low-pass
!  filters evaluated at either real or complex frequencies.
!-------------------------------------------------------------------------------------
 module butterworthFilter
     use mathConstants
     implicit none
     interface highPassButterworthFilter
         module procedure highPassComplexOmegaButterworthFilter
         module procedure highPassRealOmegaButterworthFilter
     end interface
     interface lowPassButterworthFilter
         module procedure lowPassComplexOmegaButterworthFilter
         module procedure lowPassRealOmegaButterworthFilter
     end interface
 contains
!-----------------------------------------------------------------------
!  Calculate frequency response of HIGH-pass Butterworth filter
!  at selected complex frequencies with fixed imaginary part.
!  Returns all ones if nord = 0.
!  nord:    Filter order (in)
!  nf:      Number of frequencies where transfer function is evaluated (in)
!  fc:      Corner frequency (Hz) (in)
!  fmin:    Minimum frequency (Hz) (in)
!  df:      Frequency stepping (Hz) (in) (f = fmin+(i-1)*df-i*sigma/(2*pi))
!  sigma:   Constant negative imaginary part of ANGULAR frequency (Hz) (in)
!  h:       Allocatable array with complex transfer function values (out)
!------------------------------------------------------------------------
     subroutine highPassComplexOmegaButterworthFilter(nord,nf,fc,fmin,df,sigma,h)
     integer :: nord,nf
     double precision :: fc,fmin,df,sigma
     double complex :: zfn,zfna,ep
     double complex, dimension(:), allocatable :: h
     integer :: i,j
     allocate(h(nf))
     h = 1.d0
     zfna = (fmin-mc_cid*sigma/mc_two_pid)/fc            ! first normalized complex frequency
     do j = 1,nord
        zfn = zfna
        ep = exp(mc_cid*mc_pid*(dble(j)-.5)/dble(nord))  ! unit roots
        do i = 1,nf
           h(i) = h(i)*zfn/(zfn-ep)
           zfn = zfn+df/fc
        enddo
     enddo
     end subroutine highPassComplexOmegaButterworthFilter
!-----------------------------------------------------------------------
!  Calculate frequency response of LOW-pass Butterworth filter
!  at selected complex frequencies with fixed imaginary part.
!  Returns all ones if nord = 0.
!  nord:    Filter order (in)
!  nf:      Number of frequencies where transfer function is evaluated (in)
!  fc:      Corner frequency (Hz) (in)
!  fmin:    Minimum frequency (Hz) (in)
!  df:      Frequency stepping (Hz) (in) (f = fmin+(i-1)*df-i*sigma/(2*pi))
!  sigma:   Constant negative imaginary part of ANGULAR frequency (Hz) (in)
!  h:       Allocatable array with complex transfer function values (out)
!------------------------------------------------------------------------
     subroutine lowPassComplexOmegaButterworthFilter(nord,nf,fc,fmin,df,sigma,h)
     integer :: nord,nf
     double precision :: fc,fmin,df,sigma
     double complex :: zfn,zfna,ep
     double complex, dimension(:), allocatable :: h
     integer :: i,j
     allocate(h(nf))
     h = 1.d0
     zfna = (fmin-mc_cid*sigma/mc_two_pid)/fc            ! first normalized complex frequency
     do j = 1,nord
        zfn = zfna
        ep = exp(mc_cid*mc_pid*(dble(j)-.5)/dble(nord))  ! unit roots
        do i = 1,nf
           h(i) = h(i)*(-mc_cid)/(zfn-ep)
           zfn = zfn+df/fc
        enddo
     enddo
     end subroutine lowPassComplexOmegaButterworthFilter
!-----------------------------------------------------------------------
!  Calculate frequency response of HIGH-pass Butterworth filter
!  at selected real frequencies.
!  Returns all ones if nord = 0.
!  nord:    Filter order (in)
!  nf:      Number of frequencies where transfer function is evaluated (in)
!  fc:      Corner frequency (Hz) (in)
!  fmin:    Minimum frequency (Hz) (in)
!  df:      Frequency stepping (Hz) (in) (f = fmin+(i-1)*df-i*sigma/(2*pi))
!  h:       Allocatable array with complex transfer function values (out)
!------------------------------------------------------------------------
     subroutine highPassRealOmegaButterworthFilter(nord,nf,fc,fmin,df,h)
     integer :: nord,nf
     double precision :: fc,fmin,df
     double precision :: zfn
     double complex :: ep
     double complex, dimension(:), allocatable :: h
     integer :: i,j
     allocate(h(nf))
     h = 1.d0
     do j = 1,nord
        zfn = fmin/fc                                  ! first normalized real frequency
        ep = exp(mc_cid*mc_pid*(dble(j)-.5)/dble(nord))  ! unit roots
        do i = 1,nf
           h(i) = h(i)*zfn/(zfn-ep)
           zfn = zfn+df/fc
        enddo
     enddo
     end subroutine highPassRealOmegaButterworthFilter
!-----------------------------------------------------------------------
!  Calculate frequency response of HIGH-pass Butterworth filter
!  at selected real frequencies.
!  Returns all ones if nord = 0.
!  nord:    Filter order (in)
!  nf:      Number of frequencies where transfer function is evaluated (in)
!  fc:      Corner frequency (Hz) (in)
!  fmin:    Minimum frequency (Hz) (in)
!  df:      Frequency stepping (Hz) (in) (f = fmin+(i-1)*df-i*sigma/(2*pi))
!  h:       Allocatable array with complex transfer function values (out)
!------------------------------------------------------------------------
     subroutine lowPassRealOmegaButterworthFilter(nord,nf,fc,fmin,df,h)
     integer :: nord,nf
     double precision :: fc,fmin,df
     double precision :: zfn
     double complex :: ep
     double complex, dimension(:), allocatable :: h
     integer :: i,j
     allocate(h(nf))
     h = 1.d0
     do j = 1,nord
        zfn = fmin/fc                                      ! first normalized real frequency
        ep = exp(mc_cid*mc_pid*(dble(j)-.5)/dble(nord))   ! unit roots
        do i = 1,nf
           h(i) = h(i)*(-mc_cid)/(zfn-ep)
           zfn = zfn+df/fc
        enddo
     enddo
     end subroutine lowPassRealOmegaButterworthFilter
!
 end module butterworthFilter


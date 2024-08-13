!----------------------------------------------------------------------------
!   Copyright 2016 Wolfgang Friederich (Ruhr-Universitaet Bochum, Germany)
!
!   This file is part of GEMINI_UNIFIED version 1.0.
!
!   GEMINI_UNIFIED version 1.0 are free software:
!   you can redistribute it and/or modify
!   it under the terms of the GNU General Public License as published by
!   the Free Software Foundation, either version 2 of the License, or
!   (at your option) any later version.
!
!   GEMINI_UNIFIED version 1.0 are is distributed
!   in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU General Public License for more details.
!
!   You should have received a copy of the GNU General Public License
!   along with GEMINI_UNIFIED version 1.0.
!   If not, see <http://www.gnu.org/licenses/>.
!-------------------------------------------------------------
!> \brief Routines that operate with Fourier spectra
!-------------------------------------------------------------
 module fourierSpectrum
    use fourierTransform
    use timeSeries
    implicit none
    interface dealloc; module procedure deallocFourierSpectrum; end interface
    interface operator (.fmin.); module procedure fminFourierSpectrum; end interface
    interface operator (.fmax.); module procedure fmaxFourierSpectrum; end interface
    interface operator (.df.); module procedure dfFourierSpectrum; end interface
    interface operator (.nf.); module procedure nfFourierSpectrum; end interface
    interface operator (.nfa.); module procedure nfaFourierSpectrum; end interface
    interface operator (.nfb.); module procedure nfbFourierSpectrum; end interface
    interface operator (.fs.); module procedure complexFourierSpectrum; end interface
    interface operator (.amp.); module procedure amplitudeFourierSpectrum; end interface
    interface operator (.phase.); module procedure phaseFourierSpectrum; end interface
    interface operator (.vreal.); module procedure realFourierSpectrum; end interface
    interface operator (.vimg.); module procedure imgFourierSpectrum; end interface
    type fourier_spectrum
        integer :: nf1,nf2,nf                        !< Index of first and last frequency and number of frequencies
        double precision :: df                       !< Frequency interval
        logical :: link                              !< If true data is just a pointer
        double complex, dimension(:), pointer :: y => null()   !< Pointer to data (y(1) = SP(nf1), y(nf) = SP(nf2))
    end type
!
 contains
!----------------------------------------------------------------
!> \brief Create spectrum with a link to the data
!
    subroutine createLinkFourierSpectrum(this,nf1,nf2,df,s)
    type (fourier_spectrum) :: this
    integer :: nf1,nf2
    double precision :: df
    double complex, dimension(:), target :: s
!
    this%nf1 = nf1; this%nf2 = nf2; this%df = df; this%nf = nf2-nf1+1
    this%y => s
    this%link = .true.
    end subroutine createLinkFourierSpectrum
!------------------------------------------------------------------
!> \brief Create spectrum with a true copy of the data
!
    subroutine createFromDataFourierSpectrum(this,nf1,nf2,df,s)
    type (fourier_spectrum) :: this
    integer :: nf1,nf2
    double precision :: df
    double complex, dimension(:) :: s
!
    this%nf = nf2-nf1+1
    this%nf1 = nf1; this%nf2 = nf2; this%df = df
    allocate(this%y(this%nf))
    this%y(1:this%nf) = s(1:this%nf)
    this%link = .false.
    end subroutine createFromDataFourierSpectrum
!-----------------------------------------------------------------
!   Compute Fourier spectrum of a time series
!   ns: number of samples to be used for FFT (optional, power of 2)
!
    subroutine createFromTimeSeriesFourierSpectrum(this,ts,ns)
    type (fourier_spectrum) :: this
    type (time_series) :: ts
    integer, optional :: ns
    integer :: nplog,np2,nsamp
    double precision, dimension(:), allocatable :: x,y
!
! extend number of samples to next power of 2 or prescribed power of 2
!
    nsamp = .nsamp.ts
    if (present(ns)) then
       np2 = ns
       nplog = nint(log(real(ns))/log(2.d0))
    else
       nplog = ceiling(log(real(nsamp))/log(2.))
       np2 = 2**nplog
    endif
!
! copy data to double precision and append zeros
!
    allocate(x(np2),y(np2))
    x(1:nsamp) = .trace.ts
    x(nsamp+1:np2) = 0.d0
    y = 0.d0
!
! Fourier transform
! 
    call fastFourierTransform(x,y,nplog,-1)
!
! fill fourierTransform Structure
!
    this%nf1 = 1
    this%nf2 = np2/2
    this%nf = np2/2
    this%df = 1./(np2*.dt.ts)
    this%link = .false.
    allocate(this%y(np2/2))
    this%y(1:np2/2) = dcmplx(x(1:np2/2),y(1:np2/2)) * dble(.dt.ts)
!
    deallocate(x,y)
    end subroutine createFromTimeSeriesFourierSpectrum
!------------------------------------------------------------
!> \brief Deallocate fourier spectrum
!
    subroutine deallocFourierSpectrum(this)
    type (fourier_spectrum) :: this
    if (associated(this%y)) then
        if (this%link) then; nullify(this%y); else; deallocate(this%y); endif
    endif
    end subroutine deallocFourierSpectrum
!------------------------------------------------------------
!> \brief real part of spectrum
!
    function realFourierSpectrum(this) result(vreal)
    type (fourier_spectrum), intent(in) :: this
    real, dimension(:), pointer :: vreal
!
    allocate(vreal(this%nf))
    vreal = real(dreal(this%y))
    end function realFourierSpectrum
!------------------------------------------------------------
!> \brief imaginary part of spectrum
!
    function imgFourierSpectrum(this) result(vimg)
    type (fourier_spectrum), intent(in) :: this
    real, dimension(:), pointer :: vimg
!
    allocate(vimg(this%nf))
    vimg = real(dimag(this%y))
    end function imgFourierSpectrum
!------------------------------------------------------------
!> \brief Return pointer to double complex spectrum
!
    function complexFourierSpectrum(this) result(fs)
    type (fourier_spectrum), intent(in) :: this
    double complex, dimension(:), pointer :: fs
!
    fs => this%y
    end function complexFourierSpectrum
!------------------------------------------------------------
!> \brief Amplitude spectrum
!
    function amplitudeFourierSpectrum(this) result(amp)
    type (fourier_spectrum), intent(in) :: this
    real, dimension(:), pointer :: amp
!
    allocate(amp(this%nf))
    amp = sqrt(zabs(this%y))
    end function amplitudeFourierSpectrum
!------------------------------------------------------------
!> \brief Phase spectrum
!
    function phaseFourierSpectrum(this) result(phase)
    type (fourier_spectrum), intent(in) :: this
    real, dimension(:), pointer :: phase
!
    allocate(phase(this%nf))
    phase = real(datan2(dimag(this%y),dreal(this%y)))
    end function phaseFourierSpectrum
!------------------------------------------------------------
!> \brief In place phase shift of spectrum by tau (*exp(-i*om*tau))
!
    subroutine phaseShiftFourierSpectrum(this,tau)
    type (fourier_spectrum) :: this
    double precision :: tau
    integer :: i,jf
    double precision :: om
!
    do i = 1,this%nf
       jf = i-1+this%nf1
       om = 2.d0*mc_pid*(jf-1)*this%df
       this%y(i) = this%y(i)*exp(-mc_cid*om*tau)
    enddo
    end subroutine phaseShiftFourierSpectrum
!------------------------------------------------------------
!> \brief In place scaling of spectrum by real factor
!
    subroutine scalarMultiplyFourierSpectrum(this,mag)
    type (fourier_spectrum) :: this
    double precision :: mag
!
    this%y = this%y*mag
    end subroutine scalarMultiplyFourierSpectrum
!------------------------------------------------------------
!> \brief Get df
!
    double precision function dfFourierSpectrum(this)
    type (fourier_spectrum), intent(in) :: this
    dfFourierSpectrum = this%df
    end function dfFourierSpectrum
!------------------------------------------------------------
!> \brief Get nf
!
    integer function nfFourierSpectrum(this)
    type (fourier_spectrum), intent(in) :: this
    nfFourierSpectrum = this%nf
    end function nfFourierSpectrum
!------------------------------------------------------------
!> \brief Get nf1
!
    integer function nfaFourierSpectrum(this)
    type (fourier_spectrum), intent(in) :: this
    nfaFourierSpectrum = this%nf1
    end function nfaFourierSpectrum
!------------------------------------------------------------
!> \brief Get nf2
!
    integer function nfbFourierSpectrum(this)
    type (fourier_spectrum), intent(in) :: this
    nfbFourierSpectrum = this%nf2
    end function nfbFourierSpectrum
!------------------------------------------------------------
!> \brief Get fmin
!
    real function fminFourierSpectrum(this)
    type (fourier_spectrum), intent(in) :: this
    fminFourierSpectrum = (this%nf1-1)*this%df
    end function fminFourierSpectrum
!------------------------------------------------------------
!> \brief Get fmax
!
    real function fmaxFourierSpectrum(this)
    type (fourier_spectrum), intent(in) :: this
    fmaxFourierSpectrum = (this%nf2-1)*this%df
    end function fmaxFourierSpectrum
!
 end module

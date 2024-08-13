! ===============================================================================
!  Main program for plotting Green FK spectra
! ===============================================================================
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
!-----------------------------------------------------------------------------
!   Plot frequency wavenumber spectra
!-----------------------------------------------------------------------------
program plotGreenFKSpectra
       use hdf5
    use mathConstants
    use argumentParser
    use greenFKSpectra
    use hdfWrapper
    use errorMessage
    use string
    use pgPlotImage
    use pgPlotWindow
    implicit none
    type (argument_parser) :: ap
    type (error_message) :: errmsg
    type (green_fk_spectra) :: gfk
    real, dimension(:), pointer :: fkr,fki
    real, dimension(:,:), pointer :: gfkspr,gfkspi
    real, dimension(:,:), allocatable :: fkplot
    real :: bg,fkmax,fmin,fmax,wmax,width,aspect
    integer(hid_t) :: fid
    integer :: jr,cnt,nwa,nwe,jf,isp,jsp,nf,nsp,ierr
    character (len=max_length_string) :: gfkmeta,gfkdata,ylab
    type (pgplot_image) :: pgfk
    type (pgplot_window) :: pgwin,pgps
    real :: x1,x2,y1,y2
    character (len=1) :: ch,cjsp
    character (len=18) :: myname = 'plotGreenFKSpectra'
!-----------------------------------------------------------------
    call init(ap,myname,'Plot Green frequency wavenumber spectra')
    call addPosarg(ap,'gfkmeta','sval','file name of gfk output or meta file')
    call addOption(ap,'-gfkdata',.true.,'GFK data file','sval','None')
    call addOption(ap,'-jr',.true.,'node index where fk-spectrum is evaluated','ival','1')
    call addOption(ap,'-bg',.true.,'background level for logarithmic plot','rval','-4.0')
    call addOption(ap,'-w',.true.,'plot window width','rval','8.')
    call addOption(ap,'-a',.true.,'plot window aspect ratio','rval','0.65')
    call parse(ap)
    width = ap.rval.'-w'
    aspect = ap.rval.'-a'
    bg = ap.rval.'-bg'
    jr = ap.ival.'-jr'
    gfkdata = ap.sval.'-gfkdata'
    gfkmeta = ap.sval.'gfkmeta'
    if (.level.(.errmsg.ap) == 2) then; call print(.errmsg.ap); call usage(ap); stop; endif
    call document(ap); call dealloc(ap)
!
    call new(errmsg,myname)
!-------------------------------------------------------
!  open frequency-wavenumber spectrum file
!  and read in group tree and header info
!
    call openEnvironmentHDFWrapper(errmsg)
    if (.level.errmsg == 2) goto 10
    call openFileRoHDFWrapper(gfkmeta,fid,errmsg)
    if (.level.errmsg == 2) goto 10
    call readMetaGreenFKSpectra(gfk,fid,'identity',errmsg)
    if (.level.errmsg == 2) goto 10
    call readMetaGreenFKSpectra(gfk,fid,'reals',errmsg)
    if (.level.errmsg == 2) goto 10
    call readMetaGreenFKSpectra(gfk,fid,'integers',errmsg)
    if (.level.errmsg == 2) goto 10
    call readMetaGreenFKSpectra(gfk,fid,'numberOfWavenumbersForFrequency',errmsg)
    if (.level.errmsg == 2) goto 10
    if (gfkdata.equal.'None') then
       call readMetaGreenFKSpectra(gfk,fid,'dataSpecificIntegers',errmsg)
       if (.level.errmsg == 2) goto 10
    else
       call h5fclose_f(fid,ierr)
       call openFileRoHDFWrapper(gfkdata,fid,errmsg)
       if (.not. checkIdentityGreenFKSpectra(gfk,fid,errmsg)) then
          print *,'ERROR: identity check for meta and data part of GFK file failed!'
          stop
       endif
       call readMetaGreenFKSpectra(gfk,fid,'dataSpecificIntegers',errmsg)
       if (.level.errmsg == 2) goto 10
    endif
    call printMembersGreenFKSpectra(gfk)
!
!  get frequency and wavenumber range and
!  allocate space for array to be plotted
!
    fmin = (gfk%nf1-1)*gfk%df
    fmax = (gfk%nf2-1)*gfk%df
    if (gfk%global == 1) then
       wmax = gfk%nwnmax-1
       ylab = "Harmonic degree"
    else
       wmax = (gfk%nwnmax-1)*gfk%dwn
       ylab = "Wavenumber"
    endif
    nf = gfk%nf2-gfk%nf1+1
    nsp = numberGreenFKSpectra(gfk)
    allocate(fkplot(nf,gfk%nwnmax))
!
!  read fk-spectra for node jr
!
    call readDataNodeGreenFKSpectra(gfk,fid,jr,gfkspr,gfkspi,errmsg)
    if (.level.errmsg == 2) goto 10
!
!  print out which components and jumps are available
!
    call printContentsGreenFKSpectra(gfk)
    call h5fclose_f(fid,ierr)
    call h5close_f(ierr)
!
!  open plot window
!
     call new(pgwin,width = width,aspect=aspect)
!
!  start plotting with DSV-index isp = 1 and jump-index jsp = 1
!  should always be available
!  unwrap gfkspectra    
!
    isp = 1; jsp = 1
1   cnt = positionDataGreenFKSpectra(gfk,isp,jsp)
    fkr => gfkspr(:,cnt)
    fki => gfkspi(:,cnt)
    fkmax = maxval(cabs(cmplx(fkr,fki)))
    print *,'fkmax = ',fkmax
    write(cjsp,'(i1)') jsp
!
!  initialize plot array with background value
!
4   fkplot = bg
    nwe = 0
    do jf = gfk%nf1,gfk%nf2
       nwa = nwe+1
       nwe = nwe+gfk%nwn(jf)
       fkplot(jf-gfk%nf1+1,1:gfk%nwn(jf)) = &
            max(log10(cabs(cmplx(fkr(nwa:nwe),fki(nwa:nwe)))/fkmax),bg)
    enddo
!
!  prepare plot image
!
    call new(pgfk,fkplot,x1 = fmin,xn = fmax,y1 = 0.0,yn = wmax,xlab = 'Frequency',ylab = ylab, &
          title = 'FK-Spectrum for DSV = '+names_green_fk_spectra(isp)+' and jump = '+cjsp)
!
!  plot
!
 2  call display(pgfk)
!
!  allow interacive action of user
!  'o':           go back to original extrema
!  's[isp][jsp]': select a new spectrum for plotting
!  'b':           change background value
!  'h'            produce a hardcopy PS file with name plotfk.eps
!
    print *,'Available options: A D X z Z l r u d o s[isp][jsp] b[n] h'
 3  call pgcurs(x1,y1,ch)
    select case (ch)
    case ('A','D','X','z','Z','l','r','u','d')
        call zoomBindingsPgPlotSelect(x1,y1,ch,x2,y2,.xmin.pgfk,.xmax.pgfk,.ymin.pgfk,.ymax.pgfk)
        call setExtremaPgPlotImage(pgfk,x1,x2,y1,y2)
        goto 2
    case ('o'); call setOriginalExtremaPgPlotImage(pgfk); goto 2
    case ('s'); 
       call pgcurs(x1,y1,ch); read(ch,'(i1)') isp; call pgcurs(x1,y1,ch); read(ch,'(i1)') jsp; call dealloc(pgfk)
       if (isp > 9 .or. jsp > 4 .or. isp < 1 .or. jsp < 1) then
          print *,'Invalid values for isp,jsp. Choose from 1 <= isp <= 9 and 1 <= jsp <= 4'
          goto 3
       endif
       if (.not. existGreenFKSpectra(gfk,isp,jsp)) then
          print *,'Spectrum DSV = ',names_green_fk_spectra(isp),' Jump = ',jsp,' is not available. Choose a different one'
          goto 3
       endif
       goto 1
    case ('b'); call pgcurs(x1,y1,ch); read(ch,*) bg; bg = -bg; call dealloc(pgfk); goto 4
    case ('h'); call new(pgps,plotfile = 'plotfk.eps'); call display(pgfk); call select_win(pgwin); goto 2
    end select
!
!  clean up
!
    deallocate(fkplot,gfkspr,gfkspi)
    call dealloc(pgfk)
    call dealloc(gfk)
    
10  if (.level.errmsg == 2) then
       call print(errmsg)
    endif
!    
 end program

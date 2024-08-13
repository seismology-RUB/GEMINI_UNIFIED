! ===============================================================================
!  Main program for visualizing a seismogram gather
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
!--------------------------------------------------------------------------------
!  Visualizes a seismogram gather and allows user interaction with the plot
!  such as scrolling, box, offest and time zooming.
!  Allows to overlay several seismogram gathers into one plot to
!  compare seismograms. Works as well for gathers with only one seismogram.
!--------------------------------------------------------------------------------
program plotAsciiSeismogramGather
    use interactiveBindingsPlotGather
    use argumentParser
    use pgPlotWindow
    use asciiSynseisPlotGather
    use string
    use errorMessage
!
    implicit none
    type (error_message) :: errmsg
    type (argument_parser) :: ap
    type (pgplot_window) :: pgwin,pgps
    type (ascii_synseis_plot_gather), dimension(:), allocatable :: gather
    character (len=max_length_string) :: filelist,pickfile
    character(len=max_length_string), dimension(:), pointer :: datafile
    character (len=22) :: myname = 'plotSeismogramGather'
    character (len=1) :: otherkey
    real :: width,aspect
    integer :: shift,i,nsec
    logical :: leave,txmode,hardcopy
!-------------------------------------------------------------------------------
    call init(ap,myname,'Plot a synthetic seismogram gather')
    call addPosarg(ap,'gather_files','sval','list of gather files')
    call addOption(ap,'-tx',.false.,'time axis up')
    call addOption(ap,'-w',.true.,'plot window width','rval','10.')
    call addOption(ap,'-a',.true.,'plot window aspect ratio','rval','0.65')
    call addOption(ap,'-pickfile',.true.,'name of file for picks','sval','unspec.pck')
    call parse(ap)
    if (.level.(.errmsg.ap) == 2) then; call usage(ap); call print(.errmsg.ap); stop; endif
    filelist = ap.sval.'gather_files'
    txmode = ap.optset.'-tx'
    width = ap.rval.'-w'
    aspect = ap.rval.'-a'
    pickfile = ap.sval.'-pickfile'
    call document(ap); call dealloc(ap)
!------------------------------------------------------------------------------
!  read in files
!
    datafile => getWordsString(filelist,' ',errmsg)
    if (.level.errmsg == 2) goto 10
    nsec = size(datafile)
    allocate(gather(nsec))
    do i=1,nsec
        call gather(i)%create(1,datafile(i),errmsg,ci = i,txmode = txmode,pickfile = pickfile)
        if (.level.errmsg == 2) goto 10
        print *,'Reading gather file: ',trim(datafile(i))
    enddo
!
!  open plot window and display gathers
!
    call new(pgwin,width = width, aspect = aspect)
!
!  display seismogram gathers and enable user interaction
!
    call docuInteractiveBindingsPlotGather()
    shift = 0; leave = .false.
    do while (.not. leave)
       call gather(1)%display()
       do i=2,nsec
          call gather(i)%overlay(shift)
       enddo
       call interactWithPlotGather(gather,leave,shift,hardcopy,otherkey)
       if (leave) exit
     !
     !  produce postscript output
     !
       if (hardcopy) then
          call new(pgps,plotfile = 'plotgather.ps',width=width, aspect=aspect)
          call gather(1)%display()
          do i=2,nsec
             call gather(i)%overlay(shift)
          enddo
          call dealloc(pgps); call select_win(pgwin)
       endif
    enddo
!
10  if (.level.errmsg == 2) then
       call print(errmsg)
    endif
!
end program plotAsciiSeismogramGather

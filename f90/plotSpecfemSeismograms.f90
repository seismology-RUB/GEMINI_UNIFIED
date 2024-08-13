! ===============================================================================
!  Main program for visualizing SPECEFEM seismograms
! ===============================================================================
!----------------------------------------------------------------------------
!   Copyright 2021 Wolfgang Friederich (Ruhr-Universitaet Bochum, Germany)
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
program plotSpecfemSeismograms
    use interactiveBindingsPlotGather
    use argumentParser
    use pgPlotWindow
    use specfemSynseisPlotGather
    use string
    use errorMessage
!
    implicit none
    type (error_message) :: errmsg
    type (argument_parser) :: ap
    type (pgplot_window) :: pgwin,pgps
    type (specfem_synseis_plot_gather), dimension(:), allocatable :: gather
    character (len=max_length_string) :: path,pickfile,stationfile,bandcode,seistype,evid
    character (len=22) :: myname = 'plotSpecfemSeismograms'
    character (len=1) :: otherkey
    real :: width,aspect
    integer :: shift
    logical :: leave,txmode,hardcopy,prep_evid
!-------------------------------------------------------------------------------
    call init(ap,myname,'Plot SPECFEM synthetic seismograms')
    call addPosarg(ap,'path','sval','path to Specfem seismograms incl slash')
    call addPosarg(ap,'stationfile','sval','File with list of stations')
    call addOption(ap,'-bandcode',.true.,'bandcode including component (e.g. BXX)','sval','BXZ')
    call addOption(ap,'-seistype',.true.,'type seismogram (e.g. d,v,a)','sval','d')
    call addOption(ap,'-evid',.true.,'prepend eventid to file name','sval','191207_000000')
    call addOption(ap,'-tx',.false.,'time axis up')
    call addOption(ap,'-w',.true.,'plot window width','rval','10.')
    call addOption(ap,'-a',.true.,'plot window aspect ratio','rval','0.65')
    call addOption(ap,'-pickfile',.true.,'name of file for picks','sval','unspec.pck')
    call parse(ap)
    if (.level.(.errmsg.ap) == 2) then; call usage(ap); call print(.errmsg.ap); stop; endif
    path = ap.sval.'path'
    stationfile = ap.sval.'stationfile'
    bandcode = ap.sval.'-bandcode'
    seistype = ap.sval.'-seistype'
    evid = ap.sval.'-evid'
    prep_evid = ap.optset.'-evid'
    txmode = ap.optset.'-tx'
    width = ap.rval.'-w'
    aspect = ap.rval.'-a'
    pickfile = ap.sval.'-pickfile'
    call document(ap); call dealloc(ap)
!------------------------------------------------------------------------------
!  read data
!
    allocate(gather(1))
    call gather(1)%create(path,stationfile,evid,prep_evid,bandcode,seistype,errmsg,&
                       ci = 1,txmode = txmode,pickfile = pickfile)
    if (.level.errmsg == 2) goto 10
!
!  open plot window and display gather
!
    call new(pgwin,width = width, aspect = aspect)
!
!  display seismogram gathers and enable user interaction
!
    call docuInteractiveBindingsPlotGather()
    shift = 0; leave = .false.
    do while (.not. leave)
       call gather(1)%display()
       call interactWithPlotGather(gather,leave,shift,hardcopy,otherkey)
       if (leave) exit
       if (hardcopy) then
          call new(pgps,plotfile = 'plotgather.ps')
          call gather(1)%display()
          call dealloc(pgps); call select_win(pgwin)
       endif
    enddo
!
10  if (.level.errmsg == 2) then
       call print(errmsg)
    endif
!
end program plotSpecfemSeismograms

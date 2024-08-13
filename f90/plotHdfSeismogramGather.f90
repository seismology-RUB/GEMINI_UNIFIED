! ===============================================================================
!  Main program for visualizing a HDF seismogram gather
! ===============================================================================
!----------------------------------------------------------------------------
!   Copyright 2019 Wolfgang Friederich (Ruhr-Universitaet Bochum, Germany)
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
program plotHdfSeismogramGather
    use hdf5
    use interactiveBindingsPlotGather
    use argumentParser
    use pgPlotWindow
    use hdfSynseisPlotGather
    use synseisHDF
    use string
    use errorMessage
!
    implicit none
    type (error_message) :: errmsg
    type (argument_parser) :: ap
    type (pgplot_window) :: pgwin,pgps
    type (hdf_synseis_plot_gather), dimension(:), allocatable :: gather
    character (len=max_length_string) :: hdflist,pickfile
    character (len=max_length_string), dimension(:), pointer :: datafile
    character (len=23) :: myname = 'plotHdfSeismogramGather'
    integer(hid_t), dimension(:), allocatable :: fid
    character(len=char_len_evid) :: evid
    character(len=1) :: otherkey,ch,comp
    real :: width,aspect,xa,ya,tend,tbeg
    integer :: shift,ierr,i,nsec
    logical :: leave,txmode,hardcopy,tzero
!-------------------------------------------------------------------------------
    call init(ap,myname,'Plot a synthetic seismogram gather')
    call addPosarg(ap,'hdffiles','sval','List of synseisHDF-files')
    call addOption(ap,'-c',.true.,'Component','sval','Z')
    call addOption(ap,'-tbeg',.true.,'start of time window','rval','-1.0')
    call addOption(ap,'-tend',.true.,'end of time window','rval','-1.0')
    call addOption(ap,'-tx',.false.,'time axis up')
    call addOption(ap,'-w',.true.,'plot window width','rval','10.')
    call addOption(ap,'-a',.true.,'plot window aspect ratio','rval','0.65')
    call addOption(ap,'-pickfile',.true.,'name of file for picks','sval','unspec.pck')
    call addOption(ap,'-tzero',.false.,'time axis is zero at trace start')
    call parse(ap)
    if (.level.(.errmsg.ap) == 2) then; call print(.errmsg.ap); call usage(ap); stop; endif
    hdflist = ap.sval.'hdffiles'
    comp = ap.sval.'-c'
    txmode = ap.optset.'-tx'
    width = ap.rval.'-w'
    aspect = ap.rval.'-a'
    pickfile = ap.sval.'-pickfile'
    tbeg = ap.rval.'-tbeg'
    tend = ap.rval.'-tend'
    tzero = ap.optset.'-tzero'
    call document(ap); call dealloc(ap)
    
    datafile => getWordsString(hdflist,' ',errmsg)
    if (.level.errmsg == 2) then; call print(errmsg); stop; endif
    nsec = size(datafile)
    allocate(gather(nsec),fid(nsec))
    call openEnvironmentHDFWrapper(errmsg)
    if (.level.errmsg == 2) then; call print(errmsg); stop; endif
!------------------------------------------------------------------------------
!  read in files
!
    do i = 1,nsec
       call openFileRoHDFWrapper(trim(datafile(i)),fid(i),errmsg)
       if (.level.errmsg == 2) then; call print(errmsg); stop; endif
       call gather(i)%create(fid(i),evid,comp,errmsg,tbeg = tbeg,tend = tend,ci = i,&
                             txmode = txmode,pickfile = pickfile,tzero=tzero)
       if (.level.errmsg == 2) then; call print(errmsg); stop; endif
       print *,'Reading gather file: ',trim(datafile(i)),', Event ID: ',evid,', Component: ',comp
    enddo   
!
!  open plot window and display gathers
!
    call new(pgwin,width = width, aspect = aspect)
!
!  display seismogram gather and enable user interaction
!
    call docuInteractiveBindingsPlotGather()
    shift = 0; leave = .false.
    do while (.not. leave)
      call gather(1)%display()
      do i = 2,nsec
         call gather(i)%overlay(shift)
      enddo
      call interactWithPlotGather(gather,leave,shift,hardcopy,otherkey)
      if (leave) exit
      if (otherkey == 'c') then      ! choose another component
         call pgcurs(xa,ya,ch)
         do i = 1,nsec
            call gather(i)%dealloc()
            call gather(i)%create(fid(i),evid,ch,errmsg,tbeg = tbeg,tend = tend,ci = i,&
                                  txmode = txmode,pickfile = pickfile,tzero=tzero)
            if (.level.errmsg == 2) then; call print(errmsg); stop; endif
         enddo
      endif
   !
   !  produce postscript output
   !
      if (hardcopy) then
         call new(pgps,plotfile = 'plotgather.ps')
         call gather(1)%display()
         do i = 2,nsec
            call gather(i)%overlay(shift)
         enddo
         call dealloc(pgps); call select_win(pgwin)
      endif
   enddo
!
!  clean up
!
   do i = 1,nsec
      call gather(i)%dealloc()
      call h5fclose_f(fid(i),ierr)
      if (ierr < 0)  then
         call add(errmsg,2,'Cannot close HDF file',myname); call print(errmsg)
      endif
   enddo
   call h5close_f(ierr)
   deallocate(datafile,fid)
   if (ierr < 0)  then
      call add(errmsg,2,'Cannot close HDF environment',myname); call print(errmsg)
   endif
end program plotHdfSeismogramGather

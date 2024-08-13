! ===============================================================================
!  Specific module implementing an hdfSynseis plot gather
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
!----------------------------------------------------------------------------
!  Defines extended type "hdf_synseis_plot_gather" for plotting seismogram
!  gathers in hdfSynseis format. Provides a format specific "create" routine
!  reading data and giving values to base_plot_gather type components. Could also
!  overwrite procedures of base_plot_gather type.
!-----------------------------------------------------------------------
module hdfSynseisPlotGather
       use hdf5
    use constants
    use synseisHDF
    use seismicEvent
    use seismicStation
    use errorMessage
    use basePlotGather
    implicit none
    type, extends(base_plot_gather) :: hdf_synseis_plot_gather
    contains
        procedure :: create => createHdfSynseisPlotGather
    end type hdf_synseis_plot_gather
contains
!--------------------------------------------------------------------------
!  create a pgplot_section object from a SYNSEIS HDF seismogram gather file
!  fid:  file id (assumes file is already open)
!
! optional arguments:
!  lw:            line width (def = 2)
!  ci:            colour index (def = 1)
!  ls:            line style (def = 1)
!  ch:            character height (def = 1.5)
!  xlab:          label at x-axis (def = empty)
!  ylab:          label at y-axis (def = empty)
!-----------------------------------------------------------------------------
    subroutine createHdfSynseisPlotGather(this,fid,evid,comp,errmsg,tbeg,tend,lw,ci,ls,ch,xlab,ylab,title,txmode,pickfile,tzero)
    class (hdf_synseis_plot_gather) :: this
    integer (hid_t) :: fid
    character (len=char_len_evid) :: evid
    character (len=1) :: comp
    type (error_message) :: errmsg
    integer, optional :: ls,lw,ci
    real, optional :: ch,tend,tbeg
    character (len=*), optional :: xlab,ylab,title,pickfile
    logical, optional :: txmode,tzero
    type (seismic_event) :: event
    type (seismic_station) :: station
    character (len=char_len_netcode+char_len_sta+1), dimension(:), allocatable :: stdsetname
    character(len = :), allocatable :: components
    integer(hid_t) :: evgrid
    integer :: ierr,nsta,j,poscomp,nd,nend,nbeg
    double precision :: tanf,dt
    real :: ta,te
    real, dimension(:,:), pointer :: urs
    character (len=31) :: myname = 'createHdfSynseisPlotGather'
!
    if (present(tend)) then; te = tend; else; te = -1.0; endif
    if (present(tbeg)) then; ta = tbeg; else; ta = -1.0; endif
!
    call openEventSynseisHDF(fid,evid,evgrid,event,errmsg)
    if (.level.errmsg == 2) then; call print(errmsg); stop; endif
!
!  Iterate over stations in file
!
    call  getDatasetNamesSynseisHDF(fid,stdsetname,nsta,errmsg)
    if (.level.errmsg == 2) return
    call this%allocate(nsta)
    nd = 0   
    do j = 1,nsta
       call readStationSynseisHDF(fid,trim(stdsetname(j)),station,components,tanf,dt,urs,errmsg)
       if (.level.errmsg == 2) return
       poscomp = index(components,comp)
       if (poscomp == 0) cycle                       ! desired component not present 
       nd = nd+1
       if (te > 0.0) then; nend = nint(te/dt); else; nend = size(urs,1); endif
       if (ta > 0.0) then; nbeg = nint(ta/dt); else; nbeg = 1; endif
       call createDPFromDataTimeSeries(this%seis(nd),nend-nbeg+1,tanf,dt,urs(nbeg:nend,poscomp))
       this%tag(nd) = stdsetname(j)+' '+comp
       deallocate(urs,components)
    enddo
    if (nd == 0) then
       call add(errmsg,2,'No seismograms found, component wrong?',myname)
       return
    endif
    call h5gclose_f(evgrid,ierr)
    if (ierr < 0) goto 1
!
!  remaining setup
!
    this%nseis = nd
    call this%setOffsetTraceIndices()
    call this%setup(lw,ci,ls,ch,xlab,ylab,title,txmode,pickfile,tzero)
    call dealloc(event)
    return
 1  call add(errmsg,2,'Failed in '+myname,myname)
    end subroutine createHdfSynseisPlotGather
!
end module hdfSynseisPlotGather

! ===============================================================================
!  Specific module implementing an specfemSynseis plot gather
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
!----------------------------------------------------------------------------
!  Defines extended type "ascii_synseis_plot_gather" for plotting seismogram
!  gathers in asciiSynseis format. Provides a format specific "create" routine
!  reading data and giving values to plot_gather type components. Could also
!  overwrite procedures of base_plot_gather type.
!-----------------------------------------------------------------------
module specfemSynseisPlotGather
    use asciiDataIO
    use readEventStationFile
    use errorMessage
    use basePlotGather
    implicit none
    type, extends(base_plot_gather) :: specfem_synseis_plot_gather
    contains
        procedure :: create => createSpecfemSynseisPlotGather
    end type specfem_synseis_plot_gather
contains
!-----------------------------------------------------------------------
!  create a pgplot_section object from SPECFEM synthetic
!  lu:           file unit number
!  path:         path to seismograms
!  stationfile:  list of stations
!
! optional arguments:
!  lw:            line width (def = 2)
!  ci:            colour index (def = 1)
!  ls:            line style (def = 1)
!  ch:            character height (def = 1.5)
!  xlab:          label at x-axis (def = empty)
!  ylab:          label at y-axis (def = empty)
!
  subroutine createSpecfemSynseisPlotGather(this,path,stationfile,evid,prep_evid,bandcode,seistype,errmsg,&
                                            lw,ci,ls,ch,xlab,ylab,title,txmode,pickfile)
     class (specfem_synseis_plot_gather) :: this
     character (len=*) path,stationfile,evid,bandcode,seistype
     logical :: prep_evid
     integer, optional :: ls,lw,ci
     character (len=*), optional :: xlab,ylab,title,pickfile
     type (error_message) :: errmsg
     real, optional :: ch
     logical, optional :: txmode
     type (seismic_network) :: station_list
     type (seismic_station) :: station
     integer :: nd,is
     real, dimension(:,:), pointer :: p
     character (len=max_length_string) :: staname,filename
!
!  red station list
!
    call createStationListFromStationFile(stationfile,1,'ASKI_stations',station_list,errmsg)
    if (.level.errmsg == 2) return
    nd = .nstat.station_list
   !
   !  read gather file trace by trace
   !
     allocate(this%seis(nd),this%offset(nd),this%tag(nd),this%smax(nd),this%somax(nd)) 
     do while (nextStationSeismicNetwork(station_list,station,is))
        staname = .staname.station
        if (prep_evid) then
           filename = trim(path)//trim(evid)//'.'//trim(.netcode.station)//'.'//trim(staname)//&
                      '.'//trim(bandcode)//'.sem'//trim(seistype)
        else
           filename = trim(path)//trim(.netcode.station)//'.'//trim(staname)//&
                      '.'//trim(bandcode)//'.sem'//trim(seistype)
        endif
        print *,'Reading ',trim(filename)
        call readRealMatrixAsciiDataIO(filename,1,p,errmsg,2)
        if (.level.errmsg == 2) return
        call createDPFromDataTimeSeries(this%seis(is),size(p,1),dble(p(1,1)),dble(p(2,1)-p(1,1)),p(:,2))
        this%tag(is) = (.staname.station)+'_'+trim(bandcode)
        deallocate(p)
     enddo
   !
   !  remaining setup
   !
     this%nseis = nd
     call this%setOffsetTraceIndices()
     call this%setup(lw,ci,ls,ch,xlab,ylab,title,txmode,pickfile)
  end subroutine createSpecfemSynseisPlotGather
!
end module specfemSynseisPlotGather

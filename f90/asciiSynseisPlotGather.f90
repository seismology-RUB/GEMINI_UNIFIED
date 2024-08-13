! ===============================================================================
!  Specific module implementing an asciiSynseis plot gather
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
!  Defines extended type "ascii_synseis_plot_gather" for plotting seismogram
!  gathers in asciiSynseis format. Provides a format specific "create" routine
!  reading data and giving values to plot_gather type components. Could also
!  overwrite procedures of base_plot_gather type.
!-----------------------------------------------------------------------
module asciiSynseisPlotGather
    use asciiSynseisIO
    use errorMessage
    use basePlotGather
    implicit none
    type, extends(base_plot_gather) :: ascii_synseis_plot_gather
    contains
        procedure :: create => createAsciiSynseisPlotGather
    end type ascii_synseis_plot_gather
contains
!-----------------------------------------------------------------------
!  create a pgplot_section object from a SYNSEIS ASCII seismogram gather
!  lu:        file unit number
!  filename:  name of file
!
! optional arguments:
!  lw:            line width (def = 2)
!  ci:            colour index (def = 1)
!  ls:            line style (def = 1)
!  ch:            character height (def = 1.5)
!  xlab:          label at x-axis (def = empty)
!  ylab:          label at y-axis (def = empty)
!
  subroutine createAsciiSynseisPlotGather(this,lu,filename,errmsg,lw,ci,ls,ch,xlab,ylab,title,txmode,pickfile)
     class (ascii_synseis_plot_gather) :: this
     integer :: lu
     character (len=*) filename
     integer, optional :: ls,lw,ci
     character (len=*), optional :: xlab,ylab,title,pickfile
     type (error_message) :: errmsg
     real, optional :: ch
     logical, optional :: txmode
     integer :: ios,nd,i,nsamp
     real :: dt
     double precision :: tanfdp
     logical :: do_open,do_close
     real, dimension(:), allocatable :: urs
     character (len=31) :: myname = 'createAsciiSynseisPlotGather'
     character (len=2) :: comp
     character (len=max_length_string) :: staname
   !
   !  first scan file to find out number of seismograms
   !
     nd = getTraceCountAsciiSynseisIO(lu,filename,ios)
     if (ios > 0) then
        call add(errmsg,2,'Problems accessing data file: '+filename,myname)
        return
     endif
   !
   !  read gather file trace by trace
   !
     allocate(this%seis(nd),this%offset(nd),this%tag(nd),this%smax(nd),this%somax(nd)) 
     do_open = .true.; do_close = .false.
     do i = 1,nd
        call readAsciiSynseisIO(lu,filename,nsamp,tanfdp,dt,staname,comp,urs,ios,op = do_open,cl = do_close)
        call createDPFromDataTimeSeries(this%seis(i),nsamp,tanfdp,dble(dt),urs)
        this%tag(i) = staname+'_'+comp
        deallocate(urs)
        do_open = .false.
     enddo
     close(lu)
   !
   !  remaining setup
   !
     this%nseis = nd
     call this%setOffsetTraceIndices()
     call this%setup(lw,ci,ls,ch,xlab,ylab,title,txmode,pickfile)
  end subroutine createAsciiSynseisPlotGather
!
end module asciiSynseisPlotGather

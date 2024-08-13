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
!--------------------------------------------------------
!  module to describe a GSE data block
!--------------------------------------------------------
 module gseDatablock
    implicit none
    interface new
        module procedure readGSEDatablock
    end interface
    interface dealloc
        module procedure deallocGSEDatablock
    end interface dealloc
    interface trace
        module procedure traceGSEDatablock
    end interface trace
    interface operator (.nsamp.); module procedure getNGSEDatablock; end interface
    interface operator (.dt.); module procedure getDtGSEDatablock; end interface
    interface operator (.station.); module procedure getStationGSEDatablock; end interface
    interface operator (.channel.); module procedure getChannelGSEDatablock; end interface
    interface operator (.tanf.); module procedure getTanfGSEDatablock; end interface
    interface operator (.year.); module procedure getYearGSEDatablock; end interface
    interface operator (.doy.); module procedure getDoyGSEDatablock; end interface
    type gse_datablock
        private
        character :: wid*132                 ! Identification string 
        integer :: n                         ! number of data 
        real :: dt                           ! sampling interval 
        real :: tstart                       ! t-value of first data point after midnight
        integer, dimension(:), pointer :: y           ! pointer to integer data
        real, dimension(:), pointer :: ry => null()    ! pointer to real data
        logical :: link                      ! .true. if no separate memory was allocated for y
    end type gse_datablock
!
 contains
!--------------------------------------------------------------
!  constructor
!
    subroutine readGSEDatablock(this,lu,filename,last)
    type (gse_datablock) :: this
    integer :: lu,ierr
    logical :: last
    character (len=*) :: filename
!
    open(lu,file = filename,err = 99, status = 'old')
    call gse_RWid2(lu,this%wid,this%n,this%tstart,this%dt,ierr,last)
    if(ierr /= 0) then
        print *,'Error while reading WID2 line: '
        stop
    endif
    allocate(this%y(this%n))
    this%link = .false.
    call gse_RData(lu,this%n,this%y,ierr)
    if(ierr /= 0) then
        print *,'Error while reading data block: ', this%wid(1:68)
        stop
    endif
    if(last) close(lu)
    return
 99    print *,'<readGSEDataBlock>: Error opening ',filename
    end subroutine readGSEDatablock
!-----------------------------------------------------------------------------
!  create basic GSE datablock
!
    subroutine basicGSEDatablock(this,n,dt,y,net,station,comp)
    type (gse_datablock) :: this
    integer :: n,ierr
    real :: dt
    real, dimension(:), target :: y
    character (len=*) :: net,station,comp
!
    this%n = n
    this%dt = dt
    this%ry => y
    this%link = .true.
!
    call sff_PrepWid2(n, 1./dt, station, -1,-1,-1,-1,-1,comp,net,'NSP',-1.,-1.,-1.,-1.,-1.,this%wid,ierr)
    end subroutine basicGSEDatablock
!----------------------------------------------------------------
!  destructor
!
    subroutine deallocGSEDatablock(this)
    type (gse_datablock) :: this
    if(associated(this%y)) then
        if (this%link) then; nullify(this%y); else; deallocate(this%y); endif
    endif
    end subroutine deallocGSEDatablock
!------------------------------------------------------------------------------
!  return a pointer to the data
!
    function traceGSEDatablock(this)
    type (gse_datablock), intent(in) :: this
    integer, dimension(:), pointer :: traceGSEDatablock
!
    traceGSEDatablock => this%y
    end function traceGSEDatablock
!------------------------------------------------------------------------------
!  clear the data array
!
    subroutine clearDataGSEDatablock(this)
    type (gse_datablock) :: this
!
    if(associated(this%y)) deallocate(this%y)
    end subroutine clearDataGSEDatablock
!------------------------------------------------------------------------------
!  nullify the data array pointer
!
    subroutine nullifyDataGSEDatablock(this)
    type (gse_datablock) :: this
!
    if(associated(this%y)) nullify(this%y)
    end subroutine nullifyDataGSEDatablock
!------------------------------------------------------------------------------
!  get dt of datablock
!
    function getDtGSEDatablock(this)
    type (gse_datablock), intent(in) :: this
    real :: getDtGSEDatablock
    getDtGSEDataBlock = this%dt
    end function getDtGSEDatablock
!------------------------------------------------------------------------------
!  get station name of datablock
!
    character (len=5) function getStationGSEDatablock(this)
    type (gse_datablock), intent(in) :: this
    character (len=5) :: station 
    read(this%wid(30:34),'(a5)') station
    getStationGSEDatablock = station
    end function getStationGSEDatablock
!------------------------------------------------------------------------------
!  get channel of datablock
!
    character (len=3) function getChannelGSEDatablock(this)
    type (gse_datablock), intent(in) :: this
    read(this%wid(36:38),'(a3)') getChannelGSEDatablock
    end function getChannelGSEDatablock
!------------------------------------------------------------------------------
!  get start time of data block
!
    real function getTanfGSEDataBlock(this)
    type (gse_datablock), intent(in) :: this
    getTanfGSEDataBlock = this%tstart
    end function getTanfGSEDataBlock
!------------------------------------------------------------------------------
!  get year of data block
!
    integer function getYearGSEDataBlock(this)
    type (gse_datablock), intent(in) :: this
    read(this%wid(6:9),'(i4)') getYearGSEDataBlock
    end function getYearGSEDataBlock
!------------------------------------------------------------------------------
!  get day of year of data block
!
    integer function getDoyGSEDataBlock(this)
    type (gse_datablock), intent(in) :: this
    integer year,month,day,sff_TimeGetDOY
    read(this%wid(6:9),'(i4)') year
    read(this%wid(11:12),'(i2)') month
    read(this%wid(14:15),'(i2)') day
    getDoyGSEDataBlock = sff_TimeGetDOY(year,month,day)
    end function getDoyGSEDataBlock
!------------------------------------------------------------------------------
!  get nsamp of datablock 
!
    integer function getNGSEDatablock(this)
    type (gse_datablock), intent(in) :: this
    getNGSEDatablock = this%n
    end function getNGSEDatablock
!------------------------------------------------------------------------------
!  modify tstart of datablock, also change wid2line 
!
    subroutine modifyTanfGSEDatablock(this,tstart)
    type (gse_datablock) :: this
    real :: second
    real :: tstart
    integer :: day,hour,minute
!
    this%tstart = tstart
    call sff_TimeSplit(tstart, day, hour, minute, second)
    call sff_ModWid2time(this%wid, hour, minute, second)
    end subroutine modifyTanfGSEDatablock
!------------------------------------------------------------------------------
!  modify date (year month day) in wid2line 
!
    subroutine modifyDateGSEDatablock(this,year,month,day)
    type (gse_datablock) :: this
    integer :: year,month,day
!
    call sff_ModWid2date(this%wid, year, month, day)
    end subroutine modifyDateGSEDatablock
!----------------------------------------------------------------------------
!  write gse data block to file
!
    subroutine writeGSEDatablock(this,lu,last)
    type (gse_datablock) :: this
    integer :: lu
    logical :: last
    integer :: ierr
    integer, dimension(:), allocatable :: idata
!
    allocate(idata(this%n))
    call sff_WTraceNoScale(lu, this%wid, this%n, this%ry, idata, last, ierr)
    deallocate(idata)
    end subroutine writeGSEDatablock
!
 end module gseDatablock


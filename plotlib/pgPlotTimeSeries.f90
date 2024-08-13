!-------------------------------------------------------------
!> \brief Plot a time series object
!
 module pgPlotTimeSeries
	use timeSeries
 	implicit none
	interface new
		module procedure createPgPlotTimeSeries
	end interface
	interface display; module procedure displayPgPlotTimeSeries; end interface
	interface overlay; module procedure overlayPgPlotTimeSeries; end interface
	interface dealloc; module procedure deallocPgPlotTimeSeries; end interface
	interface operator (.xmin.); module procedure getXminPgPlotTimeSeries; end interface
	interface operator (.xmax.); module procedure getXmaxPgPlotTimeSeries; end interface
	interface operator (.ymin.); module procedure getYminPgPlotTimeSeries; end interface
	interface operator (.ymax.); module procedure getYmaxPgPlotTimeSeries; end interface
!
	type pgplot_time_series
		private
		type (time_series) :: ts
		real, dimension(:), pointer :: x
		real :: xmin,xmax,ymin,ymax,ch
		integer :: ls,lw,ci,symflag
		character (len=132) :: ylab,title,xopt,yopt
	end type pgplot_time_series
 contains
!--------------------------------------------------------------------
!  an object used to plot a time series using pgline or pgpnts
!
!  ts:	           time series
!  symflag:       symbol (1) or line (0) (def = 0)
!  lw:            line width (def = 2)
!  ci:            colour index (def = 1)
!  ls:            line style (def = 1)
!  ch:            character height (def = 1.5)
!  yopt:          string describing y-axis properties
!  ylab:	  label at y-axis (def = empty)
!  title: 	  title on top (def = empty)
!
	subroutine createPgPlotTimeSeries(this,ts,symflag,lw,ci,ls,ch,xopt,yopt,ylab,title)
	type (pgplot_time_series) :: this
	type (time_series) :: ts
	character (len=*), optional :: ylab,title,xopt,yopt
	integer, optional :: ls,lw,ci,symflag
	real, optional :: ch
	real :: tanf,dt
	integer :: i
!
	call shallowCopyTimeSeries(ts,this%ts)             ! just copy link to data
	allocate(this%x(.nsamp.ts))
	dt = .dt.ts; tanf = .tanf.ts
	do i=1,.nsamp.ts; this%x(i) = tanf + (i-1)*dt; enddo
	this%xmin = tanf; this%xmax = tanf + (.nsamp.ts-1)*dt
	this%ymin = minval(.trace.ts); this%ymax = maxval(.trace.ts)
	if (present(symflag)) then; this%symflag = symflag; else; this%symflag = 0; endif
	if (present(lw)) then; this%lw = lw; else; this%lw=2; endif
	if (present(ci)) then; this%ci = ci; else; this%ci=1; endif
	if (present(ls)) then; this%ls = ls; else; this%ls=1; endif
	if (present(ch)) then; this%ch = ch; else; this%ch=1.0; endif
	if (present(xopt)) then; this%xopt = xopt; else; this%xopt = 'ZXOHBCNST'; endif
	if (present(yopt)) then; this%yopt = yopt; else; this%yopt = 'BCNTS'; endif
	if (present(ylab)) then; this%ylab = ylab; else; this%ylab=''; endif
	if (present(title)) then; this%title = title; else; this%title=' '; endif
	end subroutine createPgPlotTimeSeries
!--------------------------------------------------------------------
!  free memory
!
	subroutine deallocPgPlotTimeSeries(this)
	type (pgplot_time_series) :: this
	call dealloc(this%ts)
	if (associated(this%x)) deallocate(this%x)
	end subroutine deallocPgPlotTimeSeries
!--------------------------------------------------------------------
!  plot graph
!
	subroutine displayPgPlotTimeSeries(this)
	type (pgplot_time_series) :: this
!
	call pgslw(1);call pgsci(1);call pgscf(2)
	call pgask(.false.)
	call pgeras
	call pgsvp(0.1,0.9,0.1,0.9)
	call pgswin(this%xmin,this%xmax,this%ymin,this%ymax)
	call pgslw(3); call pgsch(this%ch)
	call pgtbox(this%xopt,0.0,0,this%yopt,0.0,0)
	call pglab('Time',this%ylab,' ')
	call pgtext(0.3*this%xmax+0.7*this%xmin,1.05*this%ymax,trim(this%title))
	call pgslw(this%lw)
	call pgsci(this%ci)
	if(this%symflag == 0) then
		call pgline(.nsamp.(this%ts),this%x,.trace.(this%ts))
	else
		call pgpt(.nsamp.(this%ts),this%x,.trace.(this%ts),this%symflag)
	endif
	end subroutine displayPgPlotTimeSeries
!----------------------------------------------------------------------
!  overlay plot
!
	subroutine overlayPgPlotTimeSeries(this)
	type (pgplot_time_series) :: this
!
	call pgslw(this%lw)
	call pgsci(this%ci)
	if(this%symflag == 0) then
		call pgline(.nsamp.(this%ts),this%x,.trace.(this%ts))
	else
		call pgpt(.nsamp.(this%ts),this%x,.trace.(this%ts),this%symflag)
	endif
	end subroutine overlayPgPlotTimeSeries
!---------------------------------------------------------------------
!  set extrema
!
	subroutine setExtremaPgPlotTimeSeries(this,x1,x2,y1,y2)
	type (pgplot_time_series) :: this
	real :: x1,x2,y1,y2
	this%xmin=x1; this%xmax=x2; this%ymin=y1; this%ymax=y2
	end subroutine setExtremaPgPlotTimeSeries
!----------------------------------------------------------------------
!  get xmin
!
	real function getXminPgPlotTimeSeries(this)
	type (pgplot_time_series), intent(in) :: this
	getXminPgPlotTimeSeries = this%xmin
	end function getXminPgPlotTimeSeries
!----------------------------------------------------------------------
!  get xmax
!
	real function getXmaxPgPlotTimeSeries(this)
	type (pgplot_time_series), intent(in) :: this
	getXmaxPgPlotTimeSeries = this%xmax
	end function getXmaxPgPlotTimeSeries
!----------------------------------------------------------------------
!  get ymin
!
	real function getYminPgPlotTimeSeries(this)
	type (pgplot_time_series), intent(in) :: this
	getYminPgPlotTimeSeries = this%ymin
	end function getYminPgPlotTimeSeries
!----------------------------------------------------------------------
!  get ymax
!
	real function getYmaxPgPlotTimeSeries(this)
	type (pgplot_time_series), intent(in) :: this
	getYmaxPgPlotTimeSeries = this%ymax
	end function getYmaxPgPlotTimeSeries
!
 end module pgPlotTimeSeries

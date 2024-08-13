module pgPlotXY
    implicit none
    interface new
        module procedure createPgPlotXY
        module procedure createPgPlotY
    end interface
    interface display; module procedure displayPgPlotXY; end interface
    interface overlay; module procedure overlayPgPlotXY; end interface
    interface dealloc; module procedure deallocPgPlotXY; end interface
    interface operator (.xmin.); module procedure getXminPgPlotXY; end interface
    interface operator (.xmax.); module procedure getXmaxPgPlotXY; end interface
    interface operator (.ymin.); module procedure getYminPgPlotXY; end interface
    interface operator (.ymax.); module procedure getYmaxPgPlotXY; end interface
!
    type pgplot_xy
        private
        real, dimension(:), pointer :: x,y
        logical :: xlink                       ! x-values just a link
        real :: xmin,xmax,ymin,ymax,ch,symh,scale
        real :: xomin,xomax,yomin,yomax
        integer :: ls,lw,ci
        integer :: eqar,symflag
        character (len=132) :: xlab,ylab,title,xopt,yopt
    end type pgplot_xy
contains
!--------------------------------------------------------------------
!  an object used to plot an xy-graph using pgline or pgpnts
!
!  x:             x-values
!  y:             y-values
!  x1,xn,y1,yn:   extremal values (def = calculate from data)
!  eqar:          equal area or not (def = 0)
!  symflag:       boxcars (99), symbol (1) or line (0) (def = 0)
!  lw:            line width (def = 2)
!  ci:            colour index (def = 1)
!  ls:            line style (def = 1)
!  ch:            character height (def = 1.0)
!  symh:          symbol height (def = 1.0)
!  xopt:          string describing x-axis properties
!  yopt:          string describing y-axis properties
!  xlab:        label at x-axis (def = empty)
!  ylab:        label at y-axis (def = empty)
!  title:         title on top (def = empty)
!
    subroutine createPgPlotXY(this,x,y,x1,xn,y1,yn,eqar,symflag,lw,ci,ls,ch,symh,xopt,yopt,xlab,ylab,title)
    type (pgplot_xy) :: this
    real, dimension(:), target :: x,y
    real, optional :: x1,xn,y1,yn
    real :: xmin,xmax,ymin,ymax,xr,yr
    character (len=*), optional :: xlab,ylab,title,xopt,yopt
    integer, optional :: eqar,ls,lw,ci,symflag
    real, optional :: ch,symh
!
    this%x => x; this%y => y
    this%xlink = .true.
    xmin = minval(x); xmax = maxval(x); ymin = minval(y); ymax = maxval(y)
    xr = xmax-xmin; yr = ymax-ymin
    if (present(x1)) then; this%xmin = x1; else; this%xmin = xmin-0.03*xr; endif
    if (present(xn)) then; this%xmax = xn; else; this%xmax = xmax+0.03*xr; endif
    if (present(y1)) then; this%ymin = y1; else; this%ymin = ymin-0.03*yr; endif
    if (present(yn)) then; this%ymax = yn; else; this%ymax = ymax+0.03*yr; endif
    if (present(eqar)) then; this%eqar = eqar; else; this%eqar = 0; endif
    if (present(symflag)) then; this%symflag = symflag; else; this%symflag = 0; endif
    if (present(lw)) then; this%lw = lw; else; this%lw=2; endif
    if (present(ci)) then; this%ci = ci; else; this%ci=1; endif
    if (present(ls)) then; this%ls = ls; else; this%ls=1; endif
    if (present(ch)) then; this%ch = ch; else; this%ch=1.0; endif
    if (present(symh)) then; this%symh = symh; else; this%symh=1.0; endif
    if (present(xopt)) then; this%xopt = xopt; else; this%xopt = 'BCNTS'; endif
    if (present(yopt)) then; this%yopt = yopt; else; this%yopt = 'BCNTS'; endif
    if (present(xlab)) then; this%xlab = xlab; else; this%xlab=''; endif
    if (present(ylab)) then; this%ylab = ylab; else; this%ylab=''; endif
    if (present(title)) then; this%title = title; else; this%title=' '; endif
    this%xomin = this%xmin; this%xomax = this%xmax
    this%yomin = this%ymin; this%yomax = this%ymax
    this%scale = 1.0
!
    end subroutine createPgPlotXY
!-----------------------------------------------------------------------------------------------------
!  create a pgplotxy object from only y-values with given dx and x1
!
    subroutine createPgPlotY(this,dx,y,x1,y1,yn,eqar,symflag,lw,ci,ls,ch,symh,xopt,yopt,xlab,ylab,title)
    type (pgplot_xy) :: this
    real, dimension(:), target :: y
    real :: dx,x1
    real, optional :: y1,yn
    character (len=*), optional :: xlab,ylab,title,xopt,yopt
    integer, optional :: eqar,ls,lw,ci,symflag
    integer :: n,i
    real, optional :: ch,symh
!
    this%y => y
    n = size(y)
    this%xlink = .false.
    allocate(this%x(n))
    do i=1,n; this%x(i) = x1+(i-1)*dx; enddo
    this%xmin = x1
    this%xmax = x1+(n-1)*dx
    if (present(y1)) then; this%ymin = y1; else; this%ymin = minval(y); endif
    if (present(yn)) then; this%ymax = yn; else; this%ymax = maxval(y); endif
    if (present(eqar)) then; this%eqar = eqar; else; this%eqar = 0; endif
    if (present(symflag)) then; this%symflag = symflag; else; this%symflag = 0; endif
    if (present(lw)) then; this%lw = lw; else; this%lw=2; endif
    if (present(ci)) then; this%ci = ci; else; this%ci=1; endif
    if (present(ls)) then; this%ls = ls; else; this%ls=1; endif
    if (present(ch)) then; this%ch = ch; else; this%ch=1.0; endif
    if (present(symh)) then; this%symh = symh; else; this%symh=1.0; endif
    if (present(xopt)) then; this%xopt = xopt; else; this%xopt = 'BCNTS'; endif
    if (present(yopt)) then; this%yopt = yopt; else; this%yopt = 'BCNTS'; endif
    if (present(xlab)) then; this%xlab = xlab; else; this%xlab=''; endif
    if (present(ylab)) then; this%ylab = ylab; else; this%ylab=''; endif
    if (present(title)) then; this%title = title; else; this%title=' '; endif
    this%xomin = this%xmin; this%xomax = this%xmax
    this%yomin = this%ymin; this%yomax = this%ymax
    this%scale = 1.0
!
    end subroutine createPgPlotY
!--------------------------------------------------------------------
!  free memory
!
    subroutine deallocPgPlotXY(this)
    type (pgplot_xy) :: this
    nullify(this%y)
    if (this%xlink) then
        nullify(this%x)
    else
        deallocate(this%x)
    endif
    end subroutine deallocPgPlotXY
!--------------------------------------------------------------------
!  plot graph
!
    subroutine displayPgPlotXY(this)
    type (pgplot_xy) :: this
    integer ::j
    real :: dxl,dxr
!
    call pgslw(1);call pgsci(1);call pgscf(2)
    call pgask(.false.)
    call pgeras
    call pgsvp(0.1,0.9,0.12,0.9)
    if (this%eqar == 0) then
        call pgswin(this%xmin,this%xmax,this%ymin,this%ymax)
    else
        call pgwnad(this%xmin,this%xmax,this%ymin,this%ymax)
    endif
    call pgslw(3); call pgsch(this%ch)
    call pgbox(this%xopt,0.0,0,this%yopt,0.0,0)
    call pglab(this%xlab,this%ylab,' ')
    call pgtext(0.3*this%xmax+0.7*this%xmin,-0.05*this%ymin+1.05*this%ymax,trim(this%title))
    call pgslw(this%lw)
    call pgsci(this%ci)
    if (this%symflag == 0) then
        call pgline(size(this%x),this%x,this%scale*this%y)
    else if (this%symflag == 99) then
        do j = 1,size(this%x)
            if (j == 1) then; dxl = 0.0; else; dxl = 0.5*(this%x(j)-this%x(j-1)); endif
            if (j == size(this%x)) then; dxr = 0.0; else; dxr = 0.5*(this%x(j+1)-this%x(j)); endif
            call pgsfs(2); call pgrect(this%x(j)-dxl,this%x(j)+dxr,0.,this%scale*this%y(j))
        enddo
    else 
        call pgsch(this%symh)
        call pgpt(size(this%x),this%x,this%scale*this%y,this%symflag)
    endif
    end subroutine displayPgPlotXY
!----------------------------------------------------------------------
!  overlay plot
!
    subroutine overlayPgPlotXY(this)
    type (pgplot_xy) :: this
    integer ::j
    real :: dxl,dxr
!
    call pgslw(this%lw)
    call pgsci(this%ci)
    if(this%symflag == 0) then
        call pgline(size(this%x),this%x,this%scale*this%y)
    else if (this%symflag == 99) then
        do j = 1,size(this%x)
            if (j == 1) then; dxl = 0.0; else; dxl = 0.5*(this%x(j)-this%x(j-1)); endif
            if (j == size(this%x)) then; dxr = 0.0; else; dxr = 0.5*(this%x(j+1)-this%x(j)); endif
            call pgsfs(2); call pgrect(this%x(j)-dxl,this%x(j)+dxr,0.,this%scale*this%y(j))
        enddo
    else
        call pgsch(this%symh)
        call pgpt(size(this%x),this%x,this%scale*this%y,this%symflag)
    endif
    end subroutine overlayPgPlotXY
!---------------------------------------------------------------------
!  set extrema
!
    subroutine setExtremaPgPlotXY(this,x1,x2,y1,y2)
    type (pgplot_xy) :: this
    real :: x1,x2,y1,y2
    this%xmin=x1; this%xmax=x2; this%ymin=y1; this%ymax=y2
    end subroutine setExtremaPgPlotXY
!----------------------------------------------------------------------
!  set original extrema
!
    subroutine setOriginalExtremaPgPlotXY(this)
    type (pgplot_xy) :: this
    this%xmin = this%xomin; this%xmax = this%xomax; this%ymin = this%yomin; this%ymax = this%yomax
    end subroutine setOriginalExtremaPgPlotXY
!---------------------------------------------------------------------
!  set scale
!
    subroutine setScalePgPlotXY(this,s)
    type (pgplot_xy) :: this
    real :: s
    this%scale = s
    end subroutine setScalePgPlotXY
!----------------------------------------------------------------------
!  get xmin
!
    real function getXminPgPlotXY(this)
    type (pgplot_xy), intent(in) :: this
    getXminPgPlotXY = this%xmin
    end function getXminPgPlotXY
!----------------------------------------------------------------------
!  get xmax
!
    real function getXmaxPgPlotXY(this)
    type (pgplot_xy), intent(in) :: this
    getXmaxPgPlotXY = this%xmax
    end function getXmaxPgPlotXY
!----------------------------------------------------------------------
!  get ymin
!
    real function getYminPgPlotXY(this)
    type (pgplot_xy), intent(in) :: this
    getYminPgPlotXY = this%ymin
    end function getYminPgPlotXY
!----------------------------------------------------------------------
!  get ymax
!
    real function getYmaxPgPlotXY(this)
    type (pgplot_xy), intent(in) :: this
    getYmaxPgPlotXY = this%ymax
    end function getYmaxPgPlotXY
!
 end module pgPlotXY

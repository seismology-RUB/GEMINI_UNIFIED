! ===============================================================================
!  Plot a node earth model
! ===============================================================================
!----------------------------------------------------------------------------
!   Copyright 2017 Wolfgang Friederich (Ruhr-Universitaet Bochum, Germany)
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
!  Plot a node earth model
!-----------------------------------------------------------------
program plotNodeEarthmodel
    use nodeEarthmodel
    use splineEarthmodelCoefficients
    use argumentParser
    use pgPlotWindow
    use pgPlotXY
    use changeAxesLimitsPgplot
    use string
    use errorMessage
    use realloc
    implicit none
    type (node_earthmodel) :: em
    type (error_message) :: errmsg
    type (argument_parser) :: ap
    type (pgplot_window) :: pgwin
    type (pgplot_xy), target :: elcon_graph,a_graph,f_graph,l_graph,n_graph,vp_graph,vs_graph,ro_graph
    type (pgplot_xy), pointer :: pgraph
    integer :: n,j,nl,nlay,nd,nlb
    real, dimension(:), pointer :: r,rho,ar,cr,fr,lr,nr,kapr,mur,vp,vs
    double precision, dimension(:), pointer :: rb
    double precision :: ro
    real :: w,asprat,rcur,dr,re,elconmax,xmin,xmax,ymin,ymax,x1,x2,y1,y2,xa,ya,rbot,vmax
    double complex, dimension(7) :: zelcon
    logical :: leave,success,stay_layer,plotvel
    character (len=1) :: ch
    character (len=max_length_string) :: emfile
    character (len=18) :: myname = 'plotNodeEarthmodel'
    data n/500/
    !-------------------------------------------------------------------------------
    call init(ap,myname,'Plot node earth model')
    call addPosarg(ap,'modelfile','sval',' Earth model file')
    call addOption(ap,'-n',.true.,'number of points to be evaluated','ival','500')
    call addOption(ap,'-w',.true.,'plot window width','rval','14.')
    call addOption(ap,'-a',.true.,'plot window aspect ratio','rval','0.65')
    call addOption(ap,'-nlb',.true.,'deepest layer considered for plot','ival','1')
    call addOption(ap,'-vel',.false.,'plot velocities')
    call parse(ap)
    emfile = ap.sval.'modelfile'
    w = ap.rval.'-w'
    asprat = ap.rval.'-a'
    n = ap.ival.'-n'
    nlb = ap.ival.'-nlb'
    plotvel = ap.optset.'-vel'
    if (.level.(.errmsg.ap) == 2) then; call usage(ap); call print(.errmsg.ap); stop; endif
    call dealloc(ap)
    !-----------------------------------------------------------------------------
    call createNodeEarthmodel(em,1,emfile,errmsg)
    if (.level.errmsg == 2) then
       call print(errmsg); stop
    endif
    call computeSplineEarthmodelCoefficients(em,.fref.em,'ELASTIC',errmsg)
    if (.level.errmsg == 2) then
       call print(errmsg); stop
    endif
    rb => getRbArrayNodeEarthmodel(em)
    nlay = .nlay.em
    nd = n+2*nlay
    if (nlb .gt. nlay) then
       call add(errmsg,2,'plot bottom layer is above top layer',myname)
       call print(errmsg); stop
    endif
    !
    !  set up plot data
    !
    allocate(r(nd),rho(nd),ar(nd),cr(nd),fr(nd),lr(nd),nr(nd),kapr(nd),mur(nd),vp(nd),vs(nd))
    re = .rearth.em
    dr = re/n
    rbot = 0.0
    if (nlb > 1) then
        dr = (re-rb(nlb-1))/n
        rbot = rb(nlb-1)
    endif
    elconmax = 0.d0
    rcur = rbot
    stay_layer = .true.
    j = 1                     ! node count
    do nl = nlb,nlay
       do while (rcur .le. rb(nl))
          call evalComplexElasticConstantsSplineEarthmodel(dble(rcur),nl,ro,zelcon)
          elconmax = max(elconmax,maxval(real(zelcon)))
          r(j) = rcur
          rho(j) = ro
          ar(j) = real(zelcon(1))
          cr(j) = real(zelcon(2))
          fr(j) = real(zelcon(3))
          lr(j) = real(zelcon(4))
          nr(j) = real(zelcon(5))
          kapr(j) = real(zelcon(6))
          mur(j) = real(zelcon(7))
          vp(j) = sqrt((kapr(j)+4d0*mur(j)/3.d0)/ro)
          vs(j) = sqrt(mur(j)/ro)
          print *,j,nl,re-r(j),r(j),ro,kapr(j),mur(j),vp(j),vs(j),vp(j)/vs(j)
          rcur = rcur+dr
          j = j+1
          if (rcur > rb(nl)) then
             rcur = rb(nl)
             if (stay_layer) then
                stay_layer = .false.
             else
                stay_layer = .true.
                exit
             endif
          endif
       enddo
    enddo
    nd = j-1
    !
    !  reallocate arrays
    !
    r => reallocate(r,nd)
    rho => reallocate(rho,nd)
    ar => reallocate(ar,nd)
    cr => reallocate(cr,nd)
    fr => reallocate(fr,nd)
    lr => reallocate(lr,nd)
    nr => reallocate(nr,nd)
    kapr => reallocate(kapr,nd)
    mur => reallocate(mur,nd)
    vp => reallocate(vp,nd)
    vs => reallocate(vs,nd)
    vmax = maxval(vp)
    !
    !  Initialize pgplotxy objects
    !
    call new(elcon_graph,cr,r,x1 = 0.0,xn = elconmax,y1 = rbot,yn = re,lw = 4,ci = 1,xlab = 'Elastic constants (GPa)',&
         & ylab = 'Radius',title = 'Earth model:'+getBehindLastSeparatorString(emfile,'/',errmsg))
    call new(a_graph,ar,r,x1 = 0.0,xn = elconmax,lw = 4,ci = 2)
    call new(f_graph,fr,r,x1 = 0.0,xn = elconmax,lw = 4,ci = 3)
    call new(l_graph,lr,r,x1 = 0.0,xn = elconmax,lw = 4,ci = 4)
    call new(n_graph,nr,r,x1 = 0.0,xn = elconmax,lw = 4,ci = 5)
    call new(vp_graph,vp,r,x1 = 0.0,xn = vmax,y1 = rbot,yn = re,lw = 4,ci = 1,xlab = 'Velocities (m/s)',&
         & ylab = 'Radius',title = 'Earth model:'+getBehindLastSeparatorString(emfile,'/',errmsg))
    call new(vs_graph,vs,r,x1 = 0.0,xn = vmax,lw = 4,ci = 2)
    call new(ro_graph,rho,r,x1 = 0.0,xn = vmax,lw = 4,ci = 3)
    !
    !  open plot window and draw figure
    !
    call new(pgwin,width = w,aspect = asprat)
    !
    !  user interaction
    !
    leave = .false.
    do while (.not. leave)
       if (.not.plotvel) then
          pgraph => elcon_graph
          call display(elcon_graph)
          call overlay(a_graph); call overlay(f_graph)
          call overlay(l_graph); call overlay(n_graph)
          xmin = .xmin.elcon_graph; xmax = .xmax.elcon_graph
          ymin = .ymin.elcon_graph; ymax = .ymax.elcon_graph
       else
          pgraph => vp_graph
          call display(vp_graph)
          call overlay(vs_graph)
          call overlay(ro_graph)
          xmin = .xmin.vp_graph; xmax = .xmax.vp_graph
          ymin = .ymin.vp_graph; ymax = .ymax.vp_graph
       endif
       call pgcurs(xa,ya,ch)
       select case (ch)
       case ('A','D','X')
          call pickAxesLimitsPgplot(ch,xa,ya,x1,x2,y1,y2,success)
          if (success) then
             if (ch == 'A') then
                call setExtremaPgPlotXY(pgraph,x1,x2,y1,y2)
             else if (ch == 'D') then
                call setExtremaPgPlotXY(pgraph,x1,x2,ymin,ymax)
             else if (ch == 'X') then
                call setExtremaPgPlotXY(pgraph,xmin,xmax,y1,y2)
             endif
          endif
       case ('z','Z','l','r','u','d')
          call shiftAxesLimitsPgplot(ch,xmin,xmax,ymin,ymax)
          call setExtremaPgPlotXY(pgraph,xmin,xmax,ymin,ymax)
       case ('o')
          call setOriginalExtremaPgPlotXY(pgraph)
       case ('x'); leave = .true.
       end select
    enddo

    !
    !  clean up
    !
    deallocate(r,rho,ar,cr,fr,lr,nr,kapr,mur,vp,vs)
    call dealloc(em)
    call dealloc(elcon_graph)
    call dealloc(vp_graph)
end program

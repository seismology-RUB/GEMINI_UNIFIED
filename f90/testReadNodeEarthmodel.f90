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
!---------------------------------------------------------------------
!  Test reading node earth model stored in a JSON file
!-----------------------------------------------------------------
program testReadNodeEarthmodel
    use argumentParser
    use nodeEarthmodel
    use errorMessage
    use checkAssertion
    implicit none
    type (node_earthmodel) :: em
    type (error_message) :: errmsg
    type (argument_parser) :: ap
    real :: p,rk32,vp32,rk18,vs18,rho43
    integer :: cnt,nlay,nk,iktop4,ikbot3,nl,top,nl5368
    logical :: res
    character (len=max_length_string) :: emfile,testmodel
    character (len=22) :: myname = 'testReadNodeEarthmodel'
    data rk32/6301.0/vp32/8.04415321/rk18/5961.0/vs18/5.09043503/rho43/2.72000003/
    data nlay/6/nk/44/iktop4/40/ikbot3/10/nl5368/1/
!-------------------------------------------------------------------------------
    call init(ap,myname,'Test reading of JSON earth model file')
    call addPosarg(ap,'modelfile','sval',' Earth model file')
    call parse(ap)
    emfile = ap.sval.'modelfile'
    if (.level.(.errmsg.ap) == 2) then; call usage(ap); call print(.errmsg.ap); stop; endif
    call dealloc(ap)
!-----------------------------------------------------------------------------
    call createNodeEarthmodel(em,1,emfile,errmsg)
    !
    !  do tests
    !
    testmodel = getBehindLastSeparatorString(emfile,'/',errmsg)
    if (.level.errmsg == 2) then
       call print(errmsg)
       stop
    endif
    cnt = 0
    if (testmodel.equal.'ak135q.json') then
       p = getPropertyNodeEarthmodel(em,'vp',32)
       call genericCheckAssertion('vp at node 32',p,vp32,1.e-7,res)
       if (res) cnt = cnt+1
       p = getPropertyNodeEarthmodel(em,'node_radius',32)
       call genericCheckAssertion('radius of node 32',p,rk32,1.e-7,res)
       if (res) cnt = cnt+1
       p = getPropertyNodeEarthmodel(em,'vs',18)
       call genericCheckAssertion('vs at node 18',p,vs18,1.e-7,res)
       if (res) cnt = cnt+1
       p = getPropertyNodeEarthmodel(em,'node_radius',18)
       call genericCheckAssertion('radius of node 18',p,rk18,1.e-7,res)
       if (res) cnt = cnt+1
       p = getPropertyNodeEarthmodel(em,'rho',43)
       call genericCheckAssertion('density at node 43',p,rho43,1.e-7,res)
       if (res) cnt = cnt+1
       call genericCheckAssertion('number of nodes',.nk.em,nk,res)
       if (res) cnt = cnt+1
       call genericCheckAssertion('number of layers',.nlay.em,nlay,res)
       if (res) cnt = cnt+1
       call genericCheckAssertion('iktop of layer 4',em.iktop.4,iktop4,res)
       if (res) cnt = cnt+1
       call genericCheckAssertion('ikbot of layer 3',em.ikbot.3,ikbot3,res)
       if (res) cnt = cnt+1
       call getLayerIndexNodeEarthmodel(em,5368.d0,nl,top)
       call genericCheckAssertion('layer index of r=5368',nl,nl5368,res)
       if (res) cnt = cnt+1
       print *,cnt,' of 10 assertions satisfied'
    endif
    call dealloc(em)
end program

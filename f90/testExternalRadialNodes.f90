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
!  Test correct setup of external radial nodes and properties at nodes
!-----------------------------------------------------------------
program testExternalRadialNodes
    use argumentParser
    use nodeEarthmodel
    use externalRadialNodes
    use errorMessage
    use checkAssertion
    implicit none
    type (node_earthmodel) :: em
    type (external_radial_nodes) :: exnod
    type (error_message) :: errmsg
    type (argument_parser) :: ap
    double precision, dimension(:), pointer :: rnod
    integer, dimension(:), pointer :: state, layer
    integer :: j
    character (len=max_length_string) :: emfile,parfile
    character (len=23) :: myname = 'testExternalRadialNodes'
!-------------------------------------------------------------------------------
    call init(ap,myname,'Test setup of external radial nodes')
    call addPosarg(ap,'modelfile','sval',' Earth model file')
    call addPosarg(ap,'parfile','sval',' Parameter file')
    call parse(ap)
    emfile = ap.sval.'modelfile'
    parfile = ap.sval.'parfile'
    if (.level.(.errmsg.ap) == 2) then; call usage(ap); call print(.errmsg.ap); stop; endif
    call dealloc(ap)
!-----------------------------------------------------------------------------
    call createNodeEarthmodel(em,1,emfile,errmsg)
    if (.level.errmsg == 2) then; call print(errmsg); stop; endif
    call createExternalRadialNodes(exnod,1,parfile,em,errmsg)
    if (.level.errmsg == 2) then; call print(errmsg); endif
!
    print *,.nnod.exnod   
    rnod => getPointerDoubleRadiiExternalRadialNodes(exnod)
    state => getPointerStateExternalRadialNodes(exnod)
    layer => getPointerLayerExternalRadialNodes(exnod)
    do j = 1,.nnod.exnod
       print *,j,rnod(j),state(j),layer(j)
    enddo
!
    call dealloc(exnod)
    call dealloc(em)
end program

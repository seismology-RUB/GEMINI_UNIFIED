! ======================================================================
!  External radial nodes
! ======================================================================
!----------------------------------------------------------------------------
!   Copyright 2016 Wolfgang Friederich (Ruhr-Universitaet Bochum, Germany)
!
!   This file is part of GEMINI_UNIFIED version 1.0.
!   This file is part of ASKI version 1.2.
!
!   GEMINI_UNIFIED version 1.0 and ASKI version 1.2 are free software:
!   you can redistribute it and/or modify
!   it under the terms of the GNU General Public License as published by
!   the Free Software Foundation, either version 2 of the License, or
!   (at your option) any later version.
!
!   GEMINI_UNIFIED version 1.0 and ASKI 1.2 are is distributed
!   in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU General Public License for more details.
!
!   You should have received a copy of the GNU General Public License
!   along with GEMINI_UNIFIED version 1.0 and ASKI version 1.2.
!   If not, see <http://www.gnu.org/licenses/>.
!-----------------------------------------------------------------
!> \brief Module dealing with external radial nodes
!
!  Definition of external radial nodes (in terms of depth):
!  For an arbitrary number of EXTERNAL_NODES_NBLOCKS blocks of nodes,
!  the vectors EXTERNAL_NODES_NNOD and EXTERNAL_NODES_DR (both of length EXTERNAL_NODES_NBLOCKS)
!  define the subdivision of each block by its number of nodes EXTERNAL_NODES_NNOD(i) contained 
!  in such a block, having constant spacing EXTERNAL_NODES_DR(i). The uppermost (last) node is always
!  slightly below the surface and is not counted in EXTERNAL_NODES_NNOD.
!  It can be shifted to greater depths using EXTERNAL_NODES_SHIFT.
!  Units of radii should be given in SI units (i.e. meter)
!
!  Example: 
!      EXTERNAL_NODES_NBLOCKS =  3
!      EXTERNAL_NODES_NNOD =  8   5   2
!      EXTERNAL_NODES_DR =  5. 10. 20.
!  this means:
!  3 blocks of NODES (one below the other starting with the first BELOW the surface) with different spacing:
!     8 nodes with spacing of 5 km between surface and 40 km depth (5,10,15,20,25,30,35,40) plus one at surface
!     5 nodes with spacing of of 10 km between 40 km and 90 km depth (50,60,70,80,90)
!     2 layers with spacing of 20 km between 90 km and 130 km depth (110,130)
!  Note: Numbering of nodes is from bottom up !!
!
!  Alternatively, external radial nodes can also be read from an HDF file.
!
!  Beside the so defined nodes, further nodes are added just below and above all interfaces of the model.
!  If a node falls on an interface, it is split into two, one below and one above.
!----------------------------------------------------------------
module externalRadialNodes
       use hdf5
    use errorMessage
    use inputParameter
    use nodeEarthmodel
    use splineEarthmodelCoefficients
    use complexElasticConstants
    use string
    use hdfWrapper
    use anyRankRealArray
    implicit none
    interface dealloc; module procedure deallocExternalRadialNodes; end interface
    interface operator (.nnod.); module procedure getNnodExternalRadialNodes; end interface
    type external_radial_nodes
       private
       double precision :: rearth                                      ! earth radius taken as surface
       double precision, dimension(:), pointer :: rnod => null()       ! nodes radii
       integer, dimension(:), pointer :: state => null()               ! fluid (1) or solid (0) material at node
       integer, dimension(:), pointer :: layer => null()               ! layer the node is in
       integer :: nnod                                                 ! number of radii
       double precision, dimension(:), allocatable :: ro               ! density
       double complex, dimension(:,:), allocatable :: elcon            ! elastic constants (A,C,F,L,N,Kap,Mue) 
       double precision, dimension(:,:), allocatable :: qinv           ! inverse quality factors (qkinv,qminv) 
     end type external_radial_nodes
!
contains
!-------------------------------------------------------------------------------
!  \brief Read out parameters from file, nodes are created on the fly upon request
!
    subroutine createExternalRadialNodes(this,lu,parfile,nem,errmsg)
    type (external_radial_nodes) :: this
    integer :: lu
    character (len=*) :: parfile
    type (node_earthmodel) :: nem
    type (error_message) :: errmsg
    integer :: nblocks                                              ! total number of nodes
    integer, dimension(:), allocatable :: nnodblocks                ! array with nodes per block
    double precision, dimension(:), allocatable :: drblocks         ! array with node distance per block
    double precision, dimension(:), allocatable :: prn              ! array with preliminary radial nodes
    double precision :: shift                                       ! shift of "surface radius"
    type (input_parameter) :: inpar    
    integer :: ios,n,ib,i,nnod
    character (len=25) :: myname = 'createExternalRadialNodes'
    character (len=80), dimension(5) :: par_keys
    data par_keys/'EXTERNAL_NODES_NBLOCKS','EXTERNAL_NODES_NNOD',&
         'EXTERNAL_NODES_DR','EXTERNAL_NODES_SHIFT','EXTERNAL_NODES_REARTH'/
!
!  read from  parameter file
!
    call addTrace(errmsg,myname)
    call createKeywordsInputParameter(inpar,par_keys)
    call readSubroutineInputParameter(inpar,lu,parfile,errmsg)
    if (.level.errmsg == 2) then; call dealloc(inpar); return; endif
    nblocks = ival(inpar,'EXTERNAL_NODES_NBLOCKS')
    allocate(nnodblocks(nblocks))
    allocate(drblocks(nblocks))
    nnodblocks = ivec(inpar,'EXTERNAL_NODES_NNOD',nblocks,ios)
    if (ios /= 0) then
       call add(errmsg,2,'problems reading out EXTERNAL_NODES_NNOD',myname)
       goto 1
    endif
    drblocks = dvec(inpar,'EXTERNAL_NODES_DR',nblocks,ios)
    if (ios /= 0) then
       call add(errmsg,2,'problems reading out EXTERNAL_NODES_DR',myname)
       goto 1
    endif
    shift = dval(inpar,'EXTERNAL_NODES_SHIFT')
    this%rearth = dval(inpar,'EXTERNAL_NODES_REARTH')
    call dealloc(inpar)
!
!  calculate node radii
!
    nnod = sum(nnodblocks)+1                            ! +1 is surface node
    allocate(prn(nnod))
    prn(nnod) = this%rearth-shift
    n = nnod
    do ib = 1,nblocks
       do i = 1,nnodblocks(ib)
          prn(n-i) = prn(n)-i*drblocks(ib)
       enddo
       n = n-nnodblocks(ib)
    enddo
!
    if (any(prn < 0.d0)) then
       call add(errmsg,2,'Negative radii, check your input',myname)
       goto 1
    endif
!
!  check for radial nodes very close to discontinuities of earth model
!  if so, split it into two, one just below and one just above
!  also generally add additional nodes just below and just above any interface
!
    call splitNodesAtDiscosExternalRadialNodes(this,nnod,prn,nem,errmsg)
!
!  compute material parameters at reference frequency at nodes
!
    call elasticConstantsExternalRadialNodes(this,nem,errmsg)
    if (.level.errmsg == 2) return
!
1   if (allocated(drblocks)) deallocate(drblocks)
    if (allocated(nnodblocks)) deallocate(nnodblocks)
    if (allocated(prn)) deallocate(prn)
    call dealloc(inpar)
    end subroutine createExternalRadialNodes
!-------------------------------------------------------------------------------
!  Create external radial nodes from radii listed in a HDF file
!  HDF environment is not opened here
!
    subroutine createFromRadiiExternalRadialNodes(this,lu,parfile,nem,errmsg)
    type (external_radial_nodes) :: this
    integer :: lu
    character (len=*) :: parfile
    type (node_earthmodel) :: nem
    type (error_message) :: errmsg
    type (input_parameter) :: inpar
    type (any_rank_real_array) :: arra
    real, dimension(:), pointer :: prn
    integer(hid_t) :: fid
    integer :: nnod,ierr
    character(len=max_length_string) :: hdffile,hdfpath
    character (len=23) :: myname = 'readExternalRadialNodes'
    character (len=80), dimension(3) :: par_keys
    data par_keys/'EXTERNAL_NODES_HDF_FILE','EXTERNAL_NODES_HDF_PATH','EXTERNAL_NODES_REARTH'/
!
!  read from  parameter file
!
    call addTrace(errmsg,myname)
    call createKeywordsInputParameter(inpar,par_keys)
    call readSubroutineInputParameter(inpar,lu,parfile,errmsg)
    if (.level.errmsg == 2) then; call dealloc(inpar); return; endif
    hdffile = inpar.sval.'EXTERNAL_NODES_HDF_FILE'
    hdfpath = inpar.sval.'EXTERNAL_NODES_HDF_PATH'
    this%rearth = dval(inpar,'EXTERNAL_NODES_REARTH')
!
!  read from HDF file
!
    call openFileRoHDFWrapper(trim(hdffile),fid,errmsg)
    if (.level.errmsg == 2) return
    call readArrayHDFWrapper(fid,trim(hdfpath),arra,errmsg)
    if (.level.errmsg == 2) return
    prn => arra%get1d()
    nnod = size(prn)
    call splitNodesAtDiscosExternalRadialNodes(this,nnod,dble(prn),nem,errmsg)
    call dealloc(inpar)
    deallocate(prn)
    call h5fclose_f(fid,ierr)
!
!  compute material parameters at reference frequency at nodes
!
    call elasticConstantsExternalRadialNodes(this,nem,errmsg)
    if (.level.errmsg == 2) return
    end subroutine createFromRadiiExternalRadialNodes
!-------------------------------------------------------------------------------
!  Write full external radial nodes instance into an existing HDF file
!  at locid
!
    subroutine writeHDFExternalRadialNodes(this,locid,errmsg)
    type (external_radial_nodes) :: this
    type (error_message) :: errmsg
    type (any_rank_real_array) :: arra
    type (any_rank_integer_array) :: aria
    integer(hid_t) :: locid
    real, dimension(:), allocatable, target :: p
    real, dimension(:,:), allocatable, target :: p2
    character (len=27) :: myname = 'writeHDFExternalRadialNodes'
!
    call addTrace(errmsg,myname)
!
    p = (/real(this%rearth)/)
    call arra%assoc1d(p)
    call writeArrayHDFWrapper(locid,'rearth',arra,errmsg)
    if (.level.errmsg == 2) return
    call arra%deassoc(); deallocate(p)
!
    allocate(p(this%nnod)); p = real(this%rnod)
    call arra%assoc1d(p)
    call writeArrayHDFWrapper(locid,'radii',arra,errmsg)
    if (.level.errmsg == 2) return
    call arra%deassoc(); deallocate(p)
!
    call aria%assoc1d(this%state)
    call writeArrayHDFWrapper(locid,'state',aria,errmsg)
    if (.level.errmsg == 2) return
!
    call aria%assoc1d(this%layer)
    call writeArrayHDFWrapper(locid,'layer',aria,errmsg)
    if (.level.errmsg == 2) return
    call aria%deassoc()
!
    allocate(p(this%nnod)); p = real(this%ro)
    call arra%assoc1d(p)
    call writeArrayHDFWrapper(locid,'density',arra,errmsg)
    if (.level.errmsg == 2) return
    call arra%deassoc(); deallocate(p)
!
    allocate(p2(2,this%nnod)); p2 = real(this%qinv)
    call arra%assoc2d(p2)
    call writeArrayHDFWrapper(locid,'inverseQ',arra,errmsg)
    if (.level.errmsg == 2) return
    call arra%deassoc(); deallocate(p2)
!
    allocate(p2(7,this%nnod)); p2 = real(this%elcon)
    call arra%assoc2d(p2)
    call writeArrayHDFWrapper(locid,'realElcon',arra,errmsg)
    if (.level.errmsg == 2) return
    call arra%deassoc()
!
    p2 = aimag(this%elcon)
    call arra%assoc2d(p2)
    call writeArrayHDFWrapper(locid,'imagElcon',arra,errmsg)
    if (.level.errmsg == 2) return
    call arra%deassoc(); deallocate(p2)
!
    end subroutine writeHDFExternalRadialNodes
!-------------------------------------------------------------------------------
!  Read a full external radial nodes instance from a HDF file
!  HDF environment is not opened here
!
    subroutine readExternalRadialNodes(this,locid,errmsg)
    type (external_radial_nodes) :: this
    integer(hid_t) :: locid
    type (error_message) :: errmsg
    type (any_rank_real_array) :: arra
    type (any_rank_integer_array) :: aria
    real, dimension(:), pointer :: p
    real, dimension(:,:), pointer :: p2
    real, dimension(:,:), pointer :: q2
!
    call readArrayHDFWrapper(locid,'rearth',arra,errmsg)
    if (.level.errmsg == 2) return
    p => arra%get1d(); this%rearth = p(1); call arra%deassoc()
!
    call readArrayHDFWrapper(locid,'radii',arra,errmsg)
    if (.level.errmsg == 2) return
    p => arra%get1d()
    this%nnod = size(p)
    allocate(this%rnod(this%nnod))
    this%rnod = p
    deallocate(p)
!
    call readArrayHDFWrapper(locid,'state',aria,errmsg)
    if (.level.errmsg == 2) return
    this%state => aria%get1d(); call aria%deassoc()
!
    call readArrayHDFWrapper(locid,'layer',aria,errmsg)
    if (.level.errmsg == 2) return
    this%layer => aria%get1d(); call aria%deassoc()
!
    call readArrayHDFWrapper(locid,'density',arra,errmsg)
    if (.level.errmsg == 2) return
    p => arra%get1d()
    allocate(this%ro(size(p)))
    this%ro = p
    deallocate(p)
!
    call readArrayHDFWrapper(locid,'inverseQ',arra,errmsg)
    if (.level.errmsg == 2) return
    p2 => arra%get2d()
    allocate(this%qinv(size(p2,1),size(p2,2)))
    this%qinv = p2
    deallocate(p2)
!
    call readArrayHDFWrapper(locid,'realElcon',arra,errmsg)
    if (.level.errmsg == 2) return
    p2 => arra%get2d()
    call readArrayHDFWrapper(locid,'imagElcon',arra,errmsg)
    if (.level.errmsg == 2) return
    q2 => arra%get2d()
    allocate(this%elcon(size(p2,1),size(p2,2)))
    this%elcon = dcmplx(p2,q2)
    deallocate(p2,q2)
    end subroutine readExternalRadialNodes
!------------------------------------------------------------------------
!  Deallocate object
!
    subroutine deallocExternalRadialNodes(this)
    type (external_radial_nodes) :: this
    if (associated(this%rnod)) deallocate(this%rnod)
    if (associated(this%state)) deallocate(this%state)
    if (associated(this%layer)) deallocate(this%layer)
    end subroutine deallocExternalRadialNodes
!------------------------------------------------------------------------
!  Split preliminary nodes at internal dicontinuities
!
    subroutine splitNodesAtDiscosExternalRadialNodes(this,nnod,prn,nem,errmsg)
    type (external_radial_nodes) :: this
    integer :: nnod
    double precision, dimension(:) :: prn
    type (node_earthmodel) :: nem
    type (error_message) :: errmsg
    double precision :: drmin
    integer :: n,ib,je,ne,iflso
    double precision, dimension(:), pointer :: rb                   ! radii of upper layer boundaries
!    
    rb => getRbArrayNodeEarthmodel(nem)
    drmin = rb(size(rb))*1.d-4
    allocate(this%rnod(nnod+2*size(rb)))
    allocate(this%state(nnod+2*size(rb)))
    allocate(this%layer(nnod+2*size(rb)))
    je = 0                                ! index of last this%rnod assigned 
    ne = 0                                ! index of last used prn
    do ib = 1,size(rb)
       n = locate(rb(ib),nnod,prn)
       if (n == 0) cycle
       iflso = nem.iflso.ib
       this%rnod(je+1:je+n-ne) = prn(ne+1:n)
       this%layer(je+1:je+n-ne) = ib 
       this%state(je+1:je+n-ne) = iflso 
       je = je+n-ne
       if (n .ge. nnod) exit
    !
       if (dabs(prn(n+1)-rb(ib)) < drmin .and. dabs(prn(n)-rb(ib)) < drmin) then
          call add(errmsg,2,'External nodes too close to each other, use bigger spacing','splitNodesAtDiscosExternalRadialnodes')
          return
       endif
    !
    !  upper node very close to disco
    !
       if (dabs(prn(n+1)-rb(ib)) < drmin) then
    !      print *,'ERN: prn(n+1) on disco (closer than drmin): ',je,prn(n),prn(n+1),rb(ib)
          this%rnod(je+1) = rb(ib)-drmin*0.5
          this%layer(je+1) = ib
          this%state(je+1) = iflso
          if (n+1 < nnod) then                     ! not at surface
             this%rnod(je+2) = rb(ib)+drmin*0.5
             this%layer(je+2) = ib+1 
             this%state(je+2) = nem.iflso.(ib+1)
             je = je+2
             ne = n+1                              ! node n+1 is replaced by rb+drmin
          else                                     ! at surface, create node just below only
             je = je+1
             ne = n
          endif
    !
    !  lower node very close to disco
    !
       else if(dabs(prn(n)-rb(ib)) < drmin) then
    !      print *,'ERN: prn(n) on disco (closer than drmin): ',je,prn(n),prn(n+1),rb(ib)
          this%rnod(je) = rb(ib)-drmin *0.5            ! overwrite already assigned node
          this%layer(je) = ib
          this%state(je) = iflso
          this%rnod(je+1) = rb(ib)+drmin*0.5           ! new node above disco
          this%layer(je+1) = ib+1 
          this%state(je+1) = nem.iflso.(ib+1)
          je = je+1
          ne = n                                       ! node n+1 not replaced
    !
    !  neither node close to disco
    !
       else
    !      print *,'ERN: neither pr-node close to disco: ',je,prn(n),prn(n+1),rb(ib)
          ne = n
          this%rnod(je+1) = rb(ib)-drmin*0.5
          this%layer(je+1) = ib 
          this%state(je+1) = iflso
          this%rnod(je+2) = rb(ib)+drmin*0.5
          this%layer(je+2) = ib+1 
          this%state(je+2) = nem.iflso.(ib+1)
          je = je+2
       endif
    enddo
    this%nnod = je
    this%rnod => reallocate(this%rnod,je)
    this%state => reallocate(this%state,je)
    this%layer => reallocate(this%layer,je)
    end subroutine splitNodesAtDiscosExternalRadialNodes
!------------------------------------------------------------------------
!  Compute density, elastic constants and invQ-factors at nodes
!  for the reference frequency and without attenuation
!
    subroutine elasticConstantsExternalRadialNodes(this,nem,errmsg)
    type (external_radial_nodes) :: this
    type (node_earthmodel) :: nem
    type (error_message) :: errmsg
    integer :: j
!
    allocate(this%ro(this%nnod),this%elcon(NELCON,this%nnod),this%qinv(2,this%nnod))
    call computeSplineEarthmodelCoefficients(nem,.fref.nem,'ELASTIC',errmsg)
    do j = 1,this%nnod
       call evalComplexElasticConstantsSplineEarthmodel(this%rnod(j),this%layer(j),this%ro(j),this%elcon(:,j),this%qinv(:,j))
    enddo
    call deallocSplineEarthmodelCoefficients
    end subroutine elasticConstantsExternalRadialNodes
!------------------------------------------------------------------------
!  Evaluate complex constants for given frequency and node
!
    subroutine getComplexConstantsSelectedExternalRadialNodes(this,fref,f,attmode,jn,zelcon)
    type (external_radial_nodes) :: this
    double precision :: fref,f
    character (len=*) attmode
    integer :: jn
    double complex, dimension(NELCON) :: zelcon
!
    call complexFromElasticConstants(f/fref,attmode,this%elcon(:,jn),this%qinv(:,jn),zelcon)
    end subroutine getComplexConstantsSelectedExternalRadialNodes
!------------------------------------------------------------------------
!  Get total number of nodes
!
    integer function getNnodExternalRadialNodes(this) result(res)
    type (external_radial_nodes), intent(in) :: this
    res = this%nnod
    end function getNnodExternalRadialNodes
!------------------------------------------------------------------------
!  Get rearth
!
    double precision function getRearthExternalRadialNodes(this) result(res)
    type (external_radial_nodes), intent(in) :: this
    res = this%rearth
    end function getRearthExternalRadialNodes
!--------------------------------------------------------------------------
!  Return allocated double pointer to radii
!
    function getDoubleRadiiExternalRadialNodes(this) result(res)
    type (external_radial_nodes), intent(in) :: this
    double precision, dimension(:), pointer :: res
   !
    allocate(res(this%nnod))
    res = this%rnod
    end function getDoubleRadiiExternalRadialNodes
!--------------------------------------------------------------------------
!  Return double precision pointer to radii
!
    function getPointerDoubleRadiiExternalRadialNodes(this) result(res)
    type (external_radial_nodes), intent(in) :: this
    double precision, dimension(:), pointer :: res
    res => this%rnod
    end function getPointerDoubleRadiiExternalRadialNodes
!--------------------------------------------------------------------------
!  Return pointer to state
!
    function getPointerStateExternalRadialNodes(this) result(res)
    type (external_radial_nodes), intent(in) :: this
    integer, dimension(:), pointer :: res
    res => this%state
    end function getPointerStateExternalRadialNodes
!--------------------------------------------------------------------------
!  Return pointer to layer
!
    function getPointerLayerExternalRadialNodes(this) result(res)
    type (external_radial_nodes), intent(in) :: this
    integer, dimension(:), pointer :: res
    res => this%layer
    end function getPointerLayerExternalRadialNodes
!--------------------------------------------------------------------------
!  Return allocated real pointer to radii
!
    function getRealRadiiExternalRadialNodes(this) result(res)
    type (external_radial_nodes), intent(in) :: this
    real, dimension(:), pointer :: res
!
    allocate(res(this%nnod))
    res = real(this%rnod)
    end function getRealRadiiExternalRadialNodes
!-------------------------------------------------------------------
!> \brief Get radius of uppermost node
!
    subroutine getDoubleRadiusTopNodeExternalRadialNodes(this,rtop)
    type (external_radial_nodes), intent(in) :: this
    double precision :: rtop
    rtop = this%rnod(this%nnod)
    end subroutine getDoubleRadiusTopNodeExternalRadialNodes
!-------------------------------------------------------------------
!> \brief Get radius of lowermost node
!
    subroutine getDoubleRadiusBottomNodeExternalRadialNodes(this,rbot)
    type (external_radial_nodes), intent(in) :: this
    double precision :: rbot
    rbot = this%rnod(1)
    end subroutine getDoubleRadiusBottomNodeExternalRadialNodes
!-------------------------------------------------------------------
!> \brief Get radius of node with index j 
!!  where validity of j was checked before
!
    function getDoubleRadiusSelectedExternalRadialNodes(this,j) result(res)
    type (external_radial_nodes), intent(in) :: this
    double precision :: res
    integer :: j
    res = this%rnod(j)
    end function getDoubleRadiusSelectedExternalRadialNodes
!-------------------------------------------------------------------
! Print external nodes to screen
!
    subroutine printExternalRadialNodes(this)
    type (external_radial_nodes) :: this
    integer :: n
    print *,'-------------------------------'
    print *,'External Radial Nodes'
    print *,'-------------------------------'
    print *,'Number of nodes: ',this%nnod
    write(6,'(a10,a15,2a10)') 'node','radius','state','layer'
    do n = 1,this%nnod
       write(6,'(i10,g15.4,2i10)') n,this%rnod(n),this%state(n),this%layer(n)
    enddo
    end subroutine printExternalRadialNodes
!--------------------------------------------------------------------
! get index of node closest to given radius
!
    function getNodeIdxFromRadiusExternalRadialNodes(this,rs) result(jr)
    type (external_radial_nodes) :: this
    double precision :: rs
    integer :: jr
!
    jr = locate(rs,this%nnod,this%rnod)              ! index of node below rs
    if (jr == 0) jr = 1                              ! use deepest node if rs is below it
    if (jr < this%nnod) then                         ! take node closest to rs
       if (abs(rs-this%rnod(jr)) > abs(rs-this%rnod(jr+1))) jr = jr+1
    endif
    end function getNodeIdxFromRadiusExternalRadialNodes
!
end module externalRadialNodes

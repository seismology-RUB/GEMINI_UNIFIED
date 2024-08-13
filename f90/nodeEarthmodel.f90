! ======================================================================
!  Data and functions for an earth model based on nodes and layers
! ======================================================================
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
!------------------------------------------------------------------
!   Module for node earthmodel object
!   Assumes use of SI units for model properties, i.e. density in kg/m^3
!   velocity in m/s and radii in m
!------------------------------------------------------------------
module nodeEarthmodel
    use fson
    use complexElasticConstants
    use cubicSpline
    use realloc
    use mathConstants
    use locatePoint
    use errorMessage
    use string
    implicit none
    interface dealloc
        module procedure deallocNodeEarthmodel
    end interface
    interface operator (.nk.); module procedure getNkNodeEarthmodel; end interface
    interface operator (.nlay.); module procedure getNlayNodeEarthmodel; end interface
    interface operator (.nic.); module procedure getInnerCoreLayerNodeEarthmodel; end interface
    interface operator (.noc.); module procedure getOuterCoreLayerNodeEarthmodel; end interface
    interface operator (.nuc.); module procedure getNlUpperCrustNodeEarthmodel; end interface
    interface operator (.rearth.); module procedure getEarthRadiusNodeEarthmodel; end interface
    interface operator (.ric.); module procedure getInnerCoreRadiusNodeEarthmodel; end interface
    interface operator (.roc.); module procedure getOuterCoreRadiusNodeEarthmodel; end interface
    interface operator (.ruc.); module procedure getOceanBottomRadiusNodeEarthmodel; end interface
    interface operator (.iktop.); module procedure getTopLayerIdxNodeEarthmodel; end interface
    interface operator (.ikbot.); module procedure getBottomLayerIdxNodeEarthmodel; end interface
    interface operator (.rhs.); module procedure getRadiusHalfspaceNodeEarthmodel; end interface
    interface operator (.aniflag.); module procedure getAnisotropyFlagNodeEarthmodel; end interface    
    interface operator (.fsflag.); module procedure getFullspaceFlagNodeEarthmodel; end interface    
    interface operator (.equidisnodes.); module procedure getEquidisnodesFlagNodeEarthmodel; end interface
    interface operator (.isfluid.); module procedure isLayerFluidNodeEarthmodel; end interface      
    interface operator (.iflso.); module procedure iflsoNodeEarthmodel; end interface      
    interface operator (.hasocean.); module procedure hasOceanNodeEarthmodel; end interface      
    interface operator (.hascore.); module procedure hasCoreNodeEarthmodel; end interface      
    interface operator (.fref.); module procedure getReferenceFrequencyNodeEarthmodel; end interface    
    type node_earthmodel
        private
        !----------------------------  FLNM PART ----------------------------------------------------------
        double precision :: fref       ! reference frequency for model paramters (Hz)
        double precision :: rearth     ! radius of model (m)
        double precision :: ric        ! radius of inner core (m)
        double precision :: roc        ! radius of outer core (m)
        integer :: dampflag            ! attenuation mode (no longer used, kept for backwards conpatibility)
        integer :: aniflag             ! 0: model is isotropic; 1: model is transversely isotropic (VTI)
        integer :: full_space_flag     ! 0: free surface, 1: halfspace on top of model
        integer :: equi_dis_nodes      ! 1: equidistant radial nodes per layer
        integer :: nk                  ! number of model knots
        double precision, dimension(:), pointer :: rk             ! radii of model knots, (m)
        double precision, dimension(:), allocatable :: ro         ! density at radial knots, (kg/m^3)
        double precision, dimension(:), allocatable :: vpv        ! vpv at radial knots, (m/s)
        double precision, dimension(:), allocatable :: vph        ! vph at radial knots, (m/s)
        double precision, dimension(:), allocatable :: vsv        ! vsv at radial knots, (m/s)
        double precision, dimension(:), allocatable :: vsh        ! vsh at radial knots, (m/s)
        double precision, dimension(:), allocatable :: eta        ! eta at radial knots
        double precision, dimension(:), allocatable :: qkinv      ! inverse qkappa at radial knots
        double precision, dimension(:), allocatable :: qminv      ! inverse qmue at radial knots (zero in liquids)
        !------------------------------- LAYER PART --------------------------------------------------------
        integer :: nlay                                           ! number of layers
        integer :: nic                                            ! layer index of inner core  (zero if there is no outer core)
        integer :: noc                                            ! layer index of outer core  (zero if there is no outer core)
        integer, dimension(:), allocatable :: iflso               ! array of ones and zeros indicating whther layer is liquid or solid
        integer, dimension(:), pointer :: iktop                   ! index of uppermost knot in layer
        integer, dimension(:), allocatable :: ikbot               ! index of lowermost knot in layer (ikbot(1) = 0)
        double precision, dimension(:), pointer :: rb             ! radii of layer boundaries (rb(n) is at top of n-th layer)
        double precision, dimension(:), pointer :: drk            ! if equidistant nodes: dr per layer, else -1.d0
    end type node_earthmodel
!
 contains
!------------------------------------------------------------------
!> \brief Read in model parameters from JSON file
!! lu:        Fortran file id
!! filename:  name of node model file
!! errmsg:    error message object
!
    subroutine createNodeEarthmodel(this,lu,filename,errmsg)
    type (node_earthmodel) :: this
    integer :: lu
    character (len=*) :: filename
    type (error_message) :: errmsg
    type (fson_value), pointer :: contents,nodeprop
    character (len=20) :: myname = 'createNodeEarthmodel'
    double precision, dimension(:), allocatable :: vals1d
    integer, dimension(:), allocatable :: ivals1d
    integer :: nl
    !
    call addTrace(errmsg,myname)
    nullify(contents)
    contents => fson_parse(filename,unit = lu)
  !
  !  extract some variables from file
  !
    call fson_get(contents,'reference_frequency_hz',this%fref)
    call fson_get(contents,'earth_radius',this%rearth)
    call fson_get(contents,'anisotropic',this%aniflag)
    call fson_get(contents,'full_space',this%full_space_flag)
    call fson_get(contents,'equidistant_nodes',this%equi_dis_nodes)
    call fson_get(contents,'num_layers',this%nlay)
    call fson_get(contents,'num_nodes',this%nk)
    call fson_get(contents,'layer_inner_core',this%nic)
    call fson_get(contents,'layer_outer_core',this%noc)
    call fson_get(contents,'radius_inner_core',this%ric)
    call fson_get(contents,'radius_outer_core',this%roc)
  !
  !  fluid flags
  !
    call fson_get(contents,'isfluid',ivals1d)
    allocate(this%iflso(this%nlay))
    this%iflso = ivals1d
    deallocate(ivals1d)
  !
  !  inhibit fullspace and ocean
  !
    if (this%iflso(this%nlay) == 1 .and. this%full_space_flag == 1) then
       call add(errmsg,2,'Full space and ocean is not allowed!!',myname)
       call fson_destroy(contents)
       deallocate(this%iflso)
       return
    endif   
  !
  !  upper layer boundaries
  !
    call fson_get(contents,'radius_upper_layer_boundary_m',vals1d)
    allocate(this%rb(this%nlay))
    this%rb = vals1d
    deallocate(vals1d)
  !
  !  upper layer boundary node index
  !
    call fson_get(contents,'index_upper_boundary',ivals1d)
    allocate(this%iktop(this%nlay))
    this%iktop = ivals1d
    allocate(this%ikbot(this%nlay))
    this%ikbot(1) = 0
    this%ikbot(2:this%nlay) = ivals1d(1:this%nlay-1)+1
    deallocate(ivals1d)
  !
  !  node radii
  !
    call fson_get(contents,'node_properties',nodeprop)
    call fson_get(nodeprop,'depth',vals1d)
    allocate(this%rk(this%nk))
    this%rk = this%rearth-vals1d
    deallocate(vals1d)
  !
  !  distance between nodes per layer (leave out halfspace layer)
  !  drk = (rk(top)-rk(bot))/nval with nval = iktop-ikbot (!)
  !  set to -1.d0 if nodes are not equidistant
  !
    allocate(this%drk(this%nlay))
    if (this%equi_dis_nodes == 1) then
       this%drk(1) = -1.d0
       forall (nl = 2:this%nlay) &
            & this%drk(nl) = (this%rk(this%iktop(nl))-this%rk(this%ikbot(nl)))/(this%iktop(nl)-this%ikbot(nl))
    else
       this%drk = -1.d0
    endif
  !
  !  properties at nodes
  !
  !  density
  !
    call fson_get(nodeprop,'density',vals1d)
    allocate(this%ro(this%nk))
    this%ro = vals1d
    deallocate(vals1d)
  !
  !  Qkappa
  !
    call fson_get(nodeprop,'Qkappa',vals1d)
    allocate(this%qkinv(this%nk))
    this%qkinv = 1./vals1d
    deallocate(vals1d)
  !
  !  Qmu
  !
    call fson_get(nodeprop,'Qmu',vals1d)
    allocate(this%qminv(this%nk))
    this%qminv = 1./vals1d
    where (this%qminv < 0.d0) this%qminv = 0.d0
    deallocate(vals1d)
    !
    if (this%aniflag == 0) then
     !
     !  P-velocity
     !
       call fson_get(nodeprop,'P-velocity',vals1d)
       allocate(this%vpv(this%nk),this%vph(this%nk))
       this%vpv = vals1d
       this%vph = vals1d
       deallocate(vals1d)
     !
     !  S-velocity
     !
       call fson_get(nodeprop,'S-velocity',vals1d)
       allocate(this%vsv(this%nk),this%vsh(this%nk))
       this%vsv = vals1d
       this%vsh = vals1d
       deallocate(vals1d)
     !
     !  Eta
     !
       allocate(this%eta(this%nk))
       this%eta = 1.0
    else
     !
     !  PV- and PH-velocity
     !
       call fson_get(nodeprop,'PV-velocity',vals1d)
       allocate(this%vpv(this%nk),this%vph(this%nk))
       this%vpv = vals1d
       deallocate(vals1d)
       call fson_get(nodeprop,'PH-velocity',vals1d)
       this%vph = vals1d
       deallocate(vals1d)
     !
     !  SV- and SH-velocity
     !
       call fson_get(nodeprop,'SV-velocity',vals1d)
       allocate(this%vsv(this%nk),this%vsh(this%nk))
       this%vsv = vals1d
       deallocate(vals1d)
       call fson_get(nodeprop,'SH-velocity',vals1d)
       this%vsh = vals1d
       deallocate(vals1d)
     !
     !  Eta
     !
       call fson_get(nodeprop,'eta',vals1d)
       allocate(this%eta(this%nk))
       this%eta = vals1d
       deallocate(vals1d)
    endif
    call fson_destroy(contents)
    end subroutine createNodeEarthmodel
!-----------------------------------------------------------------
!> \brief Deallocate node earthmodel object
!> \param this node_earthmodel object
!
    subroutine deallocNodeEarthmodel(this)
    type (node_earthmodel) :: this
    if(associated(this%rk)) deallocate(this%rk) 
    if(allocated(this%ro)) deallocate(this%ro) 
    if(allocated(this%vpv)) deallocate(this%vpv)
    if(allocated(this%vph)) deallocate(this%vph) 
    if(allocated(this%vsv)) deallocate(this%vsv) 
    if(allocated(this%vsh)) deallocate(this%vsh)
    if(allocated(this%eta)) deallocate(this%eta) 
    if(allocated(this%qkinv)) deallocate(this%qkinv) 
    if(allocated(this%qminv)) deallocate(this%qminv) 
    if(associated(this%rb)) deallocate(this%rb)
    if(allocated(this%iflso)) deallocate(this%iflso)    
    if(associated(this%iktop)) deallocate(this%iktop)    
    if(allocated(this%ikbot)) deallocate(this%ikbot)
    if(associated(this%drk)) deallocate(this%drk)    
    end subroutine deallocNodeEarthmodel
!------------------------------------------------------------------
!> \brief Compute complex elastic constants for frequency f at radial knots
!!  attmode:    attenuation mode (in)
!!  f:          frequency (in)
!!  ro:         density (out)
!!  zelcon:     elastic constants A,C,F,L,N,Kappa,Mu at knots (out)
!!  errmsg:     error message object
!
    subroutine complexElasticConstantsAtKnotsNodeEarthmodel(this,attmode,f,ro,qinv,zelcon,errmsg)
    type (node_earthmodel) :: this
    character(len=*) :: attmode
    double precision :: f
    double precision, dimension(:) :: ro
    double precision, dimension(:,:) :: qinv
    double complex, dimension(:,:) :: zelcon
    type (error_message) :: errmsg
    character(len=44) :: myname = 'complexElasticConstantsAtKnotsNodeEarthmodel'
    integer :: j
!
    call addTrace(errmsg,myname)
    if (size(ro) /= this%nk .or. size(qinv,2) /= this%nk .or. size(zelcon,2) /= this%nk) then
       call add(errmsg,2,'Inconsistent dimension of model parameters',myname)
       return
    endif
    if (attmode.equal.'ELASTIC') then
       do j = 1,this%nk
          call complexElasticVTIConstants(1.d0,this%ro(j),this%vpv(j),this%vph(j),this%vsv(j),this%vsh(j),this%eta(j),&
               & 0.d0,0.d0,zelcon(1,j),zelcon(2,j),zelcon(3,j),zelcon(4,j),zelcon(5,j),zelcon(6,j),zelcon(7,j))
       enddo
    else if (attmode.equal.'ATTENUATION_ONLY') then
       do j = 1,this%nk
          call complexElasticVTIConstants(1.d0,this%ro(j),this%vpv(j),this%vph(j),this%vsv(j),this%vsh(j),this%eta(j),&
               & this%qkinv(j),this%qminv(j),zelcon(1,j),zelcon(2,j),zelcon(3,j),zelcon(4,j),zelcon(5,j),zelcon(6,j),zelcon(7,j))
       enddo
    else if (attmode.equal.'DISPERSION_ONLY') then
       do j = 1,this%nk
          call complexElasticVTIConstants(f/this%fref,this%ro(j),this%vpv(j),this%vph(j),this%vsv(j),this%vsh(j),this%eta(j),&
               & 0.d0,0.d0,zelcon(1,j),zelcon(2,j),zelcon(3,j),zelcon(4,j),zelcon(5,j),zelcon(6,j),zelcon(7,j))
       enddo
    else if (attmode.equal.'ATTENUATION_AND_DISPERSION') then
       do j = 1,this%nk
          call complexElasticVTIConstants(f/this%fref,this%ro(j),this%vpv(j),this%vph(j),this%vsv(j),this%vsh(j),this%eta(j),&
               & this%qkinv(j),this%qminv(j),zelcon(1,j),zelcon(2,j),zelcon(3,j),zelcon(4,j),zelcon(5,j),zelcon(6,j),zelcon(7,j))
       enddo
    else
       call add(errmsg,2,'Invalid attenuation mode specified',myname)
       return
    endif
    ro = this%ro
    qinv(1,:) = this%qkinv
    qinv(2,:) = this%qminv
    end subroutine complexElasticConstantsAtKnotsNodeEarthmodel
!---------------------------------------------------------------
!> \brief Extract node radii from isotropic node model
!> \param this node_earthmodel object
!> \param nk restrict to top nk nodes
!
    function getRkArrayNodeEarthmodel(this,nk) result(rk)
    type (node_earthmodel) :: this
    integer, optional :: nk
    double precision, dimension(:), pointer :: rk
!
    if(.not.present(nk)) then
        rk => this%rk(1:this%nk)
    else
        nk=min(nk,this%nk)
        rk => this%rk(this%nk-nk+1:this%nk)
    endif
    end function getRkArrayNodeEarthmodel
!---------------------------------------------------------------
!> \brief Extract iktop array from isotropic node model
!> \param this node_earthmodel object
!
    function getIktopArrayNodeEarthmodel(this) result(res)
    type (node_earthmodel) :: this
    integer, dimension(:), pointer :: res
!
    res => this%iktop
    end function getIktopArrayNodeEarthmodel
!---------------------------------------------------------------
!> \brief Extract drk array from isotropic node model
!> \param this node_earthmodel object
!
    function getDrkArrayNodeEarthmodel(this) result(res)
    type (node_earthmodel) :: this
    double precision, dimension(:), pointer :: res
!
    res => this%drk
    end function getDrkArrayNodeEarthmodel
!----------------------------------------------------------------
!> \brief Query for anisotropy
!> \param this node_earthmodel object
!
    logical function isAnisotropicNodeEarthmodel(this)
    type (node_earthmodel) :: this
    isAnisotropicNodeEarthmodel = .false.
    if(this%aniflag == 1) isAnisotropicNodeEarthmodel = .true.
    end function isAnisotropicNodeEarthmodel
!---------------------------------------------------------------
!> \brief get anisotropy flag
!
    function getAnisotropyFlagNodeEarthmodel(this) result(res)
    type (node_earthmodel), intent(in) :: this
    integer :: res
    res = this%aniflag
    end function getAnisotropyFlagNodeEarthmodel
!----------------------------------------------------------------
!> \brief Query for full space
!> \param this node_earthmodel object
!
    function getFullspaceFlagNodeEarthmodel(this) result(res)
    type (node_earthmodel), intent(in) :: this
    integer :: res
    res = this%full_space_flag
    end function getFullspaceFlagNodeEarthmodel
!----------------------------------------------------------------
!> \brief Query for equidistant nodes flag
!> \param this node_earthmodel object
!
    function getEquidisnodesFlagNodeEarthmodel(this) result(res)
    type (node_earthmodel), intent(in) :: this
    integer :: res
    res = this%equi_dis_nodes
    end function getEquidisnodesFlagNodeEarthmodel
!---------------------------------------------------------------
!> \brief get reference frequency
!
    function getReferenceFrequencyNodeEarthmodel(this) result(res)
    type (node_earthmodel), intent(in) :: this
    double precision :: res
    res = this%fref
    end function getReferenceFrequencyNodeEarthmodel
!---------------------------------------------------------------
!> \brief Extract discontinuity radii from isotropic node model
!> \param this node_earthmodel object
!> \param nd restrict to top nd interfaces
!
    function getRbArrayNodeEarthmodel(this,nd) result(rb)
    type (node_earthmodel) :: this
    integer, optional :: nd
    double precision, dimension(:), pointer :: rb
!
    if(.not.present(nd)) then
        rb => this%rb(1:this%nlay)
    else
        rb => this%rb(this%nlay-nd+1:this%nlay)
    endif
    end function getRbArrayNodeEarthmodel
!---------------------------------------------------------------------------
!> \brief Get earth radius
!
    double precision function getEarthRadiusNodeEarthmodel(this)
    type (node_earthmodel), intent(in) :: this
    getEarthRadiusNodeEarthmodel = this%rearth
    end function getEarthRadiusNodeEarthmodel
!---------------------------------------------------------------------------
!> \brief Get radius of top of halfspace
!
    double precision function getRadiusHalfspaceNodeEarthmodel(this)
    type (node_earthmodel), intent(in) :: this
    getRadiusHalfspaceNodeEarthmodel = this%rb(1)
    end function getRadiusHalfspaceNodeEarthmodel
!---------------------------------------------------------------------------
!> \brief Get radius of inner core
!
    function getInnerCoreRadiusNodeEarthmodel(this) result(res)
    type (node_earthmodel), intent(in) :: this
    double precision :: res
    res = this%ric
    end function getInnerCoreRadiusNodeEarthmodel
!---------------------------------------------------------------------------
!> \brief Get radius of outer core
!
    function getOuterCoreRadiusNodeEarthmodel(this) result(res)
    type (node_earthmodel), intent(in) :: this
    double precision :: res
    res = this%roc
    end function getOuterCoreRadiusNodeEarthmodel
!---------------------------------------------------------------------------
!> \brief Get radius of ocean bottom
!
    function getOceanBottomRadiusNodeEarthmodel(this) result(res)
    type (node_earthmodel), intent(in) :: this
    double precision :: res
    if (this%iflso(this%nlay) == 1) then
       res = this%rb(this%nlay-1)
    else
       res = this%rearth
    endif
    end function getOceanBottomRadiusNodeEarthmodel
!---------------------------------------------------------------------------
!> \brief Get number of nodes
!
    integer function getNkNodeEarthmodel(this)
    type (node_earthmodel), intent(in) :: this
    getNkNodeEarthmodel = this%nk
    end function getNkNodeEarthmodel
!---------------------------------------------------------------------------
!> \brief Get number of layers
!
    integer function getNlayNodeEarthmodel(this)
    type (node_earthmodel), intent(in) :: this
    getNlayNodeEarthmodel = this%nlay
    end function getNlayNodeEarthmodel
!---------------------------------------------------------------------------
!> \brief Get layer index of uppermost crustal layer
!
    function getNlUpperCrustNodeEarthmodel(this) result(res)
    type (node_earthmodel), intent(in) :: this
    integer :: res    
    res = this%nlay
    if (this%iflso(this%nlay) == 1) then
       res = this%nlay -1
    endif
    end function getNlUpperCrustNodeEarthmodel
!---------------------------------------------------------------------------
!> \brief Get index of top node of given layer
!
    integer function getTopLayerIdxNodeEarthmodel(this,nl)
    type (node_earthmodel), intent(in) :: this
    integer, intent(in) :: nl
    getTopLayerIdxNodeEarthmodel = this%iktop(nl)
    end function getTopLayerIdxNodeEarthmodel
!---------------------------------------------------------------------------
!> \brief Get index of bottom node of given layer
!
    integer function getBottomLayerIdxNodeEarthmodel(this,nl)
    type (node_earthmodel), intent(in) :: this
    integer, intent(in) :: nl
    getBottomLayerIdxNodeEarthmodel = this%ikbot(nl)
    end function getBottomLayerIdxNodeEarthmodel
!---------------------------------------------------------------------------
!> \brief Get index of inner core layer
!
    function getInnerCoreLayerNodeEarthmodel(this) result(res)
    type (node_earthmodel), intent(in) :: this
    integer :: res
    res = this%nic
    end function getInnerCoreLayerNodeEarthmodel
!---------------------------------------------------------------------------
!> \brief Get index of outer core layer
!
    function getOuterCoreLayerNodeEarthmodel(this) result(res)
    type (node_earthmodel), intent(in) :: this
    integer :: res
    res = this%noc
    end function getOuterCoreLayerNodeEarthmodel
!---------------------------------------------------------------------------
!> \brief Is layer fluid?
!
    function isLayerFluidNodeEarthmodel(this,nl) result(res)
    type (node_earthmodel), intent(in) :: this
    integer, intent(in) :: nl    
    logical :: res
    if (this%iflso(nl) == 1) then
       res = .true.
    else
       res = .false.
    endif
    end function isLayerFluidNodeEarthmodel
!---------------------------------------------------------------------------
!> \brief return fluid/solid integer code
!
    function iflsoNodeEarthmodel(this,nl) result(res)
    type (node_earthmodel), intent(in) :: this
    integer, intent(in) :: nl    
    integer :: res
    res = this%iflso(nl)
    end function iflsoNodeEarthmodel
!---------------------------------------------------------------------------
!> \brief Has model a liquid core?
!
    function hasCoreNodeEarthmodel(this) result(res)
    type (node_earthmodel), intent(in) :: this
    logical :: res
    if (this%noc > 0 .and. this%iflso(this%noc) == 1) then
       res = .true.
    else
       res = .false.
    endif
    end function hasCoreNodeEarthmodel
!---------------------------------------------------------------------------
!> \brief Has model an ocean?
!
    function hasOceanNodeEarthmodel(this) result(res)
    type (node_earthmodel), intent(in) :: this
    logical :: res
    if (this%iflso(this%nlay) == 1) then
       res = .true.
    else
       res = .false.
    endif
    end function hasOceanNodeEarthmodel
!---------------------------------------------------------------------------
!> \brief Get some properties at given node index
!
    function getPropertyNodeEarthmodel(this,prop,i) result(res)
    type (node_earthmodel), intent(in) :: this
    character (len=*), intent(in) :: prop
    integer, intent(in) :: i
    real :: res
    if (prop.equal.'node_radius') then
       res = this%rk(i)
    else if (prop.equal.'vp') then
       res = this%vpv(i)
    else if (prop.equal.'vs') then
       res = this%vsv(i)
    else if (prop.equal.'rho') then
       res = this%ro(i)
    else
       res = 999999.
    endif
    end function getPropertyNodeEarthmodel
!-----------------------------------------------------------------------------
!> \brief Find layer containing given radius
!> \param this node_earthmodel object
!> \param r radius (input)
!> \param nl Layer index (output)
!> \param top Radius is equal to top boundary of layer (output)
!
    subroutine getLayerIndexNodeEarthmodel(this,r,nl,top)
    type (node_earthmodel) :: this
    double precision :: r
    integer :: nl,top
!
    top = 0
    nl = locate(r,this%nlay,this%rb)+1
    if (dabs(r/this%rb(nl)-1.d0).lt.epsilon(1.d0)) top = 1
    end subroutine getLayerIndexNodeEarthmodel
!------------------------------------------------------------------------------
!  Find top radius and index of layer just below next phase change
!  above current layer  
!
!  nl:      index of current layer
!  rpc:     top radius of layer just below phase change (out)
!  npc:     index of layer just below phase change (out)
!
    function findAbovePhaseChangeNodeEarthmodel(this,nl) result(rpc)
    type (node_earthmodel) :: this
    integer :: nl,npc
    double precision :: rpc
!
    if (this%iflso(nl) == 1) then                  ! fluid to solid
       npc = nl                                    ! I am either in core or ocean
    else                                           ! solid to fluid
       if (nl < this%noc) then                     ! I am below outer core
          npc = this%noc-1                         ! next pc is CMB
       else if (nl > this%noc) then                ! I am above outer core
          if (this%iflso(this%nlay) == 1) then     ! model has ocean
             npc = this%nlay-1                     ! next pc is ocean
          else
             npc = this%nlay                       ! next pc is surface
          endif
       else
          npc = this%nlay                          ! should never be reached, avoid compiler warning
       endif
    endif
    rpc = this%rb(npc)
    end function findAbovePhaseChangeNodeEarthmodel
!------------------------------------------------------------------------------
!  Find bottom radius and index of layer just above next phase change
!  below current layer  
!
!  nl:      index of current layer
!  rpc:     bottom radius of layer just above phase change (out)
!  npc:     index of layer just above phase change (out)
!
    function findBelowPhaseChangeNodeEarthmodel(this,nl) result(rpc)
    type (node_earthmodel) :: this
    integer :: nl,npc
    double precision :: rpc
!
    if (this%iflso(nl) == 1) then                  ! fluid to solid
       npc = nl                                    ! I am either in core or ocean
    else                                           ! solid to fluid
       if (nl < this%noc) then                     ! I am below outer core
          npc = 1                                  ! next pc is center of the earth
       else if (nl > this%noc) then                ! I am above outer core
          npc = this%noc+1                         ! next pc is CMB
       else
          npc = nl                                 ! should never be reached, avoid compiler warning
       endif
    endif
    rpc = this%rb(npc-1)
    end function findBelowPhaseChangeNodeEarthmodel
!----------------------------------------------------------------
! Print properties of node earth model to screen
!
    subroutine printNodeEarthmodel(this)
    type (node_earthmodel) :: this
    integer :: n
    print *,'---------------------------------------'
    print *,'Node earth model properties'
    print *,'---------------------------------------'
    print *,'Number of nodes: ',this%nk
    print *,'Number of layers: ', this%nlay
    print *,'Layer of inner core: ',this%nic
    print *,'Layer of outer core: ',this%noc
    print *,'Layer structure'
    write(6,'(4a10,a15)') 'layer','fluid','bottom','top','top radius'
    do n = 1,this%nlay
       write(6,'(4i10,g15.4)') n,this%iflso(n),this%ikbot(n),this%iktop(n),this%rb(n)
    enddo
    print *, 'Node structure'
    write(6,'(a10,a15,6a12,2a15)') 'node','radius','rho','vpv','vph','vsv','vsh','eta','qkinv','qminv'
    do n = 1,this%nk
       write(6,'(i10,g15.8,6f12.3,2g15.4)') n,this%rk(n),this%ro(n),this%vpv(n),this%vph(n),&
            this%vsv(n),this%vsh(n),this%eta(n),this%qkinv(n),this%qminv(n)
    enddo               
    end subroutine printNodeEarthmodel
!
 end module nodeEarthmodel

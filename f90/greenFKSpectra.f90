! ====================================================================================
!  Handling of Green functions in the fk-domain
! ====================================================================================
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
!----------------------------------------------------------------------
!  This module provides a data structure to handle Green fk-spectra,
!  provides routines for reading spectra and meta data from file,
!  and other service routines involved.
!----------------------------------------------------------------------
 module greenFKSpectra
        use hdf5
    use errorMessage
    use locatePoint
    use string
    use mathConstants
    use hdfWrapper
    use anyRankRealArray
    use anyRankIntegerArray
    use externalRadialNodes
    implicit none
    interface dealloc; module procedure deallocGreenFKSpectra; end interface
    type green_fk_spectra
       type (external_radial_nodes) :: exnod                               ! external radial nodes instance
       integer :: nf1                                                      !< index of first frequency
       integer :: nf2                                                      !< index of last frequency
       integer :: nkfmax                                                   !< size of frequency-wavenumber spectra
       integer :: nwnmax                                                   !< max number of wavenumbers
       integer :: istyp                                                    !< source type (0 = force, 1 = moment)
       integer :: ncsph, nctor                                             !< number of components spheroidal and toroidal
       integer :: njsph,njtor                                              !< number of source jumps spheroidal, toroidal
       integer :: numtasks                                                 !< number of processes by which fk-spectrum was calculated
       integer :: global                                                   !< = 1 for global calculation, = 0 else
       integer, dimension(9) :: dsvmask                                    !< mask of available DSV-components (U,R,V,S,UP,VP,W,T,WP)
       integer, dimension(:), allocatable :: nwn                           !< number of wavenumbers per frequency array
       double precision :: df                                              !< frequency spacing
       double precision :: dwn                                             !< wavenumber spacing (relevant for non-global)
       double precision :: sigma                                           !< imaginary part of omega
       double precision :: rearth                                          !< earth radius
       double precision :: rse                                             !< receiver radius OR source radius, needed outside of this module
       double precision :: fref                                            !< reference frequency of earth model
       double precision, dimension(:), allocatable :: rsnod                !< source node radii
       character (len=:), allocatable :: earth_model_name                  !< name of earth model used for calculation
       character (len=:), allocatable :: attnmode                          !< attenuation mode used in calculation
       character (len=:), allocatable :: identity                          !< unique string providing identity to meta data
    end type green_fk_spectra
!
!  Mapping of component names to component index
!
    character (len=2), dimension(9) :: names_green_fk_spectra = (/' U',' R',' V',' S','UP','VP',' W',' T','WP' /) 
!
    contains
!-----------------------------------------------------------------------
!  Read meta data from specific path from HDF locid
!
    subroutine readMetaGreenFKSpectra(this,locid,path,errmsg)
    type (green_fk_spectra) :: this
    integer(hid_t) :: locid
    character (len=*) :: path
    type (error_message) :: errmsg
    integer(hid_t) :: grpid
    type (any_rank_real_array) :: arra
    type (any_rank_integer_array) :: aria
    integer, dimension(:), pointer :: id
    real, dimension(:), pointer :: d
    integer :: dsvstep,derivflag,strlen,ierr
    character (len=max_length_string) :: cval
    character (len=23) :: myname = 'readMetaGreenFKSpectra'
!    
    if (path.equal.'identity') then
       call readStringAttributeHDFWrapper(locid,'identity',cval,strlen,errmsg)     
       this%identity = cval(1:strlen)
    else if (path.equal.'earthModelName') then
       call readStringHDFWrapper(locid,'earthModelName',cval,strlen,errmsg)
       this%earth_model_name = cval(1:strlen)
    else if (path.equal.'attenuationMode') then
       call readStringHDFWrapper(locid,'attenuationMode',cval,strlen,errmsg)
       this%attnmode = cval(1:strlen)
    else if (path.equal.'reals') then
       call readArrayHDFWrapper(locid,'reals',arra,errmsg)
       if (.level.errmsg == 2) return
       d => arra%get1d()
       this%rearth = d(1); this%rse = d(2); this%sigma = d(3)
       this%df = d(4); this%dwn = d(5); this%fref = d(6)
       deallocate(d)
    else if (path.equal.'integers') then
       call readArrayHDFWrapper(locid,'integers',aria,errmsg)
       if (.level.errmsg == 2) return
       id => aria%get1d()
       this%nf1 = id(1); this%nf2 = id(2)
       this%nkfmax = id(3); this%nwnmax = id(4)
       dsvstep = id(5); derivflag = id(6);
       this%numtasks = id(7); this%global = id(8)
       this%dsvmask = (/ 1,0,1,0,0,0,1,0,0 /)
       if (derivflag == 1) this%dsvmask((/5,6,9/)) = 1
       if (dsvstep == 1) this%dsvmask((/2,4,8/)) = 1
       deallocate(id)
    else if (path.equal.'numberOfWavenumbersForFrequency') then
       call readArrayHDFWrapper(locid,'numberOfWavenumbersForFrequency',aria,errmsg)
       if (.level.errmsg == 2) return
       this%nwn = aria%get1d(); call aria%dealloc()
    else if (path.equal.'externalNodes') then
       call h5gopen_f(locid,'externalNodes',grpid,ierr)
       call readExternalRadialNodes(this%exnod,grpid,errmsg)
       call h5gclose_f(grpid,ierr)
    else if (path.equal.'sourceNodeRadii') then
       call readArrayHDFWrapper(locid,'sourceNodeRadii',arra,errmsg)
       d => arra%get1d(); allocate(this%rsnod(size(d)))
       this%rsnod = d; deallocate(d)
    else if (path.equal.'dataSpecificIntegers') then
       call readArrayAttributeHDFWrapper(locid,'dataSpecificIntegers',aria,errmsg)
       if (.level.errmsg == 2) return
       id => aria%get1d()
       this%istyp = id(1); this%ncsph = id(2); this%nctor = id(3)
       this%njsph = id(4); this%njtor = id(5)
       deallocate(id)
       if (this%ncsph == 2) this%dsvmask((/4/)) = 0
       if (this%nctor == 0) this%dsvmask((/7,8,9/)) = 0
    else
       call add(errmsg,2,'Path_'+path+' does not exist in HDF file',myname)
    endif
    end subroutine readMetaGreenFKSpectra
!-----------------------------------------------------------------------
!  Read out wavenumber spectra for all frequencies and fixed radial node
!  Works for sequential and parallel run
!  spheroidal: (DSV = 1,2,3,4,5,6 and JUMP = 1,2,3,4)
!  toroidal:   (DSV = 1,2,3   and JUMP = 1,2)
!  jr:     external node index
!  errmsg: error message
!
    subroutine readDataNodeGreenFKSpectra(this,locid,jr,gfkr,gfki,errmsg)
    type (green_fk_spectra) :: this
    integer(hid_t) :: locid
    integer :: jr
    real, dimension(:,:), pointer :: gfkr,gfki
    type (error_message) :: errmsg
    integer :: nsp
    integer(hsize_t), dimension(3) :: offset,count
    integer(hsize_t), dimension(2) :: dimsslab
    type(any_rank_real_array) :: ara
!
    nsp = sum(this%dsvmask(1:6))*this%njsph+sum(this%dsvmask(7:9))*this%njtor
    offset = (/0,0,jr-1/)                                                             ! offset of subset to be read
    count = (/this%nkfmax,nsp,1/)                                                     ! count of data to be set
    dimsslab = (/this%nkfmax,nsp/)                                                    ! dimensions of data to be read into memory
    call readArrayHDFWrapper(locid,'GreenFKSpectraReal',ara,errmsg,&
         dimsslab = dimsslab,offset = offset,count = count)
    gfkr => ara%get2d(); call ara%deassoc()
    call readArrayHDFWrapper(locid,'GreenFKSpectraImag',ara,errmsg,&
         dimsslab = dimsslab,offset = offset,count = count)
    gfki => ara%get2d(); call ara%deassoc()
    end subroutine readDataNodeGreenFKSpectra
!------------------------------------------------------------------------
!  check identity
!
    function checkIdentityGreenFKSpectra(this,locid,errmsg) result(res)
    type (green_fk_spectra) :: this
    integer(hid_t) :: locid
    logical :: res
    type (error_message) :: errmsg
    character (len=:), allocatable :: identity
    character (len=max_length_string) :: cval
    integer :: idlen
!
    call readStringAttributeHDFWrapper(locid,'identity',cval,idlen,errmsg)
    identity = cval(1:idlen)
    if (identity.equal.this%identity) then
       res = .true.
    else
       res = .false.
    endif
    end function checkIdentityGreenFKSpectra
!------------------------------------------------------------------------
!  deallocate and close file
!
    subroutine deallocGreenFKSpectra(this)
    type (green_fk_spectra) :: this
    call dealloc(this%exnod)
    if (allocated(this%rsnod)) deallocate(this%rsnod)
    if (allocated(this%nwn)) deallocate(this%nwn)
    if (allocated(this%earth_model_name)) deallocate(this%earth_model_name)
    if (allocated(this%identity)) deallocate(this%identity)
    if (allocated(this%attnmode)) deallocate(this%attnmode)
    end subroutine deallocGreenFKSpectra
!------------------------------------------------------------------------------------
!  Calculate position of wavenumber spectrum dataset in subgroup from
!  component index and jump index
!
    function positionDataGreenFKSpectra(this,isp,jsp) result(offset)
    type (green_fk_spectra) :: this
    integer :: isp,jsp,offset
    if (isp .le. 6) then
       offset = sum(this%dsvmask(1:isp-1))*this%njsph+jsp
    else
       offset = sum(this%dsvmask(1:6))*this%njsph+sum(this%dsvmask(7:isp-1))*this%njtor+jsp
    endif
    end function positionDataGreenFKSpectra
!------------------------------------------------------------------------------------
!  Calculate number of wavenumber spectra based on available components and jumps
!
    function numberGreenFKSpectra(this) result(nsp)
    type (green_fk_spectra) :: this
    integer :: nsp
    nsp = sum(this%dsvmask(1:6))*this%njsph+sum(this%dsvmask(7:9))*this%njtor
    end function numberGreenFKSpectra
!------------------------------------------------------------------------------------
!  Does spectrum for given component and jump exist?
!
    function existGreenFKSpectra(this,isp,jsp) result(res)
    type (green_fk_spectra) :: this
    logical :: res
    integer :: isp,jsp
    res = .false.
    select case (isp)
    case (1:6)
       res = (this%dsvmask(isp) == 1) .and. (jsp .le. this%njsph) .and. (jsp > 0)
    case (7:9)
       res = (this%dsvmask(isp) == 1) .and. (jsp .le. this%njtor) .and. (jsp > 0)
    end select
    end function existGreenFKSpectra
!--------------------------------------------------------------------
! get node index from radius (in km)
!
    function getNodeIndexFromRadiusGreenFKSpectra(this,rs) result(jr)
    type (green_fk_spectra) :: this
    double precision :: rs
    integer :: jr,nnod
    double precision, dimension(:), pointer :: rnod
!
    rnod => getPointerDoubleRadiiExternalRadialNodes(this%exnod)
    nnod = .nnod.(this%exnod)
    jr = locate(rs,nnod,rnod)                        ! index of node below rs
    if (jr == 0) jr = 1                              ! use deepest node if rs is below it
    if (jr < nnod) then                              ! take node closest to rs
       if (abs(rs-rnod(jr)) > abs(rs-rnod(jr+1))) jr = jr+1
    endif
    end function getNodeIndexFromRadiusGreenFKSpectra
!--------------------------------------------------------------------
!  print table of available components and jumps
!
    subroutine printContentsGreenFKSpectra(this)
    type (green_fk_spectra) :: this
    integer :: isp,jsp
    write(6,'(a)') 'Table of contents:'
    write(6,'(a6,4a8)') 'Name','Comp','Jump','Count'
    do isp = 1,6
       do jsp = 1,this%njsph
          if (this%dsvmask(isp) == 1) then
             write(6,'(a6,3i8)') names_green_fk_spectra(isp),isp,jsp,&
                               & positionDataGreenFKSpectra(this,isp,jsp)
          endif
       enddo
    enddo
    do isp = 7,9
       do jsp = 1,this%njtor
          if (this%dsvmask(isp) == 1) then
             write(6,'(a6,3i8)') names_green_fk_spectra(isp),isp,jsp,&
                               & positionDataGreenFKSpectra(this,isp,jsp)
          endif
       enddo
    enddo
    end subroutine printContentsGreenFKSpectra
!--------------------------------------------------------------------
!  print values of member variables for information
!
    subroutine printMembersGreenFKSpectra(this)
    type (green_fk_spectra) :: this
    double precision :: rearth
!
    rearth = getRearthExternalRadialNodes(this%exnod)
    write(6,'(a,f15.3)') 'Receiver/Source node radius: ',this%rse
    write(6,'(5a15)') 'source type','fmin','fmax','df','dwn'
    write(6,'(i15,4f15.5)') this%istyp,(this%nf1-1)*this%df,(this%nf2-1)*this%df,this%df,this%dwn
    write(6,'(4a15)') 'n-Sph-Comp','n-Tor-Comps','n-Sph-Jumps','n-Tor-Jumps'
    write(6,'(4i15)') this%ncsph,this%nctor,this%njsph,this%njtor
    write(6,'(9a6)') 'U','R','V','S','UP','VP','W','T','WP'
    write(6,'(9i6)') this%dsvmask
    write(6,'(3a15)') 'numtasks','sigma','rearth'
    write(6,'(i15,2f15.3)') this%numtasks,this%sigma,this%rearth
    end subroutine printMembersGreenFKSpectra
!
end module greenFKSpectra

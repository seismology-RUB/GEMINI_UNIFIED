!--------------------------------------------------------------------------
!   Copyright 2016 Wolfgang Friederich (Ruhr-Universitaet Bochum, Germany)
!
!   This file is part of Gemini Unified
!
!   Gemini Unified is free software: you can
!   redistribute it and/or modify it under the terms of the GNU
!   General Public License as published by the Free Software
!   Foundation, either version 2 of the License, or (at your option) 
!   any later version.
!
!   Gemini Unified are distributed in the hope that they
!   will be useful, but WITHOUT ANY WARRANTY; without even the implied
!   warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!   See the GNU General Public License for more details.
!
!   You should have received a copy of the GNU General Public License
!   along with Gemini Unified.
!   If not, see <http://www.gnu.org/licenses/>.
!----------------------------------------------------------------------------
!----------------------------------------------------------------
!> \brief Handle structured stream input and output (group part)

!> \author Wolfgang Friederich

!> \par Description
!>  Support for stream access file reading and writing.
!!  File tree starts with a root group that may contain further groups
!!  and datasets. Each subgroup can contain subgroups and datasets
!!  itself. Idea is to create groups and datasets and to build up a
!!  tree under the root group which contains links to all subgroups
!!  and datasets. Datasets are written to file, file positions are stored
!!  in the group and dataset objects. Finally, the complete tree
!!  with information about the data is also written to file in a way
!!  that by reading the file the complete information and data can
!!  be retrieved.
!<---------------------------------------------------------------
module groupStreamAccess
    use fileStreamAccess
    use datasetStreamAccess
    implicit none
    interface createGroupStreamAccess
        module procedure createPostallocatedGroupStreamAccess
        module procedure createPreallocatedGroupStreamAccess
    end interface createGroupStreamAccess
    interface pushVectorGroupStreamAccess
       module procedure pushRealVectorGroupStreamAccess
       module procedure pushComplexVectorGroupStreamAccess
       module procedure pushIntegerVectorGroupStreamAccess
    end interface
    interface pullVectorGroupStreamAccess
       module procedure pullRealVectorGroupStreamAccess
       module procedure pullComplexVectorGroupStreamAccess
       module procedure pullIntegerVectorGroupStreamAccess
    end interface
    interface dealloc
        module procedure nullifyGroupStreamAccess
    end interface
    interface clearGroupTree; module procedure recursiveDeallocGroupStreamAccess; end interface
    interface operator (.ndataset.); module procedure getNdsGroupStreamAccess; end interface
    interface operator (.nsubgroup.); module procedure getNsubGroupStreamAccess; end interface
    interface operator (.groupname.); module procedure getNameGroupStreamAccess; end interface
    interface operator (.grouptag.); module procedure getTagGroupStreamAccess; end interface
!
    type group_stream_access
        private
        integer (longint) :: filepos      !< position where to find subgroup and data information in file
        integer :: tag                    !< positive integer number indexing the group in some way (root has zero!)
        character (len=132) :: name       !< name specifying group properties
        type (group_stream_access), pointer :: parent => null()  !< pointer to parent
        type (group_stream_access), dimension(:), pointer :: subgroup => null()   !< pointer to subgroups
        integer :: maxsubgroup            !< maximum number of subgroups that shall be added to this group
                                          !! optional argument in createGroupStreamAccess
                                          !! if .not.present(maxsubgroup) in createGroupStreamAccess, it is set to -1 and is not used 
                                          !<
        integer :: lastsubgroup           !< index of subgroup that was added last (internal use)
        type (data_stream_access), dimension(:), pointer :: dataset => null()     !< pointer to datasets
        integer :: maxdset                !< maximum number of datasets that shall be added to this group
                                          !! optional argument in createGroupStreamAccess
                                          !! if .not.present(maxdset) in createGroupStreamAccess, it is set to -1 and is not used
                                          !<
        integer :: lastdset               !< index of dataset that was added last (internal use)
    end type
!
    logical :: verboseStreamAccess = .false.
!
 contains
!---------------------------------------------------------------------
!> \brief Create a group. Contents to be filled in later.
!> \param this group_stream_access object
!> \param name name of the group
!> \param tag an identifier of the group
!> \par
!> Maximal number of subgroups and datasets need not to be known a priori, 
!! arrays of subgroups and datasets reallocated frequently
!>
    subroutine createPostallocatedGroupStreamAccess(this,name,tag)
    type (group_stream_access) :: this
    character (len=*) :: name
    integer :: tag
!
    this%name = name
    this%tag = tag
    this%maxsubgroup = -1
    this%lastsubgroup = -1
    this%maxdset = -1
    this%lastdset = -1
    if (verboseStreamAccess) then
        print *,'create group named ',trim(this%name),' Tag = ',this%tag
    endif
    end subroutine createPostallocatedGroupStreamAccess
!---------------------------------------------------------------------
!> \brief Create a group. Contents to be filled in later. 
!> \param this group_stream_access object
!> \param name name of the group
!> \param tag an identifier of the group
!> \param maxsubgroup number of subgroups that are accounted for in array preallocation
!> \param maxdset number of datasets that are accounted for in array preallocation
!> \par
!> Maximal number of subgroups and datasets are know a priori, 
!! arrays of subgroups and datasets are allocated in advance 
!! to avoid frequent reallocation.
!>
    subroutine createPreallocatedGroupStreamAccess(this,name,tag,maxsubgroup,maxdset)
    type (group_stream_access) :: this
    character (len=*) :: name
    integer :: tag
    integer :: maxsubgroup
    integer :: maxdset
    integer :: ierr
!
    this%name = name
    this%tag = tag
    if(maxsubgroup > 0) then
        this%maxsubgroup = maxsubgroup
        allocate(this%subgroup(1:this%maxsubgroup), stat=ierr)
        if(ierr /= 0) stop "allocate error"
    else
        this%maxsubgroup = 0
    endif
    this%lastsubgroup = 0
    if(maxdset > 0) then
        this%maxdset = maxdset
        allocate(this%dataset(1:this%maxdset), stat=ierr)
        if(ierr /= 0) stop "allocate error"
    else
        this%maxdset = 0
    endif
    this%lastdset = 0
    if (verboseStreamAccess) then
        print *,'create group named ',trim(this%name),' Tag = ',this%tag, & 
                        ' preallocated ',this%maxsubgroup,' subgroups and ',this%maxdset,' datasets'
    endif
    end subroutine createPreallocatedGroupStreamAccess
!---------------------------------------------------------------------
!> \brief Recursively deallocate group including subgroups and datasets
!> \param this group_stream_access object
!
    recursive subroutine recursiveDeallocGroupStreamAccess(this)
    type (group_stream_access) :: this
    integer :: j
    if (associated(this%subgroup)) then
        do j=1,size(this%subgroup)
            call recursiveDeallocGroupStreamAccess(this%subgroup(j))
        enddo
        deallocate(this%subgroup)
    endif
    if (associated(this%dataset)) then
        do j=1,size(this%dataset)
            call fullyDeallocDatasetStreamAccess(this%dataset(j))
        enddo
        deallocate(this%dataset)
    endif
    end subroutine recursiveDeallocGroupStreamAccess
!-----------------------------------------------------------------------
!> \brief Only nullify pointers in group object. Do not deallocate.
!> \param this group_stream_access object
!
    subroutine nullifyGroupStreamAccess(this)
    type (group_stream_access) :: this
    if (associated(this%parent)) nullify(this%parent)
    if (associated(this%subgroup)) nullify(this%subgroup)
    if (associated(this%dataset)) nullify(this%dataset)
    end subroutine nullifyGroupStreamAccess
!----------------------------------------------------------------------
!> \brief Add a subgroup to this group.
!> \param this group_stream_access object
!> \param group subgroup object to be incorporated into the tree
!> \par
!> The subgroup must be completely filled with contents.
!! Only a SHALLOW copy of the group object into the tree is performed.
!! Avoid copying the same group into different trees, because
!! clearing one tree leads to a deallocation of memory still associated
!! with pointers in other trees. Clearing these trees will then fail. 
!! After calling this routine \a group may be deallocated using dealloc. 
!<
    subroutine addSubgroupStreamAccess(this,group)
    type (group_stream_access), target :: this
    type (group_stream_access) :: group
    integer :: n,j
!
    if (associated(this%subgroup)) then
        n = size(this%subgroup)
    else
        n = 0
    endif
    if (this%maxsubgroup == -1) then  ! case of frequent reallocation
        this%subgroup => reallocateGroupStreamAccess(this%subgroup,n+1)
        this%subgroup(n+1) = group
        if (verboseStreamAccess) then
            print *,'add group named ',trim(this%subgroup(n+1)%name),' as ',n+1,' th subgroup of group ', &
                  & trim(this%name),' with tag = ',this%tag
        endif
        this%subgroup(n+1)%parent => this
    !
    !  redirect the parent member of group's subgroups from group to this%subgroup(n+1)
    !
        if (associated(group%subgroup)) then
            do j=1,size(group%subgroup)
                group%subgroup(j)%parent => this%subgroup(n+1)
            enddo
        endif
    else  ! case of preallocation, only reallocate when maxsubgroup is reached
        if (this%lastsubgroup == n) &
            this%subgroup => reallocateGroupStreamAccess(this%subgroup,n + this%maxsubgroup/10 + 1) !"+1" for case maxsubgroup < 10
        this%subgroup(this%lastsubgroup+1) = group
        if (verboseStreamAccess) then
            print *,'add group named ',trim(this%subgroup(this%lastsubgroup+1)%name),' as ',this%lastsubgroup+1, &
                  & ' th subgroup of group ',trim(this%name),' with tag = ',this%tag
        endif
        this%subgroup(this%lastsubgroup+1)%parent => this
    !
    !  redirect the parent member of group's subgroups from group to this%subgroup(lastsubgroup+1)
    !
        if (group%lastsubgroup > 0) then
            do j=1,group%lastsubgroup
                group%subgroup(j)%parent => this%subgroup(this%lastsubgroup+1)
            enddo
        endif
        this%lastsubgroup = this%lastsubgroup + 1
    endif
    end subroutine addSubgroupStreamAccess
!----------------------------------------------------------------------
!> \brief Add a dataset to some existing group
!> \param this group_stream_access object
!> \param dset dataset object
!> \par 
!> The dataset object must contain all information. This implies
!! that the data have already been written to file
!<
    subroutine addDatasetGroupStreamAccess(this,dset)
    type (group_stream_access) :: this
    type (data_stream_access) :: dset
    integer :: n
!
    if (associated(this%dataset)) then
        n = size(this%dataset)
    else
        n = 0
    endif
    if (this%maxdset == -1) then
        this%dataset => reallocateDatasetStreamAccess(this%dataset,n+1)
        call createDatasetStreamAccess(this%dataset(n+1),dset%rank,dset%dims,dset%datatype)
        this%dataset(n+1)%filepos_start = dset%filepos_start
        if (verboseStreamAccess) then
            print *,'add ',n+1,' th dataset with size ',this%dataset(n+1)%dims,' to group ',trim(this%name),' with tag ',this%tag
        endif
    else
        if (this%lastdset == n) &
            this%dataset => reallocateDatasetStreamAccess(this%dataset,n + this%maxdset/10 + 1) !"+1" for case maxdset < 10
        call createDatasetStreamAccess(this%dataset(this%lastdset+1),dset%rank,dset%dims,dset%datatype)
        this%dataset(this%lastdset+1)%filepos_start = dset%filepos_start
        if (verboseStreamAccess) then
            print *,'add ',this%lastdset+1,' th dataset with size ',this%dataset(this%lastdset+1)%dims,' to group ', &
                & trim(this%name),' with tag ',this%tag
        endif
        this%lastdset = this%lastdset + 1
    endif
    end subroutine addDatasetGroupStreamAccess
!----------------------------------------------------------------------
!> \brief Get number of datasets contained in group
!> \param this group_stream_access object
!
    integer function getNdsGroupStreamAccess(this)
    type (group_stream_access), intent(in) :: this
    if (associated(this%dataset)) then
        if (this%maxdset == -1) then
            getNdsGroupStreamAccess = size(this%dataset)
        else
            getNdsGroupStreamAccess = this%lastdset
        endif
    else
        getNdsGroupStreamAccess = 0
    endif
    end function getNdsGroupStreamAccess
!----------------------------------------------------------------------
!> \brief Get number of subgroups contains in group
!> \param this group_stream_access object
!
    integer function getNsubGroupStreamAccess(this)
    type (group_stream_access), intent(in) :: this
    if (associated(this%subgroup)) then
        if (this%maxsubgroup == -1) then
            getNsubGroupStreamAccess = size(this%subgroup)
        else
            getNsubGroupStreamAccess = this%lastsubgroup
        endif
    else
        getNsubGroupStreamAccess = 0
    endif
    end function getNsubGroupStreamAccess
!----------------------------------------------------------------------
!> \brief Get tag of group
!> \param this group_stream_access object
!
    integer function getTagGroupStreamAccess(this)
    type (group_stream_access), intent(in) :: this
    getTagGroupStreamAccess = this%tag
    end function getTagGroupStreamAccess
!----------------------------------------------------------------------
!> \brief Get name of group
!> \param this group_stream_access object
!
    character (len=132) function getNameGroupStreamAccess(this)
    type (group_stream_access), intent(in) :: this
    getNameGroupStreamAccess = this%name
    end function getNameGroupStreamAccess
!----------------------------------------------------------------------
!> \brief Get a pointer to the parent group
!> \param this group_stream_access object
!> \return a pointer to the parent group
!
    function getParentGroupStreamAccess(this) result(parent)
    type (group_stream_access) :: this
    type (group_stream_access), pointer :: parent
    parent => this%parent
    end function getParentGroupStreamAccess
!----------------------------------------------------------------------
!> \brief Get pointer to selected dataset in group
!> \param this group_stream_access object
!> \param k index of selected dataset in group
!> \return a pointer to a data_stream_access object
!
    function getSelectedDatasetGroupStreamAccess(this,k) result(dset)
    type (group_stream_access) :: this
    type (data_stream_access), pointer:: dset
    integer :: k
!
    dset => this%dataset(k)
    end function getSelectedDatasetGroupStreamAccess
!----------------------------------------------------------------------
!> \brief Recursively write the group data to file
!> \param this file_stream_access object
!> \param fda file_stream_access object
!
    recursive subroutine writeGroupStreamAccess(this,fda)
    type (group_stream_access) :: this
    type (file_stream_access) :: fda
    integer :: nsub,nds,j,lu,move
    integer (longint) :: filepos
!
    if (.nsubgroup.this > 0) then    ! first work through my subgroups
        do j=1,.nsubgroup.this
            call writeGroupStreamAccess(this%subgroup(j),fda)
        enddo
    endif
!
!  now deal with myself
!
    filepos = getCurrentFilePositionStreamAccess(fda)
    lu = getFileUnitStreamAccess(fda)
!
    this%filepos = filepos           ! position where infos about my group start
                                     ! this info will later be written to file by my parent group
!
!  if I am the root group I write now the first position of the file
!  specifiying the position where info about the root group's contents has been written
!
    if (.not. associated(this%parent)) write(lu,pos = 1) filepos
!
    nsub = getNsubGroupStreamAccess(this)
    nds = getNdsGroupStreamAccess(this)
!
    write(lu,pos=filepos) this%name,this%tag,nsub,nds   ! write name, tag etc
    if (verboseStreamAccess) then
        print *,'write group: ',trim(this%name),' with tag ',this%tag,' nsub = ',nsub,' nds = ',nds
    endif
    filepos = IncrementCurrentFilePositionStreamAccess(fda,len(this%name)+3*kind(1))
!
    if (nsub > 0) then                                  ! write subgroup record info
        write(lu,pos=filepos) (this%subgroup(j)%filepos,j=1,nsub)
        move = nsub*bit_size(filepos)/8
        filepos = IncrementCurrentFilePositionStreamAccess(fda,move)
    endif
    if (nds > 0) then                                   ! write dataset record info
        write(lu,pos=filepos) (this%dataset(j)%filepos_start,j=1,nds)
        move = nds*bit_size(filepos)/8
        filepos = IncrementCurrentFilePositionStreamAccess(fda,move)
    endif
    end subroutine writeGroupStreamAccess
!------------------------------------------------------------------------------------
!> \brief Recursively read group data from file into memory starting with root group
!> \param this group_stream_access object
!> \param fda file from which group data read
!
    recursive subroutine readGroupStreamAccess(this,fda,readInfoDset)
    type (group_stream_access), target :: this
    type (file_stream_access), target :: fda
    logical, optional :: readInfoDset
    logical :: readInfoDsetFlag
    integer :: lu,nsub,nds,j,move
    integer (longint) :: filepos
!
    if(present(readInfoDset)) then; readInfoDsetFlag=readInfoDset; else; readInfoDsetFlag = .true.;endif
!
    lu = getFileUnitStreamAccess(fda)
    filepos = getCurrentFilePositionStreamAccess(fda)
    this%filepos = filepos
    read(lu,pos=filepos) this%name,this%tag,nsub,nds      ! read first group record
    if (verboseStreamAccess) then
        print *,'read group ',trim(this%name),' with tag ',this%tag,' nsub = ',nsub,' nds = ',nds
    endif
    filepos = IncrementCurrentFilePositionStreamAccess(fda,len(this%name)+3*kind(1))
    this%maxsubgroup = nsub; this%lastsubgroup = nsub
    this%maxdset = nds ; this%lastdset = nds
!
    if (nsub > 0) then                                   ! read subgroup record info
        allocate(this%subgroup(nsub))
        do j=1,nsub
            this%subgroup(j)%parent => this
        enddo
        read(lu,pos=filepos) (this%subgroup(j)%filepos,j=1,nsub)
        move = nsub*bit_size(filepos)/8
        filepos = IncrementCurrentFilePositionStreamAccess(fda,move)
    else
        this%subgroup => null()
    endif
    if (nds > 0) then                                    ! read dataset record info
        allocate(this%dataset(nds))
        read(lu,pos=filepos) (this%dataset(j)%filepos_start,j=1,nds)
        move = nds*bit_size(filepos)/8
        filepos = IncrementCurrentFilePositionStreamAccess(fda,move)
    else
        this%dataset => null()
    endif
    if (nsub > 0) then                                   !  read in subgroups
        do j=1,nsub
            call setCurrentFilePositionStreamAccess(fda,this%subgroup(j)%filepos)
            call readGroupStreamAccess(this%subgroup(j),fda,readInfoDsetFlag)   ! read in subgroups
        enddo
    endif
    if (nds > 0 .and. readInfoDsetFlag) then
        do j = 1,nds
!            call setCurrentFilePositionStreamAccess(fda,this%dataset(j)%filepos_start)
            call readInfoDataStreamAccess(this%dataset(j),fda)
        enddo
    endif
    end subroutine readGroupStreamAccess
!----------------------------------------------------------------------
!> \brief Reallocate a group array
!> \param p pointer to an array of groups
!> \param n number of group objects to be allocated
!> \return newp pointer to the new array of groups
!> \par
!> The contents of the old array \a p is copied into the new one.
!! Afterwards, the old array \a p is deallocated.
!<
!--------------------------------------------------------------------
    function reallocateGroupStreamAccess(p, n) result(newp)
    type (group_stream_access), pointer, dimension(:) :: p, newp
    integer, intent(in) :: n
    integer :: nold, ierr
    allocate(newp(1:n), stat=ierr)
    if(ierr /= 0) stop "allocate error"
    if(.not. associated(p)) return
    nold = min(size(p), n)
    newp(1:nold) = p(1:nold)
    deallocate(p)
    end function reallocateGroupStreamAccess
!----------------------------------------------------------------------
!> \brief Traverse group tree along some given path down to some dataset.
!> \param this Path refers to this group
!> \param depth gives level of the group containing the dataset as seen from \a this
!> \param path an index array of indices specifiying the path through the group tree, 
!!  the last index identifies the data set in the group.
!<
!> \param group output pointer to group the index array points to
!> \param dset a pointer to the dataset specified by the path
!> \par
!!  If \a group does not exist, return a null pointer.
!!  If \a dset does not exist, return a null pointer.
!!  The group pointer is returned to allow fast access to other datasets in the group.
!<
    recursive subroutine traversePathGroupStreamAccess(this,depth,path,group,dset)
    type (group_stream_access), target :: this
    integer :: depth
    integer, dimension(:) :: path
    integer, dimension(:), allocatable :: path_copy
    type (group_stream_access), pointer :: group
    type (data_stream_access), pointer :: dset
    integer :: j,depth2
!
    allocate(path_copy(size(path)))
    path_copy = path
    if (depth > 0 ) then
        if (associated(this%subgroup)) then
            j = path_copy(1)                         ! take index from path
            path_copy = cshift(path_copy,1)          ! put first index at end
            depth2 = depth-1                         ! decrease depth
            call traversePathGroupStreamAccess(this%subgroup(j),depth2,path_copy,group,dset)
        else
            group => null()       ! depth can't be reached
            dset => null()
        endif
    else
        j = path_copy(1)
        group => this
        if (associated(this%dataset) .and. j > 0) then
            if (size(this%dataset) >= j) then
                dset => this%dataset(j)
            else
                dset => null()
            endif
        else
            dset => null()
        endif
    endif
    deallocate(path_copy)
    end subroutine traversePathGroupStreamAccess
!--------------------------------------------------------------------------------
!> \brief Print Info about a group
!
    subroutine printInfoGroupStreamAccess(this)
    type (group_stream_access) :: this
    write(*,'(80(1h-))')
    write(*,'(a,i4,a,i5,a,i5,a,a,a)') '|  Group: ',this%tag,'|  N-Subgroups: ',&
         .nsubgroup.this,'|  N-Datasets: ',.ndataset.this,'|  Name: ',trim(this%name),'|'
    write(*,'(80(1h-))')
    end subroutine printInfoGroupStreamAccess
!----------------------------------------------------------------------
!> \brief Recursively write the file hierarchy to screen
!> \param this file_stream_access object
!> \param fda file_stream_access object
!
    recursive subroutine printGroupStreamAccess(this)
    type (group_stream_access) :: this
    integer :: nds,j
!
!  first deal with myself
!
    call printInfoGroupStreamAccess(this)
    nds = getNdsGroupStreamAccess(this)
    if (nds > 0) then
       write(*,'(a11,2a7,a12)') '| Dataset |',' Rank |',' Type |',' Dimensions |'
       do j = 1,nds
          call printDsetInfoStreamAccess(this%dataset(j),j)
       enddo
    endif
!
    if (.nsubgroup.this > 0) then    ! work through my subgroups
        do j=1,.nsubgroup.this
            call printGroupStreamAccess(this%subgroup(j))
        enddo
    endif
    end subroutine printGroupStreamAccess
!----------------------------------------------------------------------
!  Write real vector data to the file and add dataset to group tree
!  Convenience function to create, write, add and deallocate a real vector data set
!
    subroutine pushRealVectorGroupStreamAccess(group,d,n,fda)
        type (group_stream_access) :: group
        real, dimension(:) :: d
        integer :: n
        type (file_stream_access) :: fda
        type (data_stream_access) :: dset
        call createDatasetStreamAccess(dset,1,(/n/),T_REAL)
        call writeDatasetRealVectorStreamAccess(dset,fda,d)
        call addDatasetGroupStreamAccess(group,dset)
        call dealloc(dset)
    end subroutine pushRealVectorGroupStreamAccess
!----------------------------------------------------------------------
!  Write complex vector data to the file and add dataset to group tree
!  Convenience function to create, write, add and deallocate a complex vector data set
!
    subroutine pushComplexVectorGroupStreamAccess(group,d,n,fda)
        type (group_stream_access) :: group
        complex, dimension(:) :: d
        integer :: n
        type (file_stream_access) :: fda
        type (data_stream_access) :: dset
        call createDatasetStreamAccess(dset,1,(/n/),T_COMPLEX)
        call writeDatasetComplexVectorStreamAccess(dset,fda,d)
        call addDatasetGroupStreamAccess(group,dset)
        call dealloc(dset)
    end subroutine pushComplexVectorGroupStreamAccess
!----------------------------------------------------------------------
!  Write integer vector data to the file and add dataset to group tree
!  Convenience function to create, write, add and deallocate a integer vector data set
!
    subroutine pushIntegerVectorGroupStreamAccess(group,d,n,fda)
        type (group_stream_access) :: group
        integer, dimension(:) :: d
        integer :: n
        type (file_stream_access) :: fda
        type (data_stream_access) :: dset
        call createDatasetStreamAccess(dset,1,(/n/),T_INTEGER)
        call writeDatasetIntegerVectorStreamAccess(dset,fda,d)
        call addDatasetGroupStreamAccess(group,dset)
        call dealloc(dset)
    end subroutine pushIntegerVectorGroupStreamAccess
!----------------------------------------------------------------------
!  read real vector data from file
!  Convenience function to read a real vector data set
!
    subroutine pullRealVectorGroupStreamAccess(this,fda,depth,path,d)
        type (group_stream_access) :: this
        type (file_stream_access) :: fda
        integer :: depth
        integer, dimension(:) :: path
        real, dimension(:), pointer :: d
        type (data_stream_access), pointer :: dset
        type (group_stream_access), pointer :: group
        call traversePathGroupStreamAccess(this,depth,path,group,dset)
        call readDatasetVectorStreamAccess(dset,fda,d)
    end subroutine pullRealVectorGroupStreamAccess
!----------------------------------------------------------------------
!  read integer vector data from file
!  Convenience function to read a integer vector data set
!
    subroutine pullIntegerVectorGroupStreamAccess(this,fda,depth,path,d)
        type (group_stream_access) :: this
        type (file_stream_access) :: fda
        integer :: depth
        integer, dimension(:) :: path
        integer, dimension(:), pointer :: d
        type (data_stream_access), pointer :: dset
        type (group_stream_access), pointer :: group
        call traversePathGroupStreamAccess(this,depth,path,group,dset)
        call readDatasetVectorStreamAccess(dset,fda,d)
    end subroutine pullIntegerVectorGroupStreamAccess
!----------------------------------------------------------------------
!  read complex vector data from file
!  Convenience function to read a complex vector data set
!
    subroutine pullComplexVectorGroupStreamAccess(this,fda,depth,path,d)
        type (group_stream_access) :: this
        type (file_stream_access) :: fda
        integer :: depth
        integer, dimension(:) :: path
        complex, dimension(:), pointer :: d
        type (data_stream_access), pointer :: dset
        type (group_stream_access), pointer :: group
        call traversePathGroupStreamAccess(this,depth,path,group,dset)
        call readDatasetVectorStreamAccess(dset,fda,d)
    end subroutine pullComplexVectorGroupStreamAccess
!
end module groupStreamAccess   

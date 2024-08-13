! =======================================================================
!  Drive computation of spheroidal and toroidal unit jump Green functions
! =======================================================================
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
!----------------------------------------------------------------
!! Calculate Green functions for unit vector jumps.
!! The number of unit jump vectors depends on source type and whether
!! the source is in a fluid or a solid layer. Maximum number of unit jump vectors
!! is 4 for spheroidal and 2 for toroidal motion.
!!
!! Spheroidal Green functions are stored as 3D array with indices
!! ordered according to component (U,R,V,S), jump vector and node.
!! In a fluid medium, only 2 components (U,R) are used.
!! Toroidal Green functions are stored as 3D array with indices
!! ordered according to component (W,T), jump vector and node.
!-----------------------------------------------------------------
module unitJumpGreenFunctions
    use spheroidalMinors
    use spheroidalAdjoint
    use toroidalMinors
    use geminiIntegrationEnvironment
    use splineEarthmodelCoefficients
    use string
    use errorMessage
    implicit none
    interface dealloc; module procedure deallocUnitJumpGreenFunctions; end interface
    integer, parameter :: secrecy_ujgf = 4                         ! print screen output if secrecy <= this value
    type :: unit_jump_green_functions
       double complex, dimension(:,:,:), allocatable :: ujgfsph    ! unit jump Green functions for spheroidal motion
       double complex, dimension(:,:,:), allocatable :: ujgftor    ! unit jump Green functions for toroidal motion
                                                                   ! order is (comp,jump,node)
    end type
!
contains
!-----------------------------------------------------------------------------------
!> \brief Compute unit jump Green functions from one source node to all other nodes.
!!
!!  toromin:     object holding toroidal minors
!!  spheromin:   object holding spheroidal minors
!!  jnods:       index of source node
!!  gem_intenv:  integration environment
!!  eps:         desired accuracy
!!  secrecy:     anti-debugging level
!!  errmsg:      error message
!
  subroutine computeUnitJumpGreenFunctions(this,toromin,spheromin,jnods,gem_intenv,eps,secrecy,errmsg)
     type (unit_jump_green_functions) :: this
     type (toroidal_minors) :: toromin
     type (spheroidal_minors) :: spheromin
     integer :: jnods
     type (gemini_integration_environment) :: gem_intenv
     double precision :: eps
     integer :: secrecy
     type (error_message) :: errmsg
     character(len=36) :: myname = 'computeUnitJumpGreenFunctions'
     type (spheroidal_adjoint) :: adjma,adjic,adjma_plus,adjma_minus
     integer, dimension(:), pointer :: state_exnod,layer_exnod
     double precision, dimension(:), pointer :: rnod
     double precision, dimension(:,:), allocatable :: ujvs,ujvt        ! unit jump vectors (spheroidal, toroidal)
     double complex, dimension(4) :: adjma_roc,adjma_ruc,ystart
     double complex :: d,e
     double precision :: ra,re,rstart,ruc,roc,ric
     integer :: nnod,j,n,njs,njt,i
   !-------------------------------------------------------------------
     call addTrace(errmsg,myname)
   !-------------------------------------------------------------------
   !  minors already computed ?
   !
     if (.not. .done.toromin) then
        call add(errmsg,2,'Toroidal minors not yet computed',myname)
        return
     endif
     if (.not. .done.spheromin) then
        call add(errmsg,2,'Spheroidal minors not yet computed',myname)
        return
     endif
   !-------------------------------------------------------------------
     nnod = .nnod.(gem_intenv%exnod)
     if (jnods < 1 .or. jnods > nnod) then
        call add(errmsg,2,'Invalid index for source node',myname)
        return
     endif
   !-------------------------------------------------------------------
   !  Material state at external nodes
   !
     state_exnod => getPointerStateExternalRadialNodes(gem_intenv%exnod)
     layer_exnod => getPointerLayerExternalRadialNodes(gem_intenv%exnod)
     rnod => getPointerDoubleRadiiExternalRadialNodes(gem_intenv%exnod)
   !
   !  Radii of ocean bottom/surface, CMB, ICB
   !
     ruc = .ruc.(gem_intenv%nem)
     roc = .roc.(gem_intenv%nem)
     ric = .ric.(gem_intenv%nem)
   !-------------------------------------------------------------------
   !  Number of unit jump Green functions and unit jump vectors 
   !
     if (state_exnod(jnods) == 1) then             ! source (either force or moment) in ocean
        njs = 2; njt = 0                           ! (1,0) and (0,1) to get (U,R) after reciprocity for receiver in ocean
        allocate(ujvs(2,2))
        ujvs(:,1) = (/1.d0,0.d0/)               
        ujvs(:,2) = (/0.d0,1.d0/)
     else                                          ! source (either force or moment) in solid
        njs = 4; njt = 2                           ! all unit jump GFs to later get (U,R,V,S) and (W,T) after reciprocity at receiver in solid
        allocate(ujvs(4,4),ujvt(2,2))
        ujvs(:,1) = (/1.d0,0.d0,0.d0,0.d0/)
        ujvs(:,2) = (/0.d0,1.d0,0.d0,0.d0/)
        ujvs(:,3) = (/0.d0,0.d0,1.d0,0.d0/)
        ujvs(:,4) = (/0.d0,0.d0,0.d0,1.d0/)
        ujvt(:,1) = (/1.d0,0.d0/)
        ujvt(:,2) = (/0.d0,1.d0/)
     endif
   !-----------------------------------------------------------------------
   !  allocate space for unit jump Green functions 
   !  initialize to zero
   !  if source in ocean there are no toroidal motions  
   !
     allocate(this%ujgfsph(4,njs,nnod))
     this%ujgfsph = dcmplx(0.d0,0.d0)
     if (njt > 0) then
        allocate(this%ujgftor(2,njt,nnod))
        this%ujgftor = dcmplx(0.d0,0.d0)
     endif
   ! -----------------------------------------------------------------------
   !  Loop over unit jump vectors in toroidal case.
   !  If toromin%det(jnods) is zero, then rstart was above rnod(jnods)
   !  and jnods was not reached during integration;
   !  Green function is zero then.
   !  Also zero if source is in the ocean  
   !
     if (abs(toromin%det(jnods)) > tiny(1.d0)) then
        rstart = startRadiusInitialValues(gem_intenv,eps,'tor',errmsg)
        if (secrecy <= secrecy_ujgf) print *,'UJGF: toroidal UJGF'
        do j = 1,njt
       !----------------------------------------------------------------------
       !  UJGF at nodes below source node, exclude jnods, only solid nodes, GF stays zero elsewhere
       !
           d = (ujvt(2,j)*toromin%d_tor(1,jnods)-ujvt(1,j)*toromin%d_tor(2,jnods))/toromin%det(jnods)
           do n = 1,jnods-1
              if (state_exnod(n) == 0 .and. rnod(n) > rstart) then
                 this%ujgftor(:,j,n) = d *toromin%u_tor(:,n)
              endif
           enddo
       !-----------------------------------------------------------------------------
       !  UJGF at nodes above source, include jnods, only solid nodes, GF stays zero elsewhere
       !
           e = (ujvt(2,j)*toromin%u_tor(1,jnods)-ujvt(1,j)*toromin%u_tor(2,jnods))/toromin%det(jnods)
           do n =jnods,nnod
              if (state_exnod(n) == 0) then
                 this%ujgftor(:,j,n) = e*toromin%d_tor(:,n)
              endif
           enddo
           if (secrecy <= secrecy_ujgf) then
              do n = 1,nnod
                 write(6,'(2i6,8e15.3)') j,n,this%ujgftor(1:2,j,n)
              enddo
           endif
        enddo                  ! end jumps
     endif                     ! end toroidal UJGFs
     if (allocated(ujvt)) deallocate(ujvt)
   !-----------------------------------------------------------------------------
   !  Loop over unit jump vectors in spheroidal case.
   !  If dets is zero, then rstart was above rnod(jnods)
   !  and jnods was not reached during integration;
   !  solution is zero then.
   !  If source is in the ocean:
   !      - first calculate constants e and d
   !      - calculate ujgf from already available basis solutions, e*g(r) and d*h(r)
   !      - the integrate adjoint from ocean bottom (r3) to lowermost node or starting radius or CMB (r2) 
   !        with initial condition P_4m(r3,r3) = (0,0,0,1), only once
   !      - for nodes in the outer core, use available basis solution and multiply with dP_43(r3,r2)
   !      - for nodes in the inner core, integrate adjoint from ICB downwards
   !        with initial condition P_4m(r1,r1) = (0,0,0,1) and then multiply available minors with dP_43(r3,r2)P_4m(r1,r)  
   !  If source is in the mantle:
   !      - integrate adjoint from source to ocean bottom (r3) with initial condition (eq. 48)
   !      - calculate e1 from adjoint solution at r3 and multiply basis solution in ocean by e1
   !      - integrate adjoint from source downwards and calculate b1 from solution at CMB (r2)
   !      - multiply basis solution in OC by b1
   !      - in inner core, integrate adjoint from ICB (r1) downwards with initial condition (0,0,0,1) = P_4m(r1,r1)
   !      - multiply minors by b1*P_4m(r1,r) to get solution
   !-------------------------------------------------------------------------------
   !  Case l=0, only 2 UJGFs needed, other UJGFs stay zero
   !  
     if (abs(spheromin%det(jnods)) > tiny(1.d0)) then
        if (gem_intenv%dll1 < epsilon(1.d0)) then
           if (secrecy <= secrecy_ujgf) print *,'UJGF: spheroidal UJGF, l=0'
           do j = 1,2
              d = (ujvs(2,j)*spheromin%d_minors(1,jnods)-ujvs(1,j)*spheromin%d_minors(2,jnods))/spheromin%det(jnods)
              e = (ujvs(2,j)*spheromin%u_minors(1,jnods)-ujvs(1,j)*spheromin%u_minors(2,jnods))/spheromin%det(jnods)
              do n = jnods,nnod
                 this%ujgfsph(1:2,j,n) = e*spheromin%d_minors(1:2,n)        ! UJGF above source
              enddo
              do n = jnods-1,1,-1
                 this%ujgfsph(1:2,j,n) = d*spheromin%u_minors(1:2,n)        ! UJGF below source
              enddo
              if (secrecy <= secrecy_ujgf) then
                 do n = 1,nnod
                    write(6,'(2i6,8e15.3)') j,n,this%ujgfsph(1:2,j,n)
                 enddo
              endif
           enddo
           if (allocated(ujvs)) deallocate(ujvs)
           return                                                           ! return here, rest not needed
        endif
   !                                                                        ! l > 0 case
        rstart = startRadiusInitialValues(gem_intenv,eps,'sph',errmsg)
     !
     !  propagator P_4m(r) in the inner core, integrate spheroidal adjoint once from ICB to rstart
     !  needed for both cases, only done if rstart in inner core and if there are external nodes in IC
     !
        ra = ric
        re = max(rstart,rnod(1))       ! integrate down to either start radius or bottom node (if in IC)
        if (re .lt. ra) then
           if (secrecy <= secrecy_ujgf) print *,'UJGF: propagator P_4m in inner core'
           ystart = (/0.d0,0.d0,0.d0,1.d0/)
           call computeSpheroidalAdjoint(adjic,ystart,ra,re,gem_intenv,eps,secrecy,errmsg)
           if (.level.errmsg == 2) then
              call add(errmsg,2,'Problems computing spheroidal adjoint',myname)
              call dealloc(adjic)
              return
           endif
        endif
     !
     !  source is in the ocean   
     !
        if (state_exnod(jnods) == 1 .and. layer_exnod(jnods) == .nlay.gem_intenv%nem) then
        !
        !  propagator P_4m(r) in the mantle, integrate spheroidal adjoint once from ocean bottom to CMB
        !
           if (secrecy <= secrecy_ujgf) print *,'UJGF: source is in the ocean'
           ra = ruc
           re = max(rstart,roc,rnod(1))                   ! lower limit either rstart or roc or bottom node
           adjma_roc = 0.d0                               ! zero to avoid ompiler warning
           if (re .lt. ra) then                           ! start radius in the mantle
              if (secrecy <= secrecy_ujgf) print *,'UJGF: propagator P_4m in the mantle'
              ystart = (/0.d0,0.d0,0.d0,1.d0/)
              call computeSpheroidalAdjoint(adjma,ystart,ra,re,gem_intenv,eps,secrecy,errmsg)
              if (.level.errmsg == 2) then
                 call add(errmsg,2,'Problems computing spheroidal adjoint',myname)
                 call dealloc(adjma)
                 return
              endif
              adjma_roc = ystart                          ! value at CMB (if re = roc) else not used
           endif
        !
           if (secrecy <= secrecy_ujgf) print *,'UJGF: spheroidal UJGF'
           do j = 1,njs
              d = (ujvs(2,j)*spheromin%d_minors(1,jnods)-ujvs(1,j)*spheromin%d_minors(2,jnods))/spheromin%det(jnods)
              e = (ujvs(2,j)*spheromin%u_minors(1,jnods)-ujvs(1,j)*spheromin%u_minors(2,jnods))/spheromin%det(jnods)
           !
           !  UJGF in the ocean
           !  above source
           !   
              do n = jnods,nnod
                 this%ujgfsph(1:2,j,n) = e*spheromin%d_minors(1:2,n)
              enddo
           !
           !  below source
           !
              do n = jnods-1,1,-1
                 if (state_exnod(n) == 1 .and. rnod(n) .ge. rstart) then    ! ocean and outer core
                    if (rnod(n) > ruc) then
                       this%ujgfsph(1:2,j,n) = d*spheromin%u_minors(1:2,n)
                    else if(rnod(n) < roc) then
                       this%ujgfsph(1:2,j,n) = d*adjma_roc(3)*spheromin%u_minors(1:2,n)
                    endif
                 else if(state_exnod(n) == 0 .and. rnod(n) .ge. rstart) then    ! mantle and inner core
                    if (rnod(n) > roc) then
                       this%ujgfsph(:,j,n) = d*matmul(spheromin%u_mtil(1:4,:,n),adjma%adj(:,n))
                    else if(rnod(n) < ric) then
                       this%ujgfsph(:,j,n) = d*adjma_roc(3)*matmul(spheromin%u_mtil(1:4,:,n),adjic%adj(:,n))
                    endif
                 endif
              enddo
              if (secrecy <= secrecy_ujgf) then
                 do n = 1,nnod
                    write(6,'(2i6,8e15.3)') j,n,this%ujgfsph(1:4,j,n)
                 enddo
              endif
           enddo
           call dealloc(adjma)
     !
     !  source is in the mantle
     !
        else
           if (secrecy <= secrecy_ujgf) print *,'UJGF: source is in the mantle'
           do j = 1,njs
           !
           !  adjoint solution in the mantle, above the source
           !  fm(+)(r_s) = -s_i/delta*meps(-)_mi;   "+" = above source, "-" = below source
           !  "u_meps" means upward integration, i.e. below source   
           !
              if (secrecy <= secrecy_ujgf) print *,'UJGF: adjoint solution aboce the source, jump =  ',j
              ystart = 0.d0
              do i = 1,njs
                 ystart = ystart-ujvs(i,j)*spheromin%u_meps(:,i,jnods)/spheromin%det(jnods)
              enddo
              re = ruc
              call computeSpheroidalAdjoint(adjma_plus,ystart,rnod(jnods),re,gem_intenv,eps,secrecy,errmsg)
              if (.level.errmsg == 2) then
                 call add(errmsg,2,'Problems computing spheroidal adjoint',myname)
                 call dealloc(adjma_plus)
                 return
              endif
              adjma_ruc = ystart
           !
           !  adjoint solution in the mantle, below the source
           !  fm(-)(r_s) = +s_i/delta*meps(+)_mi;   "+" = above source, "-" = below source    
           !  "d_meps" means downwward integration, i.e. above source   
           !
              if (secrecy <= secrecy_ujgf) print *,'UJGF: adjoint solution below the source, jump =  ',j
              ystart = 0.d0
              do i = 1,njs
                 ystart = ystart+ujvs(i,j)*spheromin%d_meps(:,i,jnods)/spheromin%det(jnods)
              enddo
              re = max(rstart,roc,rnod(1))             ! integrate either to rstart, roc or bottom node 
              call computeSpheroidalAdjoint(adjma_minus,ystart,rnod(jnods),re,gem_intenv,eps,secrecy,errmsg)
              if (.level.errmsg == 2) then
                 call add(errmsg,2,'Problems computing spheroidal adjoint',myname)
                 call dealloc(adjma_minus)
                 return
              endif
              adjma_roc = ystart              ! adjoint solution at CMB (if rstart, rnod(1) <= roc else not used)
           !
           !  UJGF at nodes above source (in mantle or ocean)
           !
              do n = jnods,nnod
                 if (state_exnod(n) == 0 .and. rnod(n) < ruc) then
                    this%ujgfsph(1:4,j,n) = matmul(spheromin%d_mtil(1:4,:,n),adjma_plus%adj(:,n))
                 else if(state_exnod(n) == 1 .and. rnod(n) > ruc) then
                    this%ujgfsph(1:2,j,n) = adjma_ruc(3)*spheromin%d_minors(1:2,n)
                 endif
              enddo
           !
           !  UJGF at mantle nodes below source (either mantle, outer core or inner core)
           !  if there are none, goto next jump index
           !
              do n = jnods-1,1,-1
                 if(state_exnod(n) == 0 .and. rnod(n) .ge. rstart) then      ! mantle or inner core
                    if (rnod(n) > roc) then
                       this%ujgfsph(1:4,j,n) = matmul(spheromin%u_mtil(1:4,:,n),adjma_minus%adj(:,n))
                    else if(rnod(n) < ric) then
                       this%ujgfsph(1:4,j,n) = adjma_roc(3)*matmul(spheromin%u_mtil(1:4,:,n),adjic%adj(:,n))
                    endif
                 else if(state_exnod(n) == 1 .and. rnod(n) .ge. rstart) then        ! outer core
                    this%ujgfsph(1:2,j,n) = adjma_roc(3)*spheromin%u_minors(1:2,n)
                 endif
              enddo
              call dealloc(adjma_plus)
              call dealloc(adjma_minus)
              if (secrecy <= secrecy_ujgf) then
                 do n = 1,nnod
                    write(6,'(2i6,8e15.3)') j,n,this%ujgfsph(1:4,j,n)
                 enddo
              endif
           enddo                          ! end loop over basis jumps
        endif                             ! end source in the mantle
        if (secrecy <= secrecy_ujgf) print *,'UJGF: finished'
     else
        print *,'Spheroidal determinant less than tiny(1.d0) = ',tiny(1.d0)
        print *,'l = ',sqrt(gem_intenv%dll1),' om = ',gem_intenv%omre, ' jnods = ',jnods 
     endif                                ! end spheroidal UJGFs
     if (allocated(ujvs)) deallocate(ujvs)
  end subroutine computeUnitJumpGreenFunctions
!-------------------------------------------------------------------------------------
! Deallocate unit jump Green functions object
!
  subroutine deallocUnitJumpGreenFunctions(this)
     type (unit_jump_green_functions) :: this
     if (allocated(this%ujgfsph)) deallocate(this%ujgfsph)
     if (allocated(this%ujgftor)) deallocate(this%ujgftor)
  end subroutine deallocUnitJumpGreenFunctions
!
end module unitJumpGreenFunctions

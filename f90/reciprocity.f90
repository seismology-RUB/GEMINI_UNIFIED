! ====================================================================
! Reciprocity relations for unit jump Green functions
! ====================================================================
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
!--------------------------------------------------------------------------
!  Implements reciprocity relations as described in rezifrechet.tex.
!--------------------------------------------------------------------------
module reciprocity
    use errorMessage
    use string
    implicit none
!
contains
!-----------------------------------------------------------------------------------------
!> \brief Transform spheroidal unit jump Green matrix using reciprocity principle
!! Uses the 16 explicit relations derived in the document rezifrechet.tex.
!! Convention for indices of BOTH Green matrices is:
!! -- first index:  component of solution vector (URVS)
!! -- second index: index of basis jump vector (where the delta is)
!!
!!    G(1,l)(re,rs) = +(rs/re)**2*G(i,2)(rs,re)*T(i,l)
!!    G(2,l)(re,rs) = -(rs/re)**2*G(i,1)(rs,re)*T(i,l)
!!    G(3,l)(re,rs) = +(rs/re)**2*G(i,4)(rs,re)*T(i,l)/elp1
!!    G(4,l)(re,rs) = -(rs/re)**2*G(i,3)(rs,re)*T(i,l)/elp1
!! If you want all components for one selected jump, you need to calculate
!! all 4 unit jump Green functions.
!!-----------------------------------------------------------------------------------
!!    G(1,1) = +(rs/re)**2*G(2,2)*T(2,1)      = -(rs/re)**2*G(2,2)
!!    G(2,1) = -(rs/re)**2*G(2,1)*T(2,1)      = +(rs/re)**2*G(2,1)
!!    G(3,1) = +(rs/re)**2*G(2,4)*T(2,1)/elp1 = -(rs/re)**2*G(2,4)/elp1
!!    G(4,1) = -(rs/re)**2*G(2,3)*T(2,1)/elp1 = +(rs/re)**2*G(2,3)/elp1
!!
!!    G(1,2) = +(rs/re)**2*G(1,2)*T(1,2)      = +(rs/re)**2*G(1,2)
!!    G(2,2) = -(rs/re)**2*G(1,1)*T(1,2)      = -(rs/re)**2*G(1,1)
!!    G(3,2) = +(rs/re)**2*G(1,4)*T(1,2)/elp1 = +(rs/re)**2*G(1,4)/elp1
!!    G(4,2) = -(rs/re)**2*G(1,3)*T(1,2)/elp1 = -(rs/re)**2*G(1,3)/elp1
!!    
!!    G(1,3) = +(rs/re)**2*G(4,2)*T(4,3)      = -(rs/re)**2*G(4,2)*elp1
!!    G(2,3) = -(rs/re)**2*G(4,1)*T(4,3)      = +(rs/re)**2*G(4,1)*elp1
!!    G(3,3) = +(rs/re)**2*G(4,4)*T(4,3)/elp1 = -(rs/re)**2*G(4,4)
!!    G(4,3) = -(rs/re)**2*G(4,3)*T(4,3)/elp1 = +(rs/re)**2*G(4,3)
!!
!!    G(1,4) = +(rs/re)**2*G(3,2)*T(3,4)      = +(rs/re)**2*G(3,2)*elp1
!!    G(2,4) = -(rs/re)**2*G(3,1)*T(3,4)      = -(rs/re)**2*G(3,1)*elp1
!!    G(3,4) = +(rs/re)**2*G(3,4)*T(3,4)/elp1 = +(rs/re)**2*G(3,4)
!!    G(4,4) = -(rs/re)**2*G(3,3)*T(3,4)/elp1 = -(rs/re)**2*G(3,3)
!!----------------------------------------------------------------------------------
!! Dimension of input unit jump Green functions (bgin) can be:
!!   2x2: nodes in fluid region and source in ocean
!!   4x2: nodes in solid region and source in ocean
!!   2x4: nodes in fluid region and source in mantle
!!   4x4: nodes in solid region and source in mantle
!! Output dimensions (bgout) are switched:
!!   2x2: source node in fluid region and receiver in ocean
!!   2x4: source node in solid region and receiver in ocean
!!   4x2: source node in fluid region and receiver in mantle
!!   4x4: source node in solid region and receiver in mantle
!!---------------------------------------------------------------------------------
!!  elp1:    l*(l+1) (in)
!!  re:      radius of receiver position (in)
!!  r:       node radius (in)
!!  state:   either solid (0) or fluid (1)
!!  bgin:    unit jump Green functions at nodes from source located at receiver position, G_{ik}(r,re) (in)
!!  bgout:   unit jump Green functions at receiver position for sources located at nodes, G_{nl}(re,r) (out)
!
    subroutine spheroidalReciprocity(elp1,re,r,state_node,state_src,bgin,bgout,errmsg)
    double precision :: elp1,re,r,rlp1,f
    integer :: state_node,state_src
    double complex, dimension(:,:) :: bgin                    ! (comp,jump)
    double complex, dimension(:,:), allocatable :: bgout      ! (comp,jump)
    type (error_message) :: errmsg
    character(len=21) :: myname = 'spheroidalReciprocity'
  !
    call addTrace(errmsg,myname)
    rlp1 = 1.d0/elp1
    f = (r/re)**2
    if (state_node == 1) then         ! node in fluid (2 components)
       if (state_src == 1) then       ! source in ocean (2 jumps)
          allocate(bgout(2,2))
          bgout(1,1) = -f*bgin(2,2)              !  G(1,1) = -(rs/re)**2*G(2,2)
          bgout(2,1) = +f*bgin(2,1)              !  G(2,1) = +(rs/re)**2*G(2,1)
          bgout(1,2) = +f*bgin(1,2)              !  G(1,2) = +(rs/re)**2*G(1,2)
          bgout(2,2) = -f*bgin(1,1)              !  G(2,2) = -(rs/re)**2*G(1,1)
       else if (state_src == 0) then  ! source in mantle (4 jumps)
          allocate(bgout(4,2))
          bgout = 0.d0
          bgout(1,1) = -f*bgin(2,2)              !  G(1,1) = -(rs/re)**2*G(2,2)
          bgout(2,1) = +f*bgin(2,1)              !  G(2,1) = +(rs/re)**2*G(2,1)
          if (elp1 > epsilon(1.d0)) then    ! l > 0
             bgout(3,1) = -f*bgin(2,4)*rlp1         !  G(3,1) = -(rs/re)**2*G(2,4)/elp1
             bgout(4,1) = +f*bgin(2,3)*rlp1         !  G(4,1) = +(rs/re)**2*G(2,3)/elp1
          endif
       !
          bgout(1,2) = +f*bgin(1,2)              !  G(1,2) = +(rs/re)**2*G(1,2)
          bgout(2,2) = -f*bgin(1,1)              !  G(2,2) = -(rs/re)**2*G(1,1)
          if (elp1 > epsilon(1.d0)) then    ! l > 0
             bgout(3,2) = +f*bgin(1,4)*rlp1         !  G(3,2) = +(rs/re)**2*G(1,4)/elp1
             bgout(4,2) = -f*bgin(1,3)*rlp1         !  G(4,2) = -(rs/re)**2*G(1,3)/elp1
          endif
       endif
    else if (state_node == 0) then     ! node in solid (4 components)
       if (state_src == 1) then        ! source in ocean (2 jumps)
          allocate(bgout(2,4))
          bgout(1,1) = -f*bgin(2,2)              !  G(1,1) = -(rs/re)**2*G(2,2)
          bgout(1,2) = +f*bgin(1,2)              !  G(1,2) = +(rs/re)**2*G(1,2)
          bgout(1,3) = -f*bgin(4,2)*elp1         !  G(1,3) = -(rs/re)**2*G(4,2)*elp1
          bgout(1,4) = +f*bgin(3,2)*elp1         !  G(1,4) = +(rs/re)**2*G(3,2)*elp1
       !
          bgout(2,1) = +f*bgin(2,1)              !  G(2,1) = +(rs/re)**2*G(2,1)
          bgout(2,2) = -f*bgin(1,1)              !  G(2,2) = -(rs/re)**2*G(1,1)
          bgout(2,3) = +f*bgin(4,1)*elp1         !  G(2,3) = +(rs/re)**2*G(4,1)*elp1
          bgout(2,4) = -f*bgin(3,1)*elp1         !  G(2,4) = -(rs/re)**2*G(3,1)*elp1
       else if (state_src == 0) then   ! source in mantle (4 jumps)
          allocate(bgout(4,4))
          bgout(1,1) = -f*bgin(2,2)              !  G(1,1) = -(rs/re)**2*G(2,2)
          bgout(2,1) = +f*bgin(2,1)              !  G(2,1) = +(rs/re)**2*G(2,1)
          bgout(3,1) = -f*bgin(2,4)*rlp1         !  G(3,1) = -(rs/re)**2*G(2,4)/elp1
          bgout(4,1) = +f*bgin(2,3)*rlp1         !  G(4,1) = +(rs/re)**2*G(2,3)/elp1
  !
          bgout(1,2) = +f*bgin(1,2)              !  G(1,2) = +(rs/re)**2*G(1,2)
          bgout(2,2) = -f*bgin(1,1)              !  G(2,2) = -(rs/re)**2*G(1,1)
          bgout(3,2) = +f*bgin(1,4)*rlp1         !  G(3,2) = +(rs/re)**2*G(1,4)/elp1
          bgout(4,2) = -f*bgin(1,3)*rlp1         !  G(4,2) = -(rs/re)**2*G(1,3)/elp1
  !
          bgout(1,3) = -f*bgin(4,2)*elp1         !  G(1,3) = -(rs/re)**2*G(4,2)*elp1
          bgout(2,3) = +f*bgin(4,1)*elp1         !  G(2,3) = +(rs/re)**2*G(4,1)*elp1
          bgout(3,3) = -f*bgin(4,4)              !  G(3,3) = -(rs/re)**2*G(4,4)
          bgout(4,3) = +f*bgin(4,3)              !  G(4,3) = +(rs/re)**2*G(4,3)
  !
          bgout(1,4) = +f*bgin(3,2)*elp1         !  G(1,4) = +(rs/re)**2*G(3,2)*elp1
          bgout(2,4) = -f*bgin(3,1)*elp1         !  G(2,4) = -(rs/re)**2*G(3,1)*elp1
          bgout(3,4) = +f*bgin(3,4)              !  G(3,4) = +(rs/re)**2*G(3,4)
          bgout(4,4) = -f*bgin(3,3)              !  G(4,4) = -(rs/re)**2*G(3,3)
       endif
    endif
  !
    end subroutine spheroidalReciprocity
!---------------------------------------------------------------------
!> \brief Transform toroidal basis Green matrix using reciprocity principle
!!
!! Convention for indices of BOTH Green matrices is:
!! -- first index:  component of solution vector (WT)
!! -- second index: index of basis jump vector
!!
!!  re:      radius of receiver position (in)
!!  r:       node radius (in)
!!  sourcetype: force or moment source
!!  bgin:    basis Green functions at nodes from source at receiver position (in)
!!  bgout:   basis Green functions at receiver position for sources at nodes (out)
!!
!!  only called if source or receiver are in solid medium
!!--------------------------------------------------------------------
    subroutine toroidalReciprocity(re,r,bgin,bgout)
    double precision :: re,r,f
    double complex, dimension(:,:) :: bgin                    ! (comp,jump)
    double complex, dimension(:,:), allocatable :: bgout      ! (comp,jump)
   !
    f = (r/re)**2
    allocate(bgout(2,2))
    bgout(1,1) = -f*bgin(2,2)              !  G(1,1) = -(rs/re)**2*G(2,2)
    bgout(2,1) = +f*bgin(2,1)              !  G(2,1) = +(rs/re)**2*G(2,1)
    bgout(1,2) = +f*bgin(1,2)              !  G(1,2) = +(rs/re)**2*G(1,2)
    bgout(2,2) = -f*bgin(1,1)              !  G(2,2) = -(rs/re)**2*G(1,1)
    end subroutine toroidalReciprocity
!
end module reciprocity

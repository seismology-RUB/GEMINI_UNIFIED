!--------------------------------------------------
!  test passing of array to subroutines
!  see also https://pages.mtu.edu/~shene/COURSES/cs201/NOTES/chap08/assumed.html
!
program testArrayPassing
    use testArrayPassingSubs
    implicit none
    real, dimension(0:2,-1:1) :: x
    real, dimension(:,:), allocatable :: y
    x(0,:) = [0.,3.,6.]
    x(1,:) = [1.,4.,7.]
    x(2,:) = [2.,5.,8.]
    allocate(y(0:2,-1:1))
    y(0,:) = [0.,3.,6.]
    y(1,:) = [1.,4.,7.]
    y(2,:) = [2.,5.,8.]
    call sub1(x)               ! pass explicitly dimenioned array
    call sub1(y)               ! pass allocated array
    call sub2(x(:,0))          ! pass second column of array
    call sub2(x(0,:))          ! pass first row of array
    call sub1(y(:,0:))         ! do not pass first column
end program testArrayPassing

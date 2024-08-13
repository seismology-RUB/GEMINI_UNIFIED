program testHeapsort
    use heapSort
    use mathConstants
    implicit none
    real, dimension(:,:), allocatable :: a
    integer, dimension(:), allocatable :: idx
    integer :: n1,n2,i,j
!
    n1 = 10
    n2 = 5
    allocate(a(n1,n2))
    allocate(idx(n1))
    forall (i = 1:n1) idx(i) = i
    do j = 1,n2
       do i = 1,n1
          a(i,j) = sin(mc_two_pid*dble(i-0.5)/dble(n1))*cos(mc_pid*dble(j-0.5)/dble(n2))
       enddo
       write(6,'(10f9.5)') a(:,j)
    enddo
    write(6,'(10i9)') idx
    
    call heapSort2D(a,idx,2,4)

    do j = 1,n2
       write(6,'(10f9.5)') a(:,j)
    enddo
    write(6,'(10i9)') idx
end program testHeapsort

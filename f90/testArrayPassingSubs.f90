module testArrayPassingSubs
    implicit none
contains
    
subroutine sub1(x)
    real, dimension(:,:) :: x
    integer :: i,j,dims(2)
    dims = shape(x)
    print *,'Shape of X is: ',dims
    do i = 1,dims(1)
       print *,(x(i,j),j=1,dims(2))
    enddo
end subroutine sub1

subroutine sub2(x)
    real, dimension(:) :: x
    integer :: i,dim(1)
    dim = shape(x)
    print *,'Shape of X is: ',dim
    print *,(x(i),i=1,dim(1))
end subroutine sub2
    
end module testArrayPassingSubs

program OOP_vector_operations_example
    implicit none
    
    integer(kind=4), parameter :: ngrid = 12, nmat = 2
    real(kind=8), dimension(:), allocatable :: avec, bvec
    integer(kind=4), dimension(:), allocatable :: Ivec, Jvec
    real(kind=8), dimension(:,:), allocatable :: amat, bmat
    integer(kind=4) :: igrid, jgrid 

    allocate(avec(ngrid), bvec(ngrid))
    allocate(Ivec(ngrid), Jvec(ngrid))

    do igrid = 1, ngrid
     avec(igrid) = igrid
     bvec(igrid) = 0.50d0*igrid

     Ivec(igrid) = igrid
     Jvec(igrid) = 2*igrid
    enddo

    call check_dot_product(ngrid, Ivec, Jvec)
    call vector_operation_check(ngrid, avec, bvec)
    deallocate(avec, bvec)
    
    allocate(amat(nmat, nmat), bmat(nmat, nmat))
    
    do igrid = 1, nmat
        amat(igrid, igrid) = 1.00d0
        bmat(igrid, igrid) = 2.30d0
        do jgrid = 1, igrid - 1 
            amat(igrid, jgrid) = 0.50d0
            amat(jgrid, igrid) = 0.50d0

            bmat(igrid, jgrid) = 1.50d0
            bmat(jgrid, igrid) = 1.50d0
        enddo
    enddo

    call matrix_operation_check(nmat, nmat, amat, nmat, nmat, bmat)
    deallocate(amat, bmat)
end program OOP_vector_operations_example

subroutine check_dot_product(n, avec, bvec)
    use Vector_operation_BLAS
    implicit none
    integer(kind=4), intent(in) :: n
    integer(kind=4), intent(in) :: avec(n)
    integer(kind=4), intent(in) :: bvec(n)
    integer(kind=4) :: igrd
    integer(kind=4) :: val
    val = dotProduct(n, avec, bvec, 3)
    do igrd = 1, n
        write(*,'(*(i5))') avec(igrd),bvec(igrd)
    enddo
    write(*,*)"Dot product :: ", val
end subroutine check_dot_product

subroutine vector_operation_check(ngrid, avec, bvec)
    use Vector_operation_BLAS
    implicit none

    integer(kind=4), intent(in) :: ngrid
    real(kind=8), dimension(ngrid), intent(in) :: avec, bvec

    type(Vector) :: vec1, vec2, result_vec
    real(kind=8) :: dot_result, mag_result

    ! Initialize vectors
    call vec1%initialize(ngrid, avec)
    call vec2%initialize(ngrid, bvec)
    
    ! Add two vectors
    call vec1%add(vec2, result_vec)
    print *, "Vector Addition Result:", result_vec%real_elements
    deallocate(result_vec%real_elements)
    
    ! Subtract one vector from another
    call vec1%subtract(vec2, result_vec)
    print *, "Vector Subtraction Result:", result_vec%real_elements
    
    ! Compute the dot product of two vectors
    call vec1%dotProductPre(vec2, dot_result)
    print *, "Dot Product Result:", dot_result

    ! Compute the magnitude of a vector
    mag_result = vec1%magnitude()
    print *, "Magnitude of Vector 1:", mag_result

end subroutine vector_operation_check

subroutine matrix_operation_check(nrowa, ncola, amat, nrowb, ncolb, bmat)
    use MatrixOperations
    implicit none

    integer(kind=4), intent(in) :: nrowa, ncola, nrowb, ncolb
    real(kind=8), intent(in) :: amat(nrowa, ncola)
    real(kind=8), intent(in) :: bmat(nrowb, ncolb)
    
    type(Matrix) :: rmat1, rmat2, result_rmat
    
    ! Initialize matrices
    call rmat1%init(size(amat,1), size(amat,2), amat)
    call rmat2%init(size(bmat,1), size(bmat,2), bmat)
    
    ! Perform matrix multiplication
    call rmat1%Matmul(rmat2, result_rmat)
    print *, result_rmat%real_elements
    print *, "Matrix Multiplication Result:" 
    call result_rmat%print_Matrix()
    
    ! Compute the inverse of a matrix
    call rmat1%inverse(result_rmat)
    print *, "Matrix Inverse Result:"
    print *, result_rmat%real_elements

    call rmat1%Reset_memory()
end subroutine matrix_operation_check
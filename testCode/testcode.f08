
module test_code
    implicit none
    contains

    ! Created by prash on 6/5/2024.
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

    subroutine check_rmat()
        use MatrixModule
        type(Matrix) :: matrix1, matrix2, result

        integer(kind=4) :: i

        ! Initialize the matrices
        call initialize_matrix(matrix1, 3, 3)
        call initialize_matrix(matrix2, 3, 3)

        ! Populate the matrices with some data (for demonstration purposes)
        matrix1%data = reshape([1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0], [3, 3])
        matrix2%data = reshape([9.0, 8.0, 7.0, 6.0, 5.0, 4.0, 3.0, 2.0, 1.0], [3, 3])

        print*, "Matrix 1 :"
        do i = 1, matrix1%rows
            print *, matrix1%data(i, :)
        end do

        print*, "Matrix 2 :"
        do i = 1, matrix2%rows
            print *, matrix2%data(i, :)
        end do

        ! Perform element-wise multiplication
        result = matrix1 .Multe. matrix2
        ! Display the result (for demonstration purposes)
        print *, "Result of element-wise multiplication:"
        do i = 1, result%rows
            print *, result%data(i, :)
        end do

        ! Perform Matrix multiplication
        result = matrix1 .Mult. matrix2

        ! Display the result (for demonstration purposes)
        print *, "Result of element-wise multiplication:"
        do i = 1, result%rows
            print *, result%data(i, :)
        end do
    end subroutine check_rmat

end module test_code

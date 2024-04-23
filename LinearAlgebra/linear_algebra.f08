! A module for vector and Matrix operations
! Functionality :
! DataType :: Integer, Real, Complex
!  - Dot product of two Vector
!  - Scale the vector with a constant
!  - Add vector to an another vector (scaled with constant)
!  - Matrix vector multiplication
module Vector_operation_BLAS
    implicit none
    private
    ! Define a type for a vector
    type, public :: Vector
        integer(kind=4), allocatable :: int_elements(:)
        real(kind=8), allocatable :: real_elements(:)
        complex(kind=8), allocatable :: complex_elements(:)
        integer(kind=4) :: size
    contains
        procedure :: initialize
        procedure :: add
        procedure :: subtract
        procedure :: dotProductPre
        procedure :: magnitude
    end type Vector

    interface
        ! Initialize the vector
        module subroutine initialize(this, size, values)
            implicit none
            class(Vector) :: this
            integer(kind=4), intent(in) :: size
            real(kind=8), dimension(size), intent(in) :: values
        end subroutine initialize
        ! Add two vectors
        module subroutine add(this, other, result)
            implicit none
            class(Vector), intent(in) :: this, other
            class(Vector), intent(out) :: result
        end subroutine add
        ! Subtract one vector from another
        module subroutine subtract(this, other, result)
            implicit none
            class(Vector), intent(in) :: this, other
            class(Vector), intent(out) :: result
        end subroutine subtract
        ! Dot product of two vectors
        module subroutine dotProductPre(this, other, result)
            implicit none
            class(Vector), intent(in) :: this, other
            real(kind=8), intent(out) :: result
        end subroutine dotProductPre
        ! Magnitude of a vector
        module real(kind=8) function magnitude(this)
            implicit none
            class(Vector) :: this
        end function magnitude
    end interface

    public :: dotProduct
    interface dotProduct
        procedure RealdotProduct
        procedure ComplexdotProduct
        module integer(kind=4) function IntegerdotProduct(n, vec1, vec2, increment_index) result(res)
            implicit none
            integer(kind=4), intent(in) :: n, increment_index
            integer(kind=4), intent(in) :: vec1(n)
            integer(kind=4), intent(in) :: vec2(n)
        end function IntegerdotProduct
        module integer(kind=8) function IntegerdotProduct2(n, vec1, vec2, increment_index) result(res)
            implicit none
            integer(kind=4), intent(in) :: n, increment_index
            integer(kind=8), intent(in) :: vec1(n)
            integer(kind=8), intent(in) :: vec2(n)
        end function IntegerdotProduct2
    end interface dotProduct

    public :: ScaleVector
    interface ScaleVector
        procedure Real_ScaleVector
        procedure Real_Scale_Complex_Vector
        procedure complex_Scale_Complex_Vector
        module subroutine Integer_ScaleVector(n, Vec, alpha)
            implicit none
                integer(kind=4), intent(in) :: n
                integer(kind=8), intent(in) :: alpha
                integer(kind=8), intent(inout) :: Vec(n)
        end subroutine Integer_ScaleVector
    end interface ScaleVector

    public :: Vector_Add_with_ConstVal
    interface Vector_Add_with_ConstVal
        procedure Real_VectorAdd
        procedure Real_VectorAdd2
        procedure Integer_VectorAdd
        procedure Integer_VectorAdd2
        procedure Complex_VectorAdd
        procedure Complex_VectorAdd2
    end interface Vector_Add_with_ConstVal

    public :: MatVecMul
    interface MatVecMul
        procedure Real_Matrix_Vector_Mult
        procedure Integer_Matrix_Vector_Mult
        procedure Complex_Matrix_Vector_Mult
    end interface MatVecMul
    
contains

! Dot product of two different vector
! Real(kind=8), integer(kind=4)
! integer(kind=8), complex(kind=8)
    real(kind=8) function RealdotProduct(n, vec1, vec2, increment_index) result(res)
    implicit none 
        integer(kind=4), intent(in) :: n, increment_index
        real(kind=8), intent(in) :: vec1(n)
        real(kind=8), intent(in) :: vec2(n)
        real(kind=8) :: ddot
            res = ddot(n, vec1, increment_index, vec2, increment_index)
        return
    end function RealdotProduct

    complex(kind=8) function ComplexdotProduct(n, vec1, vec2, increment_index) result(res)
    implicit none    
        integer(kind=4), intent(in) :: n, increment_index
        complex(kind=8), intent(in) :: vec1(n)
        complex(kind=8), intent(in) :: vec2(n)
        complex(kind=8) :: zdotc
            res = zdotc(n, vec1, increment_index, vec2, increment_index)
        return
    end function ComplexdotProduct

    ! Scale the vector with a constant
    ! Real double precision vector Scaled with real constant
    subroutine Real_ScaleVector(n, Vec, alpha)
    implicit none
        integer(kind=4), intent(in) :: n
        real(kind=8), intent(in) :: alpha
        real(kind=8), intent(inout) :: Vec(n)    
        call dscal(n, alpha, Vec, 1)
    end subroutine Real_ScaleVector
    ! complex double precision vector Scaled with real constant
    subroutine Real_Scale_Complex_Vector(n, Vec, alpha)
    implicit none
        integer(kind=4), intent(in) :: n
        real(kind=8), intent(in) :: alpha
        integer(kind=8), intent(inout) :: Vec(n)
        call zdscal(n, alpha, Vec, 1)    
    end subroutine Real_Scale_Complex_Vector
    ! complex double precision vector Scaled with complex constant
    subroutine complex_Scale_Complex_Vector(n, Vec, alpha)
    implicit none
        integer(kind=4), intent(in) :: n
        complex(kind=8), intent(in) :: alpha
        complex(kind=8), intent(inout) :: Vec(n)
        call zscal(n, alpha, Vec, 1)    
    end subroutine complex_Scale_Complex_Vector

! Addition of two vectors with a constant value
! Real(kind=8)
! integer(kind=8)
! complex(kind=8)
    subroutine Integer_VectorAdd(n, vec1, vec2, alpha)
    implicit none
        integer(kind=4), intent(in) :: n
        integer(kind=8), intent(in) :: alpha
        integer(kind=8), intent(in) :: vec1(n)
        integer(kind=8), intent(inout) :: vec2(n)   
    end subroutine Integer_VectorAdd

    subroutine Integer_VectorAdd2(n, vec1, vec2, vec3, alpha)
    implicit none
        integer(kind=4), intent(in) :: n
        integer(kind=8), intent(in) :: alpha
        integer(kind=8), intent(in) :: vec1(n), vec2(n)
        integer(kind=8), intent(out) :: vec3(n)
    end subroutine Integer_VectorAdd2

    subroutine Real_VectorAdd(n, vec1, vec2, alpha)
    implicit none
        integer(kind=4), intent(in) :: n
        real(kind=8), intent(in) :: alpha
        real(kind=8), intent(in) :: vec1(n)
        real(kind=8), intent(inout) :: vec2(n)
        
    end subroutine Real_VectorAdd

    subroutine Real_VectorAdd2(n, vec1, vec2, vec3, alpha)
    implicit none
        integer(kind=4), intent(in) :: n
        real(kind=8), intent(in) :: alpha
        real(kind=8), intent(in) :: vec1(n), vec2(n)
        real(kind=8), intent(out) :: vec3(n)
    end subroutine Real_VectorAdd2

    subroutine Complex_VectorAdd(n, vec1, vec2, alpha)
    implicit none
        integer(kind=4), intent(in) :: n
        complex(kind=8), intent(in) :: alpha
        complex(kind=8), intent(in) :: vec1(n)
        complex(kind=8), intent(inout) :: vec2(n)   
    end subroutine Complex_VectorAdd

    subroutine Complex_VectorAdd2(n, vec1, vec2, vec3, alpha)
    implicit none
        integer(kind=4), intent(in) :: n
        complex(kind=8), intent(in) :: alpha
        complex(kind=8), intent(in) :: vec1(n), vec2(n)
        complex(kind=8), intent(out) :: vec3(n)
    end subroutine Complex_VectorAdd2

    subroutine Integer_Matrix_Vector_Mult(n, m, mat, vec, res)
    implicit none
        integer(kind=4), intent(in) :: n, m
        integer(kind=8), intent(in) :: mat(n, m), vec(m)
        integer(kind=8), intent(out) :: res(n)
    end subroutine Integer_Matrix_Vector_Mult

    subroutine Real_Matrix_Vector_Mult(n, m, Mat, vec, res)
    implicit none
        integer(kind=4), intent(in) :: n, m
        real(kind=8), intent(in) :: mat(n, m), vec(m)
        real(kind=8), intent(out) :: res(n)
    End subroutine Real_Matrix_Vector_Mult

    subroutine Complex_Matrix_Vector_Mult(n, m, mat, vec, res)
    implicit none
        integer(kind=4), intent(in) :: n, m
        complex(kind=8), intent(in) :: mat(n, m), vec(m)
        complex(kind=8), intent(out) :: res(n)
    end subroutine Complex_Matrix_Vector_Mult

end module Vector_operation_BLAS

!- Module for the Matrix Algebra
!- Matrix Inverse
!- Matrix Multiplication (real, complex)
!- Matrix Diagonalization (real, complex)
!- Matrix singular value decomposition
!- Generalize eigen-value problem for symmetric matrix
module MatrixOperations
    implicit none

    ! Define a type for a matrix
    type :: Matrix
        real(kind=8), allocatable :: real_elements(:,:)
        complex(kind=8), allocatable :: complex_elements(:,:)
        integer(kind=4) :: rows, cols
    contains
        procedure :: real_initialize
        procedure :: complex_initialize
        procedure :: real_multiply
        procedure :: complex_multiply
        procedure :: inverse
        generic :: print_vector => print_real_vector, print_complex_vector                
        generic :: print_Matrix => print_real_matrix, print_complex_matrix
        procedure :: print_real_matrix
        procedure :: print_complex_matrix
        procedure :: print_real_vector
        procedure :: print_complex_vector
    end type Matrix

    contains

    ! Method to initialize a matrix with given dimensions and values
    subroutine real_initialize(this, rows, cols, values)
        implicit none
        class(Matrix) :: this
        integer(kind=4), intent(in) :: rows, cols
        real(kind=8), dimension(rows, cols), intent(in) :: values
        this%rows = rows
        this%cols = cols
        this%real_elements = values
    end subroutine real_initialize

    subroutine complex_initialize(this, rows, cols, values)
        implicit none
        class(Matrix) :: this
        integer(kind=4), intent(in) :: rows, cols
        complex(kind=8), dimension(rows, cols), intent(in) :: values
        this%rows = rows
        this%cols = cols
        this%complex_elements = values
    end subroutine complex_initialize

    subroutine print_real_vector(this)
        class(Matrix) :: this
        integer(kind=4) :: igrid
        write(*,*) "-----------------------------------------"
        write(*,*) " Real-Matrix "
        write(*,*) "-----------------------------------------"
        do igrid = 1, this%rows
         write(*,'(*(f10.5))') this%real_elements(igrid, 1)
        enddo
    end subroutine print_real_vector

    subroutine print_complex_vector(this, dummy)
        integer(kind=4), intent(in) :: dummy
        class(Matrix) :: this
        integer(kind=4) :: igrid
        write(*,*) "-----------------------------------------"
        write(*,*) " Real-Matrix "
        write(*,*) "-----------------------------------------"
        do igrid = 1, this%rows
         write(*,'(*(f10.5))') real(this%complex_elements(igrid, 1)), aimag(this%complex_elements(igrid, 1))
        enddo
    end subroutine print_complex_vector

    subroutine print_real_matrix(this)
        class(Matrix) :: this
        integer(kind=4) :: igrid
        write(*,*) "-----------------------------------------"
        write(*,*) " Real-Matrix "
        write(*,*) "-----------------------------------------"
        do igrid = 1, this%rows
         write(*,'(*(E10.5))') this%real_elements(igrid,:)
        enddo
    end subroutine print_real_matrix

    subroutine print_complex_matrix(this, dummy)
        integer(kind=4), intent(inout) :: dummy
        class(Matrix) :: this
        integer(kind=4) :: igrid
        dummy = 0
        write(*,*) "-----------------------------------------"
        write(*,*) " Complex-Matrix Real elements "
        write(*,*) "-----------------------------------------"
        do igrid = 1, this%rows
         write(*,'(*(E12.5))') real(this%complex_elements(igrid,:))
        enddo

      write(*,*) "-----------------------------------------"
      write(*,*) " Complex-Matrix complex elements "
      write(*,*) "-----------------------------------------"
      do igrid = 1, this%rows
        write(*,'(*(E12.5))') aimag(this%complex_elements(igrid,:))
      enddo
    end subroutine print_complex_matrix

    ! Method to perform matrix multiplication of real using LAPACK's DGEMM subroutine
    subroutine real_multiply(this, other, result)
        implicit none
        class(Matrix) :: this, other, result
        real(kind=8), dimension(:,:), allocatable :: temp
        integer(kind=4) :: m, n, k
        real(kind=8) :: alpha, beta
        
        if (this%cols /= other%rows) then
            print *, "Error: Number of columns of first matrix must equal number of rows of second matrix for multiplication"
            return
        endif
    
        m = this%rows
        n = other%cols
        k = this%cols
        alpha = 1.00d0
        beta = 0.00d0
        allocate(temp(m, n))
        call dgemm('N', 'N', m, n, k, alpha, this%real_elements, m, other%real_elements, k, beta, temp, m)
        call result%real_initialize(m, n, temp)
        deallocate(temp)
    end subroutine real_multiply

    ! Method to perform matrix multiplication complex using LAPACK's ZGEMM subroutine
    subroutine complex_multiply(this, other, result)
        implicit none
        class(Matrix) :: this, other, result
        complex(kind=8), dimension(:,:), allocatable :: temp
        integer(kind=4) :: m, n, k
        real(kind=8) :: alpha, beta
        
        if (this%cols /= other%rows) then
            print *, "Error: Number of columns of first matrix must equal number of rows of second matrix for multiplication"
            return
        endif
        
        m = this%rows
        n = other%cols
        k = this%cols
        alpha = 1.00d0
        beta = 0.00d0
        allocate(temp(m, n))
        call zgemm('N', 'N', m, n, k, alpha, this%complex_elements, m, other%complex_elements, k, beta, temp, m)
        call result%complex_initialize(m, n, temp)
        deallocate(temp)
    end subroutine complex_multiply

    ! Method to compute the inverse of a matrix using LAPACK's DGETRF and DGETRI subroutines
    subroutine inverse(this, result)
        implicit none
        class(Matrix) :: this, result
        real(kind=8), dimension(:,:), allocatable :: copy
        real(kind=8), dimension(:), allocatable :: work
        integer(kind=4), dimension(:), allocatable :: ipiv
        integer(kind=4) :: info, n, lda, lwork
        
        if (this%rows /= this%cols) then
            print *, "Error: Matrix must be square for inversion"
            return
        endif
        
        n = this%rows
        lda = n
        lwork = n*n
    allocate(copy(n,n), ipiv(n), work(lwork))
        
        copy = this%real_elements
        call dgetrf(n, n, copy, lda, ipiv, info)
        if (info /= 0) then
            print *, "Error: Matrix is singular"
            return
        endif
        
        call dgetri(n, copy, lda, ipiv, work, lwork, info)
        if (info /= 0) then
            print *, "Error: Failed to compute inverse"
            return
        endif
        
        call result%real_initialize(n, n, copy)
        
    deallocate(copy, ipiv, work)
    end subroutine inverse

    ! Method to compute selected eigen values and vectors of a real Hermitian matrix using LAPACK's dsyevx
    subroutine real_Hermitian_selected_eigenvectors(this, iu, eval, evec)
        implicit none
        class(Matrix) :: this
        character (len=1), parameter ::  jobz = "V"
        character (len=1), parameter :: range = "I"
        character (len=1), parameter ::  uplo = "L"
        
        real (kind=8), allocatable :: w(:), work(:), z(:,:)
        integer (kind=4), allocatable :: iwork(:), ifail(:)
        integer (kind=4) :: lda, il, iu
        integer (kind=4) :: ldz, lwork, info
        integer (kind=4) :: m, ix, iy
        real (kind=8) :: eval(iu), evec(this%rows, iu)
        real (kind=8) :: vl, vu, abstol, dlamch
        
        abstol = 2*dlamch('s')
        ldz = this%rows
        lda = this%rows
        lwork = 8*this%rows
        il = 1
        M = IU - IL + 1
    allocate(z(ldz, M))
    allocate(w(this%rows))
    allocate(work(lwork))
    allocate(ifail(this%rows))
    allocate(iwork(5*this%rows))
        
        call dsyevx ( JOBZ, RANGE, UPLO, this%rows, this%real_elements, lda, vl, vu, IL, IU, &
                     ABSTOL, M, w, Z, LDZ, WORK, LWORK, IWORK, IFAIL, INFO )
        
        do ix = 1, m
         eval(ix) = w(ix)
        enddo
        
        do ix = 1, this%rows
         do iy = 1, m
            evec(ix, iy) = z(ix, iy)
         enddo
        enddo    
    deallocate(w, work, z, iwork, ifail)
    end subroutine real_Hermitian_selected_eigenvectors

    ! Method to compute the eigen values and vectors of a real Hermitian matrix using LAPACK's DSYEV subroutines
    SUBROUTINE real_Hermitian_diagonalization(this, eval)
        implicit none
        class(Matrix) :: this
        real(kind=8), intent(inout) :: eval(this%rows)
        CHARACTER (LEN = 1), PARAMETER  :: JOBZ="V"
        CHARACTER (LEN = 1), PARAMETER  :: UPLO="L"
        REAL    (KIND = 8), ALLOCATABLE :: WORK(:)
        INTEGER (KIND = 4) :: NDIM,INFO
        INTEGER (KIND = 4) :: LDA, LDWORK
        
        LDA = NDIM
        LDWORK = 3*NDIM - 1
        
    ALLOCATE(WORK(LDWORK))
        CALL DSYEV (JOBZ, UPLO, NDIM, this%real_elements, LDA, eval, WORK, LDWORK, INFO)
        if (info /= 0) then
            print *, "Error: Failed to compute Eigenvalues and Eigenvectors of Real symmetric matrix"
            return
        endif
    DEALLOCATE(WORK)
    END SUBROUTINE real_Hermitian_diagonalization

    ! Method to compute the eigen values and vectors of a complex Hermitian matrix using LAPACK's zheev subroutines
    subroutine complex_Hermitian_diagonalization(this, eval)
        implicit none
        class(Matrix) :: this
        real(kind=8), intent(inout) :: eval(this%rows)
        character*1, parameter :: jobz = 'V', uplo = 'U'        
        real(kind=8), allocatable     :: w(:), rwork(:)
        complex(kind=8), allocatable :: work(:)
        integer(kind=4)   :: num, lda, lwork, info
        
        lda = num
        lwork = max(1,2*num-1)

        allocate(w(num), work(MAX(1,LWORK)), rwork(max(1, 3*num-2)))
        call zheev ( jobz, uplo, num, this%complex_elements, lda, w, work, lwork, rwork, info )
        
        if (info /= 0) then
            print *, "Error: Failed to compute Eigenvalues and Eigenvectors of complex Hermitian matrix"
            return
        endif
    eval = w
    deallocate(w, work, rwork)    
    end subroutine complex_Hermitian_diagonalization

    ! Method to compute the eigen values and vectors of a general complex matrix using LAPACK's zgeev subroutines
    SUBROUTINE general_complex_diagonalization(this, eigvals, vecsL, vecsR)
        implicit none
        character (len = 1),parameter  :: JOBVL="V"
        character (len = 1),parameter  :: JOBVR="V"
        class(Matrix) :: this
        COMPLEX (kind = 8) :: eigvals(this%rows), vecsl(this%rows,this%cols), vecsr(this%rows,this%cols)        
        integer (kind = 4) :: LDA, LDVL, LDVR, LDWORK
        integer (kind = 4) :: INFO
        COMPLEX (kind = 8), allocatable:: WORK(:),WORK2(:)
        
        LDA = this%rows
        LDVL = this%rows
        LDVR = this%rows
        LDWORK = 2*this%rows

    allocate(WORK(LDWORK),WORK2(LDWORK))
        CALL ZGEEV(JOBVL, JOBVR, this%rows, this%complex_elements, LDA, eigvals, vecsl, &
                   LDVL, vecsr, LDVR, WORK, LDWORK, WORK2, INFO)
              
        if (info /= 0) then
            print *, "Error: Failed to compute Eigenvalues and Eigenvectors of Complex Non-Hermitian Matrix"
            return
        endif
    deallocate(WORK,WORK2)
    END SUBROUTINE general_complex_diagonalization

    ! Method to solve linear equation of a real-symmetric matrix using LAPACK's dsysv subroutines
    subroutine linear_equation_real_symmetric(this, bmat, nvariables)
        implicit none
        class(Matrix) :: this
        integer(kind=4), intent(in) :: nvariables    
        character(len=1), parameter :: uplo = "U"
        real(kind=8), intent(inout) :: bmat(this%rows, nvariables)        
        integer(kind=4), allocatable :: ipiv(:)
        real(kind=8), allocatable :: work(:)
        integer(kind=4) :: lda, ldb, lwork, info
        
        lda = max(1,this%rows)
        ldb = max(1,this%rows)
        lwork = max(1,this%rows*nvariables)
        
    allocate(ipiv(this%rows), work(lwork))
        call dsysv ( uplo, this%rows, nvariables, this%real_elements, lda, ipiv, bmat, ldb, work, lwork, info )

        if (info /= 0) then
            print *, "Error: Failed to solve linear equation for real symmetic Matrix"
            return
        endif

    deallocate(ipiv, work)
    end subroutine linear_equation_real_symmetric

    ! Method to solve linear equation of a real general matrix using LAPACK's DGESV subroutines
    subroutine linear_equation_real_matrix(this, bmat, ncol)
        implicit none
        class(Matrix) :: this
        integer(kind=4), intent(in) :: ncol
        real(kind=8), intent(inout) :: bmat(this%rows, ncol)
        integer(kind=4), allocatable :: ipiv(:)
        integer(kind=4) :: lda, ldb, lwork, info
        
        lda = max(1, this%rows)
        ldb = max(1, this%rows)
        lwork = max(1, this%rows*ncol)
        
    allocate(ipiv(this%rows))
        ipiv = 0
        call DGESV (this%rows, ncol, this%real_elements, lda, ipiv, bmat, ldb, info)   
        if (info /= 0) then
            print *, "Error: Failed to solve linear equation for real General Matrix"
            return
        endif 
    deallocate(ipiv)
    end subroutine linear_equation_real_matrix

    ! Method for singular value decomposition of complex general matrix using LAPACK's zgesvd subroutines
    subroutine real_singualar_value(this, s, u, vt)
        implicit none
        class(Matrix) :: this
        character(len=1), parameter :: jobu='A', jobvt='A'
        real(kind=8), allocatable :: work(:)
      
        integer(kind=4) :: lda, ldu, ldvt, lwork, info
        real(kind=8) :: s(min(this%rows, this%cols))
        real(kind=8) :: u(this%rows, this%rows), vt(this%cols, this%cols)
      
          lda = max(1, this%rows)
          ldu = this%rows
         ldvt = this%cols
        lwork = MAX(1, 3*MIN(this%rows, this%cols) + MAX(this%rows, this%cols), &
                5*MIN(this%rows, this%cols))
      
       allocate(work(lwork))
      
        call dgesvd(jobu, jobvt, this%rows, this%cols, this%real_elements, lda, &
                    s, u, ldu, vt, ldvt, work, lwork, info)

        if (info /= 0) then
            print *, "Error: Failed to solve singular value decomposition of real matrix"
            return
        endif 
    deallocate(work)
    end subroutine real_singualar_value

    ! Method for singular value decomposition of complex general matrix using LAPACK's zgesvd subroutines
    subroutine complex_singualar_value(this, s, u, vt)
    implicit none    
        class(Matrix) :: this
        character(len=1), parameter  :: jobu='A', jobvt='A'
        complex(kind=8), allocatable :: work(:)
        real(kind=8), allocatable     :: rwork(:)
        real(kind=8)    :: s(min(this%rows, this%cols))
        complex(kind=8) :: u(this%rows, this%rows), vt(this%cols, this%cols)
        integer(kind=4) :: lda, ldu, ldvt, lwork, info
        
        lda = max(1,this%rows)
        ldu = this%rows
        ldvt = this%cols
        lwork = MAX(1,2*MIN(this%rows, this%cols) + MAX(this%rows, this%cols))
        
    allocate(work(lwork), rwork(5*min(this%rows, this%cols)))
        
        call zgesvd(jobu, jobvt, this%rows, this%cols, this%complex_elements, lda, s, u, &
                    ldu, vt, ldvt, work, lwork, rwork, info)

        if (info /= 0) then
            print *, "Error: Failed to solve singular value decomposition of real matrix"
            return
        endif

    deallocate(work, rwork)
    end subroutine complex_singualar_value

    ! Method for evaluating complex generalized symmetric-definite eigenproblem using LAPACK's dsygv subroutines
    subroutine real_geneign_sym_(test, bmat, aevl)
        implicit none
        class(Matrix) :: test
        integer(kind=4), parameter :: itype = 1
        character(len=1), parameter :: jobz = "V"
        character(len=1), parameter :: uplo = "U"
         real(kind=8), intent(in)    :: bmat(test%rows, test%cols)
         real(kind=8), intent(out)   :: aevl(test%rows)
         real(kind=8), allocatable :: work(:)
         integer(kind=4) :: info, ldb, lwork, lda
       
          lda = test%rows
          ldb = test%cols
          lwork = max(1,3*test%rows - 1)
       
        allocate(work(lwork))
       
            work = 0.00d0; aevl = 0.0d0
          call dsygv (itype, jobz, uplo, test%rows, test%real_elements, lda, bmat, &
                        ldb, aevl, work, lwork, info)
        
        if (info /= 0) then
            print *, "Error: Failed to solve generalized eigenvalue problem for real symmetric matrix"
            return
        endif
         
    deallocate(work)
    end subroutine real_geneign_sym_
end module MatrixOperations

! A fortran module for Fourier transform of the Real and complex vector
! 1, 2, 3-Dimensional vector
module FastFourierTransform
    implicit none
    private

    type, abstract :: showable

    end type showable

end module FastFourierTransform
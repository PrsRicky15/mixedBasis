submodule (Vector_operation_BLAS) Integer_vector
    implicit none
contains
    ! The Vector is initialized by the Type
    ! Method to initialize a vector with given size and values
    module procedure initialize
        implicit none
        allocate(this%real_elements(size))
        this%real_elements = values
        this%size = size
    end procedure initialize

    ! Method to add two vectors
    module procedure add
        implicit none
        integer(kind=4) :: i
        if (this%size /= other%size) then
            print *, "Error: Vectors must have the same size for addition"
            return
        endif
        allocate(result%real_elements(this%size))
        do i = 1, this%size
            result%real_elements(i) = this%real_elements(i) + other%real_elements(i)
        end do
        result%size = this%size
    end procedure add
    ! Method to subtract one vector from another
    module procedure subtract
        implicit none
        if (this%size /= other%size) then
            print *, "Error: Vectors must have the same size for subtraction"
            return
        endif
        result%real_elements = this%real_elements - other%real_elements
        result%size = this%size
    end procedure subtract
    ! Method to compute the dot product of two vectors
    module procedure dotProductPre
        implicit none
        integer(kind=4) :: i
        if (this%size /= other%size) then
            print *, "Error: Vectors must have the same size for dot product"
            return
        endif
        result = 0.0d0
        do i = 1, this%size
            result = result + this%real_elements(i) * other%real_elements(i)
        end do
    end procedure dotProductPre
    ! Method to compute the magnitude of a vector
    module procedure magnitude
        implicit none
        integer(kind=8) :: i
        magnitude = 0.0d0
        do i = 1, this%size
            magnitude = magnitude + this%real_elements(i)**2
        end do
        magnitude = sqrt(magnitude)
    end procedure magnitude

    module procedure IntegerdotProduct
        implicit none
        INTEGER(kind=4) :: incx, incy, i, ix, iy, m, mp1
        integer(kind=4) :: Itemp
        INTRINSIC mod

        incx = increment_index
        incy = increment_index
        Itemp = 0
        res = 0
        IF (n <= 0) RETURN
        IF (incx == 1 .AND. incy == 1) THEN
            m = mod(n,5)
            IF (m /= 0) THEN
                DO i = 1,m
                    Itemp = Itemp + vec1(i) * vec2(i)
                END DO
                IF (n < 5) THEN
                    res = Itemp
                    RETURN
                END IF
            END IF

            mp1 = m + 1
            DO i = mp1, n, 5
                Itemp = Itemp + vec1(i)*vec2(i) + vec1(i+1)*vec2(i+1) + &
                        vec1(i+2)*vec2(i+2) + vec1(i+3)*vec2(i+3) + vec1(i+4)*vec2(i+4)
            END DO
        else
            ix = 1
            iy = 1
            IF (incx < 0) ix = (-n+1)*incx + 1
            IF (incy < 0) iy = (-n+1)*incy + 1
            DO i = 1, n
                Itemp = Itemp + vec1(ix)*vec2(iy)
                ix = ix + incx
                iy = iy + incy
            END DO
        END IF
        res = Itemp
        return
    end procedure IntegerdotProduct

    module procedure IntegerdotProduct2
        implicit none
        INTEGER(kind=4) :: incx, incy, i, ix, iy, m, mp1
        integer(kind=8) :: stemp
        INTRINSIC mod

        incx = increment_index; incy = increment_index
        stemp = 0
        res = 0
        IF (n <= 0) RETURN
        IF (incx == 1 .AND. incy == 1) THEN
            m = mod(n,5)
            IF (m /= 0) THEN
                DO i = 1,m
                    stemp = stemp + vec1(i)*vec2(i)
                END DO
                IF (n < 5) THEN
                    res = stemp
                    RETURN
                END IF
            END IF

            mp1 = m + 1
            DO i = mp1, n, 5
                stemp = stemp + vec1(i)*vec2(i) + vec1(i+1)*vec2(i+1) + &
                        vec1(i+2)*vec2(i+2) + vec1(i+3)*vec2(i+3) + vec1(i+4)*vec2(i+4)
            END DO
        else
            ix = 1
            iy = 1
            IF (incx < 0) ix = (-n+1)*incx + 1
            IF (incy < 0) iy = (-n+1)*incy + 1
            DO i = 1,n
                stemp = stemp + vec1(ix)*vec2(iy)
                ix = ix + incx
                iy = iy + incy
            END DO
        END IF
        res = stemp
        return
    end procedure IntegerdotProduct2

    module procedure Integer_ScaleVector
        implicit none
        INTEGER(kind=4) ::  INCX
        INTEGER(kind=4) :: I, M, MP1, NINCX
        INTRINSIC MOD
        IF (N <= 0 .OR. INCX <= 0 .OR. alpha == 1) RETURN
        IF (INCX==1) THEN
            M = MOD(N,5)
            IF (M /= 0) THEN
                DO I = 1,M
                    Vec(I) = alpha * Vec(I)
                END DO
                IF (N < 5) RETURN
            END IF
            MP1 = M + 1
            DO I = MP1,N,5
                Vec(I) = alpha * Vec(I)
                Vec(I+1) = alpha * Vec(I+1)
                Vec(I+2) = alpha * Vec(I+2)
                Vec(I+3) = alpha * Vec(I+3)
                Vec(I+4) = alpha * Vec(I+4)
            END DO
        ELSE
            NINCX = N*INCX
            DO I = 1,NINCX,INCX
                Vec(I) = alpha * Vec(I)
            END DO
        END IF
        RETURN
    end procedure Integer_ScaleVector

    module procedure Real_Matrix_Vector_Mult
        implicit none
    end procedure Real_Matrix_Vector_Mult

    module procedure complex_Matrix_Vector_Mult
        implicit none
    end procedure complex_Matrix_Vector_Mult

end submodule Integer_vector
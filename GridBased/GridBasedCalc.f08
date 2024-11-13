! Created by prash on 6/6/2024.

module GridBasedCalc
    implicit none
    private
    ! Define a type for a Matrix grid
    type, public :: GridMat
        real(kind=8), public,  allocatable :: real_elements(:,:)
        complex(kind=8), public, allocatable :: complex_elements(:,:)
        integer(kind=4), public :: ndims
        integer(kind=4), public :: ngrid
    contains
        procedure, public :: init => initialize
        procedure, public :: reset => Reset_memory
        final :: Deallocate_memory
    end type GridMat

end module GridBasedCalc
module vector_module
  implicit none
  
  type :: vector
     integer, dimension(:), allocatable :: int_elements
     real, dimension(:), allocatable :: real_elements
     integer :: size
  end type vector
  
  interface operator(+)
     module procedure add_vectors_integer
     module procedure add_vectors_real
  end interface
  
  interface operator(*)
     module procedure multiply_vectors_integer
     module procedure multiply_vectors_real
  end interface
  
contains

  function add_vectors_integer(a, b) result(result_vector)
    type(vector), intent(in) :: a, b
    type(vector) :: result_vector
    
    integer :: i
    
    result_vector%size = a%size
    allocate(result_vector%int_elements(result_vector%size))
    
    do i = 1, a%size
       result_vector%int_elements(i) = a%int_elements(i) + b%int_elements(i)
    end do
    
  end function add_vectors_integer

    function add_vectors_real(a, b) result(result_vector)
    type(vector), intent(in) :: a, b
    type(vector) :: result_vector

    integer :: i

    result_vector%size = a%size
    allocate(result_vector%real_elements(result_vector%size))

    do i = 1, a%size
       result_vector%real_elements(i) = a%real_elements(i) + b%real_elements(i)
    end do

  end function add_vectors_real

  function multiply_vectors_integer(a, b) result(result_vector)
    type(vector), intent(in) :: a, b
    type(vector) :: result_vector
    
    integer :: i
    
    result_vector%size = a%size
    allocate(result_vector%int_elements(result_vector%size))
    
    do i = 1, a%size
       result_vector%int_elements(i) = a%int_elements(i) * b%int_elements(i)
    end do
    
  end function multiply_vectors_integer

    function multiply_vectors_real(a, b) result(result_vector)
    type(vector), intent(in) :: a, b
    type(vector) :: result_vector

    integer :: i

    result_vector%size = a%size
    allocate(result_vector%real_elements(result_vector%size))

    do i = 1, a%size
       result_vector%real_elements(i) = a%real_elements(i) * b%real_elements(i)
    end do

  end function multiply_vectors_real

end module vector_module

program main_program
  use vector_module
  implicit none
  
  type(vector) :: vec1, vec2, result_add, result_multiply
  integer :: i
  
  ! Initialize vectors
  vec1%size = 3
  allocate(vec1%int_elements(3))
  vec1%elements = [1, 2, 3]
  
  vec2%size = 3
  allocate(vec2%int_elements(3))
  vec2%elements = [4, 5, 6]
  
  ! Addition
  result_add = vec1 + vec2
  print *, "Addition result:"
  do i = 1, result_add%size
     print *, result_add%int_elements(i)
  end do
  
  ! Multiplication
  result_multiply = vec1 * vec2
  print *, "Multiplication result:"
  do i = 1, result_multiply%size
     print *, result_multiply%int_elements(i)
  end do
  
end program main_program


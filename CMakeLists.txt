cmake_minimum_required(VERSION 3.22)
project(mixedBasis Fortran)

enable_language(Fortran)

message(" Main source directory path :: ${CMAKE_CURRENT_SOURCE_DIR}")
message(" Build directory path :: ${CMAKE_CURRENT_BINARY_DIR}")

find_package(LAPACK REQUIRED)
if(LAPACK_FOUND)
	message("Lapack LIBRARY FOUND ${LAPACK_LIBRARIES}")
endif()

find_package(BLAS REQUIRED)
if(BLAS_FOUND)
    message("Blas LIBRARY FOUND ${BLAS_LIBRARIES}")
endif()

find_package(FFTW3 REQUIRED)
if(FFTW3_FOUND)
    message("FFT LIBRARY FOUND ${FFTW3_LIBRARIES}")
endif()

if(FFTW3_FOUND)
    set(math_library ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} ${FFTW3_LIBRARIES})
endif()

set_source_files_properties(main.f08 PROPERTIES LANGUAGE Fortran)
set_source_files_properties(LinearAlgebra/linear_algebra.f08 PROPERTIES LANGUAGE Fortran)
set_source_files_properties(TDcalc/image_time_prop.f08 PROPERTIES LANGUAGE Fortran)

add_executable(mixedBasis main.f08
                LinearAlgebra/linear_algebra.f08
                TDcalc/image_time_prop.f08)
target_link_libraries(mixedBasis PUBLIC ${math_library})

program onedim_frgrd_imagtim_prop
    use numconst
    use image_time_prop
    implicit none

    real(kind=8), allocatable :: eval(:), evec(:,:)

    real(kind=8)  :: xmin, xmax, mass, cutof,  error, dtim
    integer(kind=4) :: ix, jx, nx, iste, nste

    character (len=80) :: flerr, fleval, flevec

    open(11,file="imgtime.inp")
    read(11,*) xmin, xmax, nx
    read(11,*) mass
    read(11,*) cutof
    read(11,*) dtim, nste
    read(11,*) flerr
    close(11)

    allocate(evec(nx, nste), eval(nste))

    evec = zero
    call imag_time_1d_( xmin, xmax, mass, dtim, cutof, nx, nste, &
            evec, eval, flerr )

    open(11,file="psi.txt")
    do ix = 1, nx
        write(11,*) (evec(ix,iste),iste=1,nste)
    enddo
    close(11)

end

subroutine OOP_vector_operations_example
    use test_code
    implicit none
    
    integer(kind=4), parameter :: ngrid = 12, nmat = 2
    real(kind=8), dimension(:), allocatable :: avec, bvec
    integer(kind=4), dimension(:), allocatable :: Ivec, Jvec
    real(kind=8), dimension(:,:), allocatable :: amat, bmat
    integer(kind=4) :: igrid, jgrid

!    call check_rmat()
!    stop
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
end subroutine OOP_vector_operations_example
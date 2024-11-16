
! Constants for fftw3
module fftpar
    implicit none

    integer(kind=4), parameter :: fftw_r2hc = 0
    integer(kind=4), parameter :: fftw_hc2r = 1
    integer(kind=4), parameter :: fftw_dht = 2
    integer(kind=4), parameter :: fftw_forward = -1
    integer(kind=4), parameter :: fftw_backward = +1
    integer(kind=4), parameter :: fftw_measure = 0
    integer(kind=4), parameter :: fftw_patient = 32
    integer(kind=4), parameter :: fftw_estimate = 64

end

module numconst
    implicit none
    real(kind=8), parameter :: zero = 0.000000d0, one = 1.00000d0
    real(kind=8), parameter :: half = 0.500000d0, two = 2.00000d0
    real(kind=8), parameter :: PI = 3.1415926535897932D+00, twpi = two*pi
    complex(kind=8), parameter :: io = (0.000000d0, 1.000000d0)
    complex(kind=8), parameter :: aval = (one, zero)
    complex(kind=8), parameter :: bval = (zero, zero)
end

function pot(x)
    use numconst
    implicit none

    real(kind=8) :: x, pot

    !  pot = -half/dsqrt(x*x+0.25)
    !  pot = half*x*X
    pot = -0.800d+00*dexp(-0.150d+00*x*x)
    !  pot = zero

end

function poten_2d_(x, y) result(vxy)
    use numconst
    implicit none
    real(kind=8), intent(in)  :: x, y
    real(kind=8) :: vxy
    vxy = -two/dsqrt(x*x + one) -two/dsqrt(y*y + one) + &
            one/dsqrt((x-y)*(x-y) + one)
    return
end

module image_time_prop
    implicit none
    contains
! Subroutine:: Imaginary time-propgation
    subroutine imag_time_1d_( xmin, xmax, mass, dtim, cutof, nx, nste, evec, &
        eval, flerr )
    use numconst
    use fftpar
    implicit none

    integer(kind=4), intent(in)  :: nx, nste
    real(kind=8), intent(in)  :: xmin, xmax, mass, cutof, dtim
    real(kind=8), intent(out) :: evec(nx, nste), eval(nste)

    character(len=40), intent(in) :: flerr

    real(kind=8), allocatable :: pgrd(:), ked(:), ped(:), rcof(:)
    real(kind=8), allocatable :: gsini(:), fftgs(:), gsvec(:,:)
    real(kind=8), allocatable :: rmat(:,:), revl(:), guesn(:,:)
    real(kind=8), allocatable :: prjwf(:)
    complex(kind=8), allocatable :: cfft(:)

    real(kind=8) :: x, dx, lenx, kgrid, dpx, ppx2
    real(kind=8) :: dnrm, eng, pot, dsum, dsm2

    integer(kind=4) :: ix, indx, iste, jste, itim
    integer(kind=8) :: plan_backward, plan_forward

    character(len=39) :: dash

    dash="----------------------------------------"

    ! -----------------------------------------------------------!
    lenx = xmax - xmin             !  box length along x       !
    dpx = twpi/lenx               !  differential element px  !
    dx = twpi/(nx*dpx)           !  differential element x   !
    ppx2 = half*dpx*dpx            !  differential element kex !
    ! -----------------------------------------------------------!

    1001 format(a,f12.5,a,f12.5,a,i5)

    write(*,'(a,a)') dash, dash
    write(*,'(a40)') "variables for the calculation"
    write(*,'(a,a)') dash, dash
    write(*,1001) " | xmin |",xmin, " | xmax | ",xmax, &
            " | nx grid | ", nx
    write(*,'(a,a)') dash, dash
    write(*,'(a,e14.7)') " Cutoff error for eigenvalues :: ", cutof
    write(*,'(a,e14.7)') " Differential element of time :: ", dtim
    write(*,'(a,a)') dash, dash

    !---------- diagonal kinetic-energy along x-direction ----------

    allocate(ked(nx), ped(nx))
    allocate(pgrd(nx))

    call momentum_arry_(nx, dpx, pgrd)

    ked = zero
    do ix = 1, nx
        x = xmin + (ix-one)*dx
        ped(ix) = pot(x)
        ked(ix) = half*pgrd(ix)*pgrd(ix)
    enddo

    !--------------- define the guess vectors --------------

    allocate(gsvec(nx, nste))

    gsvec = zero
    call guess_vec_(xmin, dx, nx, nste, gsvec)

    open(11,file="guesvc.txt")
    do ix = 1, nx
        x = xmin + (ix-one)*dx
        write(11,*) x, (gsvec(ix,iste),iste=1,nste)
    enddo
    close(11)

    allocate(rcof(nx), fftgs(nx), cfft(nx/2 + 1))
    allocate(gsini(nx), prjwf(nx))

    evec = zero; evec = gsvec

    !-------------- Make plan for forward FFT and backward FFT ----------------


    dnrm = one
    open(11,file="energy.txt")
    !  do 200 while (dsqrt(dnrm) .gt. cutof)
    do itim = 1, 5000
        do iste = 1, nste
            do ix = 1, nx
                gsini(ix) = evec(ix,iste)
                fftgs(ix) = dexp(-half*ped(ix)*dtim)*gsini(ix)
            enddo

            call dfftw_plan_dft_r2c_1d_(plan_forward, nx, fftgs, cfft, FFTW_ESTIMATE)
            call dfftw_execute_dft_r2c_(plan_forward, fftgs, cfft)

            do ix = 1, nx/2 + 1
                cfft(ix) = dexp(-ked(ix)*dtim)*cfft(ix)
            enddo

            call dfftw_plan_dft_c2r_1d_(plan_backward, nx, cfft, fftgs, FFTW_ESTIMATE)
            call dfftw_execute_dft_c2r_(plan_backward, cfft, fftgs)

            do ix = 1, nx
                rcof(ix) = dexp(-half*ped(ix)*dtim)*(fftgs(ix)/(nx*one))
                gsvec(ix,iste) = rcof(ix)
            enddo
        enddo

        call guess_vecn_(xmin, dx, nx, nste, gsvec, evec)

        do iste = 1, nste
            do ix = 1, nx
                rcof(ix) = evec(ix,iste)
            enddo
            call dfftw_plan_dft_r2c_1d_(plan_forward, nx, rcof, cfft, FFTW_ESTIMATE)
            call dfftw_execute_dft_r2c_(plan_forward, rcof, cfft)
            do ix = 1, nx/2 + 1
                cfft(ix) = ked(ix)*cfft(ix)
            enddo
            call dfftw_plan_dft_c2r_1d_(plan_backward, nx, cfft, fftgs, FFTW_ESTIMATE)
            call dfftw_execute_dft_c2r_(plan_backward, cfft, fftgs)
            eng = zero
            do ix = 1, nx
                eng = eng + (rcof(ix)*fftgs(ix))/(nx*one) + &
                        rcof(ix)*ped(ix)*rcof(ix)
            enddo
            eval(iste) = eng
            write(11,'(*(a,e22.11))') " Normalization :: ", dnrm, &
                    " State energy :: ", eval(iste)
        enddo
        write(11,*)
        write(11,*)
    enddo
    close(11)
    call dfftw_destroy_plan_(plan_forward)
    call dfftw_destroy_plan_(plan_backward)
    call dfftw_cleanup_()
    end

! Evaluating eigen-values and eigen vectors for real-symmetric Hamiltonian
    subroutine real_symm_matrx_dsyev_(H, ndim, eigvals)
    IMPLICIT NONE
    CHARACTER(LEN = 1), PARAMETER  :: JOBZ="V"
    CHARACTER(LEN = 1), PARAMETER  :: UPLO="L"
    REAL(KIND = 8), ALLOCATABLE :: WORK(:)
    INTEGER(KIND = 4) :: NDIM,INFO
    INTEGER(KIND = 4) :: LDA, LDWORK

    REAL(KIND = 8) :: H(NDIM,NDIM)
    REAL(KIND = 8) :: EIGVALS(NDIM)

    LDA = NDIM
    LDWORK = 3*NDIM - 1

    ALLOCATE(WORK(LDWORK))
    CALL DSYEV(JOBZ, UPLO, NDIM, H, LDA, EIGVALS, WORK, LDWORK, INFO)
    DEALLOCATE(WORK)
    END

! Initial guess vector for Imaginary time-propagation
    subroutine guess_vec_(xmin, dx, nx, nste, gsvec)
    use numconst
    implicit none
    integer(kind=4), intent(in) :: nx, nste
    real(kind=8), intent(in) :: xmin, dx
    real(kind=8), intent(out) :: gsvec(nx,nste)
    real(kind=8), allocatable :: rmat(:,:), revl(:)
    real(kind=8), allocatable :: guesn(:,:)
    integer(kind=4) :: iste, jste, ix
    real(kind=8) :: dsum, x

    do iste = 1, nste
        do ix = 1, nx
            x = xmin + (ix-one)*dx
            gsvec(ix,iste) = (x**(iste-one))*dexp(-x*x)
        enddo
    enddo

    allocate(rmat(nste, nste), revl(nste))
    do iste = 1, nste
        do jste = 1, nste
            dsum = zero
            do ix = 1, nx
                dsum = dsum + gsvec(ix,iste)*gsvec(ix,jste)
            enddo
            rmat(iste,jste) = dsum
        enddo
    enddo

    revl = zero
    call real_symm_matrx_dsyev_(rmat, nste, revl)

    allocate(guesn(nx,nste))

    guesn = zero
    do iste = 1, nste
        do ix = 1, nx
            do jste = 1, nste
                guesn(ix,iste) = guesn(ix,iste) + (rmat(jste,iste)*&
                        gsvec(ix,jste))/dsqrt(revl(iste))
            enddo
        enddo
    enddo

    gsvec = zero
    do ix = 1, nx
        do iste = 1, nste
            gsvec(ix,iste) = guesn(ix,nste+1-iste)
        enddo
    enddo

    deallocate(rmat, revl, guesn)
end

! Evaluating the Guess-vector for Imaginary time-propgation
    subroutine guess_vecn_(xmin, dx, nx, nste, gsvec, evec)
    use numconst
    implicit none

    integer(kind=4), intent(in) :: nx, nste
    real(kind=8), intent(in) :: xmin, dx
    real(kind=8), intent(in)  :: gsvec(nx,nste)
    real(kind=8), intent(out) :: evec(nx,nste)
    real(kind=8), allocatable :: rmat(:,:), revl(:)
    real(kind=8), allocatable :: guesn(:,:)

    integer(kind=4) :: iste, jste, ix
    real(kind=8) :: dsum, x

    allocate(rmat(nste, nste), revl(nste))
    rmat = zero
    do iste = 1, nste
        do jste = 1, nste
            dsum = zero
            do ix = 1, nx
                dsum = dsum + gsvec(ix,iste)*gsvec(ix,jste)
            enddo
            rmat(iste,jste) = dsum
        enddo
    enddo

    revl = zero
    call real_symm_matrx_dsyev_(rmat, nste, revl)
    allocate(guesn(nx,nste))
    guesn = zero
    do iste = 1, nste
        do ix = 1, nx
            do jste = 1, nste
                guesn(ix,iste) = guesn(ix,iste) + (rmat(jste,iste)*&
                        gsvec(ix,jste))/dsqrt(revl(iste))
            enddo
        enddo
    enddo

    evec = zero
    do ix = 1, nx
        do iste = 1, nste
            evec(ix,iste) = guesn(ix,nste+1-iste)
        enddo
    enddo
    deallocate(rmat, revl, guesn)
end

! Momentum grid points
    subroutine momentum_arry_(nr, dpr, parry)
    use numconst
    implicit none

    integer(kind=4), intent(in) :: nr
    real(kind=8), intent(in)    :: dpr
    real(kind=8), intent(out)   :: parry(nr)
    integer(kind=4) :: ir

    if ( mod(nr,2) == 0 ) then
        parry = zero
        DO ir = 1, 1 + Nr/2
            parry(ir) = (ir-one)*dpr
        END DO
        DO Ir = 1, Nr/2 -1
            parry(Ir+1+Nr/2) = -parry(1-Ir+Nr/2)
        END DO
    else
        parry = zero
        DO Ir = 1, 1 + Nr/2
            parry(Ir) = (Ir-one)*dpr
        END DO
        DO Ir = 1, Nr/2
            parry(Ir+1+Nr/2) = -parry(-Ir+Nr/2 +2)
        END DO
    endif

    end

    subroutine imag_timeprop_2d_frgrd_( xmin, xmax, nx, ymin, ymax,  ny, mass, &
            cutof, dtim, eval, evec )
        use numconst
        use fftpar
        implicit none

        integer(kind=4), intent(in) :: nx, ny
        real(kind=8), intent(in)    :: xmin, xmax, ymin, ymax
        real(kind=8), intent(in)    :: mass, cutof, dtim

        real(kind=8), intent(out)   :: eval, evec(nx*ny)

        complex(kind=8), allocatable :: gsini(:,:), ifftgs(:,:), fftgs(:,:)
        complex(kind=8), allocatable :: kexp(:,:), vexp(:,:)
        complex(kind=8), allocatable :: vexgs(:,:), pexgs(:,:)
        complex(kind=8), allocatable :: ked(:,:), ped(:,:)

        real(kind=8), allocatable :: pgrd(:), kex(:), key(:)
        real(kind=8), allocatable :: prjmt(:), prjwf(:)

        real(kind=8) :: x, dx, y, dy, totime, poten_2d_
        real(kind=8) :: dpx, dpy, alp1, alp2, alp3

        integer(kind=4) :: ix, jx, iy, jy, ixy, iste, indx
        integer(kind=8) :: plan_backward, plan_forward

        complex(kind=8) :: eng, error

        character*40 :: dash

        dash="-----------------------------------------"

        !      -----------------------------------------------------------!
        dx = (xmax-xmin)/(nx*one)    !  differential element x   !
        dy = (ymax-ymin)/(ny*one)    !  differential element y   !
        dpx = twpi/(xmax-xmin)        !  differential element px  !
        dpy = twpi/(ymax-ymin)        !  differential element py  !
        !      -----------------------------------------------------------!

        write(*,'(a,a)') dash, dash
        write(*,'(a,2f12.7)') "DX1 e1, DX2 e2     :: ", dx, dy
        write(*,'(a,2f12.7)') "DPX1 e1, DPX2 e2   :: ", dpx, dpy
        write(*,'(a,a)') dash, dash

        allocate(kex(nx), key(ny))

        !     ---------- diagonal kinetic-energy along x-direction ----------

        allocate(pgrd(nx))

        if ( mod(nx,2) == 0 ) then

            pgrd = zero
            DO Ix = 1, 1 + NX/2
                pgrd(Ix) = (Ix-one)*dpx
            END DO

            DO Ix = 1, NX/2 -1
                pgrd(Ix+1+NX/2) = -pgrd(1-Ix+NX/2)
            END DO

        else

            pgrd = zero
            DO Ix = 1, 1 + NX/2
                pgrd(Ix) = (Ix-one)*dpx
            END DO

            DO Ix = 1, NX/2
                pgrd(Ix+1+NX/2) = -pgrd(-Ix+NX/2 +2)
            END DO

        endif

        kex = zero
        do ix = 1, nx
            kex(ix) = half*pgrd(ix)*pgrd(ix)
        enddo

        !     ---------- diagonal kinetic-energy along y-direction ----------

        deallocate(pgrd)
        allocate(pgrd(ny))

        if ( mod(ny,2) == 0 ) then

            pgrd = zero
            DO Iy = 1, 1 + Ny/2
                pgrd(Iy) = (Iy-one)*dpy
            END DO

            DO Iy = 1, Ny/2 -1
                pgrd(Iy+1+Ny/2) = -pgrd(1-Iy+Ny/2)
            END DO

        else

            pgrd = zero
            DO Iy = 1, 1 + Ny/2
                pgrd(Iy) = (Iy-one)*dpy
            END DO

            DO Iy = 1, Ny/2
                pgrd(Iy+1+Ny/2) = -pgrd(-Iy+Ny/2 +2)
            END DO

        endif

        key = zero
        do iy = 1, ny
            key(iy) = half*pgrd(iy)*pgrd(iy)
        enddo

        deallocate(pgrd)
        allocate(ked(nx,ny), ped(nx,ny))
        allocate(kexp(nx,ny), vexp(nx,ny))

        ked = (zero, zero); ped = (zero, zero)
        kexp = (zero, zero); vexp = (zero, zero)

        !     ---- exponetial and the diagonal matrix for time-propagation ---
        !     ----------- of kinetic and potential energy matrix  ------------

        open(11,file="potn.txt")
        do ix = 1, nx
            x = xmin + (ix-one)*dx
            do iy = 1, ny
                y = ymin + (iy-one)*dy

                ked(ix,iy) = kex(ix) + key(iy)
                ped(ix,iy) = poten_2d_(x,y)

                kexp(ix,iy) = cdexp(-ked(ix,iy)*dtim)
                vexp(ix,iy) = cdexp(-half*ped(ix,iy)*dtim)

                write(11,'(*(f12.6))') x, y, real(ked(ix,iy)), &
                        real(ped(ix,iy))
            enddo
            write(11,*)
        enddo
        close(11)

        deallocate(kex, key)

        !     --- imaginary time-propagation for the ground sate wavefunction --

        allocate(gsini(nx,ny), prjmt(nx*ny), prjwf(nx*ny))
        allocate(vexgs(nx,ny), pexgs(nx,ny))
        allocate(fftgs(nx,ny), ifftgs(nx,ny))

        !     ------- initial guess for the ground state wavefunction --------

        gsini = (zero, zero); eng = (zero, zero)
        do ix = 1, nx
            x = xmin + (ix-one)*dx
            do iy = 1, ny
                y = ymin + (iy-one)*dy
                gsini(ix,iy) = (x**(iste-one))*dexp(-(x*x + y*y))
                eng = eng + conjg(gsini(ix,iy))*gsini(ix,iy)
            enddo
        enddo

        gsini = gsini/(dsqrt(cdabs(eng)))

        error = one
        do while ( real(error) > cutof )

            !     -- diagonal potential expential multiplication to the guess fnc --

            vexgs = (zero, zero); ixy = 0
            do ix = 1, nx
                do iy = 1, ny
                    ixy = ixy + 1
                    vexgs(ix,iy) = vexp(ix,iy)*gsini(ix,iy)
                enddo
            enddo

            !     --------- Fourier transform of the exp(-0.5*v*dt)*psi(x,t) --------

            fftgs = (zero, zero)
            call dfftw_plan_dft_2d_ ( plan_forward, nx, ny, vexgs, fftgs, &
                    FFTW_FORWARD, FFTW_ESTIMATE )
            call dfftw_execute_ ( plan_forward )

            !     ---- diagonal kinetic exponential multiplication to the Fourier ---
            !     ------------------- transformed wavefucntion ----------------------

            pexgs = (zero, zero)
            do ix = 1, nx
                do iy = 1, ny
                    pexgs(ix,iy) = kexp(ix,iy)*fftgs(ix,iy)
                enddo
            enddo

            !     ------------------ Inverse Fourier transform of the ---------------
            !     ------------- exp(-p2/m*dt)*fft(exp(-0.5*v*dt)*psi(x,t)) ----------

            ifftgs = (zero, zero)
            call dfftw_plan_dft_2d_ ( plan_backward, nx, ny, pexgs, ifftgs, &
                    FFTW_BACKWARD, FFTW_ESTIMATE )
            call dfftw_execute_ ( plan_backward )

            !     -------- division by number of points in x and y direction --------

            ifftgs = ifftgs/(nx*ny*one)

            !     -- exp(-0.5*v*dt)ifft(exp(-p2/m*dt)*fft(exp(-0.5*v*dt)psi(x,t))) --

            fftgs = (zero, zero)
            do ix = 1, nx
                do iy = 1, ny
                    fftgs(ix,iy) = vexp(ix,iy)*ifftgs(ix,iy)
                enddo
            enddo

            eng = (zero, zero)
            do ix = 1, nx
                do iy = 1, ny
                    eng = eng + conjg(fftgs(ix,iy))*fftgs(ix,iy)
                enddo
            enddo

            fftgs = fftgs/dsqrt(cdabs(eng))

            error = (zero, zero); eng = (zero, zero)
            do ix = 1, nx
                do iy = 1, ny
                    error = error + (gsini(ix,iy) - fftgs(ix,iy))**two
                    eng = eng + conjg(fftgs(ix,iy))*fftgs(ix,iy)
                enddo
            enddo

            error = dsqrt(cdabs(error))
            gsini = (zero, zero); gsini = fftgs

        enddo

        !     --------- Fourier transform of the ground state fucntion ----------

        ifftgs = (zero, zero)
        call dfftw_plan_dft_2d_ ( plan_forward, nx, ny, fftgs, ifftgs, &
                FFTW_FORWARD, FFTW_ESTIMATE )

        call dfftw_execute_ ( plan_forward )

        !     ------- diagonal kinetic energy with Fourier transformed GS -------

        pexgs = (zero, zero)
        do ix = 1, nx
            do iy = 1, ny
                pexgs(ix,iy) = ked(ix,iy)*ifftgs(ix,iy)
            enddo
        enddo

        !     ------- inverse Fourier transform for the energy calculation ------

        ifftgs = (zero, zero)
        call dfftw_plan_dft_2d_ ( plan_backward, nx, ny, pexgs, ifftgs, &
                FFTW_BACKWARD, FFTW_ESTIMATE )

        call dfftw_execute_ ( plan_backward )

        !     -------- division by number of points in x and y direction --------

        ifftgs = ifftgs/(nx*ny*one)

        !     ------------- energy with imaginary time-propagation --------------

        eng = (zero, zero)
        do ix = 1, nx
            do iy = 1, ny
                eng = eng + gsini(ix,iy)*real(ifftgs(ix,iy)) + &
                        gsini(ix,iy)*ped(ix,iy)*gsini(ix,iy)
            enddo
        enddo

        eval = real(eng)

        ixy = 0
        do ix = 1, nx
            do iy = 1, ny
                ixy = ixy + 1
                evec(ixy) = gsini(ix,iy)
            enddo
        enddo

        call dfftw_destroy_plan_ ( plan_forward )
        call dfftw_destroy_plan_ ( plan_backward )

    end
end module image_time_prop

module time_prop_ttp

end module time_prop_ttp

module time_prop_Trotter
    contains
    !     -------------------------------------------------------------------------
    subroutine timeprop_2d_split_frgrd_( xmin, xmax, nx, ymin, ymax, &
            ny, mass, epsln, omega,     &
            phase, psi0 )
    !     -------------------------------------------------------------------------
    use numconst
    use fftpar
    implicit none

    integer(kind=4), intent(in) :: nx, ny
    real(kind=8), intent(in)    :: xmin, xmax, ymin, ymax
    real(kind=8), intent(in)    :: epsln, omega, phase
    real(kind=8), intent(in)    :: mass
    real(kind=8), intent(in)    :: psi0(nx*ny)

    integer(kind=4), parameter :: nstep = 200

    complex(kind=8), allocatable :: gsini(:,:), ifftgs(:,:), fftgs(:,:)
    complex(kind=8), allocatable :: kexp(:,:), vexp(:,:)
    complex(kind=8), allocatable :: vexgs(:,:), pexgs(:,:)
    complex(kind=8), allocatable :: ked(:,:), tmcof(:,:)

    real(kind=8), allocatable :: pgrd(:), kex(:), key(:)

    real(kind=8) :: x, dx, y, dy, poten_2d_
    real(kind=8) :: dpx, dpy, time, dtim, totime

    integer(kind=4) :: ix, jx, iy, jy, ixy, itim, ntim
    integer(kind=8) :: plan_backward, plan_forward, icof

    complex(kind=8) :: eng, dipl

    character*40 :: dash

    !      -----------------------------------------------------------!
    dx = (xmax-xmin)/(nx*one)    !  differential element x   !
    dy = (ymax-ymin)/(ny*one)    !  differential element y   !
    dpx = twpi/(xmax-xmin)        !  differential element px  !
    dpy = twpi/(ymax-ymin)        !  differential element py  !
    !      -----------------------------------------------------------!

    totime = twpi/(1.1040d0*epsln)
    ntim = int(totime*40.0d0)
    dtim = totime/(ntim-one)

    write(*,*) totime, dtim
    stop

    allocate(kex(nx), key(ny))

    !     ---------- diagonal kinetic-energy along x-direction ----------

    allocate(pgrd(nx))

    !     --- grid variable for odd and even number of grids along x ---

    if ( mod(nx,2) == 0 ) then

        pgrd = zero
        DO Ix = 1, 1 + NX/2
            pgrd(Ix) = (Ix-one)*dpx
        END DO

        DO Ix = 1, NX/2 -1
            pgrd(Ix+1+NX/2) = -pgrd(1-Ix+NX/2)
        END DO

    else

        pgrd = zero
        DO Ix = 1, 1 + NX/2
            pgrd(Ix) = (Ix-one)*dpx
        END DO

        DO Ix = 1, NX/2
            pgrd(Ix+1+NX/2) = -pgrd(-Ix+NX/2 +2)
        END DO

    endif

    kex = zero
    do ix = 1, nx
        kex(ix) = half*pgrd(ix)*pgrd(ix)
    enddo

    !     --- grid variable for odd and even number of grids along y ---

    deallocate(pgrd)
    allocate(pgrd(ny))

    if ( mod(ny,2) == 0 ) then

        pgrd = zero
        DO Iy = 1, 1 + Ny/2
            pgrd(Iy) = (Iy-one)*dpy
        END DO

        DO Iy = 1, Ny/2 -1
            pgrd(Iy+1+Ny/2) = -pgrd(1-Iy+Ny/2)
        END DO

    else

        pgrd = zero
        DO Iy = 1, 1 + Ny/2
            pgrd(Iy) = (Iy-one)*dpy
        END DO

        DO Iy = 1, Ny/2
            pgrd(Iy+1+Ny/2) = -pgrd(-Iy+Ny/2 +2)
        END DO

    endif

    key = zero
    do iy = 1, ny
        key(iy) = half*pgrd(iy)*pgrd(iy)
    enddo

    deallocate(pgrd)

    !     ---- exponetial and the diagonal matrix for time-propagation ---
    !     ----------- of kinetic and potential energy matrix  ------------

    allocate(ked(nx,ny))
    allocate(kexp(nx,ny))

    ked = (zero, zero)
    kexp = (zero, zero)

    do ix = 1, nx
        do iy = 1, ny
            ked(ix,iy) = kex(ix) + key(iy)
            kexp(ix,iy) = cdexp(-io*ked(ix,iy)*dtim)
        enddo
    enddo

    deallocate(kex, key)
    allocate(gsini(nx,ny))

    gsini = (zero, zero); ixy = 0
    do ix = 1, nx
        do iy = 1, ny
            ixy = ixy + 1
            gsini(ix,iy) = psi0(ixy)
        enddo
    enddo

    allocate(vexp(nx,ny), vexgs(nx,ny))
    allocate(fftgs(nx,ny), ifftgs(nx,ny), pexgs(nx,ny))
    allocate(tmcof(nx*ny,ntim/nstep))

    tmcof = (zero, zero)
    open(11,file="tmengdpl.txt")
    do itim = 1, ntim

        time = (itim - one)*dtim

        !     -- diagonal potential expential multiplication to the guess fnc --

        vexgs = (zero, zero)
        do ix = 1, nx
            x = xmin + (ix-one)*dx
            do iy = 1, ny
                y = ymin + (iy-one)*dy
                vexp(ix,iy) = cdexp(-io*half*(poten_2d_(x,y)+&
                        (x+y)*epsln*dcos(omega*time))*dtim)
                vexgs(ix,iy) = vexp(ix,iy)*gsini(ix,iy)
            enddo
        enddo

        !     --------- Fourier transform of the exp(-i*v*dt/2)*psi(x,t) --------

        fftgs = (zero, zero)
        call dfftw_plan_dft_2d_ ( plan_forward, nx, ny, vexgs, fftgs, &
                FFTW_FORWARD, FFTW_ESTIMATE )
        call dfftw_execute_ ( plan_forward )

        !     ---- diagonal kinetic exponential multiplication to the Fourier ---
        !     ------------------- transformed wavefucntion ----------------------

        pexgs = (zero, zero)
        do ix = 1, nx
            do iy = 1, ny
                pexgs(ix,iy) = kexp(ix,iy)*fftgs(ix,iy)
            enddo
        enddo

        !     ------------------ Inverse Fourier transform of the ---------------
        !     ------------- exp(-p2/m*dt)*fft(exp(-i*v*dt/2)*psi(x,t)) ----------

        ifftgs = (zero, zero)
        call dfftw_plan_dft_2d_ ( plan_backward, nx, ny, pexgs, ifftgs, &
                FFTW_BACKWARD, FFTW_ESTIMATE )
        call dfftw_execute_ ( plan_backward )

        !     -------- division by number of points in x and y direction --------

        ifftgs = ifftgs/(nx*ny*one)

        !     -- exp(-0.5*v*dt)ifft(exp(-p2/m*dt)*fft(exp(-0.5*v*dt)psi(x,t))) --

        fftgs = (zero, zero)
        do ix = 1, nx
            do iy = 1, ny
                fftgs(ix,iy) = vexp(ix,iy)*ifftgs(ix,iy)
            enddo
        enddo

        gsini = (zero, zero); gsini = fftgs

        !      ------ time-dependent 2d coefficients for properties ------

        if ( mod(itim ,nstep) == 0 ) then
            icof = icof + 1; ixy = 0
            do ix = 1, nx
                do iy = 1, ny
                    ixy = ixy + 1
                    tmcof(ixy, icof) = gsini(ix,iy)
                enddo
            enddo
        endif

        !     --------- Fourier transform of the ground state fucntion ----------

        ifftgs = (zero, zero)
        call dfftw_plan_dft_2d_ ( plan_forward, nx, ny, fftgs, ifftgs, &
                FFTW_FORWARD, FFTW_ESTIMATE )

        call dfftw_execute_ ( plan_forward )

        !     ------- diagonal kinetic energy with Fourier transformed GS -------

        pexgs = (zero, zero)
        do ix = 1, nx
            do iy = 1, ny
                pexgs(ix,iy) = ked(ix,iy)*ifftgs(ix,iy)
            enddo
        enddo

        !     ------- inverse Fourier transform for the energy calculation ------

        ifftgs = (zero, zero)
        call dfftw_plan_dft_2d_ ( plan_backward, nx, ny, pexgs, ifftgs, &
                FFTW_BACKWARD, FFTW_ESTIMATE )

        call dfftw_execute_ ( plan_backward )

        !     -------- division by number of points in x and y direction --------

        ifftgs = ifftgs/(nx*ny*one)

        !     ------------- energy with imaginary time-propagation --------------

        eng = (zero, zero); dipl = (zero, zero)
        do ix = 1, nx
            x = xmin + (ix-one)*dx
            do iy = 1, ny
                y = ymin + (iy-one)*dy

                eng = eng + conjg(fftgs(ix,iy))*(ifftgs(ix,iy) + &
                        poten_2d_(x,y)*fftgs(ix,iy))
                dipl = dipl + conjg(fftgs(ix,iy))*(x+y)*fftgs(ix,iy)

            enddo
        enddo

        write(11,*) time, epsln*dcos(omega*time), real(eng), &
                real(eng) + real(dipl)*epsln*dcos(omega*time), &
                real(dipl)

    enddo
    close(11)

    call dfftw_destroy_plan_ ( plan_forward )
    call dfftw_destroy_plan_ ( plan_backward )

    end

        subroutine tprop_KH_dvr2d_cart_ ( xmin, xmax, ymin, ymax, nx, ny, n2d, natm, cord,      &
                epsln, omega, phase, lmda, totime, ntim, nstrt, ncw, noff, noptc,      &
                nstep, psi0, flwfnc, fleng, fltden, fltpot )
            use numconst
            implicit none

            real(kind=8), intent(in) :: xmin, xmax, ymin, ymax, lmda
            real(kind=8), intent(in) :: epsln, omega, phase, totime
            real(kind=8), intent(in) :: psi0(n2d), cord(natm,2)
            integer(kind=4), intent(in) :: nx, ny, n2d, ntim, nstep, natm
            integer(kind=4), intent(in) :: nstrt, ncw, noff, noptc
            character(len=80), intent(in) :: flwfnc, fleng, fltden, fltpot
            real(kind=8), allocatable :: kmtx(:,:), kmty(:,:), kex(:), key(:)
            real(kind=8), allocatable :: vx(:), force(:), vcap(:), tdiag(:)
            real(kind=8), allocatable :: tmcof(:,:), tmptn(:,:)
            complex(kind=8), allocatable :: ztmtx(:,:), ztmty(:,:), vxexp(:)
            complex(kind=8), allocatable :: cofi(:), coff(:), cofp(:)
            complex(kind=8), allocatable :: cfx1(:), cfx2(:), cfy1(:), cfy2(:)
            integer(kind=4) :: ix, iy, jx, jy, i2d, nxy, itim, icof
            real(kind=8) :: x, dx, y, dy, rvl, rmax, poten, dpotn
            real(kind=8) :: tau, plstr, tcw, tof, tff, plend, tim2, epls
            real(kind=8) :: time, dtim, alpx, alpy, wfnrm
            complex(kind=8) :: hlfdt, teng, indpl

            rmax = dsqrt(xmax*xmax + ymax*ymax)
            dx = (xmax - xmin)/(nx - one)
            dy = (ymax - ymin)/(ny - one)
            dtim = totime/(ntim-one)
            tau = twpi/omega
            plstr = pi/(two*nstrt*tau)
            tcw = nstrt*tau
            tof = ncw*tau + tcw
            tff = tof + noff*tau
            plend = pi/(2.0d0*(tff-tof))

            allocate(kmtx(nx, nx), kmty(ny, ny))
            allocate(ztmtx(nx, nx), ztmty(ny, ny))
            allocate(kex(nx), key(ny))

            call Kinetic(dx, nx, kmtx)
            call CALL_DSYEV(kmtx, nx, kex)
            ztmtx = kmtx

            call Kinetic(dy, ny, kmty)
            call CALL_DSYEV(kmty, ny, key)
            ztmty = kmty
            deallocate(kmtx, kmty)
            allocate(tdiag(n2d), vcap(n2d), cofi(n2d))
            nxy = 0
            open(101,file=flwfnc)
            do ix = 1, nx
                x = xmin + (ix-one)*dx
                do iy = 1, ny
                    y = ymin + (iy-one)*dy
                    rvl = dsqrt(x*x + y*y)
                    nxy = nxy + 1
                    tdiag(nxy) = kex(ix) + key(iy)
                    cofi(nxy) = psi0(nxy)
                    vcap(nxy) = zero
                    if ( rvl > (rmax - 5.000d+00) ) &
                            vcap(nxy) = lmda*(rvl - (rmax - 5.00D+00))**2.0d0
                    write(101,*) x, y, psi0(nxy), vcap(nxy)
                enddo
                write(101,*)
            enddo
            close(101)

            deallocate(kex, key)
            allocate(vx(n2d), force(n2d), vxexp(n2d), coff(n2d), cofp(n2d))
            allocate(cfx1(nx), cfx2(nx), cfy1(ny), cfy2(ny))
            allocate(tmcof(n2d, ntim/nstep), tmptn(n2d,ntim/nstep))

            tmcof = zero; tmptn = zero
            hlfdt = half*dtim
            icof = 0
            open(11,file=fleng)
            do itim = 1, ntim
                time = (itim - one)*dtim
                epls = epsln
                tim2 = time - tof
                if(time <= tcw) epls = epls*dsin(plstr*time)**2.0D+00
                if(time >= tof) epls = epls*dcos(plend*tim2)**2.0D+00
                if(time >= tff) epls = zero

                alpx = (epls*dcos(omega*time))/(omega*omega)
                alpy = (epls*dsin(omega*time))/(omega*omega)

                nxy = 0
                do ix = 1, nx
                    x = xmin + (ix-one)*dx
                    do iy = 1, ny
                        y = ymin + (iy-one)*dy
                        nxy = nxy + 1

                        vx(nxy) = poten(natm, cord, x-alpx, y-alpy)
                        force(nxy) = dpotn(natm, cord, x-alpx, y-alpy)
                        vxexp(nxy) = cdexp(-io*hlfdt*vx(nxy))
                        coff(nxy) = vxexp(nxy)*cofi(nxy)

                    enddo
                enddo

                do iy = 1, ny
                    do ix = 1, nx
                        cfx1(ix) = coff(ny*(ix-1) + iy)
                    enddo

                    call zgemm ('T', 'N', nx, 1, nx, aval, ztmtx, nx, &
                            cfx1, nx, bval, cfx2, nx )

                    do ix = 1, nx
                        cofi(ny*(ix-1) + iy) = cfx2(ix)
                    enddo
                enddo

                 nxy = 0
                do ix = 1, nx
                    do iy = 1, ny
                        cfy1(iy) = cofi(nx*(ix-1) + iy)
                    enddo

                    call zgemm ('T', 'N', ny, 1, ny, aval, ztmty, ny, &
                            cfy1, ny, bval, cfy2, ny )

                    do iy = 1, ny
                        nxy = nxy + 1
                        coff(nxy) = cfy2(iy)
                    enddo
                enddo

                do i2d = 1, n2d
                    cofi(i2d) = cdexp(-io*tdiag(i2d)*dtim)*coff(i2d)
                enddo

                 nxy = 0
                do ix = 1, nx
                    do iy = 1, ny
                        cfy1(iy) = cofi(nx*(ix-1) + iy)
                    enddo

                    call zgemm ('N', 'N', ny, 1, ny, aval, ztmty, ny, &
                            cfy1, ny, bval, cfy2, ny )

                    do iy = 1, ny
                        nxy = nxy + 1
                        coff(nxy) = cfy2(iy)
                    enddo
                enddo

                do iy = 1, ny
                    do ix = 1, nx
                        cfx1(ix) = coff(ny*(ix-1) + iy)
                    enddo

                    call zgemm ('N', 'N', nx, 1, nx, aval, ztmtx, nx, &
                            cfx1, nx, bval, cfx2, nx )

                    do ix = 1, nx
                        cofi(ny*(ix-1) + iy) = cfx2(ix)
                    enddo
                enddo

                do i2d = 1, n2d
                    coff(i2d) = dexp(-vcap(i2d)*dtim)*vxexp(i2d)*cofi(i2d)
                enddo

                cofi = coff
                if ( mod(itim ,nstep) == 0 ) then
                    icof = icof + 1
                    do i2d = 1, n2d
                        tmcof(i2d, icof) = cdabs(conjg(cofi(i2d))*cofi(i2d))
                        tmptn(i2d, icof) = vx(i2d)
                    enddo
                endif

                do iy = 1, ny
                    do ix = 1, nx
                        cfx1(ix) = coff(ny*(ix-1) + iy)
                    enddo
                    call zgemm ('T', 'N', nx, 1, nx, aval, ztmtx, nx, &
                            cfx1, nx, bval, cfx2, nx )
                    do ix = 1, nx
                        cofp(ny*(ix-1) + iy) = cfx2(ix)
                    enddo
                enddo

                coff = (zero, zero); nxy = 0
                do ix = 1, nx
                    do iy = 1, ny
                        cfy1(iy) = cofp(nx*(ix-1) + iy)
                    enddo
                    call zgemm ('T', 'N', ny, 1, ny, aval, ztmty, ny, &
                            cfy1, ny, bval, cfy2, ny )
                    do iy = 1, ny
                        nxy = nxy + 1
                        coff(nxy) = cfy2(iy)
                    enddo
                enddo

                cofp = (zero, zero)
                do i2d = 1, n2d
                    cofp(i2d) = tdiag(i2d)*coff(i2d)
                enddo

                 nxy = 0
                do ix = 1, nx

                    do iy = 1, ny
                        cfy1(iy) = cofp(nx*(ix-1) + iy)
                    enddo

                    call zgemm ('N', 'N', ny, 1, ny, aval, ztmty, ny, &
                            cfy1, ny, bval, cfy2, ny )

                    do iy = 1, ny
                        nxy = nxy + 1
                        coff(nxy) = cfy2(iy)
                    enddo

                enddo

                cofp = (zero, zero)
                do iy = 1, ny

                    do ix = 1, nx
                        cfx1(ix) = coff(ny*(ix-1) + iy)
                    enddo

                    call zgemm ('N', 'N', nx, 1, nx, aval, ztmtx, nx, &
                            cfx1, nx, bval, cfx2, nx )

                    do ix = 1, nx
                        cofp(ny*(ix-1) + iy) = cfx2(ix)
                    enddo

                enddo

                teng = (zero, zero); wfnrm = zero; indpl = (zero, zero)
                do i2d = 1, n2d
                    teng = teng + conjg(cofi(i2d))*cofp(i2d) + &
                            conjg(cofi(i2d))*(vx(i2d)-io*vcap(i2d))*cofi(i2d)
                    wfnrm = wfnrm  + cdabs(conjg(cofi(i2d))*cofi(i2d))
                    indpl = indpl + conjg(cofi(i2d))*force(i2d)*cofi(i2d)
                enddo

                write(11,*) time, epls*dcos(omega*time), epls*dsin(omega*time),&
                        real(teng), aimag(teng), wfnrm, real(indpl)

            enddo
            close(11)

            open(12,file=fltden)
            nxy = 0
            do ix = 1, nx
                x = xmin + (ix-one)*dx
                do iy = 1, ny
                    y = ymin + (iy-one)*dy
                    nxy = nxy + 1
                    write(12,'(*(f12.7))') x, y, &
                            (tmcof(nxy,itim), itim=1,ntim/nstep)
                enddo
                write(12,*)
            enddo
            close(12)

            open(12,file=fltpot)
            nxy = 0
            do ix = 1, nx
                x = xmin + (ix-one)*dx
                do iy = 1, ny
                    y = ymin + (iy-one)*dy
                    nxy = nxy + 1
                    write(12,'(*(f12.7))') x, y, &
                            (tmptn(nxy,itim), itim=1,ntim/nstep)
                enddo
                write(12,*)
            enddo
            close(12)

            deallocate(ztmtx, ztmty, cofi, coff, cofp, cfx1, cfx2, cfy1, cfy2)
            deallocate(tdiag, vx, vxexp, vcap, force, tmcof, tmptn)

        end

    subroutine tprop_dpl_dvr2d_cart_ ( xmin, xmax, ymin, ymax, nx, ny, n2d, natm, cord, &
                epsln, omega, phase, lmda, totime, ntim, nstrt, ncw, noff, noptc,      &
                nstep, psi0, flwfnc, fleng, fltden, fltpot )
        use numconst
        implicit none

        real(kind=8), intent(in) :: xmin, xmax, ymin, ymax
        real(kind=8), intent(in) :: epsln, omega, phase
        real(kind=8), intent(in) :: lmda, totime
        real(kind=8), intent(in) :: psi0(n2d), cord(natm,2)
        integer(kind=4), intent(in) :: nx, ny, n2d, ntim, nstep, natm
        integer(kind=4), intent(in) :: nstrt, ncw, noff, noptc
        character(len=80), intent(in) :: flwfnc, fleng, fltden, fltpot

            real(kind=8), allocatable :: kmtx(:,:), kmty(:,:), kex(:), key(:)
            real(kind=8), allocatable :: vx(:), vcap(:), frce(:), tdiag(:)
            real(kind=8), allocatable :: tmcof(:,:), tmptn(:,:)
            complex(kind=8), allocatable :: ztmtx(:,:), ztmty(:,:), vxexp(:)
            complex(kind=8), allocatable :: cofi(:), coff(:), cofp(:)
            complex(kind=8), allocatable :: cfx1(:), cfx2(:), cfy1(:), cfy2(:)
            integer(kind=4) :: ix, iy, jx, jy, i2d, nxy, itim, icof
            real(kind=8) :: x, dx, y, dy, rmax, rvl, xtim, ytim, poten, dpotn
            real(kind=8) :: tau, plstr, tcw, tof, tff, plend, tim2, epls
            real(kind=8) :: time, dtim, indpl, wfnrm
            complex(kind=8) :: hlfdt, teng

            rmax = dsqrt(xmax*xmax + ymax*ymax)
            dx = (xmax - xmin)/(nx - one)
            dy = (ymax - ymin)/(ny - one)
            dtim = totime/(ntim-one)

            tau = twpi/omega
            plstr = pi/(two*nstrt*tau)
            tcw = nstrt*tau
            tof = ncw*tau + tcw
            tff = tof + noff*tau
            plend = pi/(2.0d0*(tff-tof))

            allocate(kmtx(nx, nx), kmty(ny, ny))
            allocate(ztmtx(nx, nx), ztmty(ny, ny))
            allocate(kex(nx), key(ny))

            call Kinetic(dx, nx, kmtx)
            call CALL_DSYEV(kmtx, nx, kex)
            ztmtx = kmtx

            call Kinetic(dy, ny, kmty)
            call CALL_DSYEV(kmty, ny, key)
            ztmty = kmty
            deallocate(kmtx, kmty)
            allocate(tdiag(n2d), vx(n2d), vcap(n2d), frce(n2d), cofi(n2d))
            vx = zero; vcap = zero; tdiag = zero; frce = zero
            cofi = (zero, zero)
            nxy = 0
            open(101,file=flwfnc)
            do ix = 1, nx
                x = xmin + (ix-one)*dx
                do iy = 1, ny
                    y = ymin + (iy-one)*dy
                    rvl = dsqrt(x*x + y*y)
                    nxy = nxy + 1
                    tdiag(nxy) = kex(ix) + key(iy)
                    vx(nxy) = poten(natm, cord, x, y)
                    frce(nxy) = dpotn(natm, cord, x, y)
                    cofi(nxy) = psi0(nxy)
                    vcap(nxy) = zero
                    if ( rvl > (rmax - 5.000d+00) ) &
                            vcap(nxy) = lmda*(rvl - (rmax - 5.00D+00))**2.0d0
                    write(101,*) x, y, psi0(nxy)
                enddo
                write(101,*)
            enddo
            close(101)

            deallocate(kex, key)
            allocate(vxexp(n2d), coff(n2d), cofp(n2d))
            allocate(cfx1(nx), cfx2(nx), cfy1(ny), cfy2(ny))
            allocate(tmcof(n2d, ntim/nstep), tmptn(n2d, ntim/nstep))

            hlfdt = half*dtim
            icof = 0
            open(11,file=fleng)
            do itim = 1, ntim
                time = (itim - one)*dtim
                epls = epsln
                tim2 = time - tof
                if(time <= tcw) epls = epls*dsin(plstr*time)**2.0D+00
                if(time >= tof) epls = epls*dcos(plend*tim2)**2.0D+00
                if(time >= tff) epls = zero

                nxy = 0
                do ix = 1, nx
                    x = xmin + (ix-one)*dx
                    xtim = x*dcos(omega*time)
                    do iy = 1, ny
                        y = ymin + (iy-one)*dy
                        ytim = y*dsin(omega*time)
                        nxy = nxy + 1
                        vxexp(nxy) = cdexp(-io*hlfdt*(vx(nxy) + (xtim+ytim)*epls))
                        coff(nxy) = vxexp(nxy)*cofi(nxy)
                    enddo
                enddo

                do iy = 1, ny
                    do ix = 1, nx
                        cfx1(ix) = coff(ny*(ix-1) + iy)
                    enddo
                    call zgemm ('T', 'N', nx, 1, nx, aval, ztmtx, nx, &
                            cfx1, nx, bval, cfx2, nx )
                    do ix = 1, nx
                        cofi(ny*(ix-1) + iy) = cfx2(ix)
                    enddo
                enddo

                 nxy = 0
                do ix = 1, nx

                    do iy = 1, ny
                        cfy1(iy) = cofi(nx*(ix-1) + iy)
                    enddo

                    call zgemm ('T', 'N', ny, 1, ny, aval, ztmty, ny, &
                            cfy1, ny, bval, cfy2, ny )

                    do iy = 1, ny
                        nxy = nxy + 1
                        coff(nxy) = cfy2(iy)
                    enddo

                enddo

                cofi = (zero, zero)
                do i2d = 1, n2d
                    cofi(i2d) = cdexp(-io*tdiag(i2d)*dtim)*coff(i2d)
                enddo

                 nxy = 0
                do ix = 1, nx
                    do iy = 1, ny
                        cfy1(iy) = cofi(nx*(ix-1) + iy)
                    enddo

                    call zgemm ('N', 'N', ny, 1, ny, aval, ztmty, ny, &
                            cfy1, ny, bval, cfy2, ny )

                    do iy = 1, ny
                        nxy = nxy + 1
                        coff(nxy) = cfy2(iy)
                    enddo
                enddo

                do iy = 1, ny
                    cfx1 = (zero, zero)
                    do ix = 1, nx
                        cfx1(ix) = coff(ny*(ix-1) + iy)
                    enddo

                    cfx2 = (zero, zero)
                    call zgemm ('N', 'N', nx, 1, nx, aval, ztmtx, nx, &
                            cfx1, nx, bval, cfx2, nx )

                    do ix = 1, nx
                        cofi(ny*(ix-1) + iy) = cfx2(ix)
                    enddo
                enddo

                coff = (zero, zero)
                do i2d = 1, n2d
                    coff(i2d) = dexp(-vcap(i2d)*dtim)*vxexp(i2d)*cofi(i2d)
                enddo

                cofi = (zero, zero); cofi = coff

                if ( mod(itim ,nstep) == 0 ) then
                    icof = icof + 1
                    do i2d = 1, n2d
                        tmcof(i2d, icof) = cdabs(conjg(cofi(i2d))*cofi(i2d))
                        tmptn(i2d, icof) = vx(i2d)
                    enddo
                endif

                cofp = (zero, zero)
                do iy = 1, ny

                    cfx1 = (zero, zero)
                    do ix = 1, nx
                        cfx1(ix) = coff(ny*(ix-1) + iy)
                    enddo

                    cfx2 = (zero, zero)
                    call zgemm ('T', 'N', nx, 1, nx, aval, ztmtx, nx, &
                            cfx1, nx, bval, cfx2, nx )

                    do ix = 1, nx
                        cofp(ny*(ix-1) + iy) = cfx2(ix)
                    enddo

                enddo

                coff = (zero, zero); nxy = 0
                do ix = 1, nx

                    cfy1 = (zero, zero)
                    do iy = 1, ny
                        cfy1(iy) = cofp(nx*(ix-1) + iy)
                    enddo

                    cfy2 = (zero, zero)
                    call zgemm ('T', 'N', ny, 1, ny, aval, ztmty, ny, &
                            cfy1, ny, bval, cfy2, ny )

                    do iy = 1, ny
                        nxy = nxy + 1
                        coff(nxy) = cfy2(iy)
                    enddo

                enddo

                cofp = (zero, zero)
                do i2d = 1, n2d
                    cofp(i2d) = tdiag(i2d)*coff(i2d)
                enddo

                coff = (zero, zero); nxy = 0
                do ix = 1, nx

                    cfy1 = (zero, zero)
                    do iy = 1, ny
                        cfy1(iy) = cofp(nx*(ix-1) + iy)
                    enddo

                    cfy2 = (zero, zero)
                    call zgemm ('N', 'N', ny, 1, ny, aval, ztmty, ny, &
                            cfy1, ny, bval, cfy2, ny )

                    do iy = 1, ny
                        nxy = nxy + 1
                        coff(nxy) = cfy2(iy)
                    enddo

                enddo

                cofp = (zero, zero)
                do iy = 1, ny

                    cfx1 = (zero, zero)
                    do ix = 1, nx
                        cfx1(ix) = coff(ny*(ix-1) + iy)
                    enddo

                    cfx2 = (zero, zero)
                    call zgemm ('N', 'N', nx, 1, nx, aval, ztmtx, nx, &
                            cfx1, nx, bval, cfx2, nx )

                    do ix = 1, nx
                        cofp(ny*(ix-1) + iy) = cfx2(ix)
                    enddo

                enddo

                teng = (zero, zero); wfnrm = zero; indpl = (zero, zero)
                do i2d = 1, n2d
                    teng = teng + conjg(cofi(i2d))*cofp(i2d) + &
                            conjg(cofi(i2d))*(vx(i2d)-io*vcap(i2d))*cofi(i2d)
                    wfnrm = wfnrm  + cdabs(conjg(cofi(i2d))*cofi(i2d))
                    indpl = indpl + conjg(cofi(i2d))*(frce(i2d) + (epls*( &
                            dsin(omega*time) + dcos(omega*time))))*cofi(i2d)
                enddo

                write(11,*) time, epls*dcos(omega*time), epls*dsin(omega*time), &
                        real(teng), aimag(teng), wfnrm, real(indpl)

            enddo
            close(11)

            open(12,file=fltden)
            nxy = 0
            do ix = 1, nx
                x = xmin + (ix-one)*dx
                do iy = 1, ny
                    y = ymin + (iy-one)*dy
                    nxy = nxy + 1
                    write(12,'(*(f12.7))') x, y, &
                            (tmcof(nxy,itim), itim=1,ntim/nstep)
                enddo
                write(12,*)
            enddo
            close(12)

            open(12,file=fltpot)
            nxy = 0
            do ix = 1, nx
                x = xmin + (ix-one)*dx
                do iy = 1, ny
                    y = ymin + (iy-one)*dy
                    nxy = nxy + 1
                    write(12,'(*(f12.7))') x, y, &
                            (tmptn(nxy,itim), itim=1,ntim/nstep)
                enddo
                write(12,*)
            enddo
            close(12)

            deallocate(ztmtx, ztmty, cofi, coff, cofp, cfx1, cfx2, cfy1, cfy2)
            deallocate(tdiag, vx, vxexp, vcap, frce, tmcof, tmptn)
        end
end module time_prop_Trotter

module time_prop_RK4

end module time_prop_RK4
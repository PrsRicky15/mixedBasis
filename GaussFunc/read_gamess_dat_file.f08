! Created by prash on 6/6/2024.

module read_gamess_dat_file
    implicit none
    contains

    SUBROUTINE READ_FILE_(ifnm, chrg, coratm, corbas, natm_shell, shell_func, nprm_shell, coef, expnt)
        USE const
        IMPLICIT NONE

        character(len=50), intent(in) :: ifnm
        integer(sp), intent(out) :: natm_shell(nbsc), shellfunc(nshell)
        integer(sp), intent(out) :: nprm_shell(nshell)
        REAL(dp), intent(out) :: coef(nprm), expnt(nprm)
        REAL(dp), intent(out) :: coratm(natm,3), chrg(natm)
        REAL(dp), intent(out) :: corbas(nbsc,3)
        INTEGER(sp) :: npno, ipno, IPRM
        INTEGER(sp) :: iatm, ibsc, ival, ishell, ishel_bas
        REAL(dp)    :: ANO, xx, yy, zz
        CHARACTER(LEN=50) :: VAR

90    FORMAT(A8,F7.1,3F18.10)
91    FORMAT(I6,F26.10,F12.8)
94    FORMAT(I6,F26.10,F12.8,F12.8)
93    FORMAT(A5,I7)

        OPEN(21,FILE=ifnm)
        DO WHILE (trim(adjustl(VAR)) /= "C1" )
            READ(21,'(A2)') VAR
        ENDDO
        iprm = 0; iatm = 0; ibsc = 0
        200  READ(21,90) VAR, ANO, xx, yy, zz

        if (int(ano) /= 0 ) then
            ibsc = ibsc + 1
            corbas(ibsc,1) = xx*bohr
            corbas(ibsc,2) = yy*bohr
            corbas(ibsc,3) = zz*bohr
        endif

        if ( ano > 0 ) then
            iatm = iatm + 1
            coratm(iatm,1) = xx*bohr
            coratm(iatm,2) = yy*bohr
            coratm(iatm,3) = zz*bohr
            chrg(iatm) = ano
        endif

        IF (trim(adjustl(VAR)) == "$END") THEN
            GOTO 900
        ENDIF

        ishel_bas = 0
        DO WHILE (Len_trim(VAR) /= 0)
            READ(21,93) VAR , nPNO

            select case(trim(adjustl(var)))
            case("S")
                Ishell = Ishell + 1; ishel_bas = ishel_bas + 1
                shellfunc(ishell) = 1
                nprm_shell(ishell) = npno
                DO ipno = 1, npno
                    iprm = iprm + 1
                    READ(21,91) IVAL, expnt(iprm), coef(iprm)
                enddo
            case("P")
                ishel_bas = ishel_bas + 1
                ishell = ishell + 1; shellfunc(ishell) = 3
                nprm_shell(ishell) = npno
                DO ipno = 1, npno
                    iprm = iprm + 1
                    READ(21,91) IVAL, expnt(iprm), coef(iprm)
                ENDDO
            case("L")
                ishel_bas = ishel_bas + 2
                ishell = ishell + 1; shellfunc(ishell) = 1
                nprm_shell(ishell) = npno
                DO ipno = 1, npno
                    iprm = iprm + 1
                    READ(21,94) IVAL, expnt(iprm), coef(iprm), coef(iprm+npno)
                ENDDO
                ishell = ishell + 1; shellfunc(ishell) = 3
                nprm_shell(ishell) = npno
                ival = iprm
                DO ipno = 1, npno
                    iprm = iprm + 1
                    expnt(iprm) = expnt(ival-npno+ipno)
                ENDDO
            case("SP")
                ishel_bas = ishel_bas + 2
                ishell = ishell + 1; shellfunc(ishell) = 1
                nprm_shell(ishell) = npno
                DO ipno = 1, npno
                    iprm = iprm + 1
                    READ(21,94) IVAL, expnt(iprm), coef(iprm), coef(iprm+npno)
                ENDDO
                ishell = ishell + 1; shellfunc(ishell) = 3
                nprm_shell(ishell) = npno
                ival = iprm
                DO ipno = 1, npno
                    iprm = iprm + 1
                    expnt(iprm) = expnt(ival-npno+ipno)
                ENDDO
            case("D")
                ishel_bas = ishel_bas + 1
                ishell = ishell + 1; shellfunc(ishell) = 6
                nprm_shell(ishell) = npno
                DO ipno = 1, npno
                    iprm = iprm + 1
                    READ(21,91) IVAL, expnt(iprm), COEF(iprm)
                enddo
            case("F")
                ishel_bas = ishel_bas + 1
                ishell = ishell + 1; shellfunc(ishell) = 10
                nprm_shell(ishell) = npno
                DO ipno = 1, npno
                    iprm = iprm + 1
                    READ(21,91) IVAL, expnt(iprm), COEF(iprm)
                enddo
            case("G")
                ishel_bas = ishel_bas + 1
                ishell = ishell + 1; shellfunc(ishell) = 15
                nprm_shell(ishell) = npno
                DO ipno = 1, npno
                    iprm = iprm + 1
                    READ(21,91) IVAL, expnt(iprm), COEF(iprm)
                enddo
            end select
        enddo
        natm_shell(ibsc) = ishel_bas
        GOTO 200
        900  continue
        CLOSE(21)
    END

    SUBROUTINE gaussian_basis_count(ifnm, natom, nbsc, ninpshel, nshell, ngaus, nprm, nbasis)
        USE const
        IMPLICIT NONE
        CHARACTER(len=50), INTENT(in) :: ifnm
        INTEGER(sp), INTENT(out) :: natom, nbsc, ninpshel, nshell, nprm, nbasis
        INTEGER(sp), INTENT(out) :: ngaus, nprm

        INTEGER(sp) :: iatom, ibas, iinpshel, ishell, igaus, iprm, ibasis, ipno, pno
        REAL(dp)    :: xx, yy, zz, ANO

        CHARACTER(LEN=50) :: VAR, char

        90   FORMAT(A8,F7.1,3F18.10)
        93   FORMAT(A5,I7)

        OPEN(11, FILE=ifnm)
        DO WHILE ( TRIM(VAR) /= "C1" )
            READ(11, '(A2)') VAR
        END DO

        igaus = 0; iatom = 0; ibas = 0; ishell = 0
        iprm = 0; ibasis = 0; iinpshel = 0

        DO
            READ(11, 90) VAR, ANO, xx, yy, zz

            IF ( INT(ANO) /= 0 ) ibas = ibas + 1
            IF ( INT(ANO) > 0 ) iatom = iatom + 1

            IF (TRIM(ADJUSTL(VAR)) == "$END") THEN
                natom = iatom; nbsc = ibas; nshell = ishell; nprm = iprm
                ninpshel = iinpshel; ngaus = igaus; nbasis = ibasis
                EXIT
            ENDIF

            DO WHILE (LEN_TRIM(VAR) /= 0)
                READ(11, 93) VAR, PNO
                SELECT CASE(TRIM(ADJUSTL(VAR)))
                CASE("S", "P", "L", "SP", "D", "F", "G")
                    ishell = ishell + SELECTED_INT_KIND(PNO > 0) * (1 + SELECTED_INT_KIND(TRIM(ADJUSTL(VAR)) /= "S"))
                    iinpshel = iinpshel + 1
                    igaus = igaus + ABS(PNO)
                    iprm = iprm + SELECTED_INT_KIND(TRIM(ADJUSTL(VAR)) == "S") * PNO + &
                            SELECTED_INT_KIND(TRIM(ADJUSTL(VAR)) == "P") * (3 * PNO) + &
                            SELECTED_INT_KIND(TRIM(ADJUSTL(VAR)) == "L") * (4 * PNO) + &
                            SELECTED_INT_KIND(TRIM(ADJUSTL(VAR)) == "SP") * (4 * PNO) + &
                            SELECTED_INT_KIND(TRIM(ADJUSTL(VAR)) == "D") * (6 * PNO) + &
                            SELECTED_INT_KIND(TRIM(ADJUSTL(VAR)) == "F") * (10 * PNO) + &
                            SELECTED_INT_KIND(TRIM(ADJUSTL(VAR)) == "G") * (15 * PNO)
                    ibasis = ibasis + SELECTED_INT_KIND(TRIM(ADJUSTL(VAR)) == "S") * PNO + &
                            SELECTED_INT_KIND(TRIM(ADJUSTL(VAR)) == "P") * (3 * PNO) + &
                            SELECTED_INT_KIND(TRIM(ADJUSTL(VAR)) == "L") * (4 * PNO) + &
                            SELECTED_INT_KIND(TRIM(ADJUSTL(VAR)) == "SP") * (4 * PNO) + &
                            SELECTED_INT_KIND(TRIM(ADJUSTL(VAR)) == "D") * (6 * PNO) + &
                            SELECTED_INT_KIND(TRIM(ADJUSTL(VAR)) == "F") * (10 * PNO) + &
                            SELECTED_INT_KIND(TRIM(ADJUSTL(VAR)) == "G") * (15 * PNO)
                END SELECT
            END DO
        END DO
        CLOSE(11)
    END SUBROUTINE gaussian_basis_count

end module read_gamess_dat_file
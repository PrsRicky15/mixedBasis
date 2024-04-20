module iface_help_2
    implicit none

    interface sub
        module subroutine sub1(x)
            integer, intent(inout) :: x
        end subroutine
        module subroutine sub2(x)
            real, intent(inout) :: x
        end subroutine
    end interface
end module

submodule (iface_help_2) iface_help_2_impl
    implicit none
contains
    module procedure sub1
        x = x+1
    end procedure
    module procedure sub2
        x = x+1
    end procedure
end submodule 

program main
    use iface_help_2
    real :: rx
    integer :: ix
    ix = 1
    rx = 1
    call sub(ix)
    call sub(rx)
    print *, ix, rx
end program

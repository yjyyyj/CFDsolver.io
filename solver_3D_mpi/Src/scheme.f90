module scheme_mod
    !******************************************************************
    !*     The scheme has state for Solver                            *
    !******************************************************************
    use param
    implicit none
    !
    type scheme
        procedure(),nopass,pointer :: calc_step => NULL()
        procedure(),nopass,pointer :: calc_faceQ => NULL()
        procedure(),nopass,pointer :: calc_flux => NULL()
        procedure(),nopass,pointer :: calc_visflux => NULL()
        contains
        procedure :: select_flux
        procedure :: select_visflux
        procedure :: select_muscl
        procedure :: select_step
        ! procedure :: set_q
    end type scheme

    ! interface scheme
    !     module procedure init_scheme
    ! end interface scheme

contains
    ! fluxの実装部
    ! type(scheme) function init_scheme(self, q, q0, qc)
    !     class(scheme) :: self
    !     double precision :: q, q0, qc
            
    ! end function init_scheme

    ! subroutine set_q(self,q)
    !     class(scheme) self
    !     double precision,dimension(ndmax,0:(jmax+1),0:(kmax+1)):: q
    ! end subroutine set_q

    subroutine select_flux(self,irhs)
        use flux
        class(scheme) self
        integer,intent(in) :: irhs

        select case(irhs)
        case(1)
            write(*,*) "DIVERGENCE"
            self%calc_flux => flux_div
        case(2)
            write(*,*) "UPWIND"
            self%calc_flux => flux_upwind
        case(3)
            write(*,*) "SLAU"
            self%calc_flux => flux_slau
        case(4)
            write(*,*) "KEEP_PE"
            self%calc_flux => flux_KEEP_PE
        case(5)
            write(*,*) "Proposed"
            ! self%calc_flux => flux_proposed
            self%calc_flux => flux_prodiv
        case default
            write(*,*) "[Error] invalid number of flux in main.f90: ",irhs
            stop
        end select
    end subroutine select_flux

    subroutine select_visflux(self,vflag)
        class(scheme) self
        integer,intent(in) :: vflag

        select case(vflag)
        case(0)
            write(*,*) "Euler"
            ! self%calc_visflux => NULL()
        case(1)
            write(*,*) "LAD"
            self%calc_visflux => visflux_LAD
        case(2)
            write(*,*) "DNS"
            self%calc_visflux => visflux_dns
        case default
            write(*,*) "[Error] invalid number of visflux in main.f90: ",vflag
            stop
        end select
    end subroutine select_visflux

    subroutine select_muscl(self,acc)
        class(scheme) self
        integer,intent(in) :: acc
        external firstOrderQ, muscl_va

        select case(acc)
        case(0)
            write(*,*) "1stOrder"
            self%calc_faceQ => firstOrderQ
        case(1)
            write(*,*) "MUSCL(2nd)"
            self%calc_faceQ => muscl_va
        case default
            write(*,*) "[Error] invalid number of muscl in main.f90: ",acc
            stop
        end select
    end subroutine select_muscl
    
    subroutine select_step(self,ilhs)
        class(scheme) self
        integer,intent(in) :: ilhs
        external step_euler, step_RK4

        select case(ilhs)
        case(1)
            write(*,*) "1stEuler"
            self%calc_step => step_euler
        case(2)
            write(*,*) "TVDRK3"
            self%calc_step => step_RK4
        case(3)
            write(*,*) "RK4"
            self%calc_step => step_RK4
        case default
            write(*,*) "[Error] invalid number of LHS in main.f90: ",ilhs
            stop
        end select
    end subroutine select_step
end module scheme_mod
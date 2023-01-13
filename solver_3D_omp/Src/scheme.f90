module scheme_mod
    !******************************************************************
    !*     The scheme with important state of Solver                  *
    !******************************************************************
    use param
    implicit none
    !
    type scheme
        procedure(),nopass,pointer :: calc_step => NULL()
        procedure(),nopass,pointer :: calc_faceQ => NULL()
        procedure(),nopass,pointer :: calc_flux => NULL()
        procedure(),nopass,pointer :: calc_visflux => NULL()
        procedure(),nopass,pointer :: calc_outf => NULL()
        procedure(),nopass,pointer :: calc_outf_init => NULL()
        contains
        procedure :: select_flux
        procedure :: select_visflux
        procedure :: select_muscl
        procedure :: select_step
        procedure :: select_outf
        ! procedure :: set_q
    end type scheme

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
            self%calc_flux => flux_proposed
            ! self%calc_flux => flux_prodiv
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
            self%calc_visflux => visflux_none
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
    
    subroutine select_outf(self,dim_outf)
        class(scheme) self
        integer,intent(in) :: dim_outf
        external outf_1d, outf_2d, outf_3d, outf_1d_init, outf_2d_init, outf_3d_init

        select case(dim_outf)
        case(1)
            write(*,*) "1D outflow"
            self%calc_outf => outf_1d
            self%calc_outf_init => outf_1d_init
        case(2)
            write(*,*) "2D outflow"
            self%calc_outf => outf_2d
            self%calc_outf_init => outf_2d_init
        case(3)
            write(*,*) "3D outflow"
            self%calc_outf => outf_3d
            self%calc_outf_init => outf_3d_init
        case default
            write(*,*) "[Error] invalid number of outflow dimension in main.f90: ",dim_outf
            stop
        end select
    end subroutine select_outf

    subroutine select_step(self,ilhs)
        class(scheme) self
        integer,intent(in) :: ilhs
        external step_euler, step_RK4

        select case(ilhs)
        case(1)
            write(*,*) "1stEuler"
            self%calc_step => step_euler
        case(2)
            write(*,*) "RK4"
            self%calc_step => step_RK4
        case default
            write(*,*) "[Error] invalid number of LHS in main.f90: ",ilhs
            stop
        end select
    end subroutine select_step

end module scheme_mod
program main
	!******************************************************************
  !*     3D multi-species Euler Solver                              *
  !******************************************************************
  use param_mod
  use scheme_mod

  implicit none
  integer n
  double precision,allocatable,dimension(:,:,:,:):: q, q0, qc
  type(scheme) :: myscheme

  ! set main param
  call read_stdin()
  allocate( q(ndmax,0:(jmax+1),0:(kmax+1),0:(lmax+1)) )
  allocate( q0(ndmax,jmax,kmax,lmax) )
  allocate( qc(ndmax,jmax,kmax,lmax) )

  !***** prepare calclation ******************************************************!

  !*** select scheme subroutine->  myscheme%calc_"hoge"
  call myscheme%select_step(ilhs)
  call myscheme%select_faceQ(faceAcc)
  call myscheme%select_flux(irhs)
  call myscheme%select_visflux(vflag)
  call myscheme%select_outf(dim_outf)

  call init_flow(q, q0, qc, myscheme)

  !***** main loop *****************************************************!
  write(*,*) "********* start main loop *********"

  do n=1,nmax
    call myscheme%calc_step(q, qc, myscheme)
    if(mod(n,nout)==0) then
      call sumqc(qc, n)
      call residual(n)
      call myscheme%calc_outf(q,n)
    endif
  enddo

  !********************************************************************!

  stop
end program main

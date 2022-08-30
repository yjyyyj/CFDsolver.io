program main
	!******************************************************************
  !*     3D multi-species Euler Solver                              *
  !******************************************************************
  use param
  use scheme_mod

  implicit none
  integer n,nmax,i
  double precision,allocatable,dimension(:,:,:,:):: q, q0, qc
  type(scheme) :: myscheme

  ! read param 
  write(*,*) "reading stdin"
  open(5, file='stdin') 
  read(5,*) nmax
  read(5,*) nout
  read(5,*) dt
  read(5,*) Mach
  read(5,*) Pr
  read(5,*) Re
  read(5,*) nspecies
  read(5,*) irhs
  read(5,*) acc
  read(5,*) ilhs
  read(5,*) vflag
  read(5,*) mconst
  read(5,*) bc1, bc2, bc3, bc4, bc5, bc6 
  close(5)

  write(*,*) nmax,nout,dt,Mach,Pr,Re,nspecies,irhs,acc,ilhs,vflag,bc1,bc2,bc3,bc4,bc5,bc6

  ! set param
  jmax  = 121
  kmax  = 3
  lmax  = 201
  ! jmax  = 101
  ! kmax  = 101
  ! lmax  = 101

  ndim  = 3
  ndmax = 5                ! r,ru,rv,rw,E
  ndmax = ndmax + nspecies ! + ry*N (including "r") 

  allocate( xg(jmax) )
  allocate( yg(kmax) )
  allocate( zg(lmax) )
  allocate( dx(ndim) )
  allocate( sums0(ndmax) )
  allocate( resid(ndmax) )
  allocate( qinf(ndmax) )
  allocate( q(ndmax,0:(jmax+1),0:(kmax+1),0:(lmax+1)) )
  allocate( q0(ndmax,jmax,kmax,lmax) )
  allocate( qc(ndmax,jmax,kmax,lmax) )
  allocate( gam(0:(jmax+1),0:(kmax+1),0:(lmax+1)), mw(0:(jmax+1),0:(kmax+1),0:(lmax+1)), gami(nspecies), mwi(nspecies) )

  pi = 4.0d0*atan(1.0d0)
  Ru = 8.3d0

  gami(1)  = 1.4d0        ! specific heat ratio # N2
  gami(2)  = 1.09d0        ! specific heat ratio # SF6
  ! gami(1)  = 1.4d0        ! specific heat ratio # N2
  ! gami(2)  = 1.66d0        ! specific heat ratio # He
  ! gami(3)  = 1.29d0        ! specific heat ratio # CO2
  ! gami(4)  = 1.31d0        ! specific heat ratio # CH4

  mwi(1)  = 28.8d0        ! specific heat ratio
  mwi(2)  = 146.0d0      ! specific heat ratio
  ! mwi(1)  = 28.d0        ! specific heat ratio
  ! mwi(2)  = 4.0d0        ! specific heat ratio
  ! mwi(3)  = 44.d0        ! specific heat ratio
  ! mwi(4)  = 16.d0        ! specific heat ratio

  if(mconst == 1) then 
    mwi(:)  = 1.d0        ! specific heat ratio
  endif

  do i=1,nspecies
    gami(i)  = 1d0/(gami(i) -1d0)        ! inv specific heat ratio
  enddo

  !***** main *********************************************************!

  call init(q, q0, qc)

  call myscheme%select_step(ilhs)
  call myscheme%select_muscl(acc)
  call myscheme%select_flux(irhs)

  do n=1,nmax
    call myscheme%calc_step(q, qc, myscheme)
    if(mod(n,nout)==0) then
      ! call sumdf(qc, n)
      call residual(n)
      ! call outf(q, n)
      call outf_slice(q, n)
    endif
  enddo

  !********************************************************************!

  stop
end program main

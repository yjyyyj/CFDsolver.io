subroutine read_stdin()
  !**********************************************************************
  !*     initialize                                                     *
  !**********************************************************************
  use param
  implicit none
  integer i

  !*** set param from stdin ********************

  ! read param 
  write(*,*) "reading stdin"
  open(5, file='stdin') 
    read(5,*) nmax
    read(5,*) nout
    read(5,*) restart
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

  write(*,*) nmax,nout,restart,dt,Mach,Pr,Re,nspecies,irhs,acc,ilhs,vflag,bc1,bc2,bc3,bc4,bc5,bc6

  ndim  = 3
  ndmax = 5                ! r,ru,rv,rw,E
  ndmax = ndmax + nspecies ! + N-ry (including "r" <- not solve) 

  ! allocate param 
  allocate( xg(jmax) )
  allocate( yg(kmax) )
  allocate( zg(lmax) )
  allocate( dx(ndim) )
  allocate( sums0(ndmax) )
  allocate( resid(ndmax) )
  allocate( qinf(ndmax) )
  allocate( gam(0:(jmax+1),0:(kmax+1),0:(lmax+1)), mw(0:(jmax+1),0:(kmax+1),0:(lmax+1)), gami(nspecies), mwi(nspecies) )

  ! set const param
  pi = 4.0d0*atan(1.0d0)
  Ru = 8.3d0 ! 

  ! gami(1)  = 1.4d0        ! specific heat ratio # N2
  ! gami(2)  = 1.09d0        ! specific heat ratio # SF6
  gami(1)  = 1.4d0        ! specific heat ratio # N2
  gami(2)  = 1.66d0        ! specific heat ratio # He
  ! gami(3)  = 1.29d0        ! specific heat ratio # CO2
  ! gami(4)  = 1.31d0        ! specific heat ratio # CH4

  ! mwi(1)  = 28.8d0        ! specific heat ratio
  ! mwi(2)  = 146.0d0      ! specific heat ratio
  mwi(1)  = 28.d0        ! specific heat ratio
  mwi(2)  = 4.0d0        ! specific heat ratio
  ! mwi(3)  = 44.d0        ! specific heat ratio
  ! mwi(4)  = 16.d0        ! specific heat ratio

  if(mconst == 1) then 
    mwi(:)  = 1.d0        ! specific heat ratio
  endif

  do i=1,nspecies
    gami(i)  = 1d0/(gami(i) -1d0)        ! inv specific heat ratio
  enddo

  write(*,*) "read stdin"

  return
end subroutine read_stdin


subroutine read_flow(q0)
  !**********************************************************************
  !*     initialize                                                     *
  !**********************************************************************
  use param
  implicit none
  double precision x,z,g
  integer j,k,l
  double precision, dimension(ndmax,jmax,kmax,lmax) :: q0


  q0 = 0.0d0

  !*** read restart flow ********************
  open(5, file='flow.dat')

    do l = 1, lmax
      do j = 1, jmax
        read (5, *) x, z, q0(:,j,2,l), g
        ! write(*, fmt='(11E25.15e3)',advance='No') x, z, q0(:,j,2,l), g
        ! write(*,*)
      end do
    end do

  close(5) 

  do k=1,kmax
    q0(:,:,k,:) = q0(:,:,2,:)
  enddo
  
  write(*,*) "set restart"
  
  return
end subroutine read_flow

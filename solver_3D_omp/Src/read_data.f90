subroutine read_stdin()
  !**********************************************************************
  !*     read parameters from stdin                                     *
  !**********************************************************************
  use param_mod
  implicit none
  integer i

  !*** set param from stdin ********************

  ! read param 
  write(*,*) "****** start reading stdin ******"
  open(5, file='stdin') 
    read(5,*) nmax
    read(5,*) nout
    read(5,*) restart
    read(5,*) jmax, kmax, lmax
    read(5,*) dt
    read(5,*) Mach
    read(5,*) Pr
    read(5,*) Re
    read(5,*) nspecies
    read(5,*) irhs
    read(5,*) faceAcc
    read(5,*) ilhs
    read(5,*) vflag
    read(5,*) dim_outf
    read(5,*) mconst
    read(5,*) bc1, bc2, bc3, bc4, bc5, bc6 
  close(5)

  write(*,*) nmax,nout,restart,dt
  write(*,*) jmax, kmax, lmax
  write(*,*) Mach,Pr,Re
  write(*,*) nspecies,irhs,faceAcc,ilhs,vflag,dim_outf
  write(*,*) bc1,bc2,bc3,bc4,bc5,bc6

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
  allocate( gam(0:(jmax+1),0:(kmax+1),0:(lmax+1)), mw(0:(jmax+1),0:(kmax+1),0:(lmax+1)))

  ! set const param
  pi = 4.0d0*atan(1.0d0)
  Ru = 8.3d0 ! 

  gami(1)  = 1.4d0         ! specific heat ratio # N2
  gami(2)  = 1.66d0        ! specific heat ratio # He
  gami(3)  = 1.29d0        ! specific heat ratio # CO2
  gami(4)  = 1.31d0        ! specific heat ratio # CH4
  gami(5)  = 1.09d0        ! specific heat ratio # SF6

  mwi(1)  = 28.d0        ! specific heat ratio
  mwi(2)  = 4.0d0        ! specific heat ratio
  mwi(3)  = 44.d0        ! specific heat ratio
  mwi(4)  = 16.d0        ! specific heat ratio
  mwi(5)  = 146.0d0      ! specific heat ratio

  if(mconst == 1) then 
    mwi(:)  = 1.d0        ! specific heat ratio
  endif

  do i=1,nspecies
    gami(i)  = 1d0/(gami(i) -1d0)        ! inv specific heat ratio
  enddo

  write(*,*) "****** read stdin fin *******"

  return
end subroutine read_stdin


subroutine read_flow(q0)
  !**********************************************************************
  !*     read flow data from restart file (flow.dat)                    *
  !**********************************************************************
  use param_mod
  implicit none
  double precision x,z,g
  integer j,k,l
  double precision, dimension(ndmax,jmax,kmax,lmax) :: q0


  q0 = 0.0d0

  !*** read restart flow ********************
  open(5, file='flow.dat')

    do l = 1, lmax
      do j = 1, jmax
        !** please edit this line for the format of restart data 
        read (5, *) x, z, q0(:,j,2,l), g
      end do
    end do

  close(5) 

  do k=1,kmax
    q0(:,:,k,:) = q0(:,:,2,:)
  enddo
  
  write(*,*) "set restart"
  
  return
end subroutine read_flow

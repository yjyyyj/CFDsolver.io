subroutine init(q, q0, qc)
  !**********************************************************************
  !*     initialize                                                     *
  !**********************************************************************
  use param
  implicit none
  double precision x0,x1,y0,y1,z0,z1
  double precision uvw2,eps
  integer j,k,l
  double precision, dimension(ndmax,0:(jmax+1),0:(kmax+1),0:(lmax+1)) :: q
  double precision, dimension(ndmax,jmax,kmax,lmax) :: q0, qc
  
  x0 = 0.0d0
  ! x1 = 1.0d0
  ! x1 = 4.0d0
  x1 = 12.0d0

  y0 = 0.0d0
  y1 = 1.0d0

  z0 = 0.0d0
  ! z1 = 1.0d0
  z1 = 5d0

  q = 0.0d0
  q0 = 0.0d0
  qc = 0.0d0
  sums0 = 0.0d0

  dx(1) = (x1 - x0)/ (jmax-1)
  dx(2) = (y1 - y0)/ (kmax-1)
  dx(3) = (z1 - z0)/ (lmax-1) 

  write(*,*) "dx  : ", dx(1)
  write(*,*) "dy  : ", dx(2)
  write(*,*) "dz  : ", dx(3)
  write(*,*) "dt  : ", dt
  write(*,*) "cfl : ", dt/dx(1)

  do j=1,jmax
    xg(j) = dble(j-1)*dx(1)
  end do
  do k=1,kmax
    yg(k) = dble(k-1)*dx(2)
  end do
  do l=1,lmax
    zg(l) = dble(l-1)*dx(3) + 0.5d0*dx(3)
  end do

  !*** set inital flow ********************

  ! call contact_init(q0)
  ! call contact_init_sharp(q0)
  ! call prin_init(q0)
  ! call multiprin_init(q0)
  ! call flower_init(q0)

  ! call sod_init(q0)
  ! call TGV_init(q0)
  ! call RMinstability_init(q0)
  call Blasius_init(q0)
  
  ! call nonDimtize(q0)

  q(:,1:jmax,1:kmax,1:lmax) = q0(:,1:jmax,1:kmax,1:lmax)            ! set initial q

  do l=1,lmax
    do k=1,kmax
      do j=1,jmax
        ! set Conservative Q
        mw(j,k,l)  = calmw(q(:,j,k,l))
        gam(j,k,l) = calgm(q(:,j,k,l),mw(j,k,l))
        uvw2 = q(2,j,k,l)*q(2,j,k,l) +q(3,j,k,l)*q(3,j,k,l) +q(4,j,k,l)*q(4,j,k,l)

        qc(1,j,k,l) = q(1,j,k,l)             ! rho
        qc(2,j,k,l) = q(1,j,k,l)*q(2,j,k,l)  ! rho*u
        qc(3,j,k,l) = q(1,j,k,l)*q(3,j,k,l)  ! rho*v
        qc(4,j,k,l) = q(1,j,k,l)*q(4,j,k,l)  ! rho*w
        qc(5,j,k,l) = q(5,j,k,l)*gam(j,k,l) + 0.5d0*q(1,j,k,l)*uvw2 ! E = p/(gam-1)+0.5*rho*(u^2+v^2+w^2)
        qc(6:ndmax,j,k,l) = q(6:ndmax,j,k,l) ! rho_yi
      enddo
    enddo
  enddo

  do l=1,lmax-1
    do k=1,kmax-1
      do j=1,jmax-1
        sums0(1:ndmax) = sums0(1:ndmax) + qc(1:ndmax,j,k,l)
      enddo
    enddo
  enddo

  ! zero/
  eps = 1e-10
  sums0(:) = (sums0(:)+eps)/(jmax-1)/(kmax-1)/(lmax-1) 

  ! open(150,file="sums.dat",form="formatted")
  !   write(150,fmt='(I8, 9E20.10e3)',advance='No') 0, sums0(:)-sums0(:)
  ! close(150)
  open(130,file="residual.dat",form="formatted")
    write(130,fmt='(I8, 9E20.10e3)',advance='No') 0, sums0(:)-sums0(:)
  close(130)

  !!!**** out 3d full data ********
  ! open(500,file="output_000000.dat",form="formatted")
  !   do l=1,lmax
  !     do k=1,kmax
  !       do j=1,jmax
  !         ! temp = mw0/q0(1,j,k) * q0(4,j,k)/Ru          ! T = (M/r)*(p/R)
  !         write(500,fmt='(11E25.15e3)',advance='No') xg(j),yg(k),zg(l),q0(1:ndmax,j,k,l),gam(j,k,l)
  !         write(500,*)
  !       end do
  !       write(500,*)
  !     end do
  !     write(500,*)
  !   end do
  ! close(500)

  !**** out 1d slice data ********
  open(100,file='output_000000.dat',form="formatted")
    do l=1,lmax
      do k=1,kmax
        do j=1,jmax
          if( yg(k)==0.5d0 ) then
            write(100,fmt='(11E25.15e3)',advance='No') xg(j),zg(l),q(1:ndmax,j,k,l),gam(j,k,l)
            write(100,*)
          endif
        end do
      end do
    end do
  close(100)
  
  write(*,*) "set init"
  return
end subroutine init

subroutine nonDimtize(q0)
  use param
  implicit none
  double precision, dimension(ndmax,jmax,kmax,kmax) :: q0
  double precision  gamma
  double precision  rinf,ainf,dref,pinf,uinf,vinf,winf,einf,uvw2

  gamma = gami(1)

  ! set infty parameter ( refLength = 1 )
  rinf   = 1.0d0
  ainf   = 1.0d0
  dref   = 1.0d0
  pinf   = rinf*ainf*ainf/gamma
  uinf   = Mach*ainf
  vinf   = 0d0
  winf   = 0d0
  uvw2   = uinf**2+vinf**2+winf**2
  einf   = pinf/(gamma-1)+0.5d0*rinf*uvw2

  ! compute non-dimentional parameters 
  q0(1,:,:,:) = q0(1,:,:,:)/rinf     ! rho
  q0(2,:,:,:) = q0(2,:,:,:)/ainf     ! u
  q0(3,:,:,:) = q0(3,:,:,:)/ainf     ! v
  q0(4,:,:,:) = q0(4,:,:,:)/ainf     ! v
  q0(5,:,:,:) = q0(5,:,:,:)/(rinf*ainf*ainf)    ! p
  q0(6:ndmax,:,:,:) = q0(6:ndmax,:,:,:)/rinf    ! rhoy

  return
end subroutine nonDimtize

subroutine TGV_init(q0)
  use param
  implicit none
  integer j,k,l
  double precision, dimension(ndmax,jmax,kmax,lmax) :: q0
  double precision, dimension(nspecies)   :: ry, ryw
  double precision  rinf,pinf

  Mach = 0.05d0
  ryw(1) = 1.0d0
  ryw(2) = 0.0d0

  rinf = 1d0

  do l=1,lmax
    do k=1,kmax
      do j=1,jmax

        pinf = 1d0/(1d0/gam(j,k,l)+1d0)
        ry(1) = ryw(1)*rinf
        ry(2) = ryw(2)*rinf

        q0(1,j,k,l) = ry(1) + ry(2)     ! rho
        q0(2,j,k,l) = Mach*dsin(xg(j))*dcos(yg(k))*dcos(zg(l))   ! u
        q0(3,j,k,l) = -Mach*dcos(xg(j))*dsin(yg(k))*dcos(zg(l))   ! v
        q0(4,j,k,l) = 0d0   ! w
        q0(5,j,k,l) = pinf+(rinf*Mach*Mach*(dcos(2.0d0*xg(j))+dcos(2.0d0*yg(k)))*(dcos(2.0d0*zg(l))+2d0))/16d0         ! p
        q0(6:ndmax,j,k,l) = ry(:)   ! rhoy

      end do
    end do
  end do

  return
end subroutine TGV_init

subroutine contact_init(q0)
  use param
  implicit none
  integer j,k,l
  double precision, dimension(ndmax,jmax,kmax,lmax) :: q0
  double precision, dimension(nspecies)   :: ry, ryw
  double precision  r
  double precision  xc, yc, zc, const

  r = 0d0
  xc = 0.5d0
  yc = 0.5d0
  zc = 0.5d0

  const = 15d0
  ryw(1) = 0.8d0
  ryw(2) = 0.2d0

  ! contact 
  do l=1,lmax
    do k=1,kmax
      do j=1,jmax
        r = sqrt( (xg(j)- xc)**2 )
        ry(1) = ryw(1)*(1d0 - 0.5d0*(tanh(const*( r -0.25d0)) + 1d0))
        ry(2) = ryw(2)*(0.5d0*(tanh(const*( r -0.25d0)) + 1d0))

        q0(1,j,k,l) = ry(1) + ry(2)     ! rho
        q0(2,j,k,l) = 1.0d0         ! u
        q0(3,j,k,l) = 0.0d0         ! v
        q0(4,j,k,l) = 0.0d0         ! w
        q0(5,j,k,l) = 0.9d0         ! p
        q0(6:ndmax,j,k,l) = ry(:)   ! rhoy
      end do
    end do
  end do

  return
end subroutine contact_init

subroutine contact_init_sharp(q0)
  use param
  implicit none
  integer j,k,l
  double precision, dimension(ndmax,jmax,kmax,lmax) :: q0
  double precision, dimension(nspecies)   :: ry, ryw
  double precision  r
  double precision  xc, yc, zc, rc, const

  xc = 0.5d0
  yc = 0.5d0
  zc = 0.5d0
  rc = 0.25d0

  const = 5000d0
  ! const = 50d0
  ryw(1) = 0.5d0
  ryw(2) = 0.2d0
  ! ! *** T=const  *******
  ! ryw(1) = 0.005d0*mwi(1)
  ! ryw(2) = 0.005d0*mwi(2)

  ! contact 
  do l=1,lmax
    do k=1,kmax
      do j=1,jmax
        r = sqrt( (xg(j)- xc)**2 )
        ry(1) = ryw(1)*(1d0 - 0.5d0*(tanh(const*( r - rc)) + 1d0))
        ry(2) = ryw(2)*(0.5d0*(tanh(const*( r - rc)) + 1d0))

        q0(1,j,k,l) = ry(1) + ry(2)     ! rho
        q0(2,j,k,l) = 0.0d0         ! u
        q0(3,j,k,l) = 0.0d0         ! v
        q0(4,j,k,l) = 0.0d0         ! w
        q0(5,j,k,l) = 0.9d0         ! p
        q0(6:ndmax,j,k,l) = ry(:)   ! rhoy
      end do
    end do
  end do

  return
end subroutine contact_init_sharp

subroutine sod_init(q0)
  use param
  implicit none
  integer j,k,l
  double precision, dimension(ndmax,jmax,kmax,lmax) :: q0
  double precision, dimension(nspecies)   :: ry, ryw

  ryw(1) = 1.0d0
  ryw(2) = 0.0d0

  ! contact 
  do l=1,lmax
    do k=1,kmax
      do j=1,jmax
        if(xg(j)<0.5d0) then
          ry(1) = ryw(1)*1.0d0
          ry(2) = ryw(2)*1.0d0

          q0(1,j,k,l) = ry(1) + ry(2)     ! rho
          q0(2,j,k,l) = 0.0d0         ! u
          q0(3,j,k,l) = 0.0d0         ! v
          q0(4,j,k,l) = 0.0d0         ! w
          q0(5,j,k,l) = 1.0d0         ! p
          q0(6:ndmax,j,k,l) = ry(:)   ! rhoy
        else
          ry(1) = ryw(1)*0.125d0
          ry(2) = ryw(2)*0.125d0

          q0(1,j,k,l) = ry(1) + ry(2)     ! rho
          q0(2,j,k,l) = 0.0d0         ! u
          q0(3,j,k,l) = 0.0d0         ! v
          q0(4,j,k,l) = 0.0d0         ! w
          q0(5,j,k,l) = 0.1d0         ! p
          q0(6:ndmax,j,k,l) = ry(:)   ! rhoy
        endif
      end do
    end do
  end do

  return
end subroutine sod_init

subroutine prin_init(q0)
  use param
  implicit none
  integer j,k,l
  double precision, dimension(ndmax,jmax,kmax,lmax) :: q0
  double precision, dimension(nspecies)   :: ry, ryw
  double precision  r
  double precision  xc, yc, zc, const

  r = 0d0
  xc = 0.5d0
  yc = 0.5d0
  zc = 0.5d0

  const = 15d0
  ryw(1) = 0.8d0
  ryw(2) = 0.2d0

  ! contact 
  do l=1,lmax
    do k=1,kmax
      do j=1,jmax

        r = sqrt( (xg(j)- xc)**2 + (yg(k)- yc)**2 + (zg(l)- zc)**2 )
        ry(1) = ryw(1)*(1d0 - 0.5d0*(tanh(const*( r -0.25d0)) + 1d0))
        ry(2) = ryw(2)*(0.5d0*(tanh(const*( r -0.25d0)) + 1d0))
        ! ry(1) = ryw(1)*(1d0 - 0.5d0*(tanh(const*( r -0.25d0)) + 1d0)) + 0.2d0
        ! ry(2) = 0d0

        q0(1,j,k,l) = ry(1) + ry(2)     ! rho
        q0(2,j,k,l) = 1.0d0         ! u
        q0(3,j,k,l) = 1.0d0         ! v
        q0(4,j,k,l) = 1.0d0         ! v
        q0(5,j,k,l) = 1.0d0         ! p
        q0(6:ndmax,j,k,l) = ry(:)   ! rhoy
      end do
    end do
  end do

  return
end subroutine prin_init

subroutine multiprin_init(q0)
  use param
  implicit none
  integer j,k
  double precision, dimension(ndmax,jmax,kmax) :: q0
  double precision, dimension(nspecies)   :: ry, ryw, xc, yc
  double precision  r
  double precision  const

  r = 0d0
  xc(1) = 0.5d0
  yc(1) = 0.5d0
  xc(2) = 0.0d0
  yc(2) = 0.5d0
  xc(3) = 0.5d0
  yc(3) = 0.0d0

  const = 15d0
  ryw(1) = 0.4d0
  ryw(2) = 0.2d0
  ryw(3) = 0.1d0

  ! contact 
  do k=1,kmax
    do j=1,jmax

      r = sqrt( (xg(j)- xc(1))**2 + (yg(k)- yc(1))**2 )
      ry(1) = ryw(1)*(1d0 - 0.5d0*(tanh(const*( r -0.3d0)) + 1d0))
      ry(2) = ryw(2)*(0.5d0*(tanh(const*( r -0.3d0)) + 1d0))
      ry(3) = ryw(3)*(0.5d0*(tanh(const*( r -0.2d0)) + 1d0))

      ! if( r > 0.25d0) then
      !   ry(3) = ryw(3)*(1d0 - 0.5d0*(tanh(const*( r -0.4d0)) + 1d0))
      ! else 
      !   ry(3) = ryw(3)*(0.5d0*(tanh(const*( r -0.2d0)) + 1d0))
      ! endif

      q0(1,j,k) = sum(ry(:))     ! rho
      q0(2,j,k) = 1.0d0         ! u
      q0(3,j,k) = 1.0d0         ! v
      q0(4,j,k) = 1.0d0         ! p
      q0(5:ndmax,j,k) = ry(:)   ! rhoy
    end do
  end do

  return
end subroutine multiprin_init

subroutine RMinstability_init(q0)
  use param
  implicit none
  integer j,k
  double precision, dimension(ndmax,jmax,kmax) :: q0
  double precision, dimension(nspecies)   :: y, ryw
  double precision  r
  double precision  const
  double precision  xint
  double precision  gamma
  double precision  r1,r2,u1,u2,p1,p2


  !% Richtmyer–Meshkov instability
  ! air:1, SF6:2

  mach = 1.24d0
  gamma = 1d0/gami(1)+1d0 ! gamma(air)

  ! shock condition for air
  r1 = 1.0d0 
  u1 = 0.25d0
  p1 = 1.0d0/gamma

  p2 = p1*(1.0d0+(2d0*gamma/(gamma+1d0))*(mach*mach-1d0))
  u2 = u1*((gamma+1d0)+(gamma-1d0)*(p2/p1))/((gamma-1d0)+(gamma+1d0)*(p2/p1))         ! u
  r2 = r1*((gamma-1d0)+(gamma+1d0)*(p2/p1))/((gamma+1d0)+(gamma-1d0)*(p2/p1))     ! rho

  r = 0d0
  const = 10d0

  do k=1,kmax
    do j=1,jmax
      ! shock front (M=1.24)
      ryw(1) = r1
      ryw(2) = mwi(2)/mwi(1)

      r = yg(k) - 0.5d0
      xint = 2.5d0 + 0.1d0*sin(2d0*pi*( r +0.25d0 ))
      y(2) = 1d0 - 0.5d0*((tanh(const*( xg(j) -xint)/dx(1)))+1d0)
      y(1) = 1d0 - y(2)


      q0(1,j,k) = ryw(1)*y(1) + ryw(2)*y(2)     ! rho
      q0(2,j,k) = u1         ! u
      q0(3,j,k) = 0.0d0          ! v
      q0(4,j,k) = p1    ! p
      q0(5:ndmax,j,k) = ryw(:)*y(:)   ! rhoy

      ! shock though 
      if(xg(j)>3.6d0)then 
        ryw(1) = r2
        q0(1,j,k) = ryw(1)*y(1) + ryw(2)*y(2)  ! rho
        q0(2,j,k) = u2         ! u
        q0(3,j,k) = 0.0d0      ! v
        q0(4,j,k) = p2         ! p
        q0(5:ndmax,j,k) = ryw(:)*y(:)   ! rhoy
      endif
    end do
  end do

  return
end subroutine RMinstability_init

subroutine RTI_init(q0)
  use param
  implicit none
  integer j,k
  double precision, dimension(ndmax,jmax,kmax) :: q0
  double precision, dimension(nspecies)   :: y, ryw
  double precision  r
  double precision  const
  double precision  xint
  double precision  r1,r2,u1,u2,p1,p2


  !% Rayleigh–Taylor instability
  ! tritium:1, Sn:2

  r1 = 0.520d0 
  u1 = 0.0d0
  p1 = r1*8.30d0/mwi(1)*11690

  r2 = 7.330d0 
  u2 = 0.0d0
  p2 = r2*8.30d0/mwi(2)*11690

  r = 0d0
  const = 10d0

  do k=1,kmax
    do j=1,jmax
      ! shock front (M=1.24)
      ryw(1) = r1
      ryw(2) = mwi(2)/mwi(1)

      r = yg(k) - 0.5d0
      xint = 2.5d0 + 0.1d0*sin(2d0*pi*( r +0.25d0 ))
      y(2) = 1d0 - 0.5d0*((tanh(const*( xg(j) -xint)/dx(1)))+1d0)
      y(1) = 1d0 - y(2)


      q0(1,j,k) = ryw(1)*y(1) + ryw(2)*y(2)     ! rho
      q0(2,j,k) = u1         ! u
      q0(3,j,k) = 0.0d0          ! v
      q0(4,j,k) = p1    ! p
      q0(5:ndmax,j,k) = ryw(:)*y(:)   ! rhoy

      ! shock though 
      if(xg(j)>3.6d0)then 
        ryw(1) = r2
        q0(1,j,k) = ryw(1)*y(1) + ryw(2)*y(2)  ! rho
        q0(2,j,k) = u2         ! u
        q0(3,j,k) = 0.0d0      ! v
        q0(4,j,k) = p2         ! p
        q0(5:ndmax,j,k) = ryw(:)*y(:)   ! rhoy
      endif
    end do
  end do

  return
end subroutine RTI_init

subroutine Blasius_init(q0)
  use param
  implicit none
  integer j,k,l
  double precision, dimension(ndmax,jmax,kmax,lmax) :: q0
  double precision, dimension(nspecies)   :: y, ryw
  double precision  gamma
  double precision  rinf,uinf,pinf,ainf

  !% Boundary-Layer 
  qinf = 0d0

  gamma = 1.4d0
  rinf = 1.0d0 
  ainf = 1.0d0
  uinf = Mach*ainf
  pinf = rinf*ainf*ainf/gamma

  qinf(1) = rinf
  qinf(2) = uinf
  qinf(3) = 0d0
  qinf(4) = 0d0
  qinf(5) = pinf
  qinf(6) = rinf

  ryw(1) = rinf
  y(1) = 1d0

  do l=1,lmax
    do k=1,kmax
      do j=1,jmax
        ! inlet
        q0(1,j,k,l) = qinf(1)    ! rho
        q0(2,j,k,l) = qinf(2)    ! u
        q0(3,j,k,l) = qinf(3)    ! v
        q0(4,j,k,l) = qinf(4)    ! w
        q0(5,j,k,l) = qinf(5)    ! p
        q0(6:ndmax,j,k,l) = qinf(6:ndmax)   ! rhoy
      end do
    end do
  end do

  return
end subroutine Blasius_init
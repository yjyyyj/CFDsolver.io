subroutine visflux(ql, qr, vf)
  !**********************************************************************
  !*     caluculate right-hand-side                                     *
  !**********************************************************************
  use param
  implicit none
  integer j,k,l,n
  double precision, dimension(ndmax,0:jmax,0:kmax,0:lmax,ndim) :: vf
  double precision,dimension(ndmax,0:(jmax+1),0:(kmax+1),0:(lmax+1),ndim) :: ql, qr
  double precision r1,r2,p1,p2,k1,k2,e1,e2
  double precision m1,m2,g1,g2
  double precision r_fl,kbar,hjbar
  double precision,dimension(nspecies) :: ry1,ry2,h1,h2
  double precision,dimension(nspecies) :: ry_fl,j_fl
  double precision,dimension(nspecies) :: grad_y1
  double precision,dimension(ndim) :: u1,u2,u_fl
  integer,dimension(ndim) :: dl
  double precision diff

  diff = 0.25d0
  ! diff = 0.01d0


  do n=1,ndim
    dl(:)=0
    dl(n)=1
    do l=0,lmax
      do k=0,kmax
        do j=0,jmax
    ! do l=1,lmax-1
    !   do k=1,kmax-1
    !     do j=1,jmax-1
          ! half of gradients Y_i @ j, j+1
          ! grad_y1(:) = 0.5d0*(qr(6:ndmax,j+dl(1),k+dl(2),l+dl(3),n)/qr(1,j+dl(1),k+dl(2),l+dl(3),n)  &
          !                     - ql(6:ndmax,j-dl(1),k-dl(2),l-dl(3),n)/ql(1,j-dl(1),k-dl(2),l-dl(3),n))/dx(n)
          ! grad_y2(:) = 0.5d0*(qr(6:ndmax,j+2*dl(1),k+2*dl(2),l+2*dl(3),n)/qr(1,j+2*dl(1),k+2*dl(2),l+2*dl(3),n) &
          !                     - ql(6:ndmax,j,k,l,n)/ql(1,j,k,l,n))/dx(n)
          !****************************************!
          ! left face q
          u1(1)  = ql(2,j,k,l,n)
          u1(2)  = ql(3,j,k,l,n)
          u1(3)  = ql(4,j,k,l,n)
          p1     = ql(5,j,k,l,n)
          ry1(:) = ql(6:ndmax,j,k,l,n)
          r1     = sum(ry1(:))
          k1 = r1*(u1(1)*u1(1)+u1(2)*u1(2)+u1(3)*u1(3))*0.5d0
          e1 = p1*g1
          m1 = calmw(ql(:,j,k,l,n))
          g1 = calgm(ql(:,j,k,l,n),m1)
      
          ! right face q
          u2(1)  = qr(2,j+dl(1),k+dl(2),l+dl(3),n)
          u2(2)  = qr(3,j+dl(1),k+dl(2),l+dl(3),n)
          u2(3)  = qr(4,j+dl(1),k+dl(2),l+dl(3),n)
          p2     = qr(5,j+dl(1),k+dl(2),l+dl(3),n)
          ry2(:) = qr(6:ndmax,j+dl(1),k+dl(2),l+dl(3),n)
          r2     = sum(ry2(:))
          k2 = r2*(u2(1)*u2(1)+u2(2)*u2(2)+u2(3)*u2(3))*0.5d0
          e2 = p2*g2
          m2 = calmw(qr(:,j+dl(1),k+dl(2),l+dl(3),n))
          g2 = calgm(qr(:,j+dl(1),k+dl(2),l+dl(3),n),m2)  
          
          ry_fl(:) = 0.5d0*(ry1(:) + ry2(:))
          r_fl = sum(ry_fl(:))
          u_fl(:)=0.5d0*(u1(:) + u2(:))
          kbar = 0.5d0*(k1 + k2)
          
          h1(:) = p1*(m1/r1)*(gami(:)/mwi(:) + 1d0/mwi(:))
          h2(:) = p2*(m2/r2)*(gami(:)/mwi(:) + 1d0/mwi(:))

          ! j1(:) = r1*diff*grad_y1(:) - (ry1(:)/r1)*sum(r1*diff*grad_y1(:))
          ! j2(:) = r2*diff*grad_y2(:) - (ry2(:)/r2)*sum(r2*diff*grad_y2(:))
          ! j_fl(:) = 0.5d0*(j2(:) + j1(:)) 
          ! hjbar = sum(0.5d0*(h2(:)*j2(:) + h1(:)*j1(:) ))

          grad_y1(:) = 0.5d0*(ry2(:)/r2 - ry1(:)/r1)/dx(n)  
          j_fl(:) = r_fl*diff*grad_y1(:) - 0.5d0*(ry2(:)/r2 + ry1(:)/r1)*sum(r_fl*diff*grad_y1(:))            
          hjbar = sum(0.5d0*(h2(:) + h1(:)) * j_fl(:))

          ! grad_y1(:) = 0.5d0*(m2/r2*r1/m1*ry2(:) - m1/r1*r2/m2*ry1(:))/dx(n)  
          ! j_fl(:) = diff*grad_y1(:) - 0.5d0*(ry2(:)/r2 + ry1(:)/r1)*sum(diff*grad_y1(:))       
          ! hjbar = 0.5d0*(p2*g2 - p1*g1)

          ! hjbar = sum(0.5d0*(p1*(m1/r1)*(gami(:)/mwi(:) - g1/mwi(:)) + p2*(m2/r2)*(gami(:)/mwi(:) - g2/mwi(:)))*j_fl )!** compati
          
          vf(1,j,k,l,n) = 0.d0
          vf(2,j,k,l,n) = sum(j_fl(:))*u_fl(1)
          vf(3,j,k,l,n) = sum(j_fl(:))*u_fl(2)
          vf(4,j,k,l,n) = sum(j_fl(:))*u_fl(3)
          vf(5,j,k,l,n) = sum(j_fl(:))*kbar + hjbar
          vf(6:ndmax,j,k,l,n) = j_fl(:)

          ! if(n==1 .and. k==2 .and. l==2) write(*,*) j, vf(:,j,k,l,n)
        end do
      end do
    end do
  end do

  return
end subroutine visflux

subroutine visflux_dns(ql, qr, vf, fmu)
  !**********************************************************************
  !*     caluculate right-hand-side                                     *
  !**********************************************************************
  use param
  implicit none
  integer j,k,l,n
  double precision, dimension(ndmax,0:jmax,0:kmax,0:lmax,ndim) :: vf
  double precision,dimension(ndmax,0:(jmax+1),0:(kmax+1),0:(lmax+1),ndim) :: ql, qr
  double precision r1,r2,p1,p2,k1,k2,e1,e2
  double precision m1,m2,g1,g2,c21,c22
  double precision r_fl,kbar
  double precision,dimension(nspecies) :: ry1,ry2,h1,h2
  double precision,dimension(nspecies) :: ry_fl,j_fl,hjbar
  double precision j_fl_sum
  double precision,dimension(nspecies) :: grad_y1
  double precision,dimension(ndim) :: u1,u2,u_fl
  double precision,dimension(ndim) :: dux,duy,duz,duh
  double precision,dimension(ndim) :: tau
  double precision ux,uy,uz,vx,vy,vz,wx,wy,wz,dch
  double precision txx,txy,tzx,tyy,tyz,tzz, visdiv
  double precision tau_u, qheat
  double precision,dimension(0:jmax+1,ndim,ndim) :: grad_uj
  double precision,dimension(0:kmax+1,ndim,ndim) :: grad_uk
  double precision,dimension(0:lmax+1,ndim,ndim) :: grad_ul
  integer,dimension(ndim) :: dl
  double precision diff,lamb
  double precision,dimension(0:(jmax+1),0:(kmax+1),0:(lmax+1)) :: fmu           ! viscous coefficient
  double precision fmuave

  diff = 0.0d0
  ! diff = 0.01d0
  lamb = 0.0d0

  !** loop start ****
  n = 1 
  dl(:)=0
  dl(n)=1

  do l=1,lmax
    do k=1,kmax
      ! for boundary treatment @ j=0

      j=0
      dux(1) = qr(2,j+1,k,l,n) - ql(2,j,k,l,n)
      dux(2) = qr(3,j+1,k,l,n) - ql(3,j,k,l,n)
      dux(3) = qr(4,j+1,k,l,n) - ql(4,j,k,l,n)
      duy(1) = qr(2,j,k+1,l,n) - ql(2,j,k-1,l,n)
      duy(2) = qr(3,j,k+1,l,n) - ql(3,j,k-1,l,n)
      duy(3) = qr(4,j,k+1,l,n) - ql(4,j,k-1,l,n)
      duz(1) = qr(2,j,k,l+1,n) - ql(2,j,k,l-1,n)
      duz(2) = qr(3,j,k,l+1,n) - ql(3,j,k,l-1,n)
      duz(3) = qr(4,j,k,l+1,n) - ql(4,j,k,l-1,n)

      grad_uj(j,:,1) = 0.5d0*dux(:) * int(1-dl(n))
      grad_uj(j,:,2) = 0.5d0*duy(:) * int(1-dl(n))
      grad_uj(j,:,3) = 0.5d0*duz(:) * int(1-dl(n))

      ! for boundary treatment @ j=jmax+1
      j=jmax+1
      dux(1) = qr(2,j,k,l,n) - ql(2,j-1,k,l,n)
      dux(2) = qr(3,j,k,l,n) - ql(3,j-1,k,l,n)
      dux(3) = qr(4,j,k,l,n) - ql(4,j-1,k,l,n)
      duy(1) = qr(2,j,k+1,l,n) - ql(2,j,k-1,l,n)
      duy(2) = qr(3,j,k+1,l,n) - ql(3,j,k-1,l,n)
      duy(3) = qr(4,j,k+1,l,n) - ql(4,j,k-1,l,n)
      duz(1) = qr(2,j,k,l+1,n) - ql(2,j,k,l-1,n)
      duz(2) = qr(3,j,k,l+1,n) - ql(3,j,k,l-1,n)
      duz(3) = qr(4,j,k,l+1,n) - ql(4,j,k,l-1,n)

      grad_uj(j,:,1) = 0.5d0*dux(:) * int(1-dl(n))
      grad_uj(j,:,2) = 0.5d0*duy(:) * int(1-dl(n))
      grad_uj(j,:,3) = 0.5d0*duz(:) * int(1-dl(n))

      do j=1,jmax
        ! half of gradients (u,v,w,c2) of non-horizontal 
        dux(1) = qr(2,j+1,k,l,n) - ql(2,j-1,k,l,n)
        dux(2) = qr(3,j+1,k,l,n) - ql(3,j-1,k,l,n)
        dux(3) = qr(4,j+1,k,l,n) - ql(4,j-1,k,l,n)
        duy(1) = qr(2,j,k+1,l,n) - ql(2,j,k-1,l,n)
        duy(2) = qr(3,j,k+1,l,n) - ql(3,j,k-1,l,n)
        duy(3) = qr(4,j,k+1,l,n) - ql(4,j,k-1,l,n)
        duz(1) = qr(2,j,k,l+1,n) - ql(2,j,k,l-1,n)
        duz(2) = qr(3,j,k,l+1,n) - ql(3,j,k,l-1,n)
        duz(3) = qr(4,j,k,l+1,n) - ql(4,j,k,l-1,n)

        grad_uj(j,:,1) = 0.5d0*dux(:) * int(1-dl(n))
        grad_uj(j,:,2) = 0.5d0*duy(:) * int(1-dl(n))
        grad_uj(j,:,3) = 0.5d0*duz(:) * int(1-dl(n))
      enddo 
      !
      do j=0,jmax
        !****************************************!
        ! left face q
        u1(1)  = ql(2,j,k,l,n)
        u1(2)  = ql(3,j,k,l,n)
        u1(3)  = ql(4,j,k,l,n)
        p1     = ql(5,j,k,l,n)
        ry1(:) = ql(6:ndmax,j,k,l,n)
        r1     = sum(ry1(:))
        k1 = r1*(u1(1)*u1(1)+u1(2)*u1(2)+u1(3)*u1(3))*0.5d0
        e1 = p1*g1
        m1 = calmw(ql(:,j,k,l,n))
        g1 = calgm(ql(:,j,k,l,n),m1)
        c21= p1/r1*(1d0/g1+1d0) 
    
        ! right face q
        u2(1)  = qr(2,j+dl(1),k+dl(2),l+dl(3),n)
        u2(2)  = qr(3,j+dl(1),k+dl(2),l+dl(3),n)
        u2(3)  = qr(4,j+dl(1),k+dl(2),l+dl(3),n)
        p2     = qr(5,j+dl(1),k+dl(2),l+dl(3),n)
        ry2(:) = qr(6:ndmax,j+dl(1),k+dl(2),l+dl(3),n)
        r2     = sum(ry2(:))
        k2 = r2*(u2(1)*u2(1)+u2(2)*u2(2)+u2(3)*u2(3))*0.5d0
        e2 = p2*g2
        m2 = calmw(qr(:,j+dl(1),k+dl(2),l+dl(3),n))
        g2 = calgm(qr(:,j+dl(1),k+dl(2),l+dl(3),n),m2)  
        c22= p2/r2*(1d0/g2+1d0) 

        r_fl = 0.5d0*(r2 + r1)
        ry_fl(:) = 0.5d0*(ry2(:) + ry1(:))
        u_fl(:)=0.5d0*(u1(:) + u2(:))
        kbar = 0.5d0*(k1 + k2)

        !*** diffusion flux Ji ***************************
        ! Cook 2009
        ! j1(:) = r1*diff*grad_y1(:) - (ry1(:)/r1)*sum(r1*diff*grad_y1(:))
        ! j2(:) = r2*diff*grad_y2(:) - (ry2(:)/r2)*sum(r2*diff*grad_y2(:))
        ! j_fl(:) = 0.5d0*(j2(:) + j1(:)) 

        grad_y1(:) = 0.5d0*(ry2(:)/r2 - ry1(:)/r1)/dx(n)
        j_fl(:) = r_fl*diff*grad_y1(:) - 0.5d0*(ry2(:)/r2 + ry1(:)/r1)*sum(r_fl*diff*grad_y1(:))
        j_fl_sum = sum(j_fl(:))

        !*** stress tensor flux Tau ************************
        ! gradients (u,v,w) of horizontal
        duh(:) = u2(:) - u1(:)

        ux = 0.5d0*(duh(1)*dl(1) + grad_uj(j,1,1) + grad_uj(j+1,1,1))/dx(1)
        vx = 0.5d0*(duh(2)*dl(1) + grad_uj(j,2,1) + grad_uj(j+1,2,1))/dx(1)
        wx = 0.5d0*(duh(3)*dl(1) + grad_uj(j,3,1) + grad_uj(j+1,3,1))/dx(1)
        uy = 0.5d0*(duh(1)*dl(2) + grad_uj(j,1,2) + grad_uj(j+1,1,2))/dx(2)
        vy = 0.5d0*(duh(2)*dl(2) + grad_uj(j,2,2) + grad_uj(j+1,2,2))/dx(2)
        wy = 0.5d0*(duh(3)*dl(2) + grad_uj(j,3,2) + grad_uj(j+1,3,2))/dx(2)
        uz = 0.5d0*(duh(1)*dl(3) + grad_uj(j,1,3) + grad_uj(j+1,1,3))/dx(3)
        vz = 0.5d0*(duh(2)*dl(3) + grad_uj(j,2,3) + grad_uj(j+1,2,3))/dx(3)
        wz = 0.5d0*(duh(3)*dl(3) + grad_uj(j,3,3) + grad_uj(j+1,3,3))/dx(3)

        visdiv = 2d0/3d0*(ux+vy+wz)
        txx = 2.0d0*ux - visdiv
        tyy = 2.0d0*vy - visdiv
        tzz = 2.0d0*wz - visdiv
        txy = vx + uy
        tyz = wy + vz
        tzx = uz + wx


        fmuave   = 0.5d0*(   fmu(j,k,l) +  fmu(j+dl(1),k+dl(2),l+dl(3)) )* Mach/Re
        ! diffave   = 0.5d0*(   diff(j,k,l) +  diff(j+1,k,l) )
        ! gpr2    = fmuave / (Pr*(gam-1))
        ! dre2    = Ma/Re
        ! dpe2    = Ma/Pe

        tau(1) = fmuave *(txx*dl(1) + txy*dl(2) + tzx*dl(3))
        tau(2) = fmuave *(txy*dl(1) + tyy*dl(2) + tyz*dl(3))
        tau(3) = fmuave *(tzx*dl(1) + tyz*dl(2) + tzz*dl(3))
        tau_u  = u_fl(1)*tau(1) + u_fl(2)*tau(2) + u_fl(3)*tau(3)

        !*** heat diffusion flux ***************************
        ! heat flux q
        dch = c22 - c21
        qheat  = lamb*0.5d0*(c22 - c21)/dx(n)

        ! enthalpy flux sum(hJ)
        h1(:) = p1*(m1/r1)*(gami(:)/mwi(:) + 1d0/mwi(:))
        h2(:) = p2*(m2/r2)*(gami(:)/mwi(:) + 1d0/mwi(:))
        ! hjbar(:) = 0.5d0*(h2(:)*j2(:) + h1(:)*j1(:) )
        hjbar(:) = 0.5d0*(h2(:) + h1(:))*j_fl(:)

        !!* compati
        ! hjbar(:) = 0.5d0*(p1*(m1/r1)*(gami(:)/mwi(:) - g1/mwi(:)) + p2*(m2/r2)*(gami(:)/mwi(:) - g2/mwi(:)))*j_fl 
        
        !*** construct viscous flux ***************************
        vf(1,j,k,l,n) = 0.d0
        vf(2,j,k,l,n) = tau(1) + j_fl_sum *u_fl(1)
        vf(3,j,k,l,n) = tau(2) + j_fl_sum *u_fl(2)
        vf(4,j,k,l,n) = tau(3) + j_fl_sum *u_fl(3)
        vf(5,j,k,l,n) = tau_u + qheat + sum(hjbar(:)) +j_fl_sum*kbar
        vf(6:ndmax,j,k,l,n) = j_fl(:)
      end do
    end do
  end do


  n = 2 
  dl(:)=0
  dl(n)=1

  do l=1,lmax
    do j=1,jmax
      ! for boundary treatment @ k=0
      k=0
      dux(1) = qr(2,j+1,k,l,n) - ql(2,j-1,k,l,n)
      dux(2) = qr(3,j+1,k,l,n) - ql(3,j-1,k,l,n)
      dux(3) = qr(4,j+1,k,l,n) - ql(4,j-1,k,l,n)
      duy(1) = qr(2,j,k+1,l,n) - ql(2,j,k,l,n)
      duy(2) = qr(3,j,k+1,l,n) - ql(3,j,k,l,n)
      duy(3) = qr(4,j,k+1,l,n) - ql(4,j,k,l,n)
      duz(1) = qr(2,j,k,l+1,n) - ql(2,j,k,l-1,n)
      duz(2) = qr(3,j,k,l+1,n) - ql(3,j,k,l-1,n)
      duz(3) = qr(4,j,k,l+1,n) - ql(4,j,k,l-1,n)

      grad_uk(k,:,1) = 0.5d0*dux(:) * int(1-dl(n))
      grad_uk(k,:,2) = 0.5d0*duy(:) * int(1-dl(n))
      grad_uk(k,:,3) = 0.5d0*duz(:) * int(1-dl(n))
      
      ! for boundary treatment @ k=kmax+1
      k=kmax+1
      dux(1) = qr(2,j+1,k,l,n) - ql(2,j-1,k,l,n)
      dux(2) = qr(3,j+1,k,l,n) - ql(3,j-1,k,l,n)
      dux(3) = qr(4,j+1,k,l,n) - ql(4,j-1,k,l,n)
      duy(1) = qr(2,j,k,l,n) - ql(2,j,k-1,l,n)
      duy(2) = qr(3,j,k,l,n) - ql(3,j,k-1,l,n)
      duy(3) = qr(4,j,k,l,n) - ql(4,j,k-1,l,n)
      duz(1) = qr(2,j,k,l+1,n) - ql(2,j,k,l-1,n)
      duz(2) = qr(3,j,k,l+1,n) - ql(3,j,k,l-1,n)
      duz(3) = qr(4,j,k,l+1,n) - ql(4,j,k,l-1,n)

      grad_uk(k,:,1) = 0.5d0*dux(:) * int(1-dl(n))
      grad_uk(k,:,2) = 0.5d0*duy(:) * int(1-dl(n))
      grad_uk(k,:,3) = 0.5d0*duz(:) * int(1-dl(n))

      do k=1,kmax
        ! half of gradients (u,v,w,c2) of non-horizontal 
        dux(1) = qr(2,j+1,k,l,n) - ql(2,j-1,k,l,n)
        dux(2) = qr(3,j+1,k,l,n) - ql(3,j-1,k,l,n)
        dux(3) = qr(4,j+1,k,l,n) - ql(4,j-1,k,l,n)
        duy(1) = qr(2,j,k+1,l,n) - ql(2,j,k-1,l,n)
        duy(2) = qr(3,j,k+1,l,n) - ql(3,j,k-1,l,n)
        duy(3) = qr(4,j,k+1,l,n) - ql(4,j,k-1,l,n)
        duz(1) = qr(2,j,k,l+1,n) - ql(2,j,k,l-1,n)
        duz(2) = qr(3,j,k,l+1,n) - ql(3,j,k,l-1,n)
        duz(3) = qr(4,j,k,l+1,n) - ql(4,j,k,l-1,n)

        grad_uk(k,:,1) = 0.5d0*dux(:) * int(1-dl(n))
        grad_uk(k,:,2) = 0.5d0*duy(:) * int(1-dl(n))
        grad_uk(k,:,3) = 0.5d0*duz(:) * int(1-dl(n))
      enddo 
      !
      do k=0,kmax
        !****************************************!
        ! left face q
        u1(1)  = ql(2,j,k,l,n)
        u1(2)  = ql(3,j,k,l,n)
        u1(3)  = ql(4,j,k,l,n)
        p1     = ql(5,j,k,l,n)
        ry1(:) = ql(6:ndmax,j,k,l,n)
        r1     = sum(ry1(:))
        k1 = r1*(u1(1)*u1(1)+u1(2)*u1(2)+u1(3)*u1(3))*0.5d0
        e1 = p1*g1
        m1 = calmw(ql(:,j,k,l,n))
        g1 = calgm(ql(:,j,k,l,n),m1)
        c21= p1/r1*(1d0/g1+1d0) 
    
        ! right face q
        u2(1)  = qr(2,j+dl(1),k+dl(2),l+dl(3),n)
        u2(2)  = qr(3,j+dl(1),k+dl(2),l+dl(3),n)
        u2(3)  = qr(4,j+dl(1),k+dl(2),l+dl(3),n)
        p2     = qr(5,j+dl(1),k+dl(2),l+dl(3),n)
        ry2(:) = qr(6:ndmax,j+dl(1),k+dl(2),l+dl(3),n)
        r2     = sum(ry2(:))
        k2 = r2*(u2(1)*u2(1)+u2(2)*u2(2)+u2(3)*u2(3))*0.5d0
        e2 = p2*g2
        m2 = calmw(qr(:,j+dl(1),k+dl(2),l+dl(3),n))
        g2 = calgm(qr(:,j+dl(1),k+dl(2),l+dl(3),n),m2)  
        c22= p2/r2*(1d0/g2+1d0) 

        r_fl = 0.5d0*(r2 + r1)
        ry_fl(:) = 0.5d0*(ry2(:) + ry1(:))
        u_fl(:)=0.5d0*(u1(:) + u2(:))
        kbar = 0.5d0*(k1 + k2)

        !*** diffusion flux Ji ***************************
        ! Cook 2009
        ! j1(:) = r1*diff*grad_y1(:) - (ry1(:)/r1)*sum(r1*diff*grad_y1(:))
        ! j2(:) = r2*diff*grad_y2(:) - (ry2(:)/r2)*sum(r2*diff*grad_y2(:))
        ! j_fl(:) = 0.5d0*(j2(:) + j1(:)) 

        grad_y1(:) = 0.5d0*(ry2(:)/r2 - ry1(:)/r1)/dx(n)
        j_fl(:) = r_fl*diff*grad_y1(:) - 0.5d0*(ry2(:)/r2 + ry1(:)/r1)*sum(r_fl*diff*grad_y1(:))
        j_fl_sum = sum(j_fl(:))

        !*** stress tensor flux Tau ************************
        ! gradients (u,v,w) of horizontal
        duh(:) = u2(:) - u1(:)

        ux = 0.5d0*(duh(1)*dl(1) + grad_uk(k,1,1) + grad_uk(k+1,1,1))/dx(1)
        vx = 0.5d0*(duh(2)*dl(1) + grad_uk(k,2,1) + grad_uk(k+1,2,1))/dx(1)
        wx = 0.5d0*(duh(3)*dl(1) + grad_uk(k,3,1) + grad_uk(k+1,3,1))/dx(1)
        uy = 0.5d0*(duh(1)*dl(2) + grad_uk(k,1,2) + grad_uk(k+1,1,2))/dx(2)
        vy = 0.5d0*(duh(2)*dl(2) + grad_uk(k,2,2) + grad_uk(k+1,2,2))/dx(2)
        wy = 0.5d0*(duh(3)*dl(2) + grad_uk(k,3,2) + grad_uk(k+1,3,2))/dx(2)
        uz = 0.5d0*(duh(1)*dl(3) + grad_uk(k,1,3) + grad_uk(k+1,1,3))/dx(3)
        vz = 0.5d0*(duh(2)*dl(3) + grad_uk(k,2,3) + grad_uk(k+1,2,3))/dx(3)
        wz = 0.5d0*(duh(3)*dl(3) + grad_uk(k,3,3) + grad_uk(k+1,3,3))/dx(3)

        visdiv = 2d0/3d0*(ux+vy+wz)
        txx = 2.0d0*ux - visdiv
        tyy = 2.0d0*vy - visdiv
        tzz = 2.0d0*wz - visdiv
        txy = vx + uy
        tyz = wy + vz
        tzx = uz + wx

        fmuave   = 0.5d0*( fmu(j,k,l) +  fmu(j+dl(1),k+dl(2),l+dl(3)) )* Mach/Re

        tau(1) = fmuave *(txx*dl(1) + txy*dl(2) + tzx*dl(3))
        tau(2) = fmuave *(txy*dl(1) + tyy*dl(2) + tyz*dl(3))
        tau(3) = fmuave *(tzx*dl(1) + tyz*dl(2) + tzz*dl(3))
        tau_u  = u_fl(1)*tau(1) + u_fl(2)*tau(2) + u_fl(3)*tau(3)

        !*** heat diffusion flux ***************************
        ! heat flux q
        dch = c22 - c21
        qheat  = lamb*0.5d0*(c22 - c21)/dx(n)

        ! enthalpy flux sum(hJ)
        h1(:) = p1*(m1/r1)*(gami(:)/mwi(:) + 1d0/mwi(:))
        h2(:) = p2*(m2/r2)*(gami(:)/mwi(:) + 1d0/mwi(:))
        ! hjbar(:) = 0.5d0*(h2(:)*j2(:) + h1(:)*j1(:) )
        hjbar(:) = 0.5d0*(h2(:) + h1(:))*j_fl(:)

        !!* compati
        ! hjbar(:) = 0.5d0*(p1*(m1/r1)*(gami(:)/mwi(:) - g1/mwi(:)) + p2*(m2/r2)*(gami(:)/mwi(:) - g2/mwi(:)))*j_fl 
        
        !*** construct viscous flux ***************************
        vf(1,j,k,l,n) = 0.d0
        vf(2,j,k,l,n) = tau(1) + j_fl_sum *u_fl(1)
        vf(3,j,k,l,n) = tau(2) + j_fl_sum *u_fl(2)
        vf(4,j,k,l,n) = tau(3) + j_fl_sum *u_fl(3)
        vf(5,j,k,l,n) = tau_u + qheat + sum(hjbar(:)) +j_fl_sum*kbar
        vf(6:ndmax,j,k,l,n) = j_fl(:)
      end do
    end do
  end do

  n = 3 
  dl(:)=0
  dl(n)=1

  do k=1,kmax
    do j=1,jmax
      ! for boundary treatment @ l=0
      l=0
      dux(1) = qr(2,j+1,k,l,n) - ql(2,j-1,k,l,n)
      dux(2) = qr(3,j+1,k,l,n) - ql(3,j-1,k,l,n)
      dux(3) = qr(4,j+1,k,l,n) - ql(4,j-1,k,l,n)
      duy(1) = qr(2,j,k+1,l,n) - ql(2,j,k-1,l,n)
      duy(2) = qr(3,j,k+1,l,n) - ql(3,j,k-1,l,n)
      duy(3) = qr(4,j,k+1,l,n) - ql(4,j,k-1,l,n)
      duz(1) = qr(2,j,k,l+1,n) - ql(2,j,k,l,n)
      duz(2) = qr(3,j,k,l+1,n) - ql(3,j,k,l,n)
      duz(3) = qr(4,j,k,l+1,n) - ql(4,j,k,l,n)

      grad_ul(l,:,1) = 0.5d0*dux(:) * abs(1-dl(n))
      grad_ul(l,:,2) = 0.5d0*duy(:) * abs(1-dl(n))
      grad_ul(l,:,3) = 0.5d0*duz(:) * abs(1-dl(n))

      ! for boundary treatment @ l=lmax+1
      l=lmax+1
      dux(1) = qr(2,j+1,k,l,n) - ql(2,j-1,k,l,n)
      dux(2) = qr(3,j+1,k,l,n) - ql(3,j-1,k,l,n)
      dux(3) = qr(4,j+1,k,l,n) - ql(4,j-1,k,l,n)
      duy(1) = qr(2,j,k+1,l,n) - ql(2,j,k-1,l,n)
      duy(2) = qr(3,j,k+1,l,n) - ql(3,j,k-1,l,n)
      duy(3) = qr(4,j,k+1,l,n) - ql(4,j,k-1,l,n)
      duz(1) = qr(2,j,k,l,n) - ql(2,j,k,l-1,n)
      duz(2) = qr(3,j,k,l,n) - ql(3,j,k,l-1,n)
      duz(3) = qr(4,j,k,l,n) - ql(4,j,k,l-1,n)

      grad_ul(l,:,1) = 0.5d0*dux(:) * abs(1-dl(n))
      grad_ul(l,:,2) = 0.5d0*duy(:) * abs(1-dl(n))
      grad_ul(l,:,3) = 0.5d0*duz(:) * abs(1-dl(n))

      do l=1,lmax
        ! half of gradients (u,v,w,c2) of non-horizontal 
        dux(1) = qr(2,j+1,k,l,n) - ql(2,j-1,k,l,n)
        dux(2) = qr(3,j+1,k,l,n) - ql(3,j-1,k,l,n)
        dux(3) = qr(4,j+1,k,l,n) - ql(4,j-1,k,l,n)
        duy(1) = qr(2,j,k+1,l,n) - ql(2,j,k-1,l,n)
        duy(2) = qr(3,j,k+1,l,n) - ql(3,j,k-1,l,n)
        duy(3) = qr(4,j,k+1,l,n) - ql(4,j,k-1,l,n)
        duz(1) = qr(2,j,k,l+1,n) - ql(2,j,k,l-1,n)
        duz(2) = qr(3,j,k,l+1,n) - ql(3,j,k,l-1,n)
        duz(3) = qr(4,j,k,l+1,n) - ql(4,j,k,l-1,n)

        grad_ul(l,:,1) = 0.5d0*dux(:) * abs(1-dl(n))
        grad_ul(l,:,2) = 0.5d0*duy(:) * abs(1-dl(n))
        grad_ul(l,:,3) = 0.5d0*duz(:) * abs(1-dl(n))
      enddo 
      !
      do l=0,lmax
        !****************************************!
        ! left face q
        u1(1)  = ql(2,j,k,l,n)
        u1(2)  = ql(3,j,k,l,n)
        u1(3)  = ql(4,j,k,l,n)
        p1     = ql(5,j,k,l,n)
        ry1(:) = ql(6:ndmax,j,k,l,n)
        r1     = sum(ry1(:))
        k1 = r1*(u1(1)*u1(1)+u1(2)*u1(2)+u1(3)*u1(3))*0.5d0
        e1 = p1*g1
        m1 = calmw(ql(:,j,k,l,n))
        g1 = calgm(ql(:,j,k,l,n),m1)
        c21= p1/r1*(1d0/g1+1d0) 
    
        ! right face q
        u2(1)  = qr(2,j+dl(1),k+dl(2),l+dl(3),n)
        u2(2)  = qr(3,j+dl(1),k+dl(2),l+dl(3),n)
        u2(3)  = qr(4,j+dl(1),k+dl(2),l+dl(3),n)
        p2     = qr(5,j+dl(1),k+dl(2),l+dl(3),n)
        ry2(:) = qr(6:ndmax,j+dl(1),k+dl(2),l+dl(3),n)
        r2     = sum(ry2(:))
        k2 = r2*(u2(1)*u2(1)+u2(2)*u2(2)+u2(3)*u2(3))*0.5d0
        e2 = p2*g2
        m2 = calmw(qr(:,j+dl(1),k+dl(2),l+dl(3),n))
        g2 = calgm(qr(:,j+dl(1),k+dl(2),l+dl(3),n),m2)  
        c22= p2/r2*(1d0/g2+1d0) 

        r_fl = 0.5d0*(r2 + r1)
        ry_fl(:) = 0.5d0*(ry2(:) + ry1(:))
        u_fl(:)=0.5d0*(u1(:) + u2(:))
        kbar = 0.5d0*(k1 + k2)

        !*** diffusion flux Ji ***************************
        ! Cook 2009
        ! j1(:) = r1*diff*grad_y1(:) - (ry1(:)/r1)*sum(r1*diff*grad_y1(:))
        ! j2(:) = r2*diff*grad_y2(:) - (ry2(:)/r2)*sum(r2*diff*grad_y2(:))
        ! j_fl(:) = 0.5d0*(j2(:) + j1(:)) 

        grad_y1(:) = 0.5d0*(ry2(:)/r2 - ry1(:)/r1)/dx(n)
        j_fl(:) = r_fl*diff*grad_y1(:) - 0.5d0*(ry2(:)/r2 + ry1(:)/r1)*sum(r_fl*diff*grad_y1(:))
        j_fl_sum = sum(j_fl(:))

        !*** stress tensor flux Tau ************************
        ! gradients (u,v,w) of horizontal
        duh(:) = u2(:) - u1(:)

        ux = 0.5d0*(duh(1)*dl(1) + grad_ul(l,1,1) + grad_ul(l+1,1,1))/dx(1)
        vx = 0.5d0*(duh(2)*dl(1) + grad_ul(l,2,1) + grad_ul(l+1,2,1))/dx(1)
        wx = 0.5d0*(duh(3)*dl(1) + grad_ul(l,3,1) + grad_ul(l+1,3,1))/dx(1)
        uy = 0.5d0*(duh(1)*dl(2) + grad_ul(l,1,2) + grad_ul(l+1,1,2))/dx(2)
        vy = 0.5d0*(duh(2)*dl(2) + grad_ul(l,2,2) + grad_ul(l+1,2,2))/dx(2)
        wy = 0.5d0*(duh(3)*dl(2) + grad_ul(l,3,2) + grad_ul(l+1,3,2))/dx(2)
        uz = 0.5d0*(duh(1)*dl(3) + grad_ul(l,1,3) + grad_ul(l+1,1,3))/dx(3)
        vz = 0.5d0*(duh(2)*dl(3) + grad_ul(l,2,3) + grad_ul(l+1,2,3))/dx(3)
        wz = 0.5d0*(duh(3)*dl(3) + grad_ul(l,3,3) + grad_ul(l+1,3,3))/dx(3)

        visdiv = 2d0/3d0*(ux+vy+wz)
        txx = 2.0d0*ux - visdiv
        tyy = 2.0d0*vy - visdiv
        tzz = 2.0d0*wz - visdiv
        txy = vx + uy
        tyz = wy + vz
        tzx = uz + wx

        fmuave   = 0.5d0*( fmu(j,k,l) +  fmu(j+dl(1),k+dl(2),l+dl(3)) )* Mach/Re

        tau(1) = fmuave *(txx*dl(1) + txy*dl(2) + tzx*dl(3))
        tau(2) = fmuave *(txy*dl(1) + tyy*dl(2) + tyz*dl(3))
        tau(3) = fmuave *(tzx*dl(1) + tyz*dl(2) + tzz*dl(3))
        tau_u  = u_fl(1)*tau(1) + u_fl(2)*tau(2) + u_fl(3)*tau(3)

        !*** heat diffusion flux ***************************
        ! heat flux q
        dch = c22 - c21
        qheat  = lamb*0.5d0*(c22 - c21)/dx(n)

        ! enthalpy flux sum(hJ)
        h1(:) = p1*(m1/r1)*(gami(:)/mwi(:) + 1d0/mwi(:))
        h2(:) = p2*(m2/r2)*(gami(:)/mwi(:) + 1d0/mwi(:))
        ! hjbar(:) = 0.5d0*(h2(:)*j2(:) + h1(:)*j1(:) )
        hjbar(:) = 0.5d0*(h2(:) + h1(:))*j_fl(:)

        !!* compati
        ! hjbar(:) = 0.5d0*(p1*(m1/r1)*(gami(:)/mwi(:) - g1/mwi(:)) + p2*(m2/r2)*(gami(:)/mwi(:) - g2/mwi(:)))*j_fl 
        
        !*** construct viscous flux ***************************
        vf(1,j,k,l,n) = 0.d0
        vf(2,j,k,l,n) = tau(1) + j_fl_sum *u_fl(1)
        vf(3,j,k,l,n) = tau(2) + j_fl_sum *u_fl(2)
        vf(4,j,k,l,n) = tau(3) + j_fl_sum *u_fl(3)
        vf(5,j,k,l,n) = tau_u + qheat + sum(hjbar(:)) +j_fl_sum*kbar
        vf(6:ndmax,j,k,l,n) = j_fl(:)
      end do
    end do
  end do

  return
end subroutine visflux_dns

subroutine calc_viscoefs(q, fmu)
  !**********************************************************************
  !*     caluculate right-hand-side                                     *
  !**********************************************************************
  use param
  implicit none
  integer j,k,l
  double precision, dimension(0:(jmax+1),0:(kmax+1),0:(lmax+1)) :: fmu
  double precision,dimension(ndmax,0:(jmax+1),0:(kmax+1),0:(lmax+1)) :: q
  double precision c_1,c_2,gm

  do l=1,lmax 
    do k=1,kmax 
      do j=1,jmax
        gm  = 1d0/gam(j,k,l) + 1d0
        c_2 = gm*q(5,j,k,l)/q(1,j,k,l)
        c_1 = sqrt(c_2)
        fmu(j,k,l) = vflag*c2bp*c_1**3/(c2b+c_2)
      enddo
    enddo
  enddo

  return
end subroutine calc_viscoefs
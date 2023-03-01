subroutine visflux_numerical(ql, qr, vf, fmu)
  !**********************************************************************
  !*     caluculate numerical viscuss term                              *
  !**********************************************************************
  use param_mod
  implicit none
  integer j,k,l,n
  double precision, dimension(ndmax,0:jmax,0:kmax,0:lmax,ndim) :: vf
  double precision,dimension(ndmax,0:(jmax+1),0:(kmax+1),0:(lmax+1),ndim) :: ql, qr
  double precision :: fmu           ! viscous coefficient
  double precision r1,r2,p1,p2,k1,k2,e1,e2
  double precision m1,m2,g1,g2
  double precision T1,T2
  double precision r_fl,pbar,kbar,hjbar,qheat_bar
  double precision,dimension(nspecies) :: ry1,ry2,h1,h2,phi1,phi2
  double precision,dimension(nspecies) :: rj_fl,j_fl
  double precision,dimension(nspecies) :: grad_y,diff
  double precision,dimension(ndim) :: u1,u2,u_fl
  integer,dimension(ndim) :: dl


  do n=1,ndim
    dl(:)=0
    dl(n)=1
    ! diff(:) = 0.005d0/dx(n)
    diff(:) = 0.5d0
    !$OMP parallel do default(none) &
    !$OMP & firstprivate(dl,n) &
    !$OMP & shared(vf,ql,qr,diff,jmax,kmax,lmax,ndmax) &
    !$OMP & shared(dx,gami,mwi,Ru) &
    !$OMP & private(r1,r2,p1,p2,k1,k2,e1,e2,m1,m2,g1,g2,T1,T2) & 
    !$OMP & private(r_fl,pbar,kbar,hjbar,qheat_bar) & 
    !$OMP & private(ry1,ry2,h1,h2,phi1,phi2,u1,u2,u_fl,rj_fl,j_fl,grad_y) 
    do l=0,lmax
      do k=0,kmax
        do j=0,jmax
          !****************************************!
          ! left face q
          u1(1)  = ql(2,j,k,l,n)
          u1(2)  = ql(3,j,k,l,n)
          u1(3)  = ql(4,j,k,l,n)
          p1     = ql(5,j,k,l,n)
          ry1(:) = ql(6:ndmax,j,k,l,n)
          r1     = sum(ry1(:))
          k1 = (u1(1)*u1(1)+u1(2)*u1(2)+u1(3)*u1(3))*0.5d0
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
          k2 = (u2(1)*u2(1)+u2(2)*u2(2)+u2(3)*u2(3))*0.5d0
          e2 = p2*g2
          m2 = calmw(qr(:,j+dl(1),k+dl(2),l+dl(3),n))
          g2 = calgm(qr(:,j+dl(1),k+dl(2),l+dl(3),n),m2)  
          
          phi1(:) = m1*ry1(:)/r1/mwi(:)
          phi2(:) = m2*ry2(:)/r2/mwi(:)

          r_fl = 0.5d0*(r1 + r2)
          u_fl(:)=0.5d0*(u1(:) + u2(:))
          pbar = 0.5d0*(p1 + p2)
          kbar = 0.5d0*(k1 + k2)
          
          !****** mass diff ********************************* 
          j_fl(:)  = 0.5d0*diff(:)*(phi2(:) - phi1(:)) ! nabla phi_i
          rj_fl(:) = 0.5d0*(r2/m2*mwi(:) + r1/m1*mwi(:))*j_fl(:) ! rho_i|_{j+1/2} j_i|_{j+1/2}

          ! grad_y(:) = 0.5d0*(ry2(:)/r2 - ry1(:)/r1)  
          ! rj_fl(:) = r_fl*diff(:)*grad_y(:) - 0.5d0*(ry2(:)/r2 + ry1(:)/r1)*sum(r_fl*diff*grad_y(:))  ! Cook, 2007 
          
          ! rj_fl(:) = 0.5d0*(r1+r2) * 0.5d0*diff(:)*(ry2(:)/r2 - ry1(:)/r1) ! rho|_{j+1/2} Y_i|_{j+1/2}

          ! rj_fl(:) = 0.5d0*diff(:)*(ry2(:) - ry1(:))   ! nabla rhoY_i

          r_fl = sum(rj_fl(:))
          
          !****** enthalpy diff ********************************* 
          h1(:) = p1*(m1/r1)*(gami(:)/mwi(:) + 1d0/mwi(:)) ! Y_i base
          h2(:) = p2*(m2/r2)*(gami(:)/mwi(:) + 1d0/mwi(:)) ! Y_i base
          hjbar = sum(0.5d0*(h2(:) + h1(:)) * rj_fl(:)) ! H_i|_12 = h_i|_12*J_i|_12

          ! h1(:) = p1*gami(:) + p1 ! phi_i base
          ! h2(:) = p2*gami(:) + p2 ! phi_i base
          ! hjbar = sum( 0.5d0*(h2(:) + h1(:))*j_fl(:) ) ! H_i|_12 = rho_i h_i|_12*j_i|_12
          
          ! hjbar = 0.5d0*pbar*diff(1)*(g2 - g1) ! T=const

          !****** heat diff flux ********************************* 
          T1 = p1/Ru*m1/r1
          T2 = p2/Ru*m2/r2
          qheat_bar = 0.5d0*10d0*diff(1)*(T2 - T1) ! q = d(k dT)

          !****** construct flux ********************************* 
          vf(1,j,k,l,n) = 0.d0
          vf(2,j,k,l,n) = r_fl*u_fl(1)*dl(n)
          vf(3,j,k,l,n) = r_fl*u_fl(2)*dl(n)
          vf(4,j,k,l,n) = r_fl*u_fl(3)*dl(n)
          vf(5,j,k,l,n) = r_fl*kbar*dl(n) + hjbar*dl(n) + qheat_bar*dl(n)
          vf(6:ndmax,j,k,l,n) = rj_fl(:)*dl(n)

        end do
      end do
    end do
  end do

  return
end subroutine visflux_numerical

subroutine visflux_dns(ql, qr, vf, fmu)
  !**********************************************************************
  !*     caluculate physical viscuss term  (in progress)                *
  !**********************************************************************
  use param_mod
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

  dux=0d0
  duy=0d0
  duz=0d0
  grad_uj = 0.d0
  grad_uk = 0.d0
  grad_ul = 0.d0

  !** loop start ****
  n = 1 
  dl(:)=0
  dl(n)=1

  !$OMP parallel do default(none) &
  !$OMP & firstprivate(dl,n,dux,duy,duz,grad_uj,grad_uk,grad_ul) &
  !$OMP & shared(vf,ql,qr,fmu,diff,lamb,jmax,kmax,lmax,ndmax) &
  !$OMP & shared(dx,Mach,Re,gami,mwi) &
  !$OMP & private(r1,r2,p1,p2,k1,k2,e1,e2,m1,m2,g1,g2,c21,c22) & 
  !$OMP & private(r_fl,kbar,j_fl_sum,duh,fmuave) & 
  !$OMP & private(ux,uy,uz,vx,vy,vz,wx,wy,wz,dch,txx,txy,tzx,tyy,tyz,tzz,visdiv,tau_u,qheat) & 
  !$OMP & private(ry1,ry2,h1,h2,u1,u2,u_fl,ry_fl,j_fl,hjbar,grad_y1,tau) 
  do l=1,lmax
    do k=1,kmax
      do j=1,jmax
        ! half of gradients (u,v,w,c2) of non-horizontal 
        ! dux(1) = qr(2,j+1,k,l,n) - ql(2,j-1,k,l,n)
        ! dux(2) = qr(3,j+1,k,l,n) - ql(3,j-1,k,l,n)
        ! dux(3) = qr(4,j+1,k,l,n) - ql(4,j-1,k,l,n)
        duy(1) = qr(2,j,k+1,l,n) - ql(2,j,k-1,l,n)
        duy(2) = qr(3,j,k+1,l,n) - ql(3,j,k-1,l,n)
        duy(3) = qr(4,j,k+1,l,n) - ql(4,j,k-1,l,n)
        duz(1) = qr(2,j,k,l+1,n) - ql(2,j,k,l-1,n)
        duz(2) = qr(3,j,k,l+1,n) - ql(3,j,k,l-1,n)
        duz(3) = qr(4,j,k,l+1,n) - ql(4,j,k,l-1,n)

        grad_uj(j,:,1) = 0.5d0*dux(:) * int(1-dl(1))
        grad_uj(j,:,2) = 0.5d0*duy(:) * int(1-dl(2))
        grad_uj(j,:,3) = 0.5d0*duz(:) * int(1-dl(3))
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

  !$OMP parallel do default(none) &
  !$OMP & firstprivate(dl,n,dux,duy,duz,grad_uj,grad_uk,grad_ul) &
  !$OMP & shared(vf,ql,qr,fmu,diff,lamb,jmax,kmax,lmax,ndmax) &
  !$OMP & shared(dx,Mach,Re,gami,mwi) &
  !$OMP & private(r1,r2,p1,p2,k1,k2,e1,e2,m1,m2,g1,g2,c21,c22) & 
  !$OMP & private(r_fl,kbar,j_fl_sum,duh,fmuave) & 
  !$OMP & private(ux,uy,uz,vx,vy,vz,wx,wy,wz,dch,txx,txy,tzx,tyy,tyz,tzz,visdiv,tau_u,qheat) & 
  !$OMP & private(ry1,ry2,h1,h2,u1,u2,u_fl,ry_fl,j_fl,hjbar,grad_y1,tau) 
  do l=1,lmax
    do j=1,jmax
      do k=0,kmax+1
        ! half of gradients (u,v,w,c2) of non-horizontal 
        dux(1) = qr(2,j+1,k,l,n) - ql(2,j-1,k,l,n)
        dux(2) = qr(3,j+1,k,l,n) - ql(3,j-1,k,l,n)
        dux(3) = qr(4,j+1,k,l,n) - ql(4,j-1,k,l,n)
        ! duy(1) = qr(2,j,k+1,l,n) - ql(2,j,k-1,l,n)
        ! duy(2) = qr(3,j,k+1,l,n) - ql(3,j,k-1,l,n)
        ! duy(3) = qr(4,j,k+1,l,n) - ql(4,j,k-1,l,n)
        duz(1) = qr(2,j,k,l+1,n) - ql(2,j,k,l-1,n)
        duz(2) = qr(3,j,k,l+1,n) - ql(3,j,k,l-1,n)
        duz(3) = qr(4,j,k,l+1,n) - ql(4,j,k,l-1,n)

        grad_uk(k,:,1) = 0.5d0*dux(:) * int(1-dl(1))
        grad_uk(k,:,2) = 0.5d0*duy(:) * int(1-dl(2))
        grad_uk(k,:,3) = 0.5d0*duz(:) * int(1-dl(3))
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

  !$OMP parallel do default(none) &
  !$OMP & firstprivate(dl,n,dux,duy,duz,grad_uj,grad_uk,grad_ul) &
  !$OMP & shared(vf,ql,qr,fmu,diff,lamb,jmax,kmax,lmax,ndmax) &
  !$OMP & shared(dx,Mach,Re,gami,mwi) &
  !$OMP & private(r1,r2,p1,p2,k1,k2,e1,e2,m1,m2,g1,g2,c21,c22) & 
  !$OMP & private(r_fl,kbar,j_fl_sum,duh,fmuave) & 
  !$OMP & private(ux,uy,uz,vx,vy,vz,wx,wy,wz,dch,txx,txy,tzx,tyy,tyz,tzz,visdiv,tau_u,qheat) & 
  !$OMP & private(ry1,ry2,h1,h2,u1,u2,u_fl,ry_fl,j_fl,hjbar,grad_y1,tau) 
  do k=1,kmax
    do j=1,jmax
      do l=0,lmax+1
        ! half of gradients (u,v,w,c2) of non-horizontal 
        dux(1) = qr(2,j+1,k,l,n) - ql(2,j-1,k,l,n)
        dux(2) = qr(3,j+1,k,l,n) - ql(3,j-1,k,l,n)
        dux(3) = qr(4,j+1,k,l,n) - ql(4,j-1,k,l,n)
        duy(1) = qr(2,j,k+1,l,n) - ql(2,j,k-1,l,n)
        duy(2) = qr(3,j,k+1,l,n) - ql(3,j,k-1,l,n)
        duy(3) = qr(4,j,k+1,l,n) - ql(4,j,k-1,l,n)
        ! duz(1) = qr(2,j,k,l+1,n) - ql(2,j,k,l-1,n)
        ! duz(2) = qr(3,j,k,l+1,n) - ql(3,j,k,l-1,n)
        ! duz(3) = qr(4,j,k,l+1,n) - ql(4,j,k,l-1,n)

        grad_ul(l,:,1) = 0.5d0*dux(:) * abs(1-dl(1))
        grad_ul(l,:,2) = 0.5d0*duy(:) * abs(1-dl(2))
        grad_ul(l,:,3) = 0.5d0*duz(:) * abs(1-dl(3))
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

        ! if(l==1 .and. k==2) then
        ! if(k==2) then
          ! write(*,fmt='(9E24.15e3)',advance='No') grad_ul(l,:,:)
          ! write(*,fmt='(8E24.15e3)',advance='No') vf(2:4,j,k,l,n)
          ! write(*,*)
        ! endif
      end do
    end do
  end do

  return
end subroutine visflux_dns

subroutine visflux_none(ql, qr, vf, fmu)
  !**********************************************************************
  !*     dummy subroutine for invicid flows                             *
  !**********************************************************************

  return
end subroutine visflux_none

subroutine calc_viscoefs(q, fmu)
  !**********************************************************************
  !*     caluculate the viscuss coefficients                            *
  !**********************************************************************
  use param_mod
  implicit none
  integer j,k,l
  double precision, dimension(0:(jmax+1),0:(kmax+1),0:(lmax+1)) :: fmu
  double precision,dimension(ndmax,0:(jmax+1),0:(kmax+1),0:(lmax+1)) :: q
  double precision c_1,c_2,gm

  !$omp parallel do shared(q,fmu) private(gm,c_1,c_2)
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
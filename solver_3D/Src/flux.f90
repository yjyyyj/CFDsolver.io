module flux
  implicit none
  
contains

  subroutine flux_div(ql, qr, f)
    !**********************************************************************
    !*     caluculate right-hand-side                                     *
    !**********************************************************************
    use param
    implicit none
    integer j,k,l,n
    double precision, dimension(ndmax,0:jmax,0:kmax,0:lmax,ndim) :: f
    double precision,dimension(ndmax,0:(jmax+1),0:(kmax+1),0:(lmax+1),ndim) :: ql, qr
    double precision r1,r2,p1,p2,k1,k2,e1,e2
    double precision m1,m2,g1,g2
    double precision,dimension(nspecies) :: ry1, ry2
    double precision,dimension(ndim) :: u1,u2
    integer,dimension(ndim) :: dl

    do n=1,ndim
      dl(:)=0
      dl(n)=1
      
      do l=0,lmax
        do k=0,kmax
          do j=0,jmax
          !*** x direction *************************************!
            ! left face q
            u1(1)  = ql(2,j,k,l,n)
            u1(2)  = ql(3,j,k,l,n)
            u1(3)  = ql(4,j,k,l,n)
            p1  = ql(5,j,k,l,n)
            ry1(:) = ql(6:ndmax,j,k,l,n)
            r1 = sum(ry1(:))

            m1 = calmw(ql(:,j,k,l,n))
            g1 = calgm(ql(:,j,k,l,n),m1)
        
            ! right face q
            u2(1)  = qr(2,j+dl(1),k+dl(2),l+dl(3),n)
            u2(2)  = qr(3,j+dl(1),k+dl(2),l+dl(3),n)
            u2(3)  = qr(4,j+dl(1),k+dl(2),l+dl(3),n)
            p2  = qr(5,j+dl(1),k+dl(2),l+dl(3),n)
            ry2(:) = qr(6:ndmax,j+dl(1),k+dl(2),l+dl(3),n)
            r2 = sum(ry2(:))

            m2 = calmw(qr(:,j+dl(1),k+dl(2),l+dl(3),1))
            g2 = calgm(qr(:,j+dl(1),k+dl(2),l+dl(3),1),m2)  

            k1 = r1*(u1(1)*u1(1)+u1(2)*u1(2)+u1(3)*u1(3))*0.5d0
            k2 = r2*(u2(1)*u2(1)+u2(2)*u2(2)+u2(3)*u2(3))*0.5d0

            e1 = p1*g1
            e2 = p2*g2
            
            f(1,j,k,l,n) = 0.d0
            f(2,j,k,l,n) = 0.5d0*( r1*u1(1)*u1(n) + r2*u2(1)*u2(n) + (p1 + p2)*dl(1) )
            f(3,j,k,l,n) = 0.5d0*( r1*u1(2)*u1(n) + r2*u2(2)*u2(n) + (p1 + p2)*dl(2) )
            f(4,j,k,l,n) = 0.5d0*( r1*u1(3)*u1(n) + r2*u2(3)*u2(n) + (p1 + p2)*dl(3) )
            f(5,j,k,l,n) = 0.5d0*( (k1+e1+p1)*u1(n) + (k2+e2+p2)*u2(n) )
            f(6:ndmax,j,k,l,n) = 0.5d0*( ry1(:)*u1(n) + ry2(:)*u2(n) )
          end do
        end do
      end do
    end do

    return
  end subroutine flux_div


  subroutine flux_upwind(ql, qr, f)
    !**********************************************************************
    !*     caluculate right-hand-side                                     *
    !**********************************************************************
    use param
    implicit none
    integer j,k,l,n
    double precision, dimension(ndmax,0:jmax,0:kmax,0:lmax,ndim) :: f
    double precision,dimension(ndmax,0:(jmax+1),0:(kmax+1),0:(lmax+1),ndim) :: ql, qr
    double precision r1,r2,p1,p2,k1,k2,e1,e2
    double precision m1,m2,g1,g2,c1,c2
    double precision A
    double precision,dimension(nspecies) :: ry1, ry2
    double precision,dimension(ndim) :: u1,u2
    integer,dimension(ndim) :: dl

    do n=1,ndim
      dl(:)=0
      dl(n)=1
      do l=0,lmax
        do k=0,kmax
          do j=0,jmax
          !*** x direction *************************************!
            ! left face q
            u1(1)  = ql(2,j,k,l,n)
            u1(2)  = ql(3,j,k,l,n)
            u1(3)  = ql(4,j,k,l,n)
            p1  = ql(5,j,k,l,n)
            ry1(:) = ql(6:ndmax,j,k,l,n)
            r1 = sum(ry1(:))

            m1 = calmw(ql(:,j,k,l,n))
            g1 = calgm(ql(:,j,k,l,n),m1)
        
            ! right face q
            u2(1)  = qr(2,j+dl(1),k+dl(2),l+dl(3),n)
            u2(2)  = qr(3,j+dl(1),k+dl(2),l+dl(3),n)
            u2(3)  = qr(4,j+dl(1),k+dl(2),l+dl(3),n)
            p2  = qr(5,j+dl(1),k+dl(2),l+dl(3),n)
            ry2(:) = qr(6:ndmax,j+dl(1),k+dl(2),l+dl(3),n)
            r2 = sum(ry2(:))

            m2 = calmw(qr(:,j+dl(1),k+dl(2),l+dl(3),n))
            g2 = calgm(qr(:,j+dl(1),k+dl(2),l+dl(3),n),m2)  

            k1 = r1*(u1(1)*u1(1)+u1(2)*u1(2)+u1(3)*u1(3))*0.5d0
            k2 = r2*(u2(1)*u2(1)+u2(2)*u2(2)+u2(3)*u2(3))*0.5d0
            e1 = p1*g1
            e2 = p2*g2
            c1 = sqrt(p1/r1*(1d0/g1+1d0))
            c2 = sqrt(p2/r2*(1d0/g2+1d0))

            A = max(abs(u1(n))+c1, abs(u2(n))+c2)

            f(1,j,k,l,n) = 0.d0
            f(2,j,k,l,n) = 0.5d0*( r1*u1(1)*u1(n) + r2*u2(1)*u2(n) + (p1 + p2)*dl(1) ) - 0.5d0*A*( r2*u2(1) - r1*u1(1))
            f(3,j,k,l,n) = 0.5d0*( r1*u1(2)*u1(n) + r2*u2(2)*u2(n) + (p1 + p2)*dl(2) ) - 0.5d0*A*( r2*u2(2) - r1*u1(2))
            f(4,j,k,l,n) = 0.5d0*( r1*u1(3)*u1(n) + r2*u2(3)*u2(n) + (p1 + p2)*dl(3) ) - 0.5d0*A*( r2*u2(3) - r1*u1(3))
            f(5,j,k,l,n) = 0.5d0*( (k1+e1+p1)*u1(n) + (k2+e2+p2)*u2(n) ) - 0.5d0*A*( k2+e2 -k1-e1)
            f(6:ndmax,j,k,l,n) = 0.5d0*( ry1(:)*u1(n) + ry2(:)*u2(n) ) - 0.5d0*A*(ry2(:) - ry1(:))
          end do
        end do
      end do
    end do

    return
  end subroutine flux_upwind


  subroutine flux_roe(ql, qr, f)
    !**********************************************************************
    !*     caluculate right-hand-side                                     *
    !**********************************************************************
    use param
    implicit none
    integer j,k,l,n
    double precision, dimension(ndmax,0:jmax,0:kmax,0:lmax,ndim) :: f
    double precision,dimension(ndmax,0:(jmax+1),0:(kmax+1),0:(lmax+1),ndim) :: ql, qr
    double precision r1,r2,p1,p2,k1,k2,e1,e2
    double precision m1,m2,g1,g2,c1,c2
    double precision,dimension(nspecies) :: ry1, ry2
    double precision,dimension(ndim) :: u1,u2,ut
    integer,dimension(ndim) :: dl
    double precision h1, h2, uvw2
    double precision ct,ht,gt
    double precision,dimension(ndmax,ndmax) :: A, R, RInv
    double precision,dimension(ndmax) :: tmp

    do n=1,ndim
      dl(:)=0
      dl(n)=1
      do l=0,lmax
        do k=0,kmax
          do j=0,jmax
          !*** x direction *************************************!
            ! left face q
            u1(1)  = ql(2,j,k,l,n)
            u1(2)  = ql(3,j,k,l,n)
            u1(3)  = ql(4,j,k,l,n)
            p1  = ql(5,j,k,l,n)
            ry1(:) = ql(6:ndmax,j,k,l,n)
            r1 = sum(ry1(:))

            m1 = calmw(ql(:,j,k,l,n))
            g1 = calgm(ql(:,j,k,l,n),m1)
        
            ! right face q
            u2(1)  = qr(2,j+dl(1),k+dl(2),l+dl(3),n)
            u2(2)  = qr(3,j+dl(1),k+dl(2),l+dl(3),n)
            u2(3)  = qr(4,j+dl(1),k+dl(2),l+dl(3),n)
            p2  = qr(5,j+dl(1),k+dl(2),l+dl(3),n)
            ry2(:) = qr(6:ndmax,j+dl(1),k+dl(2),l+dl(3),n)
            r2 = sum(ry2(:))

            m2 = calmw(qr(:,j+dl(1),k+dl(2),l+dl(3),n))
            g2 = calgm(qr(:,j+dl(1),k+dl(2),l+dl(3),n),m2)  

            k1 = r1*(u1(1)*u1(1)+u1(2)*u1(2)+u1(3)*u1(3))*0.5d0
            k2 = r2*(u2(1)*u2(1)+u2(2)*u2(2)+u2(3)*u2(3))*0.5d0
            e1 = p1*g1
            e2 = p2*g2
            h1 = (k1+e1+p1)/r1
            h2 = (k2+e2+p2)/r2

            c1 = sqrt(p1/r1*(1d0/g1+1d0))
            c2 = sqrt(p2/r2*(1d0/g2+1d0))

            gt = 0.5d0*(g1+g2)

            ! **Roe-averaged q*************************************************
            ut(:) = (sqrt(r2)*u2(:) + sqrt(r1)*u1(:))/(sqrt(r2) + sqrt(r1))
            uvw2  = ut(1)*ut(1) + ut(2)*ut(2) + ut(3)*ut(3)
            ht = (sqrt(r2)*h2 + sqrt(r1)*h1)/(sqrt(r2) + sqrt(r1))
            ct = sqrt((ht - 0.5d0*uvw2)/gt)
            !******************************************************************

            ! Lambda matrix
            A(1,1:3) = (/ abs(ut(n)-ct),0d0,0d0 /)
            A(2,1:3) = (/ 0d0,abs(ut(n)),0d0 /)
            A(3,1:3) = (/ 0d0,0d0,abs(ut(n)+ct) /)

            ! R matrix
            R(1,1:3) = (/ 1d0,1d0,1d0 /)
            R(2,1:3) = (/ ut(n)-ct,ut(n),ut(n)+ct /)
            R(3,1:3) = (/ ht-ct*ut,0.5d0*ut*ut,ht+ct*ut /)

            ! R^-1 matrix
            RInv(1,1:3) = (/ 0.5d0*ut*ut+ct*ut*gt,-ut-ct*gt,1d0 /)
            RInv(2,1:3) = (/ 2d0*(ht-ut*ut),2d0*ut,-2d0 /)
            RInv(3,1:3) = (/ 0.5d0*ut*ut-ct*ut*gt,-ut+ct*gt,1d0 /)

            RInv = RInv*0.5d0/gt/ct/ct

            tmp =(/ r2-r1, r2*u2-r1*u1, (k2+e2)-(k1+e1) /)
            tmp = matmul( RInv,tmp )
            tmp = matmul( A,tmp )
            tmp = matmul( R,tmp )

            ! flux
            f(1,j,k,l,n) = 0.5d0*( r1*u1(n) + r2*u2(n) - tmp(1) )
            f(2,j,k,l,n) = 0.5d0*( r1*u1(1)*u1(n) + r2*u2(1)*u2(n) + (p1 + p2)*dl(1) ) - tmp(2)
            f(3,j,k,l,n) = 0.5d0*( r1*u1(2)*u1(n) + r2*u2(2)*u2(n) + (p1 + p2)*dl(2) ) - tmp(3)
            f(4,j,k,l,n) = 0.5d0*( r1*u1(3)*u1(n) + r2*u2(3)*u2(n) + (p1 + p2)*dl(3) ) - tmp(4)
            f(5,j,k,l,n) = 0.5d0*( (k1+e1+p1)*u1(n) + (k2+e2+p2)*u2(n) ) - tmp(5)
            f(6:ndmax,j,k,l,n) = 0.5d0*( ry1(:)*u1(n) + ry2(:)*u2(n) ) - tmp(6:ndmax)
          end do
        end do
      end do
    end do

    return
  end subroutine flux_roe


  subroutine flux_slau(ql, qr, f)
    !**********************************************************************
    !*     caluculate right-hand-side                                     *
    !**********************************************************************
    use param
    implicit none
    integer j,k,l,n
    double precision, dimension(ndmax,0:jmax,0:kmax,0:lmax,ndim) :: f
    double precision,dimension(ndmax,0:(jmax+1),0:(kmax+1),0:(lmax+1),ndim) :: ql, qr
    double precision r1,r2,p1,p2,k1,k2,e1,e2
    double precision m1,m2,g1,g2,c1,c2
    double precision,dimension(nspecies) :: ry1, ry2
    double precision,dimension(ndim) :: u1,u2
    integer,dimension(ndim) :: dl
    double precision h1, h2
    double precision,dimension(ndim) :: ulm1,ulm2
    double precision xmach1,xmach2,zz,uvw21,uvw22
    double precision cbv,xm1,xm2,temp,xmh,chi,sw1,sw2,bt1,bt2,pt
    double precision g,unb,un_p,un_m,fm,fmp,fmm
    integer :: nthon

    nthon = 1

    do n=1,ndim
      dl(:)=0
      dl(n)=1
      do l=0,lmax
        do k=0,kmax
          do j=0,jmax
          !*** x direction *************************************!
            ! left face q
            u1(1)  = ql(2,j,k,l,n)
            u1(2)  = ql(3,j,k,l,n)
            u1(3)  = ql(4,j,k,l,n)
            p1     = ql(5,j,k,l,n)
            ry1(:) = ql(6:ndmax,j,k,l,n)
            r1 = sum(ry1(:))

            m1 = calmw(ql(:,j,k,l,n))
            g1 = calgm(ql(:,j,k,l,n),m1)
        
            ! right face q
            u2(1)  = qr(2,j+dl(1),k+dl(2),l+dl(3),n)
            u2(2)  = qr(3,j+dl(1),k+dl(2),l+dl(3),n)
            u2(3)  = qr(4,j+dl(1),k+dl(2),l+dl(3),n)
            p2     = qr(5,j+dl(1),k+dl(2),l+dl(3),n)
            ry2(:) = qr(6:ndmax,j+dl(1),k+dl(2),l+dl(3),n)
            r2 = sum(ry2(:))

            m2 = calmw(qr(:,j+dl(1),k+dl(2),l+dl(3),n))
            g2 = calgm(qr(:,j+dl(1),k+dl(2),l+dl(3),n),m2)  

            uvw21 = u1(1)*u1(1)+u1(2)*u1(2)+u1(3)*u1(3)
            uvw22 = u2(1)*u2(1)+u2(2)*u2(2)+u2(3)*u2(3)
            k1 = r1*uvw21*0.5d0
            k2 = r2*uvw22*0.5d0
            e1 = p1*g1
            e2 = p2*g2
            h1 = (k1+e1+p1)/r1
            h2 = (k2+e2+p2)/r2

            c1 = sqrt(p1/r1*(1d0/g1+1d0))
            c2 = sqrt(p2/r2*(1d0/g2+1d0))

            !***Thornber's coorection for low Mach number flows (i)
            if(nthon.eq.1)then
              xmach1 = sqrt(uvw21/(c1*c1))
              xmach2 = sqrt(uvw22/(c2*c2))
              zz     = min(1.0d0,max(xmach1,xmach2))
      
              ulm1(:) = 0.5d0*(u1(:)+u2(:))+zz*0.5d0*(u1(:)-u2(:))
              ulm2(:) = 0.5d0*(u1(:)+u2(:))-zz*0.5d0*(u1(:)-u2(:))
      
              u1(:)   = ulm1(:)    
              u2(:)   = ulm2(:)
            end if
            
            !*** c bar M+ M- Original SLAU
            cbv  = 2.0d0/(c1+c2)
            xm1  = u1(n)*cbv
            xm2  = u2(n)*cbv
            temp = 0.5d0*(uvw21 + uvw22)
            
            !*** M hat
            xmh  = min(1.0d0,sqrt(temp)*cbv)
            chi  = (1.0d0-xmh)*(1.0d0-xmh)
            
            !**** beta +-
            sw1  = max(0.0d0,sign(1.0d0,abs(xm1)-1.0d0))
            sw2  = max(0.0d0,sign(1.0d0,abs(xm2)-1.0d0))
        
            bt1  = (1.0d0-sw1)*0.25d0*(2.0d0-xm1)*((xm1+1.0d0)*(xm1+1.0d0)) + sw1*0.50d0*(1.0d0+sign(1.0d0,xm1))
            bt2  = (1.0d0-sw2)*0.25d0*(2.0d0+xm2)*((xm2-1.0d0)*(xm2-1.0d0)) + sw2*0.50d0*(1.0d0-sign(1.0d0,xm2))
        
            pt  = 0.5d0*( (p1+p2)+(bt1-bt2)*(p1-p2) + (p1+p2)*(1.0d0-chi)*(bt1 + bt2-1.0d0))
            
            !**** set xm_dot(fm)
            g    = -max(min(xm1,0.0d0),-1.0d0)*min(max(xm2,0.0d0),1.0d0)
            unb  = (r1*abs(u1(n)) + r2*abs(u2(n)))/(r1+r2)
            un_p = (1.0d0-g)*unb + g*abs(u1(n))
            un_m = (1.0d0-g)*unb + g*abs(u2(n))
            fm   = 0.5d0*(r1*(u1(n) + un_p) + r2*(u2(n) - un_m) - chi*(p2 - p1)*cbv  )
            fmp  = 0.5d0*( fm + abs(fm) )
            fmm  = 0.5d0*( fm - abs(fm) )
            
            !**** numerical flux (p:plus m:minus)
            ! flux
            f(1,j,k,l,n) = 0d0
            f(2,j,k,l,n) = (fmp*u1(1) + fmm*u2(1) + pt*dl(1))
            f(3,j,k,l,n) = (fmp*u1(2) + fmm*u2(2) + pt*dl(2))
            f(4,j,k,l,n) = (fmp*u1(3) + fmm*u2(3) + pt*dl(3))
            f(5,j,k,l,n) = (fmp*h1    + fmm*h2          )
            f(6,j,k,l,n) = ( fmp       + fmm            )
          end do
        end do
      end do
    end do

    return
  end subroutine flux_slau


  subroutine flux_KEEP(ql, qr, f)
    !**********************************************************************
    !*     caluculate right-hand-side                                     *
    !**********************************************************************
    use param
    implicit none
    integer j,k,l,n
    double precision, dimension(ndmax,0:jmax,0:kmax,0:lmax,ndim) :: f
    double precision,dimension(ndmax,0:(jmax+1),0:(kmax+1),0:(lmax+1),ndim) :: ql, qr
    double precision r1,r2,p1,p2
    double precision m1,m2,g1,g2,c1,c2
    double precision r_fl,cbar,pibar,kbar,ibar,pbar
    double precision,dimension(nspecies) :: ry1, ry2, ry_fl, cybar
    double precision,dimension(ndim) :: u1,u2,u_fl,mubar
    integer,dimension(ndim) :: dl

    do n=1,ndim
      dl(:)=0
      dl(n)=1
      do l=0,lmax
        do k=0,kmax
          do j=0,jmax
          !*** x direction *************************************!
            ! left face q
            u1(1)  = ql(2,j,k,l,n)
            u1(2)  = ql(3,j,k,l,n)
            u1(3)  = ql(4,j,k,l,n)
            p1  = ql(5,j,k,l,n)
            ry1(:) = ql(6:ndmax,j,k,l,n)
            r1 = sum(ry1(:))

            m1 = calmw(ql(:,j,k,l,n))
            g1 = calgm(ql(:,j,k,l,n),m1)
        
            ! right face q
            u2(1)  = qr(2,j+dl(1),k+dl(2),l+dl(3),n)
            u2(2)  = qr(3,j+dl(1),k+dl(2),l+dl(3),n)
            u2(3)  = qr(4,j+dl(1),k+dl(2),l+dl(3),n)
            p2  = qr(5,j+dl(1),k+dl(2),l+dl(3),n)
            ry2(:) = qr(6:ndmax,j+dl(1),k+dl(2),l+dl(3),n)
            r2 = sum(ry2(:))

            m2 = calmw(qr(:,j+dl(1),k+dl(2),l+dl(3),n))
            g2 = calgm(qr(:,j+dl(1),k+dl(2),l+dl(3),n),m2)  

            c1 = m1/r1*r2/m2
            c2 = m2/r2*r1/m1 

            ry_fl(:) = 0.5d0*( ry2(:) + ry1(:))
            r_fl = sum(ry_fl(:))

            u_fl(:) = 0.5d0*( u1(:) + u2(:) )
            pibar = 0.5d0*( p1 + p2 )

            cybar(:) = ry_fl(:)*u_fl(n)
            cbar = r_fl*u_fl(n)
            mubar(:) = cbar*u_fl(:)
      
            kbar = cbar*(u1(1)*u2(1)+u1(2)*u2(2)+u1(3)*u2(3))*0.5d0
            ibar = 0.5d0*(p1*g1 + p2*g2)*u_fl(n)  !KEEP_PE
            ! ibar = cbar*0.5d0*(g1+g2)*0.5d0*(p1/r1+p2/r2)   !KEEP
            pbar = ( p1*u2(n) + u1(n)*p2 )*0.5d0
        
            f(1,j,k,l,n) = 0.d0
            f(2,j,k,l,n) = mubar(1) + pibar*dl(1)
            f(3,j,k,l,n) = mubar(2) + pibar*dl(2)
            f(4,j,k,l,n) = mubar(3) + pibar*dl(3)
            f(5,j,k,l,n) = kbar + ibar + pbar
            f(6:ndmax,j,k,l,n) = cybar(:)
          end do
        end do
      end do
    end do

    return
  end subroutine flux_KEEP


  subroutine flux_KEEP_PE(ql, qr, f)
    !**********************************************************************
    !*     caluculate right-hand-side                                     *
    !**********************************************************************
    use param
    implicit none
    integer j,k,l,n
    double precision, dimension(ndmax,0:jmax,0:kmax,0:lmax,ndim) :: f
    double precision,dimension(ndmax,0:(jmax+1),0:(kmax+1),0:(lmax+1),ndim) :: ql, qr
    double precision r1,r2,p1,p2
    double precision m1,m2,g1,g2,c1,c2
    double precision r_fl,cbar,pibar,kbar,ibar,pbar
    double precision,dimension(nspecies) :: ry1, ry2, ry_fl, cybar
    double precision,dimension(ndim) :: u1,u2,u_fl,mubar
    integer,dimension(ndim) :: dl

    do n=1,ndim
      dl(:)=0
      dl(n)=1
      do l=0,lmax
        do k=0,kmax
          do j=0,jmax
          !*** x direction *************************************!
            ! left face q
            u1(1)  = ql(2,j,k,l,n)
            u1(2)  = ql(3,j,k,l,n)
            u1(3)  = ql(4,j,k,l,n)
            p1  = ql(5,j,k,l,n)
            ry1(:) = ql(6:ndmax,j,k,l,n)
            r1 = sum(ry1(:))

            m1 = calmw(ql(:,j,k,l,n))
            g1 = calgm(ql(:,j,k,l,n),m1)
        
            ! right face q
            u2(1)  = qr(2,j+dl(1),k+dl(2),l+dl(3),n)
            u2(2)  = qr(3,j+dl(1),k+dl(2),l+dl(3),n)
            u2(3)  = qr(4,j+dl(1),k+dl(2),l+dl(3),n)
            p2  = qr(5,j+dl(1),k+dl(2),l+dl(3),n)
            ry2(:) = qr(6:ndmax,j+dl(1),k+dl(2),l+dl(3),n)
            r2 = sum(ry2(:))

            m2 = calmw(qr(:,j+dl(1),k+dl(2),l+dl(3),n))
            g2 = calgm(qr(:,j+dl(1),k+dl(2),l+dl(3),n),m2)  

            c1 = m1/r1*r2/m2
            c2 = m2/r2*r1/m1 

            ry_fl(:) = 0.5d0*( ry2(:) + ry1(:))
            r_fl = sum(ry_fl(:))

            u_fl(:) = 0.5d0*( u1(:) + u2(:) )
            pibar = 0.5d0*( p1 + p2 )

            cybar(:) = ry_fl(:)*u_fl(n)
            cbar = r_fl*u_fl(n)
            mubar(:) = cbar*u_fl(:)
      
            kbar = cbar*(u1(1)*u2(1)+u1(2)*u2(2)+u1(3)*u2(3))*0.5d0
            ibar = 0.5d0*(p1*g1 + p2*g2)*u_fl(n)  !KEEP_PE
            ! ibar = cbar*0.5d0*(g1+g2)*0.5d0*(p1/r1+p2/r2)   !KEEP
            pbar = ( p1*u2(n) + u1(n)*p2 )*0.5d0
        
            f(1,j,k,l,n) = 0.d0
            f(2,j,k,l,n) = mubar(1) + pibar*dl(1)
            f(3,j,k,l,n) = mubar(2) + pibar*dl(2)
            f(4,j,k,l,n) = mubar(3) + pibar*dl(3)
            f(5,j,k,l,n) = kbar + ibar + pbar
            f(6:ndmax,j,k,l,n) = cybar(:)
          end do
        end do
      end do
    end do

    return
  end subroutine flux_KEEP_PE


  subroutine flux_proposed(ql, qr, f)
    !**********************************************************************
    !*     caluculate right-hand-side                                     *
    !**********************************************************************
    use param
    implicit none
    integer j,k,l,n
    double precision, dimension(ndmax,0:jmax,0:kmax,0:lmax,ndim) :: f
    double precision,dimension(ndmax,0:(jmax+1),0:(kmax+1),0:(lmax+1),ndim) :: ql, qr
    double precision r1,r2,p1,p2
    double precision m1,m2,g1,g2,c1,c2
    double precision r_fl,cbar,pibar,kbar,ibar,pbar
    double precision,dimension(nspecies) :: ry1, ry2, ry_fl, cybar
    double precision,dimension(ndim) :: u1,u2,u_fl,mubar
    integer,dimension(ndim) :: dl

    do n=1,ndim
      dl(:)=0
      dl(n)=1
      do l=0,lmax
        do k=0,kmax
          do j=0,jmax
          !*** x direction *************************************!
            ! left face q
            u1(1)  = ql(2,j,k,l,n)
            u1(2)  = ql(3,j,k,l,n)
            u1(3)  = ql(4,j,k,l,n)
            p1  = ql(5,j,k,l,n)
            ry1(:) = ql(6:ndmax,j,k,l,n)
            r1 = sum(ry1(:))

            m1 = calmw(ql(:,j,k,l,n))
            g1 = calgm(ql(:,j,k,l,n),m1)
        
            ! right face q
            u2(1)  = qr(2,j+dl(1),k+dl(2),l+dl(3),n)
            u2(2)  = qr(3,j+dl(1),k+dl(2),l+dl(3),n)
            u2(3)  = qr(4,j+dl(1),k+dl(2),l+dl(3),n)
            p2  = qr(5,j+dl(1),k+dl(2),l+dl(3),n)
            ry2(:) = qr(6:ndmax,j+dl(1),k+dl(2),l+dl(3),n)
            r2 = sum(ry2(:))

            m2 = calmw(qr(:,j+dl(1),k+dl(2),l+dl(3),n))
            g2 = calgm(qr(:,j+dl(1),k+dl(2),l+dl(3),n),m2)  

            c1 = m1/r1*r2/m2
            c2 = m2/r2*r1/m1 

            ry_fl(:) = 0.5d0*( c2*ry2(:) + c1*ry1(:))
            r_fl = sum(ry_fl(:))

            u_fl(:) = 0.5d0*( u1(:) + u2(:) )
            pibar = 0.5d0*( p1 + p2 )

            cybar(:) = ry_fl(:)*u_fl(n)
            cbar = r_fl*u_fl(n)
            mubar(:) = cbar*u_fl(:)
      
            kbar = cbar*(u1(1)*u2(1)+u1(2)*u2(2)+u1(3)*u2(3))*0.5d0
            ibar = 0.5d0*(p1*g1 + p2*g2)*u_fl(n)
            pbar = ( p1*u2(n) + u1(n)*p2 )*0.5d0
        
            f(1,j,k,l,n) = 0.d0
            f(2,j,k,l,n) = mubar(1) + pibar*dl(1)
            f(3,j,k,l,n) = mubar(2) + pibar*dl(2)
            f(4,j,k,l,n) = mubar(3) + pibar*dl(3)
            f(5,j,k,l,n) = kbar + ibar + pbar
            f(6:ndmax,j,k,l,n) = cybar(:)
          end do
        end do
      end do
    end do

    return
  end subroutine flux_proposed
  
end module flux

subroutine step_euler(q, qc, myscheme)
  !**********************************************************************
  !*     caluculate right-hand-side                                     *
  !**********************************************************************
  use param
  use flux
  use scheme_mod
  implicit none
  type(scheme) :: myscheme

  integer j,k,l
  double precision,dimension(ndmax,0:(jmax+1),0:(kmax+1),0:(lmax+1)) :: q       ! primitive q
  double precision,dimension(ndmax,0:(jmax+1),0:(kmax+1),0:(lmax+1),ndim) :: ql,qr   ! primitive dq
  double precision,dimension(ndmax,0:jmax,0:kmax,0:lmax,ndim) :: f       ! flux
  double precision,dimension(ndmax,0:jmax,0:kmax,0:lmax,ndim) :: vf       ! flux
  double precision,dimension(ndmax,jmax,kmax,lmax) :: s       ! rhs 
  double precision,dimension(ndmax,jmax,kmax,lmax) :: qc      ! renew q
  double precision,dimension(0:(jmax+1),0:(kmax+1),0:(lmax+1)) :: fmu           ! viscous coefficient

  f  = 0.0d0
  vf = 0.0d0
  ql = 0.0d0
  qr = 0.0d0
  s  = 0.0d0
  fmu = 0.0d0
  resid = 0.0d0

  call calc_viscoefs(q, fmu)
  call bc(q, fmu)
  call myscheme%calc_faceQ(q, ql, qr)
  call myscheme%calc_flux(ql, qr, f)
  call visflux_dns(ql, qr, vf, fmu)
  ! call visflux(ql, qr, vf)
  
  !$omp parallel do shared(q,qc,f,vf) firstprivate(dt,dx)
  do l=1,lmax
    do k=1,kmax
      do j=1,jmax
        s(1:ndmax,j,k,l)  =  -dt/dx(1) *(  f(1:ndmax,j,k,l,1) - f(1:ndmax,j-1,k,l,1) )   &
                          &  -dt/dx(2) *(  f(1:ndmax,j,k,l,2) - f(1:ndmax,j,k-1,l,2) )   &
                          &  -dt/dx(3) *(  f(1:ndmax,j,k,l,3) - f(1:ndmax,j,k,l-1,3) )   & ! f(j): ~f(j+1/2)
                          &  +dt/dx(1) *(  vf(1:ndmax,j,k,l,1) - vf(1:ndmax,j-1,k,l,1) )   & 
                          &  +dt/dx(2) *(  vf(1:ndmax,j,k,l,2) - vf(1:ndmax,j,k-1,l,2) )   & 
                          &  +dt/dx(3) *(  vf(1:ndmax,j,k,l,3) - vf(1:ndmax,j,k,l-1,3) )
        
        ! time advancement
        ! 1st Euler *************************************************************
        qc(1:ndmax,j,k,l) = qc(1:ndmax,j,k,l) + s(1:ndmax,j,k,l)

        ! renew primitive q *************************************************************
        qc(1,j,k,l) = sum(qc(6:ndmax,j,k,l))
        mw(j,k,l)  = calmw(qc(:,j,k,l))
        gam(j,k,l) = calgm(qc(:,j,k,l),mw(j,k,l))

        q(1,j,k,l) = qc(1,j,k,l) ! rho
        q(2,j,k,l) = qc(2,j,k,l)/qc(1,j,k,l) ! u
        q(3,j,k,l) = qc(3,j,k,l)/qc(1,j,k,l) ! v
        q(4,j,k,l) = qc(4,j,k,l)/qc(1,j,k,l) ! w
        q(5,j,k,l) = ( qc(5,j,k,l) - 0.5d0*( qc(2,j,k,l)**2 +qc(3,j,k,l)**2 +qc(4,j,k,l)**2 )/qc(1,j,k,l) )/gam(j,k,l) ! p
        q(6:ndmax,j,k,l) = qc(6:ndmax,j,k,l) ! rhoy

        !*****************************************************
      enddo
    enddo
  enddo

  !$omp parallel do reduction(+:resid)  
  do l=1,lmax
    do k=1,kmax
      do j=1,jmax
        ! compute l2 residual of rhs 
        resid(1) = resid(1) + s(1,j,k,l)*s(1,j,k,l) 
        resid(2) = resid(2) + s(2,j,k,l)*s(2,j,k,l) 
        resid(3) = resid(3) + s(3,j,k,l)*s(3,j,k,l) 
        resid(4) = resid(4) + s(4,j,k,l)*s(4,j,k,l) 
        resid(5) = resid(5) + s(5,j,k,l)*s(5,j,k,l)                       
        resid(6:ndmax) = resid(6:ndmax) + s(6:ndmax,j,k,l)*s(6:ndmax,j,k,l)
      enddo
    enddo
  enddo

  return
end subroutine step_euler

subroutine step_RK4(q, qc, myscheme)
  !**********************************************************************
  !*     caluculate right-hand-side                                     *
  !**********************************************************************
  use param
  use flux
  use scheme_mod
  implicit none
  type(scheme) :: myscheme

  integer j,k,l, nn
  double precision,dimension(ndmax,0:(jmax+1),0:(kmax+1),0:(lmax+1)) :: q       ! primitive q
  double precision,dimension(ndmax,0:(jmax+1),0:(kmax+1),0:(lmax+1),ndim) :: ql,qr   ! primitive dq
  double precision,dimension(ndmax,0:jmax,0:kmax,0:lmax,ndim) :: f       ! flux
  double precision,dimension(ndmax,0:jmax,0:kmax,0:lmax,ndim) :: vf       ! viscus flux
  double precision,dimension(ndmax,jmax,kmax,lmax) :: s       ! rhs 
  double precision,dimension(ndmax,jmax,kmax,lmax) :: qc      ! renew q
  double precision,dimension(ndmax,jmax,kmax,lmax) :: qcold1       ! primitive q
  double precision,dimension(ndmax,jmax,kmax,lmax) :: qcold2       ! primitive q

  double precision,parameter :: rc1(1:4) = (/1.0d0,1.0d0,1.0d0,1.0d0/)
  double precision,parameter :: rc2(1:4) = (/0.0d0,0.0d0,0.0d0,1.0d0/6.0d0/)
  double precision,parameter :: rc3(1:4) = (/1.0d0/2.0d0,1.0d0/2.0d0,1.0d0,1.0d0/6.0d0/)
  double precision,parameter :: rc4(1:4) = (/1.0d0,2.0d0,2.0d0,0.0d0/)
  double precision,dimension(0:(jmax+1),0:(kmax+1),0:(lmax+1)) :: fmu           ! viscous coefficient

  qcold1 = qc
  qcold2 = 0.d0

  do nn = 1,4

    f  = 0.0d0
    vf = 0.0d0
    ql = 0.0d0
    qr = 0.0d0
    s  = 0.0d0
    fmu = 0.0d0
    resid = 0.0d0


    call calc_viscoefs(q, fmu)
    call bc(q,fmu)
    call myscheme%calc_faceQ(q, ql, qr)
    call myscheme%calc_flux(ql, qr, f)
    call visflux_dns(ql, qr, vf,fmu)
    ! ここで落ちてる↑
    ! call dflux(ql, qr, vf)

    !$omp parallel do shared(q,qc,qcold1,qcold2,f,vf)
    do l=1,lmax
      do k=1,kmax
        do j=1,jmax
          s(1:ndmax,j,k,l)  = -dt/dx(1) *(  f(1:ndmax,j,k,l,1) - f(1:ndmax,j-1,k,l,1) )   &
                            & -dt/dx(2) *(  f(1:ndmax,j,k,l,2) - f(1:ndmax,j,k-1,l,2) )   &
                            & -dt/dx(3) *(  f(1:ndmax,j,k,l,3) - f(1:ndmax,j,k,l-1,3) )   & ! f(j): ~f(j+1/2)
                            & +dt/dx(1) *(  vf(1:ndmax,j,k,l,1) - vf(1:ndmax,j-1,k,l,1) )   & 
                            & +dt/dx(2) *(  vf(1:ndmax,j,k,l,2) - vf(1:ndmax,j,k-1,l,2) )   & 
                            & +dt/dx(3) *(  vf(1:ndmax,j,k,l,3) - vf(1:ndmax,j,k,l-1,3) )  

          ! RK4 *************************************************************
          qc(:,j,k,l) = rc1(nn)* qcold1(:,j,k,l)    &
          &           + rc2(nn)* qcold2(:,j,k,l)    &
          &           + rc3(nn)*      s(:,j,k,l)

          qcold2(:,j,k,l) = qcold2(:,j,k,l) + rc4(nn)*s(:,j,k,l) 
          ! ******************************************************************

          ! renew primitive q *********************************
          qc(1,j,k,l) = sum(qc(6:ndmax,j,k,l))
          mw(j,k,l)  = calmw(qc(:,j,k,l))
          gam(j,k,l) = calgm(qc(:,j,k,l),mw(j,k,l))

          q(1,j,k,l) = qc(1,j,k,l) ! rho
          q(2,j,k,l) = qc(2,j,k,l)/qc(1,j,k,l) ! u
          q(3,j,k,l) = qc(3,j,k,l)/qc(1,j,k,l) ! v
          q(4,j,k,l) = qc(4,j,k,l)/qc(1,j,k,l) ! w
          q(5,j,k,l) = ( qc(5,j,k,l) - 0.5d0*( qc(2,j,k,l)**2 +qc(3,j,k,l)**2 +qc(4,j,k,l)**2 )/qc(1,j,k,l) )/gam(j,k,l) ! p
          q(6:ndmax,j,k,l) = qc(6:ndmax,j,k,l) ! rhoy

          !*****************************************************
        enddo
      enddo
    enddo

    !$omp parallel do reduction(+:resid)
    do l=1,lmax
      do k=1,kmax
        do j=1,jmax
          ! compute l2 residual of rhs 
          resid(1) = resid(1) + s(1,j,k,l)*s(1,j,k,l) 
          resid(2) = resid(2) + s(2,j,k,l)*s(2,j,k,l) 
          resid(3) = resid(3) + s(3,j,k,l)*s(3,j,k,l) 
          resid(4) = resid(4) + s(4,j,k,l)*s(4,j,k,l) 
          resid(5) = resid(5) + s(5,j,k,l)*s(5,j,k,l)                       
          resid(6:ndmax) = resid(6:ndmax) + s(6:ndmax,j,k,l)*s(6:ndmax,j,k,l)
        enddo
      enddo
    enddo
  end do
  
  return
end subroutine step_RK4
subroutine firstOrderQ(q, ql, qr)
  !**********************************************************************
  !*     caluculate right-hand-side                                     *
  !**********************************************************************
  use param
  implicit none
  integer j,k,l,n
  double precision,dimension(ndmax,0:(jmax+1),0:(kmax+1),0:(lmax+1)) :: q
  double precision,dimension(ndmax,0:(jmax+1),0:(kmax+1),0:(lmax+1),ndim) :: ql, qr

  do n=1,ndim  
    do l=0,lmax+1  
      do k=0,kmax+1  
        do j=0,jmax+1  
          ql(:,j,k,l,n) = q(:,j,k,l)
          qr(:,j,k,l,n) = q(:,j,k,l)
        enddo
      enddo
    enddo
  enddo

  return
end subroutine firstOrderQ

subroutine muscl_va(q, ql, qr)
  !**********************************************************************
  !*     caluculate right-hand-side                                     *
  !**********************************************************************
  use param
  implicit none
  integer j,k,l,n,nd
  double precision dkp, dkpp1, dkpm1
  double precision dp, dm
  double precision eps, s ! limiter
  double precision,dimension(ndmax,0:(jmax+1),0:(kmax+1),0:(lmax+1)) :: q
  double precision,dimension(ndmax,0:(jmax+1),0:(kmax+1),0:(lmax+1),ndim) :: ql, qr
  integer,dimension(ndim) :: dl


  ! high order upwind
  ! dkp = -1.0d0 ! 2nd order upwind
  ! dkp = 1.0d0/3.0d0 ! 3rd order upwind
  ! dkp = 1.0d0  ! 2nd order central
  dkp = -1d0
  eps = 1.0d-6

  do nd=1,ndim  
    do l=0,lmax+1  
      do k=0,kmax+1  
        do j=0,jmax+1  
          ql(:,j,k,l,nd) = q(:,j,k,l)
          qr(:,j,k,l,nd) = q(:,j,k,l)
        enddo
      enddo
    enddo
  enddo

  do nd=1,ndim
    dl(:)=0
    dl(nd)=1
    do l=1,lmax
      do k=1,kmax
        do j=1,jmax
          do n=1,ndmax
            ! dx
            dp = q(n,j+dl(1),k+dl(2),l+dl(3)) - q(n,j,k,l)  ! dw(j+1)
            dm = q(n,j,k,l) - q(n,j-dl(1),k-dl(2),l-dl(3))  ! dw(j)

            ! set Limiter
            s = ( 2.0d0*dp*dm +eps)/( dp*dp + dm*dm +eps)   ! slope(van albada2)
            dkpp1 = 1.0d0 + s *dkp
            dkpm1 = 1.0d0 - s *dkp

            ! u_L and u_R 
            ql(n,j,k,l,nd) = q(n,j,k,l) + 0.25d0* s*( dkpm1*dm + dkpp1*dp )
            qr(n,j,k,l,nd) = q(n,j,k,l) - 0.25d0* s*( dkpp1*dm + dkpm1*dp )
          end do
        end do
      enddo 
    enddo 
  enddo 

  return
end subroutine muscl_va

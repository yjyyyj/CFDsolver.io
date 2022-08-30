subroutine outf(q, n)
  !**********************************************************************
  !*     output data                                                    *
  !**********************************************************************
  use param
  implicit none
  integer j,k,l,n
  double precision, dimension(ndmax,0:(jmax+1),0:(kmax+1),0:(lmax+1)) :: q
  character(6) loop
  
  write(loop,'(i6.6)') n
  ! write(*,*) loop
  ! fname = 
  ! write(*,*) fname

  open(100,file='output_'//loop//'.dat',form="formatted")
    do l=1,lmax
      do k=1,kmax
        do j=1,jmax
          ! mw0  = calmw(q(:,j,k)) 
          ! temp = mw0/q(1,j,k) * q(4,j,k)/Ru          ! T = (M/r)*(p/R)
          write(100,fmt='(11E25.15e3)',advance='No') xg(j),yg(k),zg(l),q(1:ndmax,j,k,l)
          write(100,*)
        end do
        write(100,*)
      end do
      write(100,*)
    end do
  close(100)

  return

end subroutine outf

subroutine outf_slice(q, n)
  !**********************************************************************
  !*     output data                                                    *
  !**********************************************************************
  use param
  implicit none
  integer j,k,l,n
  double precision, dimension(ndmax,0:(jmax+1),0:(kmax+1),0:(lmax+1)) :: q
  character(6) loop
  
  write(loop,'(i6.6)') n
  ! write(*,*) loop
  ! fname = 
  ! write(*,*) fname

  open(100,file='output_'//loop//'.dat',form="formatted")
    do l=1,lmax
      do k=1,kmax
        do j=1,jmax
          if( yg(k)==0.5d0 ) then
            write(100,fmt='(11E25.15e3)',advance='No') xg(j),zg(l),q(1:ndmax,j,k,l),gam(j,k,l)
            write(100,*)
          endif
        end do
        write(100,*)
      end do
    end do
  close(100)

  return

end subroutine outf_slice
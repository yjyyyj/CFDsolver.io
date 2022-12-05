subroutine outf(q, n)
  !**********************************************************************
  !*     output Q data in 3D                                            *
  !**********************************************************************
  use param
  implicit none
  integer j,k,l,n
  double precision, dimension(ndmax,0:(jmax+1),0:(kmax+1),0:(lmax+1)) :: q
  character(6) loop
  
  write(loop,'(i6.6)') n

  open(100,file='output_'//loop//'.dat',form="formatted")
    do l=1,lmax
      do k=1,kmax
        do j=1,jmax
          ! temp = calmw(q(:,j,k))/q(1,j,k) * q(4,j,k)/Ru          ! T = (M/r)*(p/R)
          write(100,fmt='(11E25.15e3)',advance='No') xg(j),yg(k),zg(l),q(1:ndmax,j,k,l),gam(j,k,l)
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
  !*     output Q data in 2D slice                                      *
  !**********************************************************************
  use param
  implicit none
  integer j,k,l,n
  double precision:: temp
  double precision, dimension(ndmax,0:(jmax+1),0:(kmax+1),0:(lmax+1)) :: q
  character(6) loop
  
  write(loop,'(i6.6)') n
  
  l = 2
  !**** out 2d slice data ********
  open(100,file='output_'//loop//'.dat',form="formatted")
    do k=1,kmax
      do j=1,jmax
        if( zg(l)==0.5d0 ) then
          temp = mw(j,k,l)/q(1,j,k,l) * q(5,j,k,l)/Ru           ! T = (M/r)*(p/R)
          write(100,fmt='(11E25.15e3)',advance='No') xg(j),yg(k),q(1:ndmax,j,k,l),temp,gam(j,k,l)
          write(100,*)
        endif
      end do
      write(100,*)
    end do

  close(100)

  return

end subroutine outf_slice

subroutine outf_1d(q, n)
  !**********************************************************************
  !*     output Q data in 1D                                            *
  !**********************************************************************
  use param
  implicit none
  integer j,k,l,n
  double precision, dimension(ndmax,0:(jmax+1),0:(kmax+1),0:(lmax+1)) :: q
  double precision :: temp
  character(6) loop
  
  write(loop,'(i6.6)') n

  k = 2
  l = 2

  !**** out 2d slice data ********
  open(100,file='output_'//loop//'.dat',form="formatted")
    do j=1,jmax
      temp = mw(j,k,l)/q(1,j,k,l) * q(5,j,k,l)/Ru          ! T = (M/r)*(p/R)
      write(100,fmt='(11E25.15e3)',advance='No') xg(j),q(1:ndmax,j,k,l),temp, gam(j,k,l)
      write(100,*)
    end do
  close(100)

  return

end subroutine outf_1d

subroutine outf_init(q)
  !**********************************************************************
  !*     initialize output data format                                  *
  !**********************************************************************
  use param
  implicit none
  integer j,k,l
  double precision, dimension(ndmax,0:(jmax+1),0:(kmax+1),0:(lmax+1)) :: q
  
  !!**** out 3d full data ********
  open(100,file="output_000000.dat",form="formatted")
    do l=1,lmax
      do k=1,kmax
        do j=1,jmax
          ! temp = mw0/q0(1,j,k) * q0(4,j,k)/Ru          ! T = (M/r)*(p/R)
          write(100,fmt='(11E25.15e3)',advance='No') xg(j),yg(k),zg(l),q(1:ndmax,j,k,l),gam(j,k,l)
          write(100,*)
        end do
        write(100,*)
      end do
      write(100,*)
    end do
  close(100)

  return

end subroutine outf_init

subroutine outf_slice_init(q)
  !**********************************************************************
  !*     initialize output data format                                  *
  !**********************************************************************
  use param
  implicit none
  integer j,k,l
  double precision:: temp
  double precision, dimension(ndmax,0:(jmax+1),0:(kmax+1),0:(lmax+1)) :: q
  
  l = 2
  !**** out 2d slice data ********
  open(100,file='output_000000.dat',form="formatted")
  do k=1,kmax
    do j=1,jmax
      if( zg(l)==0.5d0 ) then
        temp = mw(j,k,l)/q(1,j,k,l) * q(5,j,k,l)/Ru          ! T = (M/r)*(p/R)
        write(100,fmt='(11E25.15e3)',advance='No') xg(j),yg(k),q(1:ndmax,j,k,l),temp,gam(j,k,l)
        write(100,*)
      endif
    end do
    write(100,*)
  end do

  close(100)

  return

end subroutine outf_slice_init

subroutine outf_1d_init(q)
  !**********************************************************************
  !*     initialize output data format                                  *
  !**********************************************************************
  use param
  implicit none
  integer j,k,l
  double precision, dimension(ndmax,0:(jmax+1),0:(kmax+1),0:(lmax+1)) :: q
  double precision :: temp

  k = 2
  l = 2
  !**** out 2d slice data ********
  open(100,file='output_000000.dat',form="formatted")
    do j=1,jmax
      temp = mw(j,k,l)/q(1,j,k,l) * q(5,j,k,l)/Ru          ! T = (M/r)*(p/R)
      write(100,fmt='(11E25.15e3)',advance='No') xg(j),q(1:ndmax,j,k,l),temp,gam(j,k,l)
      write(100,*)
    end do
  close(100)

  return

end subroutine outf_1d_init

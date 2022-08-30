subroutine sumdf(qc, n)
  !**********************************************************************
  !*     output data                                                    *
  !**********************************************************************
  use param
  implicit none
  integer j,k,l,n
  double precision, dimension(ndmax) :: sums
  double precision, dimension(ndmax,jmax,kmax,lmax) :: qc

  sums(:) = 0.0d0

  do l=1,lmax-1
    do k=1,kmax-1
      do j=1,jmax-1
        sums(:) = sums(:) + qc(:,j,k,l)
      end do
    end do
  end do

  sums(:) = sums(:)/(jmax-1)/(kmax-1)/(lmax-1)
  sums(:) = (sums(:) - sums0(:)) / sums0(:)

  open(150,file="sums.dat",form="formatted",position='append')
    write(150,fmt='(I8, 9E20.10e3)',advance='No') n, sums
    write(*,fmt='(I8, 9E20.10e3)') n, sums
  close(150)

  return

end subroutine sumdf

subroutine residual(n)
  !**********************************************************************
  !*     output data                                                    *
  !**********************************************************************
  use param
  implicit none
  integer n

  resid = resid/dble((jmax)*(kmax)*(lmax) )
  resid = sqrt(resid)

  open(130,file="residual.dat",form="formatted",position='append')
    write(130,fmt='(I8, 9E20.10e3)',advance='No') n, resid
    write(*,fmt='(I8, 9E20.10e3)') n, resid
  close(130)

  return

end subroutine residual
module param
  !******************************************************************
  !*     function & parameter module                                *
  !******************************************************************
  implicit none
  
  double precision :: cfl
  double precision :: dt
  double precision :: t
  double precision :: pi
  double precision :: Ru

  double precision :: Mach
  double precision :: Re
  double precision :: Pr

  double precision,allocatable,dimension(:) :: gami
  double precision,allocatable,dimension(:) :: mwi 
  double precision,allocatable,dimension(:,:,:) :: gam
  double precision,allocatable,dimension(:,:,:) :: mw 
  
  double precision,allocatable,dimension(:) :: xg,yg,zg 
  double precision,allocatable,dimension(:) :: dx
  double precision,allocatable,dimension(:) :: sums0 
  double precision,allocatable,dimension(:) :: resid

  double precision,allocatable,dimension(:) :: qinf

  integer :: restart
  integer :: nout
  integer :: nmax
  integer :: ndim
  integer :: ndmax
  integer :: nspecies
  integer :: mconst
  integer :: jmax, kmax, lmax
  integer :: irhs
  integer :: acc
  integer :: ilhs
  integer :: vflag
  integer :: bc1, bc2, bc3, bc4, bc5, bc6

  ! for fmu calc in Sutherland’s eq.
  double precision :: tinf = 460d0
  double precision :: c2b  = 198.6d0/460d0  ! C/T0
  double precision :: c2bp = 198.6d0/460d0 + 1d0   ! (T0 + C)/T0

  contains

  double precision function calmw(qc)
    !******************************************************************
    !*     calc molar weight function                                 *
    !******************************************************************
    double precision,intent(in) :: qc(ndmax)
    double precision  :: r
    double precision  :: ry(nspecies)
    integer i

    ry(:) = qc(6:ndmax)
    r = sum(ry(:))
    calmw = 0d0
    do i = 1, nspecies
      calmw = calmw + (1d0/mwi(i))*ry(i)/r
    enddo
    calmw = 1d0/calmw
  end function 

  double precision function calgm(qc,mw_mix)
    !******************************************************************
    !*     calc G: 1/(gamma - 1) function                             *
    !******************************************************************
    double precision,intent(in) :: qc(ndmax)
    double precision,intent(in) :: mw_mix
    double precision :: r
    double precision :: ry(nspecies)
    integer i

    ry(:) = qc(6:ndmax)
    r = sum(ry(:))

    calgm = 0d0
    do i = 1, nspecies
      calgm = calgm + (gami(i)/mwi(i))*mw_mix*ry(i)/r
    enddo
  end function

end module param

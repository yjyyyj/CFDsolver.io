subroutine bc(q,fmu)
  !**********************************************************************
  !*     set boundary condition                                         *
  !**********************************************************************
  use param
  implicit none
  integer j
  double precision, dimension(ndmax,0:(jmax+1),0:(kmax+1),0:(lmax+1)) :: q
  double precision,dimension(0:(jmax+1),0:(kmax+1),0:(lmax+1)) :: fmu           

  select case( bc1 )
    case(1)
      ! outflow bc
      q(:,0,:,:) = q(:,1,:,:)
      fmu(0,:,:) = fmu(1,:,:)
    case(2)
      ! inflow bc
      q(1,0,:,:) = qinf(1)
      q(2,0,:,:) = qinf(2)
      q(3,0,:,:) = qinf(3)
      q(4,0,:,:) = qinf(4)
      q(5,0,:,:) = qinf(5)
      q(6,0,:,:) = qinf(6)
      q(7,0,:,:) = qinf(7)
      fmu(0,:,:) = fmu(1,:,:) ! tmp
    case(3)
      ! periodic bc
      q(:,0,:,:) = q(:,jmax-1,:,:) 
      fmu(0,:,:) = fmu(jmax-1,:,:)
  end select

  select case( bc2 )
    case(1)
      ! outflow bc
      if(Mach < 1d0) then ! sub-sonic
        q(:,jmax+1,:,:) = q(:,jmax,:,:)
        q(5,jmax+1,:,:) = qinf(5) ! use pinf
      else 
        q(:,jmax+1,:,:) = q(:,jmax,:,:)
      endif
      ! q(:,jmax+1,:,:) = q(:,jmax,:,:)
      fmu(jmax+1,:,:) = fmu(jmax,:,:)
    case(2)
      ! inflow bc
      q(1,jmax+1,:,:) = qinf(1)
      q(2,jmax+1,:,:) = qinf(2)
      q(3,jmax+1,:,:) = qinf(3)
      q(4,jmax+1,:,:) = qinf(4)
      q(5,jmax+1,:,:) = qinf(5)
      q(6,jmax+1,:,:) = qinf(6)
      q(7,jmax+1,:,:) = qinf(7)
      fmu(jmax+1,:,:) = fmu(jmax,:,:) ! tmp
    case(3)
      ! periodic bc
      q(:,jmax+1,:,:) = q(:,2,:,:)
      fmu(jmax+1,:,:) = fmu(2,:,:)
  end select
  
  select case( bc3 )
    case(1)
      ! outflow bc
      q(:,:,0,:) = q(:,:,1,:)
      fmu(:,0,:) = fmu(:,1,:)
    case(2)
      ! inflow bc
      ! q(:,:,0,:) = q0(:,1,:,:)
    case(3)
      ! periodic bc
      q(:,:,0,:) = q(:,:,kmax-1,:)
      fmu(:,0,:) = fmu(:,kmax-1,:)
  end select

  select case( bc4 )
    case(1)
      ! outflow bc
      q(:,:,kmax+1,:) = q(:,:,kmax,:)
      fmu(:,kmax+1,:) = fmu(:,kmax,:)
    case(2)
      ! inflow bc
      ! q(:,jmax+1) = qinf(:)
    case(3)
      ! periodic bc
    q(:,:,kmax+1,:) = q(:,:,2,:)
    fmu(:,kmax+1,:) = fmu(:,2,:)
  end select

  select case( bc5 )
    case(1)
      ! outflow bc
      q(:,:,:,0) = q(:,:,:,1)
      fmu(:,:,0) = fmu(:,:,1)
    case(2)
      ! inflow bc
      ! q(:,jmax+1) = qinf(:)
    case(3)
      ! periodic bc
      q(:,:,:,0) = q(:,:,:,lmax-1)
      fmu(:,:,0) = fmu(:,:,lmax-1)
    case(4)
      ! wall bc
      q(:,:,:,0) = q(:,:,:,1)
      fmu(:,:,0) = fmu(:,:,1)
      do j = 1, jmax
        if (xg(j) >= 2d0) then
          ! q(2:4,j,:,0) = 0d0 ! non-slip
          q(2:4,j,:,0) = -q(2:4,j,:,1) ! non-slip
        endif
      enddo
  end select

  select case( bc6 )
    case(1)
      ! outflow bc
      q(:,:,:,lmax+1) = q(:,:,:,lmax)
      fmu(:,:,lmax+1) = fmu(:,:,lmax)
    case(2)
      ! inflow bc (fixL)
      q(1,:,:,lmax+1) = qinf(1)
      q(2,:,:,lmax+1) = qinf(2)
      q(3,:,:,lmax+1) = qinf(3)
      q(4,:,:,lmax+1) = qinf(4)
      q(5,:,:,lmax+1) = qinf(5)
      q(6,:,:,lmax+1) = qinf(6)
      q(7,:,:,lmax+1) = qinf(7)
      fmu(:,:,lmax+1) = fmu(:,:,lmax) ! tmp
    case(3)
      ! periodic bc
      q(:,:,:,lmax+1) = q(:,:,:,2)
      fmu(:,:,lmax+1) = fmu(:,:,2)
  end select

  !*** termination
  return
end subroutine bc
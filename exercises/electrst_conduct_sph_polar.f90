! compile as gfortran electrst_conduct_sph_polar.f90 -llapack


! UD

! i -> r, j -> ph

program conducting_sphere
  implicit none
  integer, parameter :: ni = 30
  integer, parameter :: nj = 200
  real, parameter :: pi = 3.141592
  integer :: i,j,ka,kb,l,m
  real, dimension(1:ni*nj, 1:ni*nj) :: a
  real, dimension(1:ni, 1:ni) :: b
  real, dimension(1:ni, 1:ni) :: mi1, mi2
  real, dimension(1:ni*nj) :: c
  real, dimension(1:ni, 1:nj) :: sol
  real :: dr, dph, rmax, rmin, E0
  integer :: file1, file2
  integer :: fstat
  ! ! obstacle
  ! integer :: dkal, dkbl, dkam, dkbm
  ! real :: obst_val
  ! lapack
  integer, dimension(1:ni*nj) :: pivot
  integer :: info

  open(newunit=file1, file="out2.dat", status="replace")
  ! open(newunit=file2, file="poisson_obstacles/point.txt", status="old")

  rmax = 1. ! m
  rmin = 0.1 ! m
  E0 = 12 ! V/m ( 12 V / 6 m )

  dph = 2*pi/(nj-1)
  dr = (rmax-rmin)/(ni-1)

  a(:,:) = 0.
  b(:,:) = 0.
  mi1(:,:) = 0.
  mi2(:,:) = 0.
  c(:) = 0.

  ! ! Construction of matrix b:
  ! ! b(i,i): central cell (unknown)
  ! ! b(i,i+1): cell to the right
  ! ! b(i,i-1): cell to the left
  ! do i = 1,ni
  !   do j = 1,ni
  !     if ( j == i ) then
  !       b(i,j) = 1/(r(i)*dr) + 2/dr**2 + 2/(dph**2*r(i)**2)
  !     else if ( j == i+1 ) then
  !       b(i,j) = -1/(r(i)**2 * dph**2)
  !     else if ( j == i-1 ) then
  !       b(i,j) = -1/(r(i)**2 * dph**2)
  !     end if
  !   end do
  ! end do
  !
  !
  ! do i = 1,ni
  !   do j = 1,ni
  !     if ( j == i ) then
  !       mi1(i,j) = -(1/(r(i)*dr) + 1/dr**2 )
  !       mi2(i,j) = -( 1/dr**2 )
  !     end if
  !   end do
  ! end do
  !
  !
  ! ! Assembly of the complete stencil matrix
  ! do i = 1,nj
  !   do j = 1,nj
  !     if (j == i) then
  !       a( ((i-1)*ni+1) : i*ni, ((j-1)*ni+1): j*ni  ) = b
  !     else if ( j == i+1 ) then
  !       a( ((i-1)*ni+1) : i*ni, ((j-1)*ni+1): j*ni  ) = mi1
  !     else if ( j == i-1 ) then
  !       a( ((i-1)*ni+1) : i*ni, ((j-1)*ni+1): j*ni  ) = mi2
  !     end if
  !   end do
  ! end do

   ! call print_a(a,ni,nj)
   ! stop "Test only <<"


  ! ! Obstacle in the middle: inner boundaries
  ! ka = n/2
  ! kb = n/2
  !
  ! do
  !   read(file2,*, iostat=fstat) dkal, dkbl, dkam, dkbm, obst_val
  !   if (fstat > 0) stop "Bad file format!"
  !   if (fstat < 0) exit
  !   call get_l_from_ij(l,ka+dkal,kb+dkbl,n)
  !   call get_l_from_ij(m,ka+dkam,kb+dkbm,n)
  !   if (l == m) then
  !     a(l,:) = 0
  !     ! Warning: in the spec. file, the row of the cell where the condition
  !     ! is set (l==m) has to come first, then, the others
  !   end if
  !
  ! a(l,m) = obst_val
  ! end do


  do l = 1, ni*nj
    call get_ij_from_l(l,i,j,nj)
    a(l,l) = 1/(r(i)*dr) + 2/dr**2 + 2/(dph**2*r(i)**2)
    if (i+1<=ni) then
      call get_l_from_ij(m,i+1,j,nj)
      a(m,l) = -1/(r(i)*dr) - 1/dr**2
    end if
    if (i-1>=1) then
      call get_l_from_ij(m,i-1,j,nj)
      a(m,l) = -1/dr**2
    end if
    if (j+1<=nj) then
      call get_l_from_ij(m,i,j+1,nj)
      a(l,m) = -1/(r(i)**2*dph**2)
    end if
    if (j-1>=1) then
      call get_l_from_ij(m,i,j-1,nj)
      a(l,m) = -1/(r(i)**2*dph**2)
    end if

  end do

  ! call print_a(a,ni,nj)

  ! stop "end of test"



  do j = 1,nj
    call get_l_from_ij(l,1,j,nj)
    c(l) = 0. ! equipotential at a metal
    call get_l_from_ij(l,ni,j,nj)
    c(l) = -E0*r(ni)*cos(ph(j)) ! Constant vector E0*e_x
  end do

  do i = 1, ni
    call get_l_from_ij(l,i,1,nj)
    call get_l_from_ij(m,i,nj,nj)
    c(l) = c(m)
  end do

  ! ! Construction of the boundary vector
  ! ! a. Top and bottom boundaries (in i)
  ! do j = 1,n
  !   call get_l_from_ij(l,1,j,n)
  !   c(l) = 8./real(n)*j
  !   call get_l_from_ij(l,n,j,n)
  !   c(l) = 8./real(n)*j
  ! end do
  ! ! b. Left and right boundaries (in j)
  ! do i = 1,n
  !   call get_l_from_ij(l,i,1,n)
  !   c(l) = c(l) + 0.
  !   call get_l_from_ij(l,i,n,n)
  !   c(l) = c(l) + 8.
  ! end do
  !
  ! !call print_c(c,n)
  !
  call sgesv(ni*nj, 1, a, ni*nj, pivot, c, ni*nj, info)
  !if (info == 0) print*, "Solution successful"

  do l = 1, ni*nj
    call get_ij_from_l(l,i,j,nj)
    sol(i,j) = c(l)
  end do

  call print_sol(sol,ni,nj)

contains

  function r(i)
    integer, intent(in) :: i
    real :: r
      r = rmin + dr*(i-1)
  end function

  function ph(j)
    integer, intent(in) :: j
    real :: ph
      ph = dph*(j-1)
    end function

  subroutine get_l_from_ij(l,i,j,n)
    integer, intent(in) :: i,j,n
    integer, intent(out) :: l
    l = (i-1)*n + j
  end subroutine

  subroutine get_ij_from_l(l,i,j,n)
    integer, intent(in) :: l,n
    integer, intent(out) :: i,j
    i = l/n + 1
    j = mod(l,n)
    if (j == 0) then
      j = n
      i = i-1
    end if
  end subroutine

  subroutine print_a(a,ni,nj)
    real, intent(in) :: a(:,:)
    integer, intent(in) :: ni,nj
    integer :: i,j
    do i = 1,ni*nj
      do j = 1,ni*nj
        write(*,"(5f3.0)",advance="no") a(i,j)
        ! white space to separate blocks
        if (mod(j,nj) == 0) then
          write(*,"(a)",advance="no") " "
        end if
      end do
      write(*,*) ""
      ! white space to separate blocks
      if (mod(i,nj) == 0) then
        write(*,*) " "
      end if
    end do
  end subroutine

  subroutine print_c(c,n)
    real, intent(in) :: c(:)
    integer, intent(in) :: n
    integer :: i
    do i = 1,n*n
      write(*,"(10f8.4)") c(i)
    end do
  end subroutine

  subroutine print_sol(a,ni,nj)
    real, intent(in) :: a(:,:)
    integer, intent(in) :: ni,nj
    integer :: i,j
    do i = 1,ni
      do j = 1,nj
        write(file1,"(e12.3,t1)",advance="no") a(i,j)
      end do
      write(file1,*) ""
    end do
  end subroutine


end program

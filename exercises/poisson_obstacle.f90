! compile as gfortran poisson_obstacle.f90 -llapack

program poisson_obstacle
  implicit none
  integer, parameter :: n = 35
  integer :: i,j,ka,kb,l,m
  real, dimension(1:n*n, 1:n*n) :: a
  real, dimension(1:n, 1:n) :: b
  real, dimension(1:n, 1:n) :: id
  real, dimension(1:n*n) :: c
  real, dimension(1:n, 1:n) :: sol
  integer :: file1, file2
  integer :: fstat
  ! obstacle
  integer :: dkal, dkbl, dkam, dkbm
  real :: obst_val
  ! lapack
  integer, dimension(1:n*n) :: pivot
  integer :: info

  open(newunit=file1, file="out2.dat", status="replace")
  open(newunit=file2, file="poisson_obstacles/point.txt", status="old")

  a(:,:) = 0.
  b(:,:) = 0.
  id(:,:) = 0.
  c(:) = 0.

  ! Construction of matrix b
  do i = 1,n
    do j = 1,n
      if ( j == i ) then
        b(i,j) = 4.
      else if ( j == i+1 ) then
        b(i,j) = -1.
      else if ( j == i-1 ) then
        b(i,j) = -1.
      end if
    end do
  end do


  ! Construction of the negative identiy matrix
  do i = 1,n
    do j = 1,n
      if ( j == i ) then
        id(i,j) = -1.
      end if
    end do
  end do


  ! Assembly of the complete stencil matrix
  do i = 1,n
    do j = 1,n
      if (j == i) then
        a( ((i-1)*n+1) : i*n, ((j-1)*n+1): j*n  ) = b
      else if ( j == i+1 ) then
        a( ((i-1)*n+1) : i*n, ((j-1)*n+1): j*n  ) = id
      else if ( j == i-1 ) then
        a( ((i-1)*n+1) : i*n, ((j-1)*n+1): j*n  ) = id
      end if
    end do
  end do


  ! Obstacle in the middle: inner boundaries
  ka = n/2
  kb = n/2

  do
    read(file2,*, iostat=fstat) dkal, dkbl, dkam, dkbm, obst_val
    if (fstat > 0) stop "Bad file format!"
    if (fstat < 0) exit
    call get_l_from_ij(l,ka+dkal,kb+dkbl,n)
    call get_l_from_ij(m,ka+dkam,kb+dkbm,n)
    if (l == m) then
      a(l,:) = 0
      ! Warning: in the spec. file, the row of the cell where the condition
      ! is set (l==m) has to come first, then, the others
    end if

  a(l,m) = obst_val
  end do


  ! Example of manual obstacle handling
  ! ! Cell where the obstacle is
  ! call get_l_from_ij(l,ka,kb,n)
  ! a(l,:) = 0.
  ! a(l,l) = 1.
  ! ! Neighbor cells (cross +)
  ! call get_l_from_ij(l,ka,kb+1,n)
  ! a(l,:) = 0.
  ! a(l,l) = 1.
  ! call get_l_from_ij(m,ka,kb+2,n)
  ! a(l,m) = -1.


  !call print_a(a,n)

  ! Construction of the boundary vector
  ! a. Top and bottom boundaries (in i)
  do j = 1,n
    call get_l_from_ij(l,1,j,n)
    c(l) = 8./real(n)*j
    call get_l_from_ij(l,n,j,n)
    c(l) = 8./real(n)*j
  end do
  ! b. Left and right boundaries (in j)
  do i = 1,n
    call get_l_from_ij(l,i,1,n)
    c(l) = c(l) + 0.
    call get_l_from_ij(l,i,n,n)
    c(l) = c(l) + 8.
  end do

  !call print_c(c,n)

  call sgesv(n*n, 1, a, n*n, pivot, c, n*n, info)
  !if (info == 0) print*, "Solution successful"

  do l = 1, n*n
    call get_ij_from_l(l,i,j,n)
    sol(i,j) = c(l)
  end do

  call print_sol(sol,n)

contains

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

  subroutine print_a(a,n)
    real, intent(in) :: a(:,:)
    integer, intent(in) :: n
    integer :: i,j
    do i = 1,n*n
      do j = 1,n*n
        write(*,"(5f3.0)",advance="no") a(i,j)
        ! white space to separate blocks
        if (mod(j,n) == 0) then
          write(*,"(a)",advance="no") " "
        end if
      end do
      write(*,*) ""
      ! white space to separate blocks
      if (mod(i,n) == 0) then
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

  subroutine print_sol(a,n)
    real, intent(in) :: a(:,:)
    integer, intent(in) :: n
    integer :: i,j
    do i = 1,n
      do j = 1,n
        write(file1,"(10f8.4,t1)",advance="no") a(i,j)
      end do
      write(file1,*) ""
    end do
  end subroutine


end program

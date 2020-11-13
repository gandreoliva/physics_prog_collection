! compile as gfortran poisson_freeflow.f90 -llapack
! Poisson equation: free uniform flow
! -----------------------------------

program poisson_freeflow
  implicit none
  integer, parameter :: n = 15
  integer :: i,j,l
  real, dimension(1:n*n, 1:n*n) :: a
  real, dimension(1:n, 1:n) :: b
  real, dimension(1:n, 1:n) :: id
  real, dimension(1:n*n) :: c
  real, dimension(1:n, 1:n) :: sol
  integer :: file1
  ! lapack
  integer, dimension(1:n*n) :: pivot
  integer :: info

  open(newunit=file1, file="out1.dat", status="replace")

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

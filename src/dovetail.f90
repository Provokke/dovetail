module dovetail
  use, intrinsic :: iso_fortran_env, only: dp => real64
  implicit none
  private

  public :: dp
  public :: lu_factor, lu_solve, solve, determinant, inverse
  public :: cholesky, conjugate_gradient

contains

  subroutine lu_factor(a, pivot, info, parity)
    real(dp), intent(inout) :: a(:, :)
    integer, intent(out) :: pivot(:)
    integer, intent(out) :: info
    integer, intent(out), optional :: parity
    real(dp) :: row(size(a, 2))
    integer :: n, j, k, p, swaps

    n = size(a, 1)
    info = 0
    swaps = 0
    pivot = 0
    if (present(parity)) parity = 0

    do k = 1, n
      p = k - 1 + maxloc(abs(a(k:n, k)), dim=1)
      pivot(k) = p
      if (abs(a(p, k)) <= tiny(1.0_dp)) then
        info = k
        return
      end if
      if (p /= k) then
        row = a(k, :)
        a(k, :) = a(p, :)
        a(p, :) = row
        swaps = swaps + 1
      end if
      a(k + 1:n, k) = a(k + 1:n, k) / a(k, k)
      do j = k + 1, n
        a(k + 1:n, j) = a(k + 1:n, j) - a(k + 1:n, k) * a(k, j)
      end do
    end do

    if (present(parity)) parity = merge(1, -1, mod(swaps, 2) == 0)
  end subroutine lu_factor

  subroutine lu_solve(a, pivot, b)
    real(dp), intent(in) :: a(:, :)
    integer, intent(in) :: pivot(:)
    real(dp), intent(inout) :: b(:)
    real(dp) :: held
    integer :: n, k

    n = size(a, 1)
    do k = 1, n
      if (pivot(k) /= k) then
        held = b(k)
        b(k) = b(pivot(k))
        b(pivot(k)) = held
      end if
    end do
    do k = 2, n
      b(k) = b(k) - dot_product(a(k, 1:k - 1), b(1:k - 1))
    end do
    do k = n, 1, -1
      b(k) = (b(k) - dot_product(a(k, k + 1:n), b(k + 1:n))) / a(k, k)
    end do
  end subroutine lu_solve

  subroutine solve(a, b, x, info)
    real(dp), intent(in) :: a(:, :)
    real(dp), intent(in) :: b(:)
    real(dp), intent(out) :: x(:)
    integer, intent(out) :: info
    real(dp) :: work(size(a, 1), size(a, 2))
    integer :: pivot(size(a, 1))

    work = a
    x = b
    call lu_factor(work, pivot, info)
    if (info /= 0) return
    call lu_solve(work, pivot, x)
  end subroutine solve

  function determinant(a) result(det)
    real(dp), intent(in) :: a(:, :)
    real(dp) :: det
    real(dp) :: work(size(a, 1), size(a, 2))
    integer :: pivot(size(a, 1))
    integer :: info, parity, k

    work = a
    call lu_factor(work, pivot, info, parity)
    if (info /= 0) then
      det = 0.0_dp
      return
    end if
    det = real(parity, dp)
    do k = 1, size(a, 1)
      det = det * work(k, k)
    end do
  end function determinant

  subroutine inverse(a, ainv, info)
    real(dp), intent(in) :: a(:, :)
    real(dp), intent(out) :: ainv(:, :)
    integer, intent(out) :: info
    real(dp) :: work(size(a, 1), size(a, 2))
    integer :: pivot(size(a, 1))
    integer :: k

    work = a
    ainv = 0.0_dp
    call lu_factor(work, pivot, info)
    if (info /= 0) return
    do k = 1, size(a, 1)
      ainv(k, k) = 1.0_dp
      call lu_solve(work, pivot, ainv(:, k))
    end do
  end subroutine inverse

  subroutine cholesky(a, l, info)
    real(dp), intent(in) :: a(:, :)
    real(dp), intent(out) :: l(:, :)
    integer, intent(out) :: info
    real(dp) :: diagonal
    integer :: n, i, j

    n = size(a, 1)
    info = 0
    l = 0.0_dp
    do j = 1, n
      diagonal = a(j, j) - dot_product(l(j, 1:j - 1), l(j, 1:j - 1))
      if (diagonal <= 0.0_dp) then
        info = j
        return
      end if
      l(j, j) = sqrt(diagonal)
      do i = j + 1, n
        l(i, j) = (a(i, j) - dot_product(l(i, 1:j - 1), l(j, 1:j - 1))) / l(j, j)
      end do
    end do
  end subroutine cholesky

  subroutine conjugate_gradient(a, b, x, tolerance, max_iterations, iterations, converged)
    real(dp), intent(in) :: a(:, :)
    real(dp), intent(in) :: b(:)
    real(dp), intent(inout) :: x(:)
    real(dp), intent(in) :: tolerance
    integer, intent(in) :: max_iterations
    integer, intent(out) :: iterations
    logical, intent(out) :: converged
    real(dp) :: r(size(b)), p(size(b)), ap(size(b))
    real(dp) :: rr, rr_next, alpha, target

    r = b - matmul(a, x)
    p = r
    rr = dot_product(r, r)
    target = tolerance * max(norm2(b), 1.0_dp)
    iterations = 0
    converged = sqrt(rr) <= target

    do while (.not. converged .and. iterations < max_iterations)
      ap = matmul(a, p)
      alpha = rr / dot_product(p, ap)
      x = x + alpha * p
      r = r - alpha * ap
      rr_next = dot_product(r, r)
      iterations = iterations + 1
      converged = sqrt(rr_next) <= target
      p = r + (rr_next / rr) * p
      rr = rr_next
    end do
  end subroutine conjugate_gradient

end module dovetail


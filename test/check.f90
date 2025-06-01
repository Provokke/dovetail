program check
  use dovetail
  implicit none

  integer :: failures

  failures = 0
  call test_solve()
  call test_pivoting()
  call test_singular()
  call test_determinant()
  call test_inverse()
  call test_cholesky()
  call test_conjugate_gradient()
  call test_hilbert()

  if (failures > 0) then
    print '(i0, a)', failures, ' checks failed'
    error stop 1
  end if
  print '(a)', 'all checks passed'

contains

  subroutine expect(label, condition)
    character(*), intent(in) :: label
    logical, intent(in) :: condition

    if (condition) then
      print '(2a)', 'ok    ', label
    else
      print '(2a)', 'FAIL  ', label
      failures = failures + 1
    end if
  end subroutine expect

  subroutine test_solve()
    real(dp) :: a(3, 3), b(3), x(3)
    integer :: info

    a = reshape([2.0_dp, -3.0_dp, -2.0_dp, 1.0_dp, -1.0_dp, 1.0_dp, -1.0_dp, 2.0_dp, 2.0_dp], [3, 3])
    b = [8.0_dp, -11.0_dp, -3.0_dp]
    call solve(a, b, x, info)
    call expect('solve reports success', info == 0)
    call expect('solve finds the solution', maxval(abs(x - [2.0_dp, 3.0_dp, -1.0_dp])) < 1.0e-12_dp)
    call expect('solve leaves a small residual', maxval(abs(matmul(a, x) - b)) < 1.0e-12_dp)
  end subroutine test_solve

  subroutine test_pivoting()
    real(dp) :: a(2, 2), x(2)
    integer :: info

    a = reshape([0.0_dp, 1.0_dp, 1.0_dp, 1.0_dp], [2, 2])
    call solve(a, [3.0_dp, 5.0_dp], x, info)
    call expect('zero leading entry needs a row swap', info == 0 .and. maxval(abs(x - [2.0_dp, 3.0_dp])) < 1.0e-14_dp)
  end subroutine test_pivoting

  subroutine test_singular()
    real(dp) :: a(3, 3), x(3)
    integer :: info

    a = reshape([1.0_dp, 2.0_dp, 3.0_dp, 2.0_dp, 4.0_dp, 6.0_dp, 1.0_dp, 0.0_dp, 1.0_dp], [3, 3])
    call solve(a, [1.0_dp, 1.0_dp, 1.0_dp], x, info)
    call expect('singular matrix is reported', info /= 0)
    call expect('singular matrix has zero determinant', determinant(a) == 0.0_dp)
  end subroutine test_singular

  subroutine test_determinant()
    real(dp) :: a(3, 3), swapped(2, 2)

    a = reshape([6.0_dp, 4.0_dp, 2.0_dp, 1.0_dp, -2.0_dp, 8.0_dp, 1.0_dp, 5.0_dp, 7.0_dp], [3, 3])
    swapped = reshape([0.0_dp, 1.0_dp, 1.0_dp, 0.0_dp], [2, 2])
    call expect('determinant of a 3x3', abs(determinant(a) + 306.0_dp) < 1.0e-10_dp)
    call expect('determinant tracks row swaps', abs(determinant(swapped) + 1.0_dp) < 1.0e-14_dp)
  end subroutine test_determinant

  subroutine test_inverse()
    real(dp) :: a(4, 4), ainv(4, 4), eye(4, 4)
    integer :: info, i

    call random_seed_fixed()
    call random_number(a)
    eye = 0.0_dp
    do i = 1, 4
      a(i, i) = a(i, i) + 4.0_dp
      eye(i, i) = 1.0_dp
    end do
    call inverse(a, ainv, info)
    call expect('inverse reports success', info == 0)
    call expect('a times its inverse is the identity', maxval(abs(matmul(a, ainv) - eye)) < 1.0e-12_dp)
  end subroutine test_inverse

  subroutine test_cholesky()
    real(dp) :: a(3, 3), l(3, 3), bad(2, 2), lbad(2, 2)
    integer :: info

    a = reshape([4.0_dp, 12.0_dp, -16.0_dp, 12.0_dp, 37.0_dp, -43.0_dp, -16.0_dp, -43.0_dp, 98.0_dp], [3, 3])
    call cholesky(a, l, info)
    call expect('cholesky reports success', info == 0)
    call expect('cholesky factor matches', maxval(abs(l(:, 1) - [2.0_dp, 6.0_dp, -8.0_dp])) < 1.0e-12_dp)
    call expect('cholesky reconstructs a', maxval(abs(matmul(l, transpose(l)) - a)) < 1.0e-10_dp)
    call expect('cholesky factor is lower triangular', l(1, 2) == 0.0_dp .and. l(2, 3) == 0.0_dp)

    bad = reshape([1.0_dp, 2.0_dp, 2.0_dp, 1.0_dp], [2, 2])
    call cholesky(bad, lbad, info)
    call expect('indefinite matrix is rejected', info == 2)
  end subroutine test_cholesky

  subroutine test_conjugate_gradient()
    integer, parameter :: n = 50
    real(dp) :: a(n, n), b(n), x(n), direct(n)
    integer :: i, iterations, info
    logical :: converged

    a = 0.0_dp
    do i = 1, n
      a(i, i) = 2.0_dp
      if (i > 1) a(i, i - 1) = -1.0_dp
      if (i < n) a(i, i + 1) = -1.0_dp
    end do
    b = 1.0_dp
    x = 0.0_dp
    call conjugate_gradient(a, b, x, 1.0e-12_dp, 200, iterations, converged)
    call solve(a, b, direct, info)
    call expect('conjugate gradient converges', converged)
    call expect('conjugate gradient needs at most n steps', iterations <= n)
    call expect('conjugate gradient agrees with lu', maxval(abs(x - direct)) < 1.0e-8_dp)

    x = direct
    call conjugate_gradient(a, b, x, 1.0e-8_dp, 200, iterations, converged)
    call expect('an exact start takes no steps', converged .and. iterations == 0)

    x = 0.0_dp
    call conjugate_gradient(a, b, x, 1.0e-12_dp, 3, iterations, converged)
    call expect('iteration cap is honoured', .not. converged .and. iterations == 3)
  end subroutine test_conjugate_gradient

  subroutine test_hilbert()
    integer, parameter :: n = 6
    real(dp) :: h(n, n), ones(n), b(n), x(n)
    integer :: i, j, info

    do j = 1, n
      do i = 1, n
        h(i, j) = 1.0_dp / real(i + j - 1, dp)
      end do
    end do
    ones = 1.0_dp
    b = matmul(h, ones)
    call solve(h, b, x, info)
    call expect('ill-conditioned hilbert system', info == 0 .and. maxval(abs(x - ones)) < 1.0e-6_dp)
  end subroutine test_hilbert

  subroutine random_seed_fixed()
    integer :: length
    integer, allocatable :: seed(:)

    call random_seed(size=length)
    allocate (seed(length))
    seed = 20261
    call random_seed(put=seed)
  end subroutine random_seed_fixed

end program check

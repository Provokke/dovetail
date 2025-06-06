<div align="center">

# dovetail

Dense linear solvers in modern Fortran. One module, no BLAS or LAPACK required.

![License](https://img.shields.io/badge/LICENSE-MIT-blue?style=for-the-badge)  ![Language](https://img.shields.io/badge/FORTRAN-734f96?style=for-the-badge&logo=fortran&logoColor=white)  ![Tests](https://img.shields.io/badge/TESTS-21-brightgreen?style=for-the-badge)

Refer to [LICENSE](LICENSE) for the legal bits.

</div>

## What's in it

| Routine | Purpose |
| --- | --- |
| `lu_factor` | in-place LU factorisation with partial pivoting |
| `lu_solve` | solves with an existing factorisation |
| `solve` | one-call solve of `A x = b` |
| `determinant` | determinant from the LU factors, with the sign tracked through row swaps |
| `inverse` | explicit inverse, column by column |
| `cholesky` | lower Cholesky factor of a symmetric positive definite matrix |
| `conjugate_gradient` | iterative solve for symmetric positive definite systems |

All reals are `real64`, exported as the kind `dp`. Failures come back through an `info` argument in the LAPACK manner: zero for success, otherwise the column where the factorisation broke down.

```fortran
program demo
  use dovetail
  implicit none
  real(dp) :: a(3, 3), b(3), x(3)
  integer :: info

  a = reshape([2.0_dp, -3.0_dp, -2.0_dp, 1.0_dp, -1.0_dp, 1.0_dp, -1.0_dp, 2.0_dp, 2.0_dp], [3, 3])
  b = [8.0_dp, -11.0_dp, -3.0_dp]

  call solve(a, b, x, info)
  if (info /= 0) error stop 'singular matrix'
  print *, x
end program demo
```

Factor once and reuse it when there are several right-hand sides:

```fortran
call lu_factor(a, pivot, info)
call lu_solve(a, pivot, b1)
call lu_solve(a, pivot, b2)
```

## 🏗️ Architecture

```
 src/dovetail.f90  (module dovetail, private by default)
 ┌──────────────────────────────────────────────────────────┐
 │ lu_factor(a, pivot, info [, parity])   unblocked, in place│
 │      ├──> lu_solve(a, pivot, b)    forward / back subst.  │
 │      ├──> solve(a, b, x, info)     copy, factor, solve    │
 │      ├──> determinant(a)           product of U diagonal  │
 │      └──> inverse(a, ainv, info)   lu_solve per column    │
 │                                                          │
 │ cholesky(a, l, info)               reads lower triangle   │
 │ conjugate_gradient(a, b, x, ...)   dense matvec, no precond│
 └──────────────────────────────────────────────────────────┘
 test/check.f90   Makefile   fpm.toml
```

- These are textbook unblocked algorithms. They're easy to read and correct, and they'll be well behind an optimised LAPACK on large matrices.
- `solve`, `determinant` and `inverse` copy the matrix onto the stack before factorising. For big systems call `lu_factor` on your own array to avoid the copy.
- `conjugate_gradient` takes a dense matrix. There's no sparse storage and no preconditioner, and it reports the iteration count and whether it converged.
- `cholesky` reads only the lower triangle and doesn't check symmetry. A matrix that isn't positive definite comes back with `info` set to the failing column.

## 🚀 Quick start

```
make test
fpm test
```

`make` uses `gfortran -O2 -std=f2018 -Wall -Wextra -fcheck=all`, so bounds checking is on. I ran it with GNU Fortran 14.4. Building and running the demo above prints:

```
   2.0000000000000022        2.9999999999999991      -0.99999999999999711
```

The exact solution is `2, 3, -1`.

## 🧪 Testing

21 checks, all passing. `test/check.f90` is a plain program with an `expect` helper that prints `ok` or `FAIL` per check, and ends with `all checks passed` or exits non-zero.

- `solve` (3): success flag, the known solution and a small residual.
- Pivoting (1): a zero leading entry needs a row swap.
- Singular matrix (2): `info` is non-zero and the determinant is exactly zero.
- Determinant (2): a 3x3 with a known value and a row-swap sign check.
- Inverse (2): success, and `a` times its inverse is the identity on a fixed-seed 4x4.
- Cholesky (5): success, the known factor, `L * L^T` rebuilds `a`, the factor is lower triangular and an indefinite matrix is rejected at column 2.
- Conjugate gradient (5): converges on a 50x50 tridiagonal system, in at most `n` steps, agrees with LU, an exact starting guess takes no steps and the iteration cap is honoured.
- Hilbert (1): a 6x6 ill-conditioned system still solves to within 1e-6.

## Licence

MIT.


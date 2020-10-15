# Skeel's condition number in MATLAB

`skeel(A)` estimates Skeel's condition number, `norm(abs(inv(A))*abs(A), inf),` without computing `abs(inv(A))*abs(A)`. `skeel(A)` is always less than or
equal to `cond(A, inf)`. In practice, `skeel(A)` can be much less than `cond(A, inf)`. `skeel(A)` is invariant to row scaling.

`skeel(A, p)` directly computes Skeel's condition number in the `p`-norm: `norm(abs(inv(A))*abs(A), p)`.

## Example:

```
>> A = [1 0; 0 1e9];
>> cond(A)

ans =

   1.0000e+09

>> skeel(A)

ans =

     1
```

## References:

The esimator is based on LAPACK's [`DLA_GERCOND`](http://www.netlib.org/lapack/explore-html/dd/d9a/group__double_g_ecomputational_ga5539077fbd3a92c4d92b75bf58da5db3.html) and [`DLACN2`](http://www.netlib.org/lapack/explore-html/d8/d9b/group__double_o_t_h_e_rauxiliary_ga9b62da514b4a671acd3e3f63d018f01e.html).

[1] Robert D. Skeel, "Scaling for numerical stability in Gaussian elimination", J. ACM, 26 (1979), pp. 494-526.

[2] N.J. Higham, "FORTRAN codes for estimating the one-norm of a real or complex matrix, with applications to condition estimation", ACM Trans. Math. Soft., 14 (1988), pp. 381-396.

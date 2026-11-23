//! Small numerical helpers.

use marlu::ndarray::Array2;

/// Least-squares fit of a polynomial of order `poly_order` to the data `y` at
/// coordinates `x`, ignoring any points where `y` is `NaN`.
///
/// Returns the `poly_order + 1` coefficients `c` in the basis of `x`, lowest
/// order first, such that `y ≈ sum(c[i] * x^i)`.
///
/// The fit is done in a domain normalised to `[-1, 1]`, so that it is stable
/// for coordinates like frequencies in Hz, and the coefficients are then
/// transformed back to the basis of `x`. For this, `x` is assumed to be
/// sorted: the domain is taken from its first and last elements.
///
/// Returns `None` if:
/// - `x` and `y` are not the same length,
/// - there are not more valid points than `poly_order`,
/// - the first and last `x` are the same, or
/// - the fit is singular.
///
/// # Examples
///
/// ```rust
/// use birli::math::fit_polynomial;
///
/// // y = 1 + 2x, with a flagged point
/// let x = [0.0, 1.0, 2.0, 3.0];
/// let y = [1.0, 3.0, f32::NAN, 7.0];
/// let coeffs = fit_polynomial(&x, &y, 1).unwrap();
/// assert!((coeffs[0] - 1.0).abs() < 1e-5);
/// assert!((coeffs[1] - 2.0).abs() < 1e-5);
/// ```
pub fn fit_polynomial(x: &[f64], y: &[f32], poly_order: usize) -> Option<Vec<f32>> {
    if x.len() != y.len() {
        return None;
    }

    let mut valid_indices = Vec::with_capacity(x.len());
    let mut valid_y = Vec::with_capacity(x.len());
    for (i, &val) in y.iter().enumerate() {
        if !val.is_nan() {
            valid_indices.push(i);
            valid_y.push(val as f64);
        }
    }

    if valid_indices.len() <= poly_order {
        return None;
    }

    // Normalise x for stability during the fit: x_norm = (x - f_mid) / s_scale
    // where s_scale is half the width of the domain.
    let f_min = x[0];
    let f_max = x[x.len() - 1];
    let f_width = f_max - f_min;

    // Avoid division by zero if all x are the same
    if f_width.abs() < 1e-12 {
        return None;
    }

    let f_mid = f_min + f_width / 2.0;
    let s_scale = f_width / 2.0;

    let x_norm: Vec<f64> = valid_indices
        .iter()
        .map(|&i| (x[i] - f_mid) / s_scale)
        .collect();

    // Build Normal Equations: (X^T X) beta = X^T y
    let mut xtx = Array2::<f64>::zeros((poly_order + 1, poly_order + 1));
    let mut xty = Array2::<f64>::zeros((poly_order + 1, 1));

    for (i, &x_val) in x_norm.iter().enumerate() {
        let y_val = valid_y[i];

        // x_pow = [1, x, x^2, ...]
        let mut x_pow = vec![1.0; poly_order + 1];
        for k in 1..=poly_order {
            x_pow[k] = x_pow[k - 1] * x_val;
        }

        for j in 0..=poly_order {
            xty[[j, 0]] += x_pow[j] * y_val;
            for k in 0..=poly_order {
                xtx[[j, k]] += x_pow[j] * x_pow[k];
            }
        }
    }

    // Solve xtx * beta = xty
    let n = poly_order + 1;
    let mut aug = Array2::<f64>::zeros((n, n + 1));
    for i in 0..n {
        for j in 0..n {
            aug[[i, j]] = xtx[[i, j]];
        }
        aug[[i, n]] = xty[[i, 0]];
    }

    // Gaussian elimination with partial pivoting
    for i in 0..n {
        let mut pivot = aug[[i, i]];
        let mut pivot_row = i;
        for k in i + 1..n {
            if aug[[k, i]].abs() > pivot.abs() {
                pivot = aug[[k, i]];
                pivot_row = k;
            }
        }
        if pivot.abs() < 1e-12 {
            return None; // Singular
        }
        if pivot_row != i {
            for j in i..=n {
                let tmp = aug[[i, j]];
                aug[[i, j]] = aug[[pivot_row, j]];
                aug[[pivot_row, j]] = tmp;
            }
        }
        for j in i..=n {
            aug[[i, j]] /= pivot;
        }
        for k in 0..n {
            if k != i {
                let factor = aug[[k, i]];
                for j in i..=n {
                    aug[[k, j]] -= factor * aug[[i, j]];
                }
            }
        }
    }

    let mut beta_norm = vec![0.0; n];
    for i in 0..n {
        beta_norm[i] = aug[[i, n]];
    }

    // Convert coefficients from the normalised basis to the basis of x (f)
    // y = sum(beta_norm[j] * ((f - f_mid)/s)^j)
    // y = sum(beta_norm[j] * s^-j * sum(binom(j,k) * f^k * (-f_mid)^(j-k)))
    // coeff of f^k (c_k) is sum over j>=k of beta_norm[j] * s^-j * binom(j,k) * (-f_mid)^(j-k)

    let mut coeffs = vec![0.0; n];
    let mut binom = vec![vec![0.0; n]; n];

    // Precompute binomial coefficients
    for i in 0..n {
        binom[i][0] = 1.0;
        for j in 1..=i {
            binom[i][j] = binom[i - 1][j - 1] + binom[i - 1][j];
        }
    }

    for (k, coeff) in coeffs.iter_mut().enumerate().take(n) {
        let mut sum = 0.0;
        for j in k..n {
            let term = beta_norm[j]
                * s_scale.powi(-(j as i32))
                * binom[j][k]
                * (-f_mid).powi((j - k) as i32);
            sum += term;
        }
        *coeff = sum as f32;
    }

    Some(coeffs)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_fit_polynomial_exact() {
        // Fit y = 1 + 2x + 3x^2
        let x: Vec<f64> = (0..10).map(|i| i as f64).collect();
        let y: Vec<f32> = x
            .iter()
            .map(|&val| (1.0 + 2.0 * val + 3.0 * val * val) as f32)
            .collect();

        let coeffs = fit_polynomial(&x, &y, 2).expect("Fit failed");

        assert_eq!(coeffs.len(), 3);
        assert!((coeffs[0] - 1.0).abs() < 1e-4);
        assert!((coeffs[1] - 2.0).abs() < 1e-4);
        assert!((coeffs[2] - 3.0).abs() < 1e-4);
    }

    #[test]
    fn test_fit_polynomial_order_zero_is_the_mean() {
        let x = vec![0.0, 1.0, 2.0, 3.0];
        let y = vec![1.0, 2.0, 3.0, 6.0];
        let coeffs = fit_polynomial(&x, &y, 0).expect("Fit failed");
        assert_eq!(coeffs.len(), 1);
        assert!((coeffs[0] - 3.0).abs() < 1e-6);
    }

    #[test]
    fn test_fit_polynomial_with_nans() {
        // Fit y = 2x
        let x = vec![0.0, 1.0, 2.0, 3.0, 4.0];
        let y = vec![0.0, 2.0, f32::NAN, 6.0, 8.0];

        let coeffs = fit_polynomial(&x, &y, 1).expect("Fit failed");

        assert_eq!(coeffs.len(), 2);
        assert!(coeffs[0].abs() < 1e-4); // intercept 0
        assert!((coeffs[1] - 2.0).abs() < 1e-4); // slope 2
    }

    /// The fit is least-squares when the data aren't an exact polynomial.
    #[test]
    fn test_fit_polynomial_least_squares() {
        // the least-squares line through these is y = 0.5 + x
        let x = vec![0.0, 1.0, 2.0, 3.0];
        let y = vec![1.0, 1.0, 2.0, 4.0];
        let coeffs = fit_polynomial(&x, &y, 1).expect("Fit failed");
        assert!((coeffs[0] - 0.5).abs() < 1e-5, "{coeffs:?}");
        assert!((coeffs[1] - 1.0).abs() < 1e-5, "{coeffs:?}");
    }

    #[test]
    fn test_fit_polynomial_not_enough_points() {
        // an order 2 polynomial has 3 coefficients, so it needs 3 points.
        let x = vec![0.0, 1.0];
        let y = vec![0.0, 1.0];
        assert!(fit_polynomial(&x, &y, 2).is_none());

        // NaNs don't count
        let x = vec![0.0, 1.0, 2.0];
        let y = vec![0.0, f32::NAN, 4.0];
        assert!(fit_polynomial(&x, &y, 2).is_none());
        assert!(fit_polynomial(&x, &y, 1).is_some());

        // nor does nothing at all
        assert!(fit_polynomial(&[], &[], 0).is_none());
    }

    #[test]
    fn test_fit_polynomial_mismatched_lengths() {
        assert!(fit_polynomial(&[0.0, 1.0, 2.0], &[0.0, 1.0], 1).is_none());
    }

    #[test]
    fn test_fit_polynomial_zero_width_domain() {
        assert!(fit_polynomial(&[5.0, 5.0, 5.0], &[1.0, 2.0, 3.0], 1).is_none());
    }

    /// Two distinct x values can't constrain a parabola, however many points
    /// there are.
    #[test]
    fn test_fit_polynomial_singular() {
        let x = vec![0.0, 0.0, 1.0, 1.0];
        let y = vec![1.0, 1.0, 2.0, 2.0];
        assert!(fit_polynomial(&x, &y, 2).is_none());
        assert!(fit_polynomial(&x, &y, 1).is_some());
    }

    /// The fit is stable for large x values, like frequencies in Hz.
    #[test]
    fn test_fit_polynomial_high_values() {
        let f_min = 100e6;
        let f_step = 40e3;
        let x: Vec<f64> = (0..10).map(|i| f_min + i as f64 * f_step).collect();
        // y = x
        let y: Vec<f32> = x.iter().map(|&val| val as f32).collect();

        let coeffs = fit_polynomial(&x, &y, 1).expect("Fit failed");

        // Allow some error in the intercept due to cancellation
        assert!(coeffs[0].abs() < 1.0);
        assert!((coeffs[1] - 1.0).abs() < 1e-6);
    }
}

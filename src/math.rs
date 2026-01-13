use marlu::ndarray::Array2;

/// Fit a polynomial of order `poly_order` to the data `y` at coordinates `x`.
/// Returns the coefficients in the original basis, such that y approx sum(c_i * x^i).
/// Returns None if the fit is singular or there are not enough valid points.
///
/// The fit is performed in a normalized domain [-1, 1] for stability,
/// and then coefficients are transformed back to the original domain.
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

    // Normalize frequencies for stability during fit
    // x_norm = 2 * (f - f_min) / f_width - 1
    // x_norm = (f - f_mid) / f_half_width
    // Let x_val = (f - f_mid) / s where s = f_width / 2
    let f_min = x[0];
    let f_max = x[x.len() - 1];
    let f_width = f_max - f_min;
    
    // Avoid division by zero if all x are the same (unlikely for frequencies)
    if f_width.abs() < 1e-12 {
        return None;
    }
    
    let f_mid = f_min + f_width / 2.0;
    let s_scale = f_width / 2.0;

    let x_norm: Vec<f64> = valid_indices.iter().map(|&i| (x[i] - f_mid) / s_scale).collect();

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
    // Gaussian elimination
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

    // Convert coefficients from normalized basis to raw Hz basis
    // y = sum(beta_norm[j] * ((f - f_mid)/s)^j)
    // y = sum(beta_norm[j] * s^-j * sum(binom(j,k) * f^k * (-f_mid)^(j-k)))
    // coeff of f^k (c_k) is sum over j>=k of beta_norm[j] * s^-j * binom(j,k) * (-f_mid)^(j-k)

    let mut coeffs = vec![0.0; n];
    let mut binom = vec![vec![0.0; n]; n];
    
    // Precompute binomial coefficients
    for i in 0..n {
        binom[i][0] = 1.0;
        for j in 1..=i {
            binom[i][j] = binom[i-1][j-1] + binom[i-1][j];
        }
    }

    for k in 0..n {
        let mut sum = 0.0;
        for j in k..n {
            let term = beta_norm[j] * s_scale.powi(-(j as i32)) * binom[j][k] * (-f_mid).powi((j - k) as i32);
            sum += term;
        }
        coeffs[k] = sum as f32;
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
        let y: Vec<f32> = x.iter().map(|&val| (1.0 + 2.0 * val + 3.0 * val * val) as f32).collect();
        let order = 2;

        let coeffs = fit_polynomial(&x, &y, order).expect("Fit failed");
        
        assert_eq!(coeffs.len(), 3);
        assert!((coeffs[0] - 1.0).abs() < 1e-4);
        assert!((coeffs[1] - 2.0).abs() < 1e-4);
        assert!((coeffs[2] - 3.0).abs() < 1e-4);
    }

    #[test]
    fn test_fit_polynomial_with_nans() {
        // Fit y = 2x
        let x = vec![0.0, 1.0, 2.0, 3.0, 4.0];
        let mut y = vec![0.0, 2.0, f32::NAN, 6.0, 8.0];
        let order = 1;

        let coeffs = fit_polynomial(&x, &y, order).expect("Fit failed");
        
        assert_eq!(coeffs.len(), 2);
        assert!(coeffs[0].abs() < 1e-4); // intercept 0
        assert!((coeffs[1] - 2.0).abs() < 1e-4); // slope 2
    }

    #[test]
    fn test_fit_polynomial_not_enough_points() {
        let x = vec![0.0, 1.0];
        let y = vec![0.0, 1.0];
        let order = 2; // Need at least 3 points for order 2 (actually need order+1 points? No, order+1 params, so order+1 points)
        // fit_polynomial checks `if valid_indices.len() <= poly_order`. 
        // For order 2 (3 params), if len <= 2 (so 2 points), it returns None.
        
        let coeffs = fit_polynomial(&x, &y, order);
        assert!(coeffs.is_none());
    }
    
    #[test]
    fn test_fit_polynomial_high_values() {
        // Test with large x values (like frequencies in Hz)
        let f_min = 100e6;
        let f_step = 40e3;
        let n = 10;
        let x: Vec<f64> = (0..n).map(|i| f_min + i as f64 * f_step).collect();
        // y = x^2 (scaled down to fit in f32 range roughly)
        // Actually fit_polynomial returns f32 coeffs. 
        // If x ~ 1e8, x^2 ~ 1e16. This fits in f64 but fitting x^2 directly might have precision issues if we check f32 coeffs?
        // Let's try linear y = x
        let y: Vec<f32> = x.iter().map(|&val| val as f32).collect();
        let order = 1;
        
        let coeffs = fit_polynomial(&x, &y, order).expect("Fit failed");
        
        // y = 0 + 1*x
        assert!(coeffs[0].abs() < 1.0); // Allow some error due to large numbers cancellation
        assert!((coeffs[1] - 1.0).abs() < 1e-6);
    }
}

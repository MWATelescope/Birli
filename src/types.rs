#![allow(dead_code)]
use marlu::num_traits::{Float, Num};

pub use crate::{approx, Complex, Jones};

#[derive(Clone, Copy, Default, PartialEq)]
pub(crate) struct TestJones<F: Float + Num>(Jones<F>);

impl<F: Float> TestJones<F> {
    pub fn identity() -> TestJones<F> {
        TestJones(Jones::<F>::identity())
    }
}

impl<F: Float> From<Jones<F>> for TestJones<F> {
    #[inline]
    fn from(j: Jones<F>) -> Self {
        Self(j)
    }
}

impl From<Jones<f32>> for TestJones<f64> {
    #[inline]
    fn from(j_c32: Jones<f32>) -> Self {
        Self::from([
            Complex::new(j_c32[0].re as _, j_c32[0].im as _),
            Complex::new(j_c32[1].re as _, j_c32[1].im as _),
            Complex::new(j_c32[2].re as _, j_c32[2].im as _),
            Complex::new(j_c32[3].re as _, j_c32[3].im as _),
        ])
    }
}

impl From<Jones<f64>> for TestJones<f32> {
    #[inline]
    fn from(j_c32: Jones<f64>) -> Self {
        Self::from([
            Complex::new(j_c32[0].re as _, j_c32[0].im as _),
            Complex::new(j_c32[1].re as _, j_c32[1].im as _),
            Complex::new(j_c32[2].re as _, j_c32[2].im as _),
            Complex::new(j_c32[3].re as _, j_c32[3].im as _),
        ])
    }
}

impl<F: Float> From<[Complex<F>; 4]> for TestJones<F> {
    #[inline]
    fn from(j: [Complex<F>; 4]) -> Self {
        Self(Jones::from(j))
    }
}

impl<F: Float> std::ops::Deref for TestJones<F> {
    type Target = Jones<F>;

    #[inline]
    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl<F: Float> std::ops::DerefMut for TestJones<F> {
    #[inline]
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.0
    }
}

impl std::fmt::Display for TestJones<f32> {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(
            f,
            "[[{:e}{:+e}j, {:e}{:+e}j] [{:e}{:+e}j, {:e}{:+e}j]]",
            self[0].re,
            self[0].im,
            self[1].re,
            self[1].im,
            self[2].re,
            self[2].im,
            self[3].re,
            self[3].im,
        )
    }
}

impl std::fmt::Display for TestJones<f64> {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(
            f,
            "[[{:e}{:+e}j, {:e}{:+e}j] [{:e}{:+e}j, {:e}{:+e}j]]",
            self[0].re,
            self[0].im,
            self[1].re,
            self[1].im,
            self[2].re,
            self[2].im,
            self[3].re,
            self[3].im,
        )
    }
}

impl std::fmt::Debug for TestJones<f32> {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(
            f,
            "[[{:e}{:+e}j, {:e}{:+e}j] [{:e}{:+e}j, {:e}{:+e}j]]",
            self[0].re,
            self[0].im,
            self[1].re,
            self[1].im,
            self[2].re,
            self[2].im,
            self[3].re,
            self[3].im,
        )
    }
}

impl std::fmt::Debug for TestJones<f64> {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(
            f,
            "[[{:e}{:+e}j, {:e}{:+e}j] [{:e}{:+e}j, {:e}{:+e}j]]",
            self[0].re,
            self[0].im,
            self[1].re,
            self[1].im,
            self[2].re,
            self[2].im,
            self[3].re,
            self[3].im,
        )
    }
}

impl<F: Float + approx::AbsDiffEq> approx::AbsDiffEq for TestJones<F>
where
    F::Epsilon: Clone,
{
    type Epsilon = F::Epsilon;

    #[inline]
    fn default_epsilon() -> F::Epsilon {
        F::default_epsilon()
    }

    #[inline]
    fn abs_diff_eq(&self, other: &Self, epsilon: F::Epsilon) -> bool {
        (0..4).all(|idx| Complex::<F>::abs_diff_eq(&self[idx], &other[idx], epsilon.clone()))
    }
}

#![allow(dead_code)]
use marlu::num_traits::{Float, Num};

pub use crate::{approx, Complex, Jones};

#[derive(Clone, Copy, Default, PartialEq)]
pub(crate) struct TestJones<F: Float + Num>(Jones<F>);

impl<F: Float> TestJones<F> {
    pub fn identity() -> Self {
        Self(Jones::<F>::identity())
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
        Self::from(Jones::<f64>::from(j_c32))
    }
}

impl From<Jones<f64>> for TestJones<f32> {
    #[inline]
    fn from(j_c64: Jones<f64>) -> Self {
        Self::from(Jones::<f32>::from(j_c64))
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

#[cfg(test)]
mod tests {
    use marlu::Jones;

    use crate::types::TestJones;

    #[test]
    fn test_jones_debug_display() {
        let jones_c64 = Jones::<f64>::identity();
        let test_jones_c32: TestJones<f32> = TestJones::from(jones_c64);
        assert!(!format!("{:?}", test_jones_c32).is_empty());
        assert!(!format!("{}", test_jones_c32).is_empty());

        let jones_c32 = Jones::<f32>::identity();
        let test_jones_c64: TestJones<f64> = TestJones::from(jones_c32);
        assert!(!format!("{:?}", test_jones_c64).is_empty());
        assert!(!format!("{}", test_jones_c64).is_empty());
    }
}

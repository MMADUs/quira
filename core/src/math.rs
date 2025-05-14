use ndarray::{Array1, Array2, Array3};
use num_complex::Complex64;

/// mathematical type for complex number
pub type Complex = Complex64;

/// Mathematical structure for 1d tensor (Vector)
pub type Vector<T> = Array1<T>;

/// Mathematical structure for 2d tensor (Matrix)
pub type Matrix<T> = Array2<T>;

/// Mathematical structure for 3d tensor (Tensor)
pub type Tensor<T> = Array3<T>;

/// mathematical constant of pi = 3.14159...
pub const PI: f64 = std::f64::consts::PI;

/// Mathematical constant for π/2 (half of pi)
pub const PI_2: f64 = std::f64::consts::FRAC_PI_2;

/// Mathematical constant for π/4 (quarter of pi)
pub const PI_4: f64 = std::f64::consts::FRAC_PI_4;

/// Mathematical constant for 2π (tau, a full revolution in radians)
pub const TAU: f64 = std::f64::consts::TAU;

/// Mathematical constant for the square root of 2
pub const SQRT2: f64 = std::f64::consts::SQRT_2;

/// Mathematical constant for 1/√2 (inverse of the square root of 2)
pub const INV_SQRT2: f64 = std::f64::consts::FRAC_1_SQRT_2;

/// Mathematical constant for the natural logarithm of 2 (ln(2))
pub const LN2: f64 = std::f64::consts::LN_2;

/// Mathematical constant for Euler's number (e)
pub const E: f64 = std::f64::consts::E;

/// Smallest positive `f64` number that can be represented
/// This represents the smallest difference between two distinct `f64` values
pub const EPS: f64 = f64::EPSILON; 

/// Positive infinity in `f64`
/// This is used to represent infinite values in calculations
pub const INF: f64 = f64::INFINITY;

/// Negative infinity in `f64`
/// This represents negative infinite values in calculations
pub const NEG_INF: f64 = f64::NEG_INFINITY;

/// The imaginary unit (i), where i² = -1
/// Represented as a complex number with real part 0 and imaginary part 1
pub const I: Complex64 = Complex64::new(0.0, 1.0);

/// The complex number 1 + 0i (the real number 1 as a complex number)
pub const ONE: Complex64 = Complex64::new(1.0, 0.0);

/// The complex number 0 + 0i (the real number 0 as a complex number)
pub const ZERO: Complex64 = Complex64::new(0.0, 0.0);

/// The complex number 1/√2 + 0i (a scaled version of the inverse of sqrt(2))
pub const INV_SQRT2_C: Complex64 = Complex64::new(INV_SQRT2, 0.0);

// private sub-module defined in other files
mod array2d;
mod bounding_box;
mod point2d;
mod polynomial_regression_2d;

// exports identifiers from private sub-modules in the current module namespace
pub use self::array2d::Array2D;
pub use self::bounding_box::BoundingBox;
pub use self::point2d::Point2D;
pub use self::polynomial_regression_2d::PolynomialRegression2D;

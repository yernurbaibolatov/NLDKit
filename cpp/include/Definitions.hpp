#pragma once

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <complex>

// ===================================================
//               Type Aliases
// ===================================================

// Common Eigen types
using Vec         = Eigen::VectorXd;
using Mat         = Eigen::MatrixXd;
using RowVec      = Eigen::RowVectorXd;
using ColVec      = Eigen::VectorXd;  // alias for clarity
using Index       = Eigen::Index;     // Eigen’s default index type

// Fixed-size types
using Vec2        = Eigen::Vector2d;
using Vec3        = Eigen::Vector3d;
using Mat2        = Eigen::Matrix2d;
using Mat3        = Eigen::Matrix3d;

// Complex types
using Complex     = std::complex<double>;
using CVec        = Eigen::VectorXcd;
using CMat        = Eigen::MatrixXcd;

// Array types (for element-wise operations)
using Arr         = Eigen::ArrayXd;
using Arr2D       = Eigen::ArrayXXd;

// Sparse matrix (optional if needed)
using SparseMat   = Eigen::SparseMatrix<double>;

// ===================================================
//               Constants
// ===================================================

namespace Constants {
    // Math constants
    constexpr double PI         = 3.14159265358979323846;
    constexpr double TWO_PI     = 2.0 * PI;
    constexpr double HALF_PI    = 0.5 * PI;
    constexpr double INV_PI     = 1.0 / PI;
    constexpr double INV_2PI    = 1.0 / TWO_PI;
    constexpr double SQRT2      = 1.41421356237309504880;
    constexpr double SQRT3      = 1.73205080756887729352;

    // Degree <-> Radian conversions
    constexpr double DEG2RAD    = PI / 180.0;
    constexpr double RAD2DEG    = 180.0 / PI;

    // Numerical constants
    constexpr double EPSILON      = 1e-10;
    constexpr double TOLERANCE    = 1e-8;
    constexpr double LARGE_NUMBER = 1e10;
    constexpr double SMALL_NUMBER = 1e-10;

    // Simulation defaults (optional)
    constexpr double DEFAULT_DT   = 0.01;
    constexpr int    MAX_ITER     = 10000;
}

// No 'using namespace Constants;' here to avoid global pollution.
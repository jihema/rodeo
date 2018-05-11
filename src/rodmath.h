/*
 * rodmath.h
 */

#pragma once

#define EIGEN_INITIALIZE_MATRICES_BY_ZERO
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Geometry>

#include <stdexcept>
#include <iostream>

namespace rodeo
{

class MathException: public std::runtime_error
{
public:
    MathException(std::string const& what) :
            std::runtime_error(what)
    {
    }
};

using VecXd = Eigen::Matrix<double, Eigen::Dynamic, 1>;
using Vec3d = Eigen::Matrix<double, 3, 1>;
using MatXd = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>;
using DiagonalMatXd = Eigen::DiagonalMatrix<double, Eigen::Dynamic>;
using BandLimitedMatXd = MatXd;
// FIXME
using Mat3d = Eigen::Matrix<double, 3,3 >;
using Vec3d = Eigen::Matrix<double,3,1>;

inline Mat3d build_orthonormal_basis(Vec3d const& u)
{
    return Mat3d::Identity(); // TODO
}

inline Mat3d antisym(Vec3d const& x)
{
    Mat3d m;
    m << 0, -x[2], x[1], x[2], 0, -x[0], -x[1], x[0], 0;
    return m;
}

inline Mat3d parallel_transport_matrix(Vec3d const& a, Vec3d const& b)
{
    assert(a.squaredNorm() == 1 && b.squaredNorm() == 1);

    double const c = a.dot(b);
    if (c == 1) // a = b, no transport.
    {
        return Mat3d::Identity();
    }
    if (c == -1.) // a = -b, transport undefined.
    {
        throw MathException(__PRETTY_FUNCTION__);
    }

    Vec3d const v = a.cross(b);
    Mat3d const mv = antisym(v);
    return Mat3d::Identity() + mv + mv * mv / (1. + c);
}

inline Vec3d parallel_transport(Vec3d const& x, Vec3d const& a, Vec3d const& b)
{
    return parallel_transport_matrix(a, b) * x;
}

}

/**
 * @file ElementB3.cpp
 * Implementation of ElementB3 class (3-noded bar/membrane element).
 *
 * @author Haohang Huang
 * @date June 1, 2020
 */

#include "ElementB3.h"
#include <cmath>

ElementB3::staticMembers ElementB3::statics(3,3,1,3,3);

ElementB3::ElementB3()
{
}

ElementB3::ElementB3(const int & index, const std::vector<int> & nodeList, Node** const meshNode, Material* const material)
 : Element(index, nodeList, meshNode, material) // call the constructor of base class in the initializer list!
{
    modulusAtGaussPt = VectorXd::Constant(statics.shape->gaussianPt().size(), material_->modulus()); // 3 x 1 vector
}

ElementB3::~ElementB3()
{
}

Shape* ElementB3::shape() const
{
 return statics.shape;
}


MatrixXd ElementB3::EMatrix(const VectorXd & modulus) const
{
    // Note: to keep consistent with LinearElastic and Nonlinear Elastic scheme, here we still pass in an void modulus variable, but ignore it.
    (void)modulus;
    return material_->EMatrix();
}

MatrixXd ElementB3::BMatrix(const Vector2d & point) const
{
    MatrixXd B = MatrixXd::Zero(2, 2 * size_); // 2x6, B matrix
    // [cos * dN1/dxi  sin * dN1/dxi | cos * dN2/dxi  sin * dN2/dxi | cos * dN3/dxi  sin * dN3/dxi ]
    // [N1/r           0             | N2/r           0             | N3             0             ]

    // here B3 element is different from Q8 element
    // strain and stress of Q8 element is in global coordinates (r-z), but for B3 element they are in local coordinates (axial-hoop)
    // so we can just use local derivatives dN/dxi
    const VectorXd & N = shape()->functionVec(point); // N
    const MatrixXd & localDeriv = shape()->functionDeriv(point); // dN/dxi

    for (int n = 0; n < size_; n++) {
        B(0, 2 * n) = std::cos(_angle()) * localDeriv(0, n); // cos * dNi/dxi
        B(0, 2 * n + 1) = std::sin(_angle()) * localDeriv(0, n); // sin * dNi/dxi
        B(1, 2 * n) = N(n) / radius(point); // Ni/r 
    }

    return B;
}

MatrixXd ElementB3::_BMatrix(const int & i) const
{
    MatrixXd B = MatrixXd::Zero(2, 2 * size_); // 2x6, B matrix
    // [cos * dN1/dxi  sin * dN1/dxi | cos * dN2/dxi  sin * dN2/dxi | cos * dN3/dxi  sin * dN3/dxi ]
    // [N1/r           0             | N2/r           0             | N3             0             ]

    // here B3 element is different from Q8 element
    // strain and stress of Q8 element is in global coordinates (r-z), but for B3 element they are in local coordinates (axial-hoop)
    // so we can just use local derivatives dN/dxi
    const VectorXd & N = shape()->functionVec(i); // N
    const MatrixXd & localDeriv = shape()->functionDeriv(i); // dN/dxi

    for (int n = 0; n < size_; n++) {
        B(0, 2 * n) = std::cos(_angle()) * localDeriv(0, n); // cos * dNi/dxi
        B(0, 2 * n + 1) = std::sin(_angle()) * localDeriv(0, n); // sin * dNi/dxi
        B(1, 2 * n) = N(n) / _radius(i); // Ni/r 
    }

    return B;
}

double ElementB3::_jacobianDet(const int & i) const
{   
    // for membrane element, |J| = L/2
    (void)i;
    return _length() / 2;
}

double ElementB3::_angle() const
{
    // angle = arctan2(z2-z0/r2-r0) in [-pi, pi], radians
    double alpha = std::atan2(nodeCoord_(2,1) - nodeCoord_(0,1), nodeCoord_(2,0) - nodeCoord_(0,0));  
    return alpha;
}

double ElementB3::_length() const
{
    return std::sqrt( std::pow(nodeCoord_(2,1) - nodeCoord_(0,1), 2) + std::pow(nodeCoord_(2,0) - nodeCoord_(0,0), 2) );
}
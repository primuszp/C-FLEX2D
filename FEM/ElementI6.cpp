/**
 * @file ElementI6.cpp
 * Implementation of ElementI6 class (6-noded interface element).
 *
 * @author Haohang Huang
 * @date Dec 31, 2020
 */

// for use of M_PI
#define _USE_MATH_DEFINES
#include <cmath>
#include "ElementI6.h"

ElementI6::ElementI6()
{
}

ElementI6::ElementI6(const int & index, const std::vector<int> & nodeList, Node** const meshNode, Material* const material)
 : Element(index, nodeList, meshNode, material) // call the constructor of base class in the initializer list!
{
    // B matrix never changes, so cache it as a member variable
    B_ = MatrixXd::Zero(6, 2 * size_); // 6x12, B matrix
    // from nodal displacements to relative displacements
    
    for (int n = 0; n < size_; n++) {
        B_(n, n) = -1;
        B_(n, n + 6) = 1;
    }
}

ElementI6::~ElementI6()
{
}

Shape* ElementI6::shape() const
{   // not used but a pure virtual method, return empty to suppress error
    return NULL;
}

MatrixXd ElementI6::EMatrix(const VectorXd & modulus) const
{
    // Note: to keep consistent with LinearElastic and Nonlinear Elastic scheme, here we still pass in an void modulus variable, but ignore it.
    (void)modulus;
    
    double r_avg = (nodeCoord_(0,0) + nodeCoord_(1,0) + nodeCoord_(2,0)) / 3;
    double L = _length(); 
    double alpha = _angle();
    
    double c0 = M_PI * L / 3 * (r_avg - L/2 * std::cos(alpha));
    double c1 = M_PI * L * 4/3 * r_avg;
    double c2 = M_PI * L / 3 * (r_avg + L/2 * std::cos(alpha));
    
    double ks = material_->getInterfaceShearStiffness();
    double kn = material_->getInterfaceNormalStiffness();

    MatrixXd E = MatrixXd::Zero(6, 6); // 6x6
    E(0,0) = c0 * ks;
    E(1,1) = c0 * kn;
    E(2,2) = c1 * ks;
    E(3,3) = c1 * kn;
    E(4,4) = c2 * ks;
    E(5,5) = c2 * kn;

    return E;
}

MatrixXd ElementI6::BMatrix(const Vector2d & point) const
{
    (void)point;
    return B_; // initialized in constructor
}

MatrixXd ElementI6::_BMatrix(const int & i) const
{   // not used but a pure virtual method, return empty to suppress error
    (void)i;
    return Matrix2d::Zero();
}

double ElementI6::_angle() const
{
    // angle = arctan2(z2-z0/r2-r0) in [-pi, pi], radians
    double alpha = std::atan2(nodeCoord_(2,1) - nodeCoord_(0,1), nodeCoord_(2,0) - nodeCoord_(0,0));  
    return alpha;
}

double ElementI6::_length() const
{
    // length = sqrt((z2-z0)^2 + (r2-r0)^2)
    return std::sqrt( std::pow(nodeCoord_(2,1) - nodeCoord_(0,1), 2) + std::pow(nodeCoord_(2,0) - nodeCoord_(0,0), 2) );
}
/**
 * @file ElementQ8.cpp
 * Implementation of ElementQ8 class.
 *
 * @author Haohang Huang
 * @date Feburary 13, 2018
 */

#include "ElementQ8.h"

// static member should be initialized outside the class body because it doesn't
// depend on any instance of that class and is initialized before any instance is
// created. Initialization here does not require the static member to be public.
ElementQ8::staticMembers ElementQ8::statics(8,9,4,3,3);

ElementQ8::ElementQ8()
{
}

ElementQ8::ElementQ8(const int & index, const std::vector<int> & nodeList, Node** const meshNode, Material* const material)
    : Element(index, nodeList, meshNode, material) // call the constructor of base class in the initializer list!
{
    if (!material->anisotropy) {
        modulusAtGaussPt = VectorXd::Constant(statics.shape->gaussianPt().size(), material_->modulus()); // 9 x 1 vector, the length depends on element type, so can only be initialized in derived class
    } else {
        modulusAtGaussPt = MatrixXd(statics.shape->gaussianPt().size(), 3); // 9 x 3 Matrix. Each column: horizontal modulus, vertical modulus, shear modulus
        modulusAtGaussPt.col(0) = VectorXd::Constant(statics.shape->gaussianPt().size(), material_->modulusR());
        modulusAtGaussPt.col(1) = VectorXd::Constant(statics.shape->gaussianPt().size(), material_->modulusZ());
        modulusAtGaussPt.col(2) = VectorXd::Constant(statics.shape->gaussianPt().size(), material_->modulusG());
    }
}

ElementQ8::~ElementQ8()
{
}

Shape* ElementQ8::shape() const
{
    return statics.shape;
}

MatrixXd ElementQ8::EMatrix(const VectorXd & modulus) const
{
    if (!material_->nonlinearity)
        return material_->EMatrix();
    else
        return material_->EMatrix(modulus);
}

MatrixXd ElementQ8::BMatrix(const Vector2d & point) const
{
    MatrixXd B = MatrixXd::Zero(4, 2 * size_);
    MatrixXd globalDeriv = (shape()->functionDeriv(point) * nodeCoord_).inverse() * shape()->functionDeriv(point);

    for (int n = 0; n < size_; n++) {
        B(0, 2 * n) = globalDeriv(0, n);
        VectorXd N = shape()->functionVec(point);
        B(1, 2 * n) = N(n) / radius(point);
        B(2, 2 * n + 1) = globalDeriv(1, n);
        B(3, 2 * n) = globalDeriv(1, n);
        B(3, 2 * n + 1) = globalDeriv(0, n);
    }

    return B;
}

MatrixXd ElementQ8::_BMatrix(const int & i) const
{
    // Example dimensions are given for ElementQ8 type
    MatrixXd B = MatrixXd::Zero(4, 2 * size_); // 4x16, B matrix
    // [dN1/dr  0      | dN2/dr  0      | ... ]
    // [N1/r    0      | N2/r    0      | ... ]
    // [0       dN1/dz | 0       dN2/dz | ... ]
    // [dN1/dz  dN1/dr | dN2/dz  dN2/dr | ... ]
    // where dN/dr = J^-1 * dN/dxi, dN/dz = J^-1 * dN/deta

    // 2x8 global derivatives [dN/dr; dN/dz] = 2x2 inversed Jacobian [J^-1] * 2x8 local derivatives [dN/dxi; dN/deta]
    // where 2x2 inversed Jacobian = 2x8 local derivatives [dN/dxi; dN/deta] * 8x2 node coordinates [ri zi]
    MatrixXd globalDeriv = (shape()->functionDeriv(i) * nodeCoord_).inverse() * shape()->functionDeriv(i);

    for (int n = 0; n < size_; n++) {
        B(0, 2 * n) = globalDeriv(0, n); // dNi/dr
        const VectorXd & N = shape()->functionVec(i);
        B(1, 2 * n) = N(n) / _radius(i); // Ni/r // (shape()->functionVec(i))(n);
        B(2, 2 * n + 1) = globalDeriv(1, n); // dNi/dz
        B(3, 2 * n) = globalDeriv(1, n); // dNi/dz
        B(3, 2 * n + 1) = globalDeriv(0, n); // dNi/dr
    }

    return B;
}


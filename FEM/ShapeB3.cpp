/**
 * @file ShapeB3.cpp
 * Implementation of ShapeB3 class.
 *
 * @author Haohang Huang
 * @date June 1, 2020
 */

#include "ShapeB3.h"
#include <cmath>

ShapeB3::ShapeB3(const int & nodes, const int & gaussians, const int & edges, const int & edgeNodes, const int & edgeGaussians) :
    Shape(nodes, gaussians, edges, edgeNodes, edgeGaussians) // Call the constructor of base class
{
    // Set up shape parameters
    // Local xi coordinates of nodes
    // Note: to be consistent w/ the signature of functionVec/Deriv etc, we still make the nodeCoord and gaussianPt Vector2d by assign 0 for eta coordinate
    nodeCoord_[0] << -1, 0;
    nodeCoord_[1] << 0, 0;
    nodeCoord_[2] << 1, 0;

    // Local xi coordinates of gaussian points
    double temp = std::sqrt(0.6);
    gaussianPt_[0] << -temp, 0;
    gaussianPt_[1] << 0, 0;
    gaussianPt_[2] << temp, 0;

    // Edge Gaussian points
    edgeGaussianPt_[0] = -temp;
    edgeGaussianPt_[1] = 0;
    edgeGaussianPt_[2] = temp;

    // Gaussian weights
    double corner = 5.0 / 9.0;
    double side = 8.0 / 9.0;
    gaussianWt_[0] = corner;
    gaussianWt_[1] = side;
    gaussianWt_[2] = corner;

    // Edge Gaussian weights
    edgeGaussianWt_[0] = corner;
    edgeGaussianWt_[1] = side;
    edgeGaussianWt_[2] = corner;

    // Edge list
    std::vector<int> edge1{0,1,2};
    edgeList_[0] = edge1;

    // After setting up the above parameters, pre-cache the shape function
    _cacheShape();

}

ShapeB3::~ShapeB3()
{
}

VectorXd ShapeB3::functionVec(const Vector2d & point) const
{ // 3x1 Vector
    VectorXd result(numNodes_);
    result(0) = point(0) * (point(0) - 1) / 2;
    result(1) = 1 - point(0) * point(0);
    result(2) = point(0) * (point(0) + 1) / 2;

    return result;
}

MatrixXd ShapeB3::functionMat(const Vector2d & point) const
{ // 2x6 interleaved matrix
    // Note: to use different shape function for u and v direction, this can be modified
    MatrixXd result = MatrixXd::Zero(2, 2 * numNodes_); // 2-D, 3-Node
    
    // u direction, bar element shape functions
    double N0_u = point(0) * (point(0) - 1) / 2;
    double N1_u = 1 - point(0) * point(0);
    double N2_u = point(0) * (point(0) + 1) / 2;
    
    // v direction, beam element shape functions (see Jay Kwon 2007)
    double N0_v = (32 * std::pow(point(0),2) - 20 * std::pow(point(0),3) - 4 * std::pow(point(0),4) - 3 * std::pow(point(0),5)) / 128;
    double N1_v = (128 - 64 * std::pow(point(0),2) + 8 * std::pow(point(0),4)) / 128;
    double N2_v = (32 * std::pow(point(0),2) - 20 * std::pow(point(0),3) - 4 * std::pow(point(0),4) - 3 * std::pow(point(0),5)) / 128;
    
    result(0, 0) = N0_u;
    result(1, 1) = N0_v;
    result(0, 2) = N1_u;
    result(1, 3) = N1_v;
    result(0, 4) = N2_u;
    result(1, 5) = N2_v;

    return result;
}

MatrixXd ShapeB3::functionDeriv(const Vector2d & point) const
{ // 1x3 matrix (here we only have one xi coordinate)
    MatrixXd result(1, numNodes_);
    result(0, 0) = (2 * point(0) - 1) / 2; // dNi/dxi = (2xi-1)/2
    result(0, 1) = - 2 * point(0); // dNi/dxi = -2xi
    result(0, 2) = (2 * point(0) + 1) / 2; // dNi/dxi = (2xi+1)/2

    return result;
}

VectorXd ShapeB3::edgeFunctionVec(const double & point) const
{ // 3x1 vector
    VectorXd result(numEdgeNodes_);
    double Ni = 0;
    for (int i = 0; i < 3; i++) {
        if (i == 0) // left/bottom point
          Ni = - point * (1 - point) / 2; // Ni = -x(1-x)/2
        if (i == 2) // right/top point
          Ni = point * (1 + point) / 2; // Ni = x(1+x)/2
        if (i == 1) // mid-side point
          Ni = 1 - point * point; // Ni = (1-x^2)
        result(i) = Ni;
    }
    return result;
}

MatrixXd ShapeB3::edgeFunctionMat(const double & point) const
{ // 2x6 matrix, 3 gaussian point and 3 node at each edge
    MatrixXd result = MatrixXd::Zero(2, 2 * numEdgeNodes_);
    double Ni = 0;
    for (int i = 0; i < 3; i++) {
        if (i == 0) // left/bottom point
          Ni = - point * (1 - point) / 2; // Ni = -x(1-x)/2
        if (i == 2) // right/top point
          Ni = point * (1 + point) / 2; // Ni = x(1+x)/2
        if (i == 1) // mid-side point
          Ni = 1 - point * point; // Ni = (1-x^2)
        result(0, 2 * i) = Ni;
        result(1, 2 * i + 1) = Ni;
    }
    return result;
}

VectorXd ShapeB3::edgeFunctionDeriv(const double & point) const
{ // 3x1 vector
    VectorXd result(numEdgeNodes_);
    double Ni = 0;
    for (int i = 0; i < 3; i++) {
        if (i == 0) // left/bottom point
          Ni = (2 * point - 1) / 2; // dNi/dx = (2x-1)/2
        if (i == 2) // right/top point
          Ni = (2 * point + 1) / 2; // dNi/dx = (2x+1)/2
        if (i == 1) // mid-side point
          Ni = - 2 * point; // dNi/dx = -2x
        result(i) = Ni;
    }
    return result;
}

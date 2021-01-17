/**
 * @file ShapeB3.h
 * Derived class from Shape for the shape information of isoparametric B3 element.
 *
 * @author Haohang Huang
 * @date June 1, 2020
 */

#ifndef ShapeB3_h
#define ShapeB3_h

#include "Shape.h"

/* Derived class for the shape of isoparametric B3 element.
 * The sketch and index of the B3 element is:
 *
 * 0 -- 1 -- 2
 *
 * The sketch and index of the Gaussian integration points is:
 *
 * 0 -- 1 -- 2
 */
class ShapeB3 : public Shape
{
    public:
        /* See the documentation of base class Shape. */
        ShapeB3(const int & nodes, const int & gaussians, const int & edges, const int & edgeNodes, const int & edgeGaussians);
        ~ShapeB3();

        VectorXd functionVec(const Vector2d & point) const;
        MatrixXd functionMat(const Vector2d & point) const;
        MatrixXd functionDeriv(const Vector2d & point) const;

        VectorXd edgeFunctionVec(const double & point) const;
        MatrixXd edgeFunctionMat(const double & point) const;
        VectorXd edgeFunctionDeriv(const double & point) const;

};

#endif /* ShapeB3_h */

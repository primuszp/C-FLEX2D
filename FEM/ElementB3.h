/**
 * @file ElementB3.h
 * Derived class from Element for the isoparametric B3 element.
 *
 * @author Haohang Huang
 * @date June 1, 2020
 */

#ifndef ElementB3_h
#define ElementB3_h

#include "Element.h"

/* Derived class for the isoparametric B3 element.
*/
class ElementB3 : public Element
{
 public:
     /* See the documentation of base class Element. */
     ElementB3();
     ElementB3(const int & index, const std::vector<int> & nodeList, Node** const meshNode, Material* const material);
     ~ElementB3();

     Shape* shape() const;

     MatrixXd EMatrix(const VectorXd & modulus) const;
     MatrixXd BMatrix(const Vector2d & point) const;
     MatrixXd _BMatrix(const int & i) const;
     double _jacobianDet(const int & i) const;

 private:
     /** A static structure that manages all the static members used in this class */
     static staticMembers statics;

     /** For bar element, calculate the orientation of the element */
     double _angle() const;

     /** For bar element, calculate the element length */
     double _length() const;

};

#endif /* ElementB3_h */

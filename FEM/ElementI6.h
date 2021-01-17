/**
 * @file ElementI6.h
 * Derived class from Element for the I6 interface element.
 *
 * @author Haohang Huang
 * @date Dec 31, 2020
 */

#ifndef ElementI6_h
#define ElementI6_h

#include "Element.h"

/* Derived class for the interface I6 element.
*/
class ElementI6 : public Element
{
 public:
     /* See the documentation of base class Element. */
     ElementI6();
     ElementI6(const int & index, const std::vector<int> & nodeList, Node** const meshNode, Material* const material);
     ~ElementI6();

     Shape* shape() const;

     MatrixXd EMatrix(const VectorXd & modulus) const;
     MatrixXd BMatrix(const Vector2d & point) const;
     MatrixXd _BMatrix(const int & i) const;

 private:
     
     /** cached B matrix */
     MatrixXd B_;

     /** For interface element, calculate the orientation of the element */
     double _angle() const;

     /** For interface element, calculate the element length */
     double _length() const;

};

#endif /* ElementI6_h */

/**
 * @file ElementQ8.h
 * Derived class from Element for the isoparametric Q8 element.
 *
 * @author Haohang Huang
 * @date Feburary 13, 2018
 * @note Efficiency optimized by polymorph shape on March 26, 2018.
 * @note Efficiency optimized by storing local stiffness matrix and return-by-ref
 * on March 27, 2018
 * @note Efficiency optimized by the generalization of all element-wise operations
 * into base class Element on Apr 22, 2018.
 */

#ifndef ElementQ8_h
#define ElementQ8_h

#include "Element.h"

/* Derived class for the isoparametric Q8 element.
 */
class ElementQ8 : public Element
{
    public:
        /* See the documentation of base class Element. */
        ElementQ8();
        ElementQ8(const int & index, const std::vector<int> & nodeList, Node** const meshNode, Material* const material);
        ~ElementQ8();

        Shape* shape() const;

        MatrixXd EMatrix(const VectorXd & modulus) const;
        MatrixXd BMatrix(const Vector2d & point) const;
        MatrixXd _BMatrix(const int & i) const;

    private:
        /** A static structure that manages all the static members used in this class */
        static staticMembers statics;

};

#endif /* ElementQ8_h */

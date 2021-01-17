/**
 * @file Geosynthetic.h
 * Derived class from Material for geosynthetic materials.
 *
 * @author Haohang Huang
 * @date Dec 31, 2020
 */

#ifndef Geosynthetic_h
#define Geosynthetic_h

#include "Material.h"

class Geosynthetic : public Material
{
  public:
    /** See the documentation of base class Material. */
    Geosynthetic(const bool & anisotropy, const bool & nonlinearity, const bool & noTension, const bool & geosynthetic, const std::vector<double> & properties);
    ~Geosynthetic();
  
    // override
    double getInterfaceShearStiffness() const;
    double getInterfaceNormalStiffness() const;

  protected:
    /** Geosynthetic material properties is different from the base defintion
      * It shares the M_, v_, E_ members as the base class, and defines thickness t_, spring stiffness ks_ and kn_ 
      * @note ks and kn properties are used for I6 interface elements as well
      */
    /** membrane thickness */  
    double t_;
    
    /** shear stiffness */  
    double ks_;

    /** normal stiffness */  
    double kn_;
};

#endif /* Geosynthetic_h */

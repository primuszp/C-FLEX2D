/**
 * @file Geosynthetic.cpp
 * Implementation of Geosynthetic class.
 *
 * @author Haohang Huang
 * @date Dec 31, 2020
 */

#include "Geosynthetic.h"

Geosynthetic::Geosynthetic(const bool & anisotropy, const bool & nonlinearity, const bool & noTension, const bool & geosynthetic, const std::vector<double> & properties)
    : Material(anisotropy, nonlinearity, noTension, geosynthetic)
{   
    int i = 0;
    double M = properties[i++]; M_ = M; // initialize members of base class Material
    double v = properties[i++]; v_ = v;
    double t = properties[i++]; t_ = t;
    double ks = properties[i++]; ks_ = ks;
    double kn = properties[i++]; kn_ = kn;

    E_ = MatrixXd::Zero(2,2);
    E_ << 1, v,
          v, 1;
    E_ = E_ * M * t / (1 - v * v);

    // not used. keep dimension consistent in intergration
    bodyForce_ << 0,0;
    thermalStrain_ = VectorXd::Zero(2);
    thermalStrain_ << 0,0;
}

Geosynthetic::~Geosynthetic()
{
}

double Geosynthetic::getInterfaceShearStiffness() const
{
    return ks_;
}

double Geosynthetic::getInterfaceNormalStiffness() const
{
    return kn_;
}
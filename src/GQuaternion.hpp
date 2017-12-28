//
//  GQuaternion.hpp
//  GFC
//
//  Created by 李桢 on 2017/12/28.
//

#ifndef GQuaternion_hpp
#define GQuaternion_hpp

#include <stdio.h>
#include "GVector.hpp"



namespace gfc
{
    
    // ref: http://web.cs.iastate.edu/~cs577/handouts/quaternion.pdf
    
    class GQuaternion
    {
    public:
        
        GQuaternion();
        GQuaternion(double in_r, double in_i, double in_j, double in_k);
        GQuaternion(double real,  GVector vector);
        
        void set(double half_angle, GVector u);
        
        GVector rotate(GVector a);
        
        // Operator methods
        void conjugate(); // conjugate a quaternion
        void normalise();
        void inverse();  // inverse a quaternion
        double norm();
        double norm2();
        
        GQuaternion operator-() const;  //!< Negation
        void operator+=(GQuaternion a); //!< Addition
        void operator-=(GQuaternion a); //!< Subtraction
        
        void operator*=(GQuaternion a); //!< Hamilton product
        void operator/=(GQuaternion a); //!< Hamilton product
        void operator*=(GVector a);  //!< Hamilton product
        
        void operator*=(double k); //!< Scalar multiplication
        void operator/=(double k); //!< Scalar division
        
        
    private:
        double r = 1.0;
        GVector v;
        
    };
    
    
    GQuaternion operator+(GQuaternion a, GQuaternion b); //!< sum
    GQuaternion operator*(GQuaternion a, GQuaternion b); //!< sum
    GQuaternion operator-(GQuaternion a, GQuaternion b); //!< sum
    
    GQuaternion operator*(GQuaternion a, GVector b);  //!< Hamilton product
    GQuaternion operator*(GVector a, GQuaternion b);  //!< Hamilton product
    GQuaternion operator/(GQuaternion a, GQuaternion b);  //!< Hamilton product
    
}




#endif /* GQuaternion_hpp */

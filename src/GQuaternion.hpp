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
    class GQuaternion
    {
    public:
        
        GQuaternion();
        GQuaternion(double in_r, double in_i, double in_j, double in_k);
        GQuaternion(double real,  GVector vector);
        // Operator methods
        GQuaternion operator-() const;  //!< Negation
        void operator+=(GQuaternion a); //!< Addition
        void operator-=(GQuaternion a); //!< Subtraction
        
        void operator*=(GQuaternion a); //!< Hamilton product
        void operator*=(GVector a);  //!< Hamilton product
        
        void operator*=(double k); //!< Scalar multiplication
        void operator/=(double k); //!< Scalar division
        
    private:
        double r = 1.0;
        GVector v;
        
    };

}




#endif /* GQuaternion_hpp */

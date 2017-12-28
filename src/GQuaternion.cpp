//
//  GQuaternion.cpp
//  GFC
//
//  Created by 李桢 on 2017/12/28.
//

#include "GQuaternion.hpp"
namespace gfc
{
    GQuaternion::GQuaternion()
    {
        r = 0.0;
    }
    
    GQuaternion::GQuaternion(double in_r, double in_i, double in_j, double in_k)
    {
        r = in_r;
        v.set(in_i, in_j, in_k);
    }
    
    GQuaternion::GQuaternion(double real, GVector vector)
    {
        r = real;
        v = vector;
    }
    
    
    
}

//
//  GVector.cpp
//  GFC
//
//  Created by lizhen on 17/06/2016.
//  Copyright © 2016 lizhen. All rights reserved.
//


//============================================================================
//
//  This file is part of GFC, the GNSS FOUNDATION CLASS.
//
//  The GFC is free software; you can redistribute it and/or modify
//  it under the terms of the GNU Lesser General Public License as published
//  by the Free Software Foundation; either version 3.0 of the License, or
//  any later version.
//
//  The GFC is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public
//  License along with GFC; if not, write to the Free Software Foundation,
//  Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110, USA
//
//  Copyright 2015, lizhen
//
//============================================================================

#include "GVector.hpp"
namespace gfc
{
    // Default constructor
    GVector::GVector()
    {
        set(0.0, 0.0, 0.0);
    }
    
    // Three-argument constructor
    GVector::GVector(double in_x, double in_y, double in_z)
    {
        set(in_x, in_y, in_z);
    }
    
    void GVector::set(double in_x, double in_y, double in_z)
    { //assign specific values
        x = in_x;
        y = in_y;
        z = in_z;
    }
    
    GVector GVector::operator-() const
    {
        return GVector(-x, -y, -z);
    }
    
    void GVector::operator+=(GVector a)
    {
        x += a.x;
        y += a.y;
        z += a.z;
    }
    
    void GVector::operator-=(GVector a)
    {
        x -= a.x;
        y -= a.y;
        z -= a.z;
    }
    
    void GVector::operator-=(double k)
    {
        x -= k;
        y -= k;
        z -= k;
    }
    
    
    // Scalar multiply
    void GVector::operator*=(double k)
    {
        x *= k;
        y *= k;
        z *= k;
    }
    
    // Scalar divide
    void GVector::operator/=(double k)
    {
        x /= k;
        y /= k;
        z /= k;
    }
    
    double GVector::norm2() const
    {
        return x * x + y * y + z * z;
    }
    
    double GVector::norm() const
    {
        return sqrt(norm2());
    }
    
    // Method for scaling components so that the vector has magnitude 1.
    void GVector::normalise()
    {
        double s = norm();
        if (s > 0)
        {
            x /= s;
            y /= s;
            z /= s;
        }
        else
        {
            x = 0;
            y = 0;
            z = 0;
        }
        
//        if(fabs(x) <1.0E-16)
//        {
//            x = 0.0;
//        }
//        if(fabs(y) <1.0E-16)
//        {
//            y = 0.0;
//        }
//        if(fabs(z) <1.0E-16)
//        {
//            z = 0.0;
//        }
        
    }
    
    
    // Method for scaling components so that the vector has magnitude 1.
    GVector normalise(GVector a)
    {
        a.normalise();
        return a;
    }
    
    GVector operator+(GVector a, GVector b) //add two GVectors
    {
        a += b;
        return a;
    }
    
    GVector operator-(GVector a, GVector b) //subtract two GVectors
    {
        a -= b;
        return a;
    }
    
    GVector operator- (GVector a, double k)
    {
        a-=k;
        return a;
    }
    
    GVector operator0(double k, GVector a)
    {
        a-= k;
        return a;
    }
    
    //scalar multiply
    GVector operator*(GVector a, double b)
    {
        a *= b;
        return a;
    }
    
    //scalar multiply
    GVector operator*(double b, GVector a)
    {
        a *= b;
        return a;
    }
    
    //scalar divide
    GVector operator/(GVector a, double b)
    {
        a /= b;
        return a;
    }
    
    GVector crossproduct(const GVector &a, const GVector &b)
    {
        GVector result;
        result.x = a.y * b.z - a.z * b.y;
        result.y = a.z * b.x - a.x * b.z;
        result.z = a.x * b.y - a.y * b.x;
        
        return result;
    } //end of function crossproduct
    
    double dotproduct(const GVector &a, const GVector &b)
    {
        return a.x * b.x + a.y * b.y + a.z * b.z;
    } //end of function dotproduct

    
    
    
}  // end of namespace

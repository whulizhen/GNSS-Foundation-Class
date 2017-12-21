//
//  GVector.hpp
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


#ifndef GVector_hpp
#define GVector_hpp

#include <stdio.h>
#include <math.h>

namespace gfc
{
    class GVector
    {
    public:
        // Attributes
        double x;
        double y;
        double z;
        
        // Methods
        GVector(); //!< default constructor
        GVector(double in_x, double in_y, double in_z);
        
        
        void set(double in_x, double in_y, double in_z); //!< assignment components
        
        // Operator methods
        GVector operator-() const;  //!< negation
        void operator+=(GVector a); //!< addition
        void operator-=(GVector a); //!< subtraction
        void operator*=(double k);    //!< scalar multiplication
        void operator/=(double k);    //!< scalar division
        void operator-=(double k);
        
        double norm2() const; //!< length of vector squared
        double norm() const;  //!< length of vector
        
        void normalise(); //!< convert to unit vector.
        
        //void print() const; //!< print to screen
        
    };

    
    GVector normalise(GVector a);
    
    GVector operator+(GVector a, GVector b); //!< returns a+b
    GVector operator-(GVector a, GVector b); //!< returns a-b
    GVector operator- (GVector a, double k);
    GVector operator- (double k, GVector a);
    
    GVector operator*(GVector a, double b);
    GVector operator*(double b, GVector a);
    GVector operator/(GVector a, double b);
    
    GVector crossproduct(const GVector &a, const GVector &b);
    double dotproduct(const GVector &a, const GVector &b);
    
    
    
}  // end of namespace


#endif /* GVector_hpp */

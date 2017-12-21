//
//  GEquation.hpp
//  GFC
//
//  Created by lizhen on 28/08/2016.
//  Copyright © 2016 lizhen. All rights reserved.
//

#ifndef GEquation_hpp
#define GEquation_hpp

#include <stdio.h>

#include "GMatrix.h"

#include "GString.h"

namespace gfc
{
    
    // the variable struct,including a name and value
    struct GVariable
    {
        GString name;
        long double value;
        GVariable()
        {
            value = 0.0;
            name = "UNSET";
        }
    };
    
    // the measurement equation
    // mainly describle the relation between unknowns and observations
    struct GEquation
    {
        // the observation matrix, usually it is a column vector
        GMatrix L;
        
        GMatrix B;  // the design matrix describling the relation between unknows and observation
        
        GMatrix D;  // the covariance of observations
        
        std::vector< GString > strL;  // the description of observation
        
        std::vector< GString > strX;  // the description of unknowns
        
    };
    
    
    // the equation system, manage all the observations and construct the relations between observations and unknowns
    class GEquationSystem
    {
        
    public:
        
        
        
    };
    
    
    
    
} // end of namespace gfc




#endif /* GEquation_hpp */

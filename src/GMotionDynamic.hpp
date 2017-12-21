//
//  GMotionDynamic.hpp
//  GFC
//
//  Created by lizhen on 22/07/2016.
//  Copyright © 2016 lizhen. All rights reserved.
//

#ifndef GMotionDynamic_hpp
#define GMotionDynamic_hpp

#include <stdio.h>

namespace gfc
{
    // this class must be the father class of classes which use integrator
    class GMotionDynamic
    {
        
    public:
        
       virtual void  getDerivatives( int n, double x, double *y, double *dydx )
        {
            
        }
        
       virtual void collectStateInformation()
        {
            
        }
        
    };
    
    
}


#endif /* GMotionDynamic_hpp */

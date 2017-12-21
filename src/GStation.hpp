//
//  GStation.hpp
//  GFC
//
//  Created by lizhen on 20/08/2016.
//  Copyright © 2016 lizhen. All rights reserved.
//

#ifndef GStation_hpp
#define GStation_hpp

#include <stdio.h>

#include <list>

#include "GSourceID.h"

#include "GOBSData.h"

#include <set>

namespace gfc
{
    
    /*
     the antenna class
     storing the antenna pco and pcv information
    */
    class GAntenna
    {
        
    };
    
    
    class GStation : public GSourceID
    {
        
        
    public:
        
        
    private:
        
        // the position of the station
        double m_x;
        double m_y;
        double m_z;
        
        //the velocity of the station
        double m_vx;
        double m_vy;
        double m_vz;
        
        //std::list<GOBSData< GOBSType>  > m_obsdata;
        
    };
    
    
    
    
    
} // end of namespace gfc



#endif /* GStation_hpp */

//
//  GFMThirdBody.hpp
//  GFC
//
//  Created by lizhen on 14/06/2016.
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


#ifndef GFMThirdBody_hpp
#define GFMThirdBody_hpp

#include <stdio.h>
#include "GForceModel.hpp"
#include "GMotionState.hpp"
#include "GJPLEPH.h"
#include <vector>

namespace gfc
{
    
    
    //reference: /* Oliver Montenbruck, P69 and P248
    //
    class GFMThirdBody : public GForceModel
    {
        
    public:
        
        GFMThirdBody()
        {
            setForceName("GFMNbody"); // should be N body
            
            m_bodies.push_back(GJPLEPH::SUN);
            
            m_bodies.push_back(GJPLEPH::MOON);
            
            m_bodies.push_back(GJPLEPH::JUPITER);
            
            m_bodies.push_back(GJPLEPH::VENUS);
            
        }
        
        void setBodies(std::vector<int>& bodies);
        
        void doCompute(GVector& satpos_eci);
        
        //void bindJPLEPH(GJPLEPH* eph);
        
        //reference to Montenbruck, page 69
        GVector getAcc(double gm, GVector& posM, GVector& satpos);
        GMatrix getPartialDerivatives(double gm, GVector& posM, GVector& satpos);
        
    private:
        
        /*
         UNKNOWNBODY = 0, MERCURY = 1, VENUS = 2,
         EARTH = 3,  //add by myself, because 3 only give out the earth-moon barycenter, NOT the earth barycenter
         MARS = 4, JUPITER = 5,
         SATURN  = 6, URANUS= 7, NEPTUNE = 8, PLUTO =9, MOON = 10, SUN=11,
         SS_BARY = 12,
         EMBARY = 13,
         NUTATIONS = 14,
         LIBRATIONS = 15
         */
        static double solar_system_GM[12];
        
        std::vector<int> m_bodies;  // the bodies to be calculated, the index must be consistent with that in GJPLEPH
        
    };
    
    
} // end of namespace gfc



#endif /* GFMThirdBody_hpp */

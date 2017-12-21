//
//  GSpaceEnv.hpp
//  GFC
//
//  Created by lizhen on 15/06/2016.
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

#ifndef GSpaceEnv_hpp
#define GSpaceEnv_hpp

#include <stdio.h>

#include "GJPLEPH.h"

#include "GEarthOrientationParameter.hpp"

#include "GRadiationFlux.hpp"

namespace gfc
{
    // the class for processing the space envirnment, including the GJPLEPH , frame transformation, F10.7, Kp  et.al
    class GSpaceEnv
    {
        
    public:
        
        // here epoch time should be in UTC
        static void updateSpaceEnvironment( GTime epochUTC );
        
    public:
       
       static double F107;  // F10.7
        
       static double Kp;  // Kp index
       
       static double tsi;  //the total solar index at 1 AU, 1362w/m^2
        
       static GTime epoch_utc;  // the time in UTC
       static GTime epoch_tdb; // the time in tdb;
        
       static GJPLEPH eph;
        
       static GEarthOrientationParameter eop;
       
        
        
       static GVector planetPos_eci[12];
       //static GMatrix planetPos_eci[12];  // these are position of solar system planets referenced to earth center
       
       static GVector planetVel_eci[12];
       //static GMatrix planetVel_eci[12]; // these are velocity of solar system plantes referenced to earth center
        
       static GVector planetPos_ecef[12];
       //static GMatrix planetPos_ecef[12];
       
       static GVector planetVel_ecef[12];
       //static GMatrix planetVel_ecef[12];
        
        //these are unit vectors
        
       static GVector sunPosHat_eci;
       //static GMatrix sunPosHat_eci;
       
       static GVector sunVelHat_eci;
       //static GMatrix sunVelHat_eci;
       
       static GVector sunPosHat_ecef;
       //static GMatrix sunPosHat_ecef;
       
       static GVector sunVelHat_ecef;
       //static GMatrix sunVelHat_ecef;
       
       static double disSun_eci;
       static double disSun_sr_eci; //square root of that dis sqrt()
       static double disSun_sq_eci; // square of that dis, that is dis*dis
        
       static std::vector<int> planetsUsed;  // the planets to be calculated
       
        
        
    };

} // end of namespace


#endif /* GSpaceEnv_hpp */

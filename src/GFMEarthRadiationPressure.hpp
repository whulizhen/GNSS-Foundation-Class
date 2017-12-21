//
//  GFMEarthRadiationPressure.hpp
//  GFC
//
//  Created by lizhen on 16/4/12.
//  Copyright © 2016年 lizhen. All rights reserved.
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


#ifndef GFMEarthRadiationPressure_hpp
#define GFMEarthRadiationPressure_hpp

#include <stdio.h>

#include <math.h>

#include <fstream>

#include "GForceModel.hpp"

#include "GIcosahedron.h"

#include "GEarthOrientationParameter.hpp"

#include "GSpaceCraftAttitude.hpp"

#include "GFMSolarRadiationPressure.hpp"

#include "GRadiationFlux.hpp"

#include "GSpacecraft.hpp"

namespace gfc
{
    
   
    
    class GFMEarthRadiationPressure : public GForceModel
    {
        
    public:
        
        GFMEarthRadiationPressure()
        {
            setForceName("GFMERP");
            with_grid_on = false;
        }
        
        ~GFMEarthRadiationPressure()
        {
            
        }
        
    void doCompute( GMotionState* statePointer ,  GSpaceCraftModel* space_craft_geometry );
    
        
    GVector Solano_ERP( GVector& satpos_eci, GVector& satvel_eci );
        
    // simple box shape
    GVector erp_bus_simple( GEarthRadiationFlux& earthRadiationflux ,GSpaceCraftAttitude* svAttitude, GSpaceCraftModel* space_craft_geometry );
    
    // use the grid files from SRP, complex shape
    GVector erp_bus_grid(GEarthRadiationFlux& earthRadiationflux ,GSpaceCraftAttitude* svAttitude, GSpaceCraftModel* space_craft_geometry);
        
    GVector erp_solar_panel( GEarthRadiationFlux* earthRadiationflux, GSpaceCraftAttitude* svAttitude,GSpaceCraftModel* space_craft_geometry);
    
        
        bool with_grid_on;
        
    private:
        
        
        //int m_level;  // the discrete global grid system(DGGS) level
        
        //ofstream pf;
        
    };
    
    
    
}




#endif /* GFMEarthRadiationPressure_hpp */

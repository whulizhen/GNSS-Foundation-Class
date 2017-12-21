//
//  GFMSolarRadiationPressure.hpp
//  GFC
//
//  Created by lizhen on 18/05/2016.
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


#ifndef GFMSolarRadiationPressure_hpp
#define GFMSolarRadiationPressure_hpp

#include <stdio.h>
#include <fstream>
#include "GForceModel.hpp"
#include "GMatrix.h"
#include "GRadiationFlux.hpp"
#include "GSpacecraft.hpp"
#include "GRadiationGrid.hpp"
#include "GFMThermalRadiationForce.hpp"

namespace gfc {

/*
 this class is for the grid srp model
 */

class GFMSolarRadiationPressure : public GForceModel
{
    
 public:
  
    bool with_grid_on;  // determine whether to use grid file or not
    
    GFMSolarRadiationPressure()
    {
        setForceName("GFMSRP");
        with_grid_on = false;
    }
    
    ~GFMSolarRadiationPressure()
    {
        
    }
    
    
  GVector newBoxWing(GSolarRadiationFlux& solarRadiationflux, GSpaceCraftAttitude* svAttitude,GSpaceCraftModel* space_craft_geometry);
    
  GVector  srp_solar_panel( GSolarRadiationFlux& solarRadiationflux, GSpaceCraftAttitude* svAttitude,GSpaceCraftModel* space_craft_geometry );
  
    //box-wing model
  GVector  srp_bus_simple(  GSolarRadiationFlux& solarRadiationflux ,GSpaceCraftAttitude* svAttitude, GSpaceCraftModel* space_craft_geometry  );
  
  GVector  srp_bus_grid(GSolarRadiationFlux& solarRadiationflux ,GSpaceCraftAttitude* svAttitude, GSpaceCraftModel* space_craft_geometry );
   
    
  // this is the satellite bus model from Monternbruck
  //ref: Montenbruck 2015, Enhanced SRP modelling for Galileo satellites
  GVector srp_bus_M(double e, GSpaceCraftAttitude* svAttitude, GSpaceCraftModel* space_craft_geometry);
    
  //void doCompute(GSolarRadiationFlux& solarRadiationflux ,GSpaceCraftAttitude* svAttitude, GSpaceCraftModel* space_craft_geometry);
    
  void doCompute(GSpaceCraft* spacecraft);
    
};

}  // end of namespace

#endif /* GFMSolarRadiationPressure_hpp */

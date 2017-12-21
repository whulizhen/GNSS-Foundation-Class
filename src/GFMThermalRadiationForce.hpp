//
//  GFMThermalRadiationForce.hpp
//  GFC
//
//  Created by lizhen on 10/07/2016.
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

#ifndef GFMThermalRadiationForce_hpp
#define GFMThermalRadiationForce_hpp

#include <stdio.h>
#include "GForceModel.hpp"

#include "GSpacecraft.hpp"

namespace  gfc
{
    //this class is used to model the Thermal Reradiation Pressure for the solar panel
    class GFMThermalRadiationForce : public GForceModel
    {
        
    public:
        
        GFMThermalRadiationForce()
        {
            setForceName("GFMTRR");
        }
        
        ~GFMThermalRadiationForce() {}
        
        
        static GVector ThermalRadiationForce_MLI(gfc::GVector &n, GVector& flux_dir, double flux,double area, double alfa, double emissivity);
        
        void doCompute( GSpaceCraft* spacecraft );
        
        //get the solar panel surface temperature
        void getTemperatures( double radiation_front, double radiation_back, double power_draw,GSpaceCraftModel* space_craft_geometry);
        
    private:
        
        // this is a privat
        //void process_flux();
        
        double m_temperatures[3]; // the 3 temperatures of three surfaces of the panel
        // Stefan-Boltzmann constant
        static double sigma ; // Wm^-2K^-4;
        
    };
    
    
} // end of namespace gfc



#endif /* GFMThermalRadiationForce_hpp */

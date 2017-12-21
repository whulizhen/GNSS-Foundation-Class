//
//  GConfigure.hpp
//  GFC
//
//  Created by lizhen on 16/04/2017.
//  Copyright © 2017 lizhen. All rights reserved.
//

#ifndef GConfigure_hpp
#define GConfigure_hpp
#include "GString.h"
#include <fstream>
#include <stdio.h>
#include "GJPLEPH.h"
#include "GRadiationFlux.hpp"
#include "GForceModelMgr.hpp"

namespace gfc
{
    class GConfigure
    {
        struct cfgdata
        {
            GString spacecraftmodel; // the spacecraft models description
            GString planetEphemeris; // planet ephemeris file
            GString eopfile; // earth orientation parameter file
            GString earthGravityFile; // earth gravity coefficients file (normalised)
            GString earthRadiationDirectory; // triangular earth radiation flux
            
            int     earthGravityDegree;
            
            double stepsize;
            std::vector<GString> forcelist;
            std::vector<int>     planetlist;  // the planets list
            std::vector<bool>    forcePartial;  // withou or without partials
            bool time_variable_gravity;
            bool solid_earth_tide;
            bool ocean_earth_tide;
            bool polar_tide;
            bool y_bias;  // true or false

            int srp_option; // 0 for no srp , 1 for box-wing, 2 for grid file
            int erp_option; // 0 for no srp , 1 for box-wing, 2 for grid file
            int earth_radiation_option;  // 0 for no model, 1 for simple, 2 for ceres_trangle, 3 for ceres_original
            int max_earth_division;  // the maximum division of the TOA for CERES earth radiation flux
            
            cfgdata()
            {
                earthGravityDegree = 0;
                stepsize = 20;
                time_variable_gravity = false;
                solid_earth_tide = false;
                ocean_earth_tide = false;
                polar_tide = false;
                y_bias = false;
            }
            
        };
        
    public:
        
        void parseCfg(GString cfgfile);
        
        void configForceModelMgr(GForceModelMgr& modelMgr);
        
        cfgdata config;
        
    };
    
    
}



#endif /* GConfigure_hpp */

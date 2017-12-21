//
//  GForceModelMgr.h
//  GFC
//
//  Created by lizhen on 22/07/2016.
//  Copyright © 2016 lizhen. All rights reserved.
//

#ifndef GForceModelMgr_h
#define GForceModelMgr_h

#include <map>

#include "GTime.h"
#include "GForceModel.hpp"
#include "GFMGravity.hpp"
#include "GFMEarthGravity.hpp"
#include "GFMEarthRadiationPressure.hpp"
#include "GFMSolarRadiationPressure.hpp"
#include "GFMSolarRadiaitonPressureEM.hpp"
#include "GFMThermalRadiationForce.hpp"
#include "GFMThirdBody.hpp"
#include "GFMEarthTide.hpp"
#include "GSpacecraft.hpp"

namespace gfc
{
    // this is a class to manage all the force models
    class GForceModelMgr
    {
        
    public:
        
        GForceModelMgr()
        {
            totalNum_parameters = 0;
            gravity_order = 20;
            gravity_degree = 20;
        }
        
        // set the force model name list for the current spacecraft
        void setForceList( std::vector<GString>& forceNames ,std::vector<bool>& withPartials );
        void getDerivatives( GTime ct, int n, double x, double *y, double *dydx, GSpaceCraft* spacecraft );
        
        int gravity_order;
        int gravity_degree;
        GString gravity_file;
        
        // a string list to store all the force names acting on this spacecraft
        std::map< GString, GForceModel* > m_forceModels;
        int totalNum_parameters;  // the total number of parameters for all the force models
        
    };
    
} // end of namespace

#endif /* GForceModelMgr_h */

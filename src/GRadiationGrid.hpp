
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

//
//  GRadiationForce.hpp
//  GFC
//
//  Created by lizhen on 13/06/2016.
//  Copyright © 2016 lizhen. All rights reserved.
//

#ifndef GRadiationForce_hpp
#define GRadiationForce_hpp

#include <stdio.h>

#include "GString.h"

#include <map>

#include "GMatrix.h"
#include "GVector.hpp"

#include <fstream>


namespace gfc
{
    
    /*
     * a class to manage the solar radiation grid files
     *
     */
    class GRadiationGridMgr
    {
        
    public:
        
        static const int grid_rows = 181;
        static const int grid_cols = 361;
        static const int min_latitude = -90.0;
        static const int max_latitude = 90.0;
        static const int min_longitude = -180.0;
        static const int max_longitude = 180.0;
        
        //static double nominal_mass; //kg
        
         static double nominal_flux;  // 1368 w/m^2
        
       // static void loadGridFile( GString filename[3] , GString svType );
        
         static bool readGridFile( GString file_name,double nominal_mass, GMatrix& data );
        
         static GVector getAcc_grid( GMatrix (&griddata)[3], double longitude, double latitude);
         
    private:
        
        static GVector bilinear_interp(GMatrix (&griddata)[3], double longitude, double latitude);
        static GVector mybilinear_interp(GMatrix (&griddata)[3], double longitude, double latitude);
        
        // svType, and gridata
       // static std::map< GString, GMatrix* > griddata;
        
    };
    
    
    
} // end of namespace gfc



#endif /* GRadiationForce_hpp */

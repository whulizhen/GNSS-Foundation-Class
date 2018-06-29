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

//
//  GFMSolarRadiaitonPressureEM.hpp
//  GFC
//
//  Created by lizhen on 26/10/2016.
//  Copyright © 2016 lizhen. All rights reserved.
//

#ifndef GFMSolarRadiaitonPressureEM_hpp
#define GFMSolarRadiaitonPressureEM_hpp

#include <fstream>
#include "GForceModel.hpp"
#include <stdio.h>
#include "GMatrix.h"
#include "GRadiationFlux.hpp"
#include "GSpacecraft.hpp"
#include "GRadiationGrid.hpp"

namespace gfc
{
    class GFMSolarRadiationPressureEM : public GForceModel
    {
        
    public:
        
        int type_opt; // the name of the type, 0 for ECOM1, 1 for ECOM2, 2 for DREMT
        
        GFMSolarRadiationPressureEM();
        ~GFMSolarRadiationPressureEM();
        
        // the reduced ECOM 5 parameters model
        //ref CODE's new solar radiation pressure model for GNSS orbit determination
        GVector ECOM1( double* param,double u, double factor );
        void ECOM1_dadp(double u, double factor, GMatrix& daDYBdp );
        
        GVector ECOM2( double* param,double u, double factor );
        void ECOM2_dadp(double u, double factor, GMatrix& daDYBdp );
        
        // the reduces ECOM2
        GVector ECOM3( double* param,double u, double factor );
        void ECOM3_dadp(double u, double factor, GMatrix& daDYBdp );
        
        
        
        /*
         model description in Body Fixed System
         x =  C0 + C1*cos(phi) + C2*cos(2phi) + C3*sin(2phi)
         z = -C0 - C1*cos(phi) - C2*cos(2phi) - C3*sin(2phi)
         y =  y0
         new empirical model test
         phi: the latitude of sun in bfs
         factor: the scaling factor
         param: the parameters in this model
        */
        GVector DREM(double* param,double phi,double factor);
        void DREM_dadp(double phi,double factor,GMatrix& daXYZdp);
        
        GVector DREMT(double* param,double phi,double factor);
        void DREMT_dadp(double phi,double factor,GMatrix& daXYZdp);
        
        GVector DREM1T(double* param,double phi,double factor);
        void DREM1T_dadp(double phi,double factor,GMatrix& daXYZdp);
        
        GVector DREM1T_eps(double* param,double eps,double eta,double factor);
        void DREM1T_eps_dadp(double eps,double eta, double factor,GMatrix& daDYBdp);
        
        
        GVector DREM2T(double* param,double phi,double factor);
        void DREM2T_dadp(double phi,double factor,GMatrix& daXYZdp);
        
        // combine with DREM2T with ECOM2?
        GVector DREM3T(double* param,double phi,double beta, double eta,double factor);
        void DREM3T_dadp(double phi,double beta, double eta,double factor,GMatrix& daXYZdp);
        
        GVector Montenbruck(double* param, double eps, double factor);
        void Montenbruck_dadp(double eps,double factor,GMatrix& daXYZdp);
        
        
        GVector DREM1(double* param,double phi,double factor);
        void DREM1_dadp(double phi,double factor,GMatrix& daXYZdp);
        
        GVector DREM2(double* param,double phi,double factor);
        void DREM2_dadp(double phi,double factor,GMatrix& daXYZdp);
        
        
//        GVector empDYB(double* param,double phi, double factor);
//        void empDYB_dadp(double phi,GMatrix& daXYZdp);
//        
//        
//        void empXYZ_dadp_DYB(double phi,GMatrix& daXYZdp);
//        GVector empXYZ_DYB(double* param,double phi, double factor);
        
        void doCompute( GSpaceCraft* spacecraft );
        
   
        
        
    }; // end of class GFMSolarRadiaitonPressureEM
    
    
    
} // end of namespace





#endif /* GFMSolarRadiaitonPressureEM_hpp */

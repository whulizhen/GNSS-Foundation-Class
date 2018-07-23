
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
//  GMotionState.hpp
//  GFC
//
//  Created by lizhen on 05/06/2016.
//  Copyright © 2016 lizhen. All rights reserved.
//

#ifndef GMotionState_hpp
#define GMotionState_hpp

#include <stdio.h>
#include <vector>
#include "GTime.h"
#include "GSpaceCraftAttitude.hpp"
#include "GEarthOrientationParameter.hpp"
#include "GJPLEPH.h"

#include "GRadiationFlux.hpp"

#include "GSpaceEnv.hpp"

#include "GVector.hpp"

#include "GKPE.hpp"

namespace gfc
{
    
        
    
    
    /*
     * purpose: used to describle the geometry of satellite, planets ,and the attitude
     */
    class GMotionState
    {
        
    public:
        
        static double shadowFactor(double a, double b, GVector& sunpos_eci, GVector& satpos_eci);
        
        static int perspectiveProjection(double a, double b,GVector& sunpos_eci, GVector& satpos_eci, double& r_solar,double& S, double* EC_intersection,double* EC );
        
        static double area_ellispe(double Q1[2], double Q2[2], double Os[2], double Rs,double a, double b, bool in_out);
        static double area_hyperbola(double Q1[2], double Q2[2], double Os[2], double Rs,double a, double b, bool in_out, bool x_axis);
        
        static int myperspectiveProjection(double a, double b, GVector& sunpos_ecef, GVector& satpos_ecef, double& r_solar, double& area_bright, double& dis_boundary, double& dis_circle);
        
        static double myshadowFactor(GVector& sunpos_eci, GVector& satpos_eci);
        
        
        // if on is false, it is the montenbruck's apporach
        // if on is true, it is the modified montenbruck's apporach, very close the real solution
        static double shadowFactor_SECM(bool on, GVector& XSun, GVector& XSat);
        
        //detect all the shadow event, including begining time and end time
        void shadow_detector(GTime gpst,gfc::GVector &XSun, gfc::GVector &XSat);
        
        
        GMotionState();
        
        //time system in epoch should be in UTC
        void updateEpoch( GTime epoch );
        
        void updateState_eci(GTime epoch_utc, GVector psatpos_eci, GVector psatvel_eci );
        
        void updateState_ecef(GTime epoch_utc, GVector psatpos_ecef, GVector psatvel_ecef );
        
        //update the current state transition matrix drt/dr0
        void updateTransitionMatrix(double* phi);
        void updateTransitionMatrix(GMatrix phi);
        
        // update the current sensitivity matrix
        void updateSensitivityMatrix(double* S);
        void updateSensitivityMatrix(GMatrix S);
        
        GSpaceCraftAttitude* getAttitude();
        GSpaceCraftAttitude* getAttitude_ecef();
        
        /*dump the state to file*/
        void outputState( ofstream& file );
        
        GTime  m_epoch;  // must in UTC
        CivilTime ct;  // the year,month day in UTC
        
        // satllite positions are all in !!!! CENTER OF MASS !!!!
        
        GVector satpos_ecef;
        GVector satposHat_ecef; // the unit vector of satpos_ecef
        GVector satvel_ecef;
        GVector satvelHat_ecef; // the unit vector of satvel_ecef
        
        GVector satpos_eci;
        GVector satposHat_eci; // the unit vector of satvel_ecef
        GVector satvel_eci;
        GVector satvelHat_eci; // the unit vector of satvel_eci
        
        GKeplerianElements keplerianElement;  // the keplerian elements at this time
        
        GVector sat_sun_ecef; // the vector from sat to sun in ecef
        GVector sat_sunHat_ecef; // the unit vector from sat to sun;
        
        GVector sat_sun_eci; // the vector from sat to sun in eci
        GVector sat_sunHat_eci; // the unit vector from sat to sun in eci;
        
        double dis_sat;    // the distance from the sat to earth center
        double dis_sat_sun;  // in km
        
        GVector orbit_node_vector_eci; // the orbital node vector,ascending node
        double eps;  // the earth-probe-sun angle in radian
        //double u;    // the argument of latitude
        
        double shadow_factor;  // the shadow factor for the current state
        double dis_factor; // the distance factor for the solar flux
        
        GEarthRadiationFlux earthFlux; // the earth radiation flux for the current satellite
        
        GSolarRadiationFlux solarFlux; // the solar radiation flux for the current satellite
        
        GSpaceCraftAttitude attitude_eci; // the eci attitude
        
        GSpaceCraftAttitude attitude_ecef; // the ecef attitude
        
        GMatrix phiMatrix;  // the state transition matrix, 6 by 6
        
        // sensitivity matrix is the partial derivatives of state w.r.t force model parameters
        GMatrix senMatrix;   // the sensitivity matrix 6 by np, np is the number of other parameters
        
    };
    
    
} // end of namespace



#endif /* GMotionState_hpp */

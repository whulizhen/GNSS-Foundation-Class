//
//  GSpaceEnv.cpp
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


#include "GSpaceEnv.hpp"
namespace gfc
{
    //all the declearation for the static variables
    
    double GSpaceEnv::F107;  // F10.7
    
    double GSpaceEnv::Kp;  // Kp index
    
    double GSpaceEnv::tsi = 1362;  //the total solar index at 1 AU, 1362w/m^2
    
    GTime GSpaceEnv::epoch_utc;  // the time in UTC
    GTime GSpaceEnv::epoch_tdb; // the time in tdb;
    
    GJPLEPH GSpaceEnv::eph;
    
    GEarthOrientationParameter GSpaceEnv::eop;
    
    GVector GSpaceEnv::planetPos_eci[12];  // these are position of solar system planets referenced to earth center
    
    GVector GSpaceEnv::planetVel_eci[12]; // these are velocity of solar system plantes referenced to earth center
    
    GVector GSpaceEnv::planetPos_ecef[12];
    
    GVector GSpaceEnv::planetVel_ecef[12];
    
    //these are unit vectors
    GVector GSpaceEnv::sunPosHat_eci;
    
    GVector GSpaceEnv::sunVelHat_eci;
    
    GVector GSpaceEnv::sunPosHat_ecef;
    
    GVector GSpaceEnv::sunVelHat_ecef;
    
    double GSpaceEnv::disSun_eci;
    double GSpaceEnv::disSun_sr_eci; //square root of that dis sqrt()
    double GSpaceEnv::disSun_sq_eci; // square of that dis, that is dis*dis
    
    std::vector<int> GSpaceEnv::planetsUsed = {GJPLEPH::SUN, GJPLEPH::VENUS, GJPLEPH::MOON, GJPLEPH::JUPITER};  // the planets to be calculated
    
    
    void GSpaceEnv::updateSpaceEnvironment( GTime epochUTC )
    {
        
        if(epoch_utc != epochUTC )
        {
            epoch_utc = epochUTC;
            
            eop.setEpochTime(epoch_utc);
            
            
            F107 = 100;
            
            Kp = 100;
            
            tsi = 1362.0;
            
            JDTime jdt = GTime::GTime2JDTime( epoch_utc );
            GTime tai = GTime::UTC2TAI(epoch_utc);
            GTime tt = GTime::TAI2TT(tai);
            
            
            eph.SetUnit(true);
            eph.setEpoch(epoch_utc); ///here should be in TDB
            
            // the maximum error is 2 us for this transformation of TT to TDB
            double tdbmtt = GTime::TDBmTT( tt ,eop.getUT1mUTC() , NULL  );
            
            epoch_tdb  = GTime::TT2TDB(tt, tdbmtt);
            eph.setEpoch(epoch_tdb);
            
            double pv_eci[6] = {0.0};
            double pv_ecef[6] = {0.0};
            double p_eci[3] = {0.0};
            double v_eci[3] = {0.0};
            double p_ecef[3] = {0.0};
            double v_ecef[3] = {0.0};
            
            for( int i = 0 ; i< planetsUsed.size(); i++ )
            {
                eph.getPos(GJPLEPH::EARTH, planetsUsed[i], pv_eci);
                //here should use the apparent positions
                //eph.getPos_apparent( GJPLEPH::EARTH , planetsUsed[i], pv_eci );
                
                eop.ECI2ECEF(1, pv_eci, pv_ecef);
                
                for( int i = 0 ; i< 6; i++ )
                {
                    if( i< 3 )
                    {
                        p_eci[i] = pv_eci[i];
                        p_ecef[i] = pv_ecef[i];
                    }
                    else
                    {
                        v_eci[i-3] = pv_eci[i];
                        v_ecef[i-3]= pv_ecef[i];
                    }
                }
                
                planetPos_eci[planetsUsed[i]].set(p_eci[0], p_eci[1], p_eci[2]);
                
                planetVel_eci[planetsUsed[i]].set(v_eci[0], v_eci[1], v_eci[2]);
                
                planetPos_ecef[planetsUsed[i]].set(p_ecef[0], p_ecef[1], p_ecef[2]);
                planetVel_ecef[planetsUsed[i]].set(v_ecef[0], v_ecef[1], v_ecef[2]);
                
                if( planetsUsed[i] == GJPLEPH::SUN ) //special for sun
                {
                    
                    disSun_eci = planetPos_eci[GJPLEPH::SUN].norm(); // in km
                    
                    disSun_sr_eci = sqrt(disSun_eci); //squre root
                    
                    disSun_sq_eci = disSun_eci*disSun_eci; // sqare
                    
                    sunPosHat_eci = planetPos_eci[GJPLEPH::SUN];
                    sunPosHat_eci /= disSun_eci;
                    
                    long double disSun_ecef = planetPos_ecef[GJPLEPH::SUN].norm(); // in km
                    sunPosHat_ecef = planetPos_ecef[GJPLEPH::SUN];
                    sunPosHat_ecef /= disSun_ecef;
                    
                    double dis_vel_eci = planetVel_eci[GJPLEPH::SUN].norm();
                    double dis_vel_ecef = planetVel_ecef[GJPLEPH::SUN].norm();
                    
                    sunVelHat_eci = planetVel_eci[GJPLEPH::SUN];
                    sunVelHat_eci /= dis_vel_eci;
                    
                    sunVelHat_ecef = planetVel_ecef[GJPLEPH::SUN];
                    sunVelHat_ecef /= dis_vel_ecef;
                    
                }
                
            }
        }
        
    }
    
    
    
    
} // end of namespace

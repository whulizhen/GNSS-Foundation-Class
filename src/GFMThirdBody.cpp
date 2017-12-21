//
//  GFMThirdBody.cpp
//  GFC
//
//  Created by lizhen on 14/06/2016.
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


#include "GFMThirdBody.hpp"
namespace gfc
{
    
    //  maybe we need to get these values from GJPLEPH.constant
    /*
     UNKNOWNBODY = 0, MERCURY = 1, VENUS = 2,
     EARTH = 3,  //add by myself, because 3 only give out the earth-moon barycenter, NOT the earth barycenter
     MARS = 4, JUPITER = 5,
     SATURN  = 6, URANUS= 7, NEPTUNE = 8, PLUTO =9, MOON = 10, SUN=11,
     SS_BARY = 12,
     EMBARY = 13,
     NUTATIONS = 14,
     LIBRATIONS = 15
     */
    double GFMThirdBody::solar_system_GM[12]=
    {
        0.0,  // unknown
        22031.780000,  // mercury
        324858.592000, // venus
        398600.4415,           //earth, this is just the third body, don not consider earth
        42828.375214,  // mars
        126712764.800000, // jupiter
        37940585.200000,  // saturn
        5794548.600000,  // uranus
        6836527.100580, // neptune
        977.000000,     //pluto
        4902.800066,    // moon
        132712440041.939400  // sun
    };
    
//    void GFMThirdBody::bindJPLEPH(GJPLEPH* eph)
//    {
//        m_eph = eph;
//    }
    
    
    void GFMThirdBody::setBodies(std::vector<int>& bodies)
    {
        m_bodies = bodies;
    }
    
    /*
     * the input epoch should be in UTC
     * according to Montenbruck, the coorindates must be geocentric
     * because the distance between earth and sun, moon is very large, just ignore the light time.
     */
    void GFMThirdBody::doCompute( GVector& satpos_eci)
    {
        m_dadr.resize(3, 3);
        
        GVector acc ;  // M/s^2
        for( int i = 0 ; i< m_bodies.size() ; i++ )
        {
           acc += getAcc(solar_system_GM[m_bodies[i]], GSpaceEnv::planetPos_eci[m_bodies[i]], satpos_eci );
           
            if( m_hasPartialDerivatives == true)
            {
               m_dadr += getPartialDerivatives(solar_system_GM[m_bodies[i]],GSpaceEnv::planetPos_eci[m_bodies[i]], satpos_eci);
            }
        }
        
        setForce(acc);
        
    } // end of function doCompute
    
    
    GMatrix GFMThirdBody::getPartialDerivatives(double gm, GVector &posM, GVector &satpos)
    {
        GMatrix dadr(3,3);
        double t[9]={1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0};
        GMatrix I(t,3,3);
        
        GVector rs = satpos - posM;
        
        long double r_sc =  rs.norm2();
        
        double GMrms3 =  gm/(r_sc*sqrt(r_sc));
        
        double GMrms5 = gm/(r_sc*r_sc*sqrt(r_sc));
        
        GMatrix b(3,3);
        b(0,0) = rs.x*rs.x;b(0,1) = rs.x*rs.y;b(0,2) = rs.x*rs.z;
        b(1,0) = rs.x*rs.y;b(1,1) = rs.y*rs.y;b(1,2) = rs.y*rs.z;
        b(2,0) = rs.x*rs.z;b(2,1) = rs.z*rs.y;b(2,2) = rs.z*rs.z;
        
        dadr = 3.0*b*GMrms5 - GMrms3*I;
        
        //dadr is in km
        
        return dadr;
    }
    
    /*
     * ref: Montenbruck's book Page: 69
     * the unit of acc is in M/s^2 in ECI
     */
    GVector GFMThirdBody::getAcc(double gm, GVector& posM, GVector& satpos)
    {
        GVector acc;
        
        GVector smr = posM - satpos;
        
        long double r_sc =  smr.norm2();
        
        long double r_geo  = posM.norm2();
        
        double GMr_sc = gm / (r_sc * std::sqrt(r_sc));
        
        double GMr_geo = gm / (r_geo * std::sqrt(r_geo));
        //turn km/s^2 to m/s^2
        acc = ( GMr_sc*smr - GMr_geo*posM ) *1000.0;
        
        return acc;
        
    } // end of function getAcc
    
    
    
    
} // end of namespce

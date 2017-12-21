//
//  GFMEarthTide.cpp
//  GFC
//
//  Created by lizhen on 20/07/2016.
//  Copyright © 2016 lizhen. All rights reserved.
//

#include "GFMEarthTide.hpp"

namespace gfc
{
    void GFMEarthTide::doCompute(gfc::GVector &satpos_eci, gfc::GVector &sunpos_eci, gfc::GVector &moonpos_eci)
    {
        GVector acc;
        
        double R_Earth  =  6378.1363;  // km
        double GM_Sun   =  132712440041.939400; //km^3/s^2
        double GM_Moon = 4902.800066 ;
        
        double RSun, RMoon, RSat;
        double cosSun, cosMoon;
        double T5;
        
        //RSun = VectDot( 3, 3, SunPos, SunPos );
        //RSun = sqrt( RSun );
        RSun = sunpos_eci.norm();
        
        //RMoon = VectDot( 3, 3, MoonPos, MoonPos );
        //RMoon = sqrt( RMoon );
        RMoon = moonpos_eci.norm();
        
        //RSat = VectDot( 3, 3, Pos, Pos );
        //RSat = sqrt( RSat );
        RSat = satpos_eci.norm();
        
        //cosSun = VectDot( 3, 3, Pos, SunPos ); /*太阳、卫星的地心夹角*/
        cosSun = dotproduct(satpos_eci, sunpos_eci); /*太阳、卫星的地心夹角*/
        cosSun = cosSun / RSat / RSun;
        cosMoon = dotproduct(satpos_eci, moonpos_eci);
        cosMoon = cosMoon / RSat / RMoon;
        T5      = pow(R_Earth,5)/pow(RSat,4);
        
        
        acc = 0.15 * ( GM_Sun/pow(RSun,3) * T5 * ((3.0-15.0*cosSun*cosSun)
                                        *satpos_eci/RSat + 6.0*cosSun * sunpos_eci/RSun ))
            + 0.15 * ( GM_Moon/pow(RMoon,3) * T5 * ((3.0-15.0*cosMoon*cosMoon)
                                        *satpos_eci/RSat + 6.0*cosMoon * moonpos_eci/RMoon ));
        
        acc = acc*1000.0;
        
//        for( i=0; i<3; i++ )
//        {
//            Acc[i] = 0.15 * ( GM_Sun/pow(RSun,3) * T5 * ((3.0-15.0*cosSun*cosSun)
//                                                         *Pos[i]/RSat + 6.0*cosSun * SunPos[i]/RSun ))
//            + 0.15 * ( GM_Moon/pow(RMoon,3) * T5 * ((3.0-15.0*cosMoon*cosMoon)
//                                                    *Pos[i]/RSat + 6.0*cosMoon * MoonPos[i]/RMoon ));
//        }
        
        setForce(acc);
        
    }
    
    
}
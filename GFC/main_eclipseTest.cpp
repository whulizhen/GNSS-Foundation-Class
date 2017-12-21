//
//  main_eclipseTest.cpp
//  GFC
//
//  Created by lizhen on 10/08/2017.
//  Copyright © 2017 lizhen. All rights reserved.
//

#include <stdio.h>
#include <stdio.h>
#include <algorithm>
#include <typeinfo>


//#include "../../src/GNetcdf.h"
#include "../../src/GEarthOrientationParameter.hpp"
#include "../../src/GRungeKutta.hpp"
#include "GEarthRadiationModel.h"
#include "GFMEarthRadiationPressure.hpp"
#include "GFMThermalRadiationForce.hpp"

/*
#include <GeographicLib/Constants.hpp>
#include <GeographicLib/Geodesic.hpp>
#include <GeographicLib/Gnomonic.hpp>
#include <GeographicLib/PolygonArea.hpp>
*/


#include "GOrbitPredictor.hpp"
#include "GSpaceEnv.hpp"
#include "GOrbitFitting.hpp"
#include "GGlobalPressureTemperature.hpp"
#include "GOBSData.h"

#include "GRinex.hpp"

#include "GSLRcrd.hpp"

#include "GTemperatureTest.hpp"

#include "GSLRValidation.hpp"

#include "GMath.hpp"

#include "GAdams.hpp"

#include "GFMEarthGravity.hpp"

#include "GEllipsoidMgr.h"

#include "GConfigure.hpp"


#include "GFMThermalRadiationForce.hpp"

using namespace gfc;

int main(int argc, char* argv[])
{
    
    
    //sp3/ilrsb.orb.lageos2.160416.v35.sp3   sp3/gfz.orb.lageos2.160416.v35.sp3
    
    //    GSpaceCraftMgr::loadPreciseEphemeris("sp3/wum17860.sp3"); // wum18000.sp3,96,2014,4,6
    //    GSpaceCraftMgr::loadPreciseEphemeris("sp3/wum17861.sp3"); // wum18000.sp3,97
    //    GSpaceCraftMgr::loadPreciseEphemeris("sp3/wum17862.sp3"); // wum18000.sp3,98
    //    GSpaceCraftMgr::loadPreciseEphemeris("sp3/wum17863.sp3"); // wum18000.sp3,99
    //    GSpaceCraftMgr::loadPreciseEphemeris("sp3/wum17864.sp3"); // wum18000.sp3,100
    //    GSpaceCraftMgr::loadPreciseEphemeris("sp3/wum17865.sp3"); // wum18000.sp3,101
    //    GSpaceCraftMgr::loadPreciseEphemeris("sp3/wum17866.sp3"); // wum18000.sp3,101
    //
    //    //
    //  GSpaceCraftMgr::loadPreciseEphemeris("sp3/wum17870.sp3"); // wum18000.sp3,96,2014,4,6
    //  GSpaceCraftMgr::loadPreciseEphemeris("sp3/wum17871.sp3"); // wum18000.sp3,97
    //  GSpaceCraftMgr::loadPreciseEphemeris("sp3/wum17872.sp3"); // wum18000.sp3,98
    //  GSpaceCraftMgr::loadPreciseEphemeris("sp3/wum17873.sp3"); // wum18000.sp3,99
    //  GSpaceCraftMgr::loadPreciseEphemeris("sp3/wum17874.sp3"); // wum18000.sp3,100
    //  GSpaceCraftMgr::loadPreciseEphemeris("sp3/wum17875.sp3"); // wum18000.sp3,101
    //  GSpaceCraftMgr::loadPreciseEphemeris("sp3/wum17876.sp3"); // wum18000.sp3,102
    //
    //  GSpaceCraftMgr::loadPreciseEphemeris("sp3/wum17880.sp3"); // wum18000.sp3,103
    //  GSpaceCraftMgr::loadPreciseEphemeris("sp3/wum17881.sp3"); // wum18000.sp3,104
    //  GSpaceCraftMgr::loadPreciseEphemeris("sp3/wum17882.sp3"); // wum18000.sp3,105
    //  GSpaceCraftMgr::loadPreciseEphemeris("sp3/wum17883.sp3");
    //  GSpaceCraftMgr::loadPreciseEphemeris("sp3/wum17884.sp3");
    //  GSpaceCraftMgr::loadPreciseEphemeris("sp3/wum17885.sp3");
    //  GSpaceCraftMgr::loadPreciseEphemeris("sp3/wum17886.sp3");
    //
    //  GSpaceCraftMgr::loadPreciseEphemeris("sp3/wum17890.sp3");
    //  GSpaceCraftMgr::loadPreciseEphemeris("sp3/wum17891.sp3");
    //  GSpaceCraftMgr::loadPreciseEphemeris("sp3/wum17892.sp3");
    //  GSpaceCraftMgr::loadPreciseEphemeris("sp3/wum17893.sp3");
    //  GSpaceCraftMgr::loadPreciseEphemeris("sp3/wum17894.sp3");
    //  GSpaceCraftMgr::loadPreciseEphemeris("sp3/wum17895.sp3");
    //  GSpaceCraftMgr::loadPreciseEphemeris("sp3/wum17896.sp3");
    
    
    // choose a spacecraft, prn G11, svn46
    //these are for the gps test ()
    //   GSpaceCraftMgr::loadPreciseEphemeris("sp3/cod19142.eph");
    //   GSpaceCraftMgr::loadPreciseEphemeris("sp3/cod19143.eph");
    //   GSpaceCraftMgr::loadPreciseEphemeris("sp3/cod19144.eph");
    //   GSpaceCraftMgr::loadPreciseEphemeris("sp3/cod19145.eph");
    //    GSpaceCraftMgr::loadPreciseEphemeris("sp3/cod19146.eph");
    //    GSpaceCraftMgr::loadPreciseEphemeris("sp3/cod19150.eph");
    //    GSpaceCraftMgr::loadPreciseEphemeris("sp3/cod19151.eph");
    //    GSpaceCraftMgr::loadPreciseEphemeris("sp3/cod19152.eph");
    //    GSpaceCraftMgr::loadPreciseEphemeris("sp3/cod19153.eph");
    //    GSpaceCraftMgr::loadPreciseEphemeris("sp3/cod19154.eph");
    //
    
    //for galileo(ecllipse season)
    // GSpaceCraftMgr::loadPreciseEphemeris("sp3/codsp3-2015/com18264.sp3");
    // GSpaceCraftMgr::loadPreciseEphemeris("sp3/codsp3-2015/com18265.sp3");
    // GSpaceCraftMgr::loadPreciseEphemeris("sp3/codsp3-2015/com18266.sp3");
    // GSpaceCraftMgr::loadPreciseEphemeris("sp3/codsp3-2015/com18270.sp3");
    // GSpaceCraftMgr::loadPreciseEphemeris("sp3/codsp3-2015/com18271.sp3");
    
    
    //for galileo( non ecllipse season)
    //    GSpaceCraftMgr::loadPreciseEphemeris("sp3/codsp3-2015/com18316.sp3");  //214
    //    GSpaceCraftMgr::loadPreciseEphemeris("sp3/codsp3-2015/com18320.sp3");
    //    GSpaceCraftMgr::loadPreciseEphemeris("sp3/codsp3-2015/com18321.sp3");
    //    GSpaceCraftMgr::loadPreciseEphemeris("sp3/codsp3-2015/com18322.sp3");
    //    GSpaceCraftMgr::loadPreciseEphemeris("sp3/codsp3-2015/com18323.sp3");
    //    GSpaceCraftMgr::loadPreciseEphemeris("sp3/codsp3-2015/com18324.sp3");
    //    GSpaceCraftMgr::loadPreciseEphemeris("sp3/codsp3-2015/com18325.sp3");
    //    GSpaceCraftMgr::loadPreciseEphemeris("sp3/codsp3-2015/com18326.sp3");
    //    GSpaceCraftMgr::loadPreciseEphemeris("sp3/codsp3-2015/com18330.sp3");
    //    GSpaceCraftMgr::loadPreciseEphemeris("sp3/codsp3-2015/com18331.sp3");
    //    GSpaceCraftMgr::loadPreciseEphemeris("sp3/codsp3-2015/com18332.sp3");
    //    GSpaceCraftMgr::loadPreciseEphemeris("sp3/codsp3-2015/com18333.sp3");
    
    /*
     GSpaceCraftMgr::loadPreciseEphemeris("sp3/codsp3-2015/com18334.sp3"); // 226
     GSpaceCraftMgr::loadPreciseEphemeris("sp3/codsp3-2015/com18335.sp3");
     GSpaceCraftMgr::loadPreciseEphemeris("sp3/codsp3-2015/com18336.sp3");
     GSpaceCraftMgr::loadPreciseEphemeris("sp3/codsp3-2015/com18340.sp3");
     GSpaceCraftMgr::loadPreciseEphemeris("sp3/codsp3-2015/com18341.sp3");
     GSpaceCraftMgr::loadPreciseEphemeris("sp3/codsp3-2015/com18342.sp3");
     GSpaceCraftMgr::loadPreciseEphemeris("sp3/codsp3-2015/com18343.sp3");
     GSpaceCraftMgr::loadPreciseEphemeris("sp3/codsp3-2015/com18344.sp3");
     GSpaceCraftMgr::loadPreciseEphemeris("sp3/codsp3-2015/com18345.sp3");
     GSpaceCraftMgr::loadPreciseEphemeris("sp3/codsp3-2015/com18346.sp3");
     GSpaceCraftMgr::loadPreciseEphemeris("sp3/codsp3-2015/com18350.sp3");
     GSpaceCraftMgr::loadPreciseEphemeris("sp3/codsp3-2015/com18351.sp3");
     GSpaceCraftMgr::loadPreciseEphemeris("sp3/codsp3-2015/com18352.sp3");
     GSpaceCraftMgr::loadPreciseEphemeris("sp3/codsp3-2015/com18353.sp3");
     GSpaceCraftMgr::loadPreciseEphemeris("sp3/codsp3-2015/com18354.sp3");
     GSpaceCraftMgr::loadPreciseEphemeris("sp3/codsp3-2015/com18355.sp3");
     GSpaceCraftMgr::loadPreciseEphemeris("sp3/codsp3-2015/com18356.sp3");
     GSpaceCraftMgr::loadPreciseEphemeris("sp3/codsp3-2015/com18360.sp3");
     GSpaceCraftMgr::loadPreciseEphemeris("sp3/codsp3-2015/com18361.sp3");
     GSpaceCraftMgr::loadPreciseEphemeris("sp3/codsp3-2015/com18362.sp3");
     GSpaceCraftMgr::loadPreciseEphemeris("sp3/codsp3-2015/com18363.sp3");
     GSpaceCraftMgr::loadPreciseEphemeris("sp3/codsp3-2015/com18364.sp3");
     */
    
    //GSpaceCraftMgr::loadPreciseEphemeris("sp3/com17821.sp3");
    //GSpaceCraftMgr::loadPreciseEphemeris("sp3/com17822.sp3");
    //GSpaceCraftMgr::loadPreciseEphemeris("sp3/com17823.sp3");
    //GSpaceCraftMgr::loadPreciseEphemeris("sp3/com17824.sp3");
    
    
    //GSpaceCraftMgr::loadUCLEphemeris("sp3/ucl_svn46_01mar04.eci");
    
    // testing data
    
    argv[1] = "gfcsetup.cfg";
    argv[2] = "ssBDS";  //satellite system;
    argv[3] = "14";     //prn
    
    //2015/01/11/18:34:25.000000
    //argv[4] = "2015 1 11 17 34 43.00000000"; // start time in GPST
    //argv[5] = "2015 1 11 17 36 04.00000000"; // end time in GPST
    
    //for galileo
    //argv[4] = "2015 1 11 18 33 04.00000000"; // start time in GPST
    //argv[5] = "2015 1 11 18 34 28.00000000"; // end time in GPST
    
    //for bds
    argv[4] = "2015 01 01 00 00 00.00000000"; // start time in GPST
    argv[5] = "2015 12 31 00 00 28.00000000"; // end time in GPST
    
    argv[6]="sp3/codsp3-2015/com18264.sp3,sp3/codsp3-2015/com18265.sp3,sp3/codsp3-2015/com18266.sp3,sp3/codsp3-2015/com18270.sp3,sp3/codsp3-2015/com18271.sp3";
    
    argv[6]="sp3/whusp3-2015/wum18253.sp3,sp3/whusp3-2015/wum18254.sp3,sp3/whusp3-2015/wum18255.sp3,sp3/whusp3-2015/wum18256.sp3,sp3/whusp3-2015/wum18260.sp3";
    
    
    argv[7] = "60";  // interval
    
    
    //argv[4] = "2015 1 9 09 13 35.00000000"; // start time in GPST
    //argv[5] = "2015 1 9 10 12 32.00000000"; // end time in GPST
    
    
    
    /*
    argv[6] = "sp3/codsp3-2015/com18342.sp3,sp3/codsp3-2015/com18343.sp3,sp3/codsp3-2015/com18344.sp3,sp3/codsp3-2015/com18345.sp3,sp3/codsp3-2015/com18346.sp3,sp3/codsp3-2015/com18350.sp3,sp3/codsp3-2015/com18351.sp3,sp3/codsp3-2015/com18352.sp3,sp3/codsp3-2015/com18353.sp3,sp3/codsp3-2015/com18354.sp3,sp3/codsp3-2015/com18355.sp3";
    */
   
     
     
     
    //
   
    
    double interval = atof(argv[7]);
    
    GString configfile = argv[1];
    
    GConfigure config;
    config.parseCfg(configfile);
    
    GForceModelMgr myforceModelManager;
    config.configForceModelMgr(myforceModelManager);
    
    GSpacecraftModelMgr::initialiseModel( config.config.spacecraftmodel);
    GSpaceCraftMgr::Initializer();
    
    GJPLEPH::loadEphFile_a(config.config.planetEphemeris);
    
    GEarthOrientationParameter::loadEOP(config.config.eopfile);
    
    GString satsys = argv[2];
    int satprn = atoi(argv[3]);
    
    
    CivilTime ct0("tsGPS",argv[4]);   // 12
    CivilTime ct1("tsGPS",argv[5]); // 03, 15, 21
    
    
    GString sp3filelist = argv[6];
    
    std::vector<GString> sp3file = sp3filelist.split(',');
    for(int i = 0 ; i< sp3file.size(); i++ )
    {
        //GSpaceCraftMgr::loadPreciseEphemeris(sp3file[i]);
    }
    
    
    GSensorID myid(satsys, satprn); // C08  C06 and C14 , G11, E11(IOV101)
    GSpaceCraft& mysat = GSpaceCraftMgr::gSpacecraft[myid.getIDString()];
    
    //test the eclipse shadow function and orbit prediction
    
    GTime start_gps, end_gps;
    start_gps.SetFromCivilTime(ct0);
    end_gps.SetFromCivilTime(ct1);
    
    GTime start_utc = GTime::GPST2UTC(start_gps);
    GTime epoch_end_fit = start_utc;
    GTime epoch_gps = start_gps;
    
    
    //C02 at GPST 2015 01 01 00 00 00
    //GKeplerianElements mykpe(42165.1079788486379485,0.000307053530371378111726,0.00302151262482347904986,5.03516474942856095964,5.69774810696488764883,4.97713242928925980735);
    
    //C03 at GPST 2015 01 01 00 00 00
    //GKeplerianElements mykpe(42167.5526423147702673,0.000393354643937860229928,0.0253567136153691474743,4.88131200785190366709,0.449903337415982762515,4.62621302131688660886);
    
    //C04 at GPST 2015 01 01 00 00 00
    //GKeplerianElements mykpe(42166.2865181570873823,0.000890929908124252602566,0.0144140069034237386369,3.0385097843287671715,0.820218081949606880876,0.681508205399564470639);
    
    //C05 at GPST 2015 01 01 00 00 00
    //GKeplerianElements mykpe(42165.7107659391173193,0.000162068045336555684134,0.0182363320866512149121,3.8327760018172472023,0.371078604031111425332,4.85052711614888831804);
    
    //C06 at GPST 2015 01 01 00 00 00
    //GKeplerianElements mykpe(42162.4084826060458617,0.00353858487075050555382,0.947160032610789626161,3.44185201563725939877,3.5289183367311238726,3.15337259842582584568);
    
    //C07 at GPST 2015 01 01 00 00 00
    //GKeplerianElements mykpe(42158.2407671668466804,0.00310429849497752870929,0.944603271933532329108,3.5561210405044141325,5.61097587867011604601,0.940171369608793661143);
    
    //C08 at GPST 2015 01 01 00 00 00
    //GKeplerianElements mykpe(42170.1987404618108393,0.00250471925704411632965,0.992634997874152325003,3.30364204511180314938,1.43493153974322651401,5.32024385219189188678);
    
    //C09 at GPST 2015 01 01 00 00 00
    //GKeplerianElements mykpe(42166.3119610689637646,0.00306045314881421256017,0.952694972555786168122,3.35872720362736876254,3.56769170141733082957,2.78153374241570761072);
    
    //C10 at GPST 2015 01 01 00 00 00
    GKeplerianElements mykpe(42158.5455173419890933,0.00297678789158783898422,0.946320613861091786134,3.54313401237741487648,5.60117667874381430781,0.521560853497004761614);
    
    //C11 at GPST 2015 01 01 00 00 00
    //GKeplerianElements mykpe(27904.7761430916347916,0.00240232320729961978945,0.973095890700745203005,3.46758595575894412733,1.46975968471844375962,0.80975252861510960134);
    
    //C12
    //GKeplerianElements mykpe(27904.4660136779959139,0.00274336440899936029325,0.971964838760725591271,3.36304948820225702022,1.46073050121603542628,1.7173513865787128907);
    
    //C14
    //GKeplerianElements mykpe(27905.3627291815599456,0.00156609722616517493127,0.956695310707969204244,3.56007288522406373232,3.54991651616790226953,4.91337390378707183913);
    
    
    while( 1 )
    {
        GTime epoch_utc = GTime::GPST2UTC(epoch_gps);
        //GTime epoch_utc = epoch_end_fit;
        
        JDTime jdt = GTime::GTime2JDTime(epoch_utc);
        CivilTime ct = GTime::JDTime2CivilTime(jdt);
        
        GVector mp_ecef, mv_ecef, mp_eci, mv_eci;
        //orbit comparison
       
        GSpaceEnv::updateSpaceEnvironment(epoch_utc);
        
        
       // mysat.getEphValue_test(epoch_gps).getPV(mp_ecef,mv_ecef);
       // mysat.getStatePointer()->updateState_ecef(epoch_utc, mp_ecef, mv_ecef);
        
        
        mykpe.propagate(mykpe, (epoch_gps-start_gps).toSeconds(), mp_eci, mv_eci);
        
        
        // for the eclipse event detection
        //mysat.getStatePointer()->shadow_detector(GSpaceEnv::planetPos_ecef[GJPLEPH::SUN], mp_ecef);
        
        //for the eci cooridnates
        mysat.getStatePointer()->shadow_detector(epoch_gps,GSpaceEnv::planetPos_eci[GJPLEPH::SUN], mp_eci);
        
        
        // output the shadow function
        //double shadow_factor = GMotionState::shadowFactor_SECM(true,GSpaceEnv::planetPos_ecef[GJPLEPH::SUN] , mp_ecef);
        
        
        //printf("%s %f\n",GTime::GTime2CivilTime(epoch_gps).TimeString().c_str(), shadow_factor );
        
        
        //mysat.getStatePointer()->updateState_ecef(epoch_utc, mp_ecef, mv_ecef);
        
        if( ( epoch_gps = epoch_gps + interval) >= end_gps)
        {
            break;
        }
        
        //obscount++;
    }

    
    
    return 0;
}



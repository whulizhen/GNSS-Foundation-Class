//
//  testERP.cpp
//  GFC
//
//  Created by lizhen on 16/2/4.
//  Copyright © 2016年 lizhen. All rights reserved.
//

#include <stdio.h>
#include <algorithm>
#include <typeinfo>


#include "GEarthOrientationParameter.hpp"
#include "GRungeKutta.hpp"

#include "GFMEarthRadiationPressure.hpp"
#include "GFMThermalRadiationForce.hpp"

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

int main( int argc,char* argv[] )
{
    
    // testing data
    argv[1] = "gfcsetup.cfg";
    argv[2] = "ssGAL";  //satellite system;
    argv[3] = "11";     //prn
    
    
    //for orbit fit test,comparing with PANDA
    //argv[4] = "2016 10 14 0 0 0.00000000"; // start time in GPST
    //argv[5] = "2016 10 17 0 0 0.00000000"; // end time in GPST
    
    argv[4] = "2015 3 12 0 0 0.00000000"; // start time in GPST
    argv[5] = "2015 3 18 0 0 0.00000000"; // end time in GPST
    
    
    //argv[6] = "sp3/comsp3-2016/com19185.sp3";
    //argv[6] = "sp3/comsp3-2016/com19184.sp3,sp3/comsp3-2016/com19185.sp3,sp3/comsp3-2016/com19186.sp3,sp3/comsp3-2016/com19190.sp3,sp3/comsp3-2016/com19191.sp3,sp3/comsp3-2016/com19192.sp3,sp3/comsp3-2016/com19193.sp3";
    
    
   argv[6] = "sp3/whusp3-2015/wum18353.sp3,sp3/whusp3-2015/wum18354.sp3,sp3/whusp3-2015/wum18355.sp3,sp3/whusp3-2015/wum18356.sp3,sp3/whusp3-2015/wum18360.sp3,sp3/whusp3-2015/wum18361.sp3,sp3/whusp3-2015/wum18362.sp3,sp3/whusp3-2015/wum18363.sp3,sp3/whusp3-2015/wum18364.sp3,sp3/whusp3-2015/wum18365.sp3,sp3/whusp3-2015/wum18366.sp3,sp3/whusp3-2015/wum18370.sp3";
    
    //non eclipse test
    //argv[4] = "2015 3 8 3 0 0.00000000"; // start time in GPST
    //argv[5] = "2015 3 10 3 0.00000000"; // end time in GPST
    //argv[6] = "sp3/codsp3-2015/com18342.sp3,sp3/codsp3-2015/com18343.sp3,sp3/codsp3-2015/com18344.sp3,sp3/codsp3-2015/com18345.sp3,sp3/codsp3-2015/com18346.sp3,sp3/codsp3-2015/com18350.sp3,sp3/codsp3-2015/com18351.sp3,sp3/codsp3-2015/com18352.sp3,sp3/codsp3-2015/com18353.sp3,sp3/codsp3-2015/com18354.sp3,sp3/codsp3-2015/com18355.sp3";
    
    
    
    //2014 data
    
    //argv[4] = "2014 4 2 4 0 0.00000000";  //start time in GPST
    //argv[5] = "2014 4 3 4 0 0.00000000";  //end time in GPST
    //argv[6] = "sp3/codsp3-2014/com17862.sp3,sp3/codsp3-2014/com17863.sp3,sp3/codsp3-2014/com17864.sp3,sp3/codsp3-2014/com17865.sp3,sp3/codsp3-2014/com17866.sp3,sp3/codsp3-2014/com17870.sp3,sp3/codsp3-2014/com17871.sp3,sp3/codsp3-2014/com17872.sp3,sp3/codsp3-2014/com17873.sp3,sp3/codsp3-2014/com17874.sp3,sp3/codsp3-2014/com17875.sp3,sp3/codsp3-2014/com17876.sp3";
    
    
    
    //2004 data
    //argv[4] = "2004 2 28 4 0 0.00000000"; // start time in GPST
    //argv[5] = "2004 2 29 4 0 0.00000000"; // end time in GPST
    //argv[6] = "sp3/jplsp3-2004/jp212596.sp3,sp3/jplsp3-2004/jp212600.sp3,sp3/jplsp3-2004/jp212601.sp3,sp3/jplsp3-2004/jp212602.sp3,sp3/jplsp3-2004/jp212603.sp3,sp3/jplsp3-2004/jp212604.sp3,sp3/jplsp3-2004/jp212605.sp3,sp3/jplsp3-2004/jp212606.sp3,sp3/jplsp3-2004/jp212610.sp3,sp3/jplsp3-2004/jp212611.sp3,sp3/jplsp3-2004/jp212612.sp3,sp3/jplsp3-2004/jp212613.sp3";
    
    
    /*
    //eclipse orbit prediction test
    argv[4] = "2015 1 11 17 34 58.00000000"; // start time in GPST
    argv[5] = "2015 1 11 17 36 58.00000000"; // end time in GPST
    argv[6] ="sp3/codsp3-2015/com18264.sp3,sp3/codsp3-2015/com18265.sp3,sp3/codsp3-2015/com18266.sp3,sp3/codsp3-2015/com18270.sp3,sp3/codsp3-2015/com18271.sp3,sp3/codsp3-2015/com18272.sp3,sp3/codsp3-2015/com18273.sp3";
    */
    
    argv[7] = "1";  // 0 for without fitting, 1 for with fitting
    argv[8] = "48";  // the hours for orbit fitting measured from the start time
    argv[9] = "./";  // the directory of the log files
    
    /*
    for(int i = 1 ; i<=9; i++ )
    {
        printf("%s\n", argv[i]);
    }
    */
    
    bool with_helmert = false;
    
    GString configfile = argv[1];
    
    GConfigure config;
    config.parseCfg(configfile);
    
    
    
    GForceModelMgr myforceModelManager;
    config.configForceModelMgr(myforceModelManager);
    
    //printf("read spacecraft1 %s\n", config.config.spacecraftmodel.c_str());
    GSpacecraftModelMgr::initialiseModel( config.config.spacecraftmodel);
    //printf("read spacecraft2 %s\n", config.config.spacecraftmodel.c_str());
    GSpaceCraftMgr::Initializer();
    //printf("read spacecraft3 \n");
    
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
        GSpaceCraftMgr::loadPreciseEphemeris(sp3file[i]);
    }
    
    
    bool withfit = false;
    double orbitfittingSecond =0.0;
    if(atoi(argv[7]) != 0 )
    {
        withfit = true;
    }
    
    GString logfilepath(argv[9]);
    
    if(withfit == true)
    {
        orbitfittingSecond = atof(argv[8]);
        orbitfittingSecond = orbitfittingSecond*3600.0;
    }
    
    
//    CivilTime testct(2004,  2, 29,  21,  0,  0.00000000, "tsGPS");
//    GTime testgt;
//    testgt.SetFromCivilTime(testct);
//    GTime testutc =  GTime::GPST2UTC(testgt);
//    GEarthOrientationParameter myeop;
//    myeop.setEpochTime(testutc);
//    GTime TAI = GTime::UTC2TAI(testutc);
//    GTime TT = GTime::TAI2TT(TAI);
//    GVector mypos_eci(-11640.789627773,23502.285387444,3163.7370596617);
//    GVector mypos_ecef;
//    myeop.ECI2ECEF_pos(mypos_eci, mypos_ecef);
    
    //GSpaceCraftMgr::loadPreciseEphemeris("sp3/ilrsb.orb.lageos2.160416.v35.sp3"); // wum18000.sp3,96,2014,4,6
    //GSpaceCraftMgr::loadPreciseEphemeris("sp3/ilrsb.orb.lageos2.160416.v35.sp3"); // wum18000.sp3,96,2014,4,6
    //GSpaceCraftMgr::loadPreciseEphemeris("sp3/ilrsb.orb.lageos2.160423.v35.sp3"); // wum18000.sp3,96,2014,4,6
    //GSpaceCraftMgr::loadPreciseEphemeris("sp3/ilrsb.orb.lageos2.160430.v35.sp3"); // wum18000.sp3,96,2014,4,6
    
    
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
    
  GSensorID myid(satsys, satprn); // C08  C06 and C14 , G11, E11(IOV101)
  GSpaceCraft& mysat = GSpaceCraftMgr::gSpacecraft[myid.getIDString()];
  
    
   
    
 //non eclipse season
 // CivilTime ct0(2015, 2, 14, 20, 0, 0.00000000, "tsGPS");   // 12
 // CivilTime ct1(2015, 2, 15, 20, 0, 0.00000000, "tsGPS"); // 03, 15, 21
    
    //for eclipse season(from full phase to penumbra)
    //CivilTime ct0(2015, 1, 11, 17, 34, 43, "tsGPS");   // 12
    //CivilTime ct1(2015, 1, 11, 17, 36, 5, "tsGPS"); // 03, 15, 21
    
    // from penumbra to full phase
    //CivilTime ct0(2015, 1, 11, 18, 33, 04, "tsGPS");   // 12
    //CivilTime ct1(2015, 1, 11, 18, 34, 24, "tsGPS"); // 03, 15, 21
    
  GTime start_gps, end_gps;
  start_gps.SetFromCivilTime(ct0);
  end_gps.SetFromCivilTime(ct1);
  
    
  //ofstream dumpFile;
  //dumpFile.open("attitude.txt");
  
  GTime start_utc = GTime::GPST2UTC(start_gps);
  GTime epoch_end_fit = start_utc;
    
//    
//    double tt0 = 0;
//    GEarthOrientationParameter eop;
//    
//    CivilTime test_ct(2004,2,28,20,59,47,"tsUTC");
//    GTime test_gt = GTime::CivilTime2GTime(test_ct);
//    eop.setEpochTime(test_gt);
//    GVector testp_eci(-17755.43751	,4917.256218,	-19177.88171);
//    GVector testp_ecf(11441.354506, 14431.794210,	-19184.494766);
//    //eop.ECI2ECEF_pos(testp_eci, testp_ecf);
//    eop.ECEF2ECI_pos(testp_ecf, testp_eci);
    
//    while(1)
//    {
//        GTime test1 = start_gps + tt0;
//        GVector ps1,vs1, ps2, vs2;
//        GTime utc = GTime::GPST2UTC(test1);
//        eop.setEpochTime(utc);
//        
//        mysat.getEphValue(test1).getPV(ps1,vs1);  // these values are in ECEF
//        
//        eop.ECEF2ECI_pos(ps1, ps2);
//        eop.ECEF2ECI_vel(ps1, vs1, vs2);
//        
//        
//        
//        //printf("%s ,ecef: %12.6f %12.6f %12.6f %12.10f %12.10f %12.10f \n", GTime::GTime2CivilTime(test1).TimeString().c_str(), ps1.x, ps1.y, ps1.z, vs1.x, vs1.y, vs1.z);
//        //printf("%s ,eci : %12.6f %12.6f %12.6f %12.10f %12.10f %12.10f \n", GTime::GTime2CivilTime(utc).TimeString().c_str(), ps2.x, ps2.y, ps2.z, vs2.x, vs2.y, vs2.z);
//        
//        printf("%s ,ecef: %12.6f %12.6f %12.6f %12.6f %12.6f %12.6f \n", GTime::GTime2CivilTime(test1).TimeString().c_str(), ps1.x, ps1.y, ps1.z, ps2.x, ps2.y, ps2.z);
//        
//        tt0 = tt0 + 900;
//        
//        if( test1 >= end_gps)
//        {
//            break;
//        }
//    }
//    
    
    // these values are from sp3
    // GPS11 : 4945.8532142566091, 18404.204773754933 ,-18916.739350034703, -2.959070091997599,2.0646052197117255,1.289437698094192
    // BDS06(1320.0) : -1544.5056131903586, -26534.630773075991 , 32819.347490316803, 2.8840657859359182,0.75414467192697077,0.73305622299686313
    // BDS07 : 13612.859753915565 , -31442.781106099192, -24669.386730980459, 2.5223793718360801,-0.27256258875167894,1.7298969720692858
    // BDS08(1375.0) : 2618.7657220011865, -40902.351169073059, -9464.9895366924557, 1.7056035866900601, 0.68009850312978704, -2.4753015289771456
    // BDS11 : -9641.4958189231002,23193.605609948037, 12301.808925571615, -1.6980793230296116, -2.1030715656246803, 2.6294043824888749
    
    GVector ps,vs;
    
    //mysat.getEphValue(start_gps).getPV(ps,vs);  // these values are in ECEF
    
    mysat.getEphValue_test(start_gps).getPV(ps,vs);  // these values are in ECEF
    
    
    GSpaceEnv::updateSpaceEnvironment(start_utc);
    mysat.getStatePointer()->updateState_ecef(start_utc, ps, vs);
    
    
//    //-3821.2011680740902  -27935.652963677141 9036.4480460475643 2.0509323814961515 -1.1860677864789158 -2.8010750227933454
  //  ps.x = 1998.262166 ; ps.y = 26807.029665, ps.z = -12370.888483;
  //  vs.x = -2.2176139989; vs.y = 1.36088755712; vs.z = 2.5893515261;
    
 //   ps.x = 7000 ; ps.y = -12124, ps.z = 0.0;
 //   vs.x = 2.6679; vs.y = 4.621; vs.z = 0.0;
    
   // mysat.getStatePointer()->updateState_eci(start_utc, ps, vs);
    
 //   GVector testpp, testvv;
 //   mysat.getStatePointer()->keplerianElement.propagate(mysat.getStatePointer()->keplerianElement, 18000.0, testpp, testvv);
    
//      GKeplerianElements mykpe;
//      mykpe.PV2KP(ps, vs);
//      GVector testp, testv;
//      GKeplerianElements::propagate(mykpe,3600, testp, testv);
    
    
    DOYTime doy =  GTime::CivilTime2DOYTime(ct0);
    char tmp[1024]={0};
    sprintf(tmp, "%04d%03d",doy.m_year, doy.m_doy);

    GOrbitFitting orbfit(&mysat);
    
     //for initial orbit determination
    if(withfit == true)
    {
        
        orbfit.setForceManager(myforceModelManager);
        
        orbfit.setStepSize(config.config.stepsize);
        orbfit.setInitialValue(start_utc, mysat.getStatePointer()->satpos_eci, mysat.getStatePointer()->satvel_eci);
        epoch_end_fit = start_utc + orbitfittingSecond;  ////use the first 3 hours to adjust the inital value
        orbfit.setEndEpoch( epoch_end_fit );
        
        GString logfilename = logfilepath + "orbit_fit_state";
        logfilename += GString(tmp);
        logfilename += GString(".log");
        orbfit.setOutputFile(logfilename);
        
        //orbfit.calculateInitialError();
        
        orbfit.calculate_SRPModel();
        
        myforceModelManager = orbfit.getForceManager();
    }
    
    
   
    
    /*
    ps.x = -305.15368350834297;
    ps.y = 29515.393410085719;
    ps.z = -2199.2052776956534;
    vs.x = 0.051205519643810593;
    vs.y = 0.22453399779164063;
    vs.z = 3.002103597335366;
    */
    
    //GTime testUTC = epoch_utc + 900;
    //mysat.getStatePointer()->updateState_ecef(epoch_utc, ps, vs);
    
    //eclipse model test
    
//    epoch_gps = epoch ; // GPST
//    while(1)
//    {
//         if( epoch_gps > endTime )
//         {
//             break;
//         }
//        epoch_gps = epoch_gps + 30;
//        epoch_utc = GTime::GPST2UTC(epoch_gps);
//        
//        mysat.getEphValue(epoch_gps).getPV(ps,vs);  // these values are in ECEF
//        
//        GSpaceEnv::updateSpaceEnvironment(epoch_utc);
//        mysat.getStatePointer()->updateState_ecef(epoch_utc, ps, vs);
//        
//    }
    
    GTime epoch_gps = GTime::UTC2GPST(epoch_end_fit);
    
    //GTime epoch_gps = start_gps;
    
    bool ispredicted = false;
    
    if(epoch_gps < end_gps)
    {
        ispredicted = true;
    }
    
    if(ispredicted ==  true)
    {
        
        GOrbitPredictor op( &mysat);
        
        op.setStepsize(config.config.stepsize);
        
        if(logfilepath != "")
        {
            GString logfilename = logfilepath + "orbit_prediction_state";
            logfilename += GString(tmp);
            logfilename += GString(".log");
            op.setLogfile(logfilename);
        }
        
        
        //reset the force manager for orbit predictor
        op.setForceManager(myforceModelManager);
        
        std::vector<GVector> rms_data;
        rms_data.reserve(10000);
        
        //for helmert transformation
        std::vector<GVector> predictedOrbit;
        std::vector<GVector> preciseOrbit;
        std::vector<GVector> preciseVel;
        std::vector<CivilTime>   timeUTC;
        
        double helmertParam[7]={0.0};
        
        GVector mean_data;
        int obscount = 0;
        double interval = 900.0;  // 15 minutes
        double t0 = clock();
        //double tt = 0.0;
        //GVector mytestp, mytestv;
        while( epoch_gps < end_gps  )
        {
            GTime epoch_utc = GTime::GPST2UTC(epoch_gps);
            
            JDTime jdt = GTime::GTime2JDTime(epoch_utc);
            CivilTime ct = GTime::JDTime2CivilTime(jdt);
            
            op.PropagateTo(epoch_utc);
            
            GVector mp_ecef, mv_ecef;
            
            //orbit comparison
            mysat.getEphValue_test(epoch_gps).getPV(mp_ecef,mv_ecef);
            
            GVector mp_eci = mysat.getStatePointer()->satpos_eci;
            
            GVector mv_eci = mysat.getStatePointer()->satvel_eci;
            
            GVector p1,v1, p2,v2;
            
            GSpaceEnv::eop.ECEF2ECI_pos(mp_ecef, p1);
            
            GSpaceEnv::eop.ECEF2ECI_vel(mp_ecef, mv_ecef, v1);
            
            GSpaceEnv::eop.ECI2ECEF_pos(mp_eci, p2);
            
            GSpaceEnv::eop.ECI2ECEF_vel(mp_eci, mv_eci, v2);
            
            predictedOrbit.push_back(mp_eci);
            preciseOrbit.push_back(p1);
            preciseVel.push_back(v1);
            timeUTC.push_back(ct);
            
            //printf("%s, %.8f, %.8f, %.8f, %.8f, \n",ct.TimeString().c_str(), diff_rtn.x,diff_rtn.y,diff_rtn.z ,diff_rtn.norm());
            
            //printf("%s, %.8f, %.8f, %.8f, %.8f, %.8f, %.8f \n",ct.TimeString().c_str(),p1.x,p1.y,p1.z ,v1.x,v1.y,v1.z);
            
            //printf("%s, %.8f, %.8f, %.8f, %.8f, %.8f, %.8f \n",ct.TimeString().c_str(),mp_eci.x,mp_eci.y,mp_eci.z ,mv_eci.x,mv_eci.y,mv_eci.z);
            
            //tt+= interval;
            
            if( ( epoch_gps = epoch_gps + interval) >= end_gps)
            {
                break;
            }
            
            obscount++;
        }
        
        
        double t1 = clock();
        
        double time = (t1-t0)/((double) CLOCKS_PER_SEC);
        
        //did the helmert transformation
        if(with_helmert == true)
        {
            GMath::HelmertParameter(predictedOrbit.size(), &preciseOrbit[0], &predictedOrbit[0], helmertParam);
        }
        
        
        for(int i = 0 ; i< predictedOrbit.size(); i++ )
        {
            GVector tp = predictedOrbit[i],tv;
            
            if(with_helmert)
            {
                GMath::HelmertTransform(predictedOrbit[i], helmertParam, tp);
            }
            
            GVector diff_rtn =  op.orbitdiff( preciseOrbit[i], preciseVel[i], tp, tv);
            diff_rtn = diff_rtn*1000.0;
            
            rms_data.push_back(diff_rtn);
            
            //printf("%s, %.8f, %.8f, %.8f, %.8f, \n",timeUTC[i].TimeString().c_str(), diff_rtn.x,diff_rtn.y,diff_rtn.z ,diff_rtn.norm());
        }
        
        
        //std::cout<< "time:" << time <<" second" << std::endl;
        
        //mean_data = mean_data/rms_data.size();
        //compute the rms in RTN direction
        double rms_rtn[4] = {0.0};
        for( int i = 0 ; i< rms_data.size(); i++ )
        {
            //rms_rtn[0] += (rms_data[i]-mean_data).x * (rms_data[i]-mean_data).x;
            //rms_rtn[1] += (rms_data[i]-mean_data).y * (rms_data[i]-mean_data).y;
            //rms_rtn[2] += (rms_data[i]-mean_data).z * (rms_data[i]-mean_data).z;
            
            
            rms_rtn[0] += (rms_data[i]).x * (rms_data[i]).x;
            rms_rtn[1] += (rms_data[i]).y * (rms_data[i]).y;
            rms_rtn[2] += (rms_data[i]).z * (rms_data[i]).z;
            
            rms_rtn[3] += rms_data[i].norm()*rms_data[i].norm();
        }
        
        rms_rtn[0] = sqrt(rms_rtn[0]/rms_data.size());
        rms_rtn[1] = sqrt(rms_rtn[1]/rms_data.size());
        rms_rtn[2] = sqrt(rms_rtn[2]/rms_data.size());
        rms_rtn[3] = sqrt(rms_rtn[3]/rms_data.size());
        
        printf(" %s RMS in RAC: %.6f %.6f %.6f %.6f ",ct0.TimeString().c_str(), rms_rtn[0], rms_rtn[1], rms_rtn[2], rms_rtn[3]);
        if(withfit == true)
        {
            printf("%.3f \n", orbfit.getRMS());
        }
        else
        {
            printf("\n");
        }
        
    }
    
    
    
    
  return 0;
}

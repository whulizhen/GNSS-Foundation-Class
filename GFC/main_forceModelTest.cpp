//
//  forceModelTest.cpp
//  GFC
//
//  Created by lizhen on 06/10/2016.
//  Copyright © 2016 lizhen. All rights reserved.
//

#include <stdio.h>
#include <stdio.h>
#include <algorithm>
#include <typeinfo>
#include "GConfigure.hpp"


#include "GEarthRadiationModel.h"
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

using namespace gfc;

int main(int argc, char* argv[])
{
    
    // testing data
    argv[1] = "gfcsetup.cfg";
    argv[2] = "ssGAL";  //satellite system;
    argv[3] = "11";     //prn
    
    
    //for orbit fit test,comparing with PANDA
    //argv[4] = "2016 10 14 0 0 0.00000000"; // start time in GPST
    //argv[5] = "2016 10 17 0 0 0.00000000"; // end time in GPST
    
    argv[4] = "2015 3 12 0 0 0.00000000"; // start time in GPST
    argv[5] = "2015 3 13 0 0 0.00000000"; // end time in GPST
    
    
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
    
    //argv[7] = "1";  // 0 for without fitting, 1 for with fitting
    //argv[8] = "168";  // the hours for orbit fitting measured from the start time
    //argv[9] = "./";  // the directory of the log files
    
    /*
     for(int i = 1 ; i<=9; i++ )
     {
     printf("%s\n", argv[i]);
     }
     */
    
   
    
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
        GSpaceCraftMgr::loadPreciseEphemeris(sp3file[i]);
    }
    
    GSensorID myid(satsys, satprn); // C08  C06 and C14 , G11, E11(IOV101)
    
    GSpaceCraft& mysat = GSpaceCraftMgr::gSpacecraft[myid.getIDString()];
    
     //the following code is for force model accelerations
    
        GTime t0 = GTime::CivilTime2GTime(ct0);
        GTime t1 = GTime::CivilTime2GTime(ct1);
    
        GTime  t = t0;
        GPreciseEphemeris peph;
        GVector p, v, p_eci,v_eci;
    
        GFMEarthRadiationPressure erpForce;
        GFMSolarRadiationPressure* srpForce= ((GFMSolarRadiationPressure*)myforceModelManager.m_forceModels["GFMSRP"]);
    
    
        GVector force_eci, force_RTN, force_bfs, force_dyb, force_test;
        double interval = 30.0;
        double pi = GCONST("PI");
        FILE* testFile = fopen("forceTest.txt","w+");
    
        while( t < t1 )
        {
            peph = mysat.getEphValue_test(t);
            peph.getPV(p, v);
            GTime utc = GTime::GPST2UTC(t);
            GSpaceEnv::updateSpaceEnvironment(utc);
            mysat.getStatePointer()->updateState_ecef(utc, p, v);
            
            GSpaceEnv::eop.ECEF2ECI_pos(p, p_eci);
            GSpaceEnv::eop.ECEF2ECI_vel(p, v, v_eci);
            
            
            
            //erpForce.doCompute( mysat.getStatePointer(), mysat.getSpaceCraftGemotry());
            //GVector  lw_eci = mysat.getStatePointer()->earthFlux.totalFlux_lw;
            //GVector  sw_eci = mysat.getStatePointer()->earthFlux.totalFlux_sw;
            //force_eci = erpForce.getForce();
            
            srpForce->doCompute(&mysat);
            force_eci = srpForce->getForce();
            
            double solarflux = mysat.getStatePointer()->solarFlux.m_flux.m_shortwave+mysat.getStatePointer()->solarFlux.m_flux.m_longwave;
            
            
            double beta = mysat.getStatePointer()->attitude_eci.beta;
            
            double lat = mysat.getStatePointer()->attitude_eci.phi;
            double lon = mysat.getStatePointer()->attitude_eci.lambda;
            
            double du = mysat.getStatePointer()->attitude_eci.eta;
            
            double eps = mysat.getStatePointer()->eps;
            
            GVector R =  mysat.getStatePointer()->attitude_eci.eR;
            GVector N =  mysat.getStatePointer()->attitude_eci.eN;
            GVector T =  mysat.getStatePointer()->attitude_eci.eT;
            
            GVector D = mysat.getStatePointer()->attitude_eci.eD;
            GVector Y = mysat.getStatePointer()->attitude_eci.eY;
            GVector B = mysat.getStatePointer()->attitude_eci.eB;
            
            
            
            GVector bfs_x = mysat.getStatePointer()->attitude_eci.xhat;
            GVector bfs_y = mysat.getStatePointer()->attitude_eci.yhat;
            GVector bfs_z = mysat.getStatePointer()->attitude_eci.zhat;
            
            
            force_RTN.x = force_eci.x * R.x + force_eci.y * R.y + force_eci.z * R.z;
            force_RTN.y = force_eci.x * T.x + force_eci.y * T.y + force_eci.z * T.z;
            force_RTN.z = force_eci.x * N.x + force_eci.y * N.y + force_eci.z * N.z;
            
            force_bfs.x = force_eci.x * bfs_x.x + force_eci.y *  bfs_x.y + force_eci.z *  bfs_x.z;
            force_bfs.y = force_eci.x * bfs_y.x + force_eci.y *  bfs_y.y + force_eci.z *  bfs_y.z;
            force_bfs.z = force_eci.x * bfs_z.x + force_eci.y *  bfs_z.y + force_eci.z *  bfs_z.z;
            
            force_dyb.x = D.x*force_eci.x + D.y*force_eci.y + D.z*force_eci.z;
            force_dyb.y = Y.x*force_eci.x + Y.y*force_eci.y + Y.z*force_eci.z;
            force_dyb.z = B.x*force_eci.x + B.y*force_eci.y + B.z*force_eci.z;
            
            force_test.x = -cos(lat)*force_bfs.x + sin(lat)*force_bfs.z;
            force_test.y = -force_bfs.y;
            force_test.z = sin(lat)*force_bfs.x + cos(lat)*force_bfs.z;
            
            
            
            /*
            GVector lw_rtn, sw_rtn;
            lw_rtn.x = lw_eci.x * R.x + lw_eci.y * R.y + lw_eci.z * R.z;
            lw_rtn.y = lw_eci.x * T.x + lw_eci.y * T.y + lw_eci.z * T.z;
            lw_rtn.z = lw_eci.x * N.x + lw_eci.y * N.y + lw_eci.z * N.z;
            
            sw_rtn.x = sw_eci.x * R.x + sw_eci.y * R.y + sw_eci.z * R.z;
            sw_rtn.y = sw_eci.x * T.x + sw_eci.y * T.y + sw_eci.z * T.z;
            sw_rtn.z = sw_eci.x * N.x + sw_eci.y * N.y + sw_eci.z * N.z;
            */
            
            
            t = t + interval;
            
           
            
            //double test1 = dotproduct(mysat.getStatePointer()->attitude.phat, mysat.getStatePointer()->satposHat_eci);
            //double test2 = acos(test1)*180.0/pi;
            //double test3 = dotproduct( mysat.getStatePointer()->attitude.phat, mysat.getStatePointer()->attitude.n );
            //double test4 = acos(test3)*180.0/pi;
            
            /*
            printf("%s %.15f %.15f %.15f %.8f %.8f %.8f\n",
                   GTime::GTime2CivilTime(t).TimeString().c_str(),
                   force_RTN.x,force_RTN.y,force_RTN.z, beta, lat, lon
                  );
            */
            /*
            printf("%s %.8f %.15f %.15f %.15f %.15f %.8f %.8E %.8E %.8E\n",
                   GTime::GTime2CivilTime(t).TimeString().c_str(),solarflux,
                   beta*180/pi, du*180/pi,eps*180/pi,lat*180/pi, lon*180/pi,force_bfs.x,force_bfs.y,force_bfs.z
                   );
            */
            
            /*
            printf("%.8f %.15f %.15f %.15f %.15f %.8E %.8E %.8E\n",
                   //GTime::GTime2CivilTime(t).TimeString().c_str(),
                   solarflux,
                   beta*180/pi, du*180/pi,eps*180/pi,lat*180/pi,force_bfs.x,force_bfs.y,force_bfs.z
                   );
            */
//            
//            fprintf(testFile, "%s, %.4f,%.4f,%12.6f,%12.6f, %12.6f, %8.6f, %8.6f, %8.6f \n", GTime::GTime2CivilTime(t).TimeString().c_str(),
//                                                    mysat.getStatePointer()->beta*180/3.14159265357,
//                                                    mysat.getStatePointer()->eps*180/3.14159265357,
//                                                    p.x, p.y, p.z,
//                                                    v.x, v.y, v.z
//                                                    );
            
            
        }
    
        fclose(testFile);
    
}


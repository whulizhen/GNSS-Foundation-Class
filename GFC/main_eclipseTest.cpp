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

#include "GQuaternion.hpp"

//#include "GEarthRadiationModel.h"

using namespace gfc;


//load the grace accelerationa and attitude data
void loadGRACE_acc_sca(GString filename, GSpaceCraft& mysatA, GSpaceCraft& mysatB)
{
    std::fstream infile(filename.c_str());
    if ( infile.fail() )
    {
        std::cerr << "\nCould not open GRACE ephemeris file. Terminating...\n";
        exit(0);
    }
    const int MAX = 1024;
    char store[MAX]={0};
    long lineNumber = 0;
    int totalsat = 0;
    
    GVector mp_ecef, mv_ecef, mp_eci, mv_eci;
    
    // the body fixed frame of graceA
    GVector xbfs, ybfs, zbfs;
    
    while ( !infile.eof() )
    {
        
        infile.getline(store, MAX); //read past the file header
        GString reader(store);
        if(reader == "" )
        {
            break;
        }
        
        std::vector<GString> splitstr = reader.split();
        CivilTime ct;
        //default timesystem is GPS time
        ct.m_ts = GTimeSystem("tsGPS");
        ct.m_year = static_cast<int>(splitstr[0].asINT() ) ;
        ct.m_month = static_cast<int>( splitstr[1].asINT());
        ct.m_day = static_cast<int>(splitstr[2].asINT()) ;
        ct.m_hour = static_cast<int>(splitstr[3].asINT()) ;
        ct.m_minute = static_cast<int>(splitstr[4].asINT()) ;
        ct.m_second = splitstr[5].asDOUBLE();
        GTime epoch_gps =  GTime::CivilTime2GTime(ct);
        GTime epoch_utc = GTime::GPST2UTC(epoch_gps);
        
        GVector acc_af(splitstr[6].asDOUBLE(),splitstr[7].asDOUBLE(),splitstr[8].asDOUBLE() );
        
        double half_angle = acos(splitstr[9].asDOUBLE());
        GQuaternion q;
        q.set(half_angle, GVector(splitstr[10].asDOUBLE(),splitstr[11].asDOUBLE(),splitstr[12].asDOUBLE()  ));
        
        GVector acc_bf(acc_af.x, acc_af.y, acc_af.z);
        
        GVector acc_eci = q.rotate(acc_bf);
        
        GVector acc_RAC;
        
        GSpaceEnv::updateSpaceEnvironment(epoch_utc);
        
        
        mysatA.getEphValue_test(epoch_gps).getPV(mp_ecef,mv_ecef);
        mysatA.getStatePointer()->updateState_ecef(epoch_utc, mp_ecef, mv_ecef);
        
        mysatB.getEphValue_test(epoch_gps).getPV(mp_ecef,mv_ecef);
        mysatB.getStatePointer()->updateState_ecef(epoch_utc, mp_ecef, mv_ecef);
        
        
        //process graceA data
        mp_eci = mysatA.getStatePointer()->satpos_eci;
        mv_eci = mysatA.getStatePointer()->satvel_eci;
        
        zbfs = -mp_eci;
        zbfs.normalise();
        
        xbfs = (mysatB.getStatePointer()->satpos_eci - mysatA.getStatePointer()->satpos_eci);
        xbfs.normalise();
        
        ybfs = crossproduct(zbfs, xbfs);
        
        double RM[9];
        q.getRotationMatrix(RM);
        
        GVector acc_eci2;
        acc_eci2.x = xbfs.x * acc_bf.x + ybfs.x * acc_bf.y + zbfs.x * acc_bf.z;
        acc_eci2.y = xbfs.y * acc_bf.x + ybfs.y * acc_bf.y + zbfs.y * acc_bf.z;
        acc_eci2.z = xbfs.z * acc_bf.x + ybfs.z * acc_bf.y + zbfs.z * acc_bf.z;
        
        acc_eci = acc_eci2;
        
        GVector R = mp_eci;
        R.normalise();
        
        GVector C = crossproduct(mp_eci, mv_eci);
        C.normalise();
        
        GVector A = crossproduct(R, C);
        
        //acc_RAC.x = dotproduct(R, acc_eci);
        //acc_RAC.y = dotproduct(A, acc_eci);
        //acc_RAC.z = dotproduct(C, acc_eci);
        
        acc_RAC.x = acc_eci.x * R.x + acc_eci.y * R.y + acc_eci.z * R.z;
        acc_RAC.y = acc_eci.x * A.x + acc_eci.y * A.y + acc_eci.z * A.z;
        acc_RAC.z = acc_eci.x * C.x + acc_eci.y * C.y + acc_eci.z * C.z;
        
        printf("%s %.9E %.9E %.9E\n", ct.TimeString().c_str(), acc_RAC.x, acc_RAC.y, acc_RAC.z);
        
    }
    
    infile.close();
    int testc = 0;
}


int main(int argc, char* argv[])
{
   
    double coef[4]={-0.00000012060287539500901,0.83428284555915622,-0.00000012060287539500901,-0.18073270658108179};
    double xx[4]={0.0};
    GMath::solve_quartic(coef[0], coef[1], coef[2], coef[3], xx);
    
    //for galileo(ecllipse season)
    // GSpaceCraftMgr::loadPreciseEphemeris(sp3dir+"com18264.sp3");
    // GSpaceCraftMgr::loadPreciseEphemeris(sp3dir+"com18265.sp3");
    // GSpaceCraftMgr::loadPreciseEphemeris(sp3dir+"com18266.sp3");
    // GSpaceCraftMgr::loadPreciseEphemeris(sp3dir+"com18270.sp3");
    // GSpaceCraftMgr::loadPreciseEphemeris(sp3dir+"com18271.sp3");
    
    
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
    
    //GSpaceCraftMgr::loadUCLEphemeris("sp3/ucl_svn46_01mar04.eci");
    
    // testing data
    
    GString sp3dir = "/Users/lizhen/experiments/data/sp3/2015codesp3/";
    
    argv[1] = "/Users/lizhen/experiments/data/gfcsetup_test.cfg";
    argv[2] = "ssGAL";  //satellite system; ssGRACE
    argv[3] = "11";     //prn 1
    
    //2015/01/11/18:34:25.000000
//    argv[4] = "2015 1 11 17 34 30.00000000"; // start time in GPST
//    argv[5] = "2015 1 12 17 34 30.00000000"; // end time in GPST
    
    
//    argv[4] = "2016 01 7 15 28 22.00000000"; // start time in GPST
//    argv[5] = "2016 01 7 15 28 40.00000000"; // end time in GPST
//    argv[6] ="com18775.sp3,com18776.sp3,com18780.sp3,com18781.sp3,com18782.sp3,com18783.sp3,com18784.sp3,com18785.sp3,com18786.sp3,com18790.sp3,com18791.sp3,com18792.sp3,com18793.sp3,com18794.sp3,com18795.sp3";
    
    
    
    //for galileo
    argv[4] = "2015 1 11 00 00 00.00000000"; // start time in GPST
    argv[5] = "2015 1 12 00 00 00.00000000"; // end time in GPST
    argv[6] = "com18264.sp3,com18265.sp3,com18266.sp3,com18270.sp3,com18271.sp3,com18272.sp3,com18273.sp3,com18274.sp3,com18275.sp3,com18276.sp3,com18280.sp3,com18281.sp3,com18282.sp3,com18283.sp3";
    
    //for bds
    //argv[4] = "2015 01 01 00 00 00.00000000"; // start time in GPST
    //argv[5] = "2015 12 31 00 00 28.00000000"; // end time in GPST
    
    //for GRACE
    //argv[4] = "2008 03 25 01 47 40.000000"; // start time in GPST
    //argv[5] = "2008 03 25 01 48 30.000000"; // end time in GPST
    
   // argv[4] = "2008 03 25 01 00 00.000000"; // start time in GPST
   // argv[5] = "2008 03 25 23 00 00.000000"; // end time in GPST
    
//    argv[4] = "2007 01 20 00 01 00.000000"; // start time in GPST
//    argv[5] = "2007 01 20 23 59 00.000000"; // end time in GPST
    
    // gpst: 2007/01/20/00:40:28.000000 --> 2007/01/20/00:40:56.000000
    //argv[4] = "2007 01 20 00 29 24.000000"; // start time in GPST
    //argv[5] = "2007 01 20 00 29 52.000000"; // end time in GPST
    
    // one penumbra event grace
    //argv[4] = "2007 01 20 01 27 20.000000"; // start time in GPST
    //argv[5] = "2007 01 20 01 28 00.000000"; // end time in GPST
    //argv[4] = "2007 01 20 02 03 18.000000"; // start time in GPST
    //argv[5] = "2007 01 20 02 04 05.000000"; // end time in GPST
    
    // two penumbra event grace
    //argv[4] = "2007 01 20 06 09 20.000000"; // start time in GPST
    //argv[5] = "2007 01 20 06 10 10.000000"; // end time in GPST
    //argv[4] = "2007 01 20 06 45 15.000000"; // start time in GPST
    //argv[5] = "2007 01 20 06 45 55.000000"; // end time in GPST
    
    
    //eclipse orbit prediction test
    //argv[4] = "2015 1 9 23 18 30.00000000"; // start time in GPST
    //argv[5] = "2015 1 9 23 20 09.00000000"; // start time in GPST
    
    //argv[4] = "2015 1 10 00 16 40.00000000"; // start time in GPST
    //argv[5] = "2015 1 10 00 18 20.00000000"; // end time in GPST
    
    //argv[4] = "2015 1 10 14 22 10.00000000"; // start time in GPST
    //argv[5] = "2015 1 10 14 23 50.00000000"; // end time in GPST
    
    
    //argv[4] = "2015 1 11 17 34 30.00000000"; // start time in GPST
    //argv[5] = "2015 1 11 17 36 5.00000000"; // end time in GPST
    
    //argv[4] = "2015 1 11 18 33 5.00000000"; // start time in GPST
    //argv[5] = "2015 1 11 18 34 40.00000000"; // end time in GPST
    
    
    
   
    
    
    
//    argv[4] = "2016 01 01 01 00 00.000000"; // start time in GPST
//    argv[5] = "2016 01 01 23 00 30.000000"; // end time in GPST
    //argv[6]="com18264.sp3,com18265.sp3,com18266.sp3,com18270.sp3,com18271.sp3";
    
    //argv[6]="sp3/whusp3-2015/wum18253.sp3,sp3/whusp3-2015/wum18254.sp3,sp3/whusp3-2015/wum18255.sp3,sp3/whusp3-2015/wum18256.sp3,sp3/whusp3-2015/wum18260.sp3";
    
    
    argv[7] = "1";  // interval
    
    
    //argv[4] = "2015 1 9 09 13 35.00000000"; // start time in GPST
    //argv[5] = "2015 1 9 10 12 32.00000000"; // end time in GPST
    
    
    
    /*
    argv[6] = "sp3/codsp3-2015/com18342.sp3,sp3/codsp3-2015/com18343.sp3,sp3/codsp3-2015/com18344.sp3,sp3/codsp3-2015/com18345.sp3,sp3/codsp3-2015/com18346.sp3,sp3/codsp3-2015/com18350.sp3,sp3/codsp3-2015/com18351.sp3,sp3/codsp3-2015/com18352.sp3,sp3/codsp3-2015/com18353.sp3,sp3/codsp3-2015/com18354.sp3,sp3/codsp3-2015/com18355.sp3";
    */
    
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
    
    
    //GSpaceCraftMgr::loadGRACEEphemeris("/Users/lizhen/experiments/data/sp3/GRACE/nav_0120.data");
    
    for(int i = 0 ; i< sp3file.size(); i++ )
    {
        GSpaceCraftMgr::loadPreciseEphemeris(sp3dir + sp3file[i]);
    }
    
    GSensorID myid(satsys, satprn); // C08  C06 and C14 , G11, E11(IOV101)
    GSpaceCraft& mysat = GSpaceCraftMgr::gSpacecraft[myid.getIDString()];
    
    //GSensorID myidB("ssGRACE", 2); // C08  C06 and C14 , G11, E11(IOV101)
    //GSpaceCraft& mysatB = GSpaceCraftMgr::gSpacecraft[myidB.getIDString()];
    
    
    
    //loadGRACE_acc_sca("/Users/lizhen/experiments/data/sp3/GRACE/graceACC_SCA_0125.data", mysat,mysatB);
    
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
//    GKeplerianElements mykpe(42158.5455173419890933,0.00297678789158783898422,0.946320613861091786134,3.54313401237741487648,5.60117667874381430781,0.521560853497004761614);
    
    //C11 at GPST 2015 01 01 00 00 00
    //GKeplerianElements mykpe(27904.7761430916347916,0.00240232320729961978945,0.973095890700745203005,3.46758595575894412733,1.46975968471844375962,0.80975252861510960134);
    
    //C12
    //GKeplerianElements mykpe(27904.4660136779959139,0.00274336440899936029325,0.971964838760725591271,3.36304948820225702022,1.46073050121603542628,1.7173513865787128907);
    
    //C14
    //GKeplerianElements mykpe(27905.3627291815599456,0.00156609722616517493127,0.956695310707969204244,3.56007288522406373232,3.54991651616790226953,4.91337390378707183913);
    GKeplerianElements
mykpe(6838.48265785065597999,0.00303136072481130849909,1.55413552698389439467,1.52503752498879021184,1.10049729087647997993,3.91570627917123768803);
    
    double shadow_factor = 0.0;
    
    while( 1 )
    {
        GTime epoch_utc = GTime::GPST2UTC(epoch_gps);
        //GTime epoch_utc = epoch_end_fit;
        
        JDTime jdt = GTime::GTime2JDTime(epoch_utc);
        CivilTime ct = GTime::JDTime2CivilTime(jdt);
        
        GVector mp_ecef, mv_ecef, mp_eci, mv_eci;
        //orbit comparison
       
        GSpaceEnv::updateSpaceEnvironment(epoch_utc);
        
        mysat.getEphValue_test(epoch_gps).getPV(mp_ecef,mv_ecef);
        mysat.getStatePointer()->updateState_ecef(epoch_utc, mp_ecef, mv_ecef);
        
        //get the shadow function from motion state
       shadow_factor = mysat.getStatePointer()->shadow_factor;
//
//        printf("%.9E %.9E %.9E %.9E %.9E %.9E\n", double( mysat.getStatePointer()->keplerianElement.m_sma),
//                                            double (mysat.getStatePointer()->keplerianElement.m_ecc),
//                                            double (mysat.getStatePointer()->keplerianElement.m_inc),
//                                            double (mysat.getStatePointer()->keplerianElement.m_raan),
//                                            double(mysat.getStatePointer()->keplerianElement.m_argp),
//                                            double(mysat.getStatePointer()->keplerianElement.m_tran)  );
        
        //mykpe.propagate(mykpe, (epoch_gps-start_gps).toSeconds(), mp_eci, mv_eci);
        
        // for the eclipse event detection
        mysat.getStatePointer()->shadow_detector(epoch_gps,GSpaceEnv::planetPos_ecef[GJPLEPH::SUN], mp_ecef);
        
        //for the eci cooridnates
        //mysat.getStatePointer()->shadow_detector(epoch_gps,GSpaceEnv::planetPos_eci[GJPLEPH::SUN], mysat.getStatePointer()->satpos_eci);
        
        
        //double shadow_factor = GMotionState::shadowFactor(6371.0,6371.0,GSpaceEnv::planetPos_eci[GJPLEPH::SUN], mysat.getStatePointer()->satpos_eci);
        
        //double shadow_factor = GMotionState::shadowFactor_SECM(true,GSpaceEnv::planetPos_eci[GJPLEPH::SUN] , mysat.getStatePointer()->satpos_eci);
        
       
        //printf("%s ",GTime::GTime2CivilTime(epoch_gps).TimeString().c_str());
        //double shadow_factor = GMotionState::myshadowFactor(GSpaceEnv::planetPos_eci[GJPLEPH::SUN] ,mysat.getStatePointer()->satpos_eci);
        
        //double shadow_factor = mysat.getStatePointer()->attitude_eci.eclipse(mysat.getStatePointer()->satpos_eci,GSpaceEnv::planetPos_eci[GJPLEPH::SUN]);
        
        // output the shadow function
        //double shadow_factor = GMotionState::shadowFactor_SECM(true,GSpaceEnv::planetPos_ecef[GJPLEPH::SUN] , mp_ecef);
        
//
       // printf("%s %f %f\n",GTime::GTime2CivilTime(epoch_gps).TimeString().c_str(), shadow_factor,mysat.getStatePointer()->attitude_eci.eta*180/3.1415926 );
        
        
        //mysat.getStatePointer()->updateState_ecef(epoch_utc, mp_ecef, mv_ecef);
        
        if( ( epoch_gps = epoch_gps + interval) >= end_gps)
        {
            break;
        }
        
        //obscount++;
    }

    
    
    return 0;
}



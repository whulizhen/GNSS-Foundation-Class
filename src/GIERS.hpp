//
//  GIERS.hpp
//  GFC
//
//  Created by lizhen on 12/08/2016.
//  Copyright © 2016 lizhen. All rights reserved.
//

#ifndef GIERS_hpp
#define GIERS_hpp

#include <stdio.h>
#include <math.h>
#include <fstream>
#include <iostream>

#include "GFCCONST.h"
using namespace std;
namespace gfc
{
    
    // this class contains basic functions in IERS
    class GIERS
    {
        
    public:
        
        static void  FUNDARG( double T, double& L, double& LP, double& F, double& D, double& OM);
        static void  DOODSONARG(double JD_UT1_I, double JD_UT1_F, double JD_TT,double FNUT[5], double BETA[6]);
        
        static void  CNMTX( double MJD_TT, double* H );
        static void  ORTHO_EOP( double MJD_TT,double* EOP_correction);
        static void  PMSDNUT2( double RMJD, double* pm );
        static void UTLIBR(double RMJD, double& dut1, double& dlod);
        
        static void  RG_ZONT2(double JC_TT, double& DUT, double& DLOD, double& DOMEGA);
        
        static double  getGMST( double JD_UT1_I,double JD_UT1_F, double JD_TT);
        static double  getGAST( double JD_UT1_I, double JD_UT1_F ,double JD_TT, double dpsi );
        
        static double  getEarthRotationAngle(double JD_UT1_I, double JD_UT1_F);
        static  double getEarthRotationAngleRate(double JC_TT);
        static double  getMeanObliquity(double JC_TT);
        static double  getEECT(double JC_TT);
        static double  getEE2000( double JC_TT, double dpsi);
        
        // for solid earth tide displancement corrections
        static int  CAL2JD(int IY, int IM, int ID,double& DJM0, double& DJM);
        
        static void ST1IDIU(double* XSTA, double* XSUN, double* XMON, double FAC2SUN,double FAC2MON, double* XCORSTA);
        static void ST1L1( double* XSTA, double* XSUN, double* XMON, double FAC2SUN,double FAC2MON, double* XCORSTA);
        static void ST1ISEM ( double* XSTA, double* XSUN, double* XMON, double FAC2SUN,double FAC2MON, double* XCORSTA);
        static void STEP2DIU(double* XSTA, double FHR, double T, double* XCORSTA);
        static void STEP2LON(double* XSTA, double T, double* XCORSTA);
        //static void DEHANTTIDEINEL(double *XSTA, double JD_UTC,double leapsec, double *XSUN, double *XMON, double *DXTIDE);
        static void DEHANTTIDEINEL(double *XSTA, int YR,int MONTH, int DAY, double FHR, double leadsec, double *XSUN, double *XMON, double *DXTIDE);
        
        //function in chapter 9 FCULZD_HPA
        static void FCULZD_HPA(double latitude, double ellip_ht, double pressure, double wvp, double lambda_um, double& fcul_ztd, double& fcul_zhd,double& fcul_zwd);
        static double FCUL_A(double LATITUDE, double HEIGHT_M, double T_K, double ELEV_DEG);
        
        static void ARG2(int IYEAR, double DAY, double ANGLE[11] );
        
        static int MDAY(int iy, int im);
        
        static bool LEAP(int iy);
        
        static void TOYMD(int IT1[2], int IT2[3]);
        
    private:
        
        static double normalizeAngle(double a);
        
    };
    
    
    
    
    
    
    
    
} // end of namespace gfc









#endif /* GIERS_hpp */

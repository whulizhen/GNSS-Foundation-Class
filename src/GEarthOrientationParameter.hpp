//
//  GEarthOrientationParameter.hpp
//  GFC
//
//  Created by lizhen on 16/4/4.
//  Copyright © 2016年 lizhen. All rights reserved.
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

#ifndef GEarthOrientationParameter_hpp
#define GEarthOrientationParameter_hpp

#include <stdio.h>
#include <iostream>
#include <fstream>
#include "GException.h"
#include "GTime.h"
#include "GMath.hpp"
#include "GVector.hpp"
#include "GIERS.hpp"
#include "GEllipsoidMgr.h"
#include "GIAU.hpp"


using namespace std;
namespace gfc
{
    class GEarthOrientationParameter
    {
        
    public:
        struct EOP
        {
          public:
            double  m_polarMotionX; // polar motion x, Unit: arcsec
            double  m_polarMotionXSigma; // the error of polarmotionX
            double  m_polarMotionY; // polar motion y, Unit: arcsec
            double  m_polarMotionYSigma; // the error of polarmotionX
            double  m_UT1mUTC;  // UT1 minus UTC, Unit: second
            double  m_UT1mUTCSigma; // the error of UT1-UTC
            double  m_LOD;  // the length of day,only in BulletinA ,Unit: second
            double  m_LODSigma;
            double  m_dX;  // Unit: arcsec
            double  m_dXSigma;
            double  m_dY;  // Unit: arcsec
            double  m_dYSigma;
            double  m_dPsi;
            double  m_dPsiSigma;
            double  m_dEpsilon;
            double  m_dEpsilonSigma;
            
            double    m_dOMEGA;    // the correction to the earth rotation speed, rad/second
            
            GString  m_type; // "final" of "prediction"
            EOP()
            {
                double maxdev = 999.999;
                m_polarMotionX = 0.0;m_polarMotionXSigma = maxdev;m_polarMotionY = 0.0;m_polarMotionYSigma = maxdev;
                m_UT1mUTC =0.0; m_UT1mUTCSigma = maxdev; m_LOD =0.0; m_LODSigma = maxdev;
                m_dX =0.0; m_dXSigma= maxdev; m_dY =0.0; m_dYSigma= maxdev; m_dPsi=0.0; m_dPsiSigma = maxdev;
                m_dEpsilon =0.0; m_dEpsilonSigma= maxdev;
                
                m_dOMEGA = 0.0;
            }
        };
        
        
    public:
        
        static std::map< double, EOP> loadEOP();
        static void loadEOP(GString eopfile);
        
        GEarthOrientationParameter( GTime epochTime);
        GEarthOrientationParameter() {};
        ~GEarthOrientationParameter();
        
        void setEpochTime(GTime epochTime);
        
        double getUT1mUTC();
        double getPMX();
        double getPMY();
        
        //reference: IERS2010,technNote36, section 7.1.4
        static void getPM_mean(double mjdUTC, double& pmx, double& pmy);
        
        //get the pole tide impacts for earth gravity
        void getPoleTide( double& dC21, double& dS21);
        //get the pole tide displacement for the station position
        void getPoleTide( double yearSinceJ2k, gfc::GVector pos, gfc::GVector &displacement );
        
        
        
        
        
        void getFrameBias00(double& dpsibi,double& depsbi,double& dra);
        
        void getPrecessionRate00(double JC_TT, double& dpsipr, double& depspr);
        void getPrecessionAngle( double JC_TT, double &deltaA, double &zA, double &thetaA);
        
        /*IAU2006 model*/
        void getPrecessionAngle(double JC_TT, double & gamb, double & phib, double &psib, double& epsa);
        
        void getPrecessionMatrix( double JC_TT,double* pr);
        
        void iauNutation1980(double JC_TT,double& dpsi,double& deps);
        void iauNutation2000A(double JC_TT, double& dpsi,double& deps);
        void iauNutation2000B(double JC_TT,double& dpsi, double& deps);
        void iauNutation2006A(double JC_TT, double& dpsi, double& deps);
        //the main function to get the nutation angels
        void getNutationAngle( double JC_TT,double& dpsi,double& deps);
        
        void  getNutationMatrix( double JC_TT, double* nr);
        void  getPolarMotionMatrix( double JC_TT,double* pmr);
        
        void  getRotationBPN(double JC_TT,double* BPN);
        double getLocatorS06(double JC_TT,double x, double y);
        
        void computeRotationMatrix(double* ECI2ECEFPos, double* ECI2ECEFVel);
        
        void computeRotationMatrix1(double* ECI2ECEFPos, double* ECI2ECEFVel);
        
        double getSP2000(double JC_TT);
        
        void   ECI2ECEF(int tag, double* eci,double* ecef);
        void   ECEF2ECI(int tag, double* ecef,double* eci);
        
        void   ECEF2ECI_pos( GVector& pos_ecef, GVector& pos_eci );
        void   ECEF2ECI_vel( GVector& pos_ecef, GVector& vel_ecef, GVector& vel_eci );
        void   ECI2ECEF_pos( GVector& pos_eci, GVector& pos_ecef );
        void   ECI2ECEF_vel(GVector& pos_eci, GVector& vel_eci, GVector& vel_ecef );
        void   getECI2ECEFMatrix(double* tm);
        
    private:
        // get the eop for current time with interplation(Lagrange method)
        EOP getEOP();
        double normalizeAngle( double a );
        void matrixMultiply(int rA,int cA,double* A, int cB,double* B,double* C);
        void matrixTranspose(int r, int c, double* m);
        void Rx( double phi,double* R);
        void Ry( double phi,double* R);
        void Rz( double phi,double* R);
        
        // only be used inside this class
        static std::map< double, EOP > eopStorage;  // the eop storage parameter,the key is mjd in UTC
        
        EOP       m_eop;     // the eop for current time
        
        GTime    m_epoch;  // epoch time, usually should be in UTC time system
        double    m_eci2ecefPos[9]; // the transformation for position
        double    m_eci2ecefVel[9]; // the transformation for velocity
        
    };
    
    
    
}// end of namespace gfc


#endif /* GEarthOrientationParameter_hpp */

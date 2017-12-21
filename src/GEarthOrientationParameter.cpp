//
//  GEarthOrientationParameter.cpp
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

#include "GEarthOrientationParameter.hpp"
namespace gfc
{
    
    std::map<double,GEarthOrientationParameter::EOP> GEarthOrientationParameter::eopStorage  = GEarthOrientationParameter::loadEOP();
    
    //just for the definition of the static variable.
    std::map< double, GEarthOrientationParameter::EOP> GEarthOrientationParameter::loadEOP()
    {
        std::map<double,GEarthOrientationParameter::EOP> myeop;
        return myeop;
    }
    
    //override static function
    // read C04 earth orientation parameter, px py ut1-utc lod p
    // file source: ftp://hpiers.obspm.fr/iers/eop/eopc04/eopc04.dPsi_dEps.62-now
    // FORMAT(3(I4),I7,2(F11.6),2(F12.7),2(F11.6),2(F11.6),2(F11.7),2F12.6)
    void GEarthOrientationParameter::loadEOP(GString eopfile)
    {
        //std::map<double,GEarthOrientationParameter::EOP> eopMap;
        GTime J2000 = GTime::J2000();
        
        EOP eop;
        //start reading eop file
        ifstream inpf(eopfile.c_str());
        if(!inpf)
        {
            FileMissingException fme("Could not open C04 eop file " + eopfile);
            GFC_THROW(fme);
        }
        //first clear all the data
        eopStorage.clear();
        bool ok (true);
        GString line;
        
        // skip the first 14 lines since they are the header for description
        for(int i = 0 ; i< 14; i++ )
        {
            getline(inpf,line);
        }
        
        while(!inpf.eof() && inpf.good())
        {
            
            getline(inpf,line);
            line.stripTrailing_v('\r');
            
            if(inpf.eof()) break;
            
            // line length is actually 156
            if( inpf.bad() || line.size() < 70) { ok = false; break; }
            
            double mjd = line.substr(12,7).asDOUBLE();
            
            eop.m_polarMotionX   = line.substr(19,11).asDOUBLE();  // arcseconds
            eop.m_polarMotionY  = line.substr(30,11).asDOUBLE();  //arcseconds
            
            eop.m_UT1mUTC = line.substr(41,12).asDOUBLE();// seconds
            eop.m_LOD      = line.substr(53,12).asDOUBLE();  // seconds
            
            //eop.m_dPsi     = line.substr(65,11).asDOUBLE(); // arcsec
            //eop.m_dEpsilon  = line.substr(76,11).asDOUBLE(); //arcsec
            eop.m_dX     = line.substr(65,11).asDOUBLE(); // arcsec
            eop.m_dY  = line.substr(76,11).asDOUBLE(); //arcsec
            
            
            eop.m_polarMotionXSigma = line.substr(87,11).asDOUBLE();  // err in xp
            eop.m_polarMotionYSigma = line.substr(98,11).asDOUBLE();  // err in yp
            eop.m_UT1mUTCSigma = line.substr(109,11).asDOUBLE();  // err in yp
            eop.m_LODSigma = line.substr(120,11).asDOUBLE();  // err in yp
            
            //eop.m_dPsiSigma = line.substr(131,12).asDOUBLE();
            //eop.m_dEpsilonSigma = line.substr(143,12).asDOUBLE();
            
            eop.m_dXSigma = line.substr(131,12).asDOUBLE();
            eop.m_dYSigma = line.substr(143,12).asDOUBLE();
            
            eop.m_type = "final";
            eopStorage[mjd] = eop;
            //addEOPData(MJD(mjd,TimeSystem::UTC), EOPData(xp,yp,UT1mUTC,dPsi,dEps));
            
        };
        
        inpf.close();
        
        if( !ok )
        {
            FileMissingException fme("IERS File " + eopfile + " is corrupted or wrong format");
            GFC_THROW(fme);
        }
        
        //eopStorage = eopMap;
    }
    
    /*epochTime should in UTC time*/
    GEarthOrientationParameter::GEarthOrientationParameter(GTime epochTime)
    {
        setEpochTime(epochTime);
    }
    
    void GEarthOrientationParameter::setEpochTime(GTime epochTime)
    {
        static double  RMJD0   = 51544.5E0;
        
        TimeSystem ts =  epochTime.getTimeSystem();
        if( ts.getTimeSystemName() != "tsUTC")
        {
            printf("WARNING: the time system is not UTC in GEarthOritentationParameter class!\n");
        }
        
        m_epoch = epochTime;
        //interplate and set up the earth orientation parameter for the current time
        m_eop = getEOP();
       /*
        The long-period terms, as well as the secular variation of the libration con- tribution, are already contained in the observed polar motion and need not be added to the reported values (x,y)IERS
        */
        
        double dut =0.0, dlod =0.0, domega = 0.0;
        
        //GTime TAI =  GTime::UTC2TAI(m_epoch);
        //GTime TT  =  GTime::TAI2TT(TAI);
        
        //some corrections for the current interplated eop
        double eop_correction[3]; //dx, dy, dut1
        double pm_correction[2];
        
        double dut1=0.0;
        dlod =0.0;
        
        // The libration part of ut1
        GIERS::UTLIBR(m_epoch.toDays(), dut1, dlod);
        
        // the corrections of ocean tide of Ray model
        GIERS::ORTHO_EOP(m_epoch.toDays(), eop_correction);
        
        // the libration part
        GIERS::PMSDNUT2(m_epoch.toDays(),pm_correction);
       
        
        eop_correction[0] += pm_correction[0];
        eop_correction[1] += pm_correction[1];
        eop_correction[2] += dut1;
        
        m_eop.m_polarMotionX += (eop_correction[0]*1.0E-6);
        m_eop.m_polarMotionY += (eop_correction[1]*1.0E-6);
        m_eop.m_UT1mUTC += (eop_correction[2]*1.0E-6);
        m_eop.m_LOD += dlod*1.0E-6;
        
        // the zonal earth tides on the rotation of earth
        // reference: IERS2010, RG_ZONT2.F
        
        // how to use RZ_ZONT2 corrections ?????
        // only when eop does NOT include tide effects, add this correction
        /*
        GIERS::RG_ZONT2( (TT.toDays() - RMJD0)/36525.0, dut, dlod, domega);
        m_eop.m_UT1mUTC += dut;
        m_eop.m_LOD += dlod;
        m_eop.m_dOMEGA += domega;
        */
        
        // get the transformation matrix for the current time
        //computeRotationMatrix1(m_eci2ecefPos, m_eci2ecefVel);
        
        computeRotationMatrix(m_eci2ecefPos, m_eci2ecefVel);
        
        int testc = 0.0;
    }
    
    
    
    void GEarthOrientationParameter::matrixTranspose(int r, int c, double* m)
    {
        double *t = new double[r*c*sizeof(double)];
        memcpy(t,m,sizeof(double)*r*c);
        
        for( int i = 0 ; i< r; i++ )
        {
            for( int j = 0 ; j<c; j++ )
            {
                m[j*r+i] = t[i*c+j];
            }
        }
        
        if(t != NULL) {delete[] t, t= NULL;}
    }
    
    /*A*B=C*/
    void GEarthOrientationParameter::matrixMultiply(int rA, int cA, double *A, int cB, double *B, double *C)
    {
        //rA*cA
        //rB*cB
        for( int i = 0 ; i< rA; i++ )
        {
            for(int j = 0 ; j< cB; j++ )
            {
                C[i*cB+j] = 0.0;
                for( int k = 0 ; k< cA; k++ )
                {
                    C[i*cB+j] += A[i*cA+k]*B[k*cB+j];
                }
            }
        }
        
    } // end of the function
    
    
    GEarthOrientationParameter::~GEarthOrientationParameter()
    {
        
    }
    
    
    // get the eop for current time
    GEarthOrientationParameter::EOP GEarthOrientationParameter::getEOP()
    {
        EOP myeop;
        //check whether the current time is between the
        //NOTE the Time System of m_epoch
        TimeSystem ts;
        long mjd,sod;
        double fsod;
        double secpday = GCONST("SECPDAY");
        m_epoch.GetData( ts, mjd, sod, fsod);  // please note, m_epoch is in UTC , the time in C04 eop file is also UTC
        
        double currentTime = mjd + (sod+fsod)/secpday;
        std::map<double,EOP>::iterator it = eopStorage.begin();
        double initialTime = it->first;
        it  = eopStorage.end() ;it--;
        
        double finalTime = it->first;
        
        if( currentTime<initialTime || currentTime > finalTime)
        {
            InvalidRequest ire("GEarthOrientationParameter: currentTime is invalid,out of range\n");
            GFC_THROW(ire);
        }
        
        it = eopStorage.find(currentTime);
        if( it != eopStorage.end() )
        {
            myeop = it->second;
            return myeop;
        }
        // 考虑采用largrange 插值
        // 在ops中，采用的是采用最小二乘法每3天拟合一次系数。
        int interPoints = 5;// according to the interp.f , a window of 4 data points
        const int half = ( interPoints + 1 ) / 2;
        it = eopStorage.lower_bound(currentTime); //返回大于等于currentTime 的值
        if(currentTime > finalTime)
        {
            it = eopStorage.end();
            it--;
        }
        
        std::map<double,EOP>::const_iterator its =it;
        std::map<double,EOP>::const_iterator ite =it;
        
        if(int(eopStorage.size())> 2*half)
        {
            int ileft = half;
            for( int i = 0; i < half; i++)
            {
                if(its==eopStorage.begin()) break;
                its--;
                ileft--;
            }
            
            int iright = half-1+ileft;
            for( int i = 0; i < (half-1+ileft); i++)
            {
                ite++;
                if( ite == eopStorage.end())
                {
                    ite--;
                    break;
                }
                iright--;
            }
            
            int ileft2 = iright;
            for(int i = 0; i < iright; i++)
            {
                if(its == eopStorage.begin()) break;
                its--;
                ileft2--;
            }
            
            if( ileft2 > 0)
            {
                // the code never go here
                // just throw an exception
                InvalidRequest e("My God, it should never go here!!!");
                GFC_THROW(e);
            }
        }
        else
        {
            its = eopStorage.begin();
            ite = eopStorage.end();
            ite--;
        }
        
        const size_t N = 16; // total 16 variables need to be interplated
        // note the structure of EOP
        std::vector<double> times;
        times.reserve(10);
        std::vector<std::vector<double> > datas;;
        datas.resize(N);
        
        std::map<double,EOP>::const_iterator itrEnd = ite;
        itrEnd++;
        for( std::map<double,EOP>::const_iterator itr=its; itr!=itrEnd; itr++)
        {
           // double time = itr->first;
            EOP vd = itr->second;
            times.push_back( itr->first - its->first );
            // note the relationship of order and the variable
            datas[0].push_back(vd.m_dEpsilon);
            datas[1].push_back(vd.m_dEpsilonSigma);
            datas[2].push_back(vd.m_dPsi);
            datas[3].push_back(vd.m_dPsiSigma);
            datas[4].push_back(vd.m_dX);
            datas[5].push_back(vd.m_dXSigma);
            datas[6].push_back(vd.m_dY);
            datas[7].push_back(vd.m_dYSigma);
            datas[8].push_back(vd.m_LOD);
            datas[9].push_back(vd.m_LODSigma);
            datas[10].push_back(vd.m_polarMotionX);
            datas[11].push_back(vd.m_polarMotionXSigma);
            datas[12].push_back(vd.m_polarMotionY);
            datas[13].push_back(vd.m_polarMotionYSigma);
            datas[14].push_back(vd.m_UT1mUTC);
            datas[15].push_back(vd.m_UT1mUTCSigma);
        }
        
        std::vector<double> dd(N, 0.0);
        double dt = currentTime - its->first;
        
        for( size_t i = 0; i < N; i++ ) //interplate all the variables(16 of them)
        {
            dd[i] = GMath::SimpleLagrangeInterpolation(times,datas[i],dt);
        }
        
        //return all the interplated data;
        myeop.m_type ="interplation";
        myeop.m_dEpsilon = dd[0]; myeop.m_dEpsilonSigma = dd[1];
        myeop.m_dPsi = dd[2]; myeop.m_dPsiSigma = dd[3];
        myeop.m_dX = dd[4]; myeop.m_dXSigma = dd[5];
        myeop.m_dY = dd[6]; myeop.m_dYSigma = dd[7];
        myeop.m_LOD = dd[8]; myeop.m_LODSigma = dd[9];
        myeop.m_polarMotionX = dd[10]; myeop.m_polarMotionXSigma = dd[11];
        myeop.m_polarMotionY = dd[12]; myeop.m_polarMotionYSigma = dd[13];
        myeop.m_UT1mUTC = dd[14]; myeop.m_UT1mUTCSigma = dd[15];
        
        return myeop;
    }
    
    
    /*
     *
     *  return the ut1mutc for current epoch
     *
     */
    double GEarthOrientationParameter::getUT1mUTC()
    {
        return m_eop.m_UT1mUTC;
    }
    
    double GEarthOrientationParameter::getPMX()
    {
        return m_eop.m_polarMotionX;
    }
    
    double GEarthOrientationParameter::getPMY()
    {
        return m_eop.m_polarMotionY;
    }
    
    
    //return value: pmx and pmy in arcsec
    void GEarthOrientationParameter::getPM_mean(double mjdUTC, double& pmx, double& pmy)
    {
        static const double MJD_J2000 = 51544.5;
        static const double MJD_2010 = 55197; //the mjd of 2010.0
        
        double time_diff = (mjdUTC - MJD_J2000)/365.25;
        double time_coef[4] = { 1.0, time_diff, time_diff*time_diff, time_diff*time_diff*time_diff };
        pmx =0.0;
        pmy = 0.0;
        // before 2010.0, xp_m = 55.974,1.8243,0.18413,0.007024; yp_m = 346.346,1.7896,-0.10729,-0.000908
        // after 2010.0,  xp_m = 23.513,7.6141,0.0,0.0; yp_m = 358.891,-0.6287, 0.0, 0.0
        
        double xp1[4] = {55.974,1.8243,0.18413,0.007024};
        double yp1[4] = {346.346,1.7896,-0.10729,-0.000908};
        double xp2[4] = {23.513,7.6141,0.0,0.0};
        double yp2[4] = {358.891,-0.6287, 0.0, 0.0};
        
        double xp[4], yp[4];
        if(mjdUTC < MJD_2010)
        {
            memcpy(xp, xp1,sizeof(double)*4);
            memcpy(yp, yp1,sizeof(double)*4);
        }
        else
        {
            memcpy(xp, xp2,sizeof(double)*4);
            memcpy(yp, yp2,sizeof(double)*4);
        }
        
        for( int i = 0 ; i< 4; i++ )
        {
            pmx += time_coef[i]*xp[i];
            pmy += time_coef[i]*yp[i];
        }
        
        //convert from mas to arc sec
        pmx = pmx *1.0E-3;
        pmy = pmy *1.0E-3;
        
    }
    
    /* Solid pole tide to normalized earth potential coefficients
     *
     * @param mjdUtc   UTC in MJD
     * @param dC21     correction to normalized coefficients dC21
     * @param dS21     correction to normalized coefficients dS21
     */
    void GEarthOrientationParameter::getPoleTide( double& dC21, double& dS21)
    {
        // See IERS Conventions 2010 section 7.1.4, P84
        double xpm =0.0, ypm = 0.0;
        GEarthOrientationParameter::getPM_mean( m_epoch.toDays(), xpm, ypm);
        
        double xp = getPMX() ;  // IERS::xPole(mjdUtc);   // in arcsec
        double yp = getPMY() ;  // IERS::yPole(mjdUtc);    // in arcsec
        
        double m1 =  xp - xpm;
        double m2 = -yp + ypm;
        
        
        // Correction to normalized earth potential coefficients
        // C21 and S21 see IERS2010 section 6.4
        
        // http://iers-conventions.obspm.fr/2010/2010_official/tn36.pdf
        // this is for the solid pole tide
        dC21 += -1.333e-9 * ( m1 + 0.0115 * m2 );
        dS21 += -1.333e-9 * ( m2 - 0.0115 * m1 );
        
        
        //SEE IERS conventions 2010 section 6.5 P94
        // OCEAN POLE TIDE
        // Correction to nomalized earth potential coefficients
        dC21 += -2.1778E-10*(m1 - 0.01724*m2);
        dS21 += -1.7232E-10*(m2 - 0.03365*m1);
        
    }
    
    
    /*
    // get the station position displacement due to pole tide
    //  Displacement vector, WGS84 ECEF XYZ meters.
    // pos and displacement should be both in ECEF
     
     *  Maximum displacements because of this effect are:
     *
     *  \li Vertical:    2.5 cm
     *  \li Horizontal:  0.7 cm
     *
     *  For additional information you may consult: Wahr, J.M., 1985,
     *  "Deformation Induced by Polar Motion", Journal of Geophysical
     *  Research, Vol. 90, No B11, p. 9363-9368.
     
    */
    void GEarthOrientationParameter::getPoleTide(double yearSinceJ2k, gfc::GVector pos, gfc::GVector &displacement)
    {
        
        double lat, lon, theta;
        double PI = GCONST("PI");
        double disp[3] = {0.0};
        
        GEllipsoidMgr::GetEllipsoid("WGS84").Geocentric(pos, lat, lon);
        
        double timedif = yearSinceJ2k;
        double xpbar = 0.054 + timedif*0.00083;
        double ypbar = 0.357 + timedif*0.00395;
        
        double m1 = m_eop.m_polarMotionX - xpbar;
        double m2 = ypbar - m_eop.m_polarMotionY ;
        
        
        double sin2lat = sin(2.0*lat);
        double cos2lat = cos(2.0*lat);
        double sinlat = sin(lat);
        double sinlon = sin(lon);
        double coslon = cos(lon);
        
        theta = PI/2.0 - lat;

        // Finally, get the pole tide values, in ENU reference
        // frame and meters
        
        disp[0] =  +0.009 * sinlat  * ( m1*sinlon - m2*coslon ); // E
        disp[1] =  -0.009 * cos2lat * ( m1*coslon + m2*sinlon ); // N
        disp[2] =  -0.033 * sin2lat * ( m1*coslon + m2*sinlon ); // U
        
        displacement.x = disp[0];
        displacement.y = disp[1];
        displacement.z = disp[2];
        
//        // NEU components - IERS(1996) pg 67, eqn. 22 (in eqn 22, r==Up, theta=S, lambda=E)
//        disp[0] = -0.009 * cos(2.0*theta) * (m_eop.m_polarMotionX  * coslon - m_eop.m_polarMotionY * sinlon);  // -S = N
//        disp[1] = -0.009 * cos(theta) * (m_eop.m_polarMotionX * sinlon + m_eop.m_polarMotionY * coslon);    // E
//        disp[2] =  0.032 * sin(2.0*theta) * (m_eop.m_polarMotionX * coslon - m_eop.m_polarMotionY * sinlon);  // U
//        
//        /// transform  X=(x,y,z) into (R*X)(north,east,up)
//        //R(0,0) = -sa*co;  R(0,1) = -sa*so;  R(0,2) = ca;
//        //R(1,0) =    -so;  R(1,1) =     co;  R(1,2) = 0.;
//        //R(2,0) =  ca*co;  R(2,1) =  ca*so;  R(2,2) = sa;
//        
//        displacement.x = - sinlat*coslon*disp[0] - sinlat*sinlon*disp[1] + coslat*disp[2];
//        displacement.y =        - sinlon*disp[0] +        coslon*disp[1];
//        displacement.z =   coslat*coslon*disp[0] + coslat*sinlon*disp[1] + sinlat*disp[2];
//        
        
    } // end of function GEarthOrientationParameter::getPoleTide
    
    
    /*
     * ref: sofa iaubi00.c
     * the frame bias is a constant bias.
     */
    void GEarthOrientationParameter::getFrameBias00(double& dpsibi,double& depsbi,double& dra)
    {
        double DAS2R = GCONST("AS2R");
        /* The frame bias corrections in longitude and obliquity */
        const double DPBIAS = -0.041775  * DAS2R,
        DEBIAS = -0.0068192 * DAS2R;
        
        /* The ICRS RA of the J2000.0 equinox (Chapront et al., 2002) */
        const double DRA0 = -0.0146 * DAS2R;
        
        /* Return the results (which are fixed). */
        dpsibi = DPBIAS;
        depsbi = DEBIAS;
        dra = DRA0;
        
    }
    
    
    /*ref: sofa: iauPr00a*/
    void GEarthOrientationParameter::getPrecessionRate00(double JC_TT, double& dpsipr, double& depspr)
    {
        double T = JC_TT;
        static const double AS2R = GCONST("AS2R");
        static const double PRECOR = -0.29965 * AS2R; // Precession correction, radians per century
        static const double OBLCOR = -0.02524 * AS2R; //obliquity corrections (radians per century)
        /* Precession rate contributions with respect to IAU 1976/80. */
        dpsipr = PRECOR * T;
        depspr = OBLCOR * T;
    }
    
    void GEarthOrientationParameter::Rz(double phi,double* R)
    {
        double s = sin(phi);
        double c = cos(phi);
        R[0*3+0] = c;R[0*3+1] = s;R[0*3+2] = 0.0;
        R[1*3+0] = -s;R[1*3+1] = c;R[1*3+2] = 0.0;
        R[2*3+0] = 0.0;R[2*3+1] = 0.0;R[2*3+2] = 1.0;
    }
    
    void GEarthOrientationParameter::Ry(double phi,double* R)
    {
        double s = sin(phi);
        double c = cos(phi);
        R[0*3+0] = c;R[0*3+1] = 0.0;R[0*3+2] = -s;
        R[1*3+0] = 0.0;R[1*3+1] = 1.0;R[1*3+2] = 0.0;
        R[2*3+0] = s;R[2*3+1] = 0.0;R[2*3+2] = c;
    }
    void GEarthOrientationParameter::Rx(double phi,double* R)
    {
        double s = sin(phi);
        double c = cos(phi);
        R[0*3+0] = 1.0;R[0*3+1] = 0.0;R[0*3+2] = 0.0;
        R[1*3+0] = 0.0;R[1*3+1] = c;R[1*3+2] = s;
        R[2*3+0] = 0.0;R[2*3+1] = -s;R[2*3+2] =c;
    }
    
    /*
     * **  The CIO locator s, positioning the Celestial Intermediate Origin on
     **  the equator of the Celestial Intermediate Pole, given the CIP's X,Y
     **  coordinates.  Compatible with IAU 2006/2000A precession-nutation.
     *    JD_TT  double    TT as a 2-part Julian Date (Note 1)
     **     x,y           double    CIP coordinates (Note 3)
     * Returned (function value):
     **                   double    the CIO locator s in radians (Note 2)
     * The CIO locator s is the difference between the right ascensions
     **     of the same point in two systems:  the two systems are the GCRS
     **     and the CIP,CIO, and the point is the ascending node of the
     **     CIP equator.  The quantity s remains below 0.1 arcsecond
     **     throughout 1900-2100.
     * 3) The series used to compute s is in fact for s+XY/2, where X and Y
     **     are the x and y components of the CIP unit vector;  this series
     **     is more compact than a direct series for s would be.  This
     **     function requires X,Y to be supplied by the caller, who is
     **     responsible for providing values that are consistent with the
     **     supplied date.
     * REF:  sofa, s06.c
     */
    double GEarthOrientationParameter::getLocatorS06(double JC_TT,double x, double y)
    {
        /* Time since J2000.0, in Julian centuries */
        double T;
        double DAS2R = GCONST("AS2R");
        double D2PI = GCONST("PI")*2.0;
        /* Miscellaneous */
        int i, j;
        double a, w0, w1, w2, w3, w4, w5;
        
        /* Fundamental arguments */
        double fa[8];
        
        /* Returned value */
        double s;
        
        /* --------------------- */
        /* The series for s+XY/2 */
        /* --------------------- */
        
        typedef struct {
            int nfa[8];      /* coefficients of l,l',F,D,Om,LVe,LE,pA */
            double s, c;     /* sine and cosine coefficients */
        } TERM;
        
        /* Polynomial coefficients */
        static const double sp[] = {
            
            /* 1-6 */
            94.00e-6,
            3808.65e-6,
            -122.68e-6,
            -72574.11e-6,
            27.98e-6,
            15.62e-6
        };
        
        /* Terms of order t^0 */
        static const TERM s0[] = {
            
            /* 1-10 */
            {{ 0,  0,  0,  0,  1,  0,  0,  0}, -2640.73e-6,   0.39e-6 },
            {{ 0,  0,  0,  0,  2,  0,  0,  0},   -63.53e-6,   0.02e-6 },
            {{ 0,  0,  2, -2,  3,  0,  0,  0},   -11.75e-6,  -0.01e-6 },
            {{ 0,  0,  2, -2,  1,  0,  0,  0},   -11.21e-6,  -0.01e-6 },
            {{ 0,  0,  2, -2,  2,  0,  0,  0},     4.57e-6,   0.00e-6 },
            {{ 0,  0,  2,  0,  3,  0,  0,  0},    -2.02e-6,   0.00e-6 },
            {{ 0,  0,  2,  0,  1,  0,  0,  0},    -1.98e-6,   0.00e-6 },
            {{ 0,  0,  0,  0,  3,  0,  0,  0},     1.72e-6,   0.00e-6 },
            {{ 0,  1,  0,  0,  1,  0,  0,  0},     1.41e-6,   0.01e-6 },
            {{ 0,  1,  0,  0, -1,  0,  0,  0},     1.26e-6,   0.01e-6 },
            
            /* 11-20 */
            {{ 1,  0,  0,  0, -1,  0,  0,  0},     0.63e-6,   0.00e-6 },
            {{ 1,  0,  0,  0,  1,  0,  0,  0},     0.63e-6,   0.00e-6 },
            {{ 0,  1,  2, -2,  3,  0,  0,  0},    -0.46e-6,   0.00e-6 },
            {{ 0,  1,  2, -2,  1,  0,  0,  0},    -0.45e-6,   0.00e-6 },
            {{ 0,  0,  4, -4,  4,  0,  0,  0},    -0.36e-6,   0.00e-6 },
            {{ 0,  0,  1, -1,  1, -8, 12,  0},     0.24e-6,   0.12e-6 },
            {{ 0,  0,  2,  0,  0,  0,  0,  0},    -0.32e-6,   0.00e-6 },
            {{ 0,  0,  2,  0,  2,  0,  0,  0},    -0.28e-6,   0.00e-6 },
            {{ 1,  0,  2,  0,  3,  0,  0,  0},    -0.27e-6,   0.00e-6 },
            {{ 1,  0,  2,  0,  1,  0,  0,  0},    -0.26e-6,   0.00e-6 },
            
            /* 21-30 */
            {{ 0,  0,  2, -2,  0,  0,  0,  0},     0.21e-6,   0.00e-6 },
            {{ 0,  1, -2,  2, -3,  0,  0,  0},    -0.19e-6,   0.00e-6 },
            {{ 0,  1, -2,  2, -1,  0,  0,  0},    -0.18e-6,   0.00e-6 },
            {{ 0,  0,  0,  0,  0,  8,-13, -1},     0.10e-6,  -0.05e-6 },
            {{ 0,  0,  0,  2,  0,  0,  0,  0},    -0.15e-6,   0.00e-6 },
            {{ 2,  0, -2,  0, -1,  0,  0,  0},     0.14e-6,   0.00e-6 },
            {{ 0,  1,  2, -2,  2,  0,  0,  0},     0.14e-6,   0.00e-6 },
            {{ 1,  0,  0, -2,  1,  0,  0,  0},    -0.14e-6,   0.00e-6 },
            {{ 1,  0,  0, -2, -1,  0,  0,  0},    -0.14e-6,   0.00e-6 },
            {{ 0,  0,  4, -2,  4,  0,  0,  0},    -0.13e-6,   0.00e-6 },
            
            /* 31-33 */
            {{ 0,  0,  2, -2,  4,  0,  0,  0},     0.11e-6,   0.00e-6 },
            {{ 1,  0, -2,  0, -3,  0,  0,  0},    -0.11e-6,   0.00e-6 },
            {{ 1,  0, -2,  0, -1,  0,  0,  0},    -0.11e-6,   0.00e-6 }
        };
        
        /* Terms of order t^1 */
        static const TERM s1[] = {
            
            /* 1 - 3 */
            {{ 0,  0,  0,  0,  2,  0,  0,  0},    -0.07e-6,   3.57e-6 },
            {{ 0,  0,  0,  0,  1,  0,  0,  0},     1.73e-6,  -0.03e-6 },
            {{ 0,  0,  2, -2,  3,  0,  0,  0},     0.00e-6,   0.48e-6 }
        };
        
        /* Terms of order t^2 */
        static const TERM s2[] = {
            
            /* 1-10 */
            {{ 0,  0,  0,  0,  1,  0,  0,  0},   743.52e-6,  -0.17e-6 },
            {{ 0,  0,  2, -2,  2,  0,  0,  0},    56.91e-6,   0.06e-6 },
            {{ 0,  0,  2,  0,  2,  0,  0,  0},     9.84e-6,  -0.01e-6 },
            {{ 0,  0,  0,  0,  2,  0,  0,  0},    -8.85e-6,   0.01e-6 },
            {{ 0,  1,  0,  0,  0,  0,  0,  0},    -6.38e-6,  -0.05e-6 },
            {{ 1,  0,  0,  0,  0,  0,  0,  0},    -3.07e-6,   0.00e-6 },
            {{ 0,  1,  2, -2,  2,  0,  0,  0},     2.23e-6,   0.00e-6 },
            {{ 0,  0,  2,  0,  1,  0,  0,  0},     1.67e-6,   0.00e-6 },
            {{ 1,  0,  2,  0,  2,  0,  0,  0},     1.30e-6,   0.00e-6 },
            {{ 0,  1, -2,  2, -2,  0,  0,  0},     0.93e-6,   0.00e-6 },
            
            /* 11-20 */
            {{ 1,  0,  0, -2,  0,  0,  0,  0},     0.68e-6,   0.00e-6 },
            {{ 0,  0,  2, -2,  1,  0,  0,  0},    -0.55e-6,   0.00e-6 },
            {{ 1,  0, -2,  0, -2,  0,  0,  0},     0.53e-6,   0.00e-6 },
            {{ 0,  0,  0,  2,  0,  0,  0,  0},    -0.27e-6,   0.00e-6 },
            {{ 1,  0,  0,  0,  1,  0,  0,  0},    -0.27e-6,   0.00e-6 },
            {{ 1,  0, -2, -2, -2,  0,  0,  0},    -0.26e-6,   0.00e-6 },
            {{ 1,  0,  0,  0, -1,  0,  0,  0},    -0.25e-6,   0.00e-6 },
            {{ 1,  0,  2,  0,  1,  0,  0,  0},     0.22e-6,   0.00e-6 },
            {{ 2,  0,  0, -2,  0,  0,  0,  0},    -0.21e-6,   0.00e-6 },
            {{ 2,  0, -2,  0, -1,  0,  0,  0},     0.20e-6,   0.00e-6 },
            
            /* 21-25 */
            {{ 0,  0,  2,  2,  2,  0,  0,  0},     0.17e-6,   0.00e-6 },
            {{ 2,  0,  2,  0,  2,  0,  0,  0},     0.13e-6,   0.00e-6 },
            {{ 2,  0,  0,  0,  0,  0,  0,  0},    -0.13e-6,   0.00e-6 },
            {{ 1,  0,  2, -2,  2,  0,  0,  0},    -0.12e-6,   0.00e-6 },
            {{ 0,  0,  2,  0,  0,  0,  0,  0},    -0.11e-6,   0.00e-6 }
        };
        
        /* Terms of order t^3 */
        static const TERM s3[] = {
            
            /* 1-4 */
            {{ 0,  0,  0,  0,  1,  0,  0,  0},     0.30e-6, -23.42e-6 },
            {{ 0,  0,  2, -2,  2,  0,  0,  0},    -0.03e-6,  -1.46e-6 },
            {{ 0,  0,  2,  0,  2,  0,  0,  0},    -0.01e-6,  -0.25e-6 },
            {{ 0,  0,  0,  0,  2,  0,  0,  0},     0.00e-6,   0.23e-6 }
        };
        
        /* Terms of order t^4 */
        static const TERM s4[] = {
            
            /* 1-1 */
            {{ 0,  0,  0,  0,  1,  0,  0,  0},    -0.26e-6,  -0.01e-6 }
        };
        
        /* Number of terms in the series */
        static const int NS0 = (int) (sizeof s0 / sizeof (TERM));
        static const int NS1 = (int) (sizeof s1 / sizeof (TERM));
        static const int NS2 = (int) (sizeof s2 / sizeof (TERM));
        static const int NS3 = (int) (sizeof s3 / sizeof (TERM));
        static const int NS4 = (int) (sizeof s4 / sizeof (TERM));
        
        /*--------------------------------------------------------------------*/
        
        /* Interval between fundamental epoch J2000.0 and current date (JC). */
        T = JC_TT;
        
        /* Fundamental Arguments (from IERS Conventions 2003) */
        
        /* Mean anomaly of the Moon. */
        /* Mean anomaly of the Moon (IERS Conventions 2003). */
        fa[0] = fmod(485868.249036 +T*(1717915923.2178 + T*(31.8792 +T*(0.051635 +T * (- 0.00024470 ) ) ) ), 1296000.0 ) * DAS2R;
        
        /* Mean anomaly of the Sun. */
        /* Mean anomaly of the Sun (IERS Conventions 2003). */
        fa[1] = fmod(1287104.793048 + T*(129596581.0481 +T*(- 0.5532 +T * (0.000136 +T * (- 0.00001149 ) ) ) ), 1296000.0)*DAS2R;
        
        /* Mean longitude of the Moon minus that of the ascending node. */
        /* Mean longitude of the Moon minus that of the ascending node */
        /* (IERS Conventions 2003).                                    */
        fa[2] = fmod(335779.526232 +T * ( 1739527262.8478 +T*(- 12.7512 +T*(- 0.001037 +T*(0.00000417 ) ) ) ), 1296000.0 ) * DAS2R;
        /* Mean elongation of the Moon from the Sun. */
        /* Mean elongation of the Moon from the Sun (IERS Conventions 2003). */
        fa[3] = fmod(1072260.703692 +T*( 1602961601.2090 +T*(- 6.3706 +T*(0.006593 + T*( -0.00003169 ) ) ) ), 1296000.0 ) * DAS2R;
        /* Mean longitude of the ascending node of the Moon. */
        /* Mean longitude of the Moon's ascending node */
        /* (IERS Conventions 2003).                    */
        fa[4] = fmod(450160.398036 +T*(- 6962890.5431 +T*(7.4722 +T*(0.007702 +T*(- 0.00005939 ) ) ) ), 1296000.0 ) * DAS2R;
        /* Mean longitude of Venus. */
        /* Mean longitude of Venus (IERS Conventions 2003). */
        fa[5] = fmod(3.176146697 + 1021.3285546211 * T, D2PI);
        /* Mean longitude of Earth. */
        /* Mean longitude of Earth (IERS Conventions 2003). */
        fa[6] = fmod(1.753470314 + 628.3075849991 * T, D2PI);
        /* General precession in longitude. */
        /* General accumulated precession in longitude. */
        fa[7] = (0.024381750 + 0.00000538691 * T) * T;
        /* Evaluate s. */
        w0 = sp[0];
        w1 = sp[1];
        w2 = sp[2];
        w3 = sp[3];
        w4 = sp[4];
        w5 = sp[5];
        
        for (i = NS0-1; i >= 0; i--) {
            a = 0.0;
            for (j = 0; j < 8; j++) {
                a += (double)s0[i].nfa[j] * fa[j];
            }
            w0 += s0[i].s * sin(a) + s0[i].c * cos(a);
        }
        
        for (i = NS1-1; i >= 0; i--) {
            a = 0.0;
            for (j = 0; j < 8; j++) {
                a += (double)s1[i].nfa[j] * fa[j];
            }
            w1 += s1[i].s * sin(a) + s1[i].c * cos(a);
        }
        
        for (i = NS2-1; i >= 0; i--) {
            a = 0.0;
            for (j = 0; j < 8; j++) {
                a += (double)s2[i].nfa[j] * fa[j];
            }
            w2 += s2[i].s * sin(a) + s2[i].c * cos(a);
        }
        
        for (i = NS3-1; i >= 0; i--) {
            a = 0.0;
            for (j = 0; j < 8; j++) {
                a += (double)s3[i].nfa[j] * fa[j];
            }
            w3 += s3[i].s * sin(a) + s3[i].c * cos(a);
        }
        
        for (i = NS4-1; i >= 0; i--)
        {
            a = 0.0;
            for (j = 0; j < 8; j++)
            {
                a += (double)s4[i].nfa[j] * fa[j];
            }
            w4 += s4[i].s * sin(a) + s4[i].c * cos(a);
        }
        
        s = (w0 +
             (w1 +
              (w2 +
               (w3 +
                (w4 +
                 w5 * T) * T) * T) * T) * T) * DAS2R - x*y/2.0;
        
        return s;
        
    }
    
    
    /*get the rotation matrix with frame bias, precession and nutation
     *  ref: sofa, pnm06.c
     *
     */
    void GEarthOrientationParameter::getRotationBPN(double JC_TT,double* BPN)
    {
        
        double gamb=0.0, phib =0.0, psib =0.0, epsa =0.0, dp =0.0,de=0.0;
        
        getPrecessionAngle(JC_TT, gamb, phib, psib, epsa);
        
        iauNutation2006A(JC_TT, dp, de);
        
        double Rz_gamb[9]= {0.0}, Rx_phib[9]={0.0},Rz_psi[9]={0.0},Rx_eps[9]={0.0};
        Rz(gamb,Rz_gamb);
        Rx(phib,Rx_phib);
        Rz(-(psib+dp), Rz_psi);
        Rx(-(epsa+de),Rx_eps);
        double temp1[9]={0.0},temp2[9]={0.0};
        matrixMultiply(3, 3, Rx_phib,3, Rz_gamb, temp1);
        matrixMultiply(3, 3, Rz_psi, 3, temp1, temp2);
        matrixMultiply(3, 3, Rx_eps, 3, temp2, BPN);
        
    }
    
    
    /* the nutation 2006A model , which adds some corrections to the Nutation2000A model
     *
     *  ref: sofa iauNut06a.c
     */
    void GEarthOrientationParameter::iauNutation2006A(double JC_TT, double& dpsi, double& deps)
    {
        /* Interval between fundamental date J2000.0 and given date (JC). */
       double T = JC_TT;
        /* Factor correcting for secular variation of J2. */
       double fj2 = -2.7774e-6 * T;
        
        double dp =0.0, de =0.0;
       /* Obtain IAU 2000A nutation. */
       iauNutation2000A(JC_TT, dp, de);
       
        /* Apply P03 adjustments (Wallace & Capitaine, 2006, Eqs.5). */
        dpsi = dp + dp * (0.4697E-6 + fj2);
        deps = de + de * fj2;
        
    }
    
    /*
     ref: sofa, iauPfw06
     **  Precession angles, IAU 2006 (Fukushima-Williams 4-angle formulation).
     */
    void GEarthOrientationParameter::getPrecessionAngle(double JC_TT, double & gamb, double & phib, double &psib, double& epsa)
    {
        double T = JC_TT;
        static const double DAS2R = GCONST("AS2R");
        /* P03 bias+precession angles. */
        gamb = (    -0.052928     +
                 (    10.556378     +
                  (     0.4932044    +
                   (    -0.00031238   +
                    (    -0.000002788  +
                     (     0.0000000260 )
                     * T) * T) * T) * T) * T) * DAS2R;
        phib = ( 84381.412819     +
                 (   -46.811016     +
                  (     0.0511268    +
                   (     0.00053289   +
                    (    -0.000000440  +
                     (    -0.0000000176 )
                     * T) * T) * T) * T) * T) * DAS2R;
        psib = (    -0.041775     +
                 (  5038.481484     +
                  (     1.5584175    +
                   (    -0.00018522   +
                    (    -0.000026452  +
                     (    -0.0000000148 )
                     * T) * T) * T) * T) * T) * DAS2R;
        
        
        epsa =  GIERS::getMeanObliquity(JC_TT);
        
        
    }
    
    
    /*get the precession angle given Julian century in TDB(since J2000)*/
    /*ref: Expressions for IAU2000 precession quantities, N.Capitaine, P.T.Wallace and J.Chapront, Astronomy& Astrophysics,2003,page 572*/
    /*IAU2000A precession model*/
    /* IERS2003, tech note 32
     *  T = (TT − 2000 January 1d 12h TT)/36 525
     */
    void GEarthOrientationParameter::getPrecessionAngle( double JC_TT, double &deltaA, double &zA, double &thetaA)
    {
        double AS2R = GCONST("AS2R"); // arc second to radian
        double T = JC_TT;  // Julian Centuries in TT
        
//        //*  J2000 obliquity (Lieske et al. 1977)
        double EPS0 = 84381.448 * AS2R;
        //IAU 1980 mean obliquity of date
        //double EPSA80 = EPS0 + (  -46.8150 +(   -0.00059 + (    0.001813 ) * T ) * T ) * T * AS2R;
        //Precession rate contributions with respect to IAU 1976/80.
        //double PSIA77 =        ( 5038.7784 +(   -1.07259 +(   -0.001147 ) * T ) * T ) * T * AS2R;
        //double OMA77  = EPS0 + ((    0.05127 +(   -0.007726 ) * T ) * T ) * T * AS2R;
        //double CHIA   =        (   10.5526 +(   -2.38064 +(   -0.001125 ) * T ) * T ) * T * AS2R;
        
        //Apply IAU 2000A precession rate corrections.
        double PRECOR = -0.29965 * AS2R;
        double OBLCOR = -0.02524 * AS2R;
        //double DPSIPR = PRECOR * T;
        //double DEPSPR = OBLCOR * T;
        //double PSIA = PSIA77 + DPSIPR;
        //double OMA  = OMA77  + DEPSPR;
        //double EPSA = EPSA80 + DEPSPR;
        
        
        
//        // the other way
//        deltaA = 2.5976176 + 2306.0809506*T+0.3019015*T*T + 0.0179663*T*T*T - 0.0000327*T*T*T*T - 0.0000002*T*T*T*T*T;
//        zA  = -2.5976176+2306.0803226*T+1.0947790*T*T +0.0182273*T*T*T+0.0000470*T*T*T*T - 0.0000003*T*T*T*T*T;
//        thetaA = 2004.1917476*T - 0.4269353*T*T - 0.0418251*T*T*T - 0.0000601*T*T*T*T - 0.0000001*T*T*T*T*T;
//        
         //transfer from arcsecond to radian
        deltaA = deltaA*AS2R;
        zA     = zA*AS2R;
        thetaA = thetaA*AS2R;
        
    }
    
    /*get the precession rotation matrix */
    void GEarthOrientationParameter::getPrecessionMatrix(double JC_TT,double* pr)
    {
        memset(pr,0.0,sizeof(double)*9); // 3 rows and 3 columns
        double deltaA =0.0, zA =0.0, thetaA =0.0;
        getPrecessionAngle(JC_TT, deltaA, zA, thetaA);
        //pr = R(-zA)*R(+thetaA)*R(-deltaA);
        double rZ[9]={0.0},rT[9]={0.0},rX[9]={0.0};
        //-zA
        rZ[0*3+0]=cos(-zA);rZ[0*3+1]=sin(-zA);rZ[0*3+2]=0.0;
        rZ[1*3+0]=-sin(-zA);rZ[1*3+1]=cos(-zA);rZ[1*3+2]=0.0;
        rZ[2*3+0]=0.0;rZ[2*3+1]=0.0;rZ[2*3+2]=1.0;
        //thetaA
        rT[0*3+0]=cos(thetaA);rT[0*3+1]=0.0;rT[0*3+2]=-sin(thetaA);
        rT[1*3+0]=0.0;rT[1*3+1]=1.0;rT[1*3+2]=0.0;
        rT[2*3+0]=sin(thetaA);rT[2*3+1]=0.0;rT[2*3+2]=cos(thetaA);
        //-deltaA
        rX[0*3+0]=cos(-deltaA);rX[0*3+1]=sin(-deltaA);rX[0*3+2]=0.0;
        rX[1*3+0]=-sin(-deltaA);rX[1*3+1]=cos(-deltaA);rX[1*3+2]=0.0;
        rX[2*3+0]=0.0;rX[2*3+1]=0.0;rX[2*3+2]=1.0;
        
        double rZrT[9]={0.0}; // rZ*rT
        
        matrixMultiply(3, 3, rZ, 3, rT, rZrT);
        
        matrixMultiply(3, 3, rZrT, 3, rX, pr);
        
    } // end of the function
    
    
    
    void GEarthOrientationParameter::iauNutation2000A(double JC_TT,double& dpsi,double& deps)
    {
        short int i;
        double T, a[5], dp, de, arg, sarg, carg, factor, dpsils, depsls,
        al, alsu, af, ad, aom, alme, alve, alea, alma, alju, alsa, alur,
        alne, apa, dpsipl, depspl;
        double TWOPI = GCONST("PI")*2.0;
        double AS2R = GCONST("AS2R");
        
        /*
         Luni-Solar argument multipliers:
         L     L'    F     D     Om
         */
        
        static const short int nals_t[678][5] =
        {
            { 0,    0,    0,    0,    1},
            { 0,    0,    2,   -2,    2},
            { 0,    0,    2,    0,    2},
            { 0,    0,    0,    0,    2},
            { 0,    1,    0,    0,    0},
            { 0,    1,    2,   -2,    2},
            { 1,    0,    0,    0,    0},
            { 0,    0,    2,    0,    1},
            { 1,    0,    2,    0,    2},
            { 0,   -1,    2,   -2,    2},
            { 0,    0,    2,   -2,    1},
            {-1,    0,    2,    0,    2},
            {-1,    0,    0,    2,    0},
            { 1,    0,    0,    0,    1},
            {-1,    0,    0,    0,    1},
            {-1,    0,    2,    2,    2},
            { 1,    0,    2,    0,    1},
            {-2,    0,    2,    0,    1},
            { 0,    0,    0,    2,    0},
            { 0,    0,    2,    2,    2},
            { 0,   -2,    2,   -2,    2},
            {-2,    0,    0,    2,    0},
            { 2,    0,    2,    0,    2},
            { 1,    0,    2,   -2,    2},
            {-1,    0,    2,    0,    1},
            { 2,    0,    0,    0,    0},
            { 0,    0,    2,    0,    0},
            { 0,    1,    0,    0,    1},
            {-1,    0,    0,    2,    1},
            { 0,    2,    2,   -2,    2},
            { 0,    0,   -2,    2,    0},
            { 1,    0,    0,   -2,    1},
            { 0,   -1,    0,    0,    1},
            {-1,    0,    2,    2,    1},
            { 0,    2,    0,    0,    0},
            { 1,    0,    2,    2,    2},
            {-2,    0,    2,    0,    0},
            { 0,    1,    2,    0,    2},
            { 0,    0,    2,    2,    1},
            { 0,   -1,    2,    0,    2},
            { 0,    0,    0,    2,    1},
            { 1,    0,    2,   -2,    1},
            { 2,    0,    2,   -2,    2},
            {-2,    0,    0,    2,    1},
            { 2,    0,    2,    0,    1},
            { 0,   -1,    2,   -2,    1},
            { 0,    0,    0,   -2,    1},
            {-1,   -1,    0,    2,    0},
            { 2,    0,    0,   -2,    1},
            { 1,    0,    0,    2,    0},
            { 0,    1,    2,   -2,    1},
            { 1,   -1,    0,    0,    0},
            {-2,    0,    2,    0,    2},
            { 3,    0,    2,    0,    2},
            { 0,   -1,    0,    2,    0},
            { 1,   -1,    2,    0,    2},
            { 0,    0,    0,    1,    0},
            {-1,   -1,    2,    2,    2},
            {-1,    0,    2,    0,    0},
            { 0,   -1,    2,    2,    2},
            {-2,    0,    0,    0,    1},
            { 1,    1,    2,    0,    2},
            { 2,    0,    0,    0,    1},
            {-1,    1,    0,    1,    0},
            { 1,    1,    0,    0,    0},
            { 1,    0,    2,    0,    0},
            {-1,    0,    2,   -2,    1},
            { 1,    0,    0,    0,    2},
            {-1,    0,    0,    1,    0},
            { 0,    0,    2,    1,    2},
            {-1,    0,    2,    4,    2},
            {-1,    1,    0,    1,    1},
            { 0,   -2,    2,   -2,    1},
            { 1,    0,    2,    2,    1},
            {-2,    0,    2,    2,    2},
            {-1,    0,    0,    0,    2},
            { 1,    1,    2,   -2,    2},
            {-2,    0,    2,    4,    2},
            {-1,    0,    4,    0,    2},
            { 2,    0,    2,   -2,    1},
            { 2,    0,    2,    2,    2},
            { 1,    0,    0,    2,    1},
            { 3,    0,    0,    0,    0},
            { 3,    0,    2,   -2,    2},
            { 0,    0,    4,   -2,    2},
            { 0,    1,    2,    0,    1},
            { 0,    0,   -2,    2,    1},
            { 0,    0,    2,   -2,    3},
            {-1,    0,    0,    4,    0},
            { 2,    0,   -2,    0,    1},
            {-2,    0,    0,    4,    0},
            {-1,   -1,    0,    2,    1},
            {-1,    0,    0,    1,    1},
            { 0,    1,    0,    0,    2},
            { 0,    0,   -2,    0,    1},
            { 0,   -1,    2,    0,    1},
            { 0,    0,    2,   -1,    2},
            { 0,    0,    2,    4,    2},
            {-2,   -1,    0,    2,    0},
            { 1,    1,    0,   -2,    1},
            {-1,    1,    0,    2,    0},
            {-1,    1,    0,    1,    2},
            { 1,   -1,    0,    0,    1},
            { 1,   -1,    2,    2,    2},
            {-1,    1,    2,    2,    2},
            { 3,    0,    2,    0,    1},
            { 0,    1,   -2,    2,    0},
            {-1,    0,    0,   -2,    1},
            { 0,    1,    2,    2,    2},
            {-1,   -1,    2,    2,    1},
            { 0,   -1,    0,    0,    2},
            { 1,    0,    2,   -4,    1},
            {-1,    0,   -2,    2,    0},
            { 0,   -1,    2,    2,    1},
            { 2,   -1,    2,    0,    2},
            { 0,    0,    0,    2,    2},
            { 1,   -1,    2,    0,    1},
            {-1,    1,    2,    0,    2},
            { 0,    1,    0,    2,    0},
            { 0,   -1,   -2,    2,    0},
            { 0,    3,    2,   -2,    2},
            { 0,    0,    0,    1,    1},
            {-1,    0,    2,    2,    0},
            { 2,    1,    2,    0,    2},
            { 1,    1,    0,    0,    1},
            { 1,    1,    2,    0,    1},
            { 2,    0,    0,    2,    0},
            { 1,    0,   -2,    2,    0},
            {-1,    0,    0,    2,    2},
            { 0,    1,    0,    1,    0},
            { 0,    1,    0,   -2,    1},
            {-1,    0,    2,   -2,    2},
            { 0,    0,    0,   -1,    1},
            {-1,    1,    0,    0,    1},
            { 1,    0,    2,   -1,    2},
            { 1,   -1,    0,    2,    0},
            { 0,    0,    0,    4,    0},
            { 1,    0,    2,    1,    2},
            { 0,    0,    2,    1,    1},
            { 1,    0,    0,   -2,    2},
            {-1,    0,    2,    4,    1},
            { 1,    0,   -2,    0,    1},
            { 1,    1,    2,   -2,    1},
            { 0,    0,    2,    2,    0},
            {-1,    0,    2,   -1,    1},
            {-2,    0,    2,    2,    1},
            { 4,    0,    2,    0,    2},
            { 2,   -1,    0,    0,    0},
            { 2,    1,    2,   -2,    2},
            { 0,    1,    2,    1,    2},
            { 1,    0,    4,   -2,    2},
            {-1,   -1,    0,    0,    1},
            { 0,    1,    0,    2,    1},
            {-2,    0,    2,    4,    1},
            { 2,    0,    2,    0,    0},
            { 1,    0,    0,    1,    0},
            {-1,    0,    0,    4,    1},
            {-1,    0,    4,    0,    1},
            { 2,    0,    2,    2,    1},
            { 0,    0,    2,   -3,    2},
            {-1,   -2,    0,    2,    0},
            { 2,    1,    0,    0,    0},
            { 0,    0,    4,    0,    2},
            { 0,    0,    0,    0,    3},
            { 0,    3,    0,    0,    0},
            { 0,    0,    2,   -4,    1},
            { 0,   -1,    0,    2,    1},
            { 0,    0,    0,    4,    1},
            {-1,   -1,    2,    4,    2},
            { 1,    0,    2,    4,    2},
            {-2,    2,    0,    2,    0},
            {-2,   -1,    2,    0,    1},
            {-2,    0,    0,    2,    2},
            {-1,   -1,    2,    0,    2},
            { 0,    0,    4,   -2,    1},
            { 3,    0,    2,   -2,    1},
            {-2,   -1,    0,    2,    1},
            { 1,    0,    0,   -1,    1},
            { 0,   -2,    0,    2,    0},
            {-2,    0,    0,    4,    1},
            {-3,    0,    0,    0,    1},
            { 1,    1,    2,    2,    2},
            { 0,    0,    2,    4,    1},
            { 3,    0,    2,    2,    2},
            {-1,    1,    2,   -2,    1},
            { 2,    0,    0,   -4,    1},
            { 0,    0,    0,   -2,    2},
            { 2,    0,    2,   -4,    1},
            {-1,    1,    0,    2,    1},
            { 0,    0,    2,   -1,    1},
            { 0,   -2,    2,    2,    2},
            { 2,    0,    0,    2,    1},
            { 4,    0,    2,   -2,    2},
            { 2,    0,    0,   -2,    2},
            { 0,    2,    0,    0,    1},
            { 1,    0,    0,   -4,    1},
            { 0,    2,    2,   -2,    1},
            {-3,    0,    0,    4,    0},
            {-1,    1,    2,    0,    1},
            {-1,   -1,    0,    4,    0},
            {-1,   -2,    2,    2,    2},
            {-2,   -1,    2,    4,    2},
            { 1,   -1,    2,    2,    1},
            {-2,    1,    0,    2,    0},
            {-2,    1,    2,    0,    1},
            { 2,    1,    0,   -2,    1},
            {-3,    0,    2,    0,    1},
            {-2,    0,    2,   -2,    1},
            {-1,    1,    0,    2,    2},
            { 0,   -1,    2,   -1,    2},
            {-1,    0,    4,   -2,    2},
            { 0,   -2,    2,    0,    2},
            {-1,    0,    2,    1,    2},
            { 2,    0,    0,    0,    2},
            { 0,    0,    2,    0,    3},
            {-2,    0,    4,    0,    2},
            {-1,    0,   -2,    0,    1},
            {-1,    1,    2,    2,    1},
            { 3,    0,    0,    0,    1},
            {-1,    0,    2,    3,    2},
            { 2,   -1,    2,    0,    1},
            { 0,    1,    2,    2,    1},
            { 0,   -1,    2,    4,    2},
            { 2,   -1,    2,    2,    2},
            { 0,    2,   -2,    2,    0},
            {-1,   -1,    2,   -1,    1},
            { 0,   -2,    0,    0,    1},
            { 1,    0,    2,   -4,    2},
            { 1,   -1,    0,   -2,    1},
            {-1,   -1,    2,    0,    1},
            { 1,   -1,    2,   -2,    2},
            {-2,   -1,    0,    4,    0},
            {-1,    0,    0,    3,    0},
            {-2,   -1,    2,    2,    2},
            { 0,    2,    2,    0,    2},
            { 1,    1,    0,    2,    0},
            { 2,    0,    2,   -1,    2},
            { 1,    0,    2,    1,    1},
            { 4,    0,    0,    0,    0},
            { 2,    1,    2,    0,    1},
            { 3,   -1,    2,    0,    2},
            {-2,    2,    0,    2,    1},
            { 1,    0,    2,   -3,    1},
            { 1,    1,    2,   -4,    1},
            {-1,   -1,    2,   -2,    1},
            { 0,   -1,    0,   -1,    1},
            { 0,   -1,    0,   -2,    1},
            {-2,    0,    0,    0,    2},
            {-2,    0,   -2,    2,    0},
            {-1,    0,   -2,    4,    0},
            { 1,   -2,    0,    0,    0},
            { 0,    1,    0,    1,    1},
            {-1,    2,    0,    2,    0},
            { 1,   -1,    2,   -2,    1},
            { 1,    2,    2,   -2,    2},
            { 2,   -1,    2,   -2,    2},
            { 1,    0,    2,   -1,    1},
            { 2,    1,    2,   -2,    1},
            {-2,    0,    0,   -2,    1},
            { 1,   -2,    2,    0,    2},
            { 0,    1,    2,    1,    1},
            { 1,    0,    4,   -2,    1},
            {-2,    0,    4,    2,    2},
            { 1,    1,    2,    1,    2},
            { 1,    0,    0,    4,    0},
            { 1,    0,    2,    2,    0},
            { 2,    0,    2,    1,    2},
            { 3,    1,    2,    0,    2},
            { 4,    0,    2,    0,    1},
            {-2,   -1,    2,    0,    0},
            { 0,    1,   -2,    2,    1},
            { 1,    0,   -2,    1,    0},
            { 0,   -1,   -2,    2,    1},
            { 2,   -1,    0,   -2,    1},
            {-1,    0,    2,   -1,    2},
            { 1,    0,    2,   -3,    2},
            { 0,    1,    2,   -2,    3},
            { 0,    0,    2,   -3,    1},
            {-1,    0,   -2,    2,    1},
            { 0,    0,    2,   -4,    2},
            {-2,    1,    0,    0,    1},
            {-1,    0,    0,   -1,    1},
            { 2,    0,    2,   -4,    2},
            { 0,    0,    4,   -4,    4},
            { 0,    0,    4,   -4,    2},
            {-1,   -2,    0,    2,    1},
            {-2,    0,    0,    3,    0},
            { 1,    0,   -2,    2,    1},
            {-3,    0,    2,    2,    2},
            {-3,    0,    2,    2,    1},
            {-2,    0,    2,    2,    0},
            { 2,   -1,    0,    0,    1},
            {-2,    1,    2,    2,    2},
            { 1,    1,    0,    1,    0},
            { 0,    1,    4,   -2,    2},
            {-1,    1,    0,   -2,    1},
            { 0,    0,    0,   -4,    1},
            { 1,   -1,    0,    2,    1},
            { 1,    1,    0,    2,    1},
            {-1,    2,    2,    2,    2},
            { 3,    1,    2,   -2,    2},
            { 0,   -1,    0,    4,    0},
            { 2,   -1,    0,    2,    0},
            { 0,    0,    4,    0,    1},
            { 2,    0,    4,   -2,    2},
            {-1,   -1,    2,    4,    1},
            { 1,    0,    0,    4,    1},
            { 1,   -2,    2,    2,    2},
            { 0,    0,    2,    3,    2},
            {-1,    1,    2,    4,    2},
            { 3,    0,    0,    2,    0},
            {-1,    0,    4,    2,    2},
            { 1,    1,    2,    2,    1},
            {-2,    0,    2,    6,    2},
            { 2,    1,    2,    2,    2},
            {-1,    0,    2,    6,    2},
            { 1,    0,    2,    4,    1},
            { 2,    0,    2,    4,    2},
            { 1,    1,   -2,    1,    0},
            {-3,    1,    2,    1,    2},
            { 2,    0,   -2,    0,    2},
            {-1,    0,    0,    1,    2},
            {-4,    0,    2,    2,    1},
            {-1,   -1,    0,    1,    0},
            { 0,    0,   -2,    2,    2},
            { 1,    0,    0,   -1,    2},
            { 0,   -1,    2,   -2,    3},
            {-2,    1,    2,    0,    0},
            { 0,    0,    2,   -2,    4},
            {-2,   -2,    0,    2,    0},
            {-2,    0,   -2,    4,    0},
            { 0,   -2,   -2,    2,    0},
            { 1,    2,    0,   -2,    1},
            { 3,    0,    0,   -4,    1},
            {-1,    1,    2,   -2,    2},
            { 1,   -1,    2,   -4,    1},
            { 1,    1,    0,   -2,    2},
            {-3,    0,    2,    0,    0},
            {-3,    0,    2,    0,    2},
            {-2,    0,    0,    1,    0},
            { 0,    0,   -2,    1,    0},
            {-3,    0,    0,    2,    1},
            {-1,   -1,   -2,    2,    0},
            { 0,    1,    2,   -4,    1},
            { 2,    1,    0,   -4,    1},
            { 0,    2,    0,   -2,    1},
            { 1,    0,    0,   -3,    1},
            {-2,    0,    2,   -2,    2},
            {-2,   -1,    0,    0,    1},
            {-4,    0,    0,    2,    0},
            { 1,    1,    0,   -4,    1},
            {-1,    0,    2,   -4,    1},
            { 0,    0,    4,   -4,    1},
            { 0,    3,    2,   -2,    2},
            {-3,   -1,    0,    4,    0},
            {-3,    0,    0,    4,    1},
            { 1,   -1,   -2,    2,    0},
            {-1,   -1,    0,    2,    2},
            { 1,   -2,    0,    0,    1},
            { 1,   -1,    0,    0,    2},
            { 0,    0,    0,    1,    2},
            {-1,   -1,    2,    0,    0},
            { 1,   -2,    2,   -2,    2},
            { 0,   -1,    2,   -1,    1},
            {-1,    0,    2,    0,    3},
            { 1,    1,    0,    0,    2},
            {-1,    1,    2,    0,    0},
            { 1,    2,    0,    0,    0},
            {-1,    2,    2,    0,    2},
            {-1,    0,    4,   -2,    1},
            { 3,    0,    2,   -4,    2},
            { 1,    2,    2,   -2,    1},
            { 1,    0,    4,   -4,    2},
            {-2,   -1,    0,    4,    1},
            { 0,   -1,    0,    2,    2},
            {-2,    1,    0,    4,    0},
            {-2,   -1,    2,    2,    1},
            { 2,    0,   -2,    2,    0},
            { 1,    0,    0,    1,    1},
            { 0,    1,    0,    2,    2},
            { 1,   -1,    2,   -1,    2},
            {-2,    0,    4,    0,    1},
            { 2,    1,    0,    0,    1},
            { 0,    1,    2,    0,    0},
            { 0,   -1,    4,   -2,    2},
            { 0,    0,    4,   -2,    4},
            { 0,    2,    2,    0,    1},
            {-3,    0,    0,    6,    0},
            {-1,   -1,    0,    4,    1},
            { 1,   -2,    0,    2,    0},
            {-1,    0,    0,    4,    2},
            {-1,   -2,    2,    2,    1},
            {-1,    0,    0,   -2,    2},
            { 1,    0,   -2,   -2,    1},
            { 0,    0,   -2,   -2,    1},
            {-2,    0,   -2,    0,    1},
            { 0,    0,    0,    3,    1},
            { 0,    0,    0,    3,    0},
            {-1,    1,    0,    4,    0},
            {-1,   -1,    2,    2,    0},
            {-2,    0,    2,    3,    2},
            { 1,    0,    0,    2,    2},
            { 0,   -1,    2,    1,    2},
            { 3,   -1,    0,    0,    0},
            { 2,    0,    0,    1,    0},
            { 1,   -1,    2,    0,    0},
            { 0,    0,    2,    1,    0},
            { 1,    0,    2,    0,    3},
            { 3,    1,    0,    0,    0},
            { 3,   -1,    2,   -2,    2},
            { 2,    0,    2,   -1,    1},
            { 1,    1,    2,    0,    0},
            { 0,    0,    4,   -1,    2},
            { 1,    2,    2,    0,    2},
            {-2,    0,    0,    6,    0},
            { 0,   -1,    0,    4,    1},
            {-2,   -1,    2,    4,    1},
            { 0,   -2,    2,    2,    1},
            { 0,   -1,    2,    2,    0},
            {-1,    0,    2,    3,    1},
            {-2,    1,    2,    4,    2},
            { 2,    0,    0,    2,    2},
            { 2,   -2,    2,    0,    2},
            {-1,    1,    2,    3,    2},
            { 3,    0,    2,   -1,    2},
            { 4,    0,    2,   -2,    1},
            {-1,    0,    0,    6,    0},
            {-1,   -2,    2,    4,    2},
            {-3,    0,    2,    6,    2},
            {-1,    0,    2,    4,    0},
            { 3,    0,    0,    2,    1},
            { 3,   -1,    2,    0,    1},
            { 3,    0,    2,    0,    0},
            { 1,    0,    4,    0,    2},
            { 5,    0,    2,   -2,    2},
            { 0,   -1,    2,    4,    1},
            { 2,   -1,    2,    2,    1},
            { 0,    1,    2,    4,    2},
            { 1,   -1,    2,    4,    2},
            { 3,   -1,    2,    2,    2},
            { 3,    0,    2,    2,    1},
            { 5,    0,    2,    0,    2},
            { 0,    0,    2,    6,    2},
            { 4,    0,    2,    2,    2},
            { 0,   -1,    1,   -1,    1},
            {-1,    0,    1,    0,    3},
            { 0,   -2,    2,   -2,    3},
            { 1,    0,   -1,    0,    1},
            { 2,   -2,    0,   -2,    1},
            {-1,    0,    1,    0,    2},
            {-1,    0,    1,    0,    1},
            {-1,   -1,    2,   -1,    2},
            {-2,    2,    0,    2,    2},
            {-1,    0,    1,    0,    0},
            {-4,    1,    2,    2,    2},
            {-3,    0,    2,    1,    1},
            {-2,   -1,    2,    0,    2},
            { 1,    0,   -2,    1,    1},
            { 2,   -1,   -2,    0,    1},
            {-4,    0,    2,    2,    0},
            {-3,    1,    0,    3,    0},
            {-1,    0,   -1,    2,    0},
            { 0,   -2,    0,    0,    2},
            { 0,   -2,    0,    0,    2},
            {-3,    0,    0,    3,    0},
            {-2,   -1,    0,    2,    2},
            {-1,    0,   -2,    3,    0},
            {-4,    0,    0,    4,    0},
            { 2,    1,   -2,    0,    1},
            { 2,   -1,    0,   -2,    2},
            { 0,    0,    1,   -1,    0},
            {-1,    2,    0,    1,    0},
            {-2,    1,    2,    0,    2},
            { 1,    1,    0,   -1,    1},
            { 1,    0,    1,   -2,    1},
            { 0,    2,    0,    0,    2},
            { 1,   -1,    2,   -3,    1},
            {-1,    1,    2,   -1,    1},
            {-2,    0,    4,   -2,    2},
            {-2,    0,    4,   -2,    1},
            {-2,   -2,    0,    2,    1},
            {-2,    0,   -2,    4,    0},
            { 1,    2,    2,   -4,    1},
            { 1,    1,    2,   -4,    2},
            {-1,    2,    2,   -2,    1},
            { 2,    0,    0,   -3,    1},
            {-1,    2,    0,    0,    1},
            { 0,    0,    0,   -2,    0},
            {-1,   -1,    2,   -2,    2},
            {-1,    1,    0,    0,    2},
            { 0,    0,    0,   -1,    2},
            {-2,    1,    0,    1,    0},
            { 1,   -2,    0,   -2,    1},
            { 1,    0,   -2,    0,    2},
            {-3,    1,    0,    2,    0},
            {-1,    1,   -2,    2,    0},
            {-1,   -1,    0,    0,    2},
            {-3,    0,    0,    2,    0},
            {-3,   -1,    0,    2,    0},
            { 2,    0,    2,   -6,    1},
            { 0,    1,    2,   -4,    2},
            { 2,    0,    0,   -4,    2},
            {-2,    1,    2,   -2,    1},
            { 0,   -1,    2,   -4,    1},
            { 0,    1,    0,   -2,    2},
            {-1,    0,    0,   -2,    0},
            { 2,    0,   -2,   -2,    1},
            {-4,    0,    2,    0,    1},
            {-1,   -1,    0,   -1,    1},
            { 0,    0,   -2,    0,    2},
            {-3,    0,    0,    1,    0},
            {-1,    0,   -2,    1,    0},
            {-2,    0,   -2,    2,    1},
            { 0,    0,   -4,    2,    0},
            {-2,   -1,   -2,    2,    0},
            { 1,    0,    2,   -6,    1},
            {-1,    0,    2,   -4,    2},
            { 1,    0,    0,   -4,    2},
            { 2,    1,    2,   -4,    2},
            { 2,    1,    2,   -4,    1},
            { 0,    1,    4,   -4,    4},
            { 0,    1,    4,   -4,    2},
            {-1,   -1,   -2,    4,    0},
            {-1,   -3,    0,    2,    0},
            {-1,    0,   -2,    4,    1},
            {-2,   -1,    0,    3,    0},
            { 0,    0,   -2,    3,    0},
            {-2,    0,    0,    3,    1},
            { 0,   -1,    0,    1,    0},
            {-3,    0,    2,    2,    0},
            { 1,    1,   -2,    2,    0},
            {-1,    1,    0,    2,    2},
            { 1,   -2,    2,   -2,    1},
            { 0,    0,    1,    0,    2},
            { 0,    0,    1,    0,    1},
            { 0,    0,    1,    0,    0},
            {-1,    2,    0,    2,    1},
            { 0,    0,    2,    0,    2},
            {-2,    0,    2,    0,    2},
            { 2,    0,    0,   -1,    1},
            { 3,    0,    0,   -2,    1},
            { 1,    0,    2,   -2,    3},
            { 1,    2,    0,    0,    1},
            { 2,    0,    2,   -3,    2},
            {-1,    1,    4,   -2,    2},
            {-2,   -2,    0,    4,    0},
            { 0,   -3,    0,    2,    0},
            { 0,    0,   -2,    4,    0},
            {-1,   -1,    0,    3,    0},
            {-2,    0,    0,    4,    2},
            {-1,    0,    0,    3,    1},
            { 2,   -2,    0,    0,    0},
            { 1,   -1,    0,    1,    0},
            {-1,    0,    0,    2,    0},
            { 0,   -2,    2,    0,    1},
            {-1,    0,    1,    2,    1},
            {-1,    1,    0,    3,    0},
            {-1,   -1,    2,    1,    2},
            { 0,   -1,    2,    0,    0},
            {-2,    1,    2,    2,    1},
            { 2,   -2,    2,   -2,    2},
            { 1,    1,    0,    1,    1},
            { 1,    0,    1,    0,    1},
            { 1,    0,    1,    0,    0},
            { 0,    2,    0,    2,    0},
            { 2,   -1,    2,   -2,    1},
            { 0,   -1,    4,   -2,    1},
            { 0,    0,    4,   -2,    3},
            { 0,    1,    4,   -2,    1},
            { 4,    0,    2,   -4,    2},
            { 2,    2,    2,   -2,    2},
            { 2,    0,    4,   -4,    2},
            {-1,   -2,    0,    4,    0},
            {-1,   -3,    2,    2,    2},
            {-3,    0,    2,    4,    2},
            {-3,    0,    2,   -2,    1},
            {-1,   -1,    0,   -2,    1},
            {-3,    0,    0,    0,    2},
            {-3,    0,   -2,    2,    0},
            { 0,    1,    0,   -4,    1},
            {-2,    1,    0,   -2,    1},
            {-4,    0,    0,    0,    1},
            {-1,    0,    0,   -4,    1},
            {-3,    0,    0,   -2,    1},
            { 0,    0,    0,    3,    2},
            {-1,    1,    0,    4,    1},
            { 1,   -2,    2,    0,    1},
            { 0,    1,    0,    3,    0},
            {-1,    0,    2,    2,    3},
            { 0,    0,    2,    2,    2},
            {-2,    0,    2,    2,    2},
            {-1,    1,    2,    2,    0},
            { 3,    0,    0,    0,    2},
            { 2,    1,    0,    1,    0},
            { 2,   -1,    2,   -1,    2},
            { 0,    0,    2,    0,    1},
            { 0,    0,    3,    0,    3},
            { 0,    0,    3,    0,    2},
            {-1,    2,    2,    2,    1},
            {-1,    0,    4,    0,    0},
            { 1,    2,    2,    0,    1},
            { 3,    1,    2,   -2,    1},
            { 1,    1,    4,   -2,    2},
            {-2,   -1,    0,    6,    0},
            { 0,   -2,    0,    4,    0},
            {-2,    0,    0,    6,    1},
            {-2,   -2,    2,    4,    2},
            { 0,   -3,    2,    2,    2},
            { 0,    0,    0,    4,    2},
            {-1,   -1,    2,    3,    2},
            {-2,    0,    2,    4,    0},
            { 2,   -1,    0,    2,    1},
            { 1,    0,    0,    3,    0},
            { 0,    1,    0,    4,    1},
            { 0,    1,    0,    4,    0},
            { 1,   -1,    2,    1,    2},
            { 0,    0,    2,    2,    3},
            { 1,    0,    2,    2,    2},
            {-1,    0,    2,    2,    2},
            {-2,    0,    4,    2,    1},
            { 2,    1,    0,    2,    1},
            { 2,    1,    0,    2,    0},
            { 2,   -1,    2,    0,    0},
            { 1,    0,    2,    1,    0},
            { 0,    1,    2,    2,    0},
            { 2,    0,    2,    0,    3},
            { 3,    0,    2,    0,    2},
            { 1,    0,    2,    0,    2},
            { 1,    0,    3,    0,    3},
            { 1,    1,    2,    1,    1},
            { 0,    2,    2,    2,    2},
            { 2,    1,    2,    0,    0},
            { 2,    0,    4,   -2,    1},
            { 4,    1,    2,   -2,    2},
            {-1,   -1,    0,    6,    0},
            {-3,   -1,    2,    6,    2},
            {-1,    0,    0,    6,    1},
            {-3,    0,    2,    6,    1},
            { 1,   -1,    0,    4,    1},
            { 1,   -1,    0,    4,    0},
            {-2,    0,    2,    5,    2},
            { 1,   -2,    2,    2,    1},
            { 3,   -1,    0,    2,    0},
            { 1,   -1,    2,    2,    0},
            { 0,    0,    2,    3,    1},
            {-1,    1,    2,    4,    1},
            { 0,    1,    2,    3,    2},
            {-1,    0,    4,    2,    1},
            { 2,    0,    2,    1,    1},
            { 5,    0,    0,    0,    0},
            { 2,    1,    2,    1,    2},
            { 1,    0,    4,    0,    1},
            { 3,    1,    2,    0,    1},
            { 3,    0,    4,   -2,    2},
            {-2,   -1,    2,    6,    2},
            { 0,    0,    0,    6,    0},
            { 0,   -2,    2,    4,    2},
            {-2,    0,    2,    6,    1},
            { 2,    0,    0,    4,    1},
            { 2,    0,    0,    4,    0},
            { 2,   -2,    2,    2,    2},
            { 0,    0,    2,    4,    0},
            { 1,    0,    2,    3,    2},
            { 4,    0,    0,    2,    0},
            { 2,    0,    2,    2,    0},
            { 0,    0,    4,    2,    2},
            { 4,   -1,    2,    0,    2},
            { 3,    0,    2,    1,    2},
            { 2,    1,    2,    2,    1},
            { 4,    1,    2,    0,    2},
            {-1,   -1,    2,    6,    2},
            {-1,    0,    2,    6,    1},
            { 1,   -1,    2,    4,    1},
            { 1,    1,    2,    4,    2},
            { 3,    1,    2,    2,    2},
            { 5,    0,    2,    0,    1},
            { 2,   -1,    2,    4,    2},
            { 2,    0,    2,    4,    1}};
        
        /*
         Luni-Solar nutation coefficients, unit 1e-7 arcsec:
         longitude (sin, t*sin, cos), obliquity (cos, t*cos, sin)
         
         Each row of coefficients in 'cls_t' belongs with the corresponding
         row of fundamental-argument multipliers in 'nals_t'.
         */
        
        static const double cls_t[678][6] = {
            {-172064161.0, -174666.0,  33386.0, 92052331.0,  9086.0, 15377.0},
            { -13170906.0,   -1675.0, -13696.0,  5730336.0, -3015.0, -4587.0},
            {  -2276413.0,    -234.0,   2796.0,   978459.0,  -485.0,  1374.0},
            {   2074554.0,     207.0,   -698.0,  -897492.0,   470.0,  -291.0},
            {   1475877.0,   -3633.0,  11817.0,    73871.0,  -184.0, -1924.0},
            {   -516821.0,    1226.0,   -524.0,   224386.0,  -677.0,  -174.0},
            {    711159.0,      73.0,   -872.0,    -6750.0,     0.0,   358.0},
            {   -387298.0,    -367.0,    380.0,   200728.0,    18.0,   318.0},
            {   -301461.0,     -36.0,    816.0,   129025.0,   -63.0,   367.0},
            {    215829.0,    -494.0,    111.0,   -95929.0,   299.0,   132.0},
            {    128227.0,     137.0,    181.0,   -68982.0,    -9.0,    39.0},
            {    123457.0,      11.0,     19.0,   -53311.0,    32.0,    -4.0},
            {    156994.0,      10.0,   -168.0,    -1235.0,     0.0,    82.0},
            {     63110.0,      63.0,     27.0,   -33228.0,     0.0,    -9.0},
            {    -57976.0,     -63.0,   -189.0,    31429.0,     0.0,   -75.0},
            {    -59641.0,     -11.0,    149.0,    25543.0,   -11.0,    66.0},
            {    -51613.0,     -42.0,    129.0,    26366.0,     0.0,    78.0},
            {     45893.0,      50.0,     31.0,   -24236.0,   -10.0,    20.0},
            {     63384.0,      11.0,   -150.0,    -1220.0,     0.0,    29.0},
            {    -38571.0,      -1.0,    158.0,    16452.0,   -11.0,    68.0},
            {     32481.0,       0.0,      0.0,   -13870.0,     0.0,     0.0},
            {    -47722.0,       0.0,    -18.0,      477.0,     0.0,   -25.0},
            {    -31046.0,      -1.0,    131.0,    13238.0,   -11.0,    59.0},
            {     28593.0,       0.0,     -1.0,   -12338.0,    10.0,    -3.0},
            {     20441.0,      21.0,     10.0,   -10758.0,     0.0,    -3.0},
            {     29243.0,       0.0,    -74.0,     -609.0,     0.0,    13.0},
            {     25887.0,       0.0,    -66.0,     -550.0,     0.0,    11.0},
            {    -14053.0,     -25.0,     79.0,     8551.0,    -2.0,   -45.0},
            {     15164.0,      10.0,     11.0,    -8001.0,     0.0,    -1.0},
            {    -15794.0,      72.0,    -16.0,     6850.0,   -42.0,    -5.0},
            {     21783.0,       0.0,     13.0,     -167.0,     0.0,    13.0},
            {    -12873.0,     -10.0,    -37.0,     6953.0,     0.0,   -14.0},
            {    -12654.0,      11.0,     63.0,     6415.0,     0.0,    26.0},
            {    -10204.0,       0.0,     25.0,     5222.0,     0.0,    15.0},
            {     16707.0,     -85.0,    -10.0,      168.0,    -1.0,    10.0},
            {     -7691.0,       0.0,     44.0,     3268.0,     0.0,    19.0},
            {    -11024.0,       0.0,    -14.0,      104.0,     0.0,     2.0},
            {      7566.0,     -21.0,    -11.0,    -3250.0,     0.0,    -5.0},
            {     -6637.0,     -11.0,     25.0,     3353.0,     0.0,    14.0},
            {     -7141.0,      21.0,      8.0,     3070.0,     0.0,     4.0},
            {     -6302.0,     -11.0,      2.0,     3272.0,     0.0,     4.0},
            {      5800.0,      10.0,      2.0,    -3045.0,     0.0,    -1.0},
            {      6443.0,       0.0,     -7.0,    -2768.0,     0.0,    -4.0},
            {     -5774.0,     -11.0,    -15.0,     3041.0,     0.0,    -5.0},
            {     -5350.0,       0.0,     21.0,     2695.0,     0.0,    12.0},
            {     -4752.0,     -11.0,     -3.0,     2719.0,     0.0,    -3.0},
            {     -4940.0,     -11.0,    -21.0,     2720.0,     0.0,    -9.0},
            {      7350.0,       0.0,     -8.0,      -51.0,     0.0,     4.0},
            {      4065.0,       0.0,      6.0,    -2206.0,     0.0,     1.0},
            {      6579.0,       0.0,    -24.0,     -199.0,     0.0,     2.0},
            {      3579.0,       0.0,      5.0,    -1900.0,     0.0,     1.0},
            {      4725.0,       0.0,     -6.0,      -41.0,     0.0,     3.0},
            {     -3075.0,       0.0,     -2.0,     1313.0,     0.0,    -1.0},
            {     -2904.0,       0.0,     15.0,     1233.0,     0.0,     7.0},
            {      4348.0,       0.0,    -10.0,      -81.0,     0.0,     2.0},
            {     -2878.0,       0.0,      8.0,     1232.0,     0.0,     4.0},
            {     -4230.0,       0.0,      5.0,      -20.0,     0.0,    -2.0},
            {     -2819.0,       0.0,      7.0,     1207.0,     0.0,     3.0},
            {     -4056.0,       0.0,      5.0,       40.0,     0.0,    -2.0},
            {     -2647.0,       0.0,     11.0,     1129.0,     0.0,     5.0},
            {     -2294.0,       0.0,    -10.0,     1266.0,     0.0,    -4.0},
            {      2481.0,       0.0,     -7.0,    -1062.0,     0.0,    -3.0},
            {      2179.0,       0.0,     -2.0,    -1129.0,     0.0,    -2.0},
            {      3276.0,       0.0,      1.0,       -9.0,     0.0,     0.0},
            {     -3389.0,       0.0,      5.0,       35.0,     0.0,    -2.0},
            {      3339.0,       0.0,    -13.0,     -107.0,     0.0,     1.0},
            {     -1987.0,       0.0,     -6.0,     1073.0,     0.0,    -2.0},
            {     -1981.0,       0.0,      0.0,      854.0,     0.0,     0.0},
            {      4026.0,       0.0,   -353.0,     -553.0,     0.0,  -139.0},
            {      1660.0,       0.0,     -5.0,     -710.0,     0.0,    -2.0},
            {     -1521.0,       0.0,      9.0,      647.0,     0.0,     4.0},
            {      1314.0,       0.0,      0.0,     -700.0,     0.0,     0.0},
            {     -1283.0,       0.0,      0.0,      672.0,     0.0,     0.0},
            {     -1331.0,       0.0,      8.0,      663.0,     0.0,     4.0},
            {      1383.0,       0.0,     -2.0,     -594.0,     0.0,    -2.0},
            {      1405.0,       0.0,      4.0,     -610.0,     0.0,     2.0},
            {      1290.0,       0.0,      0.0,     -556.0,     0.0,     0.0},
            {     -1214.0,       0.0,      5.0,      518.0,     0.0,     2.0},
            {      1146.0,       0.0,     -3.0,     -490.0,     0.0,    -1.0},
            {      1019.0,       0.0,     -1.0,     -527.0,     0.0,    -1.0},
            {     -1100.0,       0.0,      9.0,      465.0,     0.0,     4.0},
            {      -970.0,       0.0,      2.0,      496.0,     0.0,     1.0},
            {      1575.0,       0.0,     -6.0,      -50.0,     0.0,     0.0},
            {       934.0,       0.0,     -3.0,     -399.0,     0.0,    -1.0},
            {       922.0,       0.0,     -1.0,     -395.0,     0.0,    -1.0},
            {       815.0,       0.0,     -1.0,     -422.0,     0.0,    -1.0},
            {       834.0,       0.0,      2.0,     -440.0,     0.0,     1.0},
            {      1248.0,       0.0,      0.0,     -170.0,     0.0,     1.0},
            {      1338.0,       0.0,     -5.0,      -39.0,     0.0,     0.0},
            {       716.0,       0.0,     -2.0,     -389.0,     0.0,    -1.0},
            {      1282.0,       0.0,     -3.0,      -23.0,     0.0,     1.0},
            {       742.0,       0.0,      1.0,     -391.0,     0.0,     0.0},
            {      1020.0,       0.0,    -25.0,     -495.0,     0.0,   -10.0},
            {       715.0,       0.0,     -4.0,     -326.0,     0.0,     2.0},
            {      -666.0,       0.0,     -3.0,      369.0,     0.0,    -1.0},
            {      -667.0,       0.0,      1.0,      346.0,     0.0,     1.0},
            {      -704.0,       0.0,      0.0,      304.0,     0.0,     0.0},
            {      -694.0,       0.0,      5.0,      294.0,     0.0,     2.0},
            {     -1014.0,       0.0,     -1.0,        4.0,     0.0,    -1.0},
            {      -585.0,       0.0,     -2.0,      316.0,     0.0,    -1.0},
            {      -949.0,       0.0,      1.0,        8.0,     0.0,    -1.0},
            {      -595.0,       0.0,      0.0,      258.0,     0.0,     0.0},
            {       528.0,       0.0,      0.0,     -279.0,     0.0,     0.0},
            {      -590.0,       0.0,      4.0,      252.0,     0.0,     2.0},
            {       570.0,       0.0,     -2.0,     -244.0,     0.0,    -1.0},
            {      -502.0,       0.0,      3.0,      250.0,     0.0,     2.0},
            {      -875.0,       0.0,      1.0,       29.0,     0.0,     0.0},
            {      -492.0,       0.0,     -3.0,      275.0,     0.0,    -1.0},
            {       535.0,       0.0,     -2.0,     -228.0,     0.0,    -1.0},
            {      -467.0,       0.0,      1.0,      240.0,     0.0,     1.0},
            {       591.0,       0.0,      0.0,     -253.0,     0.0,     0.0},
            {      -453.0,       0.0,     -1.0,      244.0,     0.0,    -1.0},
            {       766.0,       0.0,      1.0,        9.0,     0.0,     0.0},
            {      -446.0,       0.0,      2.0,      225.0,     0.0,     1.0},
            {      -488.0,       0.0,      2.0,      207.0,     0.0,     1.0},
            {      -468.0,       0.0,      0.0,      201.0,     0.0,     0.0},
            {      -421.0,       0.0,      1.0,      216.0,     0.0,     1.0},
            {       463.0,       0.0,      0.0,     -200.0,     0.0,     0.0},
            {      -673.0,       0.0,      2.0,       14.0,     0.0,     0.0},
            {       658.0,       0.0,      0.0,       -2.0,     0.0,     0.0},
            {      -438.0,       0.0,      0.0,      188.0,     0.0,     0.0},
            {      -390.0,       0.0,      0.0,      205.0,     0.0,     0.0},
            {       639.0,     -11.0,     -2.0,      -19.0,     0.0,     0.0},
            {       412.0,       0.0,     -2.0,     -176.0,     0.0,    -1.0},
            {      -361.0,       0.0,      0.0,      189.0,     0.0,     0.0},
            {       360.0,       0.0,     -1.0,     -185.0,     0.0,    -1.0},
            {       588.0,       0.0,     -3.0,      -24.0,     0.0,     0.0},
            {      -578.0,       0.0,      1.0,        5.0,     0.0,     0.0},
            {      -396.0,       0.0,      0.0,      171.0,     0.0,     0.0},
            {       565.0,       0.0,     -1.0,       -6.0,     0.0,     0.0},
            {      -335.0,       0.0,     -1.0,      184.0,     0.0,    -1.0},
            {       357.0,       0.0,      1.0,     -154.0,     0.0,     0.0},
            {       321.0,       0.0,      1.0,     -174.0,     0.0,     0.0},
            {      -301.0,       0.0,     -1.0,      162.0,     0.0,     0.0},
            {      -334.0,       0.0,      0.0,      144.0,     0.0,     0.0},
            {       493.0,       0.0,     -2.0,      -15.0,     0.0,     0.0},
            {       494.0,       0.0,     -2.0,      -19.0,     0.0,     0.0},
            {       337.0,       0.0,     -1.0,     -143.0,     0.0,    -1.0},
            {       280.0,       0.0,     -1.0,     -144.0,     0.0,     0.0},
            {       309.0,       0.0,      1.0,     -134.0,     0.0,     0.0},
            {      -263.0,       0.0,      2.0,      131.0,     0.0,     1.0},
            {       253.0,       0.0,      1.0,     -138.0,     0.0,     0.0},
            {       245.0,       0.0,      0.0,     -128.0,     0.0,     0.0},
            {       416.0,       0.0,     -2.0,      -17.0,     0.0,     0.0},
            {      -229.0,       0.0,      0.0,      128.0,     0.0,     0.0},
            {       231.0,       0.0,      0.0,     -120.0,     0.0,     0.0},
            {      -259.0,       0.0,      2.0,      109.0,     0.0,     1.0},
            {       375.0,       0.0,     -1.0,       -8.0,     0.0,     0.0},
            {       252.0,       0.0,      0.0,     -108.0,     0.0,     0.0},
            {      -245.0,       0.0,      1.0,      104.0,     0.0,     0.0},
            {       243.0,       0.0,     -1.0,     -104.0,     0.0,     0.0},
            {       208.0,       0.0,      1.0,     -112.0,     0.0,     0.0},
            {       199.0,       0.0,      0.0,     -102.0,     0.0,     0.0},
            {      -208.0,       0.0,      1.0,      105.0,     0.0,     0.0},
            {       335.0,       0.0,     -2.0,      -14.0,     0.0,     0.0},
            {      -325.0,       0.0,      1.0,        7.0,     0.0,     0.0},
            {      -187.0,       0.0,      0.0,       96.0,     0.0,     0.0},
            {       197.0,       0.0,     -1.0,     -100.0,     0.0,     0.0},
            {      -192.0,       0.0,      2.0,       94.0,     0.0,     1.0},
            {      -188.0,       0.0,      0.0,       83.0,     0.0,     0.0},
            {       276.0,       0.0,      0.0,       -2.0,     0.0,     0.0},
            {      -286.0,       0.0,      1.0,        6.0,     0.0,     0.0},
            {       186.0,       0.0,     -1.0,      -79.0,     0.0,     0.0},
            {      -219.0,       0.0,      0.0,       43.0,     0.0,     0.0},
            {       276.0,       0.0,      0.0,        2.0,     0.0,     0.0},
            {      -153.0,       0.0,     -1.0,       84.0,     0.0,     0.0},
            {      -156.0,       0.0,      0.0,       81.0,     0.0,     0.0},
            {      -154.0,       0.0,      1.0,       78.0,     0.0,     0.0},
            {      -174.0,       0.0,      1.0,       75.0,     0.0,     0.0},
            {      -163.0,       0.0,      2.0,       69.0,     0.0,     1.0},
            {      -228.0,       0.0,      0.0,        1.0,     0.0,     0.0},
            {        91.0,       0.0,     -4.0,      -54.0,     0.0,    -2.0},
            {       175.0,       0.0,      0.0,      -75.0,     0.0,     0.0},
            {      -159.0,       0.0,      0.0,       69.0,     0.0,     0.0},
            {       141.0,       0.0,      0.0,      -72.0,     0.0,     0.0},
            {       147.0,       0.0,      0.0,      -75.0,     0.0,     0.0},
            {      -132.0,       0.0,      0.0,       69.0,     0.0,     0.0},
            {       159.0,       0.0,    -28.0,      -54.0,     0.0,    11.0},
            {       213.0,       0.0,      0.0,       -4.0,     0.0,     0.0},
            {       123.0,       0.0,      0.0,      -64.0,     0.0,     0.0},
            {      -118.0,       0.0,     -1.0,       66.0,     0.0,     0.0},
            {       144.0,       0.0,     -1.0,      -61.0,     0.0,     0.0},
            {      -121.0,       0.0,      1.0,       60.0,     0.0,     0.0},
            {      -134.0,       0.0,      1.0,       56.0,     0.0,     1.0},
            {      -105.0,       0.0,      0.0,       57.0,     0.0,     0.0},
            {      -102.0,       0.0,      0.0,       56.0,     0.0,     0.0},
            {       120.0,       0.0,      0.0,      -52.0,     0.0,     0.0},
            {       101.0,       0.0,      0.0,      -54.0,     0.0,     0.0},
            {      -113.0,       0.0,      0.0,       59.0,     0.0,     0.0},
            {      -106.0,       0.0,      0.0,       61.0,     0.0,     0.0},
            {      -129.0,       0.0,      1.0,       55.0,     0.0,     0.0},
            {      -114.0,       0.0,      0.0,       57.0,     0.0,     0.0},
            {       113.0,       0.0,     -1.0,      -49.0,     0.0,     0.0},
            {      -102.0,       0.0,      0.0,       44.0,     0.0,     0.0},
            {       -94.0,       0.0,      0.0,       51.0,     0.0,     0.0},
            {      -100.0,       0.0,     -1.0,       56.0,     0.0,     0.0},
            {        87.0,       0.0,      0.0,      -47.0,     0.0,     0.0},
            {       161.0,       0.0,      0.0,       -1.0,     0.0,     0.0},
            {        96.0,       0.0,      0.0,      -50.0,     0.0,     0.0},
            {       151.0,       0.0,     -1.0,       -5.0,     0.0,     0.0},
            {      -104.0,       0.0,      0.0,       44.0,     0.0,     0.0},
            {      -110.0,       0.0,      0.0,       48.0,     0.0,     0.0},
            {      -100.0,       0.0,      1.0,       50.0,     0.0,     0.0},
            {        92.0,       0.0,     -5.0,       12.0,     0.0,    -2.0},
            {        82.0,       0.0,      0.0,      -45.0,     0.0,     0.0},
            {        82.0,       0.0,      0.0,      -45.0,     0.0,     0.0},
            {       -78.0,       0.0,      0.0,       41.0,     0.0,     0.0},
            {       -77.0,       0.0,      0.0,       43.0,     0.0,     0.0},
            {         2.0,       0.0,      0.0,       54.0,     0.0,     0.0},
            {        94.0,       0.0,      0.0,      -40.0,     0.0,     0.0},
            {       -93.0,       0.0,      0.0,       40.0,     0.0,     0.0},
            {       -83.0,       0.0,     10.0,       40.0,     0.0,    -2.0},
            {        83.0,       0.0,      0.0,      -36.0,     0.0,     0.0},
            {       -91.0,       0.0,      0.0,       39.0,     0.0,     0.0},
            {       128.0,       0.0,      0.0,       -1.0,     0.0,     0.0},
            {       -79.0,       0.0,      0.0,       34.0,     0.0,     0.0},
            {       -83.0,       0.0,      0.0,       47.0,     0.0,     0.0},
            {        84.0,       0.0,      0.0,      -44.0,     0.0,     0.0},
            {        83.0,       0.0,      0.0,      -43.0,     0.0,     0.0},
            {        91.0,       0.0,      0.0,      -39.0,     0.0,     0.0},
            {       -77.0,       0.0,      0.0,       39.0,     0.0,     0.0},
            {        84.0,       0.0,      0.0,      -43.0,     0.0,     0.0},
            {       -92.0,       0.0,      1.0,       39.0,     0.0,     0.0},
            {       -92.0,       0.0,      1.0,       39.0,     0.0,     0.0},
            {       -94.0,       0.0,      0.0,        0.0,     0.0,     0.0},
            {        68.0,       0.0,      0.0,      -36.0,     0.0,     0.0},
            {       -61.0,       0.0,      0.0,       32.0,     0.0,     0.0},
            {        71.0,       0.0,      0.0,      -31.0,     0.0,     0.0},
            {        62.0,       0.0,      0.0,      -34.0,     0.0,     0.0},
            {       -63.0,       0.0,      0.0,       33.0,     0.0,     0.0},
            {       -73.0,       0.0,      0.0,       32.0,     0.0,     0.0},
            {       115.0,       0.0,      0.0,       -2.0,     0.0,     0.0},
            {      -103.0,       0.0,      0.0,        2.0,     0.0,     0.0},
            {        63.0,       0.0,      0.0,      -28.0,     0.0,     0.0},
            {        74.0,       0.0,      0.0,      -32.0,     0.0,     0.0},
            {      -103.0,       0.0,     -3.0,        3.0,     0.0,    -1.0},
            {       -69.0,       0.0,      0.0,       30.0,     0.0,     0.0},
            {        57.0,       0.0,      0.0,      -29.0,     0.0,     0.0},
            {        94.0,       0.0,      0.0,       -4.0,     0.0,     0.0},
            {        64.0,       0.0,      0.0,      -33.0,     0.0,     0.0},
            {       -63.0,       0.0,      0.0,       26.0,     0.0,     0.0},
            {       -38.0,       0.0,      0.0,       20.0,     0.0,     0.0},
            {       -43.0,       0.0,      0.0,       24.0,     0.0,     0.0},
            {       -45.0,       0.0,      0.0,       23.0,     0.0,     0.0},
            {        47.0,       0.0,      0.0,      -24.0,     0.0,     0.0},
            {       -48.0,       0.0,      0.0,       25.0,     0.0,     0.0},
            {        45.0,       0.0,      0.0,      -26.0,     0.0,     0.0},
            {        56.0,       0.0,      0.0,      -25.0,     0.0,     0.0},
            {        88.0,       0.0,      0.0,        2.0,     0.0,     0.0},
            {       -75.0,       0.0,      0.0,        0.0,     0.0,     0.0},
            {        85.0,       0.0,      0.0,        0.0,     0.0,     0.0},
            {        49.0,       0.0,      0.0,      -26.0,     0.0,     0.0},
            {       -74.0,       0.0,     -3.0,       -1.0,     0.0,    -1.0},
            {       -39.0,       0.0,      0.0,       21.0,     0.0,     0.0},
            {        45.0,       0.0,      0.0,      -20.0,     0.0,     0.0},
            {        51.0,       0.0,      0.0,      -22.0,     0.0,     0.0},
            {       -40.0,       0.0,      0.0,       21.0,     0.0,     0.0},
            {        41.0,       0.0,      0.0,      -21.0,     0.0,     0.0},
            {       -42.0,       0.0,      0.0,       24.0,     0.0,     0.0},
            {       -51.0,       0.0,      0.0,       22.0,     0.0,     0.0},
            {       -42.0,       0.0,      0.0,       22.0,     0.0,     0.0},
            {        39.0,       0.0,      0.0,      -21.0,     0.0,     0.0},
            {        46.0,       0.0,      0.0,      -18.0,     0.0,     0.0},
            {       -53.0,       0.0,      0.0,       22.0,     0.0,     0.0},
            {        82.0,       0.0,      0.0,       -4.0,     0.0,     0.0},
            {        81.0,       0.0,     -1.0,       -4.0,     0.0,     0.0},
            {        47.0,       0.0,      0.0,      -19.0,     0.0,     0.0},
            {        53.0,       0.0,      0.0,      -23.0,     0.0,     0.0},
            {       -45.0,       0.0,      0.0,       22.0,     0.0,     0.0},
            {       -44.0,       0.0,      0.0,       -2.0,     0.0,     0.0},
            {       -33.0,       0.0,      0.0,       16.0,     0.0,     0.0},
            {       -61.0,       0.0,      0.0,        1.0,     0.0,     0.0},
            {        28.0,       0.0,      0.0,      -15.0,     0.0,     0.0},
            {       -38.0,       0.0,      0.0,       19.0,     0.0,     0.0},
            {       -33.0,       0.0,      0.0,       21.0,     0.0,     0.0},
            {       -60.0,       0.0,      0.0,        0.0,     0.0,     0.0},
            {        48.0,       0.0,      0.0,      -10.0,     0.0,     0.0},
            {        27.0,       0.0,      0.0,      -14.0,     0.0,     0.0},
            {        38.0,       0.0,      0.0,      -20.0,     0.0,     0.0},
            {        31.0,       0.0,      0.0,      -13.0,     0.0,     0.0},
            {       -29.0,       0.0,      0.0,       15.0,     0.0,     0.0},
            {        28.0,       0.0,      0.0,      -15.0,     0.0,     0.0},
            {       -32.0,       0.0,      0.0,       15.0,     0.0,     0.0},
            {        45.0,       0.0,      0.0,       -8.0,     0.0,     0.0},
            {       -44.0,       0.0,      0.0,       19.0,     0.0,     0.0},
            {        28.0,       0.0,      0.0,      -15.0,     0.0,     0.0},
            {       -51.0,       0.0,      0.0,        0.0,     0.0,     0.0},
            {       -36.0,       0.0,      0.0,       20.0,     0.0,     0.0},
            {        44.0,       0.0,      0.0,      -19.0,     0.0,     0.0},
            {        26.0,       0.0,      0.0,      -14.0,     0.0,     0.0},
            {       -60.0,       0.0,      0.0,        2.0,     0.0,     0.0},
            {        35.0,       0.0,      0.0,      -18.0,     0.0,     0.0},
            {       -27.0,       0.0,      0.0,       11.0,     0.0,     0.0},
            {        47.0,       0.0,      0.0,       -1.0,     0.0,     0.0},
            {        36.0,       0.0,      0.0,      -15.0,     0.0,     0.0},
            {       -36.0,       0.0,      0.0,       20.0,     0.0,     0.0},
            {       -35.0,       0.0,      0.0,       19.0,     0.0,     0.0},
            {       -37.0,       0.0,      0.0,       19.0,     0.0,     0.0},
            {        32.0,       0.0,      0.0,      -16.0,     0.0,     0.0},
            {        35.0,       0.0,      0.0,      -14.0,     0.0,     0.0},
            {        32.0,       0.0,      0.0,      -13.0,     0.0,     0.0},
            {        65.0,       0.0,      0.0,       -2.0,     0.0,     0.0},
            {        47.0,       0.0,      0.0,       -1.0,     0.0,     0.0},
            {        32.0,       0.0,      0.0,      -16.0,     0.0,     0.0},
            {        37.0,       0.0,      0.0,      -16.0,     0.0,     0.0},
            {       -30.0,       0.0,      0.0,       15.0,     0.0,     0.0},
            {       -32.0,       0.0,      0.0,       16.0,     0.0,     0.0},
            {       -31.0,       0.0,      0.0,       13.0,     0.0,     0.0},
            {        37.0,       0.0,      0.0,      -16.0,     0.0,     0.0},
            {        31.0,       0.0,      0.0,      -13.0,     0.0,     0.0},
            {        49.0,       0.0,      0.0,       -2.0,     0.0,     0.0},
            {        32.0,       0.0,      0.0,      -13.0,     0.0,     0.0},
            {        23.0,       0.0,      0.0,      -12.0,     0.0,     0.0},
            {       -43.0,       0.0,      0.0,       18.0,     0.0,     0.0},
            {        26.0,       0.0,      0.0,      -11.0,     0.0,     0.0},
            {       -32.0,       0.0,      0.0,       14.0,     0.0,     0.0},
            {       -29.0,       0.0,      0.0,       14.0,     0.0,     0.0},
            {       -27.0,       0.0,      0.0,       12.0,     0.0,     0.0},
            {        30.0,       0.0,      0.0,        0.0,     0.0,     0.0},
            {       -11.0,       0.0,      0.0,        5.0,     0.0,     0.0},
            {       -21.0,       0.0,      0.0,       10.0,     0.0,     0.0},
            {       -34.0,       0.0,      0.0,       15.0,     0.0,     0.0},
            {       -10.0,       0.0,      0.0,        6.0,     0.0,     0.0},
            {       -36.0,       0.0,      0.0,        0.0,     0.0,     0.0},
            {        -9.0,       0.0,      0.0,        4.0,     0.0,     0.0},
            {       -12.0,       0.0,      0.0,        5.0,     0.0,     0.0},
            {       -21.0,       0.0,      0.0,        5.0,     0.0,     0.0},
            {       -29.0,       0.0,      0.0,       -1.0,     0.0,     0.0},
            {       -15.0,       0.0,      0.0,        3.0,     0.0,     0.0},
            {       -20.0,       0.0,      0.0,        0.0,     0.0,     0.0},
            {        28.0,       0.0,      0.0,        0.0,     0.0,    -2.0},
            {        17.0,       0.0,      0.0,        0.0,     0.0,     0.0},
            {       -22.0,       0.0,      0.0,       12.0,     0.0,     0.0},
            {       -14.0,       0.0,      0.0,        7.0,     0.0,     0.0},
            {        24.0,       0.0,      0.0,      -11.0,     0.0,     0.0},
            {        11.0,       0.0,      0.0,       -6.0,     0.0,     0.0},
            {        14.0,       0.0,      0.0,       -6.0,     0.0,     0.0},
            {        24.0,       0.0,      0.0,        0.0,     0.0,     0.0},
            {        18.0,       0.0,      0.0,       -8.0,     0.0,     0.0},
            {       -38.0,       0.0,      0.0,        0.0,     0.0,     0.0},
            {       -31.0,       0.0,      0.0,        0.0,     0.0,     0.0},
            {       -16.0,       0.0,      0.0,        8.0,     0.0,     0.0},
            {        29.0,       0.0,      0.0,        0.0,     0.0,     0.0},
            {       -18.0,       0.0,      0.0,       10.0,     0.0,     0.0},
            {       -10.0,       0.0,      0.0,        5.0,     0.0,     0.0},
            {       -17.0,       0.0,      0.0,       10.0,     0.0,     0.0},
            {         9.0,       0.0,      0.0,       -4.0,     0.0,     0.0},
            {        16.0,       0.0,      0.0,       -6.0,     0.0,     0.0},
            {        22.0,       0.0,      0.0,      -12.0,     0.0,     0.0},
            {        20.0,       0.0,      0.0,        0.0,     0.0,     0.0},
            {       -13.0,       0.0,      0.0,        6.0,     0.0,     0.0},
            {       -17.0,       0.0,      0.0,        9.0,     0.0,     0.0},
            {       -14.0,       0.0,      0.0,        8.0,     0.0,     0.0},
            {         0.0,       0.0,      0.0,       -7.0,     0.0,     0.0},
            {        14.0,       0.0,      0.0,        0.0,     0.0,     0.0},
            {        19.0,       0.0,      0.0,      -10.0,     0.0,     0.0},
            {       -34.0,       0.0,      0.0,        0.0,     0.0,     0.0},
            {       -20.0,       0.0,      0.0,        8.0,     0.0,     0.0},
            {         9.0,       0.0,      0.0,       -5.0,     0.0,     0.0},
            {       -18.0,       0.0,      0.0,        7.0,     0.0,     0.0},
            {        13.0,       0.0,      0.0,       -6.0,     0.0,     0.0},
            {        17.0,       0.0,      0.0,        0.0,     0.0,     0.0},
            {       -12.0,       0.0,      0.0,        5.0,     0.0,     0.0},
            {        15.0,       0.0,      0.0,       -8.0,     0.0,     0.0},
            {       -11.0,       0.0,      0.0,        3.0,     0.0,     0.0},
            {        13.0,       0.0,      0.0,       -5.0,     0.0,     0.0},
            {       -18.0,       0.0,      0.0,        0.0,     0.0,     0.0},
            {       -35.0,       0.0,      0.0,        0.0,     0.0,     0.0},
            {         9.0,       0.0,      0.0,       -4.0,     0.0,     0.0},
            {       -19.0,       0.0,      0.0,       10.0,     0.0,     0.0},
            {       -26.0,       0.0,      0.0,       11.0,     0.0,     0.0},
            {         8.0,       0.0,      0.0,       -4.0,     0.0,     0.0},
            {       -10.0,       0.0,      0.0,        4.0,     0.0,     0.0},
            {        10.0,       0.0,      0.0,       -6.0,     0.0,     0.0},
            {       -21.0,       0.0,      0.0,        9.0,     0.0,     0.0},
            {       -15.0,       0.0,      0.0,        0.0,     0.0,     0.0},
            {         9.0,       0.0,      0.0,       -5.0,     0.0,     0.0},
            {       -29.0,       0.0,      0.0,        0.0,     0.0,     0.0},
            {       -19.0,       0.0,      0.0,       10.0,     0.0,     0.0},
            {        12.0,       0.0,      0.0,       -5.0,     0.0,     0.0},
            {        22.0,       0.0,      0.0,       -9.0,     0.0,     0.0},
            {       -10.0,       0.0,      0.0,        5.0,     0.0,     0.0},
            {       -20.0,       0.0,      0.0,       11.0,     0.0,     0.0},
            {       -20.0,       0.0,      0.0,        0.0,     0.0,     0.0},
            {       -17.0,       0.0,      0.0,        7.0,     0.0,     0.0},
            {        15.0,       0.0,      0.0,       -3.0,     0.0,     0.0},
            {         8.0,       0.0,      0.0,       -4.0,     0.0,     0.0},
            {        14.0,       0.0,      0.0,        0.0,     0.0,     0.0},
            {       -12.0,       0.0,      0.0,        6.0,     0.0,     0.0},
            {        25.0,       0.0,      0.0,        0.0,     0.0,     0.0},
            {       -13.0,       0.0,      0.0,        6.0,     0.0,     0.0},
            {       -14.0,       0.0,      0.0,        8.0,     0.0,     0.0},
            {        13.0,       0.0,      0.0,       -5.0,     0.0,     0.0},
            {       -17.0,       0.0,      0.0,        9.0,     0.0,     0.0},
            {       -12.0,       0.0,      0.0,        6.0,     0.0,     0.0},
            {       -10.0,       0.0,      0.0,        5.0,     0.0,     0.0},
            {        10.0,       0.0,      0.0,       -6.0,     0.0,     0.0},
            {       -15.0,       0.0,      0.0,        0.0,     0.0,     0.0},
            {       -22.0,       0.0,      0.0,        0.0,     0.0,     0.0},
            {        28.0,       0.0,      0.0,       -1.0,     0.0,     0.0},
            {        15.0,       0.0,      0.0,       -7.0,     0.0,     0.0},
            {        23.0,       0.0,      0.0,      -10.0,     0.0,     0.0},
            {        12.0,       0.0,      0.0,       -5.0,     0.0,     0.0},
            {        29.0,       0.0,      0.0,       -1.0,     0.0,     0.0},
            {       -25.0,       0.0,      0.0,        1.0,     0.0,     0.0},
            {        22.0,       0.0,      0.0,        0.0,     0.0,     0.0},
            {       -18.0,       0.0,      0.0,        0.0,     0.0,     0.0},
            {        15.0,       0.0,      0.0,        3.0,     0.0,     0.0},
            {       -23.0,       0.0,      0.0,        0.0,     0.0,     0.0},
            {        12.0,       0.0,      0.0,       -5.0,     0.0,     0.0},
            {        -8.0,       0.0,      0.0,        4.0,     0.0,     0.0},
            {       -19.0,       0.0,      0.0,        0.0,     0.0,     0.0},
            {       -10.0,       0.0,      0.0,        4.0,     0.0,     0.0},
            {        21.0,       0.0,      0.0,       -9.0,     0.0,     0.0},
            {        23.0,       0.0,      0.0,       -1.0,     0.0,     0.0},
            {       -16.0,       0.0,      0.0,        8.0,     0.0,     0.0},
            {       -19.0,       0.0,      0.0,        9.0,     0.0,     0.0},
            {       -22.0,       0.0,      0.0,       10.0,     0.0,     0.0},
            {        27.0,       0.0,      0.0,       -1.0,     0.0,     0.0},
            {        16.0,       0.0,      0.0,       -8.0,     0.0,     0.0},
            {        19.0,       0.0,      0.0,       -8.0,     0.0,     0.0},
            {         9.0,       0.0,      0.0,       -4.0,     0.0,     0.0},
            {        -9.0,       0.0,      0.0,        4.0,     0.0,     0.0},
            {        -9.0,       0.0,      0.0,        4.0,     0.0,     0.0},
            {        -8.0,       0.0,      0.0,        4.0,     0.0,     0.0},
            {        18.0,       0.0,      0.0,       -9.0,     0.0,     0.0},
            {        16.0,       0.0,      0.0,       -1.0,     0.0,     0.0},
            {       -10.0,       0.0,      0.0,        4.0,     0.0,     0.0},
            {       -23.0,       0.0,      0.0,        9.0,     0.0,     0.0},
            {        16.0,       0.0,      0.0,       -1.0,     0.0,     0.0},
            {       -12.0,       0.0,      0.0,        6.0,     0.0,     0.0},
            {        -8.0,       0.0,      0.0,        4.0,     0.0,     0.0},
            {        30.0,       0.0,      0.0,       -2.0,     0.0,     0.0},
            {        24.0,       0.0,      0.0,      -10.0,     0.0,     0.0},
            {        10.0,       0.0,      0.0,       -4.0,     0.0,     0.0},
            {       -16.0,       0.0,      0.0,        7.0,     0.0,     0.0},
            {       -16.0,       0.0,      0.0,        7.0,     0.0,     0.0},
            {        17.0,       0.0,      0.0,       -7.0,     0.0,     0.0},
            {       -24.0,       0.0,      0.0,       10.0,     0.0,     0.0},
            {       -12.0,       0.0,      0.0,        5.0,     0.0,     0.0},
            {       -24.0,       0.0,      0.0,       11.0,     0.0,     0.0},
            {       -23.0,       0.0,      0.0,        9.0,     0.0,     0.0},
            {       -13.0,       0.0,      0.0,        5.0,     0.0,     0.0},
            {       -15.0,       0.0,      0.0,        7.0,     0.0,     0.0},
            {         0.0,       0.0,  -1988.0,        0.0,     0.0, -1679.0},
            {         0.0,       0.0,    -63.0,        0.0,     0.0,   -27.0},
            {        -4.0,       0.0,      0.0,        0.0,     0.0,     0.0},
            {         0.0,       0.0,      5.0,        0.0,     0.0,     4.0},
            {         5.0,       0.0,      0.0,       -3.0,     0.0,     0.0},
            {         0.0,       0.0,    364.0,        0.0,     0.0,   176.0},
            {         0.0,       0.0,  -1044.0,        0.0,     0.0,  -891.0},
            {        -3.0,       0.0,      0.0,        1.0,     0.0,     0.0},
            {         4.0,       0.0,      0.0,       -2.0,     0.0,     0.0},
            {         0.0,       0.0,    330.0,        0.0,     0.0,     0.0},
            {         5.0,       0.0,      0.0,       -2.0,     0.0,     0.0},
            {         3.0,       0.0,      0.0,       -2.0,     0.0,     0.0},
            {        -3.0,       0.0,      0.0,        1.0,     0.0,     0.0},
            {        -5.0,       0.0,      0.0,        2.0,     0.0,     0.0},
            {         3.0,       0.0,      0.0,       -1.0,     0.0,     0.0},
            {         3.0,       0.0,      0.0,        0.0,     0.0,     0.0},
            {         3.0,       0.0,      0.0,        0.0,     0.0,     0.0},
            {         0.0,       0.0,      5.0,        0.0,     0.0,     0.0},
            {         0.0,       0.0,      0.0,        1.0,     0.0,     0.0},
            {         4.0,       0.0,      0.0,       -2.0,     0.0,     0.0},
            {         6.0,       0.0,      0.0,        0.0,     0.0,     0.0},
            {         5.0,       0.0,      0.0,       -2.0,     0.0,     0.0},
            {        -7.0,       0.0,      0.0,        0.0,     0.0,     0.0},
            {       -12.0,       0.0,      0.0,        0.0,     0.0,     0.0},
            {         5.0,       0.0,      0.0,       -3.0,     0.0,     0.0},
            {         3.0,       0.0,      0.0,       -1.0,     0.0,     0.0},
            {        -5.0,       0.0,      0.0,        0.0,     0.0,     0.0},
            {         3.0,       0.0,      0.0,        0.0,     0.0,     0.0},
            {        -7.0,       0.0,      0.0,        3.0,     0.0,     0.0},
            {         7.0,       0.0,      0.0,       -4.0,     0.0,     0.0},
            {         0.0,       0.0,    -12.0,        0.0,     0.0,   -10.0},
            {         4.0,       0.0,      0.0,       -2.0,     0.0,     0.0},
            {         3.0,       0.0,      0.0,       -2.0,     0.0,     0.0},
            {        -3.0,       0.0,      0.0,        2.0,     0.0,     0.0},
            {        -7.0,       0.0,      0.0,        3.0,     0.0,     0.0},
            {        -4.0,       0.0,      0.0,        2.0,     0.0,     0.0},
            {        -3.0,       0.0,      0.0,        1.0,     0.0,     0.0},
            {         0.0,       0.0,      0.0,        0.0,     0.0,     0.0},
            {        -3.0,       0.0,      0.0,        1.0,     0.0,     0.0},
            {         7.0,       0.0,      0.0,       -3.0,     0.0,     0.0},
            {        -4.0,       0.0,      0.0,        2.0,     0.0,     0.0},
            {         4.0,       0.0,      0.0,       -2.0,     0.0,     0.0},
            {        -5.0,       0.0,      0.0,        3.0,     0.0,     0.0},
            {         5.0,       0.0,      0.0,        0.0,     0.0,     0.0},
            {        -5.0,       0.0,      0.0,        2.0,     0.0,     0.0},
            {         5.0,       0.0,      0.0,       -2.0,     0.0,     0.0},
            {        -8.0,       0.0,      0.0,        3.0,     0.0,     0.0},
            {         9.0,       0.0,      0.0,        0.0,     0.0,     0.0},
            {         6.0,       0.0,      0.0,       -3.0,     0.0,     0.0},
            {        -5.0,       0.0,      0.0,        2.0,     0.0,     0.0},
            {         3.0,       0.0,      0.0,        0.0,     0.0,     0.0},
            {        -7.0,       0.0,      0.0,        0.0,     0.0,     0.0},
            {        -3.0,       0.0,      0.0,        1.0,     0.0,     0.0},
            {         5.0,       0.0,      0.0,        0.0,     0.0,     0.0},
            {         3.0,       0.0,      0.0,        0.0,     0.0,     0.0},
            {        -3.0,       0.0,      0.0,        2.0,     0.0,     0.0},
            {         4.0,       0.0,      0.0,       -2.0,     0.0,     0.0},
            {         3.0,       0.0,      0.0,       -1.0,     0.0,     0.0},
            {        -5.0,       0.0,      0.0,        2.0,     0.0,     0.0},
            {         4.0,       0.0,      0.0,       -2.0,     0.0,     0.0},
            {         9.0,       0.0,      0.0,       -3.0,     0.0,     0.0},
            {         4.0,       0.0,      0.0,        0.0,     0.0,     0.0},
            {         4.0,       0.0,      0.0,       -2.0,     0.0,     0.0},
            {        -3.0,       0.0,      0.0,        2.0,     0.0,     0.0},
            {        -4.0,       0.0,      0.0,        2.0,     0.0,     0.0},
            {         9.0,       0.0,      0.0,       -3.0,     0.0,     0.0},
            {        -4.0,       0.0,      0.0,        0.0,     0.0,     0.0},
            {        -4.0,       0.0,      0.0,        0.0,     0.0,     0.0},
            {         3.0,       0.0,      0.0,       -2.0,     0.0,     0.0},
            {         8.0,       0.0,      0.0,        0.0,     0.0,     0.0},
            {         3.0,       0.0,      0.0,        0.0,     0.0,     0.0},
            {        -3.0,       0.0,      0.0,        2.0,     0.0,     0.0},
            {         3.0,       0.0,      0.0,       -1.0,     0.0,     0.0},
            {         3.0,       0.0,      0.0,       -1.0,     0.0,     0.0},
            {        -3.0,       0.0,      0.0,        1.0,     0.0,     0.0},
            {         6.0,       0.0,      0.0,       -3.0,     0.0,     0.0},
            {         3.0,       0.0,      0.0,        0.0,     0.0,     0.0},
            {        -3.0,       0.0,      0.0,        1.0,     0.0,     0.0},
            {        -7.0,       0.0,      0.0,        0.0,     0.0,     0.0},
            {         9.0,       0.0,      0.0,        0.0,     0.0,     0.0},
            {        -3.0,       0.0,      0.0,        2.0,     0.0,     0.0},
            {        -3.0,       0.0,      0.0,        0.0,     0.0,     0.0},
            {        -4.0,       0.0,      0.0,        0.0,     0.0,     0.0},
            {        -5.0,       0.0,      0.0,        3.0,     0.0,     0.0},
            {       -13.0,       0.0,      0.0,        0.0,     0.0,     0.0},
            {        -7.0,       0.0,      0.0,        0.0,     0.0,     0.0},
            {        10.0,       0.0,      0.0,        0.0,     0.0,     0.0},
            {         3.0,       0.0,      0.0,       -1.0,     0.0,     0.0},
            {        10.0,       0.0,     13.0,        6.0,     0.0,    -5.0},
            {         0.0,       0.0,     30.0,        0.0,     0.0,    14.0},
            {         0.0,       0.0,   -162.0,        0.0,     0.0,  -138.0},
            {         0.0,       0.0,     75.0,        0.0,     0.0,     0.0},
            {        -7.0,       0.0,      0.0,        4.0,     0.0,     0.0},
            {        -4.0,       0.0,      0.0,        2.0,     0.0,     0.0},
            {         4.0,       0.0,      0.0,       -2.0,     0.0,     0.0},
            {         5.0,       0.0,      0.0,       -2.0,     0.0,     0.0},
            {         5.0,       0.0,      0.0,       -3.0,     0.0,     0.0},
            {        -3.0,       0.0,      0.0,        0.0,     0.0,     0.0},
            {        -3.0,       0.0,      0.0,        2.0,     0.0,     0.0},
            {        -4.0,       0.0,      0.0,        2.0,     0.0,     0.0},
            {        -5.0,       0.0,      0.0,        2.0,     0.0,     0.0},
            {         6.0,       0.0,      0.0,        0.0,     0.0,     0.0},
            {         9.0,       0.0,      0.0,        0.0,     0.0,     0.0},
            {         5.0,       0.0,      0.0,        0.0,     0.0,     0.0},
            {        -7.0,       0.0,      0.0,        0.0,     0.0,     0.0},
            {        -3.0,       0.0,      0.0,        1.0,     0.0,     0.0},
            {        -4.0,       0.0,      0.0,        2.0,     0.0,     0.0},
            {         7.0,       0.0,      0.0,        0.0,     0.0,     0.0},
            {        -4.0,       0.0,      0.0,        0.0,     0.0,     0.0},
            {         4.0,       0.0,      0.0,        0.0,     0.0,     0.0},
            {        -6.0,       0.0,     -3.0,        3.0,     0.0,     1.0},
            {         0.0,       0.0,     -3.0,        0.0,     0.0,    -2.0},
            {        11.0,       0.0,      0.0,        0.0,     0.0,     0.0},
            {         3.0,       0.0,      0.0,       -1.0,     0.0,     0.0},
            {        11.0,       0.0,      0.0,        0.0,     0.0,     0.0},
            {        -3.0,       0.0,      0.0,        2.0,     0.0,     0.0},
            {        -1.0,       0.0,      3.0,        3.0,     0.0,    -1.0},
            {         4.0,       0.0,      0.0,       -2.0,     0.0,     0.0},
            {         0.0,       0.0,    -13.0,        0.0,     0.0,   -11.0},
            {         3.0,       0.0,      6.0,        0.0,     0.0,     0.0},
            {        -7.0,       0.0,      0.0,        0.0,     0.0,     0.0},
            {         5.0,       0.0,      0.0,       -3.0,     0.0,     0.0},
            {        -3.0,       0.0,      0.0,        1.0,     0.0,     0.0},
            {         3.0,       0.0,      0.0,        0.0,     0.0,     0.0},
            {         5.0,       0.0,      0.0,       -3.0,     0.0,     0.0},
            {        -7.0,       0.0,      0.0,        3.0,     0.0,     0.0},
            {         8.0,       0.0,      0.0,       -3.0,     0.0,     0.0},
            {        -4.0,       0.0,      0.0,        2.0,     0.0,     0.0},
            {        11.0,       0.0,      0.0,        0.0,     0.0,     0.0},
            {        -3.0,       0.0,      0.0,        1.0,     0.0,     0.0},
            {         3.0,       0.0,      0.0,       -1.0,     0.0,     0.0},
            {        -4.0,       0.0,      0.0,        2.0,     0.0,     0.0},
            {         8.0,       0.0,      0.0,       -4.0,     0.0,     0.0},
            {         3.0,       0.0,      0.0,       -1.0,     0.0,     0.0},
            {        11.0,       0.0,      0.0,        0.0,     0.0,     0.0},
            {        -6.0,       0.0,      0.0,        3.0,     0.0,     0.0},
            {        -4.0,       0.0,      0.0,        2.0,     0.0,     0.0},
            {        -8.0,       0.0,      0.0,        4.0,     0.0,     0.0},
            {        -7.0,       0.0,      0.0,        3.0,     0.0,     0.0},
            {        -4.0,       0.0,      0.0,        2.0,     0.0,     0.0},
            {         3.0,       0.0,      0.0,       -1.0,     0.0,     0.0},
            {         6.0,       0.0,      0.0,       -3.0,     0.0,     0.0},
            {        -6.0,       0.0,      0.0,        3.0,     0.0,     0.0},
            {         6.0,       0.0,      0.0,        0.0,     0.0,     0.0},
            {         6.0,       0.0,      0.0,       -1.0,     0.0,     0.0},
            {         5.0,       0.0,      0.0,       -2.0,     0.0,     0.0},
            {        -5.0,       0.0,      0.0,        2.0,     0.0,     0.0},
            {        -4.0,       0.0,      0.0,        0.0,     0.0,     0.0},
            {        -4.0,       0.0,      0.0,        2.0,     0.0,     0.0},
            {         4.0,       0.0,      0.0,        0.0,     0.0,     0.0},
            {         6.0,       0.0,      0.0,       -3.0,     0.0,     0.0},
            {        -4.0,       0.0,      0.0,        2.0,     0.0,     0.0},
            {         0.0,       0.0,    -26.0,        0.0,     0.0,   -11.0},
            {         0.0,       0.0,    -10.0,        0.0,     0.0,    -5.0},
            {         5.0,       0.0,      0.0,       -3.0,     0.0,     0.0},
            {       -13.0,       0.0,      0.0,        0.0,     0.0,     0.0},
            {         3.0,       0.0,      0.0,       -2.0,     0.0,     0.0},
            {         4.0,       0.0,      0.0,       -2.0,     0.0,     0.0},
            {         7.0,       0.0,      0.0,       -3.0,     0.0,     0.0},
            {         4.0,       0.0,      0.0,        0.0,     0.0,     0.0},
            {         5.0,       0.0,      0.0,        0.0,     0.0,     0.0},
            {        -3.0,       0.0,      0.0,        2.0,     0.0,     0.0},
            {        -6.0,       0.0,      0.0,        2.0,     0.0,     0.0},
            {        -5.0,       0.0,      0.0,        2.0,     0.0,     0.0},
            {        -7.0,       0.0,      0.0,        3.0,     0.0,     0.0},
            {         5.0,       0.0,      0.0,       -2.0,     0.0,     0.0},
            {        13.0,       0.0,      0.0,        0.0,     0.0,     0.0},
            {        -4.0,       0.0,      0.0,        2.0,     0.0,     0.0},
            {        -3.0,       0.0,      0.0,        0.0,     0.0,     0.0},
            {         5.0,       0.0,      0.0,       -2.0,     0.0,     0.0},
            {       -11.0,       0.0,      0.0,        0.0,     0.0,     0.0},
            {         5.0,       0.0,      0.0,       -2.0,     0.0,     0.0},
            {         4.0,       0.0,      0.0,        0.0,     0.0,     0.0},
            {         4.0,       0.0,      0.0,       -2.0,     0.0,     0.0},
            {        -4.0,       0.0,      0.0,        2.0,     0.0,     0.0},
            {         6.0,       0.0,      0.0,       -3.0,     0.0,     0.0},
            {         3.0,       0.0,      0.0,       -2.0,     0.0,     0.0},
            {       -12.0,       0.0,      0.0,        0.0,     0.0,     0.0},
            {         4.0,       0.0,      0.0,        0.0,     0.0,     0.0},
            {        -3.0,       0.0,      0.0,        0.0,     0.0,     0.0},
            {        -4.0,       0.0,      0.0,        0.0,     0.0,     0.0},
            {         3.0,       0.0,      0.0,        0.0,     0.0,     0.0},
            {         3.0,       0.0,      0.0,       -1.0,     0.0,     0.0},
            {        -3.0,       0.0,      0.0,        1.0,     0.0,     0.0},
            {         0.0,       0.0,     -5.0,        0.0,     0.0,    -2.0},
            {        -7.0,       0.0,      0.0,        4.0,     0.0,     0.0},
            {         6.0,       0.0,      0.0,       -3.0,     0.0,     0.0},
            {        -3.0,       0.0,      0.0,        0.0,     0.0,     0.0},
            {         5.0,       0.0,      0.0,       -3.0,     0.0,     0.0},
            {         3.0,       0.0,      0.0,       -1.0,     0.0,     0.0},
            {         3.0,       0.0,      0.0,        0.0,     0.0,     0.0},
            {        -3.0,       0.0,      0.0,        1.0,     0.0,     0.0},
            {        -5.0,       0.0,      0.0,        3.0,     0.0,     0.0},
            {        -3.0,       0.0,      0.0,        2.0,     0.0,     0.0},
            {        -3.0,       0.0,      0.0,        2.0,     0.0,     0.0},
            {        12.0,       0.0,      0.0,        0.0,     0.0,     0.0},
            {         3.0,       0.0,      0.0,       -1.0,     0.0,     0.0},
            {        -4.0,       0.0,      0.0,        2.0,     0.0,     0.0},
            {         4.0,       0.0,      0.0,        0.0,     0.0,     0.0},
            {         6.0,       0.0,      0.0,        0.0,     0.0,     0.0},
            {         5.0,       0.0,      0.0,       -3.0,     0.0,     0.0},
            {         4.0,       0.0,      0.0,       -2.0,     0.0,     0.0},
            {        -6.0,       0.0,      0.0,        3.0,     0.0,     0.0},
            {         4.0,       0.0,      0.0,       -2.0,     0.0,     0.0},
            {         6.0,       0.0,      0.0,       -3.0,     0.0,     0.0},
            {         6.0,       0.0,      0.0,        0.0,     0.0,     0.0},
            {        -6.0,       0.0,      0.0,        3.0,     0.0,     0.0},
            {         3.0,       0.0,      0.0,       -2.0,     0.0,     0.0},
            {         7.0,       0.0,      0.0,       -4.0,     0.0,     0.0},
            {         4.0,       0.0,      0.0,       -2.0,     0.0,     0.0},
            {        -5.0,       0.0,      0.0,        2.0,     0.0,     0.0},
            {         5.0,       0.0,      0.0,        0.0,     0.0,     0.0},
            {        -6.0,       0.0,      0.0,        3.0,     0.0,     0.0},
            {        -6.0,       0.0,      0.0,        3.0,     0.0,     0.0},
            {        -4.0,       0.0,      0.0,        2.0,     0.0,     0.0},
            {        10.0,       0.0,      0.0,        0.0,     0.0,     0.0},
            {        -4.0,       0.0,      0.0,        2.0,     0.0,     0.0},
            {         7.0,       0.0,      0.0,        0.0,     0.0,     0.0},
            {         7.0,       0.0,      0.0,       -3.0,     0.0,     0.0},
            {         4.0,       0.0,      0.0,        0.0,     0.0,     0.0},
            {        11.0,       0.0,      0.0,        0.0,     0.0,     0.0},
            {         5.0,       0.0,      0.0,       -2.0,     0.0,     0.0},
            {        -6.0,       0.0,      0.0,        2.0,     0.0,     0.0},
            {         4.0,       0.0,      0.0,       -2.0,     0.0,     0.0},
            {         3.0,       0.0,      0.0,       -2.0,     0.0,     0.0},
            {         5.0,       0.0,      0.0,       -2.0,     0.0,     0.0},
            {        -4.0,       0.0,      0.0,        2.0,     0.0,     0.0},
            {        -4.0,       0.0,      0.0,        2.0,     0.0,     0.0},
            {        -3.0,       0.0,      0.0,        2.0,     0.0,     0.0},
            {         4.0,       0.0,      0.0,       -2.0,     0.0,     0.0},
            {         3.0,       0.0,      0.0,       -1.0,     0.0,     0.0},
            {        -3.0,       0.0,      0.0,        1.0,     0.0,     0.0},
            {        -3.0,       0.0,      0.0,        1.0,     0.0,     0.0},
            {        -3.0,       0.0,      0.0,        2.0,     0.0,     0.0}};
        
        /*
         Planetary argument multipliers:
         L   L'  F   D   Om  Me  Ve  E  Ma  Ju  Sa  Ur  Ne  pre
         */
        static const short int napl_t[687][14] = {
            { 0,  0,  0,  0,  0,  0,  0,  8,-16,  4,  5,  0,  0,  0},
            { 0,  0,  0,  0,  0,  0,  0, -8, 16, -4, -5,  0,  0,  2},
            { 0,  0,  0,  0,  0,  0,  0,  8,-16,  4,  5,  0,  0,  2},
            { 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, -1,  2,  2},
            { 0,  0,  0,  0,  0,  0,  0, -4,  8, -1, -5,  0,  0,  2},
            { 0,  0,  0,  0,  0,  0,  0,  4, -8,  3,  0,  0,  0,  1},
            { 0,  0,  1, -1,  1,  0,  0,  3, -8,  3,  0,  0,  0,  0},
            {-1,  0,  0,  0,  0,  0, 10, -3,  0,  0,  0,  0,  0,  0},
            { 0,  0,  0,  0,  0,  0,  0,  0,  0, -2,  6, -3,  0,  2},
            { 0,  0,  0,  0,  0,  0,  0,  4, -8,  3,  0,  0,  0,  0},
            { 0,  0,  1, -1,  1,  0,  0, -5,  8, -3,  0,  0,  0,  0},
            { 0,  0,  0,  0,  0,  0,  0, -4,  8, -3,  0,  0,  0,  1},
            { 0,  0,  0,  0,  0,  0,  0,  4, -8,  1,  5,  0,  0,  2},
            { 0,  0,  0,  0,  0,  0, -5,  6,  4,  0,  0,  0,  0,  2},
            { 0,  0,  0,  0,  0,  0,  0,  0,  0,  2, -5,  0,  0,  2},
            { 0,  0,  0,  0,  0,  0,  0,  0,  0,  2, -5,  0,  0,  1},
            { 0,  0,  1, -1,  1,  0,  0, -1,  0,  2, -5,  0,  0,  0},
            { 0,  0,  0,  0,  0,  0,  0,  0,  0,  2, -5,  0,  0,  0},
            { 0,  0,  1, -1,  1,  0,  0, -1,  0, -2,  5,  0,  0,  0},
            { 0,  0,  0,  0,  0,  0,  0,  0,  0, -2,  5,  0,  0,  1},
            { 0,  0,  0,  0,  0,  0,  0,  0,  0, -2,  5,  0,  0,  2},
            { 2,  0, -1, -1,  0,  0,  0,  3, -7,  0,  0,  0,  0,  0},
            { 1,  0,  0, -2,  0,  0, 19,-21,  3,  0,  0,  0,  0,  0},
            { 0,  0,  1, -1,  1,  0,  2, -4,  0, -3,  0,  0,  0,  0},
            { 1,  0,  0, -1,  1,  0,  0, -1,  0,  2,  0,  0,  0,  0},
            { 0,  0,  1, -1,  1,  0,  0, -1,  0, -4, 10,  0,  0,  0},
            {-2,  0,  0,  2,  1,  0,  0,  2,  0,  0, -5,  0,  0,  0},
            { 0,  0,  0,  0,  0,  0,  3, -7,  4,  0,  0,  0,  0,  0},
            { 0,  0, -1,  1,  0,  0,  0,  1,  0,  1, -1,  0,  0,  0},
            {-2,  0,  0,  2,  1,  0,  0,  2,  0, -2,  0,  0,  0,  0},
            {-1,  0,  0,  0,  0,  0, 18,-16,  0,  0,  0,  0,  0,  0},
            {-2,  0,  1,  1,  2,  0,  0,  1,  0, -2,  0,  0,  0,  0},
            {-1,  0,  1, -1,  1,  0, 18,-17,  0,  0,  0,  0,  0,  0},
            {-1,  0,  0,  1,  1,  0,  0,  2, -2,  0,  0,  0,  0,  0},
            { 0,  0,  0,  0,  0,  0, -8, 13,  0,  0,  0,  0,  0,  2},
            { 0,  0,  2, -2,  2,  0, -8, 11,  0,  0,  0,  0,  0,  0},
            { 0,  0,  0,  0,  0,  0, -8, 13,  0,  0,  0,  0,  0,  1},
            { 0,  0,  1, -1,  1,  0, -8, 12,  0,  0,  0,  0,  0,  0},
            { 0,  0,  0,  0,  0,  0,  8,-13,  0,  0,  0,  0,  0,  0},
            { 0,  0,  1, -1,  1,  0,  8,-14,  0,  0,  0,  0,  0,  0},
            { 0,  0,  0,  0,  0,  0,  8,-13,  0,  0,  0,  0,  0,  1},
            {-2,  0,  0,  2,  1,  0,  0,  2,  0, -4,  5,  0,  0,  0},
            {-2,  0,  0,  2,  2,  0,  3, -3,  0,  0,  0,  0,  0,  0},
            {-2,  0,  0,  2,  0,  0,  0,  2,  0, -3,  1,  0,  0,  0},
            { 0,  0,  0,  0,  1,  0,  3, -5,  0,  2,  0,  0,  0,  0},
            {-2,  0,  0,  2,  0,  0,  0,  2,  0, -4,  3,  0,  0,  0},
            { 0,  0, -1,  1,  0,  0,  0,  0,  2,  0,  0,  0,  0,  0},
            { 0,  0,  0,  0,  1,  0,  0, -1,  2,  0,  0,  0,  0,  0},
            { 0,  0,  1, -1,  2,  0,  0, -2,  2,  0,  0,  0,  0,  0},
            {-1,  0,  1,  0,  1,  0,  3, -5,  0,  0,  0,  0,  0,  0},
            {-1,  0,  0,  1,  0,  0,  3, -4,  0,  0,  0,  0,  0,  0},
            {-2,  0,  0,  2,  0,  0,  0,  2,  0, -2, -2,  0,  0,  0},
            {-2,  0,  2,  0,  2,  0,  0, -5,  9,  0,  0,  0,  0,  0},
            { 0,  0,  1, -1,  1,  0,  0, -1,  0,  0,  0, -1,  0,  0},
            { 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0},
            { 0,  0,  1, -1,  1,  0,  0, -1,  0,  0,  0,  0,  2,  0},
            { 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  2,  1},
            { 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  2,  2},
            {-1,  0,  0,  1,  0,  0,  0,  3, -4,  0,  0,  0,  0,  0},
            { 0,  0, -1,  1,  0,  0,  0,  1,  0,  0,  2,  0,  0,  0},
            { 0,  0,  1, -1,  2,  0,  0, -1,  0,  0,  2,  0,  0,  0},
            { 0,  0,  0,  0,  1,  0,  0, -9, 17,  0,  0,  0,  0,  0},
            { 0,  0,  0,  0,  2,  0, -3,  5,  0,  0,  0,  0,  0,  0},
            { 0,  0,  1, -1,  1,  0,  0, -1,  0, -1,  2,  0,  0,  0},
            { 0,  0,  0,  0,  0,  0,  0,  0,  0,  1, -2,  0,  0,  0},
            { 1,  0,  0, -2,  0,  0, 17,-16,  0, -2,  0,  0,  0,  0},
            { 0,  0,  1, -1,  1,  0,  0, -1,  0,  1, -3,  0,  0,  0},
            {-2,  0,  0,  2,  1,  0,  0,  5, -6,  0,  0,  0,  0,  0},
            { 0,  0, -2,  2,  0,  0,  0,  9,-13,  0,  0,  0,  0,  0},
            { 0,  0,  1, -1,  2,  0,  0, -1,  0,  0,  1,  0,  0,  0},
            { 0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  1,  0,  0,  0},
            { 0,  0, -1,  1,  0,  0,  0,  1,  0,  0,  1,  0,  0,  0},
            { 0,  0, -2,  2,  0,  0,  5, -6,  0,  0,  0,  0,  0,  0},
            { 0,  0, -1,  1,  1,  0,  5, -7,  0,  0,  0,  0,  0,  0},
            {-2,  0,  0,  2,  0,  0,  6, -8,  0,  0,  0,  0,  0,  0},
            { 2,  0,  1, -3,  1,  0, -6,  7,  0,  0,  0,  0,  0,  0},
            { 0,  0,  0,  0,  2,  0,  0,  0,  0,  1,  0,  0,  0,  0},
            { 0,  0, -1,  1,  1,  0,  0,  1,  0,  1,  0,  0,  0,  0},
            { 0,  0,  1, -1,  1,  0,  0, -1,  0,  0,  0,  2,  0,  0},
            { 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  2,  0,  1},
            { 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  2,  0,  2},
            { 0,  0,  0,  0,  0,  0,  0, -8, 15,  0,  0,  0,  0,  2},
            { 0,  0,  0,  0,  0,  0,  0, -8, 15,  0,  0,  0,  0,  1},
            { 0,  0,  1, -1,  1,  0,  0, -9, 15,  0,  0,  0,  0,  0},
            { 0,  0,  0,  0,  0,  0,  0,  8,-15,  0,  0,  0,  0,  0},
            { 1,  0, -1, -1,  0,  0,  0,  8,-15,  0,  0,  0,  0,  0},
            { 2,  0,  0, -2,  0,  0,  2, -5,  0,  0,  0,  0,  0,  0},
            {-2,  0,  0,  2,  0,  0,  0,  2,  0, -5,  5,  0,  0,  0},
            { 2,  0,  0, -2,  1,  0,  0, -6,  8,  0,  0,  0,  0,  0},
            { 2,  0,  0, -2,  1,  0,  0, -2,  0,  3,  0,  0,  0,  0},
            {-2,  0,  1,  1,  0,  0,  0,  1,  0, -3,  0,  0,  0,  0},
            {-2,  0,  1,  1,  1,  0,  0,  1,  0, -3,  0,  0,  0,  0},
            {-2,  0,  0,  2,  0,  0,  0,  2,  0, -3,  0,  0,  0,  0},
            {-2,  0,  0,  2,  0,  0,  0,  6, -8,  0,  0,  0,  0,  0},
            {-2,  0,  0,  2,  0,  0,  0,  2,  0, -1, -5,  0,  0,  0},
            {-1,  0,  0,  1,  0,  0,  0,  1,  0, -1,  0,  0,  0,  0},
            {-1,  0,  1,  1,  1,  0,-20, 20,  0,  0,  0,  0,  0,  0},
            { 1,  0,  0, -2,  0,  0, 20,-21,  0,  0,  0,  0,  0,  0},
            { 0,  0,  0,  0,  1,  0,  0,  8,-15,  0,  0,  0,  0,  0},
            { 0,  0,  2, -2,  1,  0,  0,-10, 15,  0,  0,  0,  0,  0},
            { 0,  0, -1,  1,  0,  0,  0,  1,  0,  1,  0,  0,  0,  0},
            { 0,  0,  0,  0,  1,  0,  0,  0,  0,  1,  0,  0,  0,  0},
            { 0,  0,  1, -1,  2,  0,  0, -1,  0,  1,  0,  0,  0,  0},
            { 0,  0,  1, -1,  1,  0,  0, -1,  0, -2,  4,  0,  0,  0},
            { 2,  0,  0, -2,  1,  0, -6,  8,  0,  0,  0,  0,  0,  0},
            { 0,  0, -2,  2,  1,  0,  5, -6,  0,  0,  0,  0,  0,  0},
            { 0,  0,  0,  0,  0,  0,  0,  0,  0,  0, -1,  0,  0,  1},
            { 0,  0,  1, -1,  1,  0,  0, -1,  0,  0, -1,  0,  0,  0},
            { 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0},
            { 0,  0,  1, -1,  1,  0,  0, -1,  0,  0,  1,  0,  0,  0},
            { 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  1},
            { 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  2},
            { 0,  0,  2, -2,  1,  0,  0, -9, 13,  0,  0,  0,  0,  0},
            { 0,  0,  0,  0,  1,  0,  0,  7,-13,  0,  0,  0,  0,  0},
            {-2,  0,  0,  2,  0,  0,  0,  5, -6,  0,  0,  0,  0,  0},
            { 0,  0,  0,  0,  0,  0,  0,  9,-17,  0,  0,  0,  0,  0},
            { 0,  0,  0,  0,  0,  0,  0, -9, 17,  0,  0,  0,  0,  2},
            { 1,  0,  0, -1,  1,  0,  0, -3,  4,  0,  0,  0,  0,  0},
            { 1,  0,  0, -1,  1,  0, -3,  4,  0,  0,  0,  0,  0,  0},
            { 0,  0,  0,  0,  2,  0,  0, -1,  2,  0,  0,  0,  0,  0},
            { 0,  0, -1,  1,  1,  0,  0,  0,  2,  0,  0,  0,  0,  0},
            { 0,  0, -2,  2,  0,  1,  0, -2,  0,  0,  0,  0,  0,  0},
            { 0,  0,  0,  0,  0,  0,  3, -5,  0,  2,  0,  0,  0,  0},
            {-2,  0,  0,  2,  1,  0,  0,  2,  0, -3,  1,  0,  0,  0},
            {-2,  0,  0,  2,  1,  0,  3, -3,  0,  0,  0,  0,  0,  0},
            { 0,  0,  0,  0,  1,  0,  8,-13,  0,  0,  0,  0,  0,  0},
            { 0,  0, -1,  1,  0,  0,  8,-12,  0,  0,  0,  0,  0,  0},
            { 0,  0,  2, -2,  1,  0, -8, 11,  0,  0,  0,  0,  0,  0},
            {-1,  0,  0,  1,  0,  0,  0,  2, -2,  0,  0,  0,  0,  0},
            {-1,  0,  0,  0,  1,  0, 18,-16,  0,  0,  0,  0,  0,  0},
            { 0,  0,  1, -1,  1,  0,  0, -1,  0, -1,  1,  0,  0,  0},
            { 0,  0,  0,  0,  1,  0,  3, -7,  4,  0,  0,  0,  0,  0},
            {-2,  0,  1,  1,  1,  0,  0, -3,  7,  0,  0,  0,  0,  0},
            { 0,  0,  1, -1,  2,  0,  0, -1,  0, -2,  5,  0,  0,  0},
            { 0,  0,  0,  0,  1,  0,  0,  0,  0, -2,  5,  0,  0,  0},
            { 0,  0,  0,  0,  1,  0,  0, -4,  8, -3,  0,  0,  0,  0},
            { 1,  0,  0,  0,  1,  0,-10,  3,  0,  0,  0,  0,  0,  0},
            { 0,  0,  2, -2,  1,  0,  0, -2,  0,  0,  0,  0,  0,  0},
            {-1,  0,  0,  0,  1,  0, 10, -3,  0,  0,  0,  0,  0,  0},
            { 0,  0,  0,  0,  1,  0,  0,  4, -8,  3,  0,  0,  0,  0},
            { 0,  0,  0,  0,  1,  0,  0,  0,  0,  2, -5,  0,  0,  0},
            { 0,  0, -1,  1,  0,  0,  0,  1,  0,  2, -5,  0,  0,  0},
            { 2,  0, -1, -1,  1,  0,  0,  3, -7,  0,  0,  0,  0,  0},
            {-2,  0,  0,  2,  0,  0,  0,  2,  0,  0, -5,  0,  0,  0},
            { 0,  0,  0,  0,  1,  0, -3,  7, -4,  0,  0,  0,  0,  0},
            {-2,  0,  0,  2,  0,  0,  0,  2,  0, -2,  0,  0,  0,  0},
            { 1,  0,  0,  0,  1,  0,-18, 16,  0,  0,  0,  0,  0,  0},
            {-2,  0,  1,  1,  1,  0,  0,  1,  0, -2,  0,  0,  0,  0},
            { 0,  0,  1, -1,  2,  0, -8, 12,  0,  0,  0,  0,  0,  0},
            { 0,  0,  0,  0,  1,  0, -8, 13,  0,  0,  0,  0,  0,  0},
            { 0,  0,  0,  0,  0,  0,  0,  1, -2,  0,  0,  0,  0,  1},
            { 0,  0,  1, -1,  1,  0,  0,  0, -2,  0,  0,  0,  0,  0},
            { 0,  0,  0,  0,  0,  0,  0,  1, -2,  0,  0,  0,  0,  0},
            { 0,  0,  1, -1,  1,  0,  0, -2,  2,  0,  0,  0,  0,  0},
            { 0,  0,  0,  0,  0,  0,  0, -1,  2,  0,  0,  0,  0,  1},
            {-1,  0,  0,  1,  1,  0,  3, -4,  0,  0,  0,  0,  0,  0},
            {-1,  0,  0,  1,  1,  0,  0,  3, -4,  0,  0,  0,  0,  0},
            { 0,  0,  1, -1,  1,  0,  0, -1,  0,  0, -2,  0,  0,  0},
            { 0,  0,  1, -1,  1,  0,  0, -1,  0,  0,  2,  0,  0,  0},
            { 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  2,  0,  0,  1},
            { 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  2,  0,  0,  2},
            { 0,  0,  1, -1,  0,  0,  3, -6,  0,  0,  0,  0,  0,  0},
            { 0,  0,  0,  0,  1,  0, -3,  5,  0,  0,  0,  0,  0,  0},
            { 0,  0,  1, -1,  2,  0, -3,  4,  0,  0,  0,  0,  0,  0},
            { 0,  0,  0,  0,  1,  0,  0, -2,  4,  0,  0,  0,  0,  0},
            { 0,  0,  2, -2,  1,  0, -5,  6,  0,  0,  0,  0,  0,  0},
            { 0,  0, -1,  1,  0,  0,  5, -7,  0,  0,  0,  0,  0,  0},
            { 0,  0,  0,  0,  1,  0,  5, -8,  0,  0,  0,  0,  0,  0},
            {-2,  0,  0,  2,  1,  0,  6, -8,  0,  0,  0,  0,  0,  0},
            { 0,  0,  0,  0,  1,  0,  0, -8, 15,  0,  0,  0,  0,  0},
            {-2,  0,  0,  2,  1,  0,  0,  2,  0, -3,  0,  0,  0,  0},
            {-2,  0,  0,  2,  1,  0,  0,  6, -8,  0,  0,  0,  0,  0},
            { 1,  0,  0, -1,  1,  0,  0, -1,  0,  1,  0,  0,  0,  0},
            { 0,  0,  0,  0,  0,  0,  0,  0,  0,  3, -5,  0,  0,  0},
            { 0,  0,  1, -1,  1,  0,  0, -1,  0, -1,  0,  0,  0,  0},
            { 0,  0,  0,  0,  0,  0,  0,  0,  0, -1,  0,  0,  0,  1},
            { 0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0},
            { 0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  1},
            { 0,  0,  1, -1,  1,  0,  0, -1,  0,  1,  0,  0,  0,  0},
            { 0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  1},
            { 0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  2},
            { 0,  0,  1, -1,  2,  0,  0, -1,  0,  0, -1,  0,  0,  0},
            { 0,  0,  0,  0,  1,  0,  0,  0,  0,  0, -1,  0,  0,  0},
            { 0,  0, -1,  1,  0,  0,  0,  1,  0,  0, -1,  0,  0,  0},
            { 0,  0,  0,  0,  0,  0,  0, -7, 13,  0,  0,  0,  0,  2},
            { 0,  0,  0,  0,  0,  0,  0,  7,-13,  0,  0,  0,  0,  0},
            { 2,  0,  0, -2,  1,  0,  0, -5,  6,  0,  0,  0,  0,  0},
            { 0,  0,  2, -2,  1,  0,  0, -8, 11,  0,  0,  0,  0,  0},
            { 0,  0,  2, -2,  1, -1,  0,  2,  0,  0,  0,  0,  0,  0},
            {-2,  0,  0,  2,  0,  0,  0,  4, -4,  0,  0,  0,  0,  0},
            { 0,  0,  0,  0,  0,  0,  0,  0,  0,  2, -2,  0,  0,  0},
            { 0,  0,  1, -1,  1,  0,  0, -1,  0,  0,  3,  0,  0,  0},
            { 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  3,  0,  0,  1},
            { 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  3,  0,  0,  2},
            {-2,  0,  0,  2,  0,  0,  3, -3,  0,  0,  0,  0,  0,  0},
            { 0,  0,  0,  0,  2,  0,  0, -4,  8, -3,  0,  0,  0,  0},
            { 0,  0,  0,  0,  2,  0,  0,  4, -8,  3,  0,  0,  0,  0},
            { 2,  0,  0, -2,  1,  0,  0, -2,  0,  2,  0,  0,  0,  0},
            { 0,  0,  1, -1,  2,  0,  0, -1,  0,  2,  0,  0,  0,  0},
            { 0,  0,  1, -1,  2,  0,  0,  0, -2,  0,  0,  0,  0,  0},
            { 0,  0,  0,  0,  1,  0,  0,  1, -2,  0,  0,  0,  0,  0},
            { 0,  0, -1,  1,  0,  0,  0,  2, -2,  0,  0,  0,  0,  0},
            { 0,  0, -1,  1,  0,  0,  0,  1,  0,  0, -2,  0,  0,  0},
            { 0,  0,  2, -2,  1,  0,  0, -2,  0,  0,  2,  0,  0,  0},
            { 0,  0,  1, -1,  1,  0,  3, -6,  0,  0,  0,  0,  0,  0},
            { 0,  0,  0,  0,  0,  0,  3, -5,  0,  0,  0,  0,  0,  1},
            { 0,  0,  0,  0,  0,  0,  3, -5,  0,  0,  0,  0,  0,  0},
            { 0,  0,  1, -1,  1,  0, -3,  4,  0,  0,  0,  0,  0,  0},
            { 0,  0,  0,  0,  0,  0, -3,  5,  0,  0,  0,  0,  0,  1},
            { 0,  0,  0,  0,  0,  0, -3,  5,  0,  0,  0,  0,  0,  2},
            { 0,  0,  2, -2,  2,  0, -3,  3,  0,  0,  0,  0,  0,  0},
            { 0,  0,  0,  0,  0,  0, -3,  5,  0,  0,  0,  0,  0,  2},
            { 0,  0,  0,  0,  0,  0,  0,  2, -4,  0,  0,  0,  0,  1},
            { 0,  0,  1, -1,  1,  0,  0,  1, -4,  0,  0,  0,  0,  0},
            { 0,  0,  0,  0,  0,  0,  0,  2, -4,  0,  0,  0,  0,  0},
            { 0,  0,  0,  0,  0,  0,  0, -2,  4,  0,  0,  0,  0,  1},
            { 0,  0,  1, -1,  1,  0,  0, -3,  4,  0,  0,  0,  0,  0},
            { 0,  0,  0,  0,  0,  0,  0, -2,  4,  0,  0,  0,  0,  1},
            { 0,  0,  0,  0,  0,  0,  0, -2,  4,  0,  0,  0,  0,  2},
            { 0,  0,  0,  0,  0,  0, -5,  8,  0,  0,  0,  0,  0,  2},
            { 0,  0,  2, -2,  2,  0, -5,  6,  0,  0,  0,  0,  0,  0},
            { 0,  0,  0,  0,  0,  0, -5,  8,  0,  0,  0,  0,  0,  2},
            { 0,  0,  0,  0,  0,  0, -5,  8,  0,  0,  0,  0,  0,  1},
            { 0,  0,  1, -1,  1,  0, -5,  7,  0,  0,  0,  0,  0,  0},
            { 0,  0,  0,  0,  0,  0, -5,  8,  0,  0,  0,  0,  0,  1},
            { 0,  0,  0,  0,  0,  0,  5, -8,  0,  0,  0,  0,  0,  0},
            { 0,  0,  1, -1,  2,  0,  0, -1,  0, -1,  0,  0,  0,  0},
            { 0,  0,  0,  0,  1,  0,  0,  0,  0, -1,  0,  0,  0,  0},
            { 0,  0, -1,  1,  0,  0,  0,  1,  0, -1,  0,  0,  0,  0},
            { 0,  0,  2, -2,  1,  0,  0, -2,  0,  1,  0,  0,  0,  0},
            { 0,  0,  0,  0,  0,  0,  0, -6, 11,  0,  0,  0,  0,  2},
            { 0,  0,  0,  0,  0,  0,  0,  6,-11,  0,  0,  0,  0,  0},
            { 0,  0,  0,  0,  0, -1,  0,  4,  0,  0,  0,  0,  0,  2},
            { 0,  0,  0,  0,  0,  1,  0, -4,  0,  0,  0,  0,  0,  0},
            { 2,  0,  0, -2,  1,  0, -3,  3,  0,  0,  0,  0,  0,  0},
            {-2,  0,  0,  2,  0,  0,  0,  2,  0,  0, -2,  0,  0,  0},
            { 0,  0,  2, -2,  1,  0,  0, -7,  9,  0,  0,  0,  0,  0},
            { 0,  0,  0,  0,  0,  0,  0,  0,  0,  4, -5,  0,  0,  2},
            { 0,  0,  0,  0,  0,  0,  0,  0,  0,  2,  0,  0,  0,  0},
            { 0,  0,  0,  0,  0,  0,  0,  0,  0,  2,  0,  0,  0,  1},
            { 0,  0,  1, -1,  1,  0,  0, -1,  0,  2,  0,  0,  0,  0},
            { 0,  0,  0,  0,  0,  0,  0,  0,  0,  2,  0,  0,  0,  1},
            { 0,  0,  0,  0,  0,  0,  0,  0,  0,  2,  0,  0,  0,  2},
            { 0,  0,  2, -2,  2,  0,  0, -2,  0,  2,  0,  0,  0,  0},
            { 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  5,  0,  0,  2},
            { 0,  0,  0,  0,  1,  0,  3, -5,  0,  0,  0,  0,  0,  0},
            { 0,  0, -1,  1,  0,  0,  3, -4,  0,  0,  0,  0,  0,  0},
            { 0,  0,  2, -2,  1,  0, -3,  3,  0,  0,  0,  0,  0,  0},
            { 0,  0,  0,  0,  1,  0,  0,  2, -4,  0,  0,  0,  0,  0},
            { 0,  0,  2, -2,  1,  0,  0, -4,  4,  0,  0,  0,  0,  0},
            { 0,  0,  1, -1,  2,  0, -5,  7,  0,  0,  0,  0,  0,  0},
            { 0,  0,  0,  0,  0,  0,  0,  3, -6,  0,  0,  0,  0,  0},
            { 0,  0,  0,  0,  0,  0,  0, -3,  6,  0,  0,  0,  0,  1},
            { 0,  0,  1, -1,  1,  0,  0, -4,  6,  0,  0,  0,  0,  0},
            { 0,  0,  0,  0,  0,  0,  0, -3,  6,  0,  0,  0,  0,  1},
            { 0,  0,  0,  0,  0,  0,  0, -3,  6,  0,  0,  0,  0,  2},
            { 0,  0, -1,  1,  0,  0,  2, -2,  0,  0,  0,  0,  0,  0},
            { 0,  0,  0,  0,  1,  0,  2, -3,  0,  0,  0,  0,  0,  0},
            { 0,  0,  0,  0,  0,  0,  0, -5,  9,  0,  0,  0,  0,  2},
            { 0,  0,  0,  0,  0,  0,  0, -5,  9,  0,  0,  0,  0,  1},
            { 0,  0,  0,  0,  0,  0,  0,  5, -9,  0,  0,  0,  0,  0},
            { 0,  0, -1,  1,  0,  0,  0,  1,  0, -2,  0,  0,  0,  0},
            { 0,  0,  2, -2,  1,  0,  0, -2,  0,  2,  0,  0,  0,  0},
            {-2,  0,  1,  1,  1,  0,  0,  1,  0,  0,  0,  0,  0,  0},
            { 0,  0, -2,  2,  0,  0,  3, -3,  0,  0,  0,  0,  0,  0},
            { 0,  0,  0,  0,  0,  0, -6, 10,  0,  0,  0,  0,  0,  1},
            { 0,  0,  0,  0,  0,  0, -6, 10,  0,  0,  0,  0,  0,  2},
            { 0,  0,  0,  0,  0,  0, -2,  3,  0,  0,  0,  0,  0,  2},
            { 0,  0,  0,  0,  0,  0, -2,  3,  0,  0,  0,  0,  0,  1},
            { 0,  0,  1, -1,  1,  0, -2,  2,  0,  0,  0,  0,  0,  0},
            { 0,  0,  0,  0,  0,  0,  2, -3,  0,  0,  0,  0,  0,  0},
            { 0,  0,  0,  0,  0,  0,  2, -3,  0,  0,  0,  0,  0,  1},
            { 0,  0,  0,  0,  0,  0,  0,  0,  0,  3,  0,  0,  0,  1},
            { 0,  0,  1, -1,  1,  0,  0, -1,  0,  3,  0,  0,  0,  0},
            { 0,  0,  0,  0,  0,  0,  0,  0,  0,  3,  0,  0,  0,  1},
            { 0,  0,  0,  0,  0,  0,  0,  0,  0,  3,  0,  0,  0,  2},
            { 0,  0,  0,  0,  0,  0,  0,  4, -8,  0,  0,  0,  0,  0},
            { 0,  0,  0,  0,  0,  0,  0, -4,  8,  0,  0,  0,  0,  2},
            { 0,  0, -2,  2,  0,  0,  0,  2,  0, -2,  0,  0,  0,  0},
            { 0,  0,  0,  0,  0,  0,  0, -4,  7,  0,  0,  0,  0,  2},
            { 0,  0,  0,  0,  0,  0,  0, -4,  7,  0,  0,  0,  0,  1},
            { 0,  0,  0,  0,  0,  0,  0,  4, -7,  0,  0,  0,  0,  0},
            { 0,  0,  0,  0,  1,  0, -2,  3,  0,  0,  0,  0,  0,  0},
            { 0,  0,  2, -2,  1,  0,  0, -2,  0,  3,  0,  0,  0,  0},
            { 0,  0,  0,  0,  0,  0,  0, -5, 10,  0,  0,  0,  0,  2},
            { 0,  0,  0,  0,  1,  0, -1,  2,  0,  0,  0,  0,  0,  0},
            { 0,  0,  0,  0,  0,  0,  0,  0,  0,  4,  0,  0,  0,  2},
            { 0,  0,  0,  0,  0,  0,  0, -3,  5,  0,  0,  0,  0,  2},
            { 0,  0,  0,  0,  0,  0,  0, -3,  5,  0,  0,  0,  0,  1},
            { 0,  0,  0,  0,  0,  0,  0,  3, -5,  0,  0,  0,  0,  0},
            { 0,  0,  0,  0,  0,  0,  1, -2,  0,  0,  0,  0,  0,  1},
            { 0,  0,  1, -1,  1,  0,  1, -3,  0,  0,  0,  0,  0,  0},
            { 0,  0,  0,  0,  0,  0,  1, -2,  0,  0,  0,  0,  0,  0},
            { 0,  0,  0,  0,  0,  0, -1,  2,  0,  0,  0,  0,  0,  1},
            { 0,  0,  0,  0,  0,  0, -1,  2,  0,  0,  0,  0,  0,  2},
            { 0,  0,  0,  0,  0,  0, -7, 11,  0,  0,  0,  0,  0,  2},
            { 0,  0,  0,  0,  0,  0, -7, 11,  0,  0,  0,  0,  0,  1},
            { 0,  0, -2,  2,  0,  0,  4, -4,  0,  0,  0,  0,  0,  0},
            { 0,  0,  0,  0,  0,  0,  0,  2, -3,  0,  0,  0,  0,  0},
            { 0,  0,  2, -2,  1,  0, -4,  4,  0,  0,  0,  0,  0,  0},
            { 0,  0, -1,  1,  0,  0,  4, -5,  0,  0,  0,  0,  0,  0},
            { 0,  0,  0,  0,  0,  0,  0,  1, -1,  0,  0,  0,  0,  0},
            { 0,  0,  0,  0,  0,  0, -4,  7,  0,  0,  0,  0,  0,  1},
            { 0,  0,  1, -1,  1,  0, -4,  6,  0,  0,  0,  0,  0,  0},
            { 0,  0,  0,  0,  0,  0, -4,  7,  0,  0,  0,  0,  0,  2},
            { 0,  0,  0,  0,  0,  0, -4,  6,  0,  0,  0,  0,  0,  2},
            { 0,  0,  0,  0,  0,  0, -4,  6,  0,  0,  0,  0,  0,  1},
            { 0,  0,  1, -1,  1,  0, -4,  5,  0,  0,  0,  0,  0,  0},
            { 0,  0,  0,  0,  0,  0, -4,  6,  0,  0,  0,  0,  0,  1},
            { 0,  0,  0,  0,  0,  0,  4, -6,  0,  0,  0,  0,  0,  0},
            {-2,  0,  0,  2,  0,  0,  2, -2,  0,  0,  0,  0,  0,  0},
            { 0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0},
            { 0,  0, -1,  1,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0},
            { 0,  0,  0,  0,  1,  0,  1, -1,  0,  0,  0,  0,  0,  0},
            { 0,  0,  0,  0,  0,  0,  0, -1,  0,  5,  0,  0,  0,  2},
            { 0,  0,  0,  0,  0,  0,  0,  1, -3,  0,  0,  0,  0,  0},
            { 0,  0,  0,  0,  0,  0,  0, -1,  3,  0,  0,  0,  0,  2},
            { 0,  0,  0,  0,  0,  0,  0, -7, 12,  0,  0,  0,  0,  2},
            { 0,  0,  0,  0,  0,  0, -1,  1,  0,  0,  0,  0,  0,  2},
            { 0,  0,  0,  0,  0,  0, -1,  1,  0,  0,  0,  0,  0,  1},
            { 0,  0,  1, -1,  1,  0, -1,  0,  0,  0,  0,  0,  0,  0},
            { 0,  0,  0,  0,  0,  0,  1, -1,  0,  0,  0,  0,  0,  0},
            { 0,  0,  0,  0,  0,  0,  1, -1,  0,  0,  0,  0,  0,  1},
            { 0,  0,  1, -1,  1,  0,  1, -2,  0,  0,  0,  0,  0,  0},
            { 0,  0,  0,  0,  0,  0,  0, -2,  5,  0,  0,  0,  0,  2},
            { 0,  0,  0,  0,  0,  0,  0, -1,  0,  4,  0,  0,  0,  2},
            { 0,  0,  0,  0,  0,  0,  0,  1,  0, -4,  0,  0,  0,  0},
            { 0,  0,  0,  0,  1,  0, -1,  1,  0,  0,  0,  0,  0,  0},
            { 0,  0,  0,  0,  0,  0,  0, -6, 10,  0,  0,  0,  0,  2},
            { 0,  0,  0,  0,  0,  0,  0, -6, 10,  0,  0,  0,  0,  0},
            { 0,  0,  2, -2,  1,  0,  0, -3,  0,  3,  0,  0,  0,  0},
            { 0,  0,  0,  0,  0,  0,  0, -3,  7,  0,  0,  0,  0,  2},
            {-2,  0,  0,  2,  0,  0,  4, -4,  0,  0,  0,  0,  0,  0},
            { 0,  0,  0,  0,  0,  0,  0, -5,  8,  0,  0,  0,  0,  2},
            { 0,  0,  0,  0,  0,  0,  0,  5, -8,  0,  0,  0,  0,  0},
            { 0,  0,  0,  0,  0,  0,  0, -1,  0,  3,  0,  0,  0,  2},
            { 0,  0,  0,  0,  0,  0,  0, -1,  0,  3,  0,  0,  0,  1},
            { 0,  0,  0,  0,  0,  0,  0,  1,  0, -3,  0,  0,  0,  0},
            { 0,  0,  0,  0,  0,  0,  2, -4,  0,  0,  0,  0,  0,  0},
            { 0,  0,  0,  0,  0,  0, -2,  4,  0,  0,  0,  0,  0,  1},
            { 0,  0,  1, -1,  1,  0, -2,  3,  0,  0,  0,  0,  0,  0},
            { 0,  0,  0,  0,  0,  0, -2,  4,  0,  0,  0,  0,  0,  2},
            { 0,  0,  0,  0,  0,  0, -6,  9,  0,  0,  0,  0,  0,  2},
            { 0,  0,  0,  0,  0,  0, -6,  9,  0,  0,  0,  0,  0,  1},
            { 0,  0,  0,  0,  0,  0,  6, -9,  0,  0,  0,  0,  0,  0},
            { 0,  0,  0,  0,  1,  0,  0,  1,  0, -2,  0,  0,  0,  0},
            { 0,  0,  2, -2,  1,  0, -2,  2,  0,  0,  0,  0,  0,  0},
            { 0,  0,  0,  0,  0,  0,  0, -4,  6,  0,  0,  0,  0,  2},
            { 0,  0,  0,  0,  0,  0,  0,  4, -6,  0,  0,  0,  0,  0},
            { 0,  0,  0,  0,  1,  0,  3, -4,  0,  0,  0,  0,  0,  0},
            { 0,  0,  0,  0,  0,  0,  0, -1,  0,  2,  0,  0,  0,  2},
            { 0,  0,  0,  0,  0,  0,  0,  1,  0, -2,  0,  0,  0,  0},
            { 0,  0,  0,  0,  1,  0,  0,  1,  0, -1,  0,  0,  0,  0},
            { 0,  0,  0,  0,  0,  0, -5,  9,  0,  0,  0,  0,  0,  2},
            { 0,  0,  0,  0,  0,  0,  0,  3, -4,  0,  0,  0,  0,  0},
            { 0,  0,  0,  0,  0,  0, -3,  4,  0,  0,  0,  0,  0,  2},
            { 0,  0,  0,  0,  0,  0, -3,  4,  0,  0,  0,  0,  0,  1},
            { 0,  0,  0,  0,  0,  0,  3, -4,  0,  0,  0,  0,  0,  0},
            { 0,  0,  0,  0,  0,  0,  3, -4,  0,  0,  0,  0,  0,  1},
            { 0,  0,  0,  0,  1,  0,  0,  2, -2,  0,  0,  0,  0,  0},
            { 0,  0,  0,  0,  1,  0,  0, -1,  0,  2,  0,  0,  0,  0},
            { 0,  0,  0,  0,  0,  0,  0,  1,  0,  0, -3,  0,  0,  0},
            { 0,  0,  0,  0,  0,  0,  0,  1,  0,  1, -5,  0,  0,  0},
            { 0,  0,  0,  0,  0,  0,  0, -1,  0,  1,  0,  0,  0,  1},
            { 0,  0,  0,  0,  0,  0,  0,  1,  0, -1,  0,  0,  0,  0},
            { 0,  0,  0,  0,  0,  0,  0,  1,  0, -1,  0,  0,  0,  1},
            { 0,  0,  0,  0,  0,  0,  0,  1,  0, -3,  5,  0,  0,  0},
            { 0,  0,  0,  0,  1,  0, -3,  4,  0,  0,  0,  0,  0,  0},
            { 0,  0,  0,  0,  0,  0,  0,  1,  0,  0, -2,  0,  0,  0},
            { 0,  0,  0,  0,  0,  0,  0,  2, -2,  0,  0,  0,  0,  0},
            { 0,  0,  0,  0,  0,  0,  0,  1,  0,  0, -1,  0,  0,  0},
            { 0,  0,  0,  0,  1,  0,  0, -1,  0,  1,  0,  0,  0,  0},
            { 0,  0,  0,  0,  1,  0,  0, -2,  2,  0,  0,  0,  0,  0},
            { 0,  0,  0,  0,  0,  0, -8, 14,  0,  0,  0,  0,  0,  2},
            { 0,  0,  0,  0,  0,  0,  0,  1,  0,  2, -5,  0,  0,  0},
            { 0,  0,  0,  0,  0,  0,  0,  5, -8,  3,  0,  0,  0,  0},
            { 0,  0,  0,  0,  0,  0,  0,  5, -8,  3,  0,  0,  0,  2},
            { 0,  0,  0,  0,  0,  0,  0, -1,  0,  0,  0,  0,  0,  1},
            { 0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0},
            { 0,  0,  0,  0,  0,  0,  0,  3, -8,  3,  0,  0,  0,  0},
            { 0,  0,  0,  0,  0,  0,  0, -3,  8, -3,  0,  0,  0,  2},
            { 0,  0,  0,  0,  0,  0,  0,  1,  0, -2,  5,  0,  0,  2},
            { 0,  0,  0,  0,  0,  0, -8, 12,  0,  0,  0,  0,  0,  2},
            { 0,  0,  0,  0,  0,  0, -8, 12,  0,  0,  0,  0,  0,  0},
            { 0,  0,  0,  0,  0,  0,  0,  1,  0,  1, -2,  0,  0,  0},
            { 0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  1,  0,  0,  2},
            { 0,  0,  0,  0,  0,  0,  0,  0,  2,  0,  0,  0,  0,  0},
            { 0,  0,  0,  0,  0,  0,  0,  0,  2,  0,  0,  0,  0,  2},
            { 0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  2,  0,  0,  2},
            { 0,  0,  2, -2,  1,  0, -5,  5,  0,  0,  0,  0,  0,  0},
            { 0,  0,  0,  0,  0,  0,  0,  1,  0,  1,  0,  0,  0,  0},
            { 0,  0,  0,  0,  0,  0,  0,  1,  0,  1,  0,  0,  0,  1},
            { 0,  0,  0,  0,  0,  0,  0,  1,  0,  1,  0,  0,  0,  2},
            { 0,  0,  0,  0,  0,  0,  3, -6,  0,  0,  0,  0,  0,  0},
            { 0,  0,  0,  0,  0,  0, -3,  6,  0,  0,  0,  0,  0,  1},
            { 0,  0,  0,  0,  0,  0, -3,  6,  0,  0,  0,  0,  0,  2},
            { 0,  0,  0,  0,  0,  0,  0, -1,  4,  0,  0,  0,  0,  2},
            { 0,  0,  0,  0,  0,  0, -5,  7,  0,  0,  0,  0,  0,  2},
            { 0,  0,  0,  0,  0,  0, -5,  7,  0,  0,  0,  0,  0,  1},
            { 0,  0,  1, -1,  1,  0, -5,  6,  0,  0,  0,  0,  0,  0},
            { 0,  0,  0,  0,  0,  0,  5, -7,  0,  0,  0,  0,  0,  0},
            { 0,  0,  2, -2,  1,  0,  0, -1,  0,  1,  0,  0,  0,  0},
            { 0,  0,  0,  0,  0,  0,  0, -1,  0,  1,  0,  0,  0,  0},
            { 0,  0,  0,  0,  0, -1,  0,  3,  0,  0,  0,  0,  0,  2},
            { 0,  0,  0,  0,  0,  0,  0,  1,  0,  2,  0,  0,  0,  2},
            { 0,  0,  0,  0,  0,  0,  0, -2,  6,  0,  0,  0,  0,  2},
            { 0,  0,  0,  0,  1,  0,  2, -2,  0,  0,  0,  0,  0,  0},
            { 0,  0,  0,  0,  0,  0,  0, -6,  9,  0,  0,  0,  0,  2},
            { 0,  0,  0,  0,  0,  0,  0,  6, -9,  0,  0,  0,  0,  0},
            { 0,  0,  0,  0,  0,  0, -2,  2,  0,  0,  0,  0,  0,  1},
            { 0,  0,  1, -1,  1,  0, -2,  1,  0,  0,  0,  0,  0,  0},
            { 0,  0,  0,  0,  0,  0,  2, -2,  0,  0,  0,  0,  0,  0},
            { 0,  0,  0,  0,  0,  0,  2, -2,  0,  0,  0,  0,  0,  1},
            { 0,  0,  0,  0,  0,  0,  0,  1,  0,  3,  0,  0,  0,  2},
            { 0,  0,  0,  0,  0,  0,  0, -5,  7,  0,  0,  0,  0,  2},
            { 0,  0,  0,  0,  0,  0,  0,  5, -7,  0,  0,  0,  0,  0},
            { 0,  0,  0,  0,  1,  0, -2,  2,  0,  0,  0,  0,  0,  0},
            { 0,  0,  0,  0,  0,  0,  0,  4, -5,  0,  0,  0,  0,  0},
            { 0,  0,  0,  0,  0,  0,  1, -3,  0,  0,  0,  0,  0,  0},
            { 0,  0,  0,  0,  0,  0, -1,  3,  0,  0,  0,  0,  0,  1},
            { 0,  0,  1, -1,  1,  0, -1,  2,  0,  0,  0,  0,  0,  0},
            { 0,  0,  0,  0,  0,  0, -1,  3,  0,  0,  0,  0,  0,  2},
            { 0,  0,  0,  0,  0,  0, -7, 10,  0,  0,  0,  0,  0,  2},
            { 0,  0,  0,  0,  0,  0, -7, 10,  0,  0,  0,  0,  0,  1},
            { 0,  0,  0,  0,  0,  0,  0,  3, -3,  0,  0,  0,  0,  0},
            { 0,  0,  0,  0,  0,  0, -4,  8,  0,  0,  0,  0,  0,  2},
            { 0,  0,  0,  0,  0,  0, -4,  5,  0,  0,  0,  0,  0,  2},
            { 0,  0,  0,  0,  0,  0, -4,  5,  0,  0,  0,  0,  0,  1},
            { 0,  0,  0,  0,  0,  0,  4, -5,  0,  0,  0,  0,  0,  0},
            { 0,  0,  0,  0,  0,  0,  0,  1,  1,  0,  0,  0,  0,  2},
            { 0,  0,  0,  0,  0,  0,  0, -2,  0,  5,  0,  0,  0,  2},
            { 0,  0,  0,  0,  0,  0,  0,  0,  3,  0,  0,  0,  0,  2},
            { 0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0},
            { 0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  2},
            { 0,  0,  0,  0,  0,  0, -9, 13,  0,  0,  0,  0,  0,  2},
            { 0,  0,  0,  0,  0,  0,  0, -1,  5,  0,  0,  0,  0,  2},
            { 0,  0,  0,  0,  0,  0,  0, -2,  0,  4,  0,  0,  0,  2},
            { 0,  0,  0,  0,  0,  0,  0,  2,  0, -4,  0,  0,  0,  0},
            { 0,  0,  0,  0,  0,  0,  0, -2,  7,  0,  0,  0,  0,  2},
            { 0,  0,  0,  0,  0,  0,  0,  2,  0, -3,  0,  0,  0,  0},
            { 0,  0,  0,  0,  0,  0, -2,  5,  0,  0,  0,  0,  0,  1},
            { 0,  0,  0,  0,  0,  0, -2,  5,  0,  0,  0,  0,  0,  2},
            { 0,  0,  0,  0,  0,  0, -6,  8,  0,  0,  0,  0,  0,  2},
            { 0,  0,  0,  0,  0,  0, -6,  8,  0,  0,  0,  0,  0,  1},
            { 0,  0,  0,  0,  0,  0,  6, -8,  0,  0,  0,  0,  0,  0},
            { 0,  0,  0,  0,  1,  0,  0,  2,  0, -2,  0,  0,  0,  0},
            { 0,  0,  0,  0,  0,  0,  0, -3,  9,  0,  0,  0,  0,  2},
            { 0,  0,  0,  0,  0,  0,  0,  5, -6,  0,  0,  0,  0,  0},
            { 0,  0,  0,  0,  0,  0,  0,  5, -6,  0,  0,  0,  0,  2},
            { 0,  0,  0,  0,  0,  0,  0,  2,  0, -2,  0,  0,  0,  0},
            { 0,  0,  0,  0,  0,  0,  0,  2,  0, -2,  0,  0,  0,  1},
            { 0,  0,  0,  0,  0,  0,  0,  2,  0, -2,  0,  0,  0,  2},
            { 0,  0,  0,  0,  0,  0, -5, 10,  0,  0,  0,  0,  0,  2},
            { 0,  0,  0,  0,  0,  0,  0,  4, -4,  0,  0,  0,  0,  0},
            { 0,  0,  0,  0,  0,  0,  0,  4, -4,  0,  0,  0,  0,  2},
            { 0,  0,  0,  0,  0,  0, -3,  3,  0,  0,  0,  0,  0,  1},
            { 0,  0,  0,  0,  0,  0,  3, -3,  0,  0,  0,  0,  0,  0},
            { 0,  0,  0,  0,  0,  0,  3, -3,  0,  0,  0,  0,  0,  1},
            { 0,  0,  0,  0,  0,  0,  3, -3,  0,  0,  0,  0,  0,  2},
            { 0,  0,  0,  0,  0,  0,  0,  2,  0,  0, -3,  0,  0,  0},
            { 0,  0,  0,  0,  0,  0,  0, -5, 13,  0,  0,  0,  0,  2},
            { 0,  0,  0,  0,  0,  0,  0,  2,  0, -1,  0,  0,  0,  0},
            { 0,  0,  0,  0,  0,  0,  0,  2,  0, -1,  0,  0,  0,  2},
            { 0,  0,  0,  0,  0,  0,  0,  2,  0,  0, -2,  0,  0,  0},
            { 0,  0,  0,  0,  0,  0,  0,  2,  0,  0, -2,  0,  0,  1},
            { 0,  0,  0,  0,  0,  0,  0,  3, -2,  0,  0,  0,  0,  0},
            { 0,  0,  0,  0,  0,  0,  0,  3, -2,  0,  0,  0,  0,  2},
            { 0,  0,  0,  0,  0,  0,  0,  2,  0,  0, -1,  0,  0,  2},
            { 0,  0,  0,  0,  0,  0,  0, -6, 15,  0,  0,  0,  0,  2},
            { 0,  0,  0,  0,  0,  0, -8, 15,  0,  0,  0,  0,  0,  2},
            { 0,  0,  0,  0,  0,  0, -3,  9, -4,  0,  0,  0,  0,  2},
            { 0,  0,  0,  0,  0,  0,  0,  2,  0,  2, -5,  0,  0,  2},
            { 0,  0,  0,  0,  0,  0,  0, -2,  8, -1, -5,  0,  0,  2},
            { 0,  0,  0,  0,  0,  0,  0,  6, -8,  3,  0,  0,  0,  2},
            { 0,  0,  0,  0,  0,  0,  0,  2,  0,  0,  0,  0,  0,  0},
            { 0,  0,  0,  0,  0,  0,  0,  2,  0,  0,  0,  0,  0,  0},
            { 0,  0,  0,  0,  0,  0,  0,  2,  0,  0,  0,  0,  0,  1},
            { 0,  0,  1, -1,  1,  0,  0,  1,  0,  0,  0,  0,  0,  0},
            { 0,  0,  0,  0,  0,  0,  0,  2,  0,  0,  0,  0,  0,  1},
            { 0,  0,  0,  0,  0,  0,  0,  2,  0,  0,  0,  0,  0,  2},
            { 0,  0,  0,  0,  0,  0,  0, -6, 16, -4, -5,  0,  0,  2},
            { 0,  0,  0,  0,  0,  0,  0, -2,  8, -3,  0,  0,  0,  2},
            { 0,  0,  0,  0,  0,  0,  0, -2,  8, -3,  0,  0,  0,  2},
            { 0,  0,  0,  0,  0,  0,  0,  6, -8,  1,  5,  0,  0,  2},
            { 0,  0,  0,  0,  0,  0,  0,  2,  0, -2,  5,  0,  0,  2},
            { 0,  0,  0,  0,  0,  0,  3, -5,  4,  0,  0,  0,  0,  2},
            { 0,  0,  0,  0,  0,  0, -8, 11,  0,  0,  0,  0,  0,  2},
            { 0,  0,  0,  0,  0,  0, -8, 11,  0,  0,  0,  0,  0,  1},
            { 0,  0,  0,  0,  0,  0, -8, 11,  0,  0,  0,  0,  0,  2},
            { 0,  0,  0,  0,  0,  0,  0, 11,  0,  0,  0,  0,  0,  2},
            { 0,  0,  0,  0,  0,  0,  0,  2,  0,  0,  1,  0,  0,  2},
            { 0,  0,  0,  0,  0,  0,  3, -3,  0,  2,  0,  0,  0,  2},
            { 0,  0,  2, -2,  1,  0,  0,  4, -8,  3,  0,  0,  0,  0},
            { 0,  0,  1, -1,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0},
            { 0,  0,  2, -2,  1,  0,  0, -4,  8, -3,  0,  0,  0,  0},
            { 0,  0,  0,  0,  0,  0,  0,  1,  2,  0,  0,  0,  0,  2},
            { 0,  0,  0,  0,  0,  0,  0,  2,  0,  1,  0,  0,  0,  2},
            { 0,  0,  0,  0,  0,  0, -3,  7,  0,  0,  0,  0,  0,  2},
            { 0,  0,  0,  0,  0,  0,  0,  0,  4,  0,  0,  0,  0,  2},
            { 0,  0,  0,  0,  0,  0, -5,  6,  0,  0,  0,  0,  0,  2},
            { 0,  0,  0,  0,  0,  0, -5,  6,  0,  0,  0,  0,  0,  1},
            { 0,  0,  0,  0,  0,  0,  5, -6,  0,  0,  0,  0,  0,  0},
            { 0,  0,  0,  0,  0,  0,  5, -6,  0,  0,  0,  0,  0,  2},
            { 0,  0,  0,  0,  0,  0,  0,  2,  0,  2,  0,  0,  0,  2},
            { 0,  0,  0,  0,  0,  0,  0, -1,  6,  0,  0,  0,  0,  2},
            { 0,  0,  0,  0,  0,  0,  0,  7, -9,  0,  0,  0,  0,  2},
            { 0,  0,  0,  0,  0,  0,  2, -1,  0,  0,  0,  0,  0,  0},
            { 0,  0,  0,  0,  0,  0,  2, -1,  0,  0,  0,  0,  0,  2},
            { 0,  0,  0,  0,  0,  0,  0,  6, -7,  0,  0,  0,  0,  2},
            { 0,  0,  0,  0,  0,  0,  0,  5, -5,  0,  0,  0,  0,  2},
            { 0,  0,  0,  0,  0,  0, -1,  4,  0,  0,  0,  0,  0,  1},
            { 0,  0,  0,  0,  0,  0, -1,  4,  0,  0,  0,  0,  0,  2},
            { 0,  0,  0,  0,  0,  0, -7,  9,  0,  0,  0,  0,  0,  2},
            { 0,  0,  0,  0,  0,  0, -7,  9,  0,  0,  0,  0,  0,  1},
            { 0,  0,  0,  0,  0,  0,  0,  4, -3,  0,  0,  0,  0,  2},
            { 0,  0,  0,  0,  0,  0,  0,  3, -1,  0,  0,  0,  0,  2},
            { 0,  0,  0,  0,  0,  0, -4,  4,  0,  0,  0,  0,  0,  1},
            { 0,  0,  0,  0,  0,  0,  4, -4,  0,  0,  0,  0,  0,  0},
            { 0,  0,  0,  0,  0,  0,  4, -4,  0,  0,  0,  0,  0,  1},
            { 0,  0,  0,  0,  0,  0,  4, -4,  0,  0,  0,  0,  0,  2},
            { 0,  0,  0,  0,  0,  0,  0,  2,  1,  0,  0,  0,  0,  2},
            { 0,  0,  0,  0,  0,  0,  0, -3,  0,  5,  0,  0,  0,  2},
            { 0,  0,  0,  0,  0,  0,  1,  1,  0,  0,  0,  0,  0,  0},
            { 0,  0,  0,  0,  0,  0,  1,  1,  0,  0,  0,  0,  0,  1},
            { 0,  0,  0,  0,  0,  0,  1,  1,  0,  0,  0,  0,  0,  2},
            { 0,  0,  0,  0,  0,  0, -9, 12,  0,  0,  0,  0,  0,  2},
            { 0,  0,  0,  0,  0,  0,  0,  3,  0, -4,  0,  0,  0,  0},
            { 0,  0,  2, -2,  1,  0,  1, -1,  0,  0,  0,  0,  0,  0},
            { 0,  0,  0,  0,  0,  0,  0,  7, -8,  0,  0,  0,  0,  2},
            { 0,  0,  0,  0,  0,  0,  0,  3,  0, -3,  0,  0,  0,  0},
            { 0,  0,  0,  0,  0,  0,  0,  3,  0, -3,  0,  0,  0,  2},
            { 0,  0,  0,  0,  0,  0, -2,  6,  0,  0,  0,  0,  0,  2},
            { 0,  0,  0,  0,  0,  0, -6,  7,  0,  0,  0,  0,  0,  1},
            { 0,  0,  0,  0,  0,  0,  6, -7,  0,  0,  0,  0,  0,  0},
            { 0,  0,  0,  0,  0,  0,  0,  6, -6,  0,  0,  0,  0,  2},
            { 0,  0,  0,  0,  0,  0,  0,  3,  0, -2,  0,  0,  0,  0},
            { 0,  0,  0,  0,  0,  0,  0,  3,  0, -2,  0,  0,  0,  2},
            { 0,  0,  0,  0,  0,  0,  0,  5, -4,  0,  0,  0,  0,  2},
            { 0,  0,  0,  0,  0,  0,  3, -2,  0,  0,  0,  0,  0,  0},
            { 0,  0,  0,  0,  0,  0,  3, -2,  0,  0,  0,  0,  0,  2},
            { 0,  0,  0,  0,  0,  0,  0,  3,  0, -1,  0,  0,  0,  2},
            { 0,  0,  0,  0,  0,  0,  0,  3,  0, -1,  0,  0,  0,  2},
            { 0,  0,  0,  0,  0,  0,  0,  3,  0,  0, -2,  0,  0,  2},
            { 0,  0,  0,  0,  0,  0,  0,  4, -2,  0,  0,  0,  0,  2},
            { 0,  0,  0,  0,  0,  0,  0,  3,  0,  0, -1,  0,  0,  2},
            { 0,  0,  2, -2,  1,  0,  0,  1,  0, -1,  0,  0,  0,  0},
            { 0,  0,  0,  0,  0,  0, -8, 16,  0,  0,  0,  0,  0,  2},
            { 0,  0,  0,  0,  0,  0,  0,  3,  0,  2, -5,  0,  0,  2},
            { 0,  0,  0,  0,  0,  0,  0,  7, -8,  3,  0,  0,  0,  2},
            { 0,  0,  0,  0,  0,  0,  0, -5, 16, -4, -5,  0,  0,  2},
            { 0,  0,  0,  0,  0,  0,  0,  3,  0,  0,  0,  0,  0,  2},
            { 0,  0,  0,  0,  0,  0,  0, -1,  8, -3,  0,  0,  0,  2},
            { 0,  0,  0,  0,  0,  0, -8, 10,  0,  0,  0,  0,  0,  2},
            { 0,  0,  0,  0,  0,  0, -8, 10,  0,  0,  0,  0,  0,  1},
            { 0,  0,  0,  0,  0,  0, -8, 10,  0,  0,  0,  0,  0,  2},
            { 0,  0,  0,  0,  0,  0,  0,  2,  2,  0,  0,  0,  0,  2},
            { 0,  0,  0,  0,  0,  0,  0,  3,  0,  1,  0,  0,  0,  2},
            { 0,  0,  0,  0,  0,  0, -3,  8,  0,  0,  0,  0,  0,  2},
            { 0,  0,  0,  0,  0,  0, -5,  5,  0,  0,  0,  0,  0,  1},
            { 0,  0,  0,  0,  0,  0,  5, -5,  0,  0,  0,  0,  0,  0},
            { 0,  0,  0,  0,  0,  0,  5, -5,  0,  0,  0,  0,  0,  1},
            { 0,  0,  0,  0,  0,  0,  5, -5,  0,  0,  0,  0,  0,  2},
            { 0,  0,  0,  0,  0,  0,  2,  0,  0,  0,  0,  0,  0,  0},
            { 0,  0,  0,  0,  0,  0,  2,  0,  0,  0,  0,  0,  0,  1},
            { 0,  0,  0,  0,  0,  0,  2,  0,  0,  0,  0,  0,  0,  2},
            { 0,  0,  0,  0,  0,  0,  0,  7, -7,  0,  0,  0,  0,  2},
            { 0,  0,  0,  0,  0,  0,  0,  7, -7,  0,  0,  0,  0,  2},
            { 0,  0,  0,  0,  0,  0,  0,  6, -5,  0,  0,  0,  0,  2},
            { 0,  0,  0,  0,  0,  0,  7, -8,  0,  0,  0,  0,  0,  0},
            { 0,  0,  0,  0,  0,  0,  0,  5, -3,  0,  0,  0,  0,  2},
            { 0,  0,  0,  0,  0,  0,  4, -3,  0,  0,  0,  0,  0,  2},
            { 0,  0,  0,  0,  0,  0,  1,  2,  0,  0,  0,  0,  0,  2},
            { 0,  0,  0,  0,  0,  0, -9, 11,  0,  0,  0,  0,  0,  2},
            { 0,  0,  0,  0,  0,  0, -9, 11,  0,  0,  0,  0,  0,  1},
            { 0,  0,  0,  0,  0,  0,  0,  4,  0, -4,  0,  0,  0,  2},
            { 0,  0,  0,  0,  0,  0,  0,  4,  0, -3,  0,  0,  0,  2},
            { 0,  0,  0,  0,  0,  0, -6,  6,  0,  0,  0,  0,  0,  1},
            { 0,  0,  0,  0,  0,  0,  6, -6,  0,  0,  0,  0,  0,  0},
            { 0,  0,  0,  0,  0,  0,  6, -6,  0,  0,  0,  0,  0,  1},
            { 0,  0,  0,  0,  0,  0,  0,  4,  0, -2,  0,  0,  0,  2},
            { 0,  0,  0,  0,  0,  0,  0,  6, -4,  0,  0,  0,  0,  2},
            { 0,  0,  0,  0,  0,  0,  3, -1,  0,  0,  0,  0,  0,  0},
            { 0,  0,  0,  0,  0,  0,  3, -1,  0,  0,  0,  0,  0,  1},
            { 0,  0,  0,  0,  0,  0,  3, -1,  0,  0,  0,  0,  0,  2},
            { 0,  0,  0,  0,  0,  0,  0,  4,  0, -1,  0,  0,  0,  2},
            { 0,  0,  0,  0,  0,  0,  0,  4,  0,  0, -2,  0,  0,  2},
            { 0,  0,  0,  0,  0,  0,  0,  5, -2,  0,  0,  0,  0,  2},
            { 0,  0,  0,  0,  0,  0,  0,  4,  0,  0,  0,  0,  0,  0},
            { 0,  0,  0,  0,  0,  0,  8, -9,  0,  0,  0,  0,  0,  0},
            { 0,  0,  0,  0,  0,  0,  5, -4,  0,  0,  0,  0,  0,  2},
            { 0,  0,  0,  0,  0,  0,  2,  1,  0,  0,  0,  0,  0,  2},
            { 0,  0,  0,  0,  0,  0,  2,  1,  0,  0,  0,  0,  0,  1},
            { 0,  0,  0,  0,  0,  0,  2,  1,  0,  0,  0,  0,  0,  1},
            { 0,  0,  0,  0,  0,  0, -7,  7,  0,  0,  0,  0,  0,  1},
            { 0,  0,  0,  0,  0,  0,  7, -7,  0,  0,  0,  0,  0,  0},
            { 0,  0,  0,  0,  0,  0,  4, -2,  0,  0,  0,  0,  0,  1},
            { 0,  0,  0,  0,  0,  0,  4, -2,  0,  0,  0,  0,  0,  2},
            { 0,  0,  0,  0,  0,  0,  4, -2,  0,  0,  0,  0,  0,  0},
            { 0,  0,  0,  0,  0,  0,  4, -2,  0,  0,  0,  0,  0,  0},
            { 0,  0,  0,  0,  0,  0,  0,  5,  0, -4,  0,  0,  0,  2},
            { 0,  0,  0,  0,  0,  0,  0,  5,  0, -3,  0,  0,  0,  2},
            { 0,  0,  0,  0,  0,  0,  0,  5,  0, -2,  0,  0,  0,  2},
            { 0,  0,  0,  0,  0,  0,  3,  0,  0,  0,  0,  0,  0,  2},
            { 0,  0,  0,  0,  0,  0, -8,  8,  0,  0,  0,  0,  0,  1},
            { 0,  0,  0,  0,  0,  0,  8, -8,  0,  0,  0,  0,  0,  0},
            { 0,  0,  0,  0,  0,  0,  5, -3,  0,  0,  0,  0,  0,  1},
            { 0,  0,  0,  0,  0,  0,  5, -3,  0,  0,  0,  0,  0,  2},
            { 0,  0,  0,  0,  0,  0, -9,  9,  0,  0,  0,  0,  0,  1},
            { 0,  0,  0,  0,  0,  0, -9,  9,  0,  0,  0,  0,  0,  1},
            { 0,  0,  0,  0,  0,  0, -9,  9,  0,  0,  0,  0,  0,  1},
            { 0,  0,  0,  0,  0,  0,  9, -9,  0,  0,  0,  0,  0,  0},
            { 0,  0,  0,  0,  0,  0,  6, -4,  0,  0,  0,  0,  0,  1},
            { 0,  0,  0,  0,  0,  0,  0,  6,  0,  0,  0,  0,  0,  2},
            { 0,  0,  0,  0,  0,  0,  0,  6,  0,  0,  0,  0,  0,  0},
            { 0,  0,  0,  0,  0,  0,  0,  6,  0,  0,  0,  0,  0,  0},
            { 0,  0,  0,  0,  0,  0,  0,  6,  0,  0,  0,  0,  0,  1},
            { 0,  0,  0,  0,  0,  0,  0,  6,  0,  0,  0,  0,  0,  2},
            { 0,  0,  0,  0,  0,  0,  0,  6,  0,  0,  0,  0,  0,  0},
            { 0,  0,  0,  0,  0,  0,  0,  6,  0,  0,  0,  0,  0,  1},
            { 0,  0,  0,  0,  0,  0,  0,  6,  0,  0,  0,  0,  0,  2},
            { 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  2},
            { 1,  0,  0, -2,  0,  0,  0,  2,  0, -2,  0,  0,  0,  0},
            { 1,  0,  0, -2,  0,  0,  2, -2,  0,  0,  0,  0,  0,  0},
            { 1,  0,  0, -2,  0,  0,  0,  1,  0, -1,  0,  0,  0,  0},
            { 1,  0,  0, -2,  0,  0,  1, -1,  0,  0,  0,  0,  0,  0},
            {-1,  0,  0,  0,  0,  0,  3, -3,  0,  0,  0,  0,  0,  0},
            {-1,  0,  0,  0,  0,  0,  0,  2,  0, -2,  0,  0,  0,  0},
            {-1,  0,  0,  2,  0,  0,  0,  4, -8,  3,  0,  0,  0,  0},
            { 1,  0,  0, -2,  0,  0,  0,  4, -8,  3,  0,  0,  0,  0},
            {-2,  0,  0,  2,  0,  0,  0,  4, -8,  3,  0,  0,  0,  0},
            {-1,  0,  0,  0,  0,  0,  0,  2,  0, -3,  0,  0,  0,  0},
            {-1,  0,  0,  0,  0,  0,  0,  1,  0, -1,  0,  0,  0,  0},
            {-1,  0,  0,  0,  0,  0,  1, -1,  0,  0,  0,  0,  0,  0},
            {-1,  0,  0,  2,  0,  0,  2, -2,  0,  0,  0,  0,  0,  0},
            { 1,  0, -1,  1,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0},
            {-1,  0,  0,  2,  0,  0,  0,  2,  0, -3,  0,  0,  0,  0},
            {-2,  0,  0,  0,  0,  0,  0,  2,  0, -3,  0,  0,  0,  0},
            { 1,  0,  0,  0,  0,  0,  0,  4, -8,  3,  0,  0,  0,  0},
            {-1,  0,  1, -1,  1,  0,  0, -1,  0,  0,  0,  0,  0,  0},
            { 1,  0,  1, -1,  1,  0,  0, -1,  0,  0,  0,  0,  0,  0},
            {-1,  0,  0,  0,  0,  0,  0,  4, -8,  3,  0,  0,  0,  0},
            {-1,  0,  0,  2,  1,  0,  0,  2,  0, -2,  0,  0,  0,  0},
            { 0,  0,  0,  0,  0,  0,  0,  2,  0, -2,  0,  0,  0,  0},
            {-1,  0,  0,  2,  0,  0,  0,  2,  0, -2,  0,  0,  0,  0},
            {-1,  0,  0,  2,  0,  0,  3, -3,  0,  0,  0,  0,  0,  0},
            { 1,  0,  0, -2,  1,  0,  0, -2,  0,  2,  0,  0,  0,  0},
            { 1,  0,  2, -2,  2,  0, -3,  3,  0,  0,  0,  0,  0,  0},
            { 1,  0,  2, -2,  2,  0,  0, -2,  0,  2,  0,  0,  0,  0},
            { 1,  0,  0,  0,  0,  0,  1, -1,  0,  0,  0,  0,  0,  0},
            { 1,  0,  0,  0,  0,  0,  0,  1,  0, -1,  0,  0,  0,  0},
            { 0,  0,  0, -2,  0,  0,  2, -2,  0,  0,  0,  0,  0,  0},
            { 0,  0,  0, -2,  0,  0,  0,  1,  0, -1,  0,  0,  0,  0},
            { 0,  0,  2,  0,  2,  0, -2,  2,  0,  0,  0,  0,  0,  0},
            { 0,  0,  2,  0,  2,  0,  0, -1,  0,  1,  0,  0,  0,  0},
            { 0,  0,  2,  0,  2,  0, -1,  1,  0,  0,  0,  0,  0,  0},
            { 0,  0,  2,  0,  2,  0, -2,  3,  0,  0,  0,  0,  0,  0},
            { 0,  0,  0,  2,  0,  0,  0,  2,  0, -2,  0,  0,  0,  0},
            { 0,  0,  1,  1,  2,  0,  0,  1,  0,  0,  0,  0,  0,  0},
            { 1,  0,  2,  0,  2,  0,  0,  1,  0,  0,  0,  0,  0,  0},
            {-1,  0,  2,  0,  2,  0, 10, -3,  0,  0,  0,  0,  0,  0},
            { 0,  0,  1,  1,  1,  0,  0,  1,  0,  0,  0,  0,  0,  0},
            { 1,  0,  2,  0,  2,  0,  0,  1,  0,  0,  0,  0,  0,  0},
            { 0,  0,  2,  0,  2,  0,  0,  4, -8,  3,  0,  0,  0,  0},
            { 0,  0,  2,  0,  2,  0,  0, -4,  8, -3,  0,  0,  0,  0},
            {-1,  0,  2,  0,  2,  0,  0, -4,  8, -3,  0,  0,  0,  0},
            { 2,  0,  2, -2,  2,  0,  0, -2,  0,  3,  0,  0,  0,  0},
            { 1,  0,  2,  0,  1,  0,  0, -2,  0,  3,  0,  0,  0,  0},
            { 0,  0,  1,  1,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0},
            {-1,  0,  2,  0,  1,  0,  0,  1,  0,  0,  0,  0,  0,  0},
            {-2,  0,  2,  2,  2,  0,  0,  2,  0, -2,  0,  0,  0,  0},
            { 0,  0,  2,  0,  2,  0,  2, -3,  0,  0,  0,  0,  0,  0},
            { 0,  0,  2,  0,  2,  0,  1, -1,  0,  0,  0,  0,  0,  0},
            { 0,  0,  2,  0,  2,  0,  0,  1,  0, -1,  0,  0,  0,  0},
            { 0,  0,  2,  0,  2,  0,  2, -2,  0,  0,  0,  0,  0,  0},
            {-1,  0,  2,  2,  2,  0,  0, -1,  0,  1,  0,  0,  0,  0},
            { 1,  0,  2,  0,  2,  0, -1,  1,  0,  0,  0,  0,  0,  0},
            {-1,  0,  2,  2,  2,  0,  0,  2,  0, -3,  0,  0,  0,  0},
            { 2,  0,  2,  0,  2,  0,  0,  2,  0, -3,  0,  0,  0,  0},
            { 1,  0,  2,  0,  2,  0,  0, -4,  8, -3,  0,  0,  0,  0},
            { 1,  0,  2,  0,  2,  0,  0,  4, -8,  3,  0,  0,  0,  0},
            { 1,  0,  1,  1,  1,  0,  0,  1,  0,  0,  0,  0,  0,  0},
            { 0,  0,  2,  0,  2,  0,  0,  1,  0,  0,  0,  0,  0,  0},
            { 2,  0,  2,  0,  1,  0,  0,  1,  0,  0,  0,  0,  0,  0},
            {-1,  0,  2,  2,  2,  0,  0,  2,  0, -2,  0,  0,  0,  0},
            {-1,  0,  2,  2,  2,  0,  3, -3,  0,  0,  0,  0,  0,  0},
            { 1,  0,  2,  0,  2,  0,  1, -1,  0,  0,  0,  0,  0,  0},
            { 0,  0,  2,  2,  2,  0,  0,  2,  0, -2,  0,  0,  0,  0}};
        
        /*
         Planetary nutation coefficients, unit 1e-7 arcsec:
         longitude (sin, cos), obliquity (sin, cos)
         
         Each row of coefficients in 'cpl_t' belongs with the corresponding
         row of fundamental-argument multipliers in 'napl_t'.
         */
        
        static const double cpl_t[687][4] = {
            { 1440.0,          0.0,          0.0,          0.0},
            {   56.0,       -117.0,        -42.0,        -40.0},
            {  125.0,        -43.0,          0.0,        -54.0},
            {    0.0,          5.0,          0.0,          0.0},
            {    3.0,         -7.0,         -3.0,          0.0},
            {    3.0,          0.0,          0.0,         -2.0},
            { -114.0,          0.0,          0.0,         61.0},
            { -219.0,         89.0,          0.0,          0.0},
            {   -3.0,          0.0,          0.0,          0.0},
            { -462.0,       1604.0,          0.0,          0.0},
            {   99.0,          0.0,          0.0,        -53.0},
            {   -3.0,          0.0,          0.0,          2.0},
            {    0.0,          6.0,          2.0,          0.0},
            {    3.0,          0.0,          0.0,          0.0},
            {  -12.0,          0.0,          0.0,          0.0},
            {   14.0,       -218.0,        117.0,          8.0},
            {   31.0,       -481.0,       -257.0,        -17.0},
            { -491.0,        128.0,          0.0,          0.0},
            {-3084.0,       5123.0,       2735.0,       1647.0},
            {-1444.0,       2409.0,      -1286.0,       -771.0},
            {   11.0,        -24.0,        -11.0,         -9.0},
            {   26.0,         -9.0,          0.0,          0.0},
            {  103.0,        -60.0,          0.0,          0.0},
            {    0.0,        -13.0,         -7.0,          0.0},
            {  -26.0,        -29.0,        -16.0,         14.0},
            {    9.0,        -27.0,        -14.0,         -5.0},
            {   12.0,          0.0,          0.0,         -6.0},
            {   -7.0,          0.0,          0.0,          0.0},
            {    0.0,         24.0,          0.0,          0.0},
            {  284.0,          0.0,          0.0,       -151.0},
            {  226.0,        101.0,          0.0,          0.0},
            {    0.0,         -8.0,         -2.0,          0.0},
            {    0.0,         -6.0,         -3.0,          0.0},
            {    5.0,          0.0,          0.0,         -3.0},
            {  -41.0,        175.0,         76.0,         17.0},
            {    0.0,         15.0,          6.0,          0.0},
            {  425.0,        212.0,       -133.0,        269.0},
            { 1200.0,        598.0,        319.0,       -641.0},
            {  235.0,        334.0,          0.0,          0.0},
            {   11.0,        -12.0,         -7.0,         -6.0},
            {    5.0,         -6.0,          3.0,          3.0},
            {   -5.0,          0.0,          0.0,          3.0},
            {    6.0,          0.0,          0.0,         -3.0},
            {   15.0,          0.0,          0.0,          0.0},
            {   13.0,          0.0,          0.0,         -7.0},
            {   -6.0,         -9.0,          0.0,          0.0},
            {  266.0,        -78.0,          0.0,          0.0},
            { -460.0,       -435.0,       -232.0,        246.0},
            {    0.0,         15.0,          7.0,          0.0},
            {   -3.0,          0.0,          0.0,          2.0},
            {    0.0,        131.0,          0.0,          0.0},
            {    4.0,          0.0,          0.0,          0.0},
            {    0.0,          3.0,          0.0,          0.0},
            {    0.0,          4.0,          2.0,          0.0},
            {    0.0,          3.0,          0.0,          0.0},
            {  -17.0,        -19.0,        -10.0,          9.0},
            {   -9.0,        -11.0,          6.0,         -5.0},
            {   -6.0,          0.0,          0.0,          3.0},
            {  -16.0,          8.0,          0.0,          0.0},
            {    0.0,          3.0,          0.0,          0.0},
            {   11.0,         24.0,         11.0,         -5.0},
            {   -3.0,         -4.0,         -2.0,          1.0},
            {    3.0,          0.0,          0.0,         -1.0},
            {    0.0,         -8.0,         -4.0,          0.0},
            {    0.0,          3.0,          0.0,          0.0},
            {    0.0,          5.0,          0.0,          0.0},
            {    0.0,          3.0,          2.0,          0.0},
            {   -6.0,          4.0,          2.0,          3.0},
            {   -3.0,         -5.0,          0.0,          0.0},
            {   -5.0,          0.0,          0.0,          2.0},
            {    4.0,         24.0,         13.0,         -2.0},
            {  -42.0,         20.0,          0.0,          0.0},
            {  -10.0,        233.0,          0.0,          0.0},
            {   -3.0,          0.0,          0.0,          1.0},
            {   78.0,        -18.0,          0.0,          0.0},
            {    0.0,          3.0,          1.0,          0.0},
            {    0.0,         -3.0,         -1.0,          0.0},
            {    0.0,         -4.0,         -2.0,          1.0},
            {    0.0,         -8.0,         -4.0,         -1.0},
            {    0.0,         -5.0,          3.0,          0.0},
            {   -7.0,          0.0,          0.0,          3.0},
            {  -14.0,          8.0,          3.0,          6.0},
            {    0.0,          8.0,         -4.0,          0.0},
            {    0.0,         19.0,         10.0,          0.0},
            {   45.0,        -22.0,          0.0,          0.0},
            {   -3.0,          0.0,          0.0,          0.0},
            {    0.0,         -3.0,          0.0,          0.0},
            {    0.0,          3.0,          0.0,          0.0},
            {    3.0,          5.0,          3.0,         -2.0},
            {   89.0,        -16.0,         -9.0,        -48.0},
            {    0.0,          3.0,          0.0,          0.0},
            {   -3.0,          7.0,          4.0,          2.0},
            { -349.0,        -62.0,          0.0,          0.0},
            {  -15.0,         22.0,          0.0,          0.0},
            {   -3.0,          0.0,          0.0,          0.0},
            {  -53.0,          0.0,          0.0,          0.0},
            {    5.0,          0.0,          0.0,         -3.0},
            {    0.0,         -8.0,          0.0,          0.0},
            {   15.0,         -7.0,         -4.0,         -8.0},
            {   -3.0,          0.0,          0.0,          1.0},
            {  -21.0,        -78.0,          0.0,          0.0},
            {   20.0,        -70.0,        -37.0,        -11.0},
            {    0.0,          6.0,          3.0,          0.0},
            {    5.0,          3.0,          2.0,         -2.0},
            {  -17.0,         -4.0,         -2.0,          9.0},
            {    0.0,          6.0,          3.0,          0.0},
            {   32.0,         15.0,         -8.0,         17.0},
            {  174.0,         84.0,         45.0,        -93.0},
            {   11.0,         56.0,          0.0,          0.0},
            {  -66.0,        -12.0,         -6.0,         35.0},
            {   47.0,          8.0,          4.0,        -25.0},
            {    0.0,          8.0,          4.0,          0.0},
            {   10.0,        -22.0,        -12.0,         -5.0},
            {   -3.0,          0.0,          0.0,          2.0},
            {  -24.0,         12.0,          0.0,          0.0},
            {    5.0,         -6.0,          0.0,          0.0},
            {    3.0,          0.0,          0.0,         -2.0},
            {    4.0,          3.0,          1.0,         -2.0},
            {    0.0,         29.0,         15.0,          0.0},
            {   -5.0,         -4.0,         -2.0,          2.0},
            {    8.0,         -3.0,         -1.0,         -5.0},
            {    0.0,         -3.0,          0.0,          0.0},
            {   10.0,          0.0,          0.0,          0.0},
            {    3.0,          0.0,          0.0,         -2.0},
            {   -5.0,          0.0,          0.0,          3.0},
            {   46.0,         66.0,         35.0,        -25.0},
            {  -14.0,          7.0,          0.0,          0.0},
            {    0.0,          3.0,          2.0,          0.0},
            {   -5.0,          0.0,          0.0,          0.0},
            {  -68.0,        -34.0,        -18.0,         36.0},
            {    0.0,         14.0,          7.0,          0.0},
            {   10.0,         -6.0,         -3.0,         -5.0},
            {   -5.0,         -4.0,         -2.0,          3.0},
            {   -3.0,          5.0,          2.0,          1.0},
            {   76.0,         17.0,          9.0,        -41.0},
            {   84.0,        298.0,        159.0,        -45.0},
            {    3.0,          0.0,          0.0,         -1.0},
            {   -3.0,          0.0,          0.0,          2.0},
            {   -3.0,          0.0,          0.0,          1.0},
            {  -82.0,        292.0,        156.0,         44.0},
            {  -73.0,         17.0,          9.0,         39.0},
            {   -9.0,        -16.0,          0.0,          0.0},
            {    3.0,          0.0,         -1.0,         -2.0},
            {   -3.0,          0.0,          0.0,          0.0},
            {   -9.0,         -5.0,         -3.0,          5.0},
            { -439.0,          0.0,          0.0,          0.0},
            {   57.0,        -28.0,        -15.0,        -30.0},
            {    0.0,         -6.0,         -3.0,          0.0},
            {   -4.0,          0.0,          0.0,          2.0},
            {  -40.0,         57.0,         30.0,         21.0},
            {   23.0,          7.0,          3.0,        -13.0},
            {  273.0,         80.0,         43.0,       -146.0},
            { -449.0,        430.0,          0.0,          0.0},
            {   -8.0,        -47.0,        -25.0,          4.0},
            {    6.0,         47.0,         25.0,         -3.0},
            {    0.0,         23.0,         13.0,          0.0},
            {   -3.0,          0.0,          0.0,          2.0},
            {    3.0,         -4.0,         -2.0,         -2.0},
            {  -48.0,       -110.0,        -59.0,         26.0},
            {   51.0,        114.0,         61.0,        -27.0},
            { -133.0,          0.0,          0.0,         57.0},
            {    0.0,          4.0,          0.0,          0.0},
            {  -21.0,         -6.0,         -3.0,         11.0},
            {    0.0,         -3.0,         -1.0,          0.0},
            {  -11.0,        -21.0,        -11.0,          6.0},
            {  -18.0,       -436.0,       -233.0,          9.0},
            {   35.0,         -7.0,          0.0,          0.0},
            {    0.0,          5.0,          3.0,          0.0},
            {   11.0,         -3.0,         -1.0,         -6.0},
            {   -5.0,         -3.0,         -1.0,          3.0},
            {  -53.0,         -9.0,         -5.0,         28.0},
            {    0.0,          3.0,          2.0,          1.0},
            {    4.0,          0.0,          0.0,         -2.0},
            {    0.0,         -4.0,          0.0,          0.0},
            {  -50.0,        194.0,        103.0,         27.0},
            {  -13.0,         52.0,         28.0,          7.0},
            {  -91.0,        248.0,          0.0,          0.0},
            {    6.0,         49.0,         26.0,         -3.0},
            {   -6.0,        -47.0,        -25.0,          3.0},
            {    0.0,          5.0,          3.0,          0.0},
            {   52.0,         23.0,         10.0,        -23.0},
            {   -3.0,          0.0,          0.0,          1.0},
            {    0.0,          5.0,          3.0,          0.0},
            {   -4.0,          0.0,          0.0,          0.0},
            {   -4.0,          8.0,          3.0,          2.0},
            {   10.0,          0.0,          0.0,          0.0},
            {    3.0,          0.0,          0.0,         -2.0},
            {    0.0,          8.0,          4.0,          0.0},
            {    0.0,          8.0,          4.0,          1.0},
            {   -4.0,          0.0,          0.0,          0.0},
            {   -4.0,          0.0,          0.0,          0.0},
            {   -8.0,          4.0,          2.0,          4.0},
            {    8.0,         -4.0,         -2.0,         -4.0},
            {    0.0,         15.0,          7.0,          0.0},
            { -138.0,          0.0,          0.0,          0.0},
            {    0.0,         -7.0,         -3.0,          0.0},
            {    0.0,         -7.0,         -3.0,          0.0},
            {   54.0,          0.0,          0.0,        -29.0},
            {    0.0,         10.0,          4.0,          0.0},
            {   -7.0,          0.0,          0.0,          3.0},
            {  -37.0,         35.0,         19.0,         20.0},
            {    0.0,          4.0,          0.0,          0.0},
            {   -4.0,          9.0,          0.0,          0.0},
            {    8.0,          0.0,          0.0,         -4.0},
            {   -9.0,        -14.0,         -8.0,          5.0},
            {   -3.0,         -9.0,         -5.0,          3.0},
            { -145.0,         47.0,          0.0,          0.0},
            {  -10.0,         40.0,         21.0,          5.0},
            {   11.0,        -49.0,        -26.0,         -7.0},
            {-2150.0,          0.0,          0.0,        932.0},
            {  -12.0,          0.0,          0.0,          5.0},
            {   85.0,          0.0,          0.0,        -37.0},
            {    4.0,          0.0,          0.0,         -2.0},
            {    3.0,          0.0,          0.0,         -2.0},
            {  -86.0,        153.0,          0.0,          0.0},
            {   -6.0,          9.0,          5.0,          3.0},
            {    9.0,        -13.0,         -7.0,         -5.0},
            {   -8.0,         12.0,          6.0,          4.0},
            {  -51.0,          0.0,          0.0,         22.0},
            {  -11.0,       -268.0,       -116.0,          5.0},
            {    0.0,         12.0,          5.0,          0.0},
            {    0.0,          7.0,          3.0,          0.0},
            {   31.0,          6.0,          3.0,        -17.0},
            {  140.0,         27.0,         14.0,        -75.0},
            {   57.0,         11.0,          6.0,        -30.0},
            {  -14.0,        -39.0,          0.0,          0.0},
            {    0.0,         -6.0,         -2.0,          0.0},
            {    4.0,         15.0,          8.0,         -2.0},
            {    0.0,          4.0,          0.0,          0.0},
            {   -3.0,          0.0,          0.0,          1.0},
            {    0.0,         11.0,          5.0,          0.0},
            {    9.0,          6.0,          0.0,          0.0},
            {   -4.0,         10.0,          4.0,          2.0},
            {    5.0,          3.0,          0.0,          0.0},
            {   16.0,          0.0,          0.0,         -9.0},
            {   -3.0,          0.0,          0.0,          0.0},
            {    0.0,          3.0,          2.0,         -1.0},
            {    7.0,          0.0,          0.0,         -3.0},
            {  -25.0,         22.0,          0.0,          0.0},
            {   42.0,        223.0,        119.0,        -22.0},
            {  -27.0,       -143.0,        -77.0,         14.0},
            {    9.0,         49.0,         26.0,         -5.0},
            {-1166.0,          0.0,          0.0,        505.0},
            {   -5.0,          0.0,          0.0,          2.0},
            {   -6.0,          0.0,          0.0,          3.0},
            {   -8.0,          0.0,          1.0,          4.0},
            {    0.0,         -4.0,          0.0,          0.0},
            {  117.0,          0.0,          0.0,        -63.0},
            {   -4.0,          8.0,          4.0,          2.0},
            {    3.0,          0.0,          0.0,         -2.0},
            {   -5.0,          0.0,          0.0,          2.0},
            {    0.0,         31.0,          0.0,          0.0},
            {   -5.0,          0.0,          1.0,          3.0},
            {    4.0,          0.0,          0.0,         -2.0},
            {   -4.0,          0.0,          0.0,          2.0},
            {  -24.0,        -13.0,         -6.0,         10.0},
            {    3.0,          0.0,          0.0,          0.0},
            {    0.0,        -32.0,        -17.0,          0.0},
            {    8.0,         12.0,          5.0,         -3.0},
            {    3.0,          0.0,          0.0,         -1.0},
            {    7.0,         13.0,          0.0,          0.0},
            {   -3.0,         16.0,          0.0,          0.0},
            {   50.0,          0.0,          0.0,        -27.0},
            {    0.0,         -5.0,         -3.0,          0.0},
            {   13.0,          0.0,          0.0,          0.0},
            {    0.0,          5.0,          3.0,          1.0},
            {   24.0,          5.0,          2.0,        -11.0},
            {    5.0,        -11.0,         -5.0,         -2.0},
            {   30.0,         -3.0,         -2.0,        -16.0},
            {   18.0,          0.0,          0.0,         -9.0},
            {    8.0,        614.0,          0.0,          0.0},
            {    3.0,         -3.0,         -1.0,         -2.0},
            {    6.0,         17.0,          9.0,         -3.0},
            {   -3.0,         -9.0,         -5.0,          2.0},
            {    0.0,          6.0,          3.0,         -1.0},
            { -127.0,         21.0,          9.0,         55.0},
            {    3.0,          5.0,          0.0,          0.0},
            {   -6.0,        -10.0,         -4.0,          3.0},
            {    5.0,          0.0,          0.0,          0.0},
            {   16.0,          9.0,          4.0,         -7.0},
            {    3.0,          0.0,          0.0,         -2.0},
            {    0.0,         22.0,          0.0,          0.0},
            {    0.0,         19.0,         10.0,          0.0},
            {    7.0,          0.0,          0.0,         -4.0},
            {    0.0,         -5.0,         -2.0,          0.0},
            {    0.0,          3.0,          1.0,          0.0},
            {   -9.0,          3.0,          1.0,          4.0},
            {   17.0,          0.0,          0.0,         -7.0},
            {    0.0,         -3.0,         -2.0,         -1.0},
            {  -20.0,         34.0,          0.0,          0.0},
            {  -10.0,          0.0,          1.0,          5.0},
            {   -4.0,          0.0,          0.0,          2.0},
            {   22.0,        -87.0,          0.0,          0.0},
            {   -4.0,          0.0,          0.0,          2.0},
            {   -3.0,         -6.0,         -2.0,          1.0},
            {  -16.0,         -3.0,         -1.0,          7.0},
            {    0.0,         -3.0,         -2.0,          0.0},
            {    4.0,          0.0,          0.0,          0.0},
            {  -68.0,         39.0,          0.0,          0.0},
            {   27.0,          0.0,          0.0,        -14.0},
            {    0.0,         -4.0,          0.0,          0.0},
            {  -25.0,          0.0,          0.0,          0.0},
            {  -12.0,         -3.0,         -2.0,          6.0},
            {    3.0,          0.0,          0.0,         -1.0},
            {    3.0,         66.0,         29.0,         -1.0},
            {  490.0,          0.0,          0.0,       -213.0},
            {  -22.0,         93.0,         49.0,         12.0},
            {   -7.0,         28.0,         15.0,          4.0},
            {   -3.0,         13.0,          7.0,          2.0},
            {  -46.0,         14.0,          0.0,          0.0},
            {   -5.0,          0.0,          0.0,          0.0},
            {    2.0,          1.0,          0.0,          0.0},
            {    0.0,         -3.0,          0.0,          0.0},
            {  -28.0,          0.0,          0.0,         15.0},
            {    5.0,          0.0,          0.0,         -2.0},
            {    0.0,          3.0,          0.0,          0.0},
            {  -11.0,          0.0,          0.0,          5.0},
            {    0.0,          3.0,          1.0,          0.0},
            {   -3.0,          0.0,          0.0,          1.0},
            {   25.0,        106.0,         57.0,        -13.0},
            {    5.0,         21.0,         11.0,         -3.0},
            { 1485.0,          0.0,          0.0,          0.0},
            {   -7.0,        -32.0,        -17.0,          4.0},
            {    0.0,          5.0,          3.0,          0.0},
            {   -6.0,         -3.0,         -2.0,          3.0},
            {   30.0,         -6.0,         -2.0,        -13.0},
            {   -4.0,          4.0,          0.0,          0.0},
            {  -19.0,          0.0,          0.0,         10.0},
            {    0.0,          4.0,          2.0,         -1.0},
            {    0.0,          3.0,          0.0,          0.0},
            {    4.0,          0.0,          0.0,         -2.0},
            {    0.0,         -3.0,         -1.0,          0.0},
            {   -3.0,          0.0,          0.0,          0.0},
            {    5.0,          3.0,          1.0,         -2.0},
            {    0.0,         11.0,          0.0,          0.0},
            {  118.0,          0.0,          0.0,        -52.0},
            {    0.0,         -5.0,         -3.0,          0.0},
            {  -28.0,         36.0,          0.0,          0.0},
            {    5.0,         -5.0,          0.0,          0.0},
            {   14.0,        -59.0,        -31.0,         -8.0},
            {    0.0,          9.0,          5.0,          1.0},
            { -458.0,          0.0,          0.0,        198.0},
            {    0.0,        -45.0,        -20.0,          0.0},
            {    9.0,          0.0,          0.0,         -5.0},
            {    0.0,         -3.0,          0.0,          0.0},
            {    0.0,         -4.0,         -2.0,         -1.0},
            {   11.0,          0.0,          0.0,         -6.0},
            {    6.0,          0.0,          0.0,         -2.0},
            {  -16.0,         23.0,          0.0,          0.0},
            {    0.0,         -4.0,         -2.0,          0.0},
            {   -5.0,          0.0,          0.0,          2.0},
            { -166.0,        269.0,          0.0,          0.0},
            {   15.0,          0.0,          0.0,         -8.0},
            {   10.0,          0.0,          0.0,         -4.0},
            {  -78.0,         45.0,          0.0,          0.0},
            {    0.0,         -5.0,         -2.0,          0.0},
            {    7.0,          0.0,          0.0,         -4.0},
            {   -5.0,        328.0,          0.0,          0.0},
            {    3.0,          0.0,          0.0,         -2.0},
            {    5.0,          0.0,          0.0,         -2.0},
            {    0.0,          3.0,          1.0,          0.0},
            {   -3.0,          0.0,          0.0,          0.0},
            {   -3.0,          0.0,          0.0,          0.0},
            {    0.0,         -4.0,         -2.0,          0.0},
            {-1223.0,        -26.0,          0.0,          0.0},
            {    0.0,          7.0,          3.0,          0.0},
            {    3.0,          0.0,          0.0,          0.0},
            {    0.0,          3.0,          2.0,          0.0},
            {   -6.0,         20.0,          0.0,          0.0},
            { -368.0,          0.0,          0.0,          0.0},
            {  -75.0,          0.0,          0.0,          0.0},
            {   11.0,          0.0,          0.0,         -6.0},
            {    3.0,          0.0,          0.0,         -2.0},
            {   -3.0,          0.0,          0.0,          1.0},
            {  -13.0,        -30.0,          0.0,          0.0},
            {   21.0,          3.0,          0.0,          0.0},
            {   -3.0,          0.0,          0.0,          1.0},
            {   -4.0,          0.0,          0.0,          2.0},
            {    8.0,        -27.0,          0.0,          0.0},
            {  -19.0,        -11.0,          0.0,          0.0},
            {   -4.0,          0.0,          0.0,          2.0},
            {    0.0,          5.0,          2.0,          0.0},
            {   -6.0,          0.0,          0.0,          2.0},
            {   -8.0,          0.0,          0.0,          0.0},
            {   -1.0,          0.0,          0.0,          0.0},
            {  -14.0,          0.0,          0.0,          6.0},
            {    6.0,          0.0,          0.0,          0.0},
            {  -74.0,          0.0,          0.0,         32.0},
            {    0.0,         -3.0,         -1.0,          0.0},
            {    4.0,          0.0,          0.0,         -2.0},
            {    8.0,         11.0,          0.0,          0.0},
            {    0.0,          3.0,          2.0,          0.0},
            { -262.0,          0.0,          0.0,        114.0},
            {    0.0,         -4.0,          0.0,          0.0},
            {   -7.0,          0.0,          0.0,          4.0},
            {    0.0,        -27.0,        -12.0,          0.0},
            {  -19.0,         -8.0,         -4.0,          8.0},
            {  202.0,          0.0,          0.0,        -87.0},
            {   -8.0,         35.0,         19.0,          5.0},
            {    0.0,          4.0,          2.0,          0.0},
            {   16.0,         -5.0,          0.0,          0.0},
            {    5.0,          0.0,          0.0,         -3.0},
            {    0.0,         -3.0,          0.0,          0.0},
            {    1.0,          0.0,          0.0,          0.0},
            {  -35.0,        -48.0,        -21.0,         15.0},
            {   -3.0,         -5.0,         -2.0,          1.0},
            {    6.0,          0.0,          0.0,         -3.0},
            {    3.0,          0.0,          0.0,         -1.0},
            {    0.0,         -5.0,          0.0,          0.0},
            {   12.0,         55.0,         29.0,         -6.0},
            {    0.0,          5.0,          3.0,          0.0},
            { -598.0,          0.0,          0.0,          0.0},
            {   -3.0,        -13.0,         -7.0,          1.0},
            {   -5.0,         -7.0,         -3.0,          2.0},
            {    3.0,          0.0,          0.0,         -1.0},
            {    5.0,         -7.0,          0.0,          0.0},
            {    4.0,          0.0,          0.0,         -2.0},
            {   16.0,         -6.0,          0.0,          0.0},
            {    8.0,         -3.0,          0.0,          0.0},
            {    8.0,        -31.0,        -16.0,         -4.0},
            {    0.0,          3.0,          1.0,          0.0},
            {  113.0,          0.0,          0.0,        -49.0},
            {    0.0,        -24.0,        -10.0,          0.0},
            {    4.0,          0.0,          0.0,         -2.0},
            {   27.0,          0.0,          0.0,          0.0},
            {   -3.0,          0.0,          0.0,          1.0},
            {    0.0,         -4.0,         -2.0,          0.0},
            {    5.0,          0.0,          0.0,         -2.0},
            {    0.0,         -3.0,          0.0,          0.0},
            {  -13.0,          0.0,          0.0,          6.0},
            {    5.0,          0.0,          0.0,         -2.0},
            {  -18.0,        -10.0,         -4.0,          8.0},
            {   -4.0,        -28.0,          0.0,          0.0},
            {   -5.0,          6.0,          3.0,          2.0},
            {   -3.0,          0.0,          0.0,          1.0},
            {   -5.0,         -9.0,         -4.0,          2.0},
            {   17.0,          0.0,          0.0,         -7.0},
            {   11.0,          4.0,          0.0,          0.0},
            {    0.0,         -6.0,         -2.0,          0.0},
            {   83.0,         15.0,          0.0,          0.0},
            {   -4.0,          0.0,          0.0,          2.0},
            {    0.0,       -114.0,        -49.0,          0.0},
            {  117.0,          0.0,          0.0,        -51.0},
            {   -5.0,         19.0,         10.0,          2.0},
            {   -3.0,          0.0,          0.0,          0.0},
            {   -3.0,          0.0,          0.0,          2.0},
            {    0.0,         -3.0,         -1.0,          0.0},
            {    3.0,          0.0,          0.0,          0.0},
            {    0.0,         -6.0,         -2.0,          0.0},
            {  393.0,          3.0,          0.0,          0.0},
            {   -4.0,         21.0,         11.0,          2.0},
            {   -6.0,          0.0,         -1.0,          3.0},
            {   -3.0,          8.0,          4.0,          1.0},
            {    8.0,          0.0,          0.0,          0.0},
            {   18.0,        -29.0,        -13.0,         -8.0},
            {    8.0,         34.0,         18.0,         -4.0},
            {   89.0,          0.0,          0.0,          0.0},
            {    3.0,         12.0,          6.0,         -1.0},
            {   54.0,        -15.0,         -7.0,        -24.0},
            {    0.0,          3.0,          0.0,          0.0},
            {    3.0,          0.0,          0.0,         -1.0},
            {    0.0,         35.0,          0.0,          0.0},
            { -154.0,        -30.0,        -13.0,         67.0},
            {   15.0,          0.0,          0.0,          0.0},
            {    0.0,          4.0,          2.0,          0.0},
            {    0.0,          9.0,          0.0,          0.0},
            {   80.0,        -71.0,        -31.0,        -35.0},
            {    0.0,        -20.0,         -9.0,          0.0},
            {   11.0,          5.0,          2.0,         -5.0},
            {   61.0,        -96.0,        -42.0,        -27.0},
            {   14.0,          9.0,          4.0,         -6.0},
            {  -11.0,         -6.0,         -3.0,          5.0},
            {    0.0,         -3.0,         -1.0,          0.0},
            {  123.0,       -415.0,       -180.0,        -53.0},
            {    0.0,          0.0,          0.0,        -35.0},
            {   -5.0,          0.0,          0.0,          0.0},
            {    7.0,        -32.0,        -17.0,         -4.0},
            {    0.0,         -9.0,         -5.0,          0.0},
            {    0.0,         -4.0,          2.0,          0.0},
            {  -89.0,          0.0,          0.0,         38.0},
            {    0.0,        -86.0,        -19.0,         -6.0},
            {    0.0,          0.0,        -19.0,          6.0},
            { -123.0,       -416.0,       -180.0,         53.0},
            {    0.0,         -3.0,         -1.0,          0.0},
            {   12.0,         -6.0,         -3.0,         -5.0},
            {  -13.0,          9.0,          4.0,          6.0},
            {    0.0,        -15.0,         -7.0,          0.0},
            {    3.0,          0.0,          0.0,         -1.0},
            {  -62.0,        -97.0,        -42.0,         27.0},
            {  -11.0,          5.0,          2.0,          5.0},
            {    0.0,        -19.0,         -8.0,          0.0},
            {   -3.0,          0.0,          0.0,          1.0},
            {    0.0,          4.0,          2.0,          0.0},
            {    0.0,          3.0,          0.0,          0.0},
            {    0.0,          4.0,          2.0,          0.0},
            {  -85.0,        -70.0,        -31.0,         37.0},
            {  163.0,        -12.0,         -5.0,        -72.0},
            {  -63.0,        -16.0,         -7.0,         28.0},
            {  -21.0,        -32.0,        -14.0,          9.0},
            {    0.0,         -3.0,         -1.0,          0.0},
            {    3.0,          0.0,          0.0,         -2.0},
            {    0.0,          8.0,          0.0,          0.0},
            {    3.0,         10.0,          4.0,         -1.0},
            {    3.0,          0.0,          0.0,         -1.0},
            {    0.0,         -7.0,         -3.0,          0.0},
            {    0.0,         -4.0,         -2.0,          0.0},
            {    6.0,         19.0,          0.0,          0.0},
            {    5.0,       -173.0,        -75.0,         -2.0},
            {    0.0,         -7.0,         -3.0,          0.0},
            {    7.0,        -12.0,         -5.0,         -3.0},
            {   -3.0,          0.0,          0.0,          2.0},
            {    3.0,         -4.0,         -2.0,         -1.0},
            {   74.0,          0.0,          0.0,        -32.0},
            {   -3.0,         12.0,          6.0,          2.0},
            {   26.0,        -14.0,         -6.0,        -11.0},
            {   19.0,          0.0,          0.0,         -8.0},
            {    6.0,         24.0,         13.0,         -3.0},
            {   83.0,          0.0,          0.0,          0.0},
            {    0.0,        -10.0,         -5.0,          0.0},
            {   11.0,         -3.0,         -1.0,         -5.0},
            {    3.0,          0.0,          1.0,         -1.0},
            {    3.0,          0.0,          0.0,         -1.0},
            {   -4.0,          0.0,          0.0,          0.0},
            {    5.0,        -23.0,        -12.0,         -3.0},
            { -339.0,          0.0,          0.0,        147.0},
            {    0.0,        -10.0,         -5.0,          0.0},
            {    5.0,          0.0,          0.0,          0.0},
            {    3.0,          0.0,          0.0,         -1.0},
            {    0.0,         -4.0,         -2.0,          0.0},
            {   18.0,         -3.0,          0.0,          0.0},
            {    9.0,        -11.0,         -5.0,         -4.0},
            {   -8.0,          0.0,          0.0,          4.0},
            {    3.0,          0.0,          0.0,         -1.0},
            {    0.0,          9.0,          0.0,          0.0},
            {    6.0,         -9.0,         -4.0,         -2.0},
            {   -4.0,        -12.0,          0.0,          0.0},
            {   67.0,        -91.0,        -39.0,        -29.0},
            {   30.0,        -18.0,         -8.0,        -13.0},
            {    0.0,          0.0,          0.0,          0.0},
            {    0.0,       -114.0,        -50.0,          0.0},
            {    0.0,          0.0,          0.0,         23.0},
            {  517.0,         16.0,          7.0,       -224.0},
            {    0.0,         -7.0,         -3.0,          0.0},
            {  143.0,         -3.0,         -1.0,        -62.0},
            {   29.0,          0.0,          0.0,        -13.0},
            {   -4.0,          0.0,          0.0,          2.0},
            {   -6.0,          0.0,          0.0,          3.0},
            {    5.0,         12.0,          5.0,         -2.0},
            {  -25.0,          0.0,          0.0,         11.0},
            {   -3.0,          0.0,          0.0,          1.0},
            {    0.0,          4.0,          2.0,          0.0},
            {  -22.0,         12.0,          5.0,         10.0},
            {   50.0,          0.0,          0.0,        -22.0},
            {    0.0,          7.0,          4.0,          0.0},
            {    0.0,          3.0,          1.0,          0.0},
            {   -4.0,          4.0,          2.0,          2.0},
            {   -5.0,        -11.0,         -5.0,          2.0},
            {    0.0,          4.0,          2.0,          0.0},
            {    4.0,         17.0,          9.0,         -2.0},
            {   59.0,          0.0,          0.0,          0.0},
            {    0.0,         -4.0,         -2.0,          0.0},
            {   -8.0,          0.0,          0.0,          4.0},
            {   -3.0,          0.0,          0.0,          0.0},
            {    4.0,        -15.0,         -8.0,         -2.0},
            {  370.0,         -8.0,          0.0,       -160.0},
            {    0.0,          0.0,         -3.0,          0.0},
            {    0.0,          3.0,          1.0,          0.0},
            {   -6.0,          3.0,          1.0,          3.0},
            {    0.0,          6.0,          0.0,          0.0},
            {  -10.0,          0.0,          0.0,          4.0},
            {    0.0,          9.0,          4.0,          0.0},
            {    4.0,         17.0,          7.0,         -2.0},
            {   34.0,          0.0,          0.0,        -15.0},
            {    0.0,          5.0,          3.0,          0.0},
            {   -5.0,          0.0,          0.0,          2.0},
            {  -37.0,         -7.0,         -3.0,         16.0},
            {    3.0,         13.0,          7.0,         -2.0},
            {   40.0,          0.0,          0.0,          0.0},
            {    0.0,         -3.0,         -2.0,          0.0},
            { -184.0,         -3.0,         -1.0,         80.0},
            {   -3.0,          0.0,          0.0,          1.0},
            {   -3.0,          0.0,          0.0,          0.0},
            {    0.0,        -10.0,         -6.0,         -1.0},
            {   31.0,         -6.0,          0.0,        -13.0},
            {   -3.0,        -32.0,        -14.0,          1.0},
            {   -7.0,          0.0,          0.0,          3.0},
            {    0.0,         -8.0,         -4.0,          0.0},
            {    3.0,         -4.0,          0.0,          0.0},
            {    0.0,          4.0,          0.0,          0.0},
            {    0.0,          3.0,          1.0,          0.0},
            {   19.0,        -23.0,        -10.0,          2.0},
            {    0.0,          0.0,          0.0,        -10.0},
            {    0.0,          3.0,          2.0,          0.0},
            {    0.0,          9.0,          5.0,         -1.0},
            {   28.0,          0.0,          0.0,          0.0},
            {    0.0,         -7.0,         -4.0,          0.0},
            {    8.0,         -4.0,          0.0,         -4.0},
            {    0.0,          0.0,         -2.0,          0.0},
            {    0.0,          3.0,          0.0,          0.0},
            {   -3.0,          0.0,          0.0,          1.0},
            {   -9.0,          0.0,          1.0,          4.0},
            {    3.0,         12.0,          5.0,         -1.0},
            {   17.0,         -3.0,         -1.0,          0.0},
            {    0.0,          7.0,          4.0,          0.0},
            {   19.0,          0.0,          0.0,          0.0},
            {    0.0,         -5.0,         -3.0,          0.0},
            {   14.0,         -3.0,          0.0,         -1.0},
            {    0.0,          0.0,         -1.0,          0.0},
            {    0.0,          0.0,          0.0,         -5.0},
            {    0.0,          5.0,          3.0,          0.0},
            {   13.0,          0.0,          0.0,          0.0},
            {    0.0,         -3.0,         -2.0,          0.0},
            {    2.0,          9.0,          4.0,          3.0},
            {    0.0,          0.0,          0.0,         -4.0},
            {    8.0,          0.0,          0.0,          0.0},
            {    0.0,          4.0,          2.0,          0.0},
            {    6.0,          0.0,          0.0,         -3.0},
            {    6.0,          0.0,          0.0,          0.0},
            {    0.0,          3.0,          1.0,          0.0},
            {    5.0,          0.0,          0.0,         -2.0},
            {    3.0,          0.0,          0.0,         -1.0},
            {   -3.0,          0.0,          0.0,          0.0},
            {    6.0,          0.0,          0.0,          0.0},
            {    7.0,          0.0,          0.0,          0.0},
            {   -4.0,          0.0,          0.0,          0.0},
            {    4.0,          0.0,          0.0,          0.0},
            {    6.0,          0.0,          0.0,          0.0},
            {    0.0,         -4.0,          0.0,          0.0},
            {    0.0,         -4.0,          0.0,          0.0},
            {    5.0,          0.0,          0.0,          0.0},
            {   -3.0,          0.0,          0.0,          0.0},
            {    4.0,          0.0,          0.0,          0.0},
            {   -5.0,          0.0,          0.0,          0.0},
            {    4.0,          0.0,          0.0,          0.0},
            {    0.0,          3.0,          0.0,          0.0},
            {   13.0,          0.0,          0.0,          0.0},
            {   21.0,         11.0,          0.0,          0.0},
            {    0.0,         -5.0,          0.0,          0.0},
            {    0.0,         -5.0,         -2.0,          0.0},
            {    0.0,          5.0,          3.0,          0.0},
            {    0.0,         -5.0,          0.0,          0.0},
            {   -3.0,          0.0,          0.0,          2.0},
            {   20.0,         10.0,          0.0,          0.0},
            {  -34.0,          0.0,          0.0,          0.0},
            {  -19.0,          0.0,          0.0,          0.0},
            {    3.0,          0.0,          0.0,         -2.0},
            {   -3.0,          0.0,          0.0,          1.0},
            {   -6.0,          0.0,          0.0,          3.0},
            {   -4.0,          0.0,          0.0,          0.0},
            {    3.0,          0.0,          0.0,          0.0},
            {    3.0,          0.0,          0.0,          0.0},
            {    4.0,          0.0,          0.0,          0.0},
            {    3.0,          0.0,          0.0,         -1.0},
            {    6.0,          0.0,          0.0,         -3.0},
            {   -8.0,          0.0,          0.0,          3.0},
            {    0.0,          3.0,          1.0,          0.0},
            {   -3.0,          0.0,          0.0,          0.0},
            {    0.0,         -3.0,         -2.0,          0.0},
            {  126.0,        -63.0,        -27.0,        -55.0},
            {   -5.0,          0.0,          1.0,          2.0},
            {   -3.0,         28.0,         15.0,          2.0},
            {    5.0,          0.0,          1.0,         -2.0},
            {    0.0,          9.0,          4.0,          1.0},
            {    0.0,          9.0,          4.0,         -1.0},
            { -126.0,        -63.0,        -27.0,         55.0},
            {    3.0,          0.0,          0.0,         -1.0},
            {   21.0,        -11.0,         -6.0,        -11.0},
            {    0.0,         -4.0,          0.0,          0.0},
            {  -21.0,        -11.0,         -6.0,         11.0},
            {   -3.0,          0.0,          0.0,          1.0},
            {    0.0,          3.0,          1.0,          0.0},
            {    8.0,          0.0,          0.0,         -4.0},
            {   -6.0,          0.0,          0.0,          3.0},
            {   -3.0,          0.0,          0.0,          1.0},
            {    3.0,          0.0,          0.0,         -1.0},
            {   -3.0,          0.0,          0.0,          1.0},
            {   -5.0,          0.0,          0.0,          2.0},
            {   24.0,        -12.0,         -5.0,        -11.0},
            {    0.0,          3.0,          1.0,          0.0},
            {    0.0,          3.0,          1.0,          0.0},
            {    0.0,          3.0,          2.0,          0.0},
            {  -24.0,        -12.0,         -5.0,         10.0},
            {    4.0,          0.0,         -1.0,         -2.0},
            {   13.0,          0.0,          0.0,         -6.0},
            {    7.0,          0.0,          0.0,         -3.0},
            {    3.0,          0.0,          0.0,         -1.0},
            {    3.0,          0.0,          0.0,         -1.0}};
        
        /*
         Interval between fundamental epoch J2000.0 and given date.
         */
        
        T = JC_TT;
        
        /*
         Compute fundamental arguments from Simon et al. (1994), in radians.
         */
        
        //fund_args (T, a);
        // l (mean anomaly of the Moon)
        a[0] = fmod (485868.249036+T*(1717915923.2178+T*(31.8792+T*(0.051635+T*(-0.00024470)))), 1296000.0)*AS2R;
        //l' (mean anomaly of the Sun)
        a[1] = fmod (1287104.79305+T*(129596581.0481+T*(- 0.5532 +T*(0.000136+T*(- 0.00001149)))), 1296000.0)*AS2R;
        //F (mean argument of the latitude of the Moon)
        a[2] = fmod (335779.526232+T*(1739527262.8478+T*(- 12.7512+T*(-0.001037+T* (0.00000417)))), 1296000.0)*AS2R;
        //D (mean elongation of the Moon from the Sun)
        a[3] = fmod (1072260.70369+T*(1602961601.2090+T*(- 6.3706 +T* (0.006593 +T*(- 0.00003169)))), 1296000.0)*AS2R;
        //Omega (mean longitude of the Moon's ascending node); from Simon section 3.4(b.3), precession = 5028.8200 arcsec/cy)
        a[4] = fmod (450160.398036+T*(- 6962890.5431+T*(7.4722 + T*(0.007702 + T*(- 0.00005939)))), 1296000.0)*AS2R;
        
        /*
         ** Luni-solar nutation. **
         */
        
        /*
         Initialize the nutation values.
         */
        
        dp = 0.0;
        de = 0.0;
        
        /*
         Summation of luni-solar nutation series (in reverse order).
         */
        
        for (i = 677; i >= 0; i--)
        {
            
            /*
             Argument and functions.
             */
            
            arg = fmod ((double) nals_t[i][0] * a[0]  +
                        (double) nals_t[i][1] * a[1]  +
                        (double) nals_t[i][2] * a[2]  +
                        (double) nals_t[i][3] * a[3]  +
                        (double) nals_t[i][4] * a[4], TWOPI);
            
            sarg = sin (arg);
            carg = cos (arg);
            
            /*
             Term.
             */
            
            dp += (cls_t[i][0] + cls_t[i][1] * T) * sarg
            +   cls_t[i][2] * carg;
            de += (cls_t[i][3] + cls_t[i][4] * T) * carg
            +   cls_t[i][5] * sarg;
        }
        
        /*
         Convert from 0.1 microarcsec units to radians.
         */
        
        factor = 1.0E-7 * AS2R; //ASEC2RAD;
        dpsils = dp * factor;
        depsls = de * factor;
        
        /*
         ** Planetary nutation. **
         */
        
        /*
         Mean anomaly of the Moon.
         */
        
        al = fmod (2.35555598 + 8328.6914269554 * T, TWOPI);
        
        /*
         Mean anomaly of the Sun.
         */
        
        alsu = fmod (6.24006013 + 628.301955 * T, TWOPI);
        
        /*
         Mean argument of the latitude of the Moon.
         */
        
        af = fmod (1.627905234 + 8433.466158131 * T, TWOPI);
        
        /*
         Mean elongation of the Moon from the Sun.
         */
        
        ad = fmod (5.198466741 + 7771.3771468121 * T, TWOPI);
        
        /*
         Mean longitude of the ascending node of the Moon.
         */
        
        aom = fmod (2.18243920 - 33.757045 * T, TWOPI);
        
        /*
         General accumulated precession in longitude.
         */
        
        apa = (0.02438175 + 0.00000538691 * T) * T;
        /*
         Planetary longitudes, Mercury through Neptune (Souchay et al. 1999).
         */
        
        alme = fmod (4.402608842 + 2608.7903141574 * T, TWOPI);
        alve = fmod (3.176146697 + 1021.3285546211 * T, TWOPI);
        alea = fmod (1.753470314 +  628.3075849991 * T, TWOPI);
        alma = fmod (6.203480913 +  334.0612426700 * T, TWOPI);
        alju = fmod (0.599546497 +   52.9690962641 * T, TWOPI);
        alsa = fmod (0.874016757 +   21.3299104960 * T, TWOPI);
        alur = fmod (5.481293871 +    7.4781598567 * T, TWOPI);
        alne = fmod (5.321159000 +    3.8127774000 * T, TWOPI);
        
        /*
         Initialize the nutation values.
         */
        
        dp = 0.0;
        de = 0.0;
        
        /*
         Summation of planetary nutation series (in reverse order).
         */
        
        for (i = 686; i >= 0; i--)
        {
            
            /*
             Argument and functions.
             */
            
            arg = fmod ((double) napl_t[i][ 0] * al    +
                        (double) napl_t[i][ 1] * alsu  +
                        (double) napl_t[i][ 2] * af    +
                        (double) napl_t[i][ 3] * ad    +
                        (double) napl_t[i][ 4] * aom   +
                        (double) napl_t[i][ 5] * alme  +
                        (double) napl_t[i][ 6] * alve  +
                        (double) napl_t[i][ 7] * alea  +
                        (double) napl_t[i][ 8] * alma  +
                        (double) napl_t[i][ 9] * alju  +
                        (double) napl_t[i][10] * alsa  +
                        (double) napl_t[i][11] * alur  +
                        (double) napl_t[i][12] * alne  +
                        (double) napl_t[i][13] * apa, TWOPI);
            
            sarg = sin (arg);
            carg = cos (arg);
            
            /*
             Term.
             */
            
            dp += cpl_t[i][0] * sarg + cpl_t[i][1] * carg;
            de += cpl_t[i][2] * sarg + cpl_t[i][3] * carg;
        }
        
        dpsipl = dp * factor;
        depspl = de * factor;
        
        /*
         Total: Add planetary and luni-solar components.
         */
        
        dpsi = dpsipl + dpsils;
        deps = depspl + depsls;
        
    }
    
    
    /* iau 1980 nutation model
     *
     */
    void GEarthOrientationParameter::iauNutation1980(double JC_TT, double& dpsi,double& deps)
    {
        double PI = GCONST("PI");
        // Arcseconds to radians
        const double DAS2R = GCONST("AS2R");    // 4.848136811095359935899141e-6;
        static const double nut[106][10] = {
            {   0,   0,   0,   0,   1, -6798.4, -171996, -174.2, 92025,   8.9},
            {   0,   0,   2,  -2,   2,   182.6,  -13187,   -1.6,  5736,  -3.1},
            {   0,   0,   2,   0,   2,    13.7,   -2274,   -0.2,   977,  -0.5},
            {   0,   0,   0,   0,   2, -3399.2,    2062,    0.2,  -895,   0.5},
            {   0,  -1,   0,   0,   0,  -365.3,   -1426,    3.4,    54,  -0.1},
            {   1,   0,   0,   0,   0,    27.6,     712,    0.1,    -7,   0.0},
            {   0,   1,   2,  -2,   2,   121.7,    -517,    1.2,   224,  -0.6},
            {   0,   0,   2,   0,   1,    13.6,    -386,   -0.4,   200,   0.0},
            {   1,   0,   2,   0,   2,     9.1,    -301,    0.0,   129,  -0.1},
            {   0,  -1,   2,  -2,   2,   365.2,     217,   -0.5,   -95,   0.3},
            {  -1,   0,   0,   2,   0,    31.8,     158,    0.0,    -1,   0.0},
            {   0,   0,   2,  -2,   1,   177.8,     129,    0.1,   -70,   0.0},
            {  -1,   0,   2,   0,   2,    27.1,     123,    0.0,   -53,   0.0},
            {   1,   0,   0,   0,   1,    27.7,      63,    0.1,   -33,   0.0},
            {   0,   0,   0,   2,   0,    14.8,      63,    0.0,    -2,   0.0},
            {  -1,   0,   2,   2,   2,     9.6,     -59,    0.0,    26,   0.0},
            {  -1,   0,   0,   0,   1,   -27.4,     -58,   -0.1,    32,   0.0},
            {   1,   0,   2,   0,   1,     9.1,     -51,    0.0,    27,   0.0},
            {  -2,   0,   0,   2,   0,  -205.9,     -48,    0.0,     1,   0.0},
            {  -2,   0,   2,   0,   1,  1305.5,      46,    0.0,   -24,   0.0},
            {   0,   0,   2,   2,   2,     7.1,     -38,    0.0,    16,   0.0},
            {   2,   0,   2,   0,   2,     6.9,     -31,    0.0,    13,   0.0},
            {   2,   0,   0,   0,   0,    13.8,      29,    0.0,    -1,   0.0},
            {   1,   0,   2,  -2,   2,    23.9,      29,    0.0,   -12,   0.0},
            {   0,   0,   2,   0,   0,    13.6,      26,    0.0,    -1,   0.0},
            {   0,   0,   2,  -2,   0,   173.3,     -22,    0.0,     0,   0.0},
            {  -1,   0,   2,   0,   1,    27.0,      21,    0.0,   -10,   0.0},
            {   0,   2,   0,   0,   0,   182.6,      17,   -0.1,     0,   0.0},
            {   0,   2,   2,  -2,   2,    91.3,     -16,    0.1,     7,   0.0},
            {  -1,   0,   0,   2,   1,    32.0,      16,    0.0,    -8,   0.0},
            {   0,   1,   0,   0,   1,   386.0,     -15,    0.0,     9,   0.0},
            {   1,   0,   0,  -2,   1,   -31.7,     -13,    0.0,     7,   0.0},
            {   0,  -1,   0,   0,   1,  -346.6,     -12,    0.0,     6,   0.0},
            {   2,   0,  -2,   0,   0, -1095.2,      11,    0.0,     0,   0.0},
            {  -1,   0,   2,   2,   1,     9.5,     -10,    0.0,     5,   0.0},
            {   1,   0,   2,   2,   2,     5.6,      -8,    0.0,     3,   0.0},
            {   0,  -1,   2,   0,   2,    14.2,      -7,    0.0,     3,   0.0},
            {   0,   0,   2,   2,   1,     7.1,      -7,    0.0,     3,   0.0},
            {   1,   1,   0,  -2,   0,   -34.8,      -7,    0.0,     0,   0.0},
            {   0,   1,   2,   0,   2,    13.2,       7,    0.0,    -3,   0.0},
            {  -2,   0,   0,   2,   1,  -199.8,      -6,    0.0,     3,   0.0},
            {   0,   0,   0,   2,   1,    14.8,      -6,    0.0,     3,   0.0},
            {   2,   0,   2,  -2,   2,    12.8,       6,    0.0,    -3,   0.0},
            {   1,   0,   0,   2,   0,     9.6,       6,    0.0,     0,   0.0},
            {   1,   0,   2,  -2,   1,    23.9,       6,    0.0,    -3,   0.0},
            {   0,   0,   0,  -2,   1,   -14.7,      -5,    0.0,     3,   0.0},
            {   0,  -1,   2,  -2,   1,   346.6,      -5,    0.0,     3,   0.0},
            {   2,   0,   2,   0,   1,     6.9,      -5,    0.0,     3,   0.0},
            {   1,  -1,   0,   0,   0,    29.8,       5,    0.0,     0,   0.0},
            {   1,   0,   0,  -1,   0,   411.8,      -4,    0.0,     0,   0.0},
            {   0,   0,   0,   1,   0,    29.5,      -4,    0.0,     0,   0.0},
            {   0,   1,   0,  -2,   0,   -15.4,      -4,    0.0,     0,   0.0},
            {   1,   0,  -2,   0,   0,   -26.9,       4,    0.0,     0,   0.0},
            {   2,   0,   0,  -2,   1,   212.3,       4,    0.0,    -2,   0.0},
            {   0,   1,   2,  -2,   1,   119.6,       4,    0.0,    -2,   0.0},
            {   1,   1,   0,   0,   0,    25.6,      -3,    0.0,     0,   0.0},
            {   1,  -1,   0,  -1,   0, -3232.9,      -3,    0.0,     0,   0.0},
            {  -1,  -1,   2,   2,   2,     9.8,      -3,    0.0,     1,   0.0},
            {   0,  -1,   2,   2,   2,     7.2,      -3,    0.0,     1,   0.0},
            {   1,  -1,   2,   0,   2,     9.4,      -3,    0.0,     1,   0.0},
            {   3,   0,   2,   0,   2,     5.5,      -3,    0.0,     1,   0.0},
            {  -2,   0,   2,   0,   2,  1615.7,      -3,    0.0,     1,   0.0},
            {   1,   0,   2,   0,   0,     9.1,       3,    0.0,     0,   0.0},
            {  -1,   0,   2,   4,   2,     5.8,      -2,    0.0,     1,   0.0},
            {   1,   0,   0,   0,   2,    27.8,      -2,    0.0,     1,   0.0},
            {  -1,   0,   2,  -2,   1,   -32.6,      -2,    0.0,     1,   0.0},
            {   0,  -2,   2,  -2,   1,  6786.3,      -2,    0.0,     1,   0.0},
            {  -2,   0,   0,   0,   1,   -13.7,      -2,    0.0,     1,   0.0},
            {   2,   0,   0,   0,   1,    13.8,       2,    0.0,    -1,   0.0},
            {   3,   0,   0,   0,   0,     9.2,       2,    0.0,     0,   0.0},
            {   1,   1,   2,   0,   2,     8.9,       2,    0.0,    -1,   0.0},
            {   0,   0,   2,   1,   2,     9.3,       2,    0.0,    -1,   0.0},
            {   1,   0,   0,   2,   1,     9.6,      -1,    0.0,     0,   0.0},
            {   1,   0,   2,   2,   1,     5.6,      -1,    0.0,     1,   0.0},
            {   1,   1,   0,  -2,   1,   -34.7,      -1,    0.0,     0,   0.0},
            {   0,   1,   0,   2,   0,    14.2,      -1,    0.0,     0,   0.0},
            {   0,   1,   2,  -2,   0,   117.5,      -1,    0.0,     0,   0.0},
            {   0,   1,  -2,   2,   0,  -329.8,      -1,    0.0,     0,   0.0},
            {   1,   0,  -2,   2,   0,    23.8,      -1,    0.0,     0,   0.0},
            {   1,   0,  -2,  -2,   0,    -9.5,      -1,    0.0,     0,   0.0},
            {   1,   0,   2,  -2,   0,    32.8,      -1,    0.0,     0,   0.0},
            {   1,   0,   0,  -4,   0,   -10.1,      -1,    0.0,     0,   0.0},
            {   2,   0,   0,  -4,   0,   -15.9,      -1,    0.0,     0,   0.0},
            {   0,   0,   2,   4,   2,     4.8,      -1,    0.0,     0,   0.0},
            {   0,   0,   2,  -1,   2,    25.4,      -1,    0.0,     0,   0.0},
            {  -2,   0,   2,   4,   2,     7.3,      -1,    0.0,     1,   0.0},
            {   2,   0,   2,   2,   2,     4.7,      -1,    0.0,     0,   0.0},
            {   0,  -1,   2,   0,   1,    14.2,      -1,    0.0,     0,   0.0},
            {   0,   0,  -2,   0,   1,   -13.6,      -1,    0.0,     0,   0.0},
            {   0,   0,   4,  -2,   2,    12.7,       1,    0.0,     0,   0.0},
            {   0,   1,   0,   0,   2,   409.2,       1,    0.0,     0,   0.0},
            {   1,   1,   2,  -2,   2,    22.5,       1,    0.0,    -1,   0.0},
            {   3,   0,   2,  -2,   2,     8.7,       1,    0.0,     0,   0.0},
            {  -2,   0,   2,   2,   2,    14.6,       1,    0.0,    -1,   0.0},
            {  -1,   0,   0,   0,   2,   -27.3,       1,    0.0,    -1,   0.0},
            {   0,   0,  -2,   2,   1,  -169.0,       1,    0.0,     0,   0.0},
            {   0,   1,   2,   0,   1,    13.1,       1,    0.0,     0,   0.0},
            {  -1,   0,   4,   0,   2,     9.1,       1,    0.0,     0,   0.0},
            {   2,   1,   0,  -2,   0,   131.7,       1,    0.0,     0,   0.0},
            {   2,   0,   0,   2,   0,     7.1,       1,    0.0,     0,   0.0},
            {   2,   0,   2,  -2,   1,    12.8,       1,    0.0,    -1,   0.0},
            {   2,   0,  -2,   0,   1,  -943.2,       1,    0.0,     0,   0.0},
            {   1,  -1,   0,  -2,   0,   -29.3,       1,    0.0,     0,   0.0},
            {  -1,   0,   0,   1,   1,  -388.3,       1,    0.0,     0,   0.0},
            {  -1,  -1,   0,   2,   1,    35.0,       1,    0.0,     0,   0.0},
            {   0,   1,   0,   1,   0,    27.3,       1,    0.0,     0,   0.0}
        };
        
        static const double fc[][5]={ /* coefficients for iau 1980 nutation */
            { 134.96340251, 1717915923.2178,  31.8792,  0.051635, -0.00024470},
            { 357.52910918,  129596581.0481,  -0.5532,  0.000136, -0.00001149},
            {  93.27209062, 1739527262.8478, -12.7512, -0.001037,  0.00000417},
            { 297.85019547, 1602961601.2090,  -6.3706,  0.006593, -0.00003169},
            { 125.04455501,   -6962890.2665,   7.4722,  0.007702  -0.00005939}
        };
        
        double eps = 0.0;
        dpsi = 0.0;
        deps = 0.0;
        
        // Julian cent. since J2000 ,const double T = (TT-J2000)/86400.0/36525.0;
        const double T = JC_TT;
        eps = ( 84381.448-46.8150*T-0.00059*T*T+0.001813*T*T*T)*DAS2R;  // eps
        
        double f[5]={0.0};
        {
            double tt[4]={0.0}; tt[0] = T;
            for ( int i=1; i<4; i++) tt[i]=tt[i-1]*T;
            for (int i=0; i<5; i++)
            {
                f[i]=fc[i][0]*3600.0;
                for (int j=0; j<4; j++) f[i]+=fc[i][j+1]*tt[j];
                f[i]=fmod(f[i]*DAS2R, 2.0*PI);
            }
        }
        
        for(int i = 0; i < 106; i++)
        {
            double ang(0.0);
            for(int j=0; j<5; j++) ang+=nut[i][j]*f[j];
            
            dpsi+=(nut[i][6]+nut[i][7]*T)*std::sin(ang);
            deps+=(nut[i][8]+nut[i][9]*T)*std::cos(ang);
        }
        
        dpsi *= 1E-4*DAS2R; /* 0.1 mas -> rad */
        deps *= 1E-4*DAS2R;
    }
    
    
    /*
     * get the nutation angel according to IAU2000B
     * This is the IAU2000B Nutation model,which is less precision than IAU2000A model.
     * ref: ftp://maia.usno.navy.mil/convert/conv2003/chapter5/NU2000B.f
     */
    void GEarthOrientationParameter::getNutationAngle(double JC_TT,double& dpsi,double& deps)
    {
        iauNutation2000A(JC_TT, dpsi, deps);
    }
    
    /* get the nutation rotation matrix */
    void GEarthOrientationParameter::getNutationMatrix(double JC_TT, double *nr)
    {
        double eps =0.0,deps =0.0,dpsi =0.0;
        //double PI = GCONST("PI");
        double AS2R = GCONST("AS2R");
        
        getNutationAngle( JC_TT, deps, dpsi);
        
        
        double mobl = GIERS::getMeanObliquity(JC_TT);  // mean obliquity of ecliptic
        
        double tobl = mobl + deps;  // true obliquity of ecliptic
        
        double cobm = cos(mobl), sobm = sin(mobl);
        double cobt = cos(tobl), sobt = sin(tobl);
        double cpsi = cos(dpsi), spsi = sin(dpsi);
        
        double xx = cpsi, xy = spsi * cobt,xz = spsi * sobt;
        double yx = yx = -spsi * cobm,yy = cpsi * cobm * cobt + sobm * sobt,yz = cpsi * cobm * sobt - sobm * cobt;
        double zx = spsi * cobt,zy = cpsi * sobm * cobt - cobm * sobt,zz = cpsi * sobm * sobt + cobm * cobt;
        nr[0*3+0] = xx;nr[0*3+1] = yx;nr[0*3+2] = zx;
        nr[1*3+0] = xy;nr[1*3+1] = yy;nr[1*3+2] = zy;
        nr[2*3+0] = xz;nr[2*3+1] = yz;nr[2*3+2] = zz;
        
        
    } // end of the function
    
    // Normalize angle into the range -pi <= a < +pi.
    double GEarthOrientationParameter::normalizeAngle(double a)
    {
        double D2PI = 2.0*GCONST("PI");
        double w = fmod(a, D2PI);
        if (fabs(w) >= (D2PI*0.5))
        {
            w-= ((a<0.0)?-D2PI:D2PI);
        }
        
        return w;
    }
    
    /*
     *  The quantity s', positioning the Terrestrial Ephemeris Origin on the
     *  equator of the Celestial Intermediate Pole.
     * the locator of TIP
     *  Annex to IERS Conventions 2000, Chapter 5
     * output: radians
     */
    double GEarthOrientationParameter::getSP2000( double JC_TT )
    {
        //double PI = GCONST("PI");
        double DAS2R = GCONST("AS2R");
        double T = JC_TT;
        // the more accurate way is : sp = -0.0015(ac^2/1.2 + aa^2)*t*DAS2R
        // comes from the IERS Technical Note 36.
        // using the current aa and ac, we can get the following formular.
        double sp2000 = -47E-6*T*DAS2R;
        
        return sp2000;
    }
    
    /* get the polar motion rotation matrix in current time(m_epoch) and m_eop
     * ref: pom00.c in sofa
     */
    void GEarthOrientationParameter::getPolarMotionMatrix(double JC_TT,double *pmr)
    {
        //double PI = GCONST("PI");
        double DAS2R = GCONST("AS2R"); // arcsec to radian
        
        double xp = (m_eop.m_polarMotionX )/*0.0349282*/*DAS2R;
        double yp = (m_eop.m_polarMotionY )/*0.4833163*/*DAS2R;
        
        double sp = getSP2000(JC_TT);
        
        double xm[9]={0.0}, ym[9]={0.0}, sm[9]={0.0};
        
        //xm Ry(-xp)
        xm[0*3+0] = cos(-xp);xm[0*3+1] = 0.0;xm[0*3+2] = -sin(-xp);
        xm[1*3+0] = 0.0;xm[1*3+1] = 1.0;xm[1*3+2] = 0.0;
        xm[2*3+0] = sin(-xp);xm[2*3+1] = 0.0;xm[2*3+2] = cos(-xp);
        
        //ym Rx(-yp)
        ym[0*3+0] = 1.0;ym[0*3+1] = 0.0;ym[0*3+2] = 0.0;
        ym[1*3+0] = 0.0;ym[1*3+1] = cos(-yp);ym[1*3+2] = sin(-yp);
        ym[2*3+0] = 0.0;ym[2*3+1] = -sin(-yp);ym[2*3+2] = cos(-yp);
        
        //Rz(sp)
        sm[0*3+0] = cos(sp);sm[0*3+1] = sin(sp);sm[0*3+2] = 0.0;
        sm[1*3+0] = -sin(sp);sm[1*3+1] = cos(sp);sm[1*3+2] = 0.0;
        sm[2*3+0] = 0.0;sm[2*3+1] = 0.0;sm[2*3+2] = 1.0;
        
        double tmp[9]={0.0};
        // be careful about the order of matrix multiply
        // ym*xm*sm
        matrixMultiply(3, 3, ym, 3, xm,tmp);
        matrixMultiply(3, 3, tmp, 3, sm,pmr);
        
    } // end of the function
    
    
    
    void GEarthOrientationParameter::computeRotationMatrix( double* ECI2ECEFPos, double* ECI2ECEFVel)
    {
        double DAS2R = GCONST("AS2R");
        //time since epoch J2000, please note J2000 is expressed in TT, Terrestra Time
        GTime referenceEpoch_TT = GTime::J2000();
        //GTime referenceEpoch_UTC = GTime::TAI2UTC(GTime::TT2TAI(referenceEpoch_TT));
        double secpday = 86400.0;
        TimeSystem ts;
        long mjd =0, sod = 0;
        double fsod = 0.0;
        GTime TAI = GTime::UTC2TAI( m_epoch );
        GTime tt_epoch =  GTime::TAI2TT(TAI);
        JDTime tt_jd = GTime::GTime2JDTime(tt_epoch);
        //JDTime testjd = GTime::GTime2JDTime(TT);
        //CivilTime testct = GTime::JDTime2CivilTime(testjd);
       // GTime tt_sinceJ2000 = TT - referenceEpoch_TT; // the difference between two times
        // now the time tt is in TT since J2000
        double JD_TT_I = tt_jd.m_jd ;
        double JD_TT_F = (tt_jd.m_sod+tt_jd.m_fsod)/secpday;
        double JD_TT = JD_TT_I + JD_TT_F;
        double JC_TT = JD_TT/36525.0;
        
        //double ut1mutc_epoch = m_eop.m_UT1mUTC;  // unit:s
        GTime ut1_epoch = GTime::UTC2UT1( m_epoch, m_eop.m_UT1mUTC );
        JDTime ut1_jd = GTime::GTime2JDTime(ut1_epoch);
        double JD_UT1_I = ut1_jd.m_jd ;
        double JD_UT1_F = (ut1_jd.m_sod+ut1_jd.m_fsod)/secpday;
        
        GIAU myiau;
        
        double rc2i[3][3], era, sp, rpom[3][3],temp1[9];
        
        /* Form the celestial-to-intermediate matrix for this TT. */
        myiau.iauC2i06a(JD_TT_I, JD_TT_F, rc2i);
        
        /* Predict the Earth rotation angle for this UT1. */
        era = myiau.iauEra00(JD_UT1_I, JD_UT1_F);
        
        /*earth rotation*/
        double Rera[9] = {0.0};
        Rz(era,Rera);
        
        /* Estimate s'. */
        sp = myiau.iauSp00(JD_TT_I, JD_TT_F);
        
        /* Form the polar motion matrix. */
        myiau.iauPom00(m_eop.m_polarMotionX*DAS2R, m_eop.m_polarMotionY*DAS2R, sp, rpom);
        
        
        // W*R*Q
        // W: polar motion
        // R: earth rotation angle
        // Q: precession-nutation
        matrixMultiply( 3, 3, &rpom[0][0], 3, Rera, temp1);
        matrixMultiply( 3, 3, temp1, 3, &rc2i[0][0], ECI2ECEFPos);
        
        
        //get the first derivatives
        // starting the calculation of the VELOCITY convertion matrix
        // this method ignore the changes in nutation, precession and polar motion
        double era_rate = 0.0; //getEarthRotationAngleRate(JC_TT);
        //http://www.iers.org/IERS/EN/Science/EarthRotation/UT1LOD.html
        //Omega = 72 921 151.467064 - 0.843994809 LOD , where Omega is in picoradians/s and LOD in milliseconds.
        era_rate = (72921.151467064 - 0.843994809*m_eop.m_LOD)*1.0E-9 + m_eop.m_dOMEGA;  // lod in second, era_rate in radians per second
        
        
        double dRera[9] = {0.0}; // the first derivative of ERA matrix, Rz
        double sin_era = sin(-era);
        double cos_era = cos(-era);
        dRera[0*3+0] =  sin_era*era_rate;
        dRera[0*3+1] =  cos_era*era_rate;
        dRera[1*3+0] = -dRera[0*3+1];
        dRera[1*3+1] = dRera[0*3+0];
        // W*dR*Q
        matrixMultiply( 3, 3, &rpom[0][0], 3, dRera, temp1);
        matrixMultiply( 3, 3, temp1, 3, &rc2i[0][0], ECI2ECEFVel);
        
        
    }
    
    //compute the rotation matrix between ECI and ECEF
    // CtoT
    void GEarthOrientationParameter::computeRotationMatrix1( double* ECI2ECEFPos, double* ECI2ECEFVel)
    {
        double DAS2R = GCONST("AS2R");
        //time since epoch J2000, please note J2000 is expressed in TT, Terrestra Time
        GTime referenceEpoch_TT = GTime::J2000();
        //GTime referenceEpoch_UTC = GTime::TAI2UTC(GTime::TT2TAI(referenceEpoch_TT));
        double secpday = 86400.0;
        TimeSystem ts;
        long mjd =0, sod = 0;
        double fsod = 0.0;
        GTime TAI = GTime::UTC2TAI( m_epoch );
        GTime TT =  GTime::TAI2TT(TAI);
        //JDTime testjd = GTime::GTime2JDTime(TT);
        //CivilTime testct = GTime::JDTime2CivilTime(testjd);
        GTime tt_sinceJ2000 = TT - referenceEpoch_TT; // the difference between two times
        // now the time tt is in TT since J2000
        tt_sinceJ2000.GetData(ts, mjd, sod, fsod);
        
        double JD_TT = mjd + (sod+fsod)/secpday;
        double JC_TT = JD_TT/36525.0;
        
        //double ut1mutc_epoch = m_eop.m_UT1mUTC;  // unit:s
        GTime ut1_epoch = GTime::UTC2UT1( m_epoch, m_eop.m_UT1mUTC );
        JDTime ut1_jd = GTime::GTime2JDTime(ut1_epoch);
        double JD_UT1_I = ut1_jd.m_jd ;
        double JD_UT1_F = (ut1_jd.m_sod+ut1_jd.m_fsod)/secpday;
        
        // the time JD_UT1 is very important for the accurancy of earth rotation angle.
        // be careful of the numerical problem in it, the numerical problem in it can lead to 2-3 cm in the coordinate
        //printf("JD_UT1: %18.13f\n",JD_UT1);
        
        // the CIP based new method
        double bpn[9]={0.0}; //celestial-to-true matrix (Note 1)
        double x =0.0, y=0.0,s=0.0;
        
        getRotationBPN(JC_TT, bpn);
        
        x = bpn[2*3+0];
        y = bpn[2*3+1];
        
        
        // add corrections from eop data
        //x = x + ( m_eop.m_dX )*DAS2R ;
        //y = y + ( m_eop.m_dY )*DAS2R;
        
        //get s first, then add the corrections to x and y
        s = getLocatorS06(JC_TT, x, y);
        
        
        // according to x,y and s, forming the c2i matrix(intermidiate)
        double rc2i[9]={0.0};
        double r2 = x*x + y*y;
        double e = (r2 != 0.0) ? atan2(y, x) : 0.0;
        double d = atan(sqrt(r2 / (1.0 - r2)));
        double temp1[9]={0.0};
        double Rze[9]={0.0}, Ryd[9]={0.0},Rzes[9]={0.0};
        Rz(e, Rze);
        Ry(d, Ryd);
        Rz(-(e+s), Rzes);
        matrixMultiply(3, 3, Ryd, 3, Rze, temp1);
        matrixMultiply(3, 3, Rzes, 3,temp1, rc2i);
        
        // get the earth rotation angel in radians
        double era = GIERS::getEarthRotationAngle(JD_UT1_I,JD_UT1_F);
        
        double rera[9] = {0.0};
        Rz(era,rera);
        
        // get the ploar motion rotation
        double pm[9] = {0.0};
        getPolarMotionMatrix( JC_TT, pm );
        
        // W*R*Q
        // W: polar motion
        // R: earth rotation angle
        // Q: precession-nutation
        //matrixMultiply( 3, 3, rera, 3, rc2i, temp1);
        //matrixMultiply(3, 3, pm, 3, temp1, ECI2ECEFPos);
        matrixMultiply( 3, 3, pm, 3, rera, temp1);
        matrixMultiply( 3, 3, temp1, 3, rc2i, ECI2ECEFPos);
        
        // starting the calculation of the VELOCITY convertion matrix
        // this method ignore the changes in nutation, precession and polar motion
        double era_rate = 0.0; //getEarthRotationAngleRate(JC_TT);
        //http://www.iers.org/IERS/EN/Science/EarthRotation/UT1LOD.html
        //Omega = 72 921 151.467064 - 0.843994809 LOD , where Omega is in picoradians/s and LOD in milliseconds.
        era_rate = (72921.151467064 - 0.843994809*m_eop.m_LOD)*1.0E-9 + m_eop.m_dOMEGA;  // lod in second, era_rate in radians per second
        
        
        double dRera[9] = {0.0}; // the first derivative of ERA matrix, Rz
        double sin_era = sin(-era);
        double cos_era = cos(-era);
        dRera[0*3+0] =  sin_era*era_rate;
        dRera[0*3+1] =  cos_era*era_rate;
        dRera[1*3+0] = -dRera[0*3+1];
        dRera[1*3+1] = dRera[0*3+0];
        // W*dR*Q
        matrixMultiply( 3, 3, pm, 3, dRera, temp1);
        matrixMultiply( 3, 3, temp1, 3, rc2i, ECI2ECEFVel);
        
    }
    
    /*convert from ECI to ECEF
     * tag =0 just position
     * tag =1 just velocity
     * tag =2 both position and velocity
     * eci[0], eci[1], eci[2] are positions in meters while eci[3] eci[4] eci[5] are velocity in meter/s
     */
    void GEarthOrientationParameter::ECI2ECEF(int tag, double* eci,double* ecef)
    {
        double rotationMatrixPos[9]={0.0};
        double rotationMatrixVel[9]={0.0};
        memcpy(rotationMatrixPos, m_eci2ecefPos, sizeof(double)*9);
        memcpy(rotationMatrixVel, m_eci2ecefVel, sizeof(double)*9);
        //computeRotationMatrix(rotationMatrixPos, rotationMatrixVel);
        
        if( tag == 0 ) // just position
        {
            matrixMultiply(3, 3, rotationMatrixPos, 1, eci, ecef);
        }
        else if(tag == 1)
        {
            double tmp1[3]= {0.0}, tmp2[3]= {0.0},tmp3[3]={0.0};
            tmp1[0] = eci[0];tmp1[1] = eci[1];tmp1[2] = eci[2];
            matrixMultiply(3, 3, rotationMatrixPos, 1, tmp1, tmp2);
            ecef[0] = tmp2[0];ecef[1] = tmp2[1];ecef[2] = tmp2[2];
            
            //velocity two parts:
            tmp1[0] = eci[3];tmp1[1] = eci[4];tmp1[2] = eci[5];
            matrixMultiply(3, 3, rotationMatrixPos, 1, tmp1, tmp2);
            tmp1[0] = eci[0];tmp1[1] = eci[1];tmp1[2] = eci[2];
            matrixMultiply(3, 3, rotationMatrixVel, 1, tmp1, tmp3);
            
            ecef[3] = tmp2[0]+tmp3[0];ecef[4] = tmp2[1]+tmp3[1];ecef[5] = tmp2[2]+tmp3[2];
            
        }
    }
    
    void GEarthOrientationParameter::ECI2ECEF_vel(GVector& pos_eci, GVector& vel_eci, GVector& vel_ecef )
    {
        vel_ecef.x =   m_eci2ecefPos[0]*vel_eci.x + m_eci2ecefPos[1]*vel_eci.y + m_eci2ecefPos[2]*vel_eci.z ;
        vel_ecef.y =   m_eci2ecefPos[3]*vel_eci.x + m_eci2ecefPos[4]*vel_eci.y + m_eci2ecefPos[5]*vel_eci.z ;
        vel_ecef.z =   m_eci2ecefPos[6]*vel_eci.x + m_eci2ecefPos[7]*vel_eci.y + m_eci2ecefPos[8]*vel_eci.z ;
        
        vel_ecef.x += m_eci2ecefVel[0]*pos_eci.x + m_eci2ecefVel[1]*pos_eci.y + m_eci2ecefVel[2]*pos_eci.z ;
        vel_ecef.y += m_eci2ecefVel[3]*pos_eci.x + m_eci2ecefVel[4]*pos_eci.y + m_eci2ecefVel[5]*pos_eci.z ;
        vel_ecef.z += m_eci2ecefVel[6]*pos_eci.x + m_eci2ecefVel[7]*pos_eci.y + m_eci2ecefVel[8]*pos_eci.z ;
    }
    
    
    void  GEarthOrientationParameter::ECI2ECEF_pos( GVector& pos_eci, GVector& pos_ecef )
    {
        pos_ecef.x =  m_eci2ecefPos[0]*pos_eci.x + m_eci2ecefPos[1]*pos_eci.y + m_eci2ecefPos[2]*pos_eci.z ;
        
        pos_ecef.y =  m_eci2ecefPos[3]*pos_eci.x + m_eci2ecefPos[4]*pos_eci.y + m_eci2ecefPos[5]*pos_eci.z ;
        
        pos_ecef.z =  m_eci2ecefPos[6]*pos_eci.x + m_eci2ecefPos[7]*pos_eci.y + m_eci2ecefPos[8]*pos_eci.z ;
        
    }
    
    void GEarthOrientationParameter::ECEF2ECI_pos( GVector& pos_ecef,GVector& pos_eci )
    {
        pos_eci.x =  m_eci2ecefPos[0]*pos_ecef.x + m_eci2ecefPos[3]*pos_ecef.y + m_eci2ecefPos[6]*pos_ecef.z ;
        
        pos_eci.y =  m_eci2ecefPos[1]*pos_ecef.x + m_eci2ecefPos[4]*pos_ecef.y + m_eci2ecefPos[7]*pos_ecef.z ;
        
        pos_eci.z =  m_eci2ecefPos[2]*pos_ecef.x + m_eci2ecefPos[5]*pos_ecef.y + m_eci2ecefPos[8]*pos_ecef.z ;
    }
    
    void GEarthOrientationParameter::ECEF2ECI_vel(GVector& pos_ecef, GVector& vel_ecef, GVector& vel_eci )
    {
        vel_eci.x = vel_eci.y = vel_eci.z =0.0;
        vel_eci.x =   m_eci2ecefPos[0]*vel_ecef.x + m_eci2ecefPos[3]*vel_ecef.y + m_eci2ecefPos[6]*vel_ecef.z ;
        vel_eci.y =   m_eci2ecefPos[1]*vel_ecef.x + m_eci2ecefPos[4]*vel_ecef.y + m_eci2ecefPos[7]*vel_ecef.z ;
        vel_eci.z =   m_eci2ecefPos[2]*vel_ecef.x + m_eci2ecefPos[5]*vel_ecef.y + m_eci2ecefPos[8]*vel_ecef.z ;
        
        vel_eci.x += m_eci2ecefVel[0]*pos_ecef.x + m_eci2ecefVel[3]*pos_ecef.y + m_eci2ecefVel[6]*pos_ecef.z ;
        vel_eci.y += m_eci2ecefVel[1]*pos_ecef.x + m_eci2ecefVel[4]*pos_ecef.y + m_eci2ecefVel[7]*pos_ecef.z ;
        vel_eci.z += m_eci2ecefVel[2]*pos_ecef.x + m_eci2ecefVel[5]*pos_ecef.y + m_eci2ecefVel[8]*pos_ecef.z ;
        
    }
    
    
    /*convert from ECEF to ECI
     * tag =0 just position
     * tag =1 just velocity
     * tag =2 both position and velocity
     * eci[0], eci[1], eci[2] are positions in meters while eci[3] eci[4] eci[5] are velocity in meter/s
     */
    void GEarthOrientationParameter::ECEF2ECI(int tag, double* ecef,double* eci)
    {
        double rotationMatrixPos[9]={0.0};
        double rotationMatrixVel[9]={0.0};
        memcpy(rotationMatrixPos, m_eci2ecefPos, sizeof(double)*9);
        memcpy(rotationMatrixVel, m_eci2ecefVel, sizeof(double)*9);
        
        // need to transpose the rotation matrix , because the transpose is equal to inversion
        matrixTranspose(3, 3, rotationMatrixPos);
        matrixTranspose(3, 3, rotationMatrixVel);
        
        if( tag == 0 ) // just position
        {
            matrixMultiply(3, 3, rotationMatrixPos, 1, ecef, eci);
        }
        else if(tag == 1) // position and velocity
        {
            double tmp1[3]= {0.0}, tmp2[3]= {0.0}, tmp3[3]={0.0};
            tmp1[0] = ecef[0];tmp1[1] = ecef[1];tmp1[2] = ecef[2];
            matrixMultiply(3, 3, rotationMatrixPos, 1, tmp1, tmp2);
            eci[0] = tmp2[0];eci[1] = tmp2[1];eci[2] = tmp2[2];
            
            // the velocity , two parts
            tmp1[0] = ecef[3];tmp1[1] = ecef[4];tmp1[2] = ecef[5];
            matrixMultiply(3, 3, rotationMatrixPos, 1, tmp1, tmp2);
            
            tmp1[0] = ecef[0];tmp1[1] = ecef[1];tmp1[2] = ecef[2];
            matrixMultiply(3, 3, rotationMatrixVel, 1, tmp1, tmp3);
            
            eci[3] = tmp2[0]+tmp3[0];eci[4] = tmp2[1]+tmp3[1];eci[5] = tmp2[2]+tmp3[2];
        }
    }
    
    void  GEarthOrientationParameter::getECI2ECEFMatrix(double* tm)
    {
        memcpy(tm,m_eci2ecefPos,sizeof(double)*9);
    }
    
    
    /*
     To compute the forced nutation of the non-rigid Earth based on
     the IAU 2000B precession/nutation model.
     
     REFERENCES:
     McCarthy, D. and Luzum, B. (2003). "An Abridged Model of the
     Precession & Nutation of the Celestial Pole," Celestial
     Mechanics and Dynamical Astronomy, Volume 85, Issue 1,
     Jan. 2003, p. 37. (IAU 2000B)
     IERS Conventions (2003), Chapter 5.
     
     INPUT
     ARGUMENTS: TT Julian date since J2000
     
     OUTPUT
     ARGUMENTS:
     *dpsi (double)
     Nutation (luni-solar + planetary) in longitude, in radians.
     *deps (double)
     Nutation (luni-solar + planetary) in obliquity, in radians.
     */
    void GEarthOrientationParameter::iauNutation2000B(double JC_TT,double& dpsi, double& deps)
    {
        short int i=0;
        
        /*
         Planetary nutation (arcsec).  These fixed terms account for the
         omission of the long-period planetary terms in the truncated model.
         */
        
        double dpplan = -0.000135;
        double deplan =  0.000388;
        
        double T, el, elp, f, d, om, arg, dp, de, sarg, carg, factor, dpsils,
        depsls, dpsipl, depspl;
        double AS2R = GCONST("AS2R");
        
        double TWOPI = GCONST("PI")*2.0;
        
        /*
         Luni-Solar argument multipliers:
         L     L'    F     D     Om
         */
        
        static const short int nals_t[77][5] =
        {
            { 0,    0,    0,    0,    1},
            { 0,    0,    2,   -2,    2},
            { 0,    0,    2,    0,    2},
            { 0,    0,    0,    0,    2},
            { 0,    1,    0,    0,    0},
            { 0,    1,    2,   -2,    2},
            { 1,    0,    0,    0,    0},
            { 0,    0,    2,    0,    1},
            { 1,    0,    2,    0,    2},
            { 0,   -1,    2,   -2,    2},
            { 0,    0,    2,   -2,    1},
            {-1,    0,    2,    0,    2},
            {-1,    0,    0,    2,    0},
            { 1,    0,    0,    0,    1},
            {-1,    0,    0,    0,    1},
            {-1,    0,    2,    2,    2},
            { 1,    0,    2,    0,    1},
            {-2,    0,    2,    0,    1},
            { 0,    0,    0,    2,    0},
            { 0,    0,    2,    2,    2},
            { 0,   -2,    2,   -2,    2},
            {-2,    0,    0,    2,    0},
            { 2,    0,    2,    0,    2},
            { 1,    0,    2,   -2,    2},
            {-1,    0,    2,    0,    1},
            { 2,    0,    0,    0,    0},
            { 0,    0,    2,    0,    0},
            { 0,    1,    0,    0,    1},
            {-1,    0,    0,    2,    1},
            { 0,    2,    2,   -2,    2},
            { 0,    0,   -2,    2,    0},
            { 1,    0,    0,   -2,    1},
            { 0,   -1,    0,    0,    1},
            {-1,    0,    2,    2,    1},
            { 0,    2,    0,    0,    0},
            { 1,    0,    2,    2,    2},
            {-2,    0,    2,    0,    0},
            { 0,    1,    2,    0,    2},
            { 0,    0,    2,    2,    1},
            { 0,   -1,    2,    0,    2},
            { 0,    0,    0,    2,    1},
            { 1,    0,    2,   -2,    1},
            { 2,    0,    2,   -2,    2},
            {-2,    0,    0,    2,    1},
            { 2,    0,    2,    0,    1},
            { 0,   -1,    2,   -2,    1},
            { 0,    0,    0,   -2,    1},
            {-1,   -1,    0,    2,    0},
            { 2,    0,    0,   -2,    1},
            { 1,    0,    0,    2,    0},
            { 0,    1,    2,   -2,    1},
            { 1,   -1,    0,    0,    0},
            {-2,    0,    2,    0,    2},
            { 3,    0,    2,    0,    2},
            { 0,   -1,    0,    2,    0},
            { 1,   -1,    2,    0,    2},
            { 0,    0,    0,    1,    0},
            {-1,   -1,    2,    2,    2},
            {-1,    0,    2,    0,    0},
            { 0,   -1,    2,    2,    2},
            {-2,    0,    0,    0,    1},
            { 1,    1,    2,    0,    2},
            { 2,    0,    0,    0,    1},
            {-1,    1,    0,    1,    0},
            { 1,    1,    0,    0,    0},
            { 1,    0,    2,    0,    0},
            {-1,    0,    2,   -2,    1},
            { 1,    0,    0,    0,    2},
            {-1,    0,    0,    1,    0},
            { 0,    0,    2,    1,    2},
            {-1,    0,    2,    4,    2},
            {-1,    1,    0,    1,    1},
            { 0,   -2,    2,   -2,    1},
            { 1,    0,    2,    2,    1},
            {-2,    0,    2,    2,    2},
            {-1,    0,    0,    0,    2},
            { 1,    1,    2,   -2,    2}};
        
        /*
         Luni-Solar nutation coefficients, unit 1e-7 arcsec:
         longitude (sin, t*sin, cos), obliquity (cos, t*cos, sin)
         
         Each row of coefficients in 'cls_t' belongs with the corresponding
         row of fundamental-argument multipliers in 'nals_t'.
         */
        
        static const double cls_t[77][6] = {
            {-172064161.0, -174666.0,  33386.0, 92052331.0,  9086.0, 15377.0},
            { -13170906.0,   -1675.0, -13696.0,  5730336.0, -3015.0, -4587.0},
            {  -2276413.0,    -234.0,   2796.0,   978459.0,  -485.0,  1374.0},
            {   2074554.0,     207.0,   -698.0,  -897492.0,   470.0,  -291.0},
            {   1475877.0,   -3633.0,  11817.0,    73871.0,  -184.0, -1924.0},
            {   -516821.0,    1226.0,   -524.0,   224386.0,  -677.0,  -174.0},
            {    711159.0,      73.0,   -872.0,    -6750.0,     0.0,   358.0},
            {   -387298.0,    -367.0,    380.0,   200728.0,    18.0,   318.0},
            {   -301461.0,     -36.0,    816.0,   129025.0,   -63.0,   367.0},
            {    215829.0,    -494.0,    111.0,   -95929.0,   299.0,   132.0},
            {    128227.0,     137.0,    181.0,   -68982.0,    -9.0,    39.0},
            {    123457.0,      11.0,     19.0,   -53311.0,    32.0,    -4.0},
            {    156994.0,      10.0,   -168.0,    -1235.0,     0.0,    82.0},
            {     63110.0,      63.0,     27.0,   -33228.0,     0.0,    -9.0},
            {    -57976.0,     -63.0,   -189.0,    31429.0,     0.0,   -75.0},
            {    -59641.0,     -11.0,    149.0,    25543.0,   -11.0,    66.0},
            {    -51613.0,     -42.0,    129.0,    26366.0,     0.0,    78.0},
            {     45893.0,      50.0,     31.0,   -24236.0,   -10.0,    20.0},
            {     63384.0,      11.0,   -150.0,    -1220.0,     0.0,    29.0},
            {    -38571.0,      -1.0,    158.0,    16452.0,   -11.0,    68.0},
            {     32481.0,       0.0,      0.0,   -13870.0,     0.0,     0.0},
            {    -47722.0,       0.0,    -18.0,      477.0,     0.0,   -25.0},
            {    -31046.0,      -1.0,    131.0,    13238.0,   -11.0,    59.0},
            {     28593.0,       0.0,     -1.0,   -12338.0,    10.0,    -3.0},
            {     20441.0,      21.0,     10.0,   -10758.0,     0.0,    -3.0},
            {     29243.0,       0.0,    -74.0,     -609.0,     0.0,    13.0},
            {     25887.0,       0.0,    -66.0,     -550.0,     0.0,    11.0},
            {    -14053.0,     -25.0,     79.0,     8551.0,    -2.0,   -45.0},
            {     15164.0,      10.0,     11.0,    -8001.0,     0.0,    -1.0},
            {    -15794.0,      72.0,    -16.0,     6850.0,   -42.0,    -5.0},
            {     21783.0,       0.0,     13.0,     -167.0,     0.0,    13.0},
            {    -12873.0,     -10.0,    -37.0,     6953.0,     0.0,   -14.0},
            {    -12654.0,      11.0,     63.0,     6415.0,     0.0,    26.0},
            {    -10204.0,       0.0,     25.0,     5222.0,     0.0,    15.0},
            {     16707.0,     -85.0,    -10.0,      168.0,    -1.0,    10.0},
            {     -7691.0,       0.0,     44.0,     3268.0,     0.0,    19.0},
            {    -11024.0,       0.0,    -14.0,      104.0,     0.0,     2.0},
            {      7566.0,     -21.0,    -11.0,    -3250.0,     0.0,    -5.0},
            {     -6637.0,     -11.0,     25.0,     3353.0,     0.0,    14.0},
            {     -7141.0,      21.0,      8.0,     3070.0,     0.0,     4.0},
            {     -6302.0,     -11.0,      2.0,     3272.0,     0.0,     4.0},
            {      5800.0,      10.0,      2.0,    -3045.0,     0.0,    -1.0},
            {      6443.0,       0.0,     -7.0,    -2768.0,     0.0,    -4.0},
            {     -5774.0,     -11.0,    -15.0,     3041.0,     0.0,    -5.0},
            {     -5350.0,       0.0,     21.0,     2695.0,     0.0,    12.0},
            {     -4752.0,     -11.0,     -3.0,     2719.0,     0.0,    -3.0},
            {     -4940.0,     -11.0,    -21.0,     2720.0,     0.0,    -9.0},
            {      7350.0,       0.0,     -8.0,      -51.0,     0.0,     4.0},
            {      4065.0,       0.0,      6.0,    -2206.0,     0.0,     1.0},
            {      6579.0,       0.0,    -24.0,     -199.0,     0.0,     2.0},
            {      3579.0,       0.0,      5.0,    -1900.0,     0.0,     1.0},
            {      4725.0,       0.0,     -6.0,      -41.0,     0.0,     3.0},
            {     -3075.0,       0.0,     -2.0,     1313.0,     0.0,    -1.0},
            {     -2904.0,       0.0,     15.0,     1233.0,     0.0,     7.0},
            {      4348.0,       0.0,    -10.0,      -81.0,     0.0,     2.0},
            {     -2878.0,       0.0,      8.0,     1232.0,     0.0,     4.0},
            {     -4230.0,       0.0,      5.0,      -20.0,     0.0,    -2.0},
            {     -2819.0,       0.0,      7.0,     1207.0,     0.0,     3.0},
            {     -4056.0,       0.0,      5.0,       40.0,     0.0,    -2.0},
            {     -2647.0,       0.0,     11.0,     1129.0,     0.0,     5.0},
            {     -2294.0,       0.0,    -10.0,     1266.0,     0.0,    -4.0},
            {      2481.0,       0.0,     -7.0,    -1062.0,     0.0,    -3.0},
            {      2179.0,       0.0,     -2.0,    -1129.0,     0.0,    -2.0},
            {      3276.0,       0.0,      1.0,       -9.0,     0.0,     0.0},
            {     -3389.0,       0.0,      5.0,       35.0,     0.0,    -2.0},
            {      3339.0,       0.0,    -13.0,     -107.0,     0.0,     1.0},
            {     -1987.0,       0.0,     -6.0,     1073.0,     0.0,    -2.0},
            {     -1981.0,       0.0,      0.0,      854.0,     0.0,     0.0},
            {      4026.0,       0.0,   -353.0,     -553.0,     0.0,  -139.0},
            {      1660.0,       0.0,     -5.0,     -710.0,     0.0,    -2.0},
            {     -1521.0,       0.0,      9.0,      647.0,     0.0,     4.0},
            {      1314.0,       0.0,      0.0,     -700.0,     0.0,     0.0},
            {     -1283.0,       0.0,      0.0,      672.0,     0.0,     0.0},
            {     -1331.0,       0.0,      8.0,      663.0,     0.0,     4.0},
            {      1383.0,       0.0,     -2.0,     -594.0,     0.0,    -2.0},
            {      1405.0,       0.0,      4.0,     -610.0,     0.0,     2.0},
            {      1290.0,       0.0,      0.0,     -556.0,     0.0,     0.0}};
        
        /*
         Interval between fundamental epoch J2000.0 and given date.
         */
        T = JC_TT;
        
        /*
         ** Luni-solar nutation. **
         
         Fundamental (Delaunay) arguments from Simon et al. (1994),
         in radians.
         */
        
        /*
         Mean anomaly of the Moon.
         */
        
        el  = fmod (485868.249036 +
                    T * 1717915923.2178, 1296000.0) * AS2R;
        
        /*
         Mean anomaly of the Sun.
         */
        
        elp = fmod (1287104.79305 +
                    T * 129596581.0481, 1296000.0) * AS2R;
        
        /*
         Mean argument of the latitude of the Moon.
         */
        
        f   = fmod (335779.526232 +
                    T * 1739527262.8478, 1296000.0) * AS2R;
        
        /*
         Mean elongation of the Moon from the Sun.
         */
        
        d   = fmod (1072260.70369 +
                    T * 1602961601.2090, 1296000.0) * AS2R;
        
        /*
         Mean longitude of the ascending node of the Moon.
         */
        
        om  = fmod (450160.398036 -
                    T * 6962890.5431, 1296000.0) * AS2R;
        
        /*
         Initialize the nutation values.
         */
        
        dp = 0.0;
        de = 0.0;
        
        /*
         Summation of luni-solar nutation series (in reverse order).
         */
        
        for (i = 76; i >= 0; i--)
        {
            
            /*
             Argument and functions.
             */
            
            arg = fmod ((double) nals_t[i][0] * el  +
                        (double) nals_t[i][1] * elp +
                        (double) nals_t[i][2] * f   +
                        (double) nals_t[i][3] * d   +
                        (double) nals_t[i][4] * om,   TWOPI);
            
            sarg = sin (arg);
            carg = cos (arg);
            
            /*
             Term.
             */
            
            dp += (cls_t[i][0] + cls_t[i][1] * T) * sarg
            +   cls_t[i][2] * carg;
            de += (cls_t[i][3] + cls_t[i][4] * T) * carg
            +   cls_t[i][5] * sarg;
        }
        
        /*
         Convert from 0.1 microarcsec units to radians.
         */
        
        factor = 1.0E-7 * AS2R;
        dpsils = dp * factor;
        depsls = de * factor;
        
        /*
         ** Planetary nutation. **
         
         Fixed terms to allow for long-period nutation, in radians.
         */
        
        dpsipl = dpplan * AS2R;
        depspl = deplan * AS2R;
        
        /*
         Total: Add planetary and luni-solar components.
         */
        
        dpsi = dpsipl + dpsils;
        deps = depspl + depsls;
        
    }
    

    
    
}// end of namespace


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

//
//  GMotionState.cpp
//  GFC
//
//  Created by lizhen on 05/06/2016.
//  Copyright © 2016 lizhen. All rights reserved.
//

#include "GMotionState.hpp"

namespace gfc
{
    
    GMotionState::GMotionState()
    {
        phiMatrix.resize(6, 6);
        phiMatrix(0,0) = 1.0;
        phiMatrix(1,1) = 1.0;
        phiMatrix(2,2) = 1.0;
        phiMatrix(3,3) = 1.0;
        phiMatrix(4,4) = 1.0;
        phiMatrix(5,5) = 1.0;
        
        dis_sat = 0.0;    // the distance from the sat to earth center
        dis_sat_sun = 0.0;  // in km
        eps = 0.0;  // the earth-probe-sun angle in radian
        shadow_factor = 1.0;
    }
    
    void GMotionState::updateState_ecef( GTime epoch_utc, GVector psatpos_ecef, GVector psatvel_ecef )
    {
        if( m_epoch != epoch_utc )
        {
            double PI = GCONST("PI");
            
            m_epoch = epoch_utc;
            
            JDTime jdt =  GTime::GTime2JDTime(m_epoch);
            
            ct = GTime::JDTime2CivilTime(jdt);
            
            satpos_ecef = psatpos_ecef;
            
            satvel_ecef = psatvel_ecef;
            
            GSpaceEnv::eop.ECEF2ECI_pos(satpos_ecef, satpos_eci);
            
            GSpaceEnv::eop.ECEF2ECI_vel(satpos_ecef, satvel_ecef, satvel_eci);
            
            keplerianElement.PV2KP(satpos_eci, satvel_eci);
            
            dis_sat = satpos_eci.norm();
            
            satposHat_eci = satpos_eci / dis_sat;
            
            satposHat_ecef = satpos_ecef/ dis_sat;
            
            satvelHat_eci = normalise(satvel_eci);
            
            satvelHat_ecef = normalise(satvel_ecef);
            
            sat_sun_eci = GSpaceEnv::planetPos_eci[GJPLEPH::SUN] - satpos_eci ;
            
            dis_sat_sun = sat_sun_eci.norm();
            
            sat_sunHat_eci = sat_sun_eci/dis_sat_sun;
            
            sat_sun_ecef = GSpaceEnv::planetPos_ecef[GJPLEPH::SUN] - satpos_ecef;
            
            sat_sunHat_ecef = sat_sun_ecef/ dis_sat_sun;
            
            //GVector mysunpos_eci = GSpaceEnv::planetPos_eci[GJPLEPH::SUN];
            
            // get the eps angle
            eps = acos( dotproduct( sat_sun_eci , -satpos_eci)/dis_sat/dis_sat_sun );
            GVector sunpos_u_eci = normalise( GSpaceEnv::planetPos_eci[GJPLEPH::SUN] );
            
            //double shadowfactor = 1.0;
            //shadow_factor = shadowFactor(GSpaceEnv::planetPos_ecef[GJPLEPH::SUN], satpos_ecef);
            //shadow_facor = 1.0;
            //GTime gpst = GTime::UTC2GPST(epoch_utc);
            //printf("%s %6.3f %6.3f ",GTime::GTime2CivilTime(gpst).TimeString().c_str(),beta*180.0/3.14159265357,eps*180.0/3.14159265357);
            
            //double testShadow = GSpaceCraftAttitude::eclipse( GSpaceEnv::planetPos_ecef[GJPLEPH::SUN], satpos_ecef );
            //shadow_factor = shadowFactor(GSpaceEnv::planetPos_ecef[GJPLEPH::SUN], satpos_ecef);
           
            
            //shadow_factor = shadowFunction(GSpaceEnv::planetPos_ecef[GJPLEPH::SUN], satpos_ecef);
            
            shadow_factor = shadowFactor_SECM(true,GSpaceEnv::planetPos_ecef[GJPLEPH::SUN], satpos_ecef);
            
//            double testshadow = shadowFunction(GSpaceEnv::planetPos_ecef[GJPLEPH::SUN], satpos_ecef);
//            
//            if(shadow_factor < 1.0)
//            {
//               printf("%s %.6f %.6f\n",GTime::GTime2CivilTime(GTime::UTC2GPST(epoch_utc)).TimeString().c_str(),shadow_factor, testshadow);
//            }
            
            //shadow_factor = 1.0;
            
            
            //shadow_detector(GSpaceEnv::planetPos_ecef[GJPLEPH::SUN], satpos_ecef);
            
           
            
            //printf("%.6f %.6f %.6f %.6f %.6f %.6f\n",satpos_ecef.x,satpos_ecef.y,satpos_ecef.z,
            //       GSpaceEnv::planetPos_ecef[GJPLEPH::SUN].x,
            //       GSpaceEnv::planetPos_ecef[GJPLEPH::SUN].y,
            //       GSpaceEnv::planetPos_ecef[GJPLEPH::SUN].z);

             dis_factor = (GCONST("AU")/dis_sat_sun)*(GCONST("AU")/dis_sat_sun);
            
            //update attitude
            attitude_eci.updateAttitude(satpos_eci,GSpaceEnv::planetPos_eci[GJPLEPH::SUN] , sat_sun_eci, satvel_eci,shadow_factor);
            
            // convert attitude_eci to attitude_ecef
            attitude_ecef = attitude_eci;
            attitude_ecef.convert2ECEF(GSpaceEnv::eop);
            
            
            //calculate the orbit ascending node vector
            orbit_node_vector_eci = crossproduct(GVector(0,0,1.0), attitude_eci.n);
            orbit_node_vector_eci.normalise();
            double cosu = dotproduct(orbit_node_vector_eci, satposHat_eci);
            
            GVector tmpA = crossproduct(orbit_node_vector_eci, satposHat_eci);
            double sinu = tmpA.norm();
            
            //u = atan2(sinu, cosu);
            //u = keplerianElement.m_argp + keplerianElement.m_tran;
            
            //update attitude
            if( earthFlux.m_flux.size() > 0 )
            {
                earthFlux.m_flux.clear();
            }
            
            //clock_t start, end;
            //start = clock();
            
            //update the flux value, which will use the eclipse state  of the satellite
            earthFlux.makeFlux( ct.m_month - 1, GSpaceEnv::planetPos_ecef[GJPLEPH::SUN] , satpos_ecef, dis_sat );
            
            for( int i = 0 ; i< earthFlux.m_flux.size(); i++ )
            {
                GVector t ;
                GSpaceEnv::eop.ECEF2ECI_pos(earthFlux.m_flux[i].m_dir, t);
                earthFlux.m_flux[i].m_dir = t;
            }

            //end = clock();
            
            //double time = double(end - start)/CLOCKS_PER_SEC;
            
            //update the flux value, which will use the eclipse state  of the satellite, with eci
            //earthFlux.makeFlux( ct.m_month - 1, GSpaceEnv::planetPos_eci[GJPLEPH::SUN] , satpos_eci, dis_sat );
            
            solarFlux.makeFlux( GSpaceEnv::tsi*shadow_factor*dis_factor, sat_sunHat_eci, dis_sat_sun);
            
        }
        
    }
    
    /*
     *
     * update the state
     */
    void GMotionState::updateState_eci(GTime epoch_utc, GVector psatpos_eci, GVector psatvel_eci )
    {
        
        if( epoch_utc != m_epoch )
        {
            double PI = GCONST("PI");
            
            m_epoch = epoch_utc;
            
            JDTime jdt =  GTime::GTime2JDTime(m_epoch);
            
            ct = GTime::JDTime2CivilTime(jdt);
            
            satpos_eci = psatpos_eci;
            satvel_eci = psatvel_eci;
            
            keplerianElement.PV2KP(satpos_eci, satvel_eci);
            
            GSpaceEnv::eop.ECI2ECEF_pos(satpos_eci, satpos_ecef);
            GSpaceEnv::eop.ECI2ECEF_vel(satpos_eci, satvel_eci, satvel_ecef);
            
            dis_sat = satpos_eci.norm();
            
            satposHat_eci = satpos_eci / dis_sat;
            
            satposHat_ecef = satpos_ecef/ dis_sat;
            
            satvelHat_eci = normalise(satvel_eci);
            
            satvelHat_ecef = normalise(satvel_ecef);
            //sat_sun_eci is from sat to sun
            sat_sun_eci = GSpaceEnv::planetPos_eci[GJPLEPH::SUN] - satpos_eci ;
            dis_sat_sun = sat_sun_eci.norm();
            sat_sunHat_eci = sat_sun_eci/dis_sat_sun;
            
            sat_sun_ecef = GSpaceEnv::planetPos_ecef[GJPLEPH::SUN] - satpos_ecef;
            sat_sunHat_ecef = sat_sun_ecef/ dis_sat_sun;
            
            dis_factor = (GCONST("AU")/dis_sat_sun)*(GCONST("AU")/dis_sat_sun);
            
            // test the eclipse state first,then the attitude
            //double shadowfactor = 1.0;
            // shadow_factor = shadowFactor(GSpaceEnv::planetPos_ecef[GJPLEPH::SUN], satpos_ecef);
            //shadow_factor = shadowFunction(GSpaceEnv::planetPos_ecef[GJPLEPH::SUN], satpos_ecef);
            
            //shadow_factor = 1.0;
            
            
            shadow_factor = shadowFactor_SECM(true,GSpaceEnv::planetPos_ecef[GJPLEPH::SUN], satpos_ecef);
            
            //shadow_factor = shadowFunction(GSpaceEnv::planetPos_ecef[GJPLEPH::SUN], satpos_ecef);
            
//            double testshadow = shadowFunction(GSpaceEnv::planetPos_ecef[GJPLEPH::SUN], satpos_ecef);
//            
//            if(shadow_factor <1.0)
//            {
//                printf("%s %.6f %.6f\n",GTime::GTime2CivilTime(GTime::UTC2GPST(epoch_utc)).TimeString().c_str(),shadow_factor, testshadow);
//            }
            
            //shadow_factor = 1.0;
            
            
            //shadow_detector(GSpaceEnv::planetPos_ecef[GJPLEPH::SUN], satpos_ecef);
            
    //printf("%s %6.3f ",GTime::GTime2CivilTime(GTime::UTC2GPST(epoch_utc)).TimeString().c_str(),shadow_factor);
            
            
            //printf("%10.8f\n",shadowfactor);
            
            //update attitude, always use ECI to update the attitude
            //shadow_facor = 1.0;
            attitude_eci.updateAttitude(satpos_eci,GSpaceEnv::planetPos_eci[GJPLEPH::SUN] , sat_sun_eci, satvel_eci,shadow_factor);
            
            attitude_ecef = attitude_eci;
            attitude_ecef.convert2ECEF(GSpaceEnv::eop);
            
            //attitude_ecef.updateAttitude(satpos_ecef, GSpaceEnv::planetPos_ecef[GJPLEPH::SUN], sat_sun_ecef, satvel_ecef,shadowfactor);
            
            //calculate the orbit ascending node vector
            GVector z_eci(0,0,1);
            orbit_node_vector_eci = crossproduct(z_eci, attitude_eci.n);
            orbit_node_vector_eci.normalise();
            double cosu = dotproduct(orbit_node_vector_eci, satposHat_eci);
            
            GVector tmpA = crossproduct(orbit_node_vector_eci, satposHat_eci);
            double sinu = tmpA.norm();
            
            //u = atan2( sinu, cosu );
            
            //u = keplerianElement.m_argp + keplerianElement.m_tran;
            
            // get the eps angle
            eps = acos( -dotproduct(sat_sun_eci, satpos_eci)/(dis_sat_sun*dis_sat) );
            
            //update attitude
            if( earthFlux.m_flux.size() > 0 )
            {
                earthFlux.m_flux.clear();
            }
            
            //update the flux value, which will use the eclipse state  of the satellite, with ecef, use ct.m_month
            earthFlux.makeFlux( ct.m_month - 1, GSpaceEnv::planetPos_ecef[GJPLEPH::SUN] , satpos_ecef, dis_sat );
            
            //update the flux value, which will use the eclipse state  of the satellite, with eci
            //earthFlux.makeFlux( ct.m_month - 1, GSpaceEnv::planetPos_eci[GJPLEPH::SUN] , satpos_eci, dis_sat );
            
            for( int i = 0 ; i< earthFlux.m_flux.size(); i++ )
            {
                GVector t ;
                GSpaceEnv::eop.ECEF2ECI_pos(earthFlux.m_flux[i].m_dir, t);
                earthFlux.m_flux[i].m_dir = t;
            }
            
            solarFlux.makeFlux( GSpaceEnv::tsi*shadow_factor*dis_factor, sat_sunHat_eci, dis_sat_sun);
            
        }
    }
    
    void GMotionState::updateTransitionMatrix(double *phi)
    {
        phiMatrix.setData(phi, 6, 6);
    }
    void GMotionState::updateTransitionMatrix(GMatrix phi)
    {
        phiMatrix = phi;
    }
    
    
    void GMotionState::updateSensitivityMatrix(double* S)
    {
        int np = senMatrix.getColNO();
        senMatrix.setData(S, 6, np);
    }
    
    void GMotionState::updateSensitivityMatrix(GMatrix S)
    {
        senMatrix = S;
    }
    
    GSpaceCraftAttitude* GMotionState::getAttitude()
    {
        return &attitude_eci;
    }
    
    GSpaceCraftAttitude* GMotionState::getAttitude_ecef()
    {
        GSpaceCraftAttitude atti_ecef;
        
        return &atti_ecef;
    }
    
    void GMotionState::outputState( ofstream& file )
    {
        static bool first = true;
    }
    
    void GMotionState::shadow_detector(GTime gpst,gfc::GVector &XSun, gfc::GVector &XSat)
    {
        static double pre_factor = 1.0;
        static bool first  = true;
        double myfactor = 1.0;
        if(first == true)
        {
            myfactor = GMotionState::shadowFactor(XSun, XSat);
            //myfactor = GMotionState::shadowFactor_SECM(true, XSun, XSat);
            pre_factor = myfactor;
            first = false;
            return;
        }
   
        myfactor = GMotionState::shadowFactor(XSun, XSat);
        //myfactor = GMotionState::shadowFactor_SECM(true, XSun, XSat);
        //printf("%f %f\n",myfactor1,myfactor);
        
        //GTime     gpst =  epoch_gps;
        
        double ct = gpst.toDays();
        //CivilTime ct =GTime::GTime2CivilTime(GTime::UTC2GPST(m_epoch));
        
        
        if(first == false)
        {
            //judge the eclipse event according to pre_facor and myfactor
            if(pre_factor == 1.0 && myfactor < 1.0)
            {
                printf("eclipse_start: %16.10f ", ct  );
                //printf("pre:%f -- cur:%f\n",pre_factor,myfactor);
            }
            else if(pre_factor > 0.0 && myfactor == 0.0)
            {
                printf("umbra_s: %16.10f ", ct);
                //printf("pre:%f -- cur:%f\n",pre_factor,myfactor);
            }
            else if(pre_factor == 0.0 && myfactor > 0.0)
            {
                printf("umbra_e: %16.10f ", ct);
                //printf("pre:%f -- cur:%f\n",pre_factor,myfactor);
            }
            else if(pre_factor < 1.0 && myfactor == 1.0)
            {
                printf("eclipse_end: %16.10f %s\n", ct, GTime::GTime2CivilTime(gpst).TimeString().c_str());
                //printf("pre:%f -- cur:%f\n",pre_factor,myfactor);
            }
            
            pre_factor = myfactor;
            
        }
        
        
        
    }
    
    
    double GMotionState::shadowFactor_SECM(bool on, GVector& XSun, GVector& XSat)
    {
        
        double factor = 1.0;
        
        double rad_sun = 695700.0;
        double rad_Earth= 6371.0;
        
        
        double p = 6378.137;
        double q = 6356.7523142;
        double q2 = q*q;
        double p2 = p*p;
        double q_p2 = q2/p2;
        
        
        GVector sat_sun = XSun - XSat;
        
        double xxyy = XSat.x * XSat.x + XSat.y * XSat.y;
        double zz = XSat.z * XSat.z;
        GVector n = crossproduct(XSat, sat_sun);
        
        // If M is a matrix [[q,0,0],[0,q,0],[0,0,p]]
        // A = M * M / p^2
        // B = M / sqrt(p * q)
        double rTAr = xxyy * q_p2 + zz;
        double rTAs = (XSat.x * sat_sun.x + XSat.y * sat_sun.y) * q_p2 + (XSat.z * sat_sun.z);
        double Br_cross_Bs_sqr = (n.x * n.x + n.y * n.y) + (n.z * n.z) * q_p2;
        
        // The discriminant of a quadratic used to solve for the linear combination
        // of rso and sun vectors that give the Earth edge vector.
        double disc = std::sqrt(((xxyy - p2) * q2 + zz * p2) / Br_cross_Bs_sqr);
        
        GVector edge = ((q2 - rTAs * disc) / rTAr - 1.0) * XSat + disc * sat_sun;
        
        GVector rr = edge + XSat;
        
        if(on == true)
        {
            //replace the radius of earth by the corrected rr.norm()
            rad_Earth = rr.norm();
        }
        
        GVector sat_sun_u = normalise(sat_sun);
        GVector xsat_u = normalise(XSat);
        
        double dis_sat_sun = sat_sun.norm();
        double s = XSat.norm();
        
        double a = asin(rad_sun/dis_sat_sun);
        double b = asin(rad_Earth/s);
        double c = acos( -dotproduct(xsat_u, sat_sun_u)  );
        
        if( a+b <= c )
        {
            factor = 1.0;
        }
        else if( c< b-a)
        {
            factor = 0.0;
        }
        else if(  fabs(a-b)<c && c < a+b )
        {
            double x = (c*c + a*a -b*b)/(2.0*c);
            double y = sqrt(a*a -x*x);
            double A = a*a*acos(x/a)+b*b*acos((c-x)/b)-c*y;
            factor = 1.0 - A/(M_PI*a*a);
        }
        
        
        
        return factor;
    }
    
    
    double GMotionState::shadowFunction(GVector& sunpos_ecef, GVector& satposecef)
    {
        
        
        //double mytest = satposecef.norm();
        
       
        
        
        double ratio = 1.0;
        
        double Rs = 695700.0; //km
        double a = 6378.137; //km
        double b = 6356.7523142; //km  6356.7523142
        double ab = a*b;
        double a2 = a*a;
        double b2 = b*b;
        double ab2 = ab*ab;
        
       
        GVector r = satposecef; // the position of satellite
        GVector ru= normalise(r);
        GVector rs = sunpos_ecef; // the position of sun
        
        // step1: budindg the ISF coordinate system
        GVector n = r - rs;
        double dis_sat_sun = n.norm(); // the distance from satellite to sun
        n = normalise(n);
        double f = 1000.0;  //焦距
        double t = dotproduct(n, r);
        GVector u,v;
        if(t == 0.0 )
        {
            u =r;
        }
        else
        {
            u = (1.0/dotproduct(ru,n))*ru - n;
        }
        u = normalise(u);
        
        v = crossproduct(n, u);
        v = normalise(v);
        
        //step2: calculate s1,s1,t1,t2,compute the position of A and B
        double ss = n.x*n.x*a2 + n.y*n.y*a2 + n.z*n.z*b2;
        double s1 = sqrt(1.0/ss);
        double s2 = -s1;
        //double rtn = dotproduct(r, n);
        double t1 = 1.0/s1 - t;
        double t2 = 1.0/s2 - t;
        double disAB[2],disA;
        disAB[0] = (r + t1*n -rs).norm();
        disAB[1] = (r + t2*n -rs).norm();
        disA = disAB[0]<=disAB[1]?disAB[0]:disAB[1];
        if(dis_sat_sun <= disA)
        {
            return 1.0;
        }
        
        //step3. calculate matrix M.
        double M[3][3]={0.0};
        t = (r.x*r.x/a2 + r.y*r.y/a2 + r.z*r.z/b2 - 1.0);
        M[0][0] = (r.x*r.x/a2/a2 - t/a2);
        M[0][1] = (r.x*r.y/a2/a2);
        M[0][2] = (r.x*r.z/a2/b2);
        M[1][0] = M[0][1];
        M[1][1] = (r.y*r.y/a2/a2 - t/a2);
        M[1][2] = (r.y*r.z/a2/b2);
        M[2][0] = M[0][2];
        M[2][1] = M[1][2];
        M[2][2] = (r.z*r.z/b2/b2 - t/b2);
        //calculate k parameters
        double k0,k1,k2,k3,k4,k5;
        k0 =  (u.x*(u.x*M[0][0] + u.y*M[1][0] + u.z*M[2][0])
            + u.y*(u.x*M[0][1] + u.y*M[1][1] + u.z*M[2][1])
            + u.z*(u.x*M[0][2] + u.y*M[1][2] + u.z*M[2][2]));
        
        k1 =  (v.x*(u.x*M[0][0] + u.y*M[1][0] + u.z*M[2][0])
            + v.y*(u.x*M[0][1] + u.y*M[1][1] + u.z*M[2][1])
            + v.z*(u.x*M[0][2] + u.y*M[1][2] + u.z*M[2][2]));
        
        k2 =  (v.x*(v.x*M[0][0] + v.y*M[1][0] + v.z*M[2][0])
            + v.y*(v.x*M[0][1] + v.y*M[1][1] + v.z*M[2][1])
            + v.z*(v.x*M[0][2] + v.y*M[1][2] + v.z*M[2][2]));
        
        k3 = 2.0*f*( n.x*(u.x*M[0][0] + u.y*M[1][0] + u.z*M[2][0])
            + n.y*(u.x*M[0][1] + u.y*M[1][1] + u.z*M[2][1])
            + n.z*(u.x*M[0][2] + u.y*M[1][2] + u.z*M[2][2]));
        
        k4 = 2.0*f*( n.x*(v.x*M[0][0] + v.y*M[1][0] + v.z*M[2][0])
                    + n.y*(v.x*M[0][1] + v.y*M[1][1] + v.z*M[2][1])
                    + n.z*(v.x*M[0][2] + v.y*M[1][2] + v.z*M[2][2]));
        
        k5 =  f*f*(n.x*(n.x*M[0][0] + n.y*M[1][0] + n.z*M[2][0])
                  + n.y*(n.x*M[0][1] + n.y*M[1][1] + n.z*M[2][1])
                  + n.z*(n.x*M[0][2] + n.y*M[1][2] + n.z*M[2][2]));
        
        
        //step3: determine the shape of the projection, ellipse or hyperbola
        // get the eigen values and vector for a symmetric matrix B=[k0,k1;k1,k2]
        // ref: http://courses.cs.washington.edu/courses/cse590b/02wi/eig2x2.cpp
        double lambda[2],evec1[2],evec2[2];
        double disc = sqrt((k0-k2)*(k0-k2)+4*k1*k1)/2.0;
        lambda[0] = (k0+k2)/2.0 + disc;
        lambda[1] = (k0+k2)/2.0 - disc;
        //test for the 0 eigen value
        if (fabs(lambda[0]) < 1.0E-15 || fabs(lambda[1]) < 1.0E-15)
        {
            printf("\twarning: lambda1 = %f, lambda2 = %f\n",lambda[0],lambda[1]);
            exit(0);
        }
        
        evec1[0] = -k1;  evec1[1] = k0-lambda[0];
        evec2[1] = -k1;  evec2[0] = k2-lambda[1];
        
        double v1mag = sqrt(evec1[0]*evec1[0] + evec1[1]*evec1[1]);
        double v2mag = sqrt(evec2[0]*evec2[0] + evec2[1]*evec2[1]);
        evec1[0]/= v1mag;  evec1[1]/=v1mag; // semi-major axis direction
        evec2[0]/= v2mag;  evec2[1]/=v2mag; // semi-minor axis direction
        
         // make sure we have eigenvectors (this may fail for badly-conditioned matrices)
         assert( fabs( (k0*evec1[0] + k1*evec1[1]) - lambda[0]*evec1[0]) < 1e-4);
         assert( fabs( (k1*evec1[0] + k2*evec1[1]) - lambda[0]*evec1[1]) < 1e-4);
         
         assert( fabs( (k0*evec2[0] + k1*evec2[1]) - lambda[1]*evec2[0]) < 1e-4);
         assert( fabs( (k1*evec2[0] + k2*evec2[1]) - lambda[1]*evec2[1]) < 1e-4);
        /*
        //always make the major axis at the direction of evec1
         // with out and with this process, AA and BB are different, confusing!!!
        if(fabs(lambda[0]<fabs(lambda[1])))
        {
            double testvector[2]={0.0}, testvalue=0.0;
            testvalue = lambda[1];
            lambda[1] = lambda[0];
            lambda[0] = testvalue;
            
            testvector[0] = evec2[0]; testvector[1] = evec2[1];
            evec2[0] = evec1[0];evec2[1] = evec1[1];
            evec1[0] = testvector[0];evec1[1] = testvector[1];
        }
        */
        double Y[2]= {0.0};
        Y[0] = evec1[0]*k3 + evec2[0]*k4;
        Y[1] = evec1[1]*k3 + evec2[1]*k4;
        
        //calculate phi, this conical equation is rotated by the eigen vectors.
        double phi = (  Y[0]*Y[0]/4.0/lambda[0] + Y[1]*Y[1]/4.0/lambda[1] - k5 )/lambda[0]/lambda[1];
        double AA = lambda[1]*phi; // square of the semi-major axis bind with evec2
        double BB = lambda[0]*phi; // square of the semi-minor axis bind with evec1
        double CC[2] ={ -Y[0]/2.0/lambda[0], -Y[1]/2.0/lambda[1] }; // center of the ellipse X and Y
        
        //step4:solve the quartic equation about k paramters and the apparent sun
        double XX[4]={0.0};
        double R0 = f*Rs/dis_sat_sun; // the radius of apparent sun
        double A=BB*(CC[0]+R0)*(CC[0]+R0)+AA*CC[1]*CC[1]-AA*BB;
        double B=-4*AA*R0*CC[1];
        double C=2.0*AA*CC[1]*CC[1]+2.0*BB*CC[0]*CC[0]-2.0*BB*R0*R0+4.0*AA*R0*R0-2.0*AA*BB;
        double D=-4*AA*R0*CC[1];
        double E=BB*(R0-CC[0])*(R0-CC[0])+AA*CC[1]*CC[1]-AA*BB;
        double coef[4] = {B/A,C/A,D/A,E/A};
        
        
        int num_of_solution = GMath::quarticSolver(coef,XX);
        
        
        if( num_of_solution == 3 || num_of_solution == 4)
        {
            printf("WARNING: shadowfactor, %d intersections!", num_of_solution);
            int testc = 0;
        }
        
        // the intersections, maximum number of intersections is 4
        double intersections[4][2] = {0.0};
        for(int i = 0 ; i< num_of_solution; i++ )
        {
            intersections[i][0] = (1-XX[i]*XX[i])*R0/(1+XX[i]*XX[i]);
            intersections[i][1] = 2*XX[i]*R0/(1+XX[i]*XX[i]);
        }
        
        double abc[3]={0.0};
        bool test = false;
        
        if(num_of_solution >= 2)
        {
            // the line between the two intersections
            abc[0]=intersections[1][1] - intersections[0][1];
            abc[1]=intersections[0][0] - intersections[1][0];
            abc[2]=intersections[1][0]*intersections[0][1] - intersections[0][0]*intersections[1][1];
            //test if the centre of circle and the centre of conical curve are in the same side of this line or not
            double test1 = abc[0]*CC[0] + abc[1]*CC[1] + abc[2];
            test = test1*abc[2] >0? true : false;
        }
        else  // no intersection, test the centre of the sun is inside the ellipsoid or not
        {
            //test if the centre of sun is inside the conical projection
            double test2 =  CC[0]*CC[0]/AA + CC[1]*CC[1]/BB;
            test = test2 < 1.0 ? true: false;
            
        }
        
        if( lambda[0]*lambda[1]<0.0 ) // 非椭圆
        {
            if( test == true ) // the centre of sun is inside the hyperbola
            {
                if( num_of_solution <2 )
                {
                    ratio = 1.0;   //full phase
                }
                else  // penumbra
                {
                    //printf("Hyperbola on the photo_1!\n");
                    ratio = 0.5;
                    
                    int testc =0;
                }
            }
            else // the centre of sun is ouside of hyperbola
            {
                if( num_of_solution <2 )
                {
                    ratio = 0.0;  // umbra
                }
                else  // penumbra
                {
                    // printf("Hyperbola on the photo_2!\n");
                    ratio = 0.5;
                    int testc = 0;
                }
            }
            
        }
        else  //椭圆
        {
            
            //the angle between the two intersections for the ellipse
            double v1[2]={intersections[0][0]-CC[0], intersections[0][1]-CC[1]};
            double v2[2]={intersections[1][0]-CC[0], intersections[1][1]-CC[1]};
            double len1_e = sqrt(v1[0]*v1[0]+v1[1]*v1[1]);
            double len2_e = sqrt(v2[0]*v2[0]+v2[1]*v2[1]);
            double theta_ellipse = acos((v1[0]*v2[0]+v1[1]*v2[1])/(len1_e*len2_e));
            double S_tri_e = 0.5*fabs(v1[0]*v2[1]-v1[1]*v2[0]);
            
            // the angle between the two intersections for the circle
            v1[0] = intersections[0][0];v1[1] = intersections[0][1];
            v2[0] = intersections[1][0];v2[1] = intersections[1][1];
            double len1_c =sqrt(v1[0]*v1[0]+v1[1]*v1[1]);
            double len2_c =sqrt(v2[0]*v2[0]+v2[1]*v2[1]);
            double theta_circle = acos((v1[0]*v2[0]+v1[1]*v2[1])/(len1_c*len2_c)  );
            double S_tri_c = 0.5*fabs(v1[0]*v2[1]-v1[1]*v2[0]);
            
            //ref: Table1 in "Calculating ellipse overlap areas" from Gary B.Hughes and Mochine Chraibi
            // the angle theta is the angular parameter for the ellipse, the intersections has to be in counter clockwise order
            double theta[4]={0.0};
            for( int i = 0 ; i< num_of_solution; i++ )
            {
                double x = intersections[i][0] - CC[0];
                double y = intersections[i][1] - CC[1];
                if( (x <0 && y>=0) || (x >=0 && y>=0) )
                {
                    theta[i] = acos(x/sqrt(AA));
                }
                else if( (x<0 && y<0)||(x>=0 && y<0))
                {
                    theta[i] = 2.0*M_PI - asin(x/sqrt(AA));
                }
            }
            
            if( test == true ) // the centre of sun is inside the ellipse, same side
            {
                if( num_of_solution <2 )
                {
                    ratio = 0.0;   //umbra
                }
                else  // penumbra
                {
                    double S_sector_c = 0.5*R0*R0*theta_circle;
                    double S_sector_e = 0.5*sqrt(AA*BB)*fabs(theta[1]- theta[0]);
                    double S_shadow = (S_sector_c - S_tri_c) - (S_sector_e - S_tri_e);
                    ratio = S_shadow/( M_PI*R0*R0);
                }
            }
            else // the centre of sun is ouside of ellipse, different side
            {
                if( num_of_solution <2 )
                {
                    ratio = 1.0;  // full phase
                }
                else  // penumbra
                {
                    double S_sector_c = 0.5*R0*R0*theta_circle;
                    double S_sector_e = 0.5*sqrt(AA*BB)*fabs(theta[1]- theta[0]);
                    double S_shadow = S_sector_c - S_tri_c + S_sector_e - S_tri_e;
                    ratio = 1.0 - S_shadow/(M_PI*R0*R0);
                }
            }
            
            
        }
        
        if(isnan(ratio))
        {
            int testc = 0;
        }
        
        return ratio;
    }
    
    
    // the shadow factor function to decide the scale of flux in eclipse
    /*
     * the basic idea here is put a camera on the satellite and take a photo of earth and sun,
     * and then all the calculation is done on this photo
     overlap area of ellipses
     ref:  https://arxiv.org/pdf/1106.3787.pdf
     
     // the area of right hyperbolic sector
     http://keisan.casio.com/exec/system/14770153128432
     
     */
    double GMotionState::shadowFactor(GVector& sunpos_eci, GVector& satpos_eci)
    {
        
        //the parameters for the earth ellipsoid
        double Rs = 695700.0; //km
        double a = 6378.137; //km
        double b = 6356.7523142; //km  6356.7523142
        double ab = a*b;
        double a2 = a*a;
        double b2 = b*b;
        double ab2 = ab*ab;
        
        double dis_sat_earth = satpos_eci.norm();
        
        double factor = 1.0;
        
        GVector s; // the origin of the photo coordinate system
        // these are the external orientation elements
        GVector u,v,n;
        
        n = satpos_eci - sunpos_eci;
        n.normalise();
        
        double dis_sat_sun = (satpos_eci - sunpos_eci).norm();
        
        
        // satellite and sun colinear
        
        
        
        // first test if the satellite is in the front of earth,
        // if it is in the front of earth, it is always full phase
        // otherwise,using the photogrammetry method
        int state = -1;
        
        double ss = n.x*n.x*a2 + n.y*n.y*a2 + n.z*n.z*b2;
        double s1= sqrt(1.0/ss);
        double s2= -sqrt(1.0/ss);
        
        double t1 = (1.0/s1 - dotproduct(satpos_eci, n))/dotproduct(n, n);
        double t2 = (1.0/s2 - dotproduct(satpos_eci, n))/dotproduct(n, n);
        
        GVector xA = satpos_eci + t1*n;
        GVector xB = satpos_eci + t2*n;
        
        double disA = (xA-sunpos_eci).norm();
        double disB = (xB-sunpos_eci).norm();
        double disMin = disA <= disB? disA: disB;
        
        double disMax = disA >= disB? disA: disB;
        
        
        //calculate the position of E1 and E2
        GVector E1,E2;
        E1.x = s1*a2*n.x;
        E1.y = s1*a2*n.y;
        E1.z = s1*b2*n.z;
        
        E2.x = s2*a2*n.x;
        E2.y = s2*a2*n.y;
        E2.z = s2*b2*n.z;
        
        GVector m = (xA-E1);
        m.normalise();
        
        ss = m.x*m.x*a2+m.y*m.y*a2+m.z*m.z*a2;
        s1 = sqrt(1.0/ss);
        s2 = -sqrt(1.0/ss);
        double d1 = 1.0/s1 - dotproduct(xA, m);
        double d2 = 1.0/s2 - dotproduct(xA, m);
        GVector xx; // the coordinate of the tangent point on the earth
        ss = fabs(d1)<=fabs(d2)?s1:s2;
        xx.x = ss*a2*m.x;
        xx.y = ss*a2*m.y;
        xx.z = ss*b2*m.z;
        
        double param = dotproduct(satpos_eci+xx,m);
        GVector xP = xx + param*m;
        
        double dis1 = (xP - sunpos_eci).norm();
        
        //printf("%.6f %.6f %.6f %.6f ", dis1, dis_sat_sun, disA, disB);
        
        if( dis_sat_sun <= dis1 ) // full phase
        {
            state = 0;
            return 1.0;
        }
        else if( dis_sat_sun >= dis1 && dis_sat_sun <= disMax ) // what's in the photo is hyperbola
        {
            state = 1;  // hyperbola
        }
        else if(dis_sat_sun >= disMax) // what's in the photo is ellipse
        {
            state = 2;
        }
        
        
        
        GVector r = satpos_eci ;//+ (3.5*dis_sat_earth)*n; // the coordinate of view point, it will affect the penumbra
        
        // this is the so called internal orientation elements, x0=y0=0.0 in this case
        //*****************************************************************
        double f = 1000.0; // the focal distance 1m = 0.001 km, changing f will change the solution situation which is wired !!!!
        
        s = r + f*n; // the origin of the photo coordinate system
        
        /*
         the definiton of photo coordinate system:
         the origin is projection of center of sun
         the most situation:
         u is from origin to the projection of center of earth
         n is from sun to satellite
         v form a right-hand system
         */
        
        // determine the photo-x axis vector
        //calculat the projection of eci_Z on the fundamental plane
        //u = eci_Z - dotproduct(eci_Z, n)/n.norm()/n.norm()*n;
       //u = n.norm()*n.norm()/dotproduct(eci_Z, n)*eci_Z - n;
       //u.normalise();
        
        GVector d = r;
        d.normalise();
       
        double t = 0.0;
        t = f/dotproduct(d, n);
        GVector Xorigin = r + t*d; // the coordinate of the projection of earth center
        
        u = Xorigin - s;
        u.normalise();
        
        v = crossproduct(n, u);
        v.normalise();
        
        GMatrix M(3,3);//A(3,3);
        //A(0,0) = 1.0/a2; A(1,1)=1.0/a2; A(2,2)= 1.0/b2;
        
        t = (r.x*r.x/a/a + r.y*r.y/a/a + r.z*r.z/b/b - 1.0);
        
        M(0,0) = r.x*r.x/a2/a2 - t/a2;
        M(0,1) = r.x*r.y/a2/a2;
        M(0,2) = r.x*r.z/a2/b2;
        M(1,0) = M(0,1);
        M(1,1) = r.y*r.y/a2/a2 - t/a2;
        M(1,2) = r.y*r.z/a2/b2;
        M(2,0) = M(0,2);
        M(2,1) = M(1,2);
        M(2,2) = r.z*r.z/b2/b2 - t/b2;
        
        M = M*4.0;
        
        //cout<<M<<endl;
        
        //calculate K
        double K[6] = {0.0};
        GMatrix tu(3,1), tv(3,1), ts(3,1),tr(3,1),tn(3,1);
        tu[0] = u.x;tu[1] = u.y;tu[2] = u.z;
        tv[0] = v.x;tv[1] = v.y;tv[2] = v.z;
        ts[0] = s.x;ts[1] = s.y;ts[2] = s.z;
        tr[0] = r.x;tr[1] = r.y;tr[2] = r.z;
        tn[0] = n.x;tn[1] = n.y;tn[2] = n.z;
        
        ((~tu)*M*tu).getData(K); // K[0]
        ((~tu)*M*tv).getData(K+1); // K[1]
        ((~tv)*M*tv).getData(K+2); // K[2]
        (2.0*f*(~tu*M*tn)).getData(K+3); //K[3]
        (2.0*f*(~tv*M*tn)).getData(K+4);  //K[4]
        (f*f*(~tn*M*tn)).getData(K+5); //K[5]
        
        
        for(int i =0; i< 6; i++ )
        {
            K[i] = K[i]*1.0E6;
        }
        
        if( K[0]<0.0 && K[2] <0.0 )  // B is negative definite,  make sure D is positive definite
        {
            for( int i=0; i< 6; i++ )
            {
                K[i] = -K[i];
            }
        }
        
        double scale = (r - sunpos_eci).norm(); // the distance from sun to viewpoint
        // the radius of sun on the photo, it is a circle
        double rs = f/scale*Rs;
        
        
        double XX[4]={0.0}; // the solution of quartic equation
        int num_of_solution = 0;
        
        double Kcoef[5] = {0.0};
        Kcoef[0] = K[0]*rs*rs - K[3]*rs + K[5];
        Kcoef[1] = 2*K[4]*rs - 4*K[1]*rs*rs;
        Kcoef[2] = 4.0*K[2]*rs*rs - 2.0*K[0]*rs*rs + 2.0*K[5];
        Kcoef[3] = 4.0*K[1]*rs*rs + 2.0*K[4]*rs;
        Kcoef[4] = K[0]*rs*rs + K[3]*rs + K[5];
        double Kce[4] = {Kcoef[1]/Kcoef[0],Kcoef[2]/Kcoef[0],Kcoef[3]/Kcoef[0],Kcoef[4]/Kcoef[0] };
        
        num_of_solution = GMath::quarticSolver(Kce,XX);
        if( num_of_solution == 3 || num_of_solution == 4)
        {
            printf("WARNING: shadowfactor, %d intersections!", num_of_solution);
            int testc = 0;
        }
        
        //calculate the eigen value and eigen vector of matrix B
        double tt =  (K[0]+K[2])*(K[0]+K[2]) - 4.0*(K[0]*K[2] - K[1]*K[1]) ;
        int test2 = tt>0?1:-1;
        tt = sqrt(tt);
        double lambda1 = (K[0]+K[2]+tt)/2.0;
        double lambda2 = (K[0]+K[2]-tt)/2.0;
        
        //get the eigen vector
        double  r1[2]={0.0,1.0},r2[2]={1.0,0.0};
        bool sol1 = false, sol2 =false;
        if( fabs(lambda1 - K[0])<1.0E-12)
        {
            r1[1] = 0.0;
            if( (lambda1 - K[0])*K[1] <0.0 ) // 异号
            {
                r1[0] = -1.0;
            }
            else
            {
                r1[0] = 1.0;
            }
            sol1 = true;
        }
        
        if( fabs(lambda2 - K[2])<1.0E-12)
        {
            r2[0] = 0.0;
            if( (lambda2 - K[2])*K[1] <0.0 ) // 异号
            {
                r2[1] = -1.0;
            }
            else
            {
                r2[1] = 1.0;
            }
            sol2 = true;
        }
        
        if(sol1 ==false)
        {
            r1[0] = K[1]/(lambda1 - K[0]);  // r1 for lambda1
        }
        
        if(sol2 ==false)
        {
            r2[1] = K[1]/(lambda2 - K[2]);  // r2 for lambda2
        }
       
        
        
        //get the unit vector
        double tmp = sqrt(r1[0]*r1[0]+r1[1]*r1[1]);
        r1[0] = r1[0]/tmp;r1[1] = r1[1]/tmp;
        
        tmp = sqrt(r2[0]*r2[0]+r2[1]*r2[1]);
        r2[0] = r2[0]/tmp;r2[1] = r2[1]/tmp;
        
        //the larger eigen value is lambda1
        if( lambda1 <= lambda2 )  // swap lamda1 and lamda2, with r1 and r2
        {
            double t = lambda1;
            lambda1 = lambda2;
            lambda2 = t;
            double r[2] ={r1[0],r1[1]};
            r1[0] = r2[0]; r1[1] = r2[1];
            r2[0] = r[0]; r2[1] = r[1];
        }
        
        double Y[2]= {0.0};
        Y[0] = r1[0]*K[3] + r2[0]*K[4];
        Y[1] = r1[1]*K[3] + r2[1]*K[4];
        
        
        
        //calculate phi
        double phi1 = Y[0]*Y[0]/4.0/lambda1 + Y[1]*Y[1]/4.0/lambda2 - K[5];
        
        double phi = (  Y[0]*Y[0]/4.0/lambda1 + Y[1]*Y[1]/4.0/lambda2 - K[5] )/lambda1/lambda2;
        
        double CC[2] ={ -Y[0]/2.0/lambda1, -Y[1]/2.0/lambda2 }; // center X and Y
        
        //calculate the semi-major axis and semi-minor axis
        // since the curve can be a hyperbolor, thus
        double AA = lambda2*phi; // x axis, a^2 , including sign
        double BB = lambda1*phi; // y axis, b^2, including sign
        
        
        //then solve the quartic equation
        // ref: http://users.nik.uni-obuda.hu/sanyo/gpgpu/sami2015_submission_60.pdf
        // using Ferrari method
        
        //first, forming the intersection quartic equation
        t = BB*(rs+CC[0])*(rs+CC[0]) + AA*CC[1]*CC[1] - AA*BB; // A
        double A = 1.0;
        double B = (-4.0*AA*rs*CC[1])/t;
        double C = (2.0*AA*CC[1]*CC[1] + 2.0*BB*CC[0]*CC[0] + 4.0*AA*rs*rs - 2.0*BB*rs*rs  - 2.0*AA*BB)/t;
        double D = (-4.0*AA*rs*CC[1])/t;
        double E = (BB*(rs-CC[0])*(rs-CC[0]) + AA*CC[1]*CC[1] - AA*BB)/t;
        
        // A=-15;B=0;C=376;D=0;E=-28; // for the set of x^2+y^2=9 and x^2/16+y^2/4=1
        //A = -1; B=0;C=0;D=0;E=1;
        // just need to find out all the real solutions
        
        double ce[4]={B,C,D,E};
        
        num_of_solution = GMath::quarticSolver(ce,XX);
        
        if( num_of_solution == 3 || num_of_solution == 4)
        {
            printf("WARNING: shadowfactor, %d intersections!\n", num_of_solution);
            exit(0);
        }

        //cout<<"num_of_solution:"<< num_of_solution<<endl;
        //printf("%2d %8.6f %8.6f   %8.6f %8.6f   %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f   %12.8f %12.8f %12.8f\n",
        //       num_of_solution,AA,BB,CC[0],CC[1],K[0],K[1],K[2],K[3],K[4],K[5], lambda1,lambda2,phi);
        
        //then testing, whether the center of sun is inside the earth, this is for penumbra
        
        
        if( AA*BB <=0.0 )  // hyperbola
        {
             double test =  CC[0]*CC[0]/AA + CC[1]*CC[1]/BB;
             if( test <= 1.0 )
             {
                 if( num_of_solution <2 )
                 {
                     factor = 1.0;   //full phase
                 }
                 else  // penumbra
                 {
                     //printf("Hyperbola on the photo_1!\n");
                     factor = 0.5;
                     int testc =0;
                 }
             }
            else if( test >= 1.0)
            {
                if( num_of_solution <2 )
                {
                    factor = 0.0;  // umbra
                }
                else  // penumbra
                {
                    //printf("Hyperbola on the photo_2!\n");
                    factor = 0.5;
                    int testc = 0;
                }
            }
        }
        else if( AA*BB >=0.0 )  // ellipse
        {
            double test =  CC[0]*CC[0]/AA + CC[1]*CC[1]/BB;
            
            if( test < 1.0  )  // center of sun is inside earth, different strategy for calculation or area ratio
            {
                if( num_of_solution < 2 )  //umbra
                {
                    factor = 0.0;
                }
                else //penumbra
                {
                    //calculating the area ratio
                    
                    //it seems it is impossible to get solution more than 2
                    double xc1[2]={0,0},xc2[2]={0.0}; // the coordinate of C and C'
                    xc1[0] = (1.0-XX[0]*XX[0])/(1.0+XX[0]*XX[0])*rs;
                    xc1[1] = (2.0*XX[0])/(1.0+XX[0]*XX[0])*rs;
                    
                    xc2[0] = (1.0-XX[1]*XX[1])/(1.0+XX[1]*XX[1])*rs;
                    xc2[1] = (2.0*XX[1])/(1.0+XX[1]*XX[1])*rs;
                    
                    // the area of triangle TACC' , 2d cross product
                    double TACC = 0.5*fabs( xc1[0]*xc2[1] - xc2[0]*xc1[1]);
                    double TBCC = 0.5*fabs( (xc1[0]-CC[0])*(xc2[1]-CC[1]) - (xc2[0]-CC[0])*(xc1[1]-CC[1]) );
                    
                    // the area of the circular part
                    double t1 = sqrt(xc1[0]*xc1[0] + xc1[1]*xc1[1] );
                    double t2 = sqrt(xc2[0]*xc2[0] + xc2[1]*xc2[1] );
                    double SectorACC = 0.5*acos((xc1[0]*xc2[0] + xc1[1]*xc2[1])/t1/t2)*rs*rs;
                    
                    double SectorBCC = 0.0;
                    double d[2]={0.0},e[2]={0.0};
                    //calculating the area of elliptic sector
                    if( AA >= BB )  // semi-major axis is x axis
                    {
                        //以半径为R画个圆
                        double R = sqrt(AA);
                        double y1 = CC[1]+sqrt(R*R - (xc1[0]-CC[0])*(xc1[0]-CC[0]));
                        double y2 = CC[1]-sqrt(R*R - (xc1[0]-CC[0])*(xc1[0]-CC[0]));
                        d[0] = xc1[0];
                        // get the nearest point
                        d[1] = fabs(y1 - xc1[1]) <= fabs(y2 - xc1[1]) ? y1 : y2;
                        
                        y1 = CC[1]+sqrt(R*R - (xc2[0]-CC[0])*(xc2[0]-CC[0]));
                        y2 = CC[1]-sqrt(R*R - (xc2[0]-CC[0])*(xc2[0]-CC[0]));
                        
                        e[0] = xc2[0];
                        // get the nearest point
                        e[1] = fabs(y1 - xc2[1]) <= fabs(y2 - xc2[1]) ? y1 : y2;
                        
                        //double xc1[2] = {x1[0],0.0},xc2[2]={x2[0],0.0};
                        //angle CBC', in fact, it should be Angle DBE
                        double tmp1[2] = {d[0]-CC[0],d[1]-CC[1]};
                        double tmp2[2] = {e[0]-CC[0],e[1]-CC[1]};
                        t1 = sqrt(tmp1[0]*tmp1[0]+tmp1[1]*tmp1[1]);
                        t2 = sqrt(tmp2[0]*tmp2[0]+tmp2[1]*tmp2[1]);
                        double Acbc = (tmp1[0]*tmp2[0]+tmp1[1]*tmp2[1])/t1/t2;
                        if(Acbc >1.0) Acbc =1.0;
                        
                        if(Acbc >1.0) Acbc =1.0;
                        else if(Acbc<-1.0) Acbc = -1.0;
                        
                        SectorBCC = 0.5*Acbc*sqrt(AA*BB);
                        
                    }
                    else if(AA <= BB) // semi-major axis is y axis
                    {
                        double R = sqrt(BB);
                        
                        double x1 = CC[0]+sqrt( R*R - (xc1[1]-CC[1])*(xc1[1]-CC[1]));
                        double x2 = CC[0]-sqrt( R*R - (xc1[1]-CC[1])*(xc1[1]-CC[1]));
                        d[1] = xc1[1];
                        // get the nearest point
                        d[0] = fabs( x1 - xc1[0]) <= fabs(x2 - xc1[0]) ? x1 : x2;
                        
                        x1 = CC[0]+sqrt( R*R - (xc2[1]-CC[1])*(xc2[1]-CC[1]));
                        x2 = CC[0]-sqrt( R*R - (xc2[1]-CC[1])*(xc2[1]-CC[1]));
                        e[1] = xc2[1];
                        // get the nearest point
                        e[0] = fabs(x1 - xc2[0]) <= fabs(x2 - xc2[0]) ? x1 : x2;
                        
                        //double xc1[2] = {x1[0],0.0},xc2[2]={x2[0],0.0};
                        //angle CBC', in fact, it should be Angle DBE
                        double tmp1[2] = {d[0]-CC[0],d[1]-CC[1]};
                        double tmp2[2] = {e[0]-CC[0],e[1]-CC[1]};
                        t1 = sqrt(tmp1[0]*tmp1[0]+tmp1[1]*tmp1[1]);
                        t2 = sqrt(tmp2[0]*tmp2[0]+tmp2[1]*tmp2[1]);
                        double Acbc = (tmp1[0]*tmp2[0]+tmp1[1]*tmp2[1])/t1/t2;
                        
                        if(Acbc >1.0) Acbc =1.0;
                        else if(Acbc<-1.0) Acbc = -1.0;
                        
                        Acbc = acos(Acbc);
                        SectorBCC = 0.5*Acbc*sqrt(AA*BB);
                        
                    }
                    
                    
                    
                    double S = (SectorACC - TACC) - (SectorBCC - TBCC);
                    
                    factor =  fabs(S/(3.14159265357*rs*rs));
                    
                    //factor = 0.5;
                    int testc = 0;
                    
                }
            }
            else // the center of sun is outside of the earth
            {
                if( num_of_solution < 2 )  //full phase
                {
                    
                    factor = 1.0;
                    
                }
                else  // penumbra
                {
                    //calculate the area ratio
                    
                    //it seems it is impossible to get solution more than 2
                    double xc1[2]={0,0},xc2[2]={0.0}; // the coordinate of C and C'
                    xc1[0] = (1.0-XX[0]*XX[0])/(1.0+XX[0]*XX[0])*rs;
                    xc1[1] = (2.0*XX[0])/(1.0+XX[0]*XX[0])*rs;
                    
                    xc2[0] = (1.0-XX[1]*XX[1])/(1.0+XX[1]*XX[1])*rs;
                    xc2[1] = (2.0*XX[1])/(1.0+XX[1]*XX[1])*rs;
                    
                    
                    // the area of triangle TACC' , 2d cross product
                    
                    double TACC = 0.5*fabs( xc1[0]*xc2[1] - xc2[0]*xc1[1]);
                    double TBCC = 0.5*fabs( (xc1[0]-CC[0])*(xc2[1]-CC[1]) - (xc2[0]-CC[0])*(xc1[1]-CC[1]) );
                    
                    // the area of the circular part
                    double t1 = sqrt(xc1[0]*xc1[0] + xc1[1]*xc1[1] );
                    double t2 = sqrt(xc2[0]*xc2[0] + xc2[1]*xc2[1] );
                    double SectorACC = 0.5*acos((xc1[0]*xc2[0] + xc1[1]*xc2[1])/rs/rs)*rs*rs;
                    
                    double SectorBCC = 0.0;
                    double d[2]={0.0},e[2]={0.0};
                    
                    //calculating the area of elliptic sector
                    if( AA >= BB )  // semi-major axis is x axis
                    {
                        double R = sqrt(AA);
                        
                        double y1 = CC[1]+sqrt(R*R - (xc1[0]-CC[0])*(xc1[0]-CC[0]));
                        double y2 = CC[1]-sqrt(R*R - (xc1[0]-CC[0])*(xc1[0]-CC[0]));
                        d[0] = xc1[0];
                        // get the nearest point
                        d[1] = fabs(y1 - xc1[1]) <= fabs(y2 - xc1[1]) ? y1 : y2;
                        
                        y1 = CC[1]+sqrt(R*R - (xc2[0]-CC[0])*(xc2[0]-CC[0]));
                        y2 = CC[1]-sqrt(R*R - (xc2[0]-CC[0])*(xc2[0]-CC[0]));
                        e[0] = xc2[0];
                        // get the nearest point
                        e[1] = fabs(y1 - xc2[1]) <= fabs(y2 - xc2[1]) ? y1 : y2;
                        
                        //double xc1[2] = {x1[0],0.0},xc2[2]={x2[0],0.0};
                        //angle CBC', in fact, it should be Angle DBE
                        double tmp1[2] = {d[0]-CC[0],d[1]-CC[1]};
                        double tmp2[2] = {e[0]-CC[0],e[1]-CC[1]};
                        t1 = sqrt(tmp1[0]*tmp1[0]+tmp1[1]*tmp1[1]);
                        t2 = sqrt(tmp2[0]*tmp2[0]+tmp2[1]*tmp2[1]);
                        double Acbc = (tmp1[0]*tmp2[0]+tmp1[1]*tmp2[1])/t1/t2;
                        
                        if(Acbc >1.0) Acbc =1.0;
                        else if(Acbc<-1.0) Acbc = -1.0;
                        
                        Acbc = acos(Acbc);
                        SectorBCC = 0.5*Acbc*sqrt(AA*BB);
                        
                        
                    }
                    else if(AA <= BB) // semi-major axis is y axis
                    {
                        double R = sqrt(BB);
                        
                        double x1 = CC[0]+sqrt( R*R - (xc1[1]-CC[1])*(xc1[1]-CC[1]));
                        double x2 = CC[0]-sqrt( R*R - (xc1[1]-CC[1])*(xc1[1]-CC[1]));
                        d[1] = xc1[1];
                        // get the nearest point
                        d[0] = fabs(x1 - xc1[0]) <= fabs(x2 - xc1[0]) ? x1 : x2;
                        
                        x1 = CC[0]+sqrt( R*R - (xc2[1]-CC[1])*(xc2[1]-CC[1]));
                        x2 = CC[0]-sqrt( R*R - (xc2[1]-CC[1])*(xc2[1]-CC[1]));
                        e[1] = xc2[1];
                        // get the nearest point
                        e[0] = fabs(x1 - xc2[0]) <= fabs(x2 - xc2[0]) ? x1 : x2;
                        
                        //double xc1[2] = {x1[0],0.0},xc2[2]={x2[0],0.0};
                        //angle CBC', in fact, it should be Angle DBE
                        double tmp1[2] = {d[0]-CC[0],d[1]-CC[1]};
                        double tmp2[2] = {e[0]-CC[0],e[1]-CC[1]};
                        t1 = sqrt(tmp1[0]*tmp1[0]+tmp1[1]*tmp1[1]);
                        t2 = sqrt(tmp2[0]*tmp2[0]+tmp2[1]*tmp2[1]);
                        double Acbc = (tmp1[0]*tmp2[0]+tmp1[1]*tmp2[1])/t1/t2;
                        
                        if(Acbc >1.0) Acbc =1.0;
                        else if(Acbc<-1.0) Acbc = -1.0;
                        
                        Acbc = acos(Acbc);
                        SectorBCC = 0.5*Acbc*sqrt(AA*BB);
                        
                    }
                    
                    double S = (SectorACC - TACC) + (SectorBCC - TBCC);
                    
                    factor = 1.0 - fabs(S/(3.14159265357*rs*rs));
                    
                    //factor = 0.5;
                    
                    int testc = 0;
                    
                }
                
            }
            
        }
        
        return factor;
        
    }
    
    
    
    
    
} // end of namespace

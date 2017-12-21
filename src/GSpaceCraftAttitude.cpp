//
//  GSpaceCraftAttitude.cpp
//  GFC
//
//  Created by lizhen on 16/05/2016.
//  Copyright © 2016 lizhen. All rights reserved.
//

#include "GSpaceCraftAttitude.hpp"
namespace gfc
{
    

    
    bool GSpaceCraftAttitude::initialiseAttitude(double shadowFactor,
                                                 GVector& satpos_eci,
                                                 GVector& sunpos_eci,
                                                 GVector& satsun_eci,
                                                 GVector& satvel_eci,
                                                 double yaw_angle, double pitch_angle)
    {
        bool state = false;
        if( fabs(shadowFactor-1.0)<1.0E-8 ) //full phase, the satellite can see the sun
        {
            state = true;
            
            //conduct the nominal attitude
            NorminalYawSteering_gal(satpos_eci, sunpos_eci, satsun_eci, satvel_eci);
            
        }
        else  // satellite is in eclipse and can not find the sun, need to yaw and pitch information
        {
            
            // the reason is they have different definition of Body-Fixed-Frame
            GString satsys = "sysGAL";
            
            if(satsys == "sysGPS")
            {
                
            }
            else if( satsys == "sysGAL")
            {
                //conduct the nominal attitude
                NorminalYawSteering_gal(satpos_eci, sunpos_eci, satsun_eci, satvel_eci);
            }
            else if( satsys == "sysBDS")
            {
                
            }
            
        }
        
         m_isInitialised = true; // it is true whether it is initialised properly or not
        
         return state;
    }
    
    
    
    void GSpaceCraftAttitude::convert2ECEF(GEarthOrientationParameter& eop)
    {
        GVector t;
        eop.ECI2ECEF_pos(h, t);
        h = t;
        eop.ECI2ECEF_pos(n, t);
        n = t;
        eop.ECI2ECEF_pos(xhat, t);
        xhat = t;
        eop.ECI2ECEF_pos(yhat, t);
        yhat = t;
        eop.ECI2ECEF_pos(zhat, t);
        zhat = t;
        eop.ECI2ECEF_pos(phat, t);
        phat = t;
        eop.ECI2ECEF_pos(eR, t);
        eR = t;
        eop.ECI2ECEF_pos(eT, t);
        eT = t;
        eop.ECI2ECEF_pos(eN, t);
        eN = t;
        eop.ECI2ECEF_pos(eD, t);
        eD = t;
        eop.ECI2ECEF_pos(eY, t);
        eY = t;
        eop.ECI2ECEF_pos(eB, t);
        eB = t;
        
        
    }
    
    
    /*
     u is the orbit angle measured from midnight
     beta is the sun elevation measured from orbital plane
     根据专利上的公式，主要是轨道角的定义起点不同，代码中从midnight起算，而专利中从noon起算，相隔180
     */
    void GSpaceCraftAttitude::attitudeLaw_Galileo(double eta, double beta)
    {
        //firstly, use norminal attitude control law
        yaw_model = atan2( tan(beta), -sin(eta));
        
        double cosP = cos(beta)*cos(eta);
        
        pitch_model = atan2( -sqrt(1.0-cosP*cosP),cosP);
        
    }
    
    
    void GSpaceCraftAttitude::setFromYaw_Pitch(double yaw, double pitch)
    {
        
        //旋转矩阵的定义是：与右手螺旋的方向相反
        //1, 将T绕z旋转yaw得到x
        //2, y= z X x
        //3, 将z绕y旋转pitch得到太阳帆板法向p,因为pitch是与z轴正向之间的夹角，pitch=0时，phat与zhat相同
        double cosY= cos(yaw);
        double sinY= sin(yaw);
        double cosP= cos(pitch);
        double sinP= sin(pitch);
        
        // ref: https://en.wikipedia.org/wiki/Rotation_matrix
        // rotation matrix from axis,should be replaced with quarternions for the sake of efficiency
        
        xhat.x = eT.x*( cosY + zhat.x*zhat.x*(1.0-cosY) ) + eT.y*( zhat.x*zhat.y*(1.0-cosY) - zhat.z*sinY ) + eT.z*( zhat.x*zhat.z*(1.0-cosY) + zhat.y*sinY );
        
        xhat.y = eT.x*(zhat.y*zhat.x*(1.0-cosY) + zhat.z*sinY ) + eT.y*(cosY+zhat.y*zhat.y*(1.0-cosY) ) + eT.z*(zhat.y*zhat.z*(1.0-cosY)-zhat.x*sinY );
        
        xhat.z = eT.x*( zhat.z*zhat.x*(1.0-cosY)-zhat.y*sinY ) + eT.y*( zhat.z*zhat.y*(1.0-cosY) + zhat.x*sinY) + eT.z*(cosY+zhat.z*zhat.z*(1.0-cosY) );
        
        yhat = crossproduct(zhat, xhat);
        
        xhat.normalise();
        yhat.normalise();
        
        
        phat.x = zhat.x*( cosP + yhat.x*yhat.x*(1.0-cosP) ) + zhat.y*( yhat.x*yhat.y*(1.0-cosP) - yhat.z*sinP ) + zhat.z*( yhat.x*yhat.z*(1.0-cosP) + yhat.y*sinP );
        
        phat.y = zhat.x*(yhat.y*yhat.x*(1.0-cosP) + yhat.z*sinP ) + zhat.y*(cosP+yhat.y*yhat.y*(1.0-cosP) ) + zhat.z*(yhat.y*yhat.z*(1.0-cosP)-yhat.x*sinP );
        
        phat.z = zhat.x*( yhat.z*yhat.x*(1.0-cosP)-yhat.y*sinP ) + zhat.y*( yhat.z*yhat.y*(1.0-cosP) + yhat.x*sinP) + zhat.z*(cosP+yhat.z*yhat.z*(1.0-cosP) );
        
        phat.normalise();
        
        
    }
    
    //convet the coordinate into target coordinate system
    GVector GSpaceCraftAttitude::convert2Target(GVector target)
    {
        
        // force is in BFS, transform from bfs to eci
        //force =  force_bfs.x * svAttitude->xhat + force_bfs.y * svAttitude->yhat +force_bfs.z * svAttitude->zhat ;
        
       GVector p = target.x*xhat + target.y*yhat + target.z*zhat;
        
//        GVector p;
//        
//        p.x = xhat.x * target.x + yhat.x * target.y + zhat.x * target.z;
//        p.y = xhat.y * target.x + yhat.y * target.y + zhat.y * target.z;
//        p.z = xhat.z * target.x + yhat.z * target.y + zhat.z * target.z;
//        
        return p;
    }
    
    void GSpaceCraftAttitude::updateAttitude(GVector& satpos_eci, GVector& sunpos_eci, GVector& satsun_eci, GVector& satvel_eci, double shadowfactor)
    {
        //first choose spacecraft, then update the eclipse state, finally, update the attitude
        //m_eclipse = eclipse(sunpos, satposvel);
        //m_eclipse = 0;
        //start the control part
        
        
        
        double PI = GCONST("PI");
        
        h = crossproduct(satpos_eci, satvel_eci);  // angular momentum
        n = normalise(h);  // orbital normal
        
        // the local orbital frame RTN
        // T is not the velocity direction
        eR = satpos_eci;
        eR.normalise();
        eN = n;
        // here T is not real the orbital tangent
        // T should point to velocity direction
        eT = crossproduct(eN, eR);
        eT.normalise();
        
        //this definition is from Arnold's paper, new SRP model for GNSS orbit dermination
        //from satellite to sun
        eD = satsun_eci;
        eD.normalise();
        
        eY = crossproduct(eD, eR);
        eY.normalise();
        
        eB = crossproduct(eD, eY);
        eB.normalise();
        
        // this is always the case
        zhat = -eR;
        
        GVector satsun_u = normalise(satsun_eci);
        GVector midnight_u, noon_u;
        
        beta = asin( dotproduct(satsun_u, n )/( satsun_u.norm()*n.norm() ) );
        
        midnight_u = sin(beta)*n - satsun_u;
        noon_u = - midnight_u;
        midnight_u.normalise();
        noon_u.normalise();
        
        int d = 1;
        //calculating the orbital angle measured from midnight
        GVector t = crossproduct(midnight_u, eR) ;
        if( ( normalise(t) - n).norm() >1.0E-8 )
        {
            d = -1;
        }
        
        eta = atan2( t.norm()*d, dotproduct(midnight_u, eR) );
        
        if(eta < 0.0) {eta += M_PI*2;}
        
        if( m_isInitialised == false )
        {
            initialiseAttitude(shadowfactor, satpos_eci, sunpos_eci, satsun_eci, satvel_eci);
        }
        else
        {
            // this is for testing
            //NorminalYawSteering_gal(satpos_eci, sunpos_eci, satsun_eci, satvel_eci);
            //NorminalYawSteering_gps(satpos_eci, sunpos_eci, satsun_eci, satvel_eci);
//            
//            GVector xx = xhat;
//            GVector yy = yhat;
//            GVector zz = zhat;
//            GVector pp = phat;
//            
//            GVector t1 = crossproduct(eT, xx) ;
//            GVector t2 = crossproduct(zz, satsun_u);
//            double mYaw = acos(dotproduct(eT, xx));
//            double myPitch = acos(dotproduct(zz, satsun_u));
//            
//            d = 1;
//            if( ( normalise(t1) - zz).norm() >1.0E-8) // t is the same as n
//            {
//                d = -1;
//            }
//            
//            mYaw = atan2( t1.norm()*d, dotproduct(eT, xx) );
//            d = 1;
//            if( ( normalise(t2) - yy).norm() >1.0E-8) // t is the same as n
//            {
//                d =-1;
//            }
//            
//            myPitch = atan2( t2.norm()*d, dotproduct(zz, satsun_u) );
//
            //****the accuracy is about 1E-5, which is basically wrong *****
            attitudeLaw_Galileo(eta, beta);
            setFromYaw_Pitch(yaw_model, pitch_model);
            
            //conduct the nominal attitude
            //NorminalYawSteering_gal(satpos_eci, sunpos_eci, satsun_eci, satvel_eci);
            
            //test EPS angle
            
            double eps = acos( dotproduct( satsun_eci , -satpos_eci)/satpos_eci.norm()/satsun_eci.norm() );
            
//            double testX = (xx - xhat).norm();
//            double testY = (yy - yhat).norm();
//            double testZ = (zz - zhat).norm();
//            double testP = (pp - phat).norm();
//            
//            if( testX > 1.0E-8 )
//            {
//                printf("X error!\n");
//            }
//            
//            if( testY > 1.0E-8 )
//            {
//                printf("Y error!\n");
//            }
//            if( testZ > 1.0E-8 )
//            {
//                printf("Z error!\n");
//            }
//            if( testP > 1.0E-8 )
//            {
//                printf("P error!\n");
//            }
//            
//            printf("%.6f %.6f %.6f %.6f\n",beta*180/PI,eta*180/PI,yaw_model*180/PI,pitch_model*180/PI);
//////
            
            
        }
        
        
        
        GVector body_frame_sun, satsunvec_eci;
        satsunvec_eci = satsun_eci;
        satsunvec_eci.normalise();
        body_frame_sun.x =  dotproduct(xhat, satsunvec_eci);
        
        body_frame_sun.y =  dotproduct(yhat, satsunvec_eci);
        
        body_frame_sun.z =  dotproduct(zhat, satsunvec_eci);
        
        body_frame_sun.normalise();
        
        //longitude -pi to pi  , body_frame_sun.y 的正负决定-180和+180，光压差别大
        lambda = atan2(body_frame_sun.y, body_frame_sun.x );
        
        //latitude  -pi/2 to pi/2
        phi = 0.0;
        
        if( fabs(body_frame_sun.z - 1.0)<1.0E-14 )
        {
            phi = PI/2.0;
        }
        else if(fabs(body_frame_sun.z + 1.0)<1.0E-14)
        {
            phi = -PI/2.0;
        }
        else
        {
            phi = std::asin(body_frame_sun.z) ;
        }
        
        
        int testc = 0;
        
        
        
    }
    
    
    
    void GSpaceCraftAttitude::NorminalYawSteering_gps(GVector& satpos_eci, GVector& sunpos_eci, GVector& satsun_eci, GVector& satvel_eci)
    {
        zhat = -satpos_eci;
        // the direction vector of solar panel
        phat = satsun_eci;
        
        //yaw-steering attitude
        // the IGS like BFS, +x panel is illuminated by sun
        // but, for Galileo, +x panel is oriented away from sun,
        
        // for GPSIIA, IGS like attitude
        yhat = crossproduct(zhat,phat);
        
        //for Galileo
        //yhat = crossproduct(phat, zhat);
        
        xhat = crossproduct(yhat, zhat);
        
        
        //GMath::crossproduct( m_yhat, m_zhat, m_xhat);
        xhat.normalise();
        yhat.normalise();
        zhat.normalise();
        phat.normalise();
        
        //calculating yaw angle, pitch angle, yaw rate ,yaw acceleration, pitch rate
        //double testp = dotproduct(xhat, phat);
        
        //int testc =0;

    }
    
    
    /*
     *
     * the yaw-steering attitude mode, which is the nomial attitude for most of the GNSS satellites
     * note that: the error in z direction will affect the y and x direction, Is there any method to adjust this error ??
     *
     */
    void GSpaceCraftAttitude::NorminalYawSteering_gal( GVector& satpos_eci, GVector& sunpos_eci, GVector& satsun_eci, GVector& satvel_eci )
    {
        
        zhat = -satpos_eci;
        // the direction vector of solar panel
        phat = satsun_eci;
        
        //yaw-steering attitude
        // the IGS like BFS, +x panel is illuminated by sun
        // but, for Galileo, +x panel is oriented away from sun,
        
        // for GPSIIA, IGS like attitude
        //yhat = crossproduct(zhat,phat);
        
        //for Galileo
        yhat = crossproduct(phat,zhat);
        
        xhat = crossproduct(yhat, zhat);
        
        
        
        //GMath::crossproduct( m_yhat, m_zhat, m_xhat);
        xhat.normalise();
        yhat.normalise();
        zhat.normalise();
        phat.normalise();
        
        //GVector zt = crossproduct(xhat, yhat);
        
        
        //calculating yaw angle, pitch angle, yaw rate ,yaw acceleration, pitch rate
        //double testp = dotproduct(xhat, phat);
        
        //int testc =0;
    }
    
    
    /*
     * this is a special attitude mode for BeiDou satellite, GEO always use this attitude
     * IGSO/MEO will change from Yaw-steering to Yaw-fixed on some conditions
     * Yaw-fixed means the yaw-angle is zero
     */
    void GSpaceCraftAttitude::YawFixed( GVector& satpos_eci, GVector& satvel_eci )
    {
        // because the yaw-angle is zero, namely, the x and velocity have the same direction
        // according the montenbruck's paper: GNSS satellite geometry and attitude models
        //m_zhat[0] = - satpos_u[0]; m_zhat[1] = - satpos_u[1];m_zhat[2] = - satpos_u[2];
        zhat = - satpos_eci;
        yhat = crossproduct(satvel_eci, satpos_eci);
        xhat = crossproduct(yhat, zhat);
        
        xhat.normalise();
        yhat.normalise();
        zhat.normalise();
        
        //GMath::crossproduct( m_yhat, m_zhat, m_xhat);
       // m_yawangle =0.0;
        
    }
    
    
    
    
    // Returns true if there is an intersection, false otherwise
    bool GSpaceCraftAttitude::sun_edge_earth_intersection(GVector rso, GVector edge)
    {
        constexpr double p = 6378.137; //sgnlOPS::wgs84_equatorial_radius; km
        
        constexpr double q = 6356.7523142; //sgnlOPS::wgs84_polar_radius; 6356.7523142 km
        
        double sun_radius = 695700.0;  //km
        constexpr double p2 = p*p;
        constexpr double q2 = q*q;
        
        GVector b = rso - edge;
        
        // Temporary variables used to evaluate B and C
        double temp1 = b.x * rso.y - b.y * rso.x;
        double temp2 = b.x * rso.z - b.z * rso.x;
        
        double B = b.y * q2 * temp1 + b.z * p2 * temp2; // Missing multiply by 2
        
        // We only need the squares of b components for A and C calculations
        b.x *= b.x;
        b.y *= b.y;
        b.z *= b.z;
        
        double A = (b.x + b.y) * q2 + b.z * p2;
        
        double C = (temp1 * temp1 - b.x * p2) * q2 + temp2 * temp2 * p2;
        
        // Check for intersection (discriminant greater than zero)
        // Should be 4AC, but don't include 4, because we didn't multiply B by 2
        return (B * B - A * C > 0.0);
    }
    
    /* @fn eclipse
     * @author Santosh Bhattarai (adapted from routines by S. Adhya & A. Sibthorpe)
     * @date 12 March 2015
     * @return int (value specifies eclipse state)
     * @cite<Adhya04> Sima Adhya, Anthony Sibthorpe, Marek Ziebart and Paul Cross
     *     "Oblate Earth Eclipse State Algorithm for Low-Earth-OrbitingSatellites",
     *     Journal of Spacecraft and Rockets, Vol. 41, No. 1 (2004), pp. 157-159.
     *     doi: 10.2514/1.1485
     *
     * The first version of this function was written by Sima Adhya and Ant
     * Sibthorpe of the Department of Geomatic Engineering, UCL in 2002.
     *
     * This function determines eclipse states based on an (oblate) spheroidal
     * Earth model. Returns int meaning; 0 full phase, 1 penumbra, 2 umbra.
     * sun_eci in ECI
     *
     * the unit of sun_eci and posvel_eci should be in KM
     */
    double GSpaceCraftAttitude::eclipse( GVector& sunpos_eci, GVector& satpos_eci )
    {
        double eclipse_state = 1.0;
        double sun_radius = 695700.0;  //km
        double rso_sun_distance = (sunpos_eci-satpos_eci).norm();
        double eci_sun_distance = sunpos_eci.norm();
        
        // Only run eclipse check if spacecraft is further from sun than Earth is
        if (rso_sun_distance > eci_sun_distance)
        {
            GVector rso = satpos_eci;
            GVector sun = sunpos_eci;
            
            // sun_radius * unit vector in sun, earth, sc plane which is perp.
            // to earth sun vector, defines edges of sun (sun +/- sun_perp)
            GVector sun_perp = sun * dotproduct(sun, rso) - rso * sun.norm2();
            
            // Normalise sun_perp and multiply by sun radius
            sun_perp *= sun_radius  / sun_perp.norm();
            
            // Test sun-edge 1 to spacecraft vector for Earth intersection
            if (sun_edge_earth_intersection(rso, (sun + sun_perp))) {
                eclipse_state -= 0.5;
            }
            
            // Test sun-edge 2 to spacecraft vector for Earth intersection
            if (sun_edge_earth_intersection(rso, (sun - sun_perp))) {
                eclipse_state -= 0.5;
            }
        }
        
        return eclipse_state;
    }

    
    
    
    
} // enf of the namespace

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
        
        GVector sun_u = normalise(sunpos_eci);
        GVector midnight_u, noon_u;
        
        beta = asin( dotproduct(sun_u, n )/( sun_u.norm()*n.norm() ) );
        
        noon_u = sun_u - sin(beta)*n;
        
        //midnight_u = - noon_u;
        //midnight_u.normalise();
        noon_u.normalise();
        
        int d = 1;
        //calculating the orbital angle measured from midnight
        GVector t = crossproduct(noon_u, eR) ;
        if( ( normalise(t) - n).norm() >1.0E-8 )
        {
            d = -1;
        }
        
        // tan(u) = sin(u) / cos(u)
        eta = atan2( t.norm()*d, dotproduct(noon_u, eR) );
        
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
            //attitudeLaw_Galileo(eta, beta);
            //setFromYaw_Pitch(yaw_model, pitch_model);
            
            //conduct the nominal attitude
            NorminalYawSteering_gal(satpos_eci, sunpos_eci, satsun_eci, satvel_eci);
            
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
    bool GSpaceCraftAttitude::sun_edge_earth_intersection(GVector rso, GVector b)
    {
        
        constexpr double a2 = 6378.137*6378.137;
        constexpr double a_b2 = 1.00336408982098*1.00336408982098;
        
        // Matrix A = [[1,0,0],[0,1,0],[0,0,(a/b)^2]]
        double rTAr = (rso.x * rso.x + rso.y * rso.y) + (rso.z * rso.z) * a_b2;
        double rTAb = (rso.x * b.x + rso.y * b.y) + (rso.z * b.z) * a_b2;
        double bTAb = (b.x * b.x + b.y * b.y) + (b.z * b.z) * a_b2;
        
        return (rTAb * rTAb - bTAb * (rTAr - a2));
        
    }
    
    
    // Returns a vector from the rso that satisfies 3 conditions:
    // 1) vector is tangential to the surface of the earth ellipsoid
    // 2) vector is in the same plane as the earth and sun vectors
    // 3) since there are 2 solutions to above, return the one closest to the sun
    // Conditions 2 and 3 are easily satisfied as the solution is made from a
    // linear combination of the rso and sun vectors.
    GVector GSpaceCraftAttitude::get_earth_edge(GVector rso, GVector sun)
    {
         constexpr double a2 = 6378.137*6378.137;
         constexpr double a_b2 = 1.00336408982098*1.00336408982098;
        
        // Matrix A = [[1,0,0],[0,1,0],[0,0,(a/b)^2]]
        double rTAr = (rso.x * rso.x + rso.y * rso.y) + (rso.z * rso.z) * a_b2;
        double rTAs = (rso.x * sun.x + rso.y * sun.y) + (rso.z * sun.z) * a_b2;
        double sTAs = (sun.x * sun.x + sun.y * sun.y) + (sun.z * sun.z) * a_b2;
        
        // Use this to solve for the linear combination of rso and sun vectors
        // that give the Earth edge vectors.
        double D = std::sqrt((rTAr - a2) * a2 / (rTAr * sTAs - rTAs * rTAs));
        
        // Earth edge closest to sun
        GVector edge = ((a2 - rTAs * D) / rTAr - 1.0) * rso + D * sun;
        
        // Testing:
        // This is the second edge, always further away from the sun:
        // Cartesian edge2 = ((a2 + rTAs * D) / rTAr - 1.0) * rso - D * sun;
        
        // double inter1 = sun_edge_earth_intersection(rso, edge);
        // double inter2 = sun_edge_earth_intersection(rso, edge2);
        // std::cout << "intersection1: " << inter1 << std::endl;
        // std::cout << "intersection2: " << inter2 << std::endl;
        // std::cout << "------------------------------" << std::endl;
        
        return edge;
    }
    
    
    // Adapted from conical shadow function from Montenbruck 3.4.2, pages 82 - 83
    // Angles a, b and c are calculated differently.
    // In particular angle b is calculated based on the Earth "radius" at the point
    // of eclipse from the WGS84 ellipsoidal model.
    double GSpaceCraftAttitude::penumbra_flux_scale(GVector rso, GVector earth_edge,
                                              GVector sun, GVector sun_edge)
    {
        const double sun2 = sun.norm2();
        const double rso2 = rso.norm2();
        
        const double a = std::acos(dotproduct(sun_edge, sun) /
                                   std::sqrt(sun_edge.norm2() * sun2));
        const double b = std::acos(-dotproduct(earth_edge, rso) /
                                   std::sqrt(earth_edge.norm2() * rso2));
        const double c = std::acos(-dotproduct(rso, sun) / std::sqrt(rso2 * sun2));
        
        const double a2 = a * a;
        const double b2 = b * b;
        
        double factor;
        
        if (c >= a + b) {
            // No eclipse
            factor = 1.0;
        } else if (b - a >= c) {
            // Sun completely occluded
            factor = 0.0;
        } else if (a - b >= c) {
            // Annular eclipse, ie: sun is visible around the occluding body
            factor = 1.0 - b2 / a2;
        } else {
            double x = (c * c + a2 - b2) / (2.0 * c);
            double y = std::sqrt(a2 - x * x);
            double A = a2 * std::acos(x / a) + b2 * std::acos((c - x) / b) - c * y;
            factor = 1.0 - A / (GCONST("PI") * a2);
        }
        
        factor = std::min(factor, 1.0);
        factor = std::max(factor, 0.0);
        
        return factor;
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
    double GSpaceCraftAttitude::eclipse( GVector rso, GVector sun )
    {
        
        double rso_sun_distance = (rso - sun).norm();
        double eci_sun_distance = sun.norm();
        
        constexpr double R2 = 695700.0 * 695700.0;
        
        double eclipse_state = 1.0;
        
        // Only run eclipse check if spacecraft is further from sun than Earth is.
        if (rso_sun_distance > eci_sun_distance) {
            double sTs = sun.norm2();
            double sTr = dotproduct(sun, rso);
            double rTr = rso.norm2();
            
            // Use this to solve for the linear combination of rso and sun vectors
            // that give the sun edge vectors.
            double D = std::sqrt((sTs - R2) * R2 / (sTs * rTr - sTr * sTr));
            GVector Drso = D * rso;
            
            // Create the sun edge vector closest to Earth.
            GVector sun_edge1 = (1.0 + (sTr * D - R2) / sTs) * sun - Drso;
            
            if (sun_edge_earth_intersection(rso, sun_edge1) > 0.0) {
                
                // Create the sun edge vector furthest from Earth.
                GVector sun_edge2 = (1.0 - (sTr * D + R2) / sTs) * sun + Drso;
                
                if (sun_edge_earth_intersection(rso, sun_edge2) >= 0.0) {
                    // Both sun edges are blocked, so no flux.
                    eclipse_state = 0.0;
                } else {
                    // Find the Earth edge vector that lies between the sun edges
                    GVector earth_edge = get_earth_edge(rso, sun);
                    
                    // Only one edge of the sun is blocked, find the penumbral flux.
                    // eclipse_state =
                    //     penumbra_flux_scale(sun_edge1, earth_edge, sun_edge2);
                    eclipse_state =
                    penumbra_flux_scale(rso, earth_edge, sun, sun_edge1);
                }
            }
        }
        
        return eclipse_state;
        
        
    }

    
    
    
    
} // enf of the namespace

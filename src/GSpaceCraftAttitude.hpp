//
//  GSpaceCraftAttitude.hpp
//  GFC
//
//  Created by lizhen on 16/05/2016.
//  Copyright © 2016 lizhen. All rights reserved.
//

#ifndef GSpaceCraftAttitude_hpp
#define GSpaceCraftAttitude_hpp

#include <stdio.h>
#include "GEarthOrientationParameter.hpp"
#include "GTime.h"
#include "GMath.hpp"
#include "GMatrix.h"
#include "GVector.hpp"

namespace gfc
{
    /*
     * usually, the attitude of the spacecraft is describled in ECI, that means to determined the three coordinate axis of BFS in ECI coordinate
     * generally, these two coordinate system have the same origin, thus only need to determine the three oure angles
     * different spacecraft have different laws to contrl the attitude
     */
   class GSpaceCraftAttitude
    {
        
    public:
        
        static double eclipse( GVector& sunpos_eci, GVector& satpos_eci );  //units: km
        static bool sun_edge_earth_intersection(GVector rso, GVector edge);
        
        GSpaceCraftAttitude()
        {
            m_isInitialised = false;
        }
        
        void convert2ECEF(GEarthOrientationParameter& eop);
        
        /*
         acquire the attitude the first time
         if it is eclipse, the satellite can NOT know where sun is thus return false
        */
        bool initialiseAttitude(double shadowFactor, GVector& satpos_eci,
                                GVector& sunpos_eci,
                                GVector& satsun_eci,
                                GVector& satvel_eci,double yaw_angle = 0.0, double pitch_angle = 0.0);
        
        
        /*before calling this function, z vector and T vector is known*/
        void setFromYaw_Pitch( double yaw, double pitch );
        
        //convet the coordinate into target coordinate system
        GVector convert2Target(GVector target);
        
        void updateAttitude( GVector& satpos_eci, GVector& sunpos_eci, GVector& satsun_eci, GVector& satvel_eci, double shadowfactor);
        
        
       /*
         * function: attitude control law of the spacecraft, namely the attitude dynamics
         * sunpos is the position(center of mass ) of sun in ECI
         * satpos is the position(center of mass ) of sat in ECI
         * satvel is the velocity of sat in ECI
         *
        */
        //void attitudeLaw_BDS_IGSO( GVector& satpos_eci, GVector& sunpos_eci, GVector& satsun_eci, GVector& satvel_eci );
        
        void attitudeLaw_Galileo(double eta, double beta);
        
        /*
          *
          * the nominal yaw-steering attitude mode, which is the nomial attitude for most of the GNSS satellites
          * note that: the error in z direction will affect the y and x direction, Is there any method to adjust this error ??
          *
        */
        void NorminalYawSteering_gps(  GVector& satpos_eci, GVector& sunpos_eci, GVector& satsun_eci, GVector& satvel_eci );
        
        void NorminalYawSteering_gal(  GVector& satpos_eci, GVector& sunpos_eci, GVector& satsun_eci, GVector& satvel_eci );
        
         /*
          * this is a special attitude mode for BeiDou satellite, GEO always use this attitude
          * IGSO/MEO will change from Yaw-steering to Yaw-fixed on some conditions
          * Yaw-fixed means the yaw-angle is zero
          */
        void YawFixed( GVector& satpos_eci, GVector& satvel_eci );
        
        GVector xhat;  // the BFS x in ECI
        GVector yhat;  // the BFS y in ECI
        GVector zhat;  // the BFS z in ECI
        GVector phat;  // the solar panel normal in ECI
        
        GVector h;  // the angular momentum of orbital plane
        GVector n;  // normal vector of orbital plane
        
        GVector eD;  // the D vector of DYB coordinate system in ECI
        GVector eY;  // the Y vector of DYB coordinate system in ECI
        GVector eB;  // the B vector of DYB coordinate system in ECI
        
        // local orbital frame
        GVector eR; // radial direction
        GVector eT; // perpendicular to radial direction (not exactly the tangent)
        GVector eN; // normal direction
        
        
        double beta;  // the sun elevation
        double eta;   // the orbital angle measured from midnight
        double phi;    // the latitude of sun in bfs
        double lambda; // the longitude of sun in bfs
        

        /* 
         
         ref:http://www.sciencedirect.com/science/article/pii/S0273117715004378
         yaw angle is the angle between Trtn and Xbfs
         
        */
        
        double yaw_nomial;  // the nominal yaw angle
        double yaw_model;    // yaw angle form attitude model
        double yaw_rate_model; // yaw angular rate
        double yaw_acc_model; // yaw angular acceleration
        // SADM(solar array drive mechanism), solar panel rotation angle
        double pitch_nominal; // the pitch angle from nomial
        double pitch_model;  // the pitch angle from attitude model
        
    private:
        
        bool m_isInitialised; // whether the attitude is initialised
        
    };
    
    
} // end of the namespace



#endif /* GSpaceCraftAttitude_hpp */

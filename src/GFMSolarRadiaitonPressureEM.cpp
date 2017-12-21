//
//  GFMSolarRadiaitonPressureEM.cpp
//  GFC
//
//  Created by lizhen on 26/10/2016.
//  Copyright © 2016 lizhen. All rights reserved.
//

#include "GFMSolarRadiaitonPressureEM.hpp"
namespace gfc
{
    GFMSolarRadiationPressureEM::GFMSolarRadiationPressureEM()
    {
        setForceName("GFMEMP");
        
        num_param = 5;  //10 for +-z , 7 for the common z panel, 5 for ecom
        
        parameters.resize(num_param);
        m_dadp.resize(3, num_param);
        
        //unit: nano meter/s^2
        // the d0 parameter for 0.1 solar panel reflectivity of E11
        //parameters[0] = -15.35E-9;
        // ybias test
        //parameters[1] = 1.65;
        
        
    }
    
    //param[5]: D0, Y0,B0, Bc and Bs, u is the satellite's argument of latitude
    // here factor should be distance factor times shadow_factor
    GVector GFMSolarRadiationPressureEM::reducedECOM(double* param,double u, double factor )
    {
        GVector acc_dyb,acc_xyz;
        //D
        acc_dyb.x = param[0];
        //Y
        acc_dyb.y = param[1];
        //B
        acc_dyb.z = param[2] + param[3]*cos(u) + param[4]*sin(u);
        
        acc_dyb = acc_dyb*factor;
        
        return acc_dyb*1.0E-9;
        
    }
    
    /*  da_eci/dp = dT/dp*a_dyb + T* da_dyb/dp. T is the transformation matrix from DYB to ECI
     *  the parameters in ECOM model is independent from T
     *  thus, dT/dp = 0, i.e da_eci/dp = T*da_dyb/dp
     *
     *   da_dby     | 1  0  0  0     0    |
     *   ------ =   | 0  1  0  0     0    |
     *     dp       | 0  0  1  cosu  sinu |
     *
     */
    
    void GFMSolarRadiationPressureEM::reducedECOM_dadp(double u, GMatrix& daDYBdp )
    {
        daDYBdp(0,0) = 1.0; daDYBdp(1,1) = 1.0;
        daDYBdp(2,2) = 1.0;
        daDYBdp(2,3) = cos(u);
        daDYBdp(2,4) = sin(u);
        
        daDYBdp = daDYBdp*1.0E-9;
        
    }
    
    
    
    /*
     the Body Fixed Frame version of ECOM
     */
    GVector GFMSolarRadiationPressureEM::reducedECOM_XYZ(double* param,double phi, double factor )
    {
        GVector acc_xyz;
        //X
        acc_xyz.x = -param[0]*cos(phi) + param[2]*sin(phi) + 0.5*param[3]*(1.0-cos(2.0*phi)) + 0.5*param[4]*sin(2.0*phi);
        //Y
        acc_xyz.y = -param[1];
        //Z
        acc_xyz.z = param[0]*sin(phi) + param[2]*cos(phi) + 0.5*param[3]*sin(2.0*phi) + 0.5*param[4]*(cos(2.0*phi)+1.0);
        
        acc_xyz = acc_xyz*factor;
        
        return acc_xyz*1.0E-9;
        
    }
    
    void GFMSolarRadiationPressureEM::reducedECOM_dadp_XYZ(double phi, GMatrix& daXYZdp )
    {
        daXYZdp(0,0) = -cos(phi); daXYZdp(0,1) = 0.0; daXYZdp(0,2) = sin(phi); daXYZdp(0,3) = 0.5*(1.0-cos(2.0*phi));
        daXYZdp(0,4) = 0.5*sin(2.0*phi);
        daXYZdp(1,1) = -1.0;
        daXYZdp(2,0) = sin(phi);daXYZdp(2,1) = 0.0;daXYZdp(2,2) = cos(phi);daXYZdp(2,3) = 0.5*sin(2.0*phi);daXYZdp(2,4) = 0.5*(cos(2.0*phi)+1.0);
        
        daXYZdp = daXYZdp*1.0E-9;
        
    }

    
    GVector GFMSolarRadiationPressureEM::constantDYB(double *param, double u, double factor)
    {
        GVector acc_dyb;
        //D
        acc_dyb.x = param[0];
        //Y
        acc_dyb.y = param[1];
        //B
        acc_dyb.z = param[2];
        
        acc_dyb = acc_dyb*factor;
        
        return acc_dyb*1.0E-9;
    }
    
    void GFMSolarRadiationPressureEM::constantDYB_dadp(GMatrix& daDYBdp)
    {
        daDYBdp(0,0) = 1.0;
        daDYBdp(1,1) = 1.0;
        daDYBdp(2,2) = 1.0;
        daDYBdp = daDYBdp*1.0E-9;
    }
    
    /*
     
     
     */
    GVector GFMSolarRadiationPressureEM::empXYZ_DYB(double* param,double phi, double factor)
    {
        double delta = 0.0;
        if( phi >= 0.0 )
        {
            delta = 1.0;
        }
        else
        {
            delta = -1.0;
        }
        
        GVector acc_dyb;
        
        //D
        acc_dyb.x = -param[0]*cos(phi)*cos(phi) - param[1]*cos(phi)*(1.0+cos(2*phi)) - param[2]*delta*cos(phi)*sin(2.0*phi) - param[3]*sin(phi)*sin(phi) + param[4]*delta*sin(phi)*(cos(2.0*phi)-1.0) - param[5]*sin(phi)*sin(2.0*phi);
        
        //Y
        acc_dyb.y = -param[6];
        //B
        acc_dyb.z = param[0]*sin(phi)*cos(phi) + param[1]*sin(phi)*(1.0+cos(2.0*phi)) + delta*param[2]*sin(phi)*sin(2.0*phi) - param[3]*sin(phi)*cos(phi) + param[4]*delta*cos(phi)*(cos(2.0*phi)-1.0) - param[5]*cos(phi)*sin(2.0*phi) ;

         acc_dyb = acc_dyb*factor;
        
        return acc_dyb*1.0E-9;
        
    }
    
    void GFMSolarRadiationPressureEM::empXYZ_dadp_DYB(double phi,GMatrix& daXYZdp)
    {
        daXYZdp.clear();
        
        double delta = 0.0;
        if(phi >= 0.0)
        {
            delta = 1.0;
        }
        else
        {
            delta = -1.0;
        }
        
        daXYZdp(0,0) = -cos(phi)*cos(phi);
        daXYZdp(0,1) = -cos(phi)*(1.0+cos(2.0*phi));
        daXYZdp(0,2) = -delta*sin(2.0*phi)*cos(phi);
        daXYZdp(0,3) = -sin(phi)*sin(phi);
        daXYZdp(0,4) = delta*sin(phi)*(cos(2.0*phi) - 1.0);
        daXYZdp(0,5) = -sin(phi)*sin(2.0*phi);
        
        daXYZdp(1,6) = -1.0;
        
        daXYZdp(2,0) = sin(phi)*cos(phi);
        daXYZdp(2,1) = sin(phi)*(1.0+cos(2.0*phi));
        daXYZdp(2,2) = delta*sin(2.0*phi)*sin(phi);
        daXYZdp(2,3) = -sin(phi)*cos(phi);
        daXYZdp(2,4) = delta*cos(phi)*(cos(2.0*phi) - 1.0);
        daXYZdp(2,5) = -cos(phi)*sin(2.0*phi);
        
        daXYZdp = daXYZdp*1.0E-9;
        
    }
    
    /*
     model description in Body Fixed System
     x =  P1*cos(phi) + P2*cos(phi)*cos(phi) - d2*P3*sin(2phi) + d1*P4*sin(2phi)
     z =  d2*P5*sin(phi) + d2*P7*sin(phi)*sin(phi) - P9*sin(2phi) + d1*P6*sin(phi)- d1*P8*sin(phi)*sin(phi) ;
     y =  P10
     
     if phi >= 0.0, d1=1, d2 = 0 ; if phi <0.0, d1=0, d2=1;
     
     new empirical model test
     phi: the latitude of sun in bfs, -90 : 90
     factor: the scaling factor
     param: the parameters in this model
     */

    GVector GFMSolarRadiationPressureEM::empXYZ(double* param,double phi, double factor)
    {
        GVector acc_xyz;
        double d1 = 0, d2 = 0, mass = 696.815;
        double W_C = 1368/299792458.0;
        double cp = cos(phi);
        double cp2 = cp*cp;
        double c2p = cos(2.0*phi);
        
        double sp  = sin(phi);
        double sp2 = sp*sp;
        double s2p = sin(2.0*phi);
        
        if (phi >= 0)
        {
            d1 = 1;
            d2 = 0;
        }
        else
        {
            d1 = 0.0;
            d2 = 1.0;
        }
        
        acc_xyz.x =  param[0]*cp + param[1]*cp2 - d2*param[2]*s2p + d1*param[3]*s2p;
        acc_xyz.y =  param[9];  // y bias
        acc_xyz.z = d2*param[4]*sp + d1*param[5]*sp + d2*param[6]*sp2 - d1*param[7]*sp2 - param[8]*s2p;
        
        // y bias does not apply the factor
        acc_xyz.x = acc_xyz.x*factor;
        acc_xyz.z = acc_xyz.z*factor;
        
        return acc_xyz*1.0E-9;
        
    }
    
    void GFMSolarRadiationPressureEM::empXYZ_dadp(double phi,GMatrix& daXYZdp)
    {
        double d1 = 0, d2 = 0,mass = 696.815;
        double cp = cos(phi);
        double cp2 = cp*cp;
        double c2p = cos(2.0*phi);
        
        double sp  = sin(phi);
        double sp2 = sp*sp;
        double s2p = sin(2.0*phi);

        if (phi >= 0)
        {
            d1 = 1;
            d2 = 0;
        }
        else
        {
            d1 = 0.0;
            d2 = 1.0;
        }
        
        double W_C = 1368/299792458.0;
        
        // x component
        daXYZdp(0,0) = cp;
        daXYZdp(0,1) = cp2;
        daXYZdp(0,2) = -d2*s2p;
        daXYZdp(0,3) = d1*s2p;
        // y component
        daXYZdp(1,9) = 1.0;
        // z component
        daXYZdp(2,4) = d2*sp;
        daXYZdp(2,5) = d1*sp;
        daXYZdp(2,6) = d2*sp2;
        daXYZdp(2,7) = -d1*sp2;
        daXYZdp(2,8) = -s2p;
        
        daXYZdp = daXYZdp;
        
        //cout<<daXYZdp;
        
        daXYZdp = daXYZdp*1.0E-9;
        
    }
    
    /*
     This is to validate the simplified empirial model in DYB system
     D = D0
     B = (D0 + P1)*cot(phi) + P2*cos(phi)*cos(phi)/sin(phi) - 2*d2*P3*cos(phi) + 2*d1*P4*cos(phi)
     it is still under developing
     */
    GVector GFMSolarRadiationPressureEM::empDYB(double* param,double phi, double factor)
    {
        double d1 = 0, d2 = 0,mass = 696.815;
        
        if (phi >= 0)
        {
            d1 = 1;
            d2 = 0;
        }
        else
        {
            d1 = 0.0;
            d2 = 1.0;
        }
        
        double W_C = 1368/299792458.0/mass;
        
        GVector acc_dyb;
        
        acc_dyb.x = param[0];
        acc_dyb.y = param[5];
        acc_dyb.y = param[1]/tan(phi) + param[2]*cos(phi)/tan(phi) - 2*d2*param[3]*cos(phi) + 2*d1*param[4]*cos(phi);
        
        // y bias does not to apply the factor
        acc_dyb.x = acc_dyb.x*factor*W_C;
        acc_dyb.z = acc_dyb.z*factor*W_C;
        
        return acc_dyb;
        
    }
    
    void GFMSolarRadiationPressureEM::empDYB_dadp(double phi,GMatrix& daXYZdp)
    {
        double delta = 0.0;
        if(phi >= 0.0)
        {
            delta = 1.0;
        }
        else
        {
            delta = -1.0;
        }

        
        daXYZdp.clear(); // set it to zeros
        
        
        
    }
    
    
    
    
    //must implement doCompute function !!!!
    void GFMSolarRadiationPressureEM::doCompute(GSpaceCraft* spacecraft )
    {
        double dis_factor = spacecraft->getStatePointer()->dis_factor;
        
        double shadow_factor = spacecraft->getStatePointer()->shadow_factor;
        
        double factor = dis_factor*shadow_factor;
        
        double eta = spacecraft->getStatePointer()->attitude_eci.eta;
        
        //double u = spacecraft->getStatePointer()->u;
        
        
        double eps = spacecraft->getStatePointer()->eps;
        
        double phi = spacecraft->getStatePointer()->attitude_eci.phi;
        
        //if(phi >=0.0) phi = M_PI_2 - phi;
        //if(phi<0.0)   phi = M_PI_2 + phi;
        
        GVector acc_xyz, acc_dyb, acc_bfs;
        
        double PI = GCONST("PI");
        
        //this 5 parameter ECOM has to use u as the angle
        acc_dyb = reducedECOM(&parameters[0], eta , factor);
        
        //acc_bfs  = reducedECOM_XYZ(&parameters[0],phi,factor);
        
        //acc_dyb = constantDYB(&parameters[0], u, factor);
        
        //acc_dyb = empXYZ_DYB(&parameters[0], phi, factor);
        
        //acc_bfs = empXYZ(&parameters[0],  phi, factor);
        
        //acc_bfs = empXYZ2(&parameters[0], phi, factor);
        
        
        // transform from DYB coordinate system to ECI system
        GVector& ed = spacecraft->getStatePointer()->attitude_eci.eD;
        GVector& ey = spacecraft->getStatePointer()->attitude_eci.eY;
        GVector& eb = spacecraft->getStatePointer()->attitude_eci.eB;
        GMatrix T(3,3);
        T(0,0) = ed.x ; T(0,1) = ey.x; T(0,2) = eb.x;
        T(1,0) = ed.y ; T(1,1) = ey.y; T(1,2) = eb.y;
        T(2,0) = ed.z ; T(2,1) = ey.z; T(2,2) = eb.z;
        //they are in m/s^2
        acc_xyz.x = T(0,0) * acc_dyb.x + T(0,1)*acc_dyb.y + T(0,2)*acc_dyb.z;
        acc_xyz.y = T(1,0) * acc_dyb.x + T(1,1)*acc_dyb.y + T(1,2)*acc_dyb.z;
        acc_xyz.z = T(2,0) * acc_dyb.x + T(2,1)*acc_dyb.y + T(2,2)*acc_dyb.z;
        
        
        /*
        // transform from bfs coordinate system to ECI system
        GVector& eb = spacecraft->getStatePointer()->attitude_eci.xhat;
        GVector& ef = spacecraft->getStatePointer()->attitude_eci.yhat;
        GVector& es = spacecraft->getStatePointer()->attitude_eci.zhat;
        GMatrix T(3,3);
        T(0,0) = eb.x ; T(0,1) = ef.x; T(0,2) = es.x;
        T(1,0) = eb.y ; T(1,1) = ef.y; T(1,2) = es.y;
        T(2,0) = eb.z ; T(2,1) = ef.z; T(2,2) = es.z;
        //acc_xyz = acc_bfs.x* eb + acc_bfs.y*ef + acc_bfs.z*es;
        //they are in m/s^2
        acc_xyz.x = T(0,0) * acc_bfs.x + T(0,1)*acc_bfs.y + T(0,2)*acc_bfs.z;
        acc_xyz.y = T(1,0) * acc_bfs.x + T(1,1)*acc_bfs.y + T(1,2)*acc_bfs.z;
        acc_xyz.z = T(2,0) * acc_bfs.x + T(2,1)*acc_bfs.y + T(2,2)*acc_bfs.z;
        */


        
        
        
        // transform from acceleration to force, it is in unit of Newton
        acc_xyz = acc_xyz * spacecraft->getSpaceCraftGemotry()->m_mass;
        
        setForce(acc_xyz); // this is force in unit Newton
        
        //set m_dadr, m_dadv and m_dadp for this force model
        
        if( m_hasPartialDerivatives == true )
        {
            
            GMatrix mydadp(3,num_param); // in local system
            
            reducedECOM_dadp( eta , mydadp);
            
            //reducedECOM_dadp_XYZ(phi, mydadp);
            
            //constantDYB_dadp( mydadp);
            
            //empXYZ_dadp_DYB(phi,mydadp);
            
            
            
           // empXYZ_dadp(phi, mydadp);
            
           // empXYZ2_dadp(phi, mydadp);
            
            
            //cout<<mydadp;
            
            //set the m_dadp, according to the law of matrix differential
            // dAX/dX = AT
            m_dadp = factor*(T*mydadp); // m_dadp in ECI
            
            //cout<< m_dadp;
            
        }
        
    }
    
    
    
} // end of namespace

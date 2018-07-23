//
//  GForceModel.hpp
//  GFC
//
//  Created by lizhen on 16/4/11.
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


#ifndef GForceModel_hpp
#define GForceModel_hpp

#include <stdio.h>
#include "GFCCONST.h"
#include "GString.h"
#include "GVector.hpp"
#include "GSpaceEnv.hpp"

namespace gfc {

// all the force should be describled in ECI
class GForceModel
    {
        
 public:
        
    GForceModel()
        {
            m_dadr.resize(3, 3);
            m_dadv.resize(3, 3);
            m_hasPartialDerivatives = false;
            num_param = 0;
        }
        
    ~GForceModel() {}
        
  // this function must be inheried by the children force model classes
    virtual void doCompute(){};
    
    void hasPartialDerivatives(bool pd) {m_hasPartialDerivatives = pd;}
     
    void setModelParameters(int num, double* param)
        {
            if( num != 0 )
            {
                parameters.resize(num);
                memcpy(&parameters[0],param,sizeof(double)*num);
                m_dadp.resize(3, num);
            }
        }
        
    // get the force in three directions
    GVector getForce() { return m_force;}
    
    void setForce(GVector force) { m_force = force; }
    
    void setForceName(GString name) { m_name = name; }
    
    GString getForceName() { return m_name; }
        
    GMatrix m_dadr; // 3 by 3, should be in ECI
    GMatrix m_dadv; // 3 by 3, should be in ECI
    
     // should be in ECI
    GMatrix m_dadp; // the derivatives of acceleartion w.r.t force model parameters, 3 by np, depends on different force models
    
    std::vector<double> parameters; // the initial parameters for this force model
    int num_param; // the number of parameters for this force model
        
    bool  m_hasPartialDerivatives; // it is false by default
        
 private:
  
  GString m_name;     // the name of forceModel
  
  
        
  GVector m_force;  // the acceleration in 3 directions,unit m/s2
  
};

    


    
// the class for antenna thrust
// basicly, this force is descrilbed in Body-Fixed-System
class GFMAntennaThrust : public GForceModel {
 public:
  GFMAntennaThrust() {
    setForceName("GFMANT");
  }

  ~GFMAntennaThrust(){};

  /*
   set the power in three axis directions
   */
  void setPower( GVector power )
    {
        m_powerRate = power;
    }
    
  // inheried from GForceModel class
  void doCompute()
    {
    
        double CLIGHT = GCONST("CLIGHT");
        
        GVector force;
        // the area of antenna disk
        //double r = 0.8;
        //double area = 3.14159265357*r*r;
        force = -m_powerRate/CLIGHT;
        
        setForce(force);
    }

 private:
  // power rate in the body fixed system direction, should be in watts
  GVector m_powerRate;
    
};

   class GFMGeneralRelativity : public GForceModel
    {
        
    public:
        
        GFMGeneralRelativity()
        {
            setForceName("GFMGR");
        }
        
        void doCompute(double GM,  GVector& satpos_eci, GVector& satvel_eci)
        {
            GVector force;
            
            force =generalRelativityCorrection( GM,  satpos_eci, satvel_eci);
            
            setForce(force);
        }
        
        // the general relativity correction
        GVector generalRelativityCorrection(double GM,  GVector& satpos_eci, GVector& satvel_eci)
        {
            GVector force_ecef, force_eci;
            double c = GCONST("CLIGHT")/1000.0; //unit: km
            double C2 = c*c;
            double GMc2 = GM/C2;  // the unit of m_GM is km^3/s^2 ; the unit of minusGMc2 is meter
            double r2 =  satpos_eci.norm2();//  pow(pv_ecef[0],2.0) + pow(pv_ecef[1],2.0) + pow(pv_ecef[2],2.0);
            double v2 =  satvel_eci.norm2(); // pow(pv_ecef[3],2.0) + pow(pv_ecef[4],2.0) + pow(pv_ecef[5],2.0);
            double r = sqrt(r2);
            
            GVector er = satpos_eci;
            er.normalise();
            GVector ev = satvel_eci;
            ev.normalise();
            double rdv = dotproduct(er, ev);
            
            //here acc should be in km/s^2
            force_eci = -GM/r2*( (4.0*GMc2/r - v2/C2)*er + 4.0*v2/C2*rdv*ev   );
            
            /*
            double minusGMc2r3 = minusGMc2 / (r * r2);
            double pvd = 0.0; // unit: km^2
            
            pvd = dotproduct(satpos_eci, satvel_eci);
            
            double r_coef = minusGMc2r3 * (4.0 * GM / r - v2); //unit: 1.0/s^4
            double v_coef = minusGMc2r3 * 4.0 * pvd;
            
            //double a[3] = {0.0}; // in ECEF, expected to be in km/s^2
            //acceleration in m/s^2
            force_eci.x = (r_coef * satpos_eci.x + v_coef * satvel_eci.x )*1000.0;
            force_eci.y = (r_coef * satpos_eci.y + v_coef * satvel_eci.y )*1000.0;
            force_eci.z = (r_coef * satpos_eci.z + v_coef * satvel_eci.z )*1000.0;
            */
            
            //turn km/s^2 to m/s^2
            
            return force_eci*1000.0;
        }
        
    };

    
    
    //Ybias force
    class GFMBFSbias : public GForceModel
    {
        
    public:
        GFMBFSbias()
        {
            setForceName("GFMBFSbias");
        }
        
        // ybias is in nm s^{-2}
        void setParam(GVector bias)
        {
            m_bias = bias*1.0E-9;
        }
        
        void doCompute(GVector& xhat, GVector& yhat, GVector& zhat)
        {
            //force in ECI
            GVector force = m_bias.x*xhat + m_bias.y*yhat + m_bias.z*zhat;
            setForce(force);
        }
        
    private:
        GVector m_bias;  // unit: Newton
        
    };
    
    
    
    

}  // end of namespace gfc

#endif /* GForceModel_hpp */

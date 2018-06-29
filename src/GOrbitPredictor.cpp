
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
//  GOrbitPredictor.cpp
//  GFC
//
//  Created by lizhen on 21/05/2016.
//  Copyright © 2016 lizhen. All rights reserved.
//

#include "GOrbitPredictor.hpp"
#include "GSpacecraft.hpp"
#include "GRungeKutta.hpp"

namespace gfc
{
    
    GOrbitPredictor::GOrbitPredictor(GSpaceCraft* spacecraft)
    {
        m_integrator = NULL;
        
        m_integrator = new GRungeKuttaFehlberg();
        
        //m_integrator = new GAdams(13);
        
        m_spaceCraft = spacecraft;
        
        //m_orbitLog.open("orbit_prediction_state.log");
        m_logfilename = "orbit_prediction_state.log";
        log_on = false;
    }
    
    
    GOrbitPredictor::~GOrbitPredictor()
    {
        if( m_integrator != NULL )
        {
            //note that, must NOT use delete[] , just use delete
            delete m_integrator;
            m_integrator = NULL;
        }
        
        m_orbitLog.close();
    }
    
    
    void GOrbitPredictor::setLogfile(gfc::GString logfilename)
    {
        //m_orbitLog.open(logfilename);
        m_logfilename = logfilename;
        if(m_logfilename != "")
        {
            log_on = true;
        }
    }
    
    void GOrbitPredictor::setStepsize(double stepsize)
    {
        m_integrator->setStepsize(stepsize);
    }
    
    
    void GOrbitPredictor::setSpaceCraft(GSpaceCraft* spacecraft)
    {
        m_spaceCraft = spacecraft;
    }
    
    //return the orbit diff in R T N directions
    // all the inputs must be in ECI
    // p1 and v1 are the reference orbit
    GVector GOrbitPredictor::orbitdiff(GVector& p1, GVector& v1, GVector& p2, GVector& v2)
    {
        GVector diff;
        
        GVector R = normalise(p1);
        GVector C = normalise(crossproduct(p1, v1));
        GVector A = crossproduct(C, R);
        A.normalise();
        
        //force = force_bfs.x * svAttitude->xhat + force_bfs.y * svAttitude->yhat +force_bfs.z * svAttitude->zhat ;
        GVector d = p2 - p1;
        
        diff.x = d.x * R.x + d.y * R.y + d.z * R.z;
        diff.y = d.x * A.x + d.y * A.y + d.z * A.z;
        diff.z = d.x * C.x + d.y * C.y + d.z * C.z;
        
        //diff = p1 - p2;
        
        return  diff;
    }
    
    GForceModelMgr GOrbitPredictor::getForceManager()
    {
        return m_forceManager;
    }
    
    void GOrbitPredictor::setForceManager( GForceModelMgr& manager)
    {
        m_forceManager = manager;
        
        //set hasPartialDerivatives FALSEE
        //m_forceManager.m_forceModels.size();
        std::map< GString, GForceModel* >::iterator   it = m_forceManager.m_forceModels.begin();
        for( ;it != m_forceManager.m_forceModels.end(); ++it )
        {
            it->second->hasPartialDerivatives(false);
        }
        
        // as a predictor, there is no need to estimate any parameters
        if(m_forceManager.totalNum_parameters!= 0 )
        {
            m_forceManager.totalNum_parameters = 0;
        }
        
        
        if(log_on == true)
        {
            //open log file
            m_orbitLog.open(m_logfilename);
            
            m_orbitLog << getLogHeader();
        }
        
        
    }
    
    //this function is used in the class GRungekutta for providing the acceleration
    // usually x is time in unit of seconds
    void  GOrbitPredictor::getDerivatives( int n, double x, double *y, double *dydx )
    {
        //        // need to reset the time each step for the calculation of gravity force
        GVector p, v;
        p.x = y[0]; p.y = y[1]; p.z = y[2];
        v.x = y[3]; v.y = y[4]; v.z = y[5];
        
        //need to reset the time each step for the calculation of gravity force
        GTime ct = m_t0 + x;  // for each step, should be in UTC
        //        //need to update the space environment and the motion state at the same time
        GSpaceEnv::updateSpaceEnvironment(ct);  //update the
        //need to update the spacecraft state in every step
        m_spaceCraft->getStatePointer()->updateState_eci(ct,p,v);
        
        // without partial derivatives
        m_forceManager.getDerivatives(ct,n, x, y, dydx, m_spaceCraft);
        
    }
    
    
    GString GOrbitPredictor::getLogHeader()
    {
        GString res = "timeUTC,ECI_PX(m),ECI_PY(m),ECI_PZ(m),ECI_VX(m/s),ECI_VY(m/s),ECI_VZ(m/s),";
        
        res += "EPS,";
        
        res += "beta,";
        
        res += "eta,"; // the angle of satellite measured from the midnight in orbital plane
        
        res += "phi,"; // the phi angle used in radiation pressure calculation
        
        res += "semi-major-axis,eccentricity,inclination,argument_of_perigee,longitude-of-ascending-node,true-anomaly,";
        std::map< GString, GForceModel* >::iterator   it = m_forceManager.m_forceModels.begin();
        for( ;it != m_forceManager.m_forceModels.end(); ++it )
        {
            GString forceName = it->second->getForceName();
            res += forceName +"_X,";
            res += forceName +"_Y,";
            res += forceName +"_Z,";
        }
        
        res += "SolarFlux_X(w/m^2),";
        res += "SolarFlux_Y(w/m^2),";
        res += "SolarFlux_Z(w/m^2),";
        
        res += "EarthFlux(longwave)_X(w/m^2),";
        res += "EarthFlux(longwave)_Y(w/m^2),";
        res += "EarthFlux(longwave)_Z(w/m^2),";
        
        res += "EarthFlux(shortwave)_X(w/m^2),";
        res += "EarthFlux(shortwave)_Y(w/m^2),";
        res += "EarthFlux(shortwave)_Z(w/m^2),";
        
        res += "refpos_X,";
        res += "refpos_Y,";
        res += "refpos_Z,";
        res += "refvel_X,";
        res += "refvel_Y,";
        res += "refvel_Z,";
        
        //res.strip_v(',');
        res += "\n";
        return res;
    }
    
    
   // add more information to the m_logStr
    void GOrbitPredictor::addStateInformation(GVector& refpos, GVector& refvel)
    {
        if(log_on == true)
        {
            m_logStr = m_logStr.substr(0, m_logStr.length() - 1);
            m_logStr +=
            ( GString(refpos.x*1000.0,6) + ","
             + GString(refpos.y*1000.0,6) + ","
             + GString(refpos.z*1000.0,6) + ","
             + GString(refvel.x*1000.0,6) + ","
             + GString(refvel.y*1000.0,6) + ","
             + GString(refvel.z*1000.0,6) + ","
             );
            
            m_logStr+="\n";
        }
    }
    
    void GOrbitPredictor::outputLog()
    {
       m_orbitLog << m_logStr;
    }
    
    void GOrbitPredictor::collectStateInformation()
    {
        if(log_on == true)
        {
            double R2D = GCONST("R2D");
            
            m_logStr = "";
            
            m_logStr += (m_spaceCraft->getStatePointer()->ct.TimeString()+",");
            
            m_logStr +=
            ( GString(m_spaceCraft->getStatePointer()->satpos_eci.x*1000.0,6) + ","
             + GString(m_spaceCraft->getStatePointer()->satpos_eci.y*1000.0,6) + ","
             + GString(m_spaceCraft->getStatePointer()->satpos_eci.z*1000.0,6) + ","
             + GString(m_spaceCraft->getStatePointer()->satvel_eci.x*1000.0,6) + ","
             + GString(m_spaceCraft->getStatePointer()->satvel_eci.y*1000.0,6) + ","
             + GString(m_spaceCraft->getStatePointer()->satvel_eci.z*1000.0,6) + ","
             
             + GString(m_spaceCraft->getStatePointer()->eps*R2D) + ","
             
             + GString(m_spaceCraft->getStatePointer()->attitude_eci.beta*R2D) + ","
             
             + GString(m_spaceCraft->getStatePointer()->attitude_eci.eta*R2D) + ","
             
             + GString(m_spaceCraft->getStatePointer()->attitude_eci.phi*R2D) + ","
             
             + GString(m_spaceCraft->getStatePointer()->keplerianElement.m_sma) + ","
             + GString(m_spaceCraft->getStatePointer()->keplerianElement.m_ecc) + ","
             + GString(m_spaceCraft->getStatePointer()->keplerianElement.m_inc*R2D) + ","
             + GString(m_spaceCraft->getStatePointer()->keplerianElement.m_argp*R2D) + ","
             + GString(m_spaceCraft->getStatePointer()->keplerianElement.m_raan*R2D) + ","
             + GString(m_spaceCraft->getStatePointer()->keplerianElement.m_tran*R2D) + ","
             );
            
            std::map< GString, GForceModel* >::iterator   it = m_forceManager.m_forceModels.begin();
            for( ;it != m_forceManager.m_forceModels.end(); ++it )
            {
                
                GVector force =  it->second->getForce();
                
                //GString forceName = it->second->getForceName();
                m_logStr += ( GString(force.x) + "," + GString(force.y) + "," + GString(force.z) + ","  );
            }
            
            // m_logStr.stripTrailing_v(',');
            
            GVector totalSolarFlux = m_spaceCraft->getStatePointer()->solarFlux.m_flux.m_dir*
            (m_spaceCraft->getStatePointer()->solarFlux.m_flux.m_longwave
             +m_spaceCraft->getStatePointer()->solarFlux.m_flux.m_shortwave);
            
            m_logStr += GString(totalSolarFlux.x) +"," + GString(totalSolarFlux.y) +"," + GString(totalSolarFlux.z) + ",";
            
            m_logStr +=  GString(m_spaceCraft->getStatePointer()->earthFlux.totalFlux_lw.x) +","
            + GString(m_spaceCraft->getStatePointer()->earthFlux.totalFlux_lw.y) +","
            + GString(m_spaceCraft->getStatePointer()->earthFlux.totalFlux_lw.z) +","
            + GString(m_spaceCraft->getStatePointer()->earthFlux.totalFlux_sw.x) +","
            + GString(m_spaceCraft->getStatePointer()->earthFlux.totalFlux_sw.y) +","
            + GString(m_spaceCraft->getStatePointer()->earthFlux.totalFlux_sw.z) +"," ;
            
            //m_logStr.strip_v(',');
            m_logStr += "\n";
            
            
        }
        
    }
    
       
    void GOrbitPredictor::PropagateTo( GTime t )
    {
        
        m_t0 = m_spaceCraft->getStatePointer()->m_epoch;
        
        int ndim = 6;
        
        gfc::TimeSystem ts;
        
        long mjd = 0,sod = 0;
        
        double fsod = 0.0;
        
        (t - m_t0).GetData( ts, mjd, sod, fsod );
        
        double sec_end = mjd*GCONST("SECPDAY") + sod + fsod;
        
        GVector p = m_spaceCraft->getStatePointer()->satpos_eci;
        GVector v = m_spaceCraft->getStatePointer()->satvel_eci;
        double y0[6] = {p.x,p.y,p.z,v.x,v.y,v.z};
        
        double yend[6] = {0.0};
        
        
        m_integrator->IntegrateTo( this , ndim, 0.0, y0, sec_end, yend);
        
        //GKeplerianElements::propagate(m_spaceCraft->getStatePointer()->keplerianElement, <#double t#>, <#gfc::GVector &p#>, <#gfc::GVector &v#>)
        
        // considering collect state information here, not in the RungeKutta intergrator
        collectStateInformation();
        
        p.x = yend[0];p.y = yend[1];p.z = yend[2];
        v.x = yend[3];v.y = yend[4];v.z = yend[5];
        
        m_spaceCraft->getStatePointer()->updateState_eci(t, p, v);
        
        int testc = 0;
        
        
        
        
    }
    
    
} // end of namespace

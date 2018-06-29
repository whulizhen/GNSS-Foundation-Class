//
//  GOrbitPredictor.hpp
//  GFC
//
//  Created by lizhen on 21/05/2016.
//  Copyright © 2016 lizhen. All rights reserved.
//

#ifndef GOrbitPredictor_hpp
#define GOrbitPredictor_hpp

#include <stdio.h>
#include "GRungeKutta.hpp"
#include "GAdams.hpp"

#include "GMatrix.h"
#include "GTime.h"
#include "GForceModel.hpp"
#include "GSpacecraft.hpp"
#include "GFMGravity.hpp"
#include "GFMEarthRadiationPressure.hpp"
#include "GFMSolarRadiationPressure.hpp"
#include "GFMThirdBody.hpp"
#include "GFMEarthTide.hpp"
#include "GMotionDynamic.hpp"

#include "GForceModelMgr.hpp"

namespace gfc
{
    
    // this external function is used as the friend of the class GOrbitPredictor
    //void getDerivatives( int n, double x, double *y, double *dydx);
    
    class GOrbitPredictor : public GMotionDynamic
    {
        
    public:
        
        void getDerivatives( int n, double x, double *y, double *dydx);
        
        GOrbitPredictor( GSpaceCraft* spacecraft);
        
        virtual ~GOrbitPredictor();
        
        void setLogfile(GString logfilename);
        
        void setStepsize( double stepsize);
        
        void setSpaceCraft( GSpaceCraft* spacecraft);
        
        void setForceManager( GForceModelMgr& manager);
        
        GForceModelMgr getForceManager();
        
        
        //integrate from t0 (with the state) to t
        void PropagateTo( GTime t );
        
        //return the orbit diff in R T N directions
        static GVector orbitdiff( GVector& p1, GVector& v1, GVector& p2, GVector& v2);
        
        GString getLogHeader();
        //collect the state information
        void collectStateInformation();
        void addStateInformation(GVector& refpos, GVector& refvel);
        void outputLog();
        
    private:
        
        GTime m_t0;
        
        GString m_logStr; // the string for storing the current log record
        GString m_logfilename; // the logging file name
        ofstream m_orbitLog; // the log file for the whole orbit prediction process
        bool log_on;
        GForceModelMgr m_forceManager;
        
        // a string list to store all the force names acting on this spacecraft
        //std::map< GString, GForceModel* > m_forceModels;
        
        /// current state
        // r        3
        // v        3
        // dr_dr0   3*3
        // dr_dv0   3*3
        // dr_dp0   3*np
        // dv_dr0   3*3
        // dv_dv0   3*3
        // dv_dp0   3*np
        
        // the GRungeKuttaFehlbert should be replaced later by the father class of integrator
        GRungeKuttaFehlberg* m_integrator;
        
        // Adams Bashforth and Moulton intergrator
        //GAdams* m_integrator;
        
        GSpaceCraft* m_spaceCraft; // it should be a pointer to a static variable
        
        //GTime m_epoch;
        
        //GMatrix m_state;  // the current state
        
        //GMatrix m_phi;  // the state transmission matrix
        
        //GMatrix m_dphi; // the first derivative of phi, obtained from variation equation
        
        /// the sensitivity matrix
        //GMatrix m_s;         // 6*np
        
    };
    
    
    
    
    
}



#endif /* GOrbitPredictor_hpp */

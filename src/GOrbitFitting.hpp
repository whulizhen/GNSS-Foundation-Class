
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
//  GInitialOrbit.hpp
//  GFC
//
//  Created by lizhen on 10/07/2016.
//  Copyright © 2016 lizhen. All rights reserved.
//
// this class is to determine the initial value using sp3 as observation
// basically, this class fit the sp3 orbit with dynamic force models

#ifndef GInitialOrbit_hpp
#define GInitialOrbit_hpp

#include <stdio.h>
#include <map>
#include "GString.h"
#include "GForceModel.hpp"
#include "GRungeKutta.hpp"
#include "GSpacecraft.hpp"
#include "GMotionDynamic.hpp"
#include "GForceModelMgr.hpp"

namespace gfc
{

class GOrbitFitting : public GMotionDynamic
{
    
public:
    
    GOrbitFitting();
    
    GOrbitFitting(GSpaceCraft* spacecraft);
    
    void setInitialValue(GTime t, GVector p, GVector v);
    
    void setEndEpoch(GTime& t);
    void setOutputFile(GString output)
    {
        m_outputFile = output;
    }
    
    GVector getPosition() {return m_initialP;}
    GVector getVelocity() {return m_initialV;}
    double getRMS() {return m_fittingRMS;}
    void setSpaceCraft( GSpaceCraft* spacecraft);
    
    GForceModelMgr getForceManager();
    void setForceManager(GForceModelMgr manager);
    
     //this should include state vector and the variation equation
     // here, n may be 6 + 6*6 + 6*np, np is the number of force model parameters
    void getDerivatives( int n, double x, double *y, double *dydx );
    
    //integrate from t0 (with the state) to t
    void PropagateTo( GTime t );
    
    // the main function to calculate the initial value error
    void calculateInitialError();
    void calculate_SRPModel();
    
    // the main function to filter the corresponding parameters
    void fitting();
    
    void setStepSize(double stepsize);
    
private:
    
    GVector m_initialP;  // initial position
    
    GVector m_initialV;  // initial velocity
    
    GTime   m_initialEpoch; // the time for the initial value
    
    GTime   m_endEpoch;  // the end time for the orbit fitting
    
    GString m_outputFile;  // the output filename
    
    GTime m_t0; // temp variable
    
    double m_fittingRMS;
    
    GForceModelMgr m_forceManager;
      // the GRungeKuttaFehlbert should be replaced later by the father class of integrator
    GRungeKuttaFehlberg* m_integrator;
    
    GSpaceCraft* m_spaceCraft; // it should be a pointer to a static variable, spacecraft include the motion state
    
};

    
} // end of namespace


#endif /* GInitialOrbit_hpp */

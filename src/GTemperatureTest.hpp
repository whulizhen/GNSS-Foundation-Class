//
//  GTemperatureTest.hpp
//  GFC
//
//  Created by lizhen on 25/08/2016.
//  Copyright © 2016 lizhen. All rights reserved.
//

#ifndef GTemperatureTest_hpp
#define GTemperatureTest_hpp

#include <stdio.h>

#include "GMotionDynamic.hpp"

#include "GRungeKutta.hpp"

namespace gfc
{
    class GTemperatureTest : public GMotionDynamic
    {
        
    public:
        
        GTemperatureTest()
        {
            m_integrator.setStepsize(30);
        }
        
        //considering the
        void getDerivatives( int n, double x, double *y, double *dydx )
        {
            //dt1
            dydx[0] = 1.0/500*( 1550.0 + 3.8E-3*(y[1]-y[0]) + 12.8E-9*(pow(y[1],4.0) -pow(y[0],4.0) ) - 51.3E-9*pow(y[0],4.0) - 64.1E-9*pow(y[0],4.0) );
            dydx[1] = 1.0/15000.0*( 3.8E-3*(y[1]-y[0]) - 12.8E-9*( pow(y[1],4) - pow(y[0],4.0) ) - 166E-9*pow(y[1],4.0) );
            
            printf("y: %f %f\n",y[0],y[1]);
            //printf("dydx: %f %f\n",dydx[0],dydx[1]);
        }
        
        // to solve the temperature with time
        void PropagateTo(double t)
        {
            double y0[2] = {300,300};
            double yend[2]  = {0.0};
            int ndim = 2;
            t= 30000;
            m_integrator.IntegrateTo( this , 2, 0.0, y0, t, yend);
            
            printf("yend: %f %f\n", yend[0],yend[1]);
        }
        
        
        double dT1()
        {
            double dt1 = 0.0;
            
            
            
            return dt1;
        }
        
        double dT2()
        {
            double dt2 = 0.0;
            
            
            
            return dt2;
        }
        
        
    private:
        
        GRungeKuttaFehlberg m_integrator;
        
    };
    
    
    
} // end of namespace




#endif /* GTemperatureTest_hpp */

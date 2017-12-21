//
//  GAdamsCowell.hpp
//  GFC
//
//  Created by lizhen on 17/12/2016.
//  Copyright © 2016 lizhen. All rights reserved.
//

#ifndef GAdamsCowell_hpp
#define GAdamsCowell_hpp
#include <vector>
#include <queue>
#include <stdio.h>
#include "GRungeKutta.hpp"
#include "GMotionDynamic.hpp"
#include "GMatrix.h"

namespace  gfc
{
    
    // class for multi-step intergration
    // Adams-Bashforth (Predictor) and Adams-Moulton (Corrector)
    /*
     reference:
     http://drum.lib.umd.edu/bitstream/handle/1903/2202/2004-berry-healy-jas.pdf?sequence=7&isAllowed=y
     http://www.sml.ece.upatras.gr/UploadedFiles/BOOK-CK/04-MultistepIntegrationMethods.pdf
     
     Bashforth:
     Yn+1 = Yn + h*( a1*Y'n + a2*Y'n-1 + a3*Y'n-2 + ... + ak*Y'n-k+1 )
     
     Moulton:
     Yn+1 = Yn + h*( b1*Y'n+1 + b2*Y'n + b3*Y'n-1 + ... + bk*Y'n-k+2 )
     
    */
    class GAdams
    {
        
    public:
       
       GAdams(int order = 9 );
       
       void setStepsize(double stepsize);
        
       void singleStep( GMotionDynamic* porb, int n, double x0, double *y0, double h, double *y);
       
       void IntegrateTo(GMotionDynamic* porb,  int n , double start, double* ystart, double end,double* yend );
        
        //clear the startup information, when change the intergrate directions
       void clearStartup();
        
    private:
        
        //void getBashforthCoeff();
        //void getMoultonCoeff();
        
        void getCoef();
        
        // RFK as a start up intergrator
        GRungeKuttaFehlberg* m_rkf;
        double m_stepsize;
        int m_order; // the order of adams-cowell method
        std::vector<double> m_bC;  //bashforth coefficients
        std::vector<double> m_mC;  // moulton coefficients
        
        std::deque<GMatrix>  m_startup; // the start up , including
        
    }; // end of class Adams-Cowell
    
    
}  // end of namespace








#endif /* GAdamsCowell_hpp */

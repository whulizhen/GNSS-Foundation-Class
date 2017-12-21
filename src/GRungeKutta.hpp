//
//  GRungeKutta.hpp
//  GFC
//
//  Created by lizhen on 16/3/9.
//  Copyright © 2016年 lizhen. All rights reserved.
//


//============================================================================
//
//  This file is part of GFC, the GNSS FOUNDATION CLASS.
//
//  The GFC is free software; you can redistribute it and/or modify
//  it under the terms of the GNU Lesser General Public License as published
//  by the Free Software Foundation; either version 2.1 of the License, or
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


#ifndef GRungeKutta_hpp
#define GRungeKutta_hpp

#include <stdio.h>
#include <vector>
#include <math.h>
#include "GMotionDynamic.hpp"

using namespace std;

namespace gfc
{
    
    // RKF7(8)
    class GRungeKuttaFehlberg
    {
        
    public:
        //定义一个函数指针类型，类型名为DYDX
        //typedef void (*pDYDX)(int, double , double* , double* );
        
        GRungeKuttaFehlberg();
        virtual ~GRungeKuttaFehlberg();
        
        GRungeKuttaFehlberg& operator= (const GRungeKuttaFehlberg& right); // 赋值重载
        GRungeKuttaFehlberg( const GRungeKuttaFehlberg& right );   //copy construction function
        
        
        //void IntegrateTo( GOrbitPredictor* porb,  int n , double start, double* ystart, double end,double* yend );
        void IntegrateTo( GMotionDynamic* porb,  int n , double start, double* ystart, double end,double* yend );
        
        
        void setStepsize( double stepsize);
        
        //int  singleStep( GOrbitPredictor* porb, int n, double x0, double *y0,double& h, double *y, double *yerr);
        int  singleStep( GMotionDynamic* porb, int n, double x0, double *y0,double& h, double *y, double *yerr);
        
//        void setDerivativeFunction(GOrbitPredictor* porb)
//        {
//            //m_porb = porb;
//        }
        
        /*the test function providing the acceleration and the */
        /*
         *  in the runtime, use member setDerivativeFunction to provide this function
            this one is just for test.
         */
        void static testacceleration(int n, double x0, double* y0, double* dydx )
        {
            //int n, double x, double* y, double* dydx
            if( fabs(y0[0])< GRungeKuttaFehlberg::EPS )
            {
                y0[0] = GRungeKuttaFehlberg::EPS;
            }
            
            dydx[0] = 4.0*x0+3.0;//y0[0] - 2.0*x0/y0[0];
        }
        
        
    private:
        
        //a test function
        //void getDerivatives(GOrbitPredictor* porb, int n, double x, double* y, double* dydx );
        
        double A(int i)  {return m_a[i];}
        double B(int i,int j) { return m_b[i*12+j]; }
        double C1(int i)   {return m_c1[i];}
        double C2(int i)  {return m_c2[i];}
        double DC(int i)  {return (m_c1[i] - m_c2[i]);}
        
        static const double EPS;  // to control the precision
        
        void getCoefficients();
        
        //函数指针，指向了函数dydx,其函数原型为 void dydx(double, double*)
        //GOrbitPredictor*  m_porb;
        
        bool m_isAdaptive;
        
        int m_order;  // the order of the integrator
        
        double m_stepLen;   // the stepsize of the integretor
        //int m_num_coef;  // the number of coefficient
        //int m_bCol ;
        //int m_bRow ;
        
        double m_a[13];
        double m_b[13*12];
        double m_c1[13];
        double m_c2[13];
        
        
        /*
        double m_a[6];
        double m_b[6*5];
        double m_c1[6];
        double m_c2[6];
        */
        
        
    };


    
}  // end of namespace gfc




#endif /* GRungeKutta_hpp */

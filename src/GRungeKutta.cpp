//
//  GRungeKutta.cpp
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


#include "GRungeKutta.hpp"
#include "GOrbitPredictor.hpp"

namespace gfc
{
    
    const double GRungeKuttaFehlberg::EPS  = 1.0-30;
    
    
    GRungeKuttaFehlberg::GRungeKuttaFehlberg()
    {
        m_order = 13;
        m_isAdaptive = false;
        //determine the coefficients
        getCoefficients();
        
    }
    
    
    GRungeKuttaFehlberg& GRungeKuttaFehlberg::operator= (const GRungeKuttaFehlberg& right) // 赋值重载
    {
        m_isAdaptive = right.m_isAdaptive;
        m_order = right.m_order;
        m_stepLen = right.m_stepLen;
        
        getCoefficients();
        
        return *this;
    }
    
    GRungeKuttaFehlberg::GRungeKuttaFehlberg( const GRungeKuttaFehlberg& right )   //copy construction function
    {
        m_isAdaptive = right.m_isAdaptive;
        m_order = right.m_order;
        m_stepLen = right.m_stepLen;
        
        getCoefficients();
        
    }
    
    
    
    GRungeKuttaFehlberg::~GRungeKuttaFehlberg()
    {
        
    }
    
    /*private function to get the coefficients for the certain order
      it is said that is very difficult to do so
     */
    void GRungeKuttaFehlberg::getCoefficients( )
    {
        
        
        double a[13]={0.0, 2.0 / 27.0, 1.0 / 9.0, 1.0 / 6.0, 5.0 / 12.0, 1.0 / 2.0, 5.0 / 6.0, 1.0 / 6.0, 2.0 / 3.0, 1.0 / 3.0,1.0,0.0,1.0};
        double b[13*12] = {
            0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0     ,
            2.0/27.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0 ,
            1.0/36.0,1.0/12.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0 ,
            1.0/24.0,0.0,1.0/8.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0 ,
            5.0/12.0,0.0,-25.0/16.0,25.0/16.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
            1.0/20.0,0.0,0.0,1.0/4.0,1.0/5.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0 ,
            -25.0/108.0,0.0,0.0,125.0/108.0,-65.0/27.0,125.0/54.0,0.0,0.0,0.0,0.0,0.0,0.0,
            31.0/300.0,0.0,0.0,0.0,61.0/225.0,-2.0/9.0,13.0/900.0,0.0,0.0,0.0,0.0,0.0 ,
            2.0,0.0,0.0,-53.0/6.0,704.0/45.0,-107.0/9.0,67.0/90.0,3.0,0.0,0.0,0.0,0.0 ,
            -91.0/108.0,0.0,0.0,23.0/108.0,-976.0/135.0,311.0/54.0,-19.0/60.0,17.0/6.0,-1.0/12.0,0.0,0.0,0.0,
            2383.0/4100.0,0.0,0.0,-341.0/164.0,4496.0/1025.0,-301.0/82.0,2133.0/4100,45.0/82.0,45.0/164.0,18.0/41.0,0.0,0.0,
            3.0/205.0,0.0,0.0,0.0,0.0,-6.0/41.0,-3.0/205.0,-3.0/41.0,3.0/41.0,6.0/41.0,0.0,0.0,
            -1777.0/4100.0,0.0,0.0,-341.0/164.0,4496.0/1025.0,-289.0/82.0,2193.0/4100.0,51.0/82.0,33.0/164.0,12.0/41.0,0.0,1.0
        };
        double c1[13]={41.0/840,0.0,0.0,0.0,0.0,34.0/105.0,9.0/35.0,9.0/35.0,9.0/280.0,9.0/280.0,41.0/840.0,0.0,0.0};
        double c2[13]={ 0.0,0.0,0.0,0.0,0.0,34.0/105.0,9.0/35.0,9.0/35.0,9.0/280.0,9.0/280.0,0.0,41.0/840.0,41.0/840.0};
        
         memcpy(m_a,a,sizeof(double)*13);
         memcpy(m_b,b,sizeof(double)*13*12);
         memcpy(m_c1, c1, sizeof(double)*13);
         memcpy(m_c2,c2,sizeof(double)*13);
         
         
        
        /*
        double a[6] = {0.0, 1.0/4.0, 3.0/8.0, 12.0/13.0, 1.0, 1.0/2.0 };
        double b[6*5] = {  0.0, 0.0, 0.0 ,0.0, 0.0,
                           1.0/4.0, 0.0, 0.0 ,0.0, 0.0,
                           3.0/32.0, 9.0/32.0, 0.0 ,0.0, 0.0,
                           1932.0/2197.0, -7200.0/2197.0, 7296.0/2197.0, 0.0, 0.0 ,
                           439.0/216.0, -8.0, 3680.0/513.0, -845.0/4104.0, 0.0,
                           -8.0/27.0, 2.0, -3544.0/2565.0, 1859.0/4104.0, -11.0/40.0
                        };
        
        double c1[6] = { 25.0/216.0, 0.0, 1408.0/2565.0, 2197.0/4101.0, -1.0/5.0, 0.0 };
        double c2[6] = { 16.0/135.0, 0.0, 6656.0/12825.0, 28561.0/56430.0, -9.0/50.0, 2.0/55.0 };
   
        memcpy(m_a,a,sizeof(double)*6);
        memcpy(m_b,b,sizeof(double)*6*5);
        memcpy(m_c1, c1, sizeof(double)*6);
        memcpy(m_c2,c2,sizeof(double)*6);
        */
    }
    
    void GRungeKuttaFehlberg::setStepsize(double stepsize)
    {
        m_stepLen = stepsize;
    }
    
    
    /* single step,  the step can be adapted acoording to yerr 
         n: number of variable
         x0:  the variable, usually, it is the time
         y0:   the initail value
         y:    the integrated value
         h: stepsize
     */
    
    //int GRungeKuttaFehlberg::singleStep(GOrbitPredictor* porb, int n, double x0, double *y0,double& h, double *y, double *yerr)
    int GRungeKuttaFehlberg::singleStep( GMotionDynamic* porb, int n, double x0, double *y0,double& h, double *y, double *yerr)
    {
        
        double* dydx = new double[n];
        memset(dydx,0,sizeof(double)*n);
        
        double** K = new double*[m_order];
        for( int i = 0;i<m_order; i++ )
        {
            K[i] = new double[n];
            memset(K[i],0,sizeof(double)*n);
        }
        
        double* yk = new double[n];
        memcpy(yk,y0,sizeof(double)*n);
        
        
        //this is for index ==0
        porb->getDerivatives(n, x0, yk, dydx);
        for( int i = 0; i< n; i++ )
        {
            K[0][i] = dydx[i]*h;
        }
        
        double coef = 0.0;
        for( int index = 0 ; index < m_order; index++ )
        {
            for(int i = 0 ; i< n; i++ )
            {
                coef = 0.0;
                for(int j = 0; j< index; j++ )
                {
                    coef += m_b[index*(m_order-1)+j]*K[j][i];
                }
                yk[i] = y0[i] + coef;
            }
            
            //
            porb->getDerivatives(n, x0+m_a[index]*h, yk, dydx);
            
            for( int i = 0; i< n; i++ )
            {
                K[index][i] = dydx[i]*h;
            }
            
        }
        
        
        for(int i = 0 ; i< n; i++ )
        {
            double coef =0.0;
            for( int index = 0; index < m_order; index++ )
            {
                coef += m_c2[index]*K[index][i];
            }
            //this y[i] is for error check
            //获取最后的y[i]
            
            y[i] = y0[i] + coef;
        }
        
        //update all the derivatives finally, including spacecraft state,
        //porb->getDerivatives(n, x0+h, y, dydx);
        
        
        if(yk != NULL ) { delete[] yk; yk = NULL; }
        if(dydx != NULL ) { delete[] dydx; dydx = NULL; }
        
        for(int i = 0 ; i< m_order; i++ )
        {
            if(K[i]!= NULL ) {delete[] K[i];K[i]=NULL;}
        }
        if(K != NULL) {delete[] K; K = NULL;}

        
//        
//        ///DYDX is m_bRow row and n colomn
//        double** dydx = new double*[13];
//        for( int i = 0;i<13; i++ )
//        {
//            dydx[i] = new double[n];
//            memset(dydx[i],0,sizeof(double)*n);
//        }
//        
//        
//        double* ytmp = new double[n];
//        memcpy(ytmp, y0, sizeof(double)*n);
//        
//       // porb->getDerivatives(n, x0, ytmp, dydx[0]);
//        
//        for( int i = 0  ;i< 13 ; i++ )
//        {
//            for( int k = 0; k< n; k++ )  // all the dimentions
//            {
//                double coff = 0.0;
//                for( int j = 0 ; j< i; j++ )
//                {
//                    coff += m_b[i*12+j]*dydx[j][k];
//                }
//                
//                ytmp[k] = y0[k] + h*coff;
//                //ytmp[k] = y0[k] + coff;
//            }
//            
//            // get the derivatives
//            porb->getDerivatives(n, x0+m_a[i]*h, ytmp, dydx[i]);
//            
//        }
//        
//        
//        for( int i = 0 ; i< n ; i++ )  //// all the dimentions
//        {
//            double coff = 0.0;
//            for( int j = 0 ; j< 13; j++ )
//            {
//                coff = coff + m_c2[j] * dydx[j][i] ;
//                //coff = coff + h* (m_c2[j] * dydx[j][i]) ;
//                
//            }
//            
//            y[i] =  y0[i] + h*coff;
//            
//            yerr[i] = h*m_c1[0]*(dydx[11][i] + dydx[12][i] - dydx[0][i] - dydx[10][i]);
//        }
//        
//        if(ytmp != NULL ) { delete[] ytmp; ytmp = NULL; }
//        for(int i = 0 ; i< 13; i++ )
//        {
//            if(dydx[i]!= NULL ) {delete[] dydx[i];dydx[i]=NULL;}
//        }
//        if(dydx != NULL) {delete[] dydx; dydx = NULL;}
//        
        
        /*
        double* dydx = new double[n];
        porb->getDerivatives(n, x0, y0, dydx); //ak1
        
        
        double* ytemp = new double[n];
        for( int i = 0 ; i< n; i++ )
        {
            ytemp[i] = y0[i] + h*B(1,0)*dydx[i];
        }
        //ak2
        double* ak2 = new double[n];
        porb->getDerivatives(n,x0+A(1)*h, ytemp, ak2);
        
        for( int i = 0 ; i< n; i++ )
        {
            ytemp[i] = y0[i] + h*(B(2,0)*dydx[i] + B(2,1)*ak2[i] );
        }
        //ak3
        double* ak3 = new double[n];
        porb->getDerivatives(n, x0+A(2)*h, ytemp, ak3);
        
        for( int i = 0 ; i< n; i++ )
        {
            ytemp[i] = y0[i] + h*(B(3,0)*dydx[i] + B(3,1)*ak2[i] + B(3,2)*ak3[i] );
        }
        //ak4
        double* ak4 = new double[n];
        porb->getDerivatives(n, x0+A(3)*h,ytemp, ak4);
        
        for( int i = 0 ; i< n; i++ )
        {
            ytemp[i] = y0[i] + h*(B(4,0)*dydx[i] + B(4,1)*ak2[i] + B(4,2)*ak3[i] + B(4,3)*ak4[i] );
        }
        //ak5
        double* ak5 = new double[n];
        porb->getDerivatives(n, x0+A(4)*h,ytemp, ak5);
        
        for( int i = 0 ; i< n; i++ )
        {
            ytemp[i] = y0[i] + h*(B(5,0)*dydx[i] + B(5,1)*ak2[i] + B(5,2)*ak3[i] + B(5,3)*ak4[i]+B(5,4)*ak5[i] );
        }
        //ak6
        double* ak6 = new double[n];
        porb->getDerivatives(n, x0+A(5)*h, ytemp, ak6);
        
        for( int i = 0 ; i< n; i++ )
        {
            ytemp[i] = y0[i] + h*(B(6,0)*dydx[i] + B(6,1)*ak2[i] + B(6,2)*ak3[i] + B(6,3)*ak4[i]+B(6,4)*ak5[i]+B(6,5)*ak6[i] );
        }
        //ak7
        double* ak7 = new double[n];
        porb->getDerivatives(n, x0+A(6)*h, ytemp, ak7);
        
        for( int i = 0 ; i< n; i++ )
        {
            ytemp[i] = y0[i] + h*(B(7,0)*dydx[i] + B(7,1)*ak2[i] + B(7,2)*ak3[i] +
                                         B(7,3)*ak4[i]+B(7,4)*ak5[i]+B(7,5)*ak6[i]+B(7,6)*ak7[i] );
        }
        //ak8
        double* ak8 = new double[n];
        porb->getDerivatives(n, x0+A(7)*h, ytemp, ak8 );
        
        for( int i = 0 ; i< n; i++ )
        {
            ytemp[i] = y0[i] + h*(B(8,0)*dydx[i] + B(8,1)*ak2[i] + B(8,2)*ak3[i]
                                         +B(8,3)*ak4[i]+B(8,4)*ak5[i]+B(8,5)*ak6[i]+B(8,6)*ak7[i]
                                         +B(8,7)*ak8[i]);
        }
        //ak9
        double* ak9 = new double[n];
        porb->getDerivatives(n, x0+A(8)*h, ytemp, ak9 );
        
        for( int i = 0 ; i< n; i++ )
        {
            ytemp[i] = y0[i] + h*(B(9,0)*dydx[i] + B(9,1)*ak2[i] + B(9,2)*ak3[i]
                                         +B(9,3)*ak4[i]+B(9,4)*ak5[i]+B(9,5)*ak6[i]+B(9,6)*ak7[i]
                                         +B(9,7)*ak8[i]+B(9,8)*ak9[i]);
        }
        //ak10
        double* ak10 = new double[n];
        porb->getDerivatives(n, x0+A(9)*h, ytemp, ak10 );
        
        for( int i = 0 ; i< n; i++ )
        {
            ytemp[i] = y0[i] + h*(B(10,0)*dydx[i] + B(10,1)*ak2[i] + B(10,2)*ak3[i]
                                         +B(10,3)*ak4[i]+B(10,4)*ak5[i]+B(10,5)*ak6[i]+B(10,6)*ak7[i]
                                         +B(10,7)*ak8[i]+B(10,8)*ak9[i]+B(10,9)*ak10[i]);
        }
        //ak11
        double* ak11 = new double[n];
        porb->getDerivatives(n, x0+A(10)*h, ytemp, ak11 );
        
        
        for( int i = 0 ; i< n; i++ )
        {
            ytemp[i] = y0[i] + h*(B(11,0)*dydx[i] + B(11,1)*ak2[i] + B(11,2)*ak3[i]
                                         +B(11,3)*ak4[i]+B(11,4)*ak5[i]+B(11,5)*ak6[i]+B(11,6)*ak7[i]
                                         +B(11,7)*ak8[i]+B(11,8)*ak9[i]+B(11,9)*ak10[i]+B(11,10)*ak11[i]);
        }
        //ak12
        double* ak12 = new double[n];
        porb->getDerivatives(n, x0+A(11)*h, ytemp, ak12 );
        
        
        for( int i = 0 ; i< n; i++ )
        {
            ytemp[i] = y0[i] + h*(B(12,0)*dydx[i] + B(12,1)*ak2[i] + B(12,2)*ak3[i]
                                         +B(12,3)*ak4[i]+B(12,4)*ak5[i]+B(12,5)*ak6[i]+B(12,6)*ak7[i]
                                         +B(12,7)*ak8[i]+B(12,8)*ak9[i]+B(12,9)*ak10[i]+B(12,10)*ak11[i]
                                         +B(12,11)*ak12[i]);
        }
        //ak13
        double* ak13 = new double[n];
        porb->getDerivatives(n, x0+A(12)*h, ytemp, ak13 );
        
        memset(y,0,sizeof(double)*n);
        memset(yerr,0,sizeof(double)*n);
        
        for(int i = 0 ; i< n ; i++ )
        {
            //the 8th
            y[i] =  y0[i] + h * (C2(0) * dydx[i] + C2(1) * ak2[i] + C2(2) * ak3[i] + C2(3) * ak4[i]
                                + C2(4) * ak5[i] + C2(5) * ak6[i] + C2(6) * ak7[i] + C2(7) * ak8[i] + C2(8) * ak9[i]
                                + C2(9) * ak10[i] + C2(10) * ak11[i] + C2(11) * ak12[i] + C2(12) * ak13[i]) ;
        
        
            yerr[i] = h*C1(0)*(ak12[i] + ak13[i] - dydx[i] - ak11[i]);
        }
        
        
        if(dydx != NULL) {delete[] dydx, dydx = NULL;}
        if(ytemp != NULL) {delete[] ytemp, ytemp = NULL;}
        if(ak2 != NULL) {delete[] ak2, ak2 = NULL;}
        if(ak3 != NULL) {delete[] ak3, ak3 = NULL;}
        if(ak4 != NULL) {delete[] ak4, ak4 = NULL;}
        if(ak5 != NULL) {delete[] ak5, ak5 = NULL;}
        if(ak6 != NULL) {delete[] ak6, ak6 = NULL;}
        if(ak7 != NULL) {delete[] ak7, ak7 = NULL;}
        if(ak8 != NULL) {delete[] ak8, ak8 = NULL;}
        if(ak9 != NULL) {delete[] ak9, ak9 = NULL;}
        if(ak10 != NULL) {delete[] ak10, ak10 = NULL;}
        if(ak11 != NULL) {delete[] ak11, ak11 = NULL;}
        if(ak12 != NULL) {delete[] ak12, ak12 = NULL;}
        if(ak13 != NULL) {delete[] ak13, ak13 = NULL;}
         
         */
//
//
        
        //here collect the state and acceleration information
       // porb->collectStateInformation();
        
        return 0;
    }
    
    
    //void GRungeKuttaFehlberg::IntegrateTo(GOrbitPredictor* porb,  int n , double start, double* ystart, double end,double* yend )
    void GRungeKuttaFehlberg::IntegrateTo(GMotionDynamic* porb,  int n , double start, double* ystart, double end,double* yend )
    {
        bool rec = false;
        int sign = -1; // determine the direction of intergration
        
        if( end >= start)
        {
            sign = 1;
        }
        
        double* yerr = new double[n];
        double* yold = new double[n];
        memcpy(yold,ystart,sizeof(double)*n);
        memcpy(yend,ystart,sizeof(double)*n);
        memset(yerr,0,sizeof(double)*n);
        
        double dt = m_stepLen*sign;
        double tt = start;
        
        /*
        while( tt<end )
        {
            if( (tt + dt) >= end ) // here must include the equal sign
            {
                dt = end - tt;
                break;
            }
            
            singleStep(porb,n, tt, yold, dt,yend, yerr);
            
            memcpy(yold, yend, sizeof(double)*n);
            tt+=dt;
        }
        
        //process the last step, in case of the last step is not a fully step.
        //dt = end - tt;
        if( fabs(dt) > 1.0E-10 )
        {
            singleStep(porb,n, tt, yold,dt, yend, yerr);
        }
        */
        
        
        
        while( (end - tt)*sign > 0.0 )
        {
            if( ( (tt + m_stepLen*sign) - end )*sign >= 0.0 ) // here must include the equal sign
            {
                dt = end - tt;
                rec = true;
                break;
            }
            
            singleStep( porb,n, tt, yold, dt ,yend, yerr);
            memcpy(yold, yend, sizeof(double)*n);
            tt+= m_stepLen*sign;
        }
        
         if(rec == true && fabs(dt) > 1.0E-10)
         {
             singleStep(porb,n, tt, yold,dt, yend, yerr);
         }
        
        
        
        
        if(yerr != NULL )
        {
            delete[] yerr;
            yerr = NULL;
        }
        if(yold != NULL )
        {
            delete[] yold;
            yold = NULL;
        }
       
    }
    
    
    
}





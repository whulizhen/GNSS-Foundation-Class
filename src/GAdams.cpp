//
//  GAdamsCowell.cpp
//  GFC
//
//  Created by lizhen on 17/12/2016.
//  Copyright © 2016 lizhen. All rights reserved.
//

#include "GAdams.hpp"
#include "GMath.hpp"
namespace gfc
{
    
    
    GAdams::GAdams(int order)
    {
        m_order = order;
        m_bC.resize(m_order);
        m_mC.resize(m_order);
        
        getCoef();
        
        m_rkf = new GRungeKuttaFehlberg();
        m_stepsize = 10.0;
    }
    
    void GAdams::setStepsize(double stepsize)
    {
        m_stepsize = stepsize;
        m_rkf->setStepsize(m_stepsize);
    }
    
    void GAdams::clearStartup()
    {
        m_startup.clear();
    }
    
    
//    // this method: bad precision
//    // get the adams-bashforth coefficients
//    // http://www.udsspace.uds.edu.gh/bitstream/123456789/294/1/A%20MATRIX%20SYSTEM%20FOR%20COMPUTING%20THE%20COEFFICIENTS%20OF%20THE%20ADAMS%20BASHFORTH-MOULTON%20PREDICTOR-CORRECTOR%20FORMULAE.pdf
//    void GAdams::getBashforthCoeff()
//    {
//        
//        GMatrix A(m_order,m_order);
//        GMatrix L(m_order,1);
//        L(0,0) = 1.0;
//        for( int i = 0 ; i< m_order; i++ )
//        {
//            A(0,i) = 1.0;
//        }
//        int t = 1;
//        for( int i =1; i< m_order; i++)
//        {
//            for( int j = 1; j< m_order; j++)
//            {
//                A(i,j) = pow(j,i*1.0);
//            }
//            
//            t = t*(-1);
//            L(i,0) = t/(i+1.0);
//        }
//        
//        // solve the coefficients
//        GMatrix X = (!A)*L;
//        
//        //cout<< A;
//        //cout<< ~L <<endl;
//        
//        X.getData(&m_bC[0]);
//        
//       // cout<< ~X;
//        
//        //
////        double b = 1440.0;
////        m_bC[0] = 4277.0/b;
////        m_bC[1] = -7923.0/b;
////        m_bC[2] = 9482.0/b;
////        m_bC[3] = -6798.0/b;
////        m_bC[4] = 2627.0/b;
////        m_bC[5] = -475.0/b;
//        
//        int testc = 0;
//        
//    }
//    
//    // this method: bad precision
//    void GAdams::getMoultonCoeff()
//    {
//        
//        GMatrix A(m_order,m_order);
//        GMatrix L(m_order,1);
//        for( int i = 0 ; i< m_order ; i++ )
//        {
//            for( int j = 0 ; j< m_order; j++)
//            {
//                A(i,j) = pow((j-1),i);
//            }
//            L(i,0) = pow(-1,i)/(i+1);
//        }
//        
//        // solve the coefficients
//        GMatrix X = (!A)*L;
//        X.getData(&m_mC[0]);
//        
////        double b = 1440.0;
////        m_mC[0] = 475.0/b;
////        m_mC[1] = 1427.0/b;
////        m_mC[2] = -798.0/b;
////        m_mC[3] = 482.0/b;
////        m_mC[4] = -173.0/b;
////        m_mC[5] = 27.0/b;
//        cout << ~X;
//        int testc = 0;
//        
//    }
//    
    
    void GAdams::getCoef()
    {
        double *c = new double[m_order];
        double *g = new double[m_order];
        c[0] = 1.0;
        for( int n = 1 ; n< m_order; n++ )
        {
            double s = 0.0;
            for( int i =0 ; i<= n-1; i++)
            {
                s += c[i]/(n+1-i);
            }
            c[n] = -s;
        }
        
        for( int n = 0 ; n< m_order ; n++)
        {
            double s = 0.0;
            for(int k = 0 ; k<=n; k++)
            {
                s += c[k];
            }
            g[n] = s;
        }
        
        
        //starting the bashforth coefficients
        for( int n = m_order-1 ;n >= 0 ; n-- )
        {
            int sign = -1;
            for( int i = 0 ; i<=n; i++ ) // m_order
            {
                sign *= -1;
                m_bC[i] += sign*g[n]*GMath::nchoosek(n, i);
            }
        }
        
        //starting the moulton coefficients
        for( int n = m_order-1 ;n >= 0 ; n-- )
        {
            int sign = -1;
            for( int i = 0 ; i<=n; i++ ) // m_order
            {
                sign *= -1;
                m_mC[i] += sign*c[n]*GMath::nchoosek(n, i);
            }
        }
        
        if( c!= NULL) {delete[] c; c = NULL;}
        if( g!= NULL) {delete[] g; g = NULL;}
        
        int testc = 0;
        
    }
    
    
    void GAdams::singleStep(gfc::GMotionDynamic *porb, int n, double x0, double *y0, double h, double *y)
    {
        
        double dt = -h;
        double tstep = dt/m_order; // the stepsize of RFK for start up information
        
        double tt = x0;
        int count = 0;
        double* yp = new double[n];
        memset(yp,0,sizeof(double)*n);
        GMatrix mp(n,1);
        
        if( m_startup.size() == 0 || m_startup.size() < m_order ) // the start up process
        {
            
            double* yerr = new double[n];
            double* yold = new double[n];
            memcpy(yold,y0,sizeof(double)*n);
            
            //reset the step size of RKF
            m_rkf->setStepsize(fabs(tstep));
            
            //get the start-up information
            for( int i = 0 ; i< m_order; i++ )
            {
                
                if( i == 0 )
                {
                    porb->getDerivatives(n, x0, y0, yp);
                    mp.setData(yp, n, 1);
                    m_startup.push_back(mp);
                }
                else
                {
                    
                    m_rkf->IntegrateTo(porb, n, tt, yold, tt+dt, y);
                    //m_rkf->singleStep(porb, n, tt, yold, dt, y, yerr);
                    
                    porb->getDerivatives(n, tt+dt, y, yp);
                    mp.setData(yp, n, 1);
                    m_startup.push_back(mp);
                    
                    memcpy(yold, y, sizeof(double)*n);
                    
                    tt += dt;
                }
                
            }
            
            if(yerr != NULL) {delete[] yerr; yerr = NULL;}
            if(yold != NULL) {delete[] yold; yold = NULL;}
            
        } // end of if
        
        // execuate the multi-step process
        /*
         Bashforth:
         Yn+1 = Yn + h*( a1*Y'n + a2*Y'n-1 + a3*Y'n-2 + ... + ak*Y'n-k+1 )
         */
        std::deque<GMatrix> YP = m_startup;
        
        double* Ynp1 = new double[n];
        memset(Ynp1, 0, sizeof(double)*n);
        
        double* ynp1 = new double[n];
        memset(ynp1, 0, sizeof(double)*n);
        
        
        //GMatrix p;
        for(int k = 0 ; k< m_order ; k++ )
        {
            mp = YP.front();
            YP.pop_front();
            for( int i = 0 ; i< n; i++ )
            {
                Ynp1[i] += m_bC[k]*mp[i];
            }
        }
        
        for( int i = 0 ; i< n; i++ )
        {
            Ynp1[i] = y0[i] + h*Ynp1[i];
        }
        
        memcpy(ynp1, Ynp1, sizeof(double)*n);
        
        /*
         Moulton:
         Yn+1 = Yn + h*( b1*Y'n+1 + b2*Y'n + b3*Y'n-1 + ... + bk*Y'n-k+2 )
         */
        //update the m_startup firstly
        porb->getDerivatives(n, x0 + h , ynp1, yp);
        mp.setData(yp, n, 1);
        m_startup.pop_back();
        m_startup.push_front(mp);
        
        
        YP = m_startup;
        
        memset(Ynp1, 0, sizeof(double)*n);
        //GMatrix p;
        for( int k = 0 ; k< m_order ; k++ )
        {
            mp = YP.front();
            YP.pop_front();
            for( int i = 0 ; i< n; i++ )
            {
                Ynp1[i] += m_mC[k]*mp[i];
            }
        }
        
        for( int i = 0 ; i< n; i++ )
        {
            Ynp1[i] = y0[i] + h*Ynp1[i];
        }
        
        
        //update the m_startup firstly
        porb->getDerivatives(n, x0 + h , Ynp1, yp);
        mp.setData(yp, n, 1);
        m_startup.pop_front();
        m_startup.push_front(mp);
        
        //after that, obtain the final Ynp1
        memcpy(y,Ynp1,sizeof(double)*n);
        
        if(yp!= NULL ) {delete[] yp; yp = NULL;}
        if(ynp1!= NULL ) {delete[] ynp1; ynp1 = NULL;}
        if(Ynp1!= NULL ) {delete[] Ynp1; Ynp1 = NULL;}
        
    }
    
    
    void GAdams::IntegrateTo(GMotionDynamic* porb,  int n , double start, double* ystart, double end,double* yend )
    {
        
        bool rec = false;
        
        int sign = -1; // determine the direction of intergration
        
        if( end >= start )
        {
            sign = 1;
        }
        
        double* yold = new double[n];
        memcpy(yold,ystart,sizeof(double)*n);
        
        double dt = m_stepsize*sign;
        double tt = start;
        
        while( (end - tt)*sign > 0.0 )
        {
            if( ( (tt + m_stepsize*sign) - end )*sign >= 0.0 ) // here must include the equal sign
            {
                dt = end - tt;
                rec = true;
                break;
            }
            
            singleStep( porb,n, tt, yold, dt ,yend);
            memcpy(yold, yend, sizeof(double)*n);
            tt+= m_stepsize*sign;
        }
        
        // how to deal with the last step ???
        // the current method with big error !!!, can't reduce the step size
        
        /*
         Solution:
         transform the history state to Nordsieck vector by a matrix T
         then change the step size at this Nordsieck vector
         at last, transform back to get the new history state with new step size
         */
        
        if( rec == true && fabs(dt) > 1.0E-10)
        {
            singleStep(porb,n, tt, yold,dt, yend);
        }
        
        
        if( yold != NULL )
        {
            delete[] yold;
            yold = NULL;
        }
        
    
        
    }  // end of function intergrateTo
    
    
    
    
    
}  // end of namespace

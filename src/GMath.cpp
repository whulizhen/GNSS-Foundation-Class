
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
//  GMath.cpp
//  GFC
//
//  Created by lizhen on 15/10/26.
//  Copyright © 2015年 lizhen. All rights reserved.
//

#include "GMath.hpp"



namespace gfc
{
    
    
    std::vector< std::vector<double> > GMath::AllanDeviation(std::vector<double>& phase, double tau0 ) throw(GException)
    {
        
        std::vector< std::vector<double> > result;
        
        int N = (int)(phase.size() - 1);
        int numGaps = 0;
        if( N < 1 )
        {
            GException e("Need more than 2 point to compute a meaningful allan variance.");
            GFC_THROW(e);
        }
        
        // Actual Overlapping Allan Deviation Calculation is done here
        // The Overlapping Allan Deviation is calculated as follows
        //  Sigma^2(Tau) = 1 / (2*(N-2*m)*Tau^2) * Sum(X[i+2*m]-2*X[i+m]+X[i], i=1, i=N-2*m)
        //  Where Tau is the averaging time, N is the total number of points, and Tau = m*Tau0
        //  Where Tau0 is the basic measurement interval
        double sum , sigma = 0.0 ;
        for(int m = 1; m <= (N-1)/2; m++)
        {
            double tau = m*tau0;
            sigma = 0;
            
            for(int i = 0; i < (N-2*m); i++)
            {
                sum = 0;
                if((phase[i+2*m]==0 ||  phase[i+m]==0 || phase[i]==0)
                   && i!=0 && i!=(N-2*m-1))
                    numGaps++;
                else
                    sum = phase[i+2*m] - 2*phase[i+m] + phase[i];
                    sigma += sum * sum;
                    }
            
            sigma = sigma / (2.0*((double)N-(double)numGaps-0-2.0*(double)m)*tau*tau);
            sigma = sqrt(sigma);
            
            std::vector<double> tmpV;
            
            tmpV.push_back(tau);
            tmpV.push_back(sigma);
            
            result.push_back(tmpV);
            
            //deviation.push_back(sigma);
            //time.push_back(tau);
        }
        
        return result;
    } // end of function AllanDeviation
    
    
    /* Computes the Gamma function using a simple Lanczos approximation.
     *
     * This implementation typically gives 15 correct decimal places, and
     * it is adapted from free Python code found in:
     *
     * http://en.wikipedia.org/wiki/Lanczos_approximation
     *
     * \warning Be aware that Gamma function is not defined for 0, -1, -2,...
     */
    double GMath::gamma(double val)
    {
        
        double inf( 9.0e+99 );
        
        if( val == 0.0 )
        {
            return inf;
        }
        
        if ( ( val < 0.0 ) &&
            ( floor(val) == val ) )
        {
            return inf;
        }
        
        // Set the number of coefficients being used
        int g(7);
        const double lanczos_coef[] = { 0.99999999999980993,
            676.5203681218851,
            -1259.1392167224028,
            771.32342877765313,
            -176.61502916214059,
            12.507343278686905,
            -0.13857109526572012,
            9.9843695780195716e-6,
            1.5056327351493116e-7 };
        
        if( val < 0.5 )
        {
            return ( GCONST("PI") / (sin( GCONST("PI")*val)*gamma(1.0-val)) );
        }
        else
        {
            val -= 1.0;
            
            double x( lanczos_coef[0] );
            
            for(int i = 1; i<g+2; i++)
            {
                x += lanczos_coef[i]/(val+(double)i);
            }
            
            double t (val + static_cast<double>(g) + 0.5);
            
            return ( 2.5066282746310002 * pow( t, (val+0.5) ) * exp(-t) * x );
            
        }
        
    }  // End of function 'gamma()'
    
    
    
    /* Computes the natural logarithm of Gamma function
     * using the Lanczos approximation.
     *
     * \warning This version does not work for values <= 0.0
     */
    double  GMath::lngamma(double val)
    {
        
        double inf( 9.0e+99 );
        
        if( val <= 0.0 )
        {
            return inf;
        }
        
        // Set the number of coefficients being used
        int g(7);
        const double lanczos_coef[] = { 0.99999999999980993,
            676.5203681218851,
            -1259.1392167224028,
            771.32342877765313,
            -176.61502916214059,
            12.507343278686905,
            -0.13857109526572012,
            9.9843695780195716e-6,
            1.5056327351493116e-7 };
        
        if(val < 0.5)
        {
            return ( 1.1447298858494002 - (log(sin(GCONST("PI")*val)) + lngamma(1.0-val)) );
        }
        else
        {
            val -= 1.0;
            
            double x( lanczos_coef[0] );
            
            for( int i = 1; i<g+2; i++ )
            {
                x += lanczos_coef[i]/(val+(double)i);
            }
            
            double t (val + static_cast<double>(g) + 0.5);
            
            return ( 0.918938533204672741781 + (val+0.5)*log(t) + (-t) + log(x) );
        }
        
        
    }  // End of function 'lngamma()'
    
    
    
    // Auxiliar Kummer function.
   // double kummerFunc(const double a, const double z);
    
    // We compute the lower incomplete gamma function g(a,z) using the
    // formula:
    //
    //        g(a,z) = z**a * exp(-z) * S(a,z) / a
    //
    // where:
    //
    //                  oo
    //                 ___            k
    //                \              z
    //   S(a,z) = 1 +  )     ------------------.
    //                /___   (a+1)(a+2)...(a+k)
    //                k = 1
    //
    // S(a,z) is computed with the Kummer function "kummerFunc()".
    double GMath::kummerFunc(const double aval, const double zval)
    {
        
        double eps(1.0e-15);    // Small threshold controling precision
        
        double z( abs(zval) );    // We only allow positive values of 'z'
        double a( abs(aval) );    // We only allow positive values of 'a'
        
        double den(a);          // Variable to store denominator
        
        double s(1.0);          // Result will be stored here
        
        double coef(1.0);       // Initialize coefficient
        
        while( abs(coef) > eps)
        {
            coef = coef*z;   // Compute numerator
            den += 1.0;
            coef = coef/den;  // Compute coefficient
            s += coef;        // Add new coefficient to result
        }
        
        return s;
        
    }  // End of function 'kummerFunc()'
    
    
    
    // Lower incomplete gamma function.
    double GMath::lower_gamma(const double a, const double z)
    {
        
        double zp( abs(z) );    // We only allow positive values of 'z'
        double ap( abs(a) );    // We only allow positive values of 'a'
        
        double s( kummerFunc(ap, zp) );
        
        return exp(log(zp)*ap) * exp(-zp) * s / ap;
        
    }  // End of function 'lower_gamma()'
    
    
    
    // Upper incomplete gamma function.
    double GMath::upper_gamma(const double a, const double z)
    {
        
        return ( gamma(a) - lower_gamma(a, z) );
        
    }  // End of function 'upper_gamma()'
    
    // Lower incomplete regularized gamma function P(a,z).
    double GMath::gammaP(const double a, const double z)
    {
        return ( lower_gamma(a,z) / gamma(a) );
    }
    
    
    
    // Upper incomplete regularized gamma function Q(a,z).
    double GMath::gammaQ(const double a, const double z)
    {
        return ( 1.0 - gammaP(a,z) );
    }
    
    
    /* Computes factorial of integer number n.
     *
     * This implementation typically gives 15 correct decimal places, and
     * returns the result as double.
     */
    double GMath::factorial(const int n)
    {
        
        // Return 0 if n<0
        if( n < 0 ) return 0.0;
        
        double facttable[] = { 1.0, 1.0, 2.0, 6.0, 24.0, 120.0, 720.0, 5040.0,
            40320.0, 362880.0, 3628800.0, 39916800.0,
            479001600.0, 6227020800.0, 87178291200.0,
            1307674368000.0 };
        
        // Use a look-up table for small values of n
        if( (n >= 0) && (n <= 15) )
        {
            return facttable[n];
        }
        
        // Use gamma function properties for big values of n
        return gamma(n+1.0);
    }
    
    /* Computes factorial of double number d.
     */
    double GMath::factorial(const double d)
    {
        // copy the input
        double n(d);
        
        double returnValue(0.0);
        
        if (n < 0.0)
        {
            return returnValue;
        }
        else if ((n == 0.0)||(n == 1.0))
        {
            returnValue = 1.0;
            
        } else
        {
            returnValue = n * factorial(n - 1.0);
        }
        
        return returnValue;
    }
    
    
    //calculate combination : Cnm = n!/m!/(n-m)!

    double GMath::nchoosek(int n, int m)
    {
        double res = 0.0;
        
        if( m == 0 )
        {
            return 1;
        }
        
        if( m == n )
        {
            return 1;
        }
        
        if( m > n )
        {
            return 0.0;
        }
        
        if ( m > n/2.0 )
        {
            m = n - m ;
        }
        double s1 =0.0, s2 =0.0;
        for(int i = m+ 1; i<=n; i++)
        {
            s1 += log((double)i);
        }
        
        for( int  i = 2; i<= n-m ; i++)
        {
            s2 += log((double)i);
        }
        
        res = exp(s1-s2);
        
        if( res < (numeric_limits<long>::max)() )
        {
           res =  static_cast<long>(res + 0.5);
        }
        
        return res;
    }
    
    
    /* Error function.
     *  Gauss error function
     *  https://en.wikipedia.org/wiki/Error_function
     *  This is a C++ implementation of the free Python code found in:
     *
     *   http://code.activestate.com/recipes/576391/
     *
     * Such code was based in a C code base with OpenBSD license from:
     *
     * ====================================================
     * Copyright (C) 1993 by Sun Microsystems, Inc. All rights reserved.
     *
     * Developed at SunPro, a Sun Microsystems, Inc. business.
     * Permission to use, copy, modify, and distribute this
     * software is freely granted, provided that this notice
     * is preserved.
     * ====================================================
     */
    double GMath::erf(const double x)
    {
        
        /*
         *                            |x
         *                    2       |\
         *     erf(x)  =  ---------   |   exp(-t*t)dt
         *                 sqrt(pi)  \|
         *                            |0
         *
         *     erfc(x) =  1-erf(x)
         *  Note that
         *              erf(-x) = -erf(x)
         *             erfc(-x) = 2 - erfc(x)
         *
         * Method:
         *      1. For |x| in [0, 0.84375]
         *          erf(x)  = x + x*R(x^2)
         *          erfc(x) = 1 - erf(x)           if x in [-.84375,0.25]
         *                  = 0.5 + ((0.5-x)-x*R)  if x in [0.25,0.84375]
         *
         *         where R = P/Q where P is an odd poly of degree 8 and
         *         Q is an odd poly of degree 10.
         *
         *                      | R - (erf(x)-x)/x | <= 2^(-57.90)
         *
         *
         *         Remark. The formula is derived by noting
         *            erf(x) = (2/sqrt(pi))*(x - x^3/3 + x^5/10 - x^7/42 + ....)
         *         and that
         *            2/sqrt(pi) = 1.128379167095512573896158903121545171688
         *         is close to one. The interval is chosen because the fix
         *         point of erf(x) is near 0.6174 (i.e., erf(x)=x when x is
         *         near 0.6174), and by some experiment, 0.84375 is chosen to
         *         guarantee the error is less than one ulp for erf.
         *
         *      2. For |x| in [0.84375,1.25], let s = |x| - 1, and
         *         c = 0.84506291151 rounded to single (24 bits)
         *              erf(x)  = sign(x) * (c  + P1(s)/Q1(s))
         *              erfc(x) = (1-c)  - P1(s)/Q1(s) if x > 0
         *                        1+(c+P1(s)/Q1(s))    if x < 0
         *              |P1/Q1 - (erf(|x|)-c)| <= 2^(-59.06)
         *
         *         Remark: here we use the taylor series expansion at x=1.
         *              erf(1+s) = erf(1) + s*Poly(s)
         *                       = 0.845.. + P1(s)/Q1(s)
         *         That is, we use rational approximation to approximate
         *              erf(1+s) - (c = (single)0.84506291151)
         *         Note that |P1/Q1|< 0.078 for x in [0.84375,1.25]
         *         where
         *              P1(s) = degree 6 poly in s
         *              Q1(s) = degree 6 poly in s
         *
         *      3. For x in [1.25,1/0.35(~2.857143)],
         *              erfc(x) = (1/x)*exp(-x*x-0.5625+R1/S1)
         *              erf(x)  = 1 - erfc(x)
         *         where
         *              R1(z) = degree 7 poly in z, (z=1/x^2)
         *              S1(z) = degree 8 poly in z
         *
         *      4. For x in [1/0.35,28]
         *              erfc(x) = (1/x)*exp(-x*x-0.5625+R2/S2) if x > 0
         *                      = 2.0 - (1/x)*exp(-x*x-0.5625+R2/S2) if -6<x<0
         *                      = 2.0 - tiny      (if x <= -6)
         *              erf(x)  = sign(x)*(1.0 - erfc(x)) if x < 6, else
         *              erf(x)  = sign(x)*(1.0 - tiny)
         *         where
         *              R2(z) = degree 6 poly in z, (z=1/x^2)
         *              S2(z) = degree 7 poly in z
         *
         *         Note1:
         *            To compute exp(-x*x-0.5625+R/S), let s be a single
         *            precision number and s := x; then
         *              -x*x = -s*s + (s-x)*(s+x)
         *              exp(-x*x-0.5626+R/S) =
         *                               exp(-s*s-0.5625)*exp((s-x)*(s+x)+R/S);
         *
         *         Note2:
         *            Here 4 and 5 make use of the asymptotic series
         *                         exp(-x*x)
         *              erfc(x) ~ ---------- * ( 1 + Poly(1/x^2) )
         *                        x*sqrt(pi)
         *            We use rational approximation to approximate
         *              g(s)=f(1/x^2) = log(erfc(x)*x) - x*x + 0.5625
         *            Here is the error bound for R1/S1 and R2/S2
         *              |R1/S1 - f(x)|  < 2**(-62.57)
         *              |R2/S2 - f(x)|  < 2**(-61.52)
         *
         *      5. For inf > x >= 28
         *              erf(x)  = sign(x) *(1 - tiny)  (raise inexact)
         *              erfc(x) = tiny*tiny (raise underflow) if x > 0
         *                      = 2 - tiny if x<0
         *
         *      7. Special case:
         *              erf(0)  = 0, erf(inf)  = 1, erf(-inf) = -1,
         *              erfc(0) = 1, erfc(inf) = 0, erfc(-inf) = 2,
         *              erfc/erf(NaN) is NaN
         */
        
        double inf( 9e+99 );
        
        if (x >= inf)
        {
            return 1.0;
        }
        
        if (x <= -inf)
        {
            return -1.0;
        }
        
        if ( fabs(x) < 0.84375 )
        {
            return erf1(x);
        }
        else if ( ( 0.84375 <= fabs(x) ) &&
                 ( fabs(x) < 1.25 ) )
        {
            return erf2(x);
        }
        else if ( ( 1.25 <= fabs(x) ) &&
                 ( fabs(x) < 2.857142) )
        {
            return erf3(x);
        }
        else if ( ( 2.857142 <= fabs(x) ) &&
                 ( fabs(x) < 6.0 ) )
        {
            return erf4(x);
        }
        
        // If we got here, it means that ( fabs(x) >= 6.0 )
        return erf5(x);
        
    }  // End of function 'erf()'
    
    
    
    // Auxiliar error function #1. erf(x) for x in [0,0.84375].
    double GMath::erf1(const double x)
    {
        
        int i;
        
        // Get the base-2 exponent of input 'x'
        frexp (x , &i);
        
        if ( fabs( static_cast<double>(i) ) > 28.0 )
        {
            if( fabs( static_cast<double>(i) ) > 57.0 )
            {
                return ( 1.128379167095512586316e+0 * x );
            }
            
            return ( x * ( 2.28379167095512586316e-01 ) );
        }
        
        double z( x*x );
        double r( 1.28379167095512558561e-01
                 + z * ( -3.25042107247001499370e-01
                        + z * ( -2.84817495755985104766e-02
                               + z * ( -5.77027029648944159157e-03
                                      + z * -2.37630166566501626084e-05 ) ) ) );
        
        
        double s( 1.0 + z * ( 3.97917223959155352819e-01
                             + z * ( 6.50222499887672944485e-02
                                    + z * ( 5.08130628187576562776e-03
                                           + z * ( 1.32494738004321644526e-04
                                                  + z * -3.96022827877536812320e-06 ) ) ) ) );
        
        return ( x * (1.0 + r/s ) );
        
    }  // End of function 'erf1()'
    
    
    
    // Auxiliar error function #2. erf(x) for x in [0.84375,1.25]
    double GMath::erf2(const double x)
    {
        
        double s( fabs(x) - 1.0 );
        
        double P( -2.36211856075265944077e-03
                 + s * ( 4.14856118683748331666e-01
                        + s * ( -3.72207876035701323847e-01
                               + s * ( 3.18346619901161753674e-01
                                      + s * ( -1.10894694282396677476e-01
                                             + s * ( 3.54783043256182359371e-02
                                                    + s * -2.16637559486879084300e-03 ) ) ) ) ) );
        
        double Q( 1.0 + s * ( 1.06420880400844228286e-01
                             + s * ( 5.40397917702171048937e-01
                                    + s * ( 7.18286544141962662868e-02
                                           + s * ( 1.26171219808761642112e-01
                                                  + s * ( 1.36370839120290507362e-02
                                                         + s * 1.19844998467991074170e-02 ) ) ) ) ) );
        
        if( x >= 0.0 )
        {
            return ( 8.45062911510467529297e-01 + P/Q );
        }
        else
        {
            return ( -8.45062911510467529297e-01 - P/Q );
        }
        
    }  // End of function 'erf2()'
    
    
    
    // Auxiliar error function #3. erf(x) for x in [1.25,2.857142]
    double GMath::erf3(const double xval)
    {
        
        double x0(xval);
        double x( fabs(xval) );
        double s( 1.0/(x*x) );
        
        double R( -9.86494403484714822705e-03
                 + s * ( -6.93858572707181764372e-01
                        + s * ( -1.05586262253232909814e+01
                               + s * ( -6.23753324503260060396e+01
                                      + s * ( -1.62396669462573470355e+02
                                             + s * ( -1.84605092906711035994e+02
                                                    + s * ( -8.12874355063065934246e+01
                                                           + s * -9.81432934416914548592e+00 ) ) ) ) ) ) );
        
        double S( 1.0 + s * ( 1.96512716674392571292e+01
                             + s * ( 1.37657754143519042600e+02
                                    + s * ( 4.34565877475229228821e+02
                                           + s * ( 6.45387271733267880336e+02
                                                  + s * ( 4.29008140027567833386e+02
                                                         + s * ( 1.08635005541779435134e+02
                                                                + s * ( 6.57024977031928170135e+00
                                                                       + s * -6.04244152148580987438e-02 ) ) ) ) ) ) ) );
        
        double r( exp(-x0*x0-0.5625) * exp( (x0-x)*(x0+x)+R/S) );
        
        if( x0 >= 0.0 )
        {
            return ( 1.0 - r/x );
        }
        else
        {
            return ( r/x - 1.0 );
        }
        
    }  // End of function 'erf3()'
    
    
    
    // Auxiliar error function #4. erf(x) for x in [ 2.857142, 6.0 ]
    double GMath::erf4(const double xval)
    {
        
        double x0(xval);
        double x( fabs(xval) );
        double s( 1.0/(x*x) );
        
        double R( -9.86494292470009928597e-03
                 + s * ( -7.99283237680523006574e-01
                        + s * ( -1.77579549177547519889e+01
                               + s * ( -1.60636384855821916062e+02
                                      + s * ( -6.37566443368389627722e+02
                                             + s * ( -1.02509513161107724954e+03
                                                    + s * -4.83519191608651397019e+02 ) ) ) ) ) );
        
        double S( 1.0 + s * ( 3.03380607434824582924e+01
                             + s * ( 3.25792512996573918826e+02
                                    + s * ( 1.53672958608443695994e+03
                                           + s * ( 3.19985821950859553908e+03
                                                  + s * ( 2.55305040643316442583e+03
                                                         + s * ( 4.74528541206955367215e+02
                                                                + s * -2.24409524465858183362e+01 ) ) ) ) ) ) );
        
        double r( exp( -x0 * x0 - 0.5625 ) * exp( (x0-x)*(x0+x) + R/S ) );
        
        if( x0 >= 0.0 )
        {
            return ( 1.0 - r/x );
        }
        else
        {
            return ( r/x - 1.0 );
        }
        
    }  // End of function 'erf4()'
    
    
    
    // Auxiliar error function #5. erf(x) for x in [ 6.0, inf ]
    double GMath::erf5(const double x)
    {
        
        double tiny( 1e-99 );
        
        if ( x > 0.0 )
        {
            return ( 1.0 - tiny );
        }
        else
        {
            return ( tiny - 1.0 );
        }
        
    }  // End of function 'erf5()'
    
    
    
    // Complementary error function.
    double GMath::erfc(const double x)
    {
        
        return ( 1.0 - erf(x) );
        
    }  // End of function 'erfc()'
    
    
    
    /* Inverse of error function.
     *
     * \ warning Value "z" must be in the range (-1, 1)
     */
    double GMath::inverf(const double z)
    {
        
        double inf( 9.0e+99 );
        
        // Check limits
        if( z >= 1.0  ) return inf;
        if( z <= -1.0 ) return -inf;
        
        double z2PI( z*z*GCONST("PI") );
        
        // Compute the approximation found in:
        //    http://en.wikipedia.org/wiki/Error_function
        // The module of this approximation is always under fabs(inverf())
        double x( 0.88622692545275794 * z
                 * ( 1.0 + z2PI * ( 0.083333333333333333
                                   + z2PI * ( 0.014583333333333334
                                             + z2PI * ( 0.0031498015873015874
                                                       + z2PI * ( 0.00075248704805996468
                                                                 + z2PI * ( 0.0001907475361251403
                                                                           + z2PI * ( 1.8780048076923078e-5
                                                                                     + z2PI * ( 1.3623642420578133e-5 ) ) ) ) ) ) ) ) );
        
        double delta(1.0);
        double threshold( 1.0e-10 );
        int iter( 0 );
        
        while ( ( fabs(delta) > threshold ) &&
               ( iter < 100 ) )
        {
            
            // Use the Newton-Raphson method. Denominator is d(erf(x))/dx
            delta = (z - erf(x)) / (1.1283791670955126*std::exp(-x*x));
            
            x = x + delta;
            
            ++iter;
        }
        
        return x;
        
    }  // End of function 'inverf()'
    
    
    
    /* Beta function.
     *
     * \warning This version may not work for values > 130.0
     */
    double GMath::beta(const double x, const double y )
    {
        
        return ( gamma(x) * gamma(y) / gamma( x + y ) );
        
    }  // End of function 'beta()'
    
    
    /* Computes the natural logarithm of Beta function
     * logE(Beta)
     * \warning This version does not work for values <= 0.0
     */
    double GMath::lnbeta(double x, double y)
    {
        
        x = fabs(x);
        y = fabs(y);
        
        return ( lngamma(x) + lngamma(y) - lngamma( x + y ) );
        
    }  // End of function 'lbeta()'
    
    
    
    /* Auxiliar function to compute incomplete beta function.
     * Power series for incomplete beta integral.
     * Use when b*x is small and x not too close to 1.
     */
    double GMath::incompletebetaps(const double x, const double a, const double b)
    {
        
        // Declare working variables
        double epsilon(1.0e-30);
        double big(1.0e99);
        double small(1.0e-99);
        double maxgam(171.0);   // Approximate maximum input for gamma function
        double ai( 1.0/a );
        double u( (1.0-b) * x );
        double v( u / (a+1.0) );
        double t1(v);
        double t(u);
        double n(2.0);
        double s(0.0);
        double z( epsilon*ai );
        
        while( std::fabs(v) > z )
        {
            u = (n-b) * x / n;
            t = t * u;
            v = t/(a+n);
            s = s+v;
            n = n+1.0;
        }
        
        s = s + t1 + ai;
        u = a * std::log(x);
        
        // Check if we can compute gamma values
        if( (a+b < maxgam) &&
           ( std::fabs(u) < std::log(big) ) )
        {
            
            t = gamma(a+b) / ( gamma(a) * gamma(b) );
            s = s * t * pow(x, a);
            
        }
        else
        {
            
            t = lngamma(a+b) - lngamma(a)- lngamma(b) + u + std::log(s);
            
            if( t < std::log(small) )
            {
                s = 0.0;
            }
            else
            {
                s = std::exp(t);
            }
            
        }  // End of block 'if( (a+b < maxgam) && ( std::fabs(u) ...'
        
        return s;
        
    }  // End of function 'incompletebetaps()'
    
    
    
    // Auxiliar function to compute incomplete beta function.
    double GMath::incompletebetafe(const double x, const double a, const double b)
    {
        
        // Declare auxiliar parameters
        double epsilon(1.0e-30);
        double big( 1.0e16 );
        double biginv( 1.0/big );
        
        double k1( a );
        double k2( a+b );
        double k3( a );
        double k4( a+1.0 );
        double k5( 1.0 );
        double k6( b-1.0 );
        double k7( k4 );
        double k8( a+2.0 );
        
        double pkm1( 1.0 );
        double pkm2( 0.0 );
        double qkm1( 1.0 );
        double qkm2( 1.0 );
        
        double ans( 1.0 );
        double r( 1.0 );
        double thresh( 3.0 * epsilon );
        
        int n( 0 );
        
        // Start iteration
        do
        {
            
            double xk( -x * k1 * k2 / (k3*k4) );
            double pk( pkm1 + pkm2 * xk );
            double qk( qkm1 + qkm2 * xk );
            
            pkm2 = pkm1;
            pkm1 = pk;
            qkm2 = qkm1;
            qkm1 = qk;
            
            xk = x * k5 * k6 / (k7*k8);
            pk = pkm1 + pkm2 * xk;
            qk = qkm1 + qkm2 * xk;
            
            pkm2 = pkm1;
            pkm1 = pk;
            qkm2 = qkm1;
            qkm1 = qk;
            
            if( qk != 0.0 )
            {
                r = pk/qk;
            }
            
            double t;
            
            if( r != 0.0 )
            {
                t = std::fabs( (ans-r)/r );
                ans = r;
            }
            else
            {
                t = 1.0;
            }
            
            // Exit if we are under threshold
            if( t < thresh )
            {
                break;
            }
            
            // Recompute auxiliar parameters
            k1 = k1 + 1.0;
            k2 = k2 + 1.0;
            k3 = k3 + 2.0;
            k4 = k4 + 2.0;
            k5 = k5 + 1.0;
            k6 = k6 - 1.0;
            k7 = k7 + 2.0;
            k8 = k8 + 2.0;
            
            // Check if qk, pk are too big
            if( (std::fabs(qk) + std::fabs(pk)) > big )
            {
                pkm2 = pkm2 * biginv;
                pkm1 = pkm1 * biginv;
                qkm2 = qkm2 * biginv;
                qkm1 = qkm1 * biginv;
            }
            
            // Check if qk, pk are too small
            if( (std::fabs(qk) < biginv) ||
               (std::fabs(pk) < biginv ) )
            {
                pkm2 = pkm2 * big;
                pkm1 = pkm1 * big;
                qkm2 = qkm2 * big;
                qkm1 = qkm1 * big;
            }
            
            // Increase n
            n = n+1;
            
        }
        while( n < 300 );
        
        // Return result
        return ans;
        
    }  // End of function 'incompletebetafe()'
    
    
    
    // Auxiliar function to compute incomplete beta function.
    double GMath::incompletebetafe2(const double x, const double a, const double b)
    {
        
        // Declare auxiliar parameters
        double epsilon(1.0e-30);
        double big( 1.0e16 );
        double biginv( 1.0/big );
        
        double k1( a );
        double k2( b-1.0 );
        double k3( a );
        double k4( a+1.0 );
        double k5( 1.0 );
        double k6( a+b );
        double k7( a+1.0 );
        double k8( a+2.0 );
        
        double pkm1( 1.0 );
        double pkm2( 0.0 );
        double qkm1( 1.0 );
        double qkm2( 1.0 );
        
        double z( x / (1.0-x) );
        double ans( 1.0 );
        double r( 1.0 );
        double thresh( 3.0 * epsilon );
        
        int n( 0 );
        
        // Start iteration
        do
        {
            
            double xk( -z * k1 * k2 / (k3*k4) );
            double pk( pkm1 + pkm2 * xk );
            double qk( qkm1 + qkm2 * xk );
            
            pkm2 = pkm1;
            pkm1 = pk;
            qkm2 = qkm1;
            qkm1 = qk;
            
            xk = z * k5 * k6 / (k7*k8);
            pk = pkm1 + pkm2 * xk;
            qk = qkm1 + qkm2 * xk;
            
            pkm2 = pkm1;
            pkm1 = pk;
            qkm2 = qkm1;
            qkm1 = qk;
            
            if( qk != 0.0 )
            {
                r = pk/qk;
            }
            
            double t;
            
            if( r != 0.0 )
            {
                t = std::fabs( (ans-r) / r );
                ans = r;
            }
            else
            {
                t = 1.0;
            }
            
            // Exit if we are under threshold
            if( t < thresh )
            {
                break;
            }
            
            
            // Recompute auxiliar parameters
            k1 = k1 + 1.0;
            k2 = k2 - 1.0;
            k3 = k3 + 2.0;
            k4 = k4 + 2.0;
            k5 = k5 + 1.0;
            k6 = k6 + 1.0;
            k7 = k7 + 2.0;
            k8 = k8 + 2.0;
            
            // Check if qk, pk are too big
            if( (std::fabs(qk) + std::fabs(pk)) > big )
            {
                pkm2 = pkm2 * biginv;
                pkm1 = pkm1 * biginv;
                qkm2 = qkm2 * biginv;
                qkm1 = qkm1 * biginv;
            }
            
            // Check if qk, pk are too small
            if( (std::fabs(qk) < biginv) ||
               (std::fabs(pk) < biginv ) )
            {
                pkm2 = pkm2 * big;
                pkm1 = pkm1 * big;
                qkm2 = qkm2 * big;
                qkm1 = qkm1 * big;
            }
            
            // Increase n
            n = n+1;
            
        }
        while( n < 300 );
        
        // Return result
        return ans;
        
    }  // End of function 'incompletebetafe2()'
    
    
    
    /* Computes the regularized incomplete Beta function Ix(a,b).
     *
     * This code is a C++ implementation and adaptation from code found
     * in Cephes Math Library Release 2.8, copyright by Stephen L. Moshier,
     * released under a BSD license.
     */
    double GMath::regIncompleteBeta(const double x, const double a, const double b)
    throw(InvalidParameter)
    {
        
        // Let's do some basic checks
        if( (a <= 0.0) || (b <= 0.0) )
        {
            InvalidParameter e("Function 'regIncompleteBeta()': 'a' and 'b' must \
                               be greater than zero." );
            GFC_THROW(e);
        }
        
        if( (x < 0.0) || (x > 1.0) )
        {
            InvalidParameter e("Function 'regIncompleteBeta()': 'x' must be \
                               within the interval [0,1]." );
            GFC_THROW(e);
        }
        
        if (x == 0.0) return 0.0;
        if (x == 1.0) return 1.0;
        
        // Declare working parameters
        int flag(0);
        double epsilon(1.0e-30);
        double xx(x);
        double aa(a);
        double bb(b);
        
        
        // Check if bb*xx is small and xx not too close to 1.
        if( (bb*xx <= 1.0) &&
           (xx <= 0.95) )
        {
            return incompletebetaps(xx, aa, bb);
        }
        
        double w( 1.0-xx );
        double t, xc;
        
        if( xx > (aa/(aa+bb)) )
        {
            flag = 1;
            t = aa;
            aa = bb;
            bb = t;
            xc = xx;
            xx = w;
        }
        else
        {
            xc = w;
        }
        
        double result(0.0);
        
        if( (flag==1)      &&
           (bb*xx <= 1.0) &&
           (xx <= 0.95) )
        {
            
            // bb*xx is small and xx not too close to 1.
            t = incompletebetaps(xx, aa, bb);
            
            if( t <= epsilon )
            {
                result = 1.0 - epsilon;
            }
            else
            {
                result = 1.0 - t;
            }
            
            return result;
            
        }
        
        double y( xx * (aa+bb-2.0) - (aa-1.0) );
        
        if( y < 0.0 )
        {
            w = incompletebetafe(xx, aa, bb);
        }
        else
        {
            w = incompletebetafe2(xx, aa, bb) / xc;
        }
        
        t = std::pow(xc, bb);
        t = t * std::pow(xx, aa);
        t = t/aa;
        t = t*w;
        t = t * ( gamma(aa+bb) / ( gamma(aa) * gamma(bb) ) );
        
        if( flag==1 )
        {
            if( t <= epsilon )
            {
                result = 1.0 - epsilon;
            }
            else
            {
                result = 1.0 - t;
            }
        }
        else
        {
            result = t;
        }
        
        return result;
        
    }  // End of function 'incompleteBeta()'
    
    
    // calculate the 7 Helmert transformation parameters
    // transformation from p to q, p-->q
    // q = k*R*p + dx
    void GMath::HelmertParameter(int npoints, GVector*q , GVector* p,double param[])
    {
//        test data
//        with dx = 0.15 dy=-0.41 dz=1.7 m=0.0000245 a = 30deg b= 30deg c = 60deg
//        int npoints = 5;
//        double param[7] = {0.0};
//        GVector p[5],q[5];
//        p[0].x = 242.3244 ;p[0].y = -98721.4320 ;p[0].z = 765392.1258;
//        p[1].x = -66532.4250 ;p[1].y = 4356.2336 ;p[1].z = -3258.2130;
//        p[2].x = 32434.8774 ;p[2].y = -21238.3434 ;p[2].z = 32394.1974;
//        p[3].x = 43435.2324 ;p[3].y = -2345.2356 ;p[3].z = 76634.3230;
//        p[4].x = -3235.9873 ;p[4].y = 43765.4464 ;p[4].z = -7653.6780;
//        
//        q[0].x = 227600.6206 ;q[0].y = -691004.3367 ;q[0].z = 257512.6095;
//        q[1].x = -53328.5451 ;q[1].y = -25016.3194 ;q[1].z = -31408.9545;
//        q[2].x = 44614.9152 ;q[2].y = -18899.3181 ;q[2].z = 14317.8308;
//        q[3].x = 50635.4047 ;q[3].y = -48756.3882 ;q[3].z = 53145.3267;
//        q[4].x = -31438.0846 ;q[4].y = 14771.1936 ;q[4].z = 27894.3347;
        
        GMatrix R(3,3),dR_a(3,3),dR_b(3,3),dR_c(3,3);
        double a = 0.0, b =0.0, c =0.0, k =1.0, dx=0.0,dy=0.0,dz=0.0;
        memset(param, 0, sizeof(double)*7);
        
        GMatrix B(npoints*3,7);  // the design matrix
        GMatrix L(npoints*3,1);  // the observation
        GMatrix X(7,1);        // the solution
        GMatrix V(npoints,1); // the residual
        dx = X[0];
        dy = X[1];
        dz = X[2];
        k =  1.0+X[3];
        a =  X[4];
        b =  X[5];
        c =  X[6];
        
        while(1)
        {
            double sina=sin(a),cosa=cos(a);
            double sinb=sin(b),cosb=cos(b);
            double sinc=sin(c),cosc=cos(c);
            
            // rotaiton matrix R
            R(0,0) = cosa*cosb;
            R(0,1) = -sina*cosc-cosa*sinb*sinc;
            R(0,2) = sina*sinc - cosa*sinb*cosc;
            R(1,0) = sina*cosb;
            R(1,1) = cosa*cosc - sina*sinb*sinc;
            R(1,2) = -cosa*sinc - sina*sinb*cosc;
            R(2,0) = sinb;
            R(2,1) = cosb*sinc;
            R(2,2) = cosb*cosc;
            
            // part a of the dR
            dR_a(0,0) = -sina*cosb;
            dR_a(0,1) = -cosa*cosc+sina*sinb*sinc;
            dR_a(0,2) = cosa*sinc +sina*sinb*cosc;
            dR_a(1,0) = cosb*cosa;
            dR_a(1,1) = -sina*cosc-cosa*sinb*sinc;
            dR_a(1,2) = sina*sinc-cosa*sinb*cosc;
            dR_a(2,0) = 0.0;
            dR_a(2,1) = 0.0;
            dR_a(2,2) = 0.0;
            
            // part b of the dR
            dR_b(0,0) = -cosa*sinb;
            dR_b(0,1) = -cosa*cosb*sinc;
            dR_b(0,2) = -cosa*cosb*cosc;
            dR_b(1,0) = -sina*sinb;
            dR_b(1,1) = -sina*cosb*sinc;
            dR_b(1,2) = -sina*cosb*cosc;
            dR_b(2,0) = cosb;
            dR_b(2,1) = -sinb*sinc;
            dR_b(2,2) = -sinb*cosc;
            
            // part c of the dR
            dR_c(0,0) = 0.0;
            dR_c(0,1) = sina*sinc-cosa*sinb*cosc;
            dR_c(0,2) = sina*cosc+cosa*sinb*sinc;
            dR_c(1,0) = 0.0;
            dR_c(1,1) = -cosa*sinc-sina*sinb*cosc;
            dR_c(1,2) = -cosa*cosc+sina*sinb*sinc;
            dR_c(2,0) = 0.0;
            dR_c(2,1) = cosb*cosc;
            dR_c(2,2) = -sinc*cosb;
            
            for( int i =0; i< npoints; i++ )
            {
                
                double ll[3] = {0.0};
                ll[0] = k*(R(0,0)*p[i].x + R(0,1)*p[i].y + R(0,2)*p[i].z);
                ll[1] = k*(R(1,0)*p[i].x + R(1,1)*p[i].y + R(1,2)*p[i].z);
                ll[2] = k*(R(2,0)*p[i].x + R(2,1)*p[i].y + R(2,2)*p[i].z);
                
                L[i*3+0] = q[i].x - (ll[0] + dx);
                L[i*3+1] = q[i].y - (ll[1] + dy);
                L[i*3+2] = q[i].z - (ll[2] + dz);
                
                //form matrix B
                B(i*3+0,0) = 1.0; B(i*3+1,1) = 1.0; B(i*3+2,2) = 1.0;
                // dk
                B(i*3+0,3) = R(0,0)*p[i].x + R(0,1)*p[i].y + R(0,2)*p[i].z;
                B(i*3+1,3) = R(1,0)*p[i].x + R(1,1)*p[i].y + R(1,2)*p[i].z;
                B(i*3+2,3) = R(2,0)*p[i].x + R(2,1)*p[i].y + R(2,2)*p[i].z;
                
                // da
                B(i*3+0,4) = k*(dR_a(0,0)*p[i].x + dR_a(0,1)*p[i].y + dR_a(0,2)*p[i].z);
                B(i*3+1,4) = k*(dR_a(1,0)*p[i].x + dR_a(1,1)*p[i].y + dR_a(1,2)*p[i].z);
                B(i*3+2,4) = k*(dR_a(2,0)*p[i].x + dR_a(2,1)*p[i].y + dR_a(2,2)*p[i].z);
                
                // db
                B(i*3+0,5) = k*(dR_b(0,0)*p[i].x + dR_b(0,1)*p[i].y + dR_b(0,2)*p[i].z);
                B(i*3+1,5) = k*(dR_b(1,0)*p[i].x + dR_b(1,1)*p[i].y + dR_b(1,2)*p[i].z);
                B(i*3+2,5) = k*(dR_b(2,0)*p[i].x + dR_b(2,1)*p[i].y + dR_b(2,2)*p[i].z);
                
                // dc
                B(i*3+0,6) = k*(dR_c(0,0)*p[i].x + dR_c(0,1)*p[i].y + dR_c(0,2)*p[i].z);
                B(i*3+1,6) = k*(dR_c(1,0)*p[i].x + dR_c(1,1)*p[i].y + dR_c(1,2)*p[i].z);
                B(i*3+2,6) = k*(dR_c(2,0)*p[i].x + dR_c(2,1)*p[i].y + dR_c(2,2)*p[i].z);
                
            }
            
            X = !(~B*B)*(~B)*L;
            
            if(  fabs(X[0])<1.0E-10
               &&fabs(X[1])<1.0E-10
               &&fabs(X[2])<1.0E-10
               &&fabs(X[3])<1.0E-10
               &&fabs(X[4])<1.0E-9/20/3.14
               &&fabs(X[5])<1.0E-9/20/3.14
               &&fabs(X[6])<1.0E-9/20/3.14
               )
            {
                V = L;
                param[0] = dx;
                param[1] = dy;
                param[2] = dz;
                param[3] = k-1.0;
                param[4] = a;
                param[5] = b;
                param[6] = c;
                //printf("%f %f %f %f %f %f %f\n",dx,dy,dz,k,a,b,c);
                // (~V).print();
                break;
            }
            else
            {
                dx+= X[0];
                dy+= X[1];
                dz+= X[2];
                k += X[3];
                a += X[4];
                b += X[5];
                c += X[6];
            }
            
        }
        
    }  // end of function HelmertParameter
    
    // do the Helmert transformation: q = DX + k*R*p, where k = 1+m
    void GMath::HelmertTransform(gfc::GVector &p, double *param, gfc::GVector &q)
    {
        GMatrix R(3,3);
        
        double a = param[4], b= param[5], c = param[6];
        
        double sina=sin(a),cosa=cos(a);
        double sinb=sin(b),cosb=cos(b);
        double sinc=sin(c),cosc=cos(c);
        
        // rotaiton matrix R
        R(0,0) = cosa*cosb;
        R(0,1) = -sina*cosc-cosa*sinb*sinc;
        R(0,2) = sina*sinc - cosa*sinb*cosc;
        R(1,0) = sina*cosb;
        R(1,1) = cosa*cosc - sina*sinb*sinc;
        R(1,2) = -cosa*sinc - sina*sinb*cosc;
        R(2,0) = sinb;
        R(2,1) = cosb*sinc;
        R(2,2) = cosb*cosc;
        
        q.x = param[0] + (1.0+param[3])*( R(0,0)*p.x + R(0,1)*p.y +R(0,2)*p.z );
        q.y = param[1] + (1.0+param[3])*( R(1,0)*p.x + R(1,1)*p.y +R(1,2)*p.z );
        q.z = param[2] + (1.0+param[3])*( R(2,0)*p.x + R(2,1)*p.y +R(2,2)*p.z );
        
        
    }
    
    
    
    //https://github.com/dxli/quarticSolver
    
    // find real roots for quartic equations with real coefficients
    // Author: Dongxu Li
    // Written for the LibreCAD project: http://librecad.org/
    //
    
    // ce is a pointer to an array of equation coefficients
    // root is a pointer to roots to be stored,
    // there's no attempt verify validity of the argument pointer
    //
    unsigned int GMath::quadraticSolver(double * ce,  double * roots)
    //quadratic solver for
    // x^2 + ce[0] x + ce[1] =0
    {
        double discriminant=0.25*ce[0]*ce[0]-ce[1];
        if (discriminant < 0.0 ) return 0;
        roots[0]= -0.5*ce[0] + sqrt(discriminant);
        roots[1]= -ce[0] - roots[0];
        return 2;
    }
    
    
    unsigned int GMath::cubicSolver(double * ce, double *roots)
    //cubic equation solver
    // x^3 + ce[0] x^2 + ce[1] x + ce[2] = 0
    {
        // depressed cubic, Tschirnhaus transformation, x= t - b/(3a)
        // t^3 + p t +q =0
        unsigned int ret=0;
        double shift=(1./3)*ce[0];
        double p=ce[1] -shift*ce[0];
        double q=ce[0]*( (2./27)*ce[0]*ce[0]-(1./3)*ce[1])+ce[2];
        //Cardano's method,
        //	t=u+v
        //	u^3 + v^3 + ( 3 uv + p ) (u+v) + q =0
        //	select 3uv + p =0, then,
        //	u^3 + v^3 = -q
        //	u^3 v^3 = - p^3/27
        //	so, u^3 and v^3 are roots of equation,
        //	z^2 + q z - p^3/27 = 0
        //	and u^3,v^3 are,
        //		-q/2 \pm sqrt(q^2/4 + p^3/27)
        //	discriminant= q^2/4 + p^3/27
        //std::cout<<"cubicSolver:: p="<<p<<"\tq="<<q<<std::endl;
        
        double discriminant= (1./27)*p*p*p+(1./4)*q*q;
        if ( fabs(p)< 1.0e-75)
        {
            ret=1;
            *roots=(q>0)?-pow(q,(1./3)):pow(-q,(1./3));
            *roots -= shift;
            return ret;
        }
        
        //std::cout<<"cubicSolver:: discriminant="<<discriminant<<std::endl;
        
        if( discriminant>0.0 )
        {
            double ce2[2]= {q, -1./27*p*p*p},u3[2];
            ret=quadraticSolver(ce2,u3);
            if (! ret )
            { //should not happen
                std::cerr<<"cubicSolver::Error cubicSolver("<<ce[0]<<' '<<ce[1]<<' '<<ce[2]<<")\n";
            }
            ret=1;
            double u,v;
            u= (q<=0) ? pow(u3[0], 1./3): -pow(-u3[1],1./3);
            //u=(q<=0)?pow(-0.5*q+sqrt(discriminant),1./3):-pow(0.5*q+sqrt(discriminant),1./3);
            v=(-1./3)*p/u;
            
            //std::cout<<"cubicSolver:: u="<<u<<"\tv="<<v<<std::endl;
            //std::cout<<"cubicSolver:: u^3="<<u*u*u<<"\tv^3="<<v*v*v<<std::endl;
            
            *roots=u+v - shift;
            return ret;
        }
        ret=3;
        std::complex<double> u(q,0),rt[3];
        u=pow(-0.5*u-sqrt(0.25*u*u+p*p*p/27),1./3);
        rt[0]=u-p/(3.*u)-shift;
        std::complex<double> w(-0.5,sqrt(3.)/2);
        rt[1]=u*w-p/(3.*u*w)-shift;
        rt[2]=u/w-p*w/(3.*u)-shift;
        //	std::cout<<"Roots:\n";
        //	std::cout<<rt[0]<<std::endl;
        //	std::cout<<rt[1]<<std::endl;
        //	std::cout<<rt[2]<<std::endl;
        
        roots[0]=rt[0].real();
        roots[1]=rt[1].real();
        roots[2]=rt[2].real();
        return ret;
    }
    
    unsigned int GMath::quarticSolver(double * ce, double *roots)
    //quartic solver
    // x^4 + ce[0] x^3 + ce[1] x^2 + ce[2] x + ce[3] = 0
    {
        // x^4 + a x^3 + b x^2 +c x + d = 0
        // depressed quartic, x= t - a/4
        // t^4 + ( b - 3/8 a^2 ) t^2 + (c - a b/2 + a^3/8) t + d - a c /4 + a^2 b/16 - 3 a^4/256 =0
        // t^4 + p t^2 + q t + r =0
        // p= b - (3./8)*a*a;
        // q= c - 0.5*a*b+(1./8)*a*a*a;
        // r= d - 0.25*a*c+(1./16)*a*a*b-(3./256)*a^4
        unsigned int ret=0;
        double shift=0.25*ce[0];
        double shift2=shift*shift;
        double a2=ce[0]*ce[0];
        double p= ce[1] - (3./8)*a2;
        double q= ce[2] + ce[0]*((1./8)*a2 - 0.5*ce[1]);
        double r= ce[3] - shift*ce[2] + (ce[1] - 3.*shift2)*shift2;
        double eps = 1.0E-16;
        //std::cout<<"quarticSolver:: p="<<p<<"\tq="<<q<<"\tr="<<r<<std::endl;
        
        if (fabs(q) <= eps) // q == 0
        {
            // Biquadratic equations
            double discriminant= 0.25*p*p -r;
            
            if(discriminant < 0.)
            {
                return 0;
            }
            
            double t2[2];
            t2[0]=-0.5*p-sqrt(discriminant);  // t2[0] > t2[1]
            t2[1]= -p - t2[0];
            
            if( fabs(discriminant) < eps )
            {
                roots[0]=sqrt(t2[0])-shift;
                roots[1]= -sqrt(t2[0])-shift;
                return 2;
            }
            
            if ( t2[0] >= 0. )
            {// four real roots
                roots[0]=sqrt(t2[0])-shift;
                roots[1]= -sqrt(t2[0])-shift;
                roots[2]=sqrt(t2[1])-shift;
                roots[3]= -sqrt(t2[1])-shift;
                return 4;
            }
            if ( t2[1] >= 0.)
            { // two real roots
                roots[0]=sqrt(t2[1])-shift;
                roots[1]= -roots[0]-shift;
                return 2;
            }
            return 0;
        }
        if ( fabs(r)< eps )
        {
            double cubic[3]= {0.,p,q};
            roots[0]=0.;
            ret=1+cubicSolver(cubic,roots+1);
            for(unsigned int i=0; i<ret; i++) roots[i] -= shift;
            return ret;
        }
        // depressed quartic to two quadratic equations
        // t^4 + p t^2 + q t + r = ( t^2 + u t + v) ( t^2 - u t + w)
        // so,
        // 	p = w+v-u^2
        // 	q = (w-v)*u
        // 	r = wv
        // so when u != 0, here, u can't be 0
        //  (p+u^2)^2 - (q/u)^2 = 4 r
        //  y=u^2,
        //  y^3 + 2 p y^2 + ( p^2 - 4 r) y - q^2 =0
        //
        double cubic[3]= {2.*p,p*p-4.*r,-q*q},croots[3]={0.0};
        ret = cubicSolver(cubic,croots);
        
        //std::cout<<"quarticSolver:: real roots from cubic: "<<ret<<std::endl;
//        for( unsigned int i=0; i<ret; i++ )
//        {
//            std::cout<<"cubic["<<i<<"]="<<cubic[i]<<" x= "<<croots[i]<<std::endl;
//        }
        
        int rootNum = 0; // this is for roots
        
        for( int i = 0 ; i< ret ; i++ )
        {
            int num_root = 0;
            double myroots[4] = {0.0};
            double droots[2] = {0.0};
            double ce2[2];
            
            if( croots[i] < 0.0 ) // only can be positive
            {
                continue;
            }
            
            //double u2 = croots[i];
            double u = sqrt(croots[i]);
            double v = 0.5*(p+croots[i])-0.5*q/u ;
            double w = 0.5*(p+croots[i])+0.5*q/u ;
            ce2[0]=	u;
            ce2[1]= v;
            int ret1=quadraticSolver(ce2,droots);
            for( int j = 0 ; j< ret1 ; j++ )
            {
                if( num_root == 0 )
                {
                    myroots[num_root++] = droots[0];
                }
                else
                {
                    bool exist = false;
                    for(int k = 0 ; k< num_root; k++ )
                    {
                        if( myroots[k] == droots[j])
                        {
                            exist = true;
                            break;
                        }
                    }
                    
                    if(exist == false)
                    {
                        myroots[num_root++] = droots[j];
                    }
                    
                }
                
                
            }
            
            ce2[0]=	-u;
            ce2[1]= w;
            int ret2=quadraticSolver(ce2,droots);
            for( int j = 0 ; j< ret2 ; j++ )
            {
                if( num_root == 0 )
                {
                    myroots[num_root++] = droots[0];
                }
                else
                {
                    bool exist = false;
                    for(int k = 0 ; k< num_root; k++ )
                    {
                        if( myroots[k] == droots[j])
                        {
                            exist = true;
                            break;
                        }
                    }
                    
                    if(exist == false)
                    {
                        myroots[num_root++] = droots[j];
                    }
                    
                }
            }
            
            //comparing myroots with roots
            for( int j = 0 ; j< num_root ; j++ )
            {
                if( rootNum == 0 )
                {
                    roots[rootNum++] = myroots[0];
                }
                else
                {
                    bool exist = false;
                    for(int k = 0 ; k< rootNum; k++ )
                    {
                        if( roots[k] == myroots[j])
                        {
                            exist = true;
                            break;
                        }
                    }
                    
                    if(exist == false)
                    {
                        roots[rootNum++] = myroots[j];
                    }
                }
            }
        }
        
        //do the shift
        for(int i = 0 ; i< rootNum; i++ )
        {
            roots[i] -= shift;
        }
        return rootNum;
        
        
//        
//        if( ret==1 )
//        { //one real root from cubic
//            if (croots[0]< 0.0)
//            {//this should not happen
//                std::cerr<<"Quartic Error:: Found one real root for cubic, but negative\n";
//                return 0;
//            }
//            double sqrtz0=sqrt(croots[0]);
//            double ce2[2];
//            ce2[0]=	-sqrtz0;
//            ce2[1]=0.5*(p+croots[0])+0.5*q/sqrtz0;
//            ret=quadraticSolver(ce2,roots);
//            if (! ret )
//            {
//                ce2[0]=	sqrtz0;
//                ce2[1]=0.5*(p+croots[0])-0.5*q/sqrtz0;
//                ret=quadraticSolver(ce2,roots);
//            }
//            ret=2;
//            for(unsigned int i=0; i<ret; i++) roots[i] -= shift;
//            return ret;
//        }
//        if ( croots[0]> 0. && croots[1] > 0. )
//        {
//            double sqrtz0=sqrt(croots[0]);
//            double ce2[2];
//            ce2[0]=	-sqrtz0;
//            ce2[1]=0.5*(p+croots[0])+0.5*q/sqrtz0;
//            ret=quadraticSolver(ce2,roots);
//            ce2[0]=	sqrtz0;
//            ce2[1]=0.5*(p+croots[0])-0.5*q/sqrtz0;
//            ret=quadraticSolver(ce2,roots+2);
//            // two equal solutions
//            if(ret == 2  && roots[0] == roots[2] )
//            {
//                ret = 1;
//            }
//            
//            for(unsigned int i=0; i<ret; i++) roots[i] -= shift;
//            return ret;
//        }
//        return 0;
        
        
    }
    
    // the code for cubic spline interplation
    
   void GSpline::set_boundary(gfc::GSpline::bd_type left, double left_value,
                              gfc::GSpline::bd_type right, double right_value,
                              bool force_linear_extrapolation)
    {
        assert(m_x.size()==0);          // set_points() must not have happened yet
        m_left=left;
        m_right=right;
        m_left_value=left_value;
        m_right_value=right_value;
        m_force_linear_extrapolation=force_linear_extrapolation;
    }
    
    void GSpline::set_points(const std::vector<double>& x,
                            const std::vector<double>& y, bool cubic_spline)
    {
        assert(x.size()==y.size());
        assert(x.size()>2);
        m_x=x;
        m_y=y;
        int   n=x.size();
        // TODO: maybe sort x and y, rather than returning an error
        for( int i=0; i<n-1; i++ )
        {
            assert(m_x[i]<m_x[i+1]);
        }
        
        if(cubic_spline==true)
        {
            // cubic spline interpolation
            // setting up the matrix and right hand side of the equation system
            // for the parameters b[]
            //band_matrix A(n,1,1);
            GMatrix A(n,n);
            
            //std::vector<double>  rhs(n);
            GMatrix rhs(n,1);
            
            for(int i=1; i<n-1; i++)
            {
                A(i,i-1)=1.0/3.0*(x[i]-x[i-1]);
                A(i,i)=2.0/3.0*(x[i+1]-x[i-1]);
                A(i,i+1)=1.0/3.0*(x[i+1]-x[i]);
                rhs[i]=(y[i+1]-y[i])/(x[i+1]-x[i]) - (y[i]-y[i-1])/(x[i]-x[i-1]);
            }
            // boundary conditions
            if(m_left == gfc::GSpline::second_deriv )
            {
                // 2*b[0] = f''
                A(0,0)=2.0;
                A(0,1)=0.0;
                rhs[0]=m_left_value;
            }
            else if(m_left == gfc::GSpline::first_deriv)
            {
                // c[0] = f', needs to be re-expressed in terms of b:
                // (2b[0]+b[1])(x[1]-x[0]) = 3 ((y[1]-y[0])/(x[1]-x[0]) - f')
                A(0,0)=2.0*(x[1]-x[0]);
                A(0,1)=1.0*(x[1]-x[0]);
                rhs[0]=3.0*((y[1]-y[0])/(x[1]-x[0])-m_left_value);
            }
            else
            {
                assert(false);
            }
            if(m_right == gfc::GSpline::second_deriv)
            {
                // 2*b[n-1] = f''
                A(n-1,n-1)=2.0;
                A(n-1,n-2)=0.0;
                rhs[n-1]=m_right_value;
            }
            else if(m_right == gfc::GSpline::first_deriv)
            {
                // c[n-1] = f', needs to be re-expressed in terms of b:
                // (b[n-2]+2b[n-1])(x[n-1]-x[n-2])
                // = 3 (f' - (y[n-1]-y[n-2])/(x[n-1]-x[n-2]))
                A(n-1,n-1)=2.0*(x[n-1]-x[n-2]);
                A(n-1,n-2)=1.0*(x[n-1]-x[n-2]);
                rhs[n-1]=3.0*(m_right_value-(y[n-1]-y[n-2])/(x[n-1]-x[n-2]));
            }
            else
            {
                assert(false);
            }
            
            // solve the equation system to obtain the parameters b[]
            //m_b=A.lu_solve(rhs);
            //cout<<A;
            m_b.resize(n,0.0);
            ((!A)*rhs).getData( &m_b[0] );
            
            // calculate parameters a[] and c[] based on b[]
            m_a.resize(n);
            m_c.resize(n);
            for(int i=0; i<n-1; i++ )
            {
                m_a[i]=1.0/3.0*(m_b[i+1]-m_b[i])/(x[i+1]-x[i]);
                m_c[i]=(y[i+1]-y[i])/(x[i+1]-x[i])
                - 1.0/3.0*(2.0*m_b[i]+m_b[i+1])*(x[i+1]-x[i]);
            }
        }
        else
        {
            // linear interpolation
            m_a.resize(n);
            m_b.resize(n);
            m_c.resize(n);
            for( int i=0; i<n-1; i++ )
            {
                m_a[i]=0.0;
                m_b[i]=0.0;
                m_c[i]=(m_y[i+1]-m_y[i])/(m_x[i+1]-m_x[i]);
            }
        }
        
        // for left extrapolation coefficients
        m_b0 = (m_force_linear_extrapolation==false) ? m_b[0] : 0.0;
        m_c0 = m_c[0];
        
        // for the right extrapolation coefficients
        // f_{n-1}(x) = b*(x-x_{n-1})^2 + c*(x-x_{n-1}) + y_{n-1}
        double h=x[n-1]-x[n-2];
        // m_b[n-1] is determined by the boundary condition
        m_a[n-1]=0.0;
        m_c[n-1]=3.0*m_a[n-2]*h*h+2.0*m_b[n-2]*h+m_c[n-2];   // = f'_{n-2}(x_{n-1})
        if(m_force_linear_extrapolation==true)
            m_b[n-1]=0.0;
    }
    
    double GSpline::operator() (double x) const
    {
        size_t n=m_x.size();
        // find the closest point m_x[idx] < x, idx=0 even if x<m_x[0]
        std::vector<double>::const_iterator it;
        it=std::lower_bound(m_x.begin(),m_x.end(),x);
        int idx=std::max( int(it-m_x.begin())-1, 0);
        
        double h=x-m_x[idx];
        double interpol;
        if(x<m_x[0])
        {
            // extrapolation to the left
            interpol=(m_b0*h + m_c0)*h + m_y[0];
        }
        else if(x>m_x[n-1])
        {
            // extrapolation to the right
            interpol=(m_b[n-1]*h + m_c[n-1])*h + m_y[n-1];
        }
        else
        {
            // interpolation
            interpol=((m_a[idx]*h + m_b[idx])*h + m_c[idx])*h + m_y[idx];
        }
        
        return interpol;
    }
    
    double GSpline::deriv(int order, double x) const
    {
        assert(order>0);
        
        size_t n=m_x.size();
        // find the closest point m_x[idx] < x, idx=0 even if x<m_x[0]
        std::vector<double>::const_iterator it;
        it=std::lower_bound(m_x.begin(),m_x.end(),x);
        int idx=std::max( int(it-m_x.begin())-1, 0);
        
        double h=x-m_x[idx];
        double interpol;
        if(x<m_x[0])
        {
            // extrapolation to the left
            switch(order)
            {
                case 1:
                    interpol=2.0*m_b0*h + m_c0;
                    break;
                case 2:
                    interpol=2.0*m_b0*h;
                    break;
                default:
                    interpol=0.0;
                    break;
            }
        }
        else if(x>m_x[n-1])
        {
            // extrapolation to the right
            switch(order)
            {
                case 1:
                    interpol=2.0*m_b[n-1]*h + m_c[n-1];
                    break;
                case 2:
                    interpol=2.0*m_b[n-1];
                    break;
                default:
                    interpol=0.0;
                    break;
            }
        }
        else
        {
            // interpolation
            switch(order)
            {
                case 1:
                    interpol=(3.0*m_a[idx]*h + 2.0*m_b[idx])*h + m_c[idx];
                    break;
                case 2:
                    interpol=6.0*m_a[idx]*h + 2.0*m_b[idx];
                    break;
                case 3:
                    interpol=6.0*m_a[idx];
                    break;
                default:
                    interpol=0.0;
                    break;
            }
        }
        
        return interpol;
    }
    
    
    
    
    
} // end of namespace gfc




//
//  GMath.hpp
//  GFC
//
//  Created by lizhen on 15/10/25.
//  Copyright © 2015年 lizhen. All rights reserved.
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


#ifndef GFC_GMath_hpp
#define GFC_GMath_hpp

#include <stdio.h>
#include <assert.h>

#include "GException.h"
#include "GFCCONST.h"
#include "GVector.hpp"
#include "GMatrix.h"

#include <cmath>
using namespace std;


//@fix MSVC doesnt like std::sqrt or std::abs, and disabling
//extensions allows abs(double) to be used instead of fabs()
#ifdef _MSC_VER
#undef _MSC_EXTENSIONS
#endif
#include <cmath>
#ifdef _MSC_VER
#define _MSC_EXTENSIONS
#endif


namespace gfc
{
    
    
    //Only provide the static methods and this class cannot be inheried
    // 静态函数的性能需要考虑，考虑将此类改为名字空间
    class GMath
    {
        
    public:
        
    template < class T >
    static  T SQRT( T x, T y )
    {
       return std::sqrt(x,y);
    }
        
    template < class T >
    static  T ABS( T x, T y )
    {
        return std::abs(x,y);
    }
        
    template < class T >
    static  T MAX( T x, T y )
    {
        if(x > y)
        {
            return x;
        }
        else
        {
            return y;
        }
        
        //return std::max(x,y);
    }
      
   template <class T>
   static  T MAX( T *items, size_t NumItems )
    {
        T MaxFound;
        size_t i;
        for( i=0; i<NumItems; i++ )
        {
            if( i == 0 ) MaxFound = items[i];
            else if( MaxFound < items[i] ) MaxFound = items[i];
        }
            
        return MaxFound;
    }
        
    template < class T >
    static  T MIN( T x, T y )
    {
        return std::min(x,y);
    }
        
    template <class T>
    static  T MIN( T *items, size_t NumItems )
    {
        T MinFound;
        size_t i;
        for( i=0; i<NumItems; i++ )
        {
            if( i == 0 ) MinFound = items[i];
            else if( MinFound > items[i] ) MinFound = items[i];
        }
        return MinFound;
    }
    
    template <class T>
    static T ROUND( T x )
    {
       return double(std::floor(x + 0.5) );
    }
        
     
        /// This is a straightforward version of Lagrange Interpolation.
        /// Y must have size at least as large as X, and X.size() must be >= 2;
        /// x should lie within the range of X.
        template <class T>
        static  T SimpleLagrangeInterpolation(const std::vector<T>& X, const std::vector<T>& Y,
                                      const T x) throw(GException)
        {
            if(Y.size() < X.size())
            {
                GFC_THROW(GException("Input vectors must be of same size"));
            }
            size_t i,j;
            T Yx(0);
            for(i=0; i<X.size(); i++)
            {
                if(x==X[i]) return Y[i];
                
                T Li(1);
                for(j=0; j<X.size(); j++)
                {
                    if(i!=j)
                    {
                        Li = Li*((x-X[j])/(X[i]-X[j]));
                    }
                }
                
                Yx += Li*Y[i];
            }
            return Yx;
        }  // end T LagrangeInterpolation(const vector, const vector, const T)
        
       
        
        /// Lagrange interpolation on data (X[i],Y[i]), i=0,N-1 to compute Y(x).
        /// Also return an estimate of the estimation error in 'err'.
        /// This routine assumes that N=X.size() is even and that x is centered on the
        /// interval, that is X[N/2-1] <= x <= X[N/2].
        /// NB This routine will work for N as small as 4, however tests with satellite
        /// ephemerides have shown that N=4 yields m-level errors, N=6 cm-level,
        /// N=8 ~0.1mm level and N=10 ~numerical noise errors; best to use N>=8.
        template <class T>
  static T LagrangeInterpolation(const std::vector<T>& X, const std::vector<T>& Y,
                                const T& x, T& err) throw(GException)
        {
            if( Y.size() < X.size() || X.size() < 4 )
            {
                GFC_THROW( GException("Input vectors must be of same length, at least 4"));
            }
            
            std::size_t i,j,k;
            T y,del;
            std::vector<T> D,Q;
            
            err = T(0);
            k = X.size()/2;
            if(x == X[k]) return Y[k];
            if(x == X[k-1]) return Y[k-1];
            
            if( std::abs( x-X[k-1] ) < std::abs(x-X[k]) ) k=k-1;
            
            for(i=0; i<X.size(); i++ )
            {
                Q.push_back(Y[i]);
                D.push_back(Y[i]);
            }
            y = Y[k--];
            for(j=1; j<X.size(); j++) {
                for(i=0; i<X.size()-j; i++) {
                    del = (Q[i+1]-D[i])/(X[i]-X[i+j]);
                    D[i] = (X[i+j]-x)*del;
                    Q[i] = (X[i]-x)*del;
                }
                err = (2*(k+1) < X.size()-j ? Q[k+1] : D[k--]);    // NOT 2*k
                y += err;
            }
            return y;
        }  // end T LagrangeInterpolation(vector, vector, const T, T&)
        
        // The following is a
        // Straightforward implementation of Lagrange polynomial and its derivative
        // { all sums are over index=0,N-1; Xi is short for X[i]; Lp is dL/dx;
        //   y(x) is the function being approximated. }
        // y(x) = SUM[Li(x)*Yi]
        // Li(x) = PROD(j!=i)[x-Xj] / PROD(j!=i)[Xi-Xj]
        // dy(x)/dx = SUM[Lpi(x)*Yi]
        // Lpi(x) = SUM(k!=i){PROD(j!=i,j!=k)[x-Xj]} / PROD(j!=i)[Xi-Xj]
        // Define Pi = PROD(j!=i)[x-Xj], Di = PROD(j!=i)[Xi-Xj],
        // Qij = PROD(k!=i,k!=j)[x-Xk] and Si = SUM(j!=i)Qij.
        // then Li(x) = Pi/Di, and Lpi(x) = Si/Di.
        // Qij is symmetric, there are only N(N+1)/2 - N of them, so store them
        // in a vector of length N(N+1)/2, where Qij==Q[i+j*(j+1)/2] (ignore i=j).
        
        /// Perform Lagrange interpolation on the data (X[i],Y[i]), i=1,N (N=X.size()),
        /// returning the value of Y(x) and dY(x)/dX.
        /// Assumes that x is between X[k-1] and X[k], where k=N/2 and N > 2;
        /// Warning: for use with the precise (SP3) ephemeris only when velocity is not
        /// available; estimates of velocity, and especially clock drift, not as accurate.
        template <class T>
       static void LagrangeInterpolation(const std::vector<T>& X, const std::vector<T>& Y,
                                   const T& x, T& y, T& dydx) throw(GException)
        {
            if(Y.size() < X.size() || X.size() < 4)
            {
                GFC_THROW( GException("Input vectors must be of same length, at least 4"));
            }
            
            std::size_t i,j,k,N=X.size(),M;
            M = (N*(N+1))/2;
            std::vector<T> P(N,T(1)),Q(M,T(1)),D(N,T(1));
            for(i=0; i<N; i++) {
                for(j=0; j<N; j++) {
                    if(i != j) {
                        P[i] *= x-X[j];
                        D[i] *= X[i]-X[j];
                        if(i < j) {
                            for(k=0; k<N; k++) {
                                if(k == i || k == j) continue;
                                Q[i+(j*(j+1))/2] *= (x-X[k]);
                            }
                        }
                    }
                }
            }
            y = dydx = T(0);
            for(i=0; i<N; i++) {
                y += Y[i]*(P[i]/D[i]);
                T S(0);
                for(k=0; k<N; k++) if(i != k)
                {
                    if(k<i) S += Q[k+(i*(i+1))/2]/D[i];
                    else    S += Q[i+(k*(k+1))/2]/D[i];
                }
                dydx += Y[i]*S;
            }
        }  // end void LagrangeInterpolation(vector, vector, const T, T&, T&)
        
        /// Returns the second derivative of Lagrange interpolation.
        template <class T>
        static T LagrangeInterpolating2ndDerivative(const std::vector<T>& pos,
                                                    const std::vector<T>& val,
                                                    const T desiredPos)
        {
            int degree(pos.size());
            int i,j,m,n;
            
            // First, compute interpolation factors
            typedef std::vector< T > vectorType;
            std::vector< vectorType > delta(degree, vectorType(degree, 0.0));
            
            for(i=0; i < degree; ++i) {
                for(j=0; j < degree; ++j) {
                    if(j != i) {
                        delta[i][j] = ((desiredPos - pos[j])/(pos[i] - pos[j]));
                    }
                }
            }
            
            double retVal(0.0);
            for(i=0; i < degree; ++i) {
                double sum(0.0);
                
                for(m=0; m < degree; ++m) {
                    if(m != i) {
                        double weight1(1.0/(pos[i]-pos[m]));
                        double sum2(0.0);
                        
                        for(j=0; j < degree; ++j) {
                            if((j != i) && (j != m)) {
                                double weight2(1.0/(pos[i]-pos[j]));
                                for(n=0; n < degree; ++n) {
                                    if((n != j) && (n != m) && (n != i)) {
                                        weight2 *= delta[i][n];
                                    }
                                }
                                sum2 += weight2;
                            }
                        }
                        sum += sum2*weight1;
                    }
                }
                retVal += val[i] * sum;
            }
            
            return retVal;
            
        }  // End of 'lagrangeInterpolating2ndDerivative()'
        
        
        
        
        /** Computes the Gamma function using a simple Lanczos approximation.
         *
         * This implementation typically gives 15 correct decimal places, and
         * it is adapted from free Python code found in:
         *
         * http://en.wikipedia.org/wiki/Lanczos_approximation
         *
         * \warning Be aware that Gamma function is not defined for 0, -1, -2,...
         */
      static  double gamma(const double val);
        
        
        /** Computes the natural logarithm of Gamma function
         *  using the Lanczos approximation.
         *
         * \warning This version does not work for values <= 0.0
         */
      static  double lngamma(double val);
        
        
        /// Lower incomplete gamma function.
      static double lower_gamma(const double a, const double z);
        
        
        /// Upper incomplete gamma function.
       static double upper_gamma(const double a, const double z);
        
        
        /// Lower incomplete regularized gamma function P(a,z).
       static double gammaP(const double a, const double z);
        
        
        /// Upper incomplete regularized gamma function Q(a,z).
       static double gammaQ(const double a, const double z);
       
       static double kummerFunc(const double aval, const double zval);
        
        /** Computes factorial of integer number n.
         *
         * This implementation typically gives 15 correct decimal places, and
         * returns the result as double.
         */
       static double factorial(const int n);
        
        /** Computes factorial of double number n.
         *  d < 360 has been tested
         */
       static double factorial(const double d);
        
        
        /** Error function.
         *
         * This is a C++ implementation of the free Python code found in:
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
       static double erf(const double x);
        
        
        /// Complementary error function.
       static double erfc(const double x);
        
        
        /** Inverse of error function.
         *
         * \ warning Value "z" must be in the range (-1, 1)
         */
       static double inverf(const double z);
        
        
        /** Beta function.
         *
         * \warning This version may not work for values > 130.0
         */
       static double beta(const double x, const double y);
        
        
        /** Computes the natural logarithm of Beta function
         *
         * \warning This version does not work for values <= 0.0
         */
       static double lnbeta(double x, double y);
        
        
        /** Computes the regularized incomplete Beta function Ix(a,b).
         *
         * This code is a C++ implementation and adaptation from code found
         * in Cephes Math Library Release 2.8, copyright by Stephen L. Moshier,
         * released under a BSD license.
         */
       static double regIncompleteBeta(const double x, const double a, const double b)
        throw(InvalidParameter);
        
    static std::vector< std::vector<double> > AllanDeviation(std::vector<double>& phase, double tau0 ) throw(GException);
 
        // Auxiliar error function #1. erf(x) for x in [ 0, 0.84375 ]
    static double erf1(const double x);
        // Auxiliar error function #2. erf(x) for x in [ 0.84375, 1.25 ]
        static double erf2(const double x);
        
        
        // Auxiliar error function #3. erf(x) for x in [ 1.25, 2.857142 ]
        static double erf3(const double x);
        
        
        // Auxiliar error function #4. erf(x) for x in [ 2.857142, 6.0 ]
        static double erf4(const double x);
        
        
        // Auxiliar error function #5. erf(x) for x in [ 6.0, inf ]
        static double erf5(const double x);
        
        // calculate the helmert 7 parameters
        /*
         
         */
        static void HelmertParameter(int npoints, GVector*q , GVector* p,double param[]);
        
        static void HelmertTransform(GVector& p, double param[7], GVector& q);
        
        
        /* quadratic, cubic and quartic equation solution */
        
        static unsigned int solveP3(double *x,double a,double b,double c);
        static unsigned int solve_quartic(double a, double b, double c, double d,double x[4] );
        
        static unsigned int quadraticSolver(double * ce,  double * roots);
        static unsigned int cubicSolver(double * ce, double *roots);
        static unsigned int quarticSolver(double * ce, double *roots);
        
        static unsigned int mycubicSolver(double* ce, double* rt);
        
        static double nchoosek(int n, int m);
        
    private:
        GMath();
        ~GMath();
        
       
       static double incompletebetaps(const double x, const double a, const double b);
       
       static double incompletebetafe(const double x, const double a, const double b);
        
       static double incompletebetafe2(const double x, const double a, const double b);
        
        
    };
    
    
    // cubic spline class
    class GSpline
    {
        
    public:
        enum bd_type
        {
            first_deriv =  1,
            second_deriv = 2
        };
        
        // set default boundary condition to be zero curvature at both ends
        GSpline(): m_left(second_deriv), m_right(second_deriv),
        m_left_value(0.0), m_right_value(0.0),
        m_force_linear_extrapolation(false)
        {
            
        }
        
        // optional, but if called it has to come be before set_points()
        void set_boundary(bd_type left, double left_value,
                          bd_type right, double right_value,
                          bool force_linear_extrapolation=false);
        void set_points(const std::vector<double>& x,
                        const std::vector<double>& y, bool cubic_spline=true);
        double operator() (double x) const;
        double deriv(int order, double x) const;
        
        
    private:
        std::vector<double> m_x,m_y;            // x,y coordinates of points
        // interpolation parameters
        // f(x) = a*(x-x_i)^3 + b*(x-x_i)^2 + c*(x-x_i) + y_i
        std::vector<double> m_a,m_b,m_c;        // spline coefficients
        
        double  m_b0, m_c0;                     // for left extrapol
        
        bd_type m_left, m_right;
        double  m_left_value, m_right_value;
        bool    m_force_linear_extrapolation;
        
        
    };
    
    
    
    
}  // end of namespace gfc





#endif /* GMath_hpp */

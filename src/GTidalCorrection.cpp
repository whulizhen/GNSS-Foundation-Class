//
//  GEarthTide.cpp
//  GFC
//
//  Created by lizhen on 22/07/2016.
//  Copyright © 2016 lizhen. All rights reserved.
//

#include "GTidalCorrection.hpp"

namespace gfc
{
    // For dC21 and dS21
    // The coefficients we choose are in-phase(ip) amplitudes and out-of-phase amplitudes of the
    // corrections for frequency dependence, and multipliers of the Delaunay variables
    // Refers to Table 6.5a in IERS2010
    // the multipliers for 5 fundamental arguments
    const double GTidalCorrection::Argu_C21[48][13]=
    {
    //   IP,      OP,      l,    l', F,  D, OMEG,τ,     s,  h,   p,   N′,   ps
        {-0.1,    0,       2,    0,  2,  0,  2 , 1 ,	-3,	0	, 2	, 0,	0   },
        {-0.1,    0,       0,    0,  2,  2,  2 , 1 ,	-3,	2	, 0	, 0,	0   },
        {-0.1,    0,       1,    0,  2,  0,  1 , 1 ,	-2,	0	, 1	,-1,	0   },
        {-0.7,    0.1,     1,    0,  2,  0,  2 , 1 ,	-2,	0	, 1	, 0,	0   },
        {-0.1,    0,      -1,    0,  2,  2,  2 , 1 ,	-2,	2	,-1	, 0,	0   },
        {-1.3,    0.1,     0,    0,  2,  0,  1 , 1 ,	-1,	0	, 0	,-1,	0   },
        {-6.8,    0.6,     0,    0,  2,  0,  2 , 1 ,	-1,	0	, 0	, 0,	0   },
        {0.1,    0,       0,    0,  0,  2,   0 , 1 ,	-1,	2	, 0	, 0,	0   },
        {0.1,    0,       1,    0,  2, -2,   2 , 1 ,	0	,-2	, 1	, 0,	0   },
        {0.1,    0,      -1,    0,  2,  0,   1 , 1 ,	0	, 0	,-1	,-1,	0   },
        {0.4,    0,      -1,    0,  2,  0,   2 , 1 ,	0	, 0	,-1	, 0,	0   },
        {1.3,   -0.1,     1,    0,  0,  0,   0 , 1 ,	0	, 0	, 1	, 0,	0   },
        {0.3,    0,       1,    0,  0,  0,   1 , 1 ,	0	, 0	, 1	, 1,	0   },
        {0.3,    0,      -1,    0,  0,  2,   0 , 1 ,	0	, 2	,-1	, 0,	0   },
        {0.1,    0,      -1,    0,  0,  2,   1 , 1 ,	0	, 2	,-1	, 1,	0   },
        {-1.9,    0.1,     0,    1,  2, -2,  2 , 1 ,	1	,-3	, 0	, 0,	1   },
        {0.5,    0,       0,    0,  2, -2,   1 , 1 ,	1	,-2	, 0	,-1,	0   },
        {-43.4,    2.9,    0,    0,  2, -2,  2 , 1 ,	1	,-2	, 0	, 0,	0   },
        {0.6,    0,       0,   -1,  2, -2,   2 , 1 ,	1	,-1	, 0	, 0,	-1  },
        {1.6,   -0.1,     0,    1,  0,  0,   0 , 1 ,	1	,-1	, 0	, 0,	1   },
        {0.1,    0,      -2,    0,  2,  0,   1 , 1 ,	1	, 0	,-2	,-1,	0   },
        {0.1,    0,       0,    0,  0,  0,  -2 , 1 ,	1	, 0	, 0	,-2,	0   },
        {-8.8,    0.5,     0,    0,  0,  0, -1 , 1 ,	1	, 0	, 0	,-1,	0   },
        {470.9, -30.2,    0,    0,  0,  0,   0 , 1 ,	1	, 0	, 0	, 0,	0   },
        {68.1,  -4.6,     0,    0,  0,  0,   1 , 1 ,	1	, 0	, 0	, 1,	0   },
        {-1.6,    0.1,     0,    0,  0,  0,  2 , 1 ,	1	, 0	, 0	, 2,	0   },
        {0.1,    0,      -1,    0,  0,  1,   0 , 1 ,	1	, 1	,-1	, 0,	0   },
        {-0.1,    0,       0,   -1,  0,  0, -1 , 1 ,	1	, 1	, 0	,-1,	-1  },
        {-20.6,  -0.3,     0,   -1,  0,  0,  0 , 1 ,	1	, 1	, 0	, 0,	-1  },
        {0.3,    0,       0,    1, -2,  2,  -2 , 1 ,	1	, 1	, 0	, 0,	1   },
        {-0.3,    0,       0,   -1,  0,  0,  1 , 1 ,	1	, 1	, 0	, 1,	-1  },
        {-0.2,    0,      -2,    0,  0,  2,  0 , 1 ,	1	, 2	,-2	, 0,	0   },
        {-0.1,    0,      -2,    0,  0,  2,  1 , 1 ,	1	, 2	,-2	, 1,	0   },
        {-5.0,    0.3,     0,    0, -2,  2, -2 , 1 ,	1	, 2	, 0	, 0,	0   },
        {0.2,    0,       0,    0, -2,  2,  -1 , 1 ,	1	, 2	, 0	, 1,	0   },
        {-0.2,    0,       0,   -1, -2,  2, -2 , 1 ,	1	, 3	, 0	, 0,	-1  },
        {-0.5,    0,       1,    0,  0, -2,  0 , 1 ,	2	,-2	, 1	, 0,	0   },
        {-0.1,    0,       1,    0,  0, -2,  1 , 1 ,	2	,-2	, 1	, 1,	0   },
        {0.1,    0,      -1,    0,  0,  0,  -1 , 1 ,	2	, 0	,-1	,-1,	0   },
        {-2.1,    0.1,    -1,    0,  0,  0,  0 , 1 ,	2	, 0	,-1	, 0,	0   },
        {-0.4,    0,      -1,    0,  0,  0,  1 , 1 ,	2	, 0	,-1	, 1,	0   },
        {-0.2,    0,       0,    0,  0, -2,  0 , 1 ,	3	,-2	, 0	, 0,	0   },
        {-0.1,    0,      -2,    0,  0,  0,  0 , 1 ,	3	, 0	,-2	, 0,	0   },
        {-0.6,    0,       0,    0, -2,  0, -2 , 1 ,	3	, 0	, 0	, 0,	0   },
        {-0.4,    0,       0,    0, -2,  0, -1 , 1 ,	3	, 0	, 0	, 1,	0   },
        {-0.1,    0,       0,    0, -2,  0,  0 , 1 ,	3	, 0	, 0	, 2,	0   },
        {-0.1,    0,      -1,    0, -2,  0, -2 , 1 ,	4	, 0	,-1	, 0,	0   },
        {-0.1,    0,      -1,    0, -2,  0, -1 , 1 ,	4	, 0	,-1	, 1,	0   }
        
    };
    
    // For dC22 and dS22
    // Refer to Table 6.5c in IERS2010
    // (0.30102 . i 0.00130).
    // the multipliers for 5 fundamental arguments
    const double GTidalCorrection::Argu_C22[2][13] =
    {
    //    amp  dkfR    l  l' F  D OMEG t   s  h  p  N' ps
        {-0.3, 0.0,    1, 0, 2, 0, 2,  2, -1, 0, 1, 0, 0 },
        {-1.2, 0.0,    0, 0, 2, 0, 2,  2,  0, 0, 0, 0, 0 }
    };
    
    // For dC20
    // Refer to Table 6.5b in IERS2010
    // The nominal value k20 for the zonal tides is taken as 0.30190
    // used with formula 6.8a
    // the multipliers for 5 fundamental arguments
    const double GTidalCorrection::Argu_C20[21][13]=
    {
       //IP     OP,   l,  l', F,  D,  OMEG, t, s, h, p, N', ps
        {16.6, -6.7,  0,  0,  0,  0,  1 , 0,	0 ,	0	, 0 ,	1	, 0 },
        {-0.1,   0.1,  0,  0,  0,  0,  2, 0,	0 ,	0	, 0 ,	2	, 0 },
        {-1.2,   0.8,  0, -1,  0,  0,  0, 0,	0 ,	1	, 0 ,	0	,-1 },
        {-5.5,   4.3,  0,  0, -2,  2, -2, 0,	0 ,	2	, 0 ,	0	, 0 },
        {0.1,  -0.1,  0,  0, -2,  2, -1 , 0,	0 ,	2	, 0 ,	1	, 0 },
        {-0.3,   0.2,  0, -1, -2,  2, -2, 0,	0 ,	3	, 0 ,	0	,-1 },
        {-0.3,   0.7,  1,  0,  0, -2,  0, 0,	1 ,-2   , 1 ,	0	, 0 },
        {0.1,  -0.2, -1,  0,  0,  0, -1 , 0,	1 ,	0	,-1 ,  -1	, 0 },
        {-1.2,   3.7, -1,  0,  0,  0,  0, 0,	1 ,	0	,-1 ,	0	, 0 },
        {0.1,  -0.2, -1,  0,  0,  0,  1 , 0,	1 ,	0	,-1 ,	1	, 0 },
        {0.1,  -0.2,  1,  0, -2,  0, -2 , 0,	1 ,	0	, 1 ,	0	, 0 },
        {0,     0.6,  0,  0,  0, -2,  0 , 0,	2 ,-2	, 0 ,	0	, 0 },
        {0,     0.3, -2,  0,  0,  0,  0 , 0,	2 ,	0	,-2 ,	0	, 0 },
        {0.6,   6.3,  0,  0, -2,  0, -2 , 0,	2 ,	0	, 0 ,	0	, 0 },
        {0.2,   2.6,  0,  0, -2,  0, -1 , 0,	2 ,	0	, 0 ,	1	, 0 },
        {0,     0.2,  0,  0, -2,  0,  0 , 0,	2 ,	0	, 0 ,	2	, 0 },
        {0.1,   0.2,  1,  0, -2, -2, -2 , 0,	3 ,-2	, 1 ,	0	, 0 },
        {0.4,   1.1, -1,  0, -2,  0, -2 , 0,	3 ,	0	,-1 ,	0	, 0 },
        {0.2,   0.5, -1,  0, -2,  0, -1 , 0,	3 ,	0	,-1 ,	1	, 0 },
        {0.1,   0.2,  0,  0, -2, -2, -2 , 0,	4 ,-2	, 0 ,	0	, 0 },
        {0.1,   0.1, -2,  0, -2,  0, -2 , 0,	4 ,	0	,-2 ,	0	, 0 }
    };
    
    //  Legendre polynomial
    double GTidalCorrection::legendrePoly(int n, int m, double u)
    {
        // reference:Satellite Orbits Montenbruck. P66
        if(0==n && 0==m)
        {
            return 1.0;
        }
        else if(m==n)
        {
            return (2.0*m-1.0)*sqrt(1.0-u*u)*legendrePoly(n-1,m-1,u);
        }
        else if(n==m+1)
        {
            return (2.0*m+1)*u*legendrePoly(m,m,u);
        }
        else
        {
            return ((2.0*n-1.0)*u*legendrePoly(n-1,m,u)-(n+m-1.0)*legendrePoly(n-2,m,u))/(n-m);
        }
        
    }  // End of method 'EarthSolidTide::legendrePoly()'
    
    
    // Nnm IERS2003 P60
    double GTidalCorrection::normFactor(int n, int m)
    {
        // The input should be n >= m >= 0
        
        double fac(1.0);
        for( int i = (n-m+1); i <= (n+m); i++)
        {
            fac = fac * double(i);
        }
        
        double delta  = (m == 0) ? 1.0 : 0.0;
        
        double num = (2.0 * n + 1.0) * (2.0 - delta);
        
        // We should make sure fac!=0, but it won't happen on the case,
        // so we just skip handling it
        double out = sqrt(num/fac);
        
        return out;
        
    }  // End of method 'EarthSolidTide::normFactor'
    
    
    
    /*
     corrections for:
     C20, C21, C22, C30, C31, C32, C33, C40, C41, C42
     
     */
    void GTidalCorrection::getSolidEarthTideCorrection(double JD_UT1_I,double JD_UT1_F,double JD_TT, gfc::GVector &sunpos_ecef, gfc::GVector &moonpos_ecef, double dC[], double dS[])
    {
        
        //Table 6.3: Nominal values of solid Earth tide external potential Love numbers
        double LoveNumber[][5] = {  /* Knm,Knm+ (Elastic) | Re Knm , Im Knm,Knm+  (Anelastic)*/
                                    { 0.29525, -0.00087, 0.30190, -0.00000, -0.00089},  // n = 2, m = 0
                                    { 0.29470, -0.00079, 0.29830, -0.00144, -0.00080},  // n = 2, m = 1
                                    { 0.29801, -0.00057, 0.30102, -0.00130, -0.00057},  // n = 2, m = 2
                                    { 0.09300, -0.00000, 0.00000, -0.00000, -0.00000},  // n = 3, m = 0
                                    { 0.09300, -0.00000, 0.00000, -0.00000, -0.00000},  // n = 3, m = 1
                                    { 0.09300, -0.00000, 0.00000, -0.00000, -0.00000},  // n = 3, m = 2
                                    { 0.09400, -0.00000, 0.00000, -0.00000, -0.00000},  // n = 3, m = 3
            
                                 };
        
        static const double GM_Sun      = 1.3271250E11; //1.3271220e+20;    // [km^3/s^2]; STK
        static const double GM_Moon     = 4.9027890E3;         // // [km^3/s^2]; STK
        static const double R_Earth     = 6378.1370;   // kilometer
        static const double GM_Earth   = 3.986004415E5;

        static const double PI = GCONST("PI");
        
        
        double GMi[2] = {GM_Sun, GM_Moon};
        GVector planetPos[2] = {sunpos_ecef, moonpos_ecef};
        double lambda[2], phi[2];
        //N00, N10, N11, N20, N21, N22, N30, N31, N32, N33
        double Pnm[10]={0.0};
        for( int i = 0 ; i< 2; i++) // only two plannet, sun and moon
        {
            double gm = GMi[i]/GM_Earth;
            double reri = R_Earth/planetPos[i].norm();
            phi[i] = asin(planetPos[i].z/planetPos[i].norm());
            lambda[i] = atan2(planetPos[i].y,planetPos[i].x);
            //if(lambda[i]<0.0)  lambda[i]+= 2.0*PI;
                
            for(int n =2; n<=3; n++)
            {
                for(int m = 0 ; m<=n; m++)
                {
                    int index = n*(n+1)/2+m - 3;
                    double knm = LoveNumber[index][0];  //knm
                    double Nnm = normFactor( n, m );        //normalization coefficents of degree n and order m
                    // Pnm: normalized Legendre polynomials of degress n and order m
                    // 0 for sun and 1 for lunar each
                    double legendre =  Nnm*legendrePoly(n, m, sin(phi[i]));
                    double factor = 1.0/(2*n+1)*gm*pow(reri,n+1)*legendre;
                    
                    dC[index] += knm*factor*cos(m*lambda[i]);
                    dS[index] += knm*factor*sin(m*lambda[i]);
                    
                    if(n==2)
                    {
                        double knmp = LoveNumber[index][1];  //knm+
                        dC[7+m]+= knmp*factor*cos(m*lambda[i]); // delta C4m
                        dS[7+m]+= knmp*factor*sin(m*lambda[i]); // delta S4m
                    }
                }
                
            }
            
        }
        
        //step2:
        double BETA[6] = {0.0};
        double Dela[5] = {0.0};
        GIERS::DOODSONARG(JD_UT1_I, JD_UT1_F, JD_TT, Dela, BETA);
        
        //ref: IERS2010, equation 6.8a and table 6.5b
        //calculate theta_f
        for(int i = 0 ; i<21; i++ )
        {
            double theta_f = BETA[0]*Argu_C20[i][7] + BETA[1]*Argu_C20[i][8] + BETA[2]*Argu_C20[i][9]
                            +BETA[3]*Argu_C20[i][10] +BETA[4]*Argu_C20[i][11] +BETA[5]*Argu_C20[i][12];
            
            //thet_f = thet_f*PI/180.0;
            double t_s = std::sin(theta_f);
            double t_c = std::cos(theta_f);
            
            // Resulted from formula 6.8a in chapter 6.1 in IERS2010
            dC[0]+=((Argu_C20[i][0]*t_c-Argu_C20[i][1]*t_s)*1e-12);
        }
        
        //corrections for C21 and S21
        //ref: IERS2010, formula 6.8b and table 6.5a
        for( int i=0;i<48;i++)
        {
            // Computation of thet_f
            double theta_f = BETA[0]*Argu_C21[i][7] + BETA[1]*Argu_C21[i][8] + BETA[2]*Argu_C21[i][9]
                            +BETA[3]*Argu_C21[i][10] +BETA[4]*Argu_C21[i][11] +BETA[5]*Argu_C21[i][12];
            
            double t_s = std::sin(theta_f);
            double t_c = std::cos(theta_f);
            
            // Resulted from formula 5b in chapter 6.1
            dC[1] += ((Argu_C21[i][0]*t_s+Argu_C21[i][1]*t_c )*1e-12);
            dS[1] += ((Argu_C21[i][0]*t_c-Argu_C21[i][1]*t_s )*1e-12);
        }
        
        //corrections for C22, the fomular is exactly the same as C21, both equation 6.8b
        // but the table is different, using table 6.5c
        
        for( int i=0;i<2;i++)
        {
            // Input the computation of thet_f
            // Computation of thet_f
            double theta_f = BETA[0]*Argu_C22[i][7] + BETA[1]*Argu_C22[i][8] + BETA[2]*Argu_C22[i][9]
                            +BETA[3]*Argu_C22[i][10] +BETA[4]*Argu_C22[i][11] +BETA[5]*Argu_C22[i][12];
            
            
            //thet_f = thet_f*PI/180.0;
            double t_s = std::sin(theta_f);
            double t_c = std::cos(theta_f);
            
            // Resulted from formula 5b in chapter 6.1
            // The corrections are only to the real part.
            dC[2] += ((Argu_C22[i][0]*t_c)*1e-12 );
            dS[2] += -((Argu_C22[i][0]*t_s)*1e-12 );
        }

        //the third step
        
        
        
    }
    
    
    /*
       C20 C21 C22 C30 C31 C32 C33 C40 C41 C42
       double dC[10], double dS[10]
     * @param mjdUtc  UTC in MJD
     * @param dC      correction to normalized coefficients dC
     * @param dS      correction to normalized coefficients dS
     */
    void GTidalCorrection::getSolidEarthTideCorrection_1( double JD_UT1_I,double JD_UT1_F,double JD_TT, GVector& sunpos_ecef,GVector& moonpos_ecef, double dC[], double dS[])
    {
        
        static const double GM_Sun      = 1.3271250E11; //1.3271220e+20;    // [km^3/s^2]; STK
        static const double GM_Moon     = 4.9027890E3;         // // [km^3/s^2]; STK
        static const double R_Earth     = 6378.1370;   // kilometer
        static const double GM_Earth   = 3.986004415E5;
        static const double PI = GCONST("PI");
        
        double r_sun = sunpos_ecef.norm();
        double r_lunar = moonpos_ecef.norm();
        double GMi[2] ={GM_Sun, GM_Moon};
        double ri[2]  ={r_sun, r_lunar};
        
        
        
        double phi_sun =0.0, lamda_sun =0.0, phi_lunar =0.0, lamda_lunar =0.0;
        
        phi_sun = asin( sunpos_ecef.z / r_sun );
        lamda_sun = atan2(sunpos_ecef.y, sunpos_ecef.x);
        
        phi_lunar = asin(moonpos_ecef.z / r_lunar);
        lamda_lunar = atan2(moonpos_ecef.y, moonpos_ecef.x);
               
        // Euler's formula:  exp(ix) = cosx + isinx
        
        // reference bern 5 TIDPT2.f
        /*
         PERTURBING ACCELERATION DUE TO TIDES CORRESPONDING TO IERS STANDARDS 2003.
         STEP 1 CORRECTIONS OF SOLID EARTH TIDES INCLUDED,
         STEP 2 ONLY TERM DUE TO K1. SOLID EARTH POLE TIDES INCLUDED
         OCEAN TIDE TERMS UP TO N=M=4 INCLUDED
         */
        
        /*       IERS2003,  P60
         Elastic Earth           Anelastic Earth
         n m    knm     k+nm    Reknm   Imknm    k+nm
         2 0 0.29525 .0.00087 0.30190 .0.00000 .0.00089
         2 1 0.29470 .0.00079 0.29830 .0.00144 .0.00080
         2 2 0.29801 .0.00057 0.30102 .0.00130 .0.00057
         3 0 0.093 ?°Ë ?°Ë ?°Ë
         3 1 0.093 ?°Ë ?°Ë ?°Ë
         3 2 0.093 ?°Ë ?°Ë ?°Ë
         3 3 0.094 ?°Ë ?°Ë ?°Ë
         */
        complex<double> k[10] =      // Anelastic Earth
        {
            complex<double >(0.30190, 0.0),          // 20
            complex<double >(0.29830,-0.00144),      // 21
            complex<double >(0.30102,-0.00130),      // 22
            complex<double >(0.093, 0.0),            // 30
            complex<double >(0.093, 0.0),            // 31
            complex<double >(0.093, 0.0),            // 32
            complex<double >(0.094, 0.0),            // 33
            complex<double >(-0.00089, 0.0),         // k+ 20
            complex<double >(-0.00080, 0.0),         // k+ 21
            complex<double >(-0.00057, 0.0)          // k+ 22
        };
        
        complex<double> res[7];

        //----------------------------------------------------------------------
        // The first step of the computation ,refer to "IERS conventions 2003" P59
        // Each iteration for dC[n,m] and dS[n,m]
        for( int n=2;n<=3;n++ )
        {
            for( int m=0;m<=n;m++ )
            {
                int index = n * n - 2 * n + m;          //index in the returning value array
                
                double Nnm = normFactor( n, m );        //normalization coefficents of degree n and order m
                
                // Pnm: normalized Legendre polynomials of degress n and order m
                // 0 for sun and 1 for lunar each
                double sunPnm  = Nnm * legendrePoly( n, m, std::sin( phi_sun) );
                double moonPnm  = Nnm * legendrePoly( n, m, std::sin( phi_lunar) );
                
                double sunTemp = (GM_Sun/GM_Earth)*std::pow(R_Earth/r_sun,n+1) * sunPnm;
                double moonTemp = (GM_Moon/GM_Earth)*std::pow(R_Earth/r_lunar,n+1)*moonPnm;
                
                // Exp(-m*lamda*i) for sun and lunar each
                complex<double> c_sun   = complex<double>( std::cos( - m * lamda_sun ), std::sin( - m * lamda_sun ) );
                complex<double> c_lunar = complex<double>( std::cos( - m * lamda_lunar ), std::sin( - m * lamda_lunar ) );
                
                res[index] =  sunTemp * c_sun + moonTemp * c_lunar;
                
                dC[index]  =  (k[index]*res[index]).real()/(2.0*n+1.0);
                dS[index]  = -(k[index]*res[index]).imag()/(2.0*n+1.0);
                
            }  // 'for(int m=0;m<=n;m++)'
            
        }  // 'for(int n=2;n<=3;n++)'
        
        // The correction of dC[4,i] and dS[4,i](i=0,1,2) produced by degree 2 tide
        // The only difference from the above dC[2,i] and dS[2,i] is value of k replaced by k+
        for( int n = 0; n <= 2; n++ )
        {
            int index   = 2 * 2 - 2 * 2 + n;                     // liuwk
            complex<double> c_temp   = k[n+7 ] * res[ index ];   // liuwk
            dC[7+n] = c_temp.real() / 5.0;
            dS[7+n] =-c_temp.imag() / 5.0;
        }

        
       
       
        //   COMPUTE DOODSON'S FUNDAMENTAL ARGUMENTS (BETA)
        double BETA[6] = {0.0};
        double Dela[5] = {0.0};
        
        GIERS::DOODSONARG(JD_UT1_I, JD_UT1_F, JD_TT, Dela, BETA);
        double GMST = GIERS::getGMST(JD_UT1_I, JD_UT1_F, JD_TT);
        //here GMST should be in degres,
        
        for( int i=0;i<21;i++ )
        {
            // Input the computation of thet_f
            double thet_f = -( Argu_C20[i][2]*Dela[0]+Argu_C20[i][3]*Dela[1]+Argu_C20[i][4]*Dela[2]
                              + Argu_C20[i][5]*Dela[3]+Argu_C20[i][6]*Dela[4]);
            
            //thet_f = thet_f*PI/180.0;
            double t_s = std::sin(thet_f);
            double t_c = std::cos(thet_f);
            
            // Resulted from formula 5a in chapter 6.1
            // Modified, 05.12.2009
            //         dC[0] += ( ( Argu_C20[i][0] * t_c - Argu_C20[i][1] * t_s ) * 1e-12 );
            dC[0]+=((Argu_C20[i][0]*t_c-Argu_C20[i][1]*t_s)*1e-12);
        }

        
        
        //-------------------------------------------------------------
        // The second step
        // corrections for C21, S21 and C22 S22
        
        for( int i=0;i<48;i++)
        {
            // Computation of thet_f
            double thet_f = (GMST+PI)-(Argu_C21[i][2]*Dela[0]+Argu_C21[i][3]*Dela[1]+Argu_C21[i][4]*Dela[2]
                                                   + Argu_C21[i][5]*Dela[3]+Argu_C21[i][6]*Dela[4]);
            
            
            //thet_f = thet_f*PI/180.0;
            
            double t_s = std::sin(thet_f);
            double t_c = std::cos(thet_f);
            
            // Resulted from formula 5b in chapter 6.1
            dC[1] += ((Argu_C21[i][0]*t_s+Argu_C21[i][1]*t_c )*1e-12);
            dS[1] += ((Argu_C21[i][0]*t_c-Argu_C21[i][1]*t_s )*1e-12);
        }
        
//        printf("dC:\n");
//        for(int i = 0 ; i< 10; i++ )
//        {
//            printf("%E ", dC[i]);
//        }
//
//        printf("\ndS:\n");
//        for(int i = 0 ; i< 10; i++ )
//        {
//            printf("%E ", dS[i]);
//        }
        
        
        for( int i=0;i<2;i++)
        {
            // Input the computation of thet_f
            double thet_f = 2.0*(GMST+PI)-(Argu_C22[i][1]*Dela[0]+Argu_C22[i][2]*Dela[1]+Argu_C22[i][3]*Dela[2]
                                                     + Argu_C22[i][4]*Dela[3]+Argu_C22[i][5]*Dela[4]);
            //thet_f = thet_f*PI/180.0;
            double t_s = std::sin(thet_f);
            double t_c = std::cos(thet_f);
            
            // Resulted from formula 5b in chapter 6.1
            // The corrections are only to the real part.
            dC[2] += ((Argu_C22[i][0]*t_c)*1e-12 );
            dS[2] += ((-Argu_C22[i][0]*t_s)*1e-12 );
        }
        
        
        
        //--------------------------------------------------------------------------------
        // the third step
        //
        //C REMOVE PREMANENT TIDE FROM C02 (FOR JGM-3 NOT FOR GEM-T3)
        //C ---------------------------------------------------------
        
        
        
        
        
    } // end of function GTidalCorrection::getSolidEarthTideCorrection
    
    
    
    
    
    
    
} // end of namespace

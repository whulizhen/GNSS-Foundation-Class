//
//  GFMEarthGravity.cpp
//  GFC
//
//  Created by lizhen on 09/04/2017.
//  Copyright © 2017 lizhen. All rights reserved.
//

#include "GFMEarthGravity.hpp"

namespace gfc
{
    
    
    /*
     calculate the normalization facor for C and S
     
     Fnm = sqrt( (2n+1)(n-m)!/(n+m)!   )
     递推关系:
     Fn,m = Fn-1,m * sqrt( (n-m)(2n+1)/(n+m)/(2n-1)   )
     Fm,m = Fm-1,m-1 *sqrt( 2m+1/(2m-1)/(2m)/(2m-1)  )
     */

    void GFMEarthGravity::getNormalFactor()
    {
        
        int nn = m_desiredDegree + 2;
        int mm = m_desiredOrder  + 2;
        int size = (nn+1)*(nn+2)/2;
        double delta = 1.0;  // sqrt(2-delta0,m)
        N[0] = 1.0; // f00 = 1.0
        double factor = 0.0;
        
        //the case for m==0
        for( int n = 1; n<= m_nn; n++ )
        {
            factor = sqrt( (2.0*n+1.0)/(2.0*n-1.0)  );
            N[n*(n+1)/2] = N[(n-1)*n/2]*factor;
        }
        
        for(int m = 1; m<= m_mm; m++ )
        {
            factor = sqrt((2.0*m+1.0)/(2.0*m-1.0)/(2.0*m-1.0)/(2.0*m)  );
            N[m*(m+1)/2+m] = N[m*(m-1)/2+m-1]*factor;
            
            for( int n = m+1 ; n<= m_nn; n++ )
            {
                factor = sqrt( (2.0*n+1.0)*(n-m)/(2.0*n-1.0)/(n+m)  );
                N[n*(n+1)/2+m] = N[(n-1)*n/2+m]*factor;
            }
        }
        
        for(int m = 1; m<= m_mm; m++ )
        {
            for( int n = m ; n<= m_nn; n++ )
            {
                N[n*(n+1)/2+m] = N[n*(n+1)/2+m]*sqrt(2.0);
            }
        }
        
    }
    
    //denormalise the C and S
    void GFMEarthGravity::denormaliseCS(std::vector<double>& Cn, std::vector<double>& Sn)
    {
        getNormalFactor();
        
        for(int i = 0 ; i< N.size(); i++ )
        {
            C[i] = Cn[i]*N[i];
            S[i] = Sn[i]*N[i];
        }
        
        //copy C and S
        Cp = C;
        Sp = S;
        
    }
    
    
    void GFMEarthGravity::loadJGM3(gfc::GString filename)
    {
        
        
        
    }
    
    
    void GFMEarthGravity::loadEGM2008(gfc::GString filename)
    {
        //m_gravityCoefFilePath = filename;
        
        R = 6378.1363;  //km
        GM = 398600.4415;  // Km^3/s^2; 0.3986004418D15 value on website in m^3/s^2
        //update the GM constant
        
        
        //refMJD = 46431.0;
        //should be corrected to normalized C and S
        //dotC20 = 1.1627553400E-11;  // per year
        //dotC30 = 0.490000E-11;      // per year
        //dotC40 = 0.470000E-11;      // per year
        
        //https://www.iers.org/SharedDocs/Publikationen/EN/IERS/Publications/tn/TechnNote36/tn36.pdf?__blob=publicationFile&v=1
        // C21 and S21 should be calculated according to C20, C22 and S22 and mean polar motion
        //dotS21 = 1.6200000000E-11; // per year
        //dotC21 = -0.3200000000E-11;  // per year
        
        const int MAX = 1024;
        char store[MAX]={0};
        std::vector<double> C_norm(size,0.0 );
        std::vector<double> S_norm(size,0.0 );
        C_norm[0] = 1.0;
        
        //in the file, missing C10, C11 and C00 terms, they maybe zero ?
        
        int n =0; // degree
        int m =0; //order
        double cnm_normalised, snm_normalised;
        double dummy;
        std::ifstream infile(filename.c_str());
        if (infile.fail())
        {
            std::cerr << "\nCould not open gravity file. Terminating...\n";
            exit(0);
        }
        bool recorded = false;
        while ( !infile.eof() )
        {
            infile >> n; //degree of coefficient
            infile >> m; //order of coefficient
            infile >> cnm_normalised;
            infile >> snm_normalised;
            infile >> dummy; //not used
            infile >> dummy; // not used
            
            if ( (n <= m_nn) && (m <= m_mm) )
            {
                C_norm[n*(n+1)/2+m] = cnm_normalised;
                S_norm[n*(n+1)/2+m] = snm_normalised;
                recorded = true;
            }
            
            if( recorded == false ) // exit reading file
            {
                break;
            }
            
            recorded = false;
            
        } //end of loop reading over gravity field coefficients in file
        
        infile.close();
        
        denormaliseCS(C_norm, S_norm);
        
    }
    
    
    void GFMEarthGravity::generateVW(GVector& pos_ecef)
    {
        double r = pos_ecef.norm();
        double r2 = pos_ecef.norm2();
        

        
        V[0] = R/r;
        W[0] = 0.0;
        
        double xR = pos_ecef.x*R/r2;
        double yR = pos_ecef.y*R/r2;
        double zR = pos_ecef.z*R/r2;
        double rR = R*R/r2;
        
        //V[1] = zR * V[0];
        //W[1] = 0.0;
        
        //the case for m==0, n starts from 2
        for( int n = 1; n<= m_nn; n++ )
        {
            V[n*(n+1)/2] = ((2.0*n - 1.0) * zR * V[(n-1)*n/2] - (n - 1) * rR * V[ (n-2)*(n-1)/2 ]) /n;
            //W[n*(n+1)/2] = ((2.0*n - 1.0) * zR * W[(n-1)*n/2] - (n - 1) * rR * W[ (n-2)*(n-1)/2 ]) /n;
            W[n*(n+1)/2] = 0.0;
        }
        
        // Calculate tesseral and sectorial terms
        for (int m = 1; m <= m_mm; m++ )
        {
            // Calculate V(m,m) .. V(n_max+1,m)
            
            V[m*(m+1)/2+m] = (2 * m - 1) * ( xR * V[m*(m-1)/2+m-1] - yR * W[m*(m-1)/2+m-1] );
            W[m*(m+1)/2+m] = (2 * m - 1) * ( xR * W[m*(m-1)/2+m-1] + yR * V[m*(m-1)/2+m-1] );
            
            // Calculate V(m+1,m), when n=m+1
            int n = m+1;
            if( n<= m_nn )
            {
                V[n*(n+1)/2+m] = (2 * n - 1)*zR*V[n*(n-1)/2+m];
                W[n*(n+1)/2+m] = (2 * n - 1)*zR*W[n*(n-1)/2+m];
            }
            
            for( int n = m+2; n<= m_nn; n++ ) // start from m+2
            {
                V[n*(n+1)/2+m] = ((2*n-1.0)*zR*V[n*(n-1)/2+m] - (n+m-1.0)*rR*V[(n-2)*(n-1)/2+m])/(n-m);
                W[n*(n+1)/2+m] = ((2*n-1.0)*zR*W[n*(n-1)/2+m] - (n+m-1.0)*rR*W[(n-2)*(n-1)/2+m])/(n-m);
            }
            
        }  // End 'for (int m = 1; m <= (desiredOrder + 2); m++) '
        
    }
    
    /*
     
     calculate the earth gravity acceleration
     
     */
    GVector GFMEarthGravity::calculateAcc(gfc::GVector &pos_ecef)
    {
        GVector a_ecef;
        
        for (int m = 0; m <= m_desiredOrder; m++)
        {
            for (int n = m; n <= m_desiredDegree; n++)
            {
                
                double Cnm = C[n*(n+1)/2 + m];  // = C_n,m
                double Snm = S[n*(n+1)/2 + m]; // = S_n,m
                
                if (m==0)
                {
                    double Cn0 = C[n*(n+1)/2]; // = C_n,0
                    
                    a_ecef.x -= Cn0*V[(n+1)*(n+2)/2+1];
                    a_ecef.y -= Cn0*W[(n+1)*(n+2)/2+1];
                }
                else
                {
                    double Fac =  (n-m+1) * (n-m+2);
                    a_ecef.x += 0.5*((-Cnm*V[(n+1)*(n+2)/2+m+1] - Snm*W[(n+1)*(n+2)/2+m+1] ) + Fac*(Cnm*V[(n+1)*(n+2)/2+m-1] + Snm*W[(n+1)*(n+2)/2+m-1]));
                    a_ecef.y += 0.5*((-Cnm*W[(n+2)*(n+1)/2+m+1] + Snm*V[(n+2)*(n+1)/2+m+1]) + Fac*(-Cnm*W[(n+2)*(n+1)/2+m-1] + Snm*V[(n+1)*(n+2)/2+m-1]));
                }
                
                a_ecef.z += (n-m+1.0)*(-Cnm*V[(n+1)*(n+2)/2+m] - Snm*W[(n+1)*(n+2)/2+m]  );
                
            }  // End of 'for (int n = m; n <= (desiredDegree+1) ; n++)'
            
        }  // End of 'for (int m = 0; m <= (desiredOrder+1); m++)'
        
        
        a_ecef = a_ecef * ( GM / (R*R) )*1000.0;  // unit: m/s^2
        
        return a_ecef;
        
    }
    
    
    void GFMEarthGravity::getPartialDerivatives(GVector& pos_ecef)
    {
        double GMR3 = GM/(R*R*R);
        
        m_dadr.resize(3, 3);
        double paxpx =0.0, paxpy =0.0, paxpz =0.0, paypz =0.0, pazpz =0.0, paypy;
        int n =0, m = 0;
        int Inm=0,In1m=0, In1m1,In1m2, In1m_1, In1m_2, In0, In1, In12, In10, In11, In13;
        
        for ( n = m_desiredDegree; n >=0 ; --n )  // n==0 is a special situation
        {
            for ( m = n; m >= 0; --m )
            {
                
                
                In0 = n*(n+1)/2;
                In10 = (n+2)*(n+3)/2;
                In1 = In0+1;
                In11 = In10+1;
                In12 = In10+2;
                In13 = In10+3;
                
                Inm = In0+m;
                In1m =  In10+m;
                In1m1 = In1m+1;
                In1m2 = In1m1+1;
                In1m_1 = In1m-1;
                In1m_2 = In1m-2;
                
                
                //calculate paz/pz
                pazpz +=  (n-m+2)*(n-m+1)*( C[Inm]*V[In1m] + S[Inm]*W[In1m] );
                
                if( m == 0 )
                {
                    paxpx += 0.5*( C[In0]*V[In12] - (n+2)*(n+1)*C[In0]*V[In10] );
                    
                    //calculate pax/py
                    paxpy += 0.5*(C[In0]*W[In12]);
                    
                    // calculate pax/pz
                    paxpz += (n+1)*C[In0]*V[In11];
                    
                    //calculate pay/pz
                    paypz += (n+1)*C[In0]*W[In11];
                    
                }
                else if( m > 0 )
                {
                    paxpz +=  (n-m+1.0)*0.5*(C[Inm]*V[In1m1] + S[Inm]*W[In1m1])
                                      -0.5*(n-m+3)*(n-m+2)*(n-m+1)*( C[Inm]*V[In1m_1] + S[Inm]*W[In1m_1]);
                    
                    paypz +=   0.5*(n-m+1)*(C[Inm]*W[In1m1] - S[Inm]*V[In1m1])
                                      + 0.5*(n-m+3)*(n-m+2)*(n-m+1)*( C[Inm]*W[In1m_1] - S[Inm]*V[In1m_1] );
                    
                    if( m==1 )
                    {
                        paxpx +=  0.25*(C[In1]*V[In13] + S[In1]*W[In13])
                                        - n*(n+1)*( 3.0*C[In1]*V[In11] + S[In1]*W[In11]);
                        
                        paxpy += 0.25*(C[In1]*W[In13] - S[In1]*V[In13] )
                                        - n*(n+1)*( C[In1]*W[In11] + S[In1]*V[In11]);
                        
                    }
                    else if( m > 1)
                    {
                        paxpx +=  0.25*( C[Inm]*V[In1m2] + S[Inm]*W[In1m2] )
                                  -2.0*(n-m+2)*(n-m+1)*( C[Inm]*V[In1m] + S[Inm]*W[In1m] )
                                  +(n-m+4)*(n-m+3)*(n-m+2)*(n-m+1)*(C[Inm]*V[In1m_2] + S[Inm]*W[In1m_2]);
                        
                        paxpy += 0.25*( C[Inm]*W[In1m2] - S[Inm]*V[In1m2])
                                 +(n-m+4)*(n-m+3)*(n-m+2)*(n-m+1)*( -C[Inm]*W[In1m_2]+ S[Inm]*V[In1m_2]);
                    }
                }
            }
        }
        
        //format the matrix and then convert it to ECI frame
        paypy = -paxpx - pazpz;
        m_dadr(0,0) = paxpx;m_dadr(0,1) = paxpy;m_dadr(0,2) = paxpz;
        m_dadr(1,0) = paxpy;m_dadr(1,1) = paypy;m_dadr(1,2) = paypz;
        m_dadr(2,0) = paxpz;m_dadr(2,1) = paypz;m_dadr(2,2) = pazpz;
        m_dadr = m_dadr*GMR3;
        
        double tm[9] = {0.0};
        GSpaceEnv::eop.getECI2ECEFMatrix(tm);
        GMatrix U(tm,3,3);
        
        m_dadr = (~U)*m_dadr*(U);
        
        //dadr is in km
    }
    
    
    void GFMEarthGravity::correctTideCS(gfc::GTime ctUTC, double ut1mutc)
    {
        GTime UT1 =  GTime::UTC2UT1(ctUTC, ut1mutc);
        GTime TAI = GTime::UTC2TAI(ctUTC);
        GTime TT = GTime::TAI2TT(TAI);
        JDTime jdUT1  = GTime::GTime2JDTime(UT1);
        JDTime jdTT  = GTime::GTime2JDTime(TT);
        double JD_UT1_I = int(jdUT1.jdt());
        double JD_UT1_F = jdUT1.jdt() - JD_UT1_I;
        
        double N20 = GTidalCorrection::normFactor(2, 0);
        double N21 = GTidalCorrection::normFactor(2, 1);
        double N22 = GTidalCorrection::normFactor(2, 2);
        double N30 = GTidalCorrection::normFactor(3, 0);
        double N31 = GTidalCorrection::normFactor(3, 1);
        double N32 = GTidalCorrection::normFactor(3, 2);
        double N33 = GTidalCorrection::normFactor(3, 3);
        double N40 = GTidalCorrection::normFactor(4, 0);
        double N41 = GTidalCorrection::normFactor(4, 1);
        double N42 = GTidalCorrection::normFactor(4, 2);
        
        //corrections for the unnormalized coefficients
        double detC20 =0.0,detC21=0.0,detC22 =0.0,detC30=0.0, detC31=0.0, detC32=0.0, detC33=0.0,detC40=0.0,detC41=0.0,detC42=0.0;
        double detS20 =0.0,detS21=0.0,detS22 =0.0,detS30=0.0, detS31=0.0, detS32=0.0, detS33=0.0,detS40=0.0,detS41=0.0,detS42=0.0;
        
        //******** different correction form !!!!!!
        //C[2*(2+1)/2+1] = Cp[2*(2+1)/2+1] ; //C21, Cnm , C[n*(n+1)/2 + m]
        //S[2*(2+1)/2+1] = Sp[2*(2+1)/2+1] ;  //S21, Snm , S[n*(n+1)/2 + m]
        //unnormalized coefficients
        double C21 =  0.0, S21 = 0.0;
        
        if(time_variable == true)
        {
            double Cn20 =  Cp[2*(2+1)/2+0]/N20, Cn22 =  Cp[2*(2+1)/2+2]/N22, Sn22 = Sp[2*(2+1)/2+2]/N22;
            
            //时变重力场???
            //ref : https://www.iers.org/SharedDocs/Publikationen/EN/IERS/Publications/tn/TechnNote36/tn36.pdf?__blob=publicationFile&v=1
            // formular: 6.5
            double refMJD =  51544.5; // reference MJD J2000.0 in UTC
            double dotCn20 = 11.6E-12;
            double dotCn30 = 4.9E-12;
            double dotCn40 = 4.7E-12;
            
            //double dotC21 = -0.33700E-11;
            //double dotS21 = 1.60600E-11;
            
            double PI = GCONST("PI");
            //add the corrections to the denormalized C and S coefficients
            double leapYears = ( ctUTC.toDays() - refMJD )/365.25;
            double pxm,pym;
            GEarthOrientationParameter::getPM_mean(ctUTC.toDays(), pxm, pym);
            pxm = pxm * PI/(180.0*3600); // from arc sec to radian
            pym = pym * PI/(180.0*3600); // from arc sec to radian
            
            detC20 += N20*leapYears*dotCn20;
            detC30 += N30*leapYears*dotCn30;
            detC40 += N40*leapYears*dotCn40;
            
            //detC21 += N21*leapYears*dotC21;
            //detS21 += N21*leapYears*dotS21;
            
            //C[2*(2+1)/2+0] = Cp[2*(2+1)/2+0] +  detC20;  //C20, Cnm , C[n*(n+1)/2 + m]
            //C[3*(3+1)/2+0] = Cp[3*(3+1)/2+0] +  detC30;  //C30, Cnm , C[n*(n+1)/2 + m]
            //C[4*(4+1)/2+0] = Cp[4*(4+1)/2+0] +  detC40;  //C40, Cnm , C[n*(n+1)/2 + m]
            
           detC21 += (( sqrt(3.0)*pxm*Cn20- pxm*Cn22 + pym*Sn22 )*N21 - Cp[2*(2+1)/2+1]);
           detS21 += (( -sqrt(3.0)*pym*Cn20 - pym*Cn22 - pxm*Sn22 )*N21- Sp[2*(2+1)/2+1]);
            
            
        }
        
        GTidalCorrection tidecorrection;
        
        if( solid_tide == true )
        {
            // normalized coefficients C20 C21 C22 C30 C31 C32 C33 C40 C41 C42
            double dS[10]={0.0},dC[10] = {0.0}; //corrections for normalized C and S
            tidecorrection.getSolidEarthTideCorrection( JD_UT1_I, JD_UT1_F, jdTT.jdt(),
                                                       GSpaceEnv::planetPos_ecef[GJPLEPH::SUN],
                                                       GSpaceEnv::planetPos_ecef[GJPLEPH::MOON], dC, dS);
           
            detC20 +=N20*dC[0];
            detC21 +=N21*dC[1];
            detC22 +=N22*dC[2];
            detC30 +=N30*dC[3];
            detC31 +=N31*dC[4];
            detC32 +=N32*dC[5];
            detC33 +=N33*dC[6];
            detC40 +=N40*dC[7];
            detC41 +=N41*dC[8];
            detC42 +=N42*dC[9];
            
            detS20 +=N20*dS[0];
            detS21 +=N21*dS[1];
            detS22 +=N22*dS[2];
            detS30 +=N30*dS[3];
            detS31 +=N31*dS[4];
            detS32 +=N32*dS[5];
            detS33 +=N33*dS[6];
            detS40 +=N40*dS[7];
            detS41 +=N41*dS[8];
            detS42 +=N42*dS[9];
            
        }
        
        if( polar_tide )  // pole motion tide correction
        {
            double dC21 =0.0, dS21 = 0.0;
            
            GSpaceEnv::eop.getPoleTide(dC21, dS21);
            
            detC21 += N21*dC21;
            detS21 += N21*dS21;
            
        }
        
        
            C[2*(2+1)/2+0] = Cp[2*(2+1)/2+0] + detC20;
            C[2*(2+1)/2+1] = Cp[2*(2+1)/2+1] + detC21;
            C[2*(2+1)/2+2] = Cp[2*(2+1)/2+2] + detC22;
            C[3*(3+1)/2+0] = Cp[3*(3+1)/2+0] + detC30;
            C[3*(3+1)/2+1] = Cp[3*(3+1)/2+1] + detC31;
            C[3*(3+1)/2+2] = Cp[3*(3+1)/2+2] + detC32;
            C[3*(3+1)/2+3] = Cp[3*(3+1)/2+3] + detC33;
            C[4*(4+1)/2+0] = Cp[4*(4+1)/2+0] + detC40;
            C[4*(4+1)/2+1] = Cp[4*(4+1)/2+1] + detC41;
            C[4*(4+1)/2+2] = Cp[4*(4+1)/2+2] + detC42;
        
            S[2*(2+1)/2+0] = Sp[2*(2+1)/2+0] + detS20; // S20
            S[2*(2+1)/2+1] = Sp[2*(2+1)/2+1] + detS21;
        
            S[2*(2+1)/2+2] = Sp[2*(2+1)/2+2] + detS22;
            S[3*(3+1)/2+0] = Sp[3*(3+1)/2+0] + detS30; // S30
            S[3*(3+1)/2+1] = Sp[3*(3+1)/2+1] + detS31;
            S[3*(3+1)/2+2] = Sp[3*(3+1)/2+2] + detS32;
            S[3*(3+1)/2+3] = Sp[3*(3+1)/2+3] + detS33;
            S[4*(4+1)/2+0] = Sp[4*(4+1)/2+0] + detS40; // S40
            S[4*(4+1)/2+1] = Sp[4*(4+1)/2+1] + detS41;
            S[4*(4+1)/2+2] = Sp[4*(4+1)/2+2] + detS42;
        
        
    }
    
    // must implement this function
    void GFMEarthGravity::doCompute(GTime ctUTC,double ut1mutc, GVector& pos_ecef  )
    {
        GVector a_ecef, a_eci;
        
        
        if(m_desiredDegree >=4)
        {
            correctTideCS(ctUTC, ut1mutc);
        }
        
        generateVW(pos_ecef);
        
        a_ecef = calculateAcc(pos_ecef);
        
        GSpaceEnv::eop.ECEF2ECI_pos(a_ecef, a_eci);
        
        setForce(a_eci);
        
        if( m_hasPartialDerivatives == true)
        {
            getPartialDerivatives(pos_ecef);
        }
        
    }
    
    
    
    
} // end of namespace

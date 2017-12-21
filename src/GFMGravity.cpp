/*! @file Force_earth_gravity.cpp
	@author David Harrison
	@date 16 May 2016
	@brief SGNL OPS implementation of Earth gravity.
 */

#include "GFMGravity.hpp"
namespace gfc
{
    
    
    void GFMGravity::loadGRACE_ggm03c(GString filename)
    {
        
        m_R = 6378.1363;   // 6378136.30 value in file in m
        m_GM = 398600.4415; // 398600.44150E+09 value in file in m^3/s^2;  km^3/s^2
        
        m_gravityCoefFilePath = filename;
        const int MAX = 100;
        char store[MAX]={0};
        std::vector<std::vector<double> > C_norm( m_N + 1, std::vector<double>(m_M + 1, 0.0));
        std::vector<std::vector<double> > S_norm( m_N + 1, std::vector<double>(m_M + 1, 0.0));
        C_norm[0][0] = 1.0;
        bool recorded = false;
        int count = 0, i = 0;
        int n =0; // degree
        int m =0; //order
        double cnm_normalised = 0.0, snm_normalised = 0.0;
        double dummy = 0.0;
        char temp = 0;
        std::ifstream infile(m_gravityCoefFilePath.c_str());
        if (infile.fail())
        {
            std::cerr << "\nCould not open gravity file. Terminating...\n";
            exit(0);
        }

        //std::cout << __FILE__ << "\t" << __LINE__ << std::endl;
        
        //read past the file header
        for (count = 1; count <= 3; count++)
        {
            infile.getline(store, MAX);
        }
        //next line - the normalised gravity field coefficients start
        while (!infile.eof())
        {
            //get past the first 6 characters
            for (count = 1; count <= 6; count++) infile >> temp;
            infile >> n; //degree of coefficient
            
            if (n > 360) //fix to deal with degree/order looking like 360100
            {
                i = n;
                n /= 1000;
                m = i - n * 1000;
            }
            else
            {
                infile >> m; //order of coefficient
            }
            
            infile >> cnm_normalised;
            infile >> snm_normalised;
            infile >> dummy;
            infile >> dummy;
            infile >> dummy;
            
            if ((n <= m_N) && (m <= m_M))
            {
                C_norm[n][m] = cnm_normalised;
                S_norm[n][m] = snm_normalised;
                recorded = true;
            }
            
            
            if( recorded == false ) // exit reading file
            {
                break;
            }
            
            recorded = false;
            
        } //end of loop reading over gravity field coefficients in file
        
        
        infile.close();
        
        denorm_coef(C_norm, S_norm);
        
        generate_h_coefs();
        
    }
    
    
    
    void GFMGravity::loadEGM2008(gfc::GString filename)
    {
        m_gravityCoefFilePath = filename;
        
        m_R = 6378.1363;  //km
        m_GM = 398600.4415;  // Km^3/s^2; 0.3986004418D15 value on website in m^3/s^2
        
        refMJD = 46431.0;
        //should be corrected to normalized C and S
        dotC20 = 1.1627553400E-11;
        dotC30 = 0.490000E-11;
        dotC40 = 0.470000E-11;
        
        dotS21 = 1.6200000000E-11;
        dotC21 = -0.3200000000E-11;
        
        //m_R = 6378.1363;  //km
        //m_GM = 398600.4415;  // Km^3/s^2; 0.3986004418D15 value on website in m^3/s^2
        const int MAX = 1024;
        char store[MAX]={0};
        std::vector<std::vector<double> > C_norm( m_N + 1, std::vector<double>(m_M + 1, 0.0));
        std::vector<std::vector<double> > S_norm( m_N + 1, std::vector<double>(m_M + 1, 0.0));
        C_norm[0][0] = 1.0;
        int n =0; // degree
        int m =0; //order
        double cnm_normalised, snm_normalised;
        double dummy;
        std::ifstream infile(m_gravityCoefFilePath.c_str());
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
            
            if ( (n <= m_N) && (m <= m_M) )
            {
                C_norm[n][m] = cnm_normalised;
                S_norm[n][m] = snm_normalised;
                recorded = true;
            }
            
            if( recorded == false ) // exit reading file
            {
                break;
            }
            
            recorded = false;
            
        } //end of loop reading over gravity field coefficients in file
        
        infile.close();
        
        denorm_coef(C_norm, S_norm);
        
        generate_h_coefs();
        
    }
    
     // http://earth-info.nga.mil/GandG/wgs84/gravitymod/egm96/egm96.html
    void GFMGravity::loadEGM96(gfc::GString filename)
    {
        
        m_gravityCoefFilePath = filename;
        
        m_R = 6378.1363;  //km
        m_GM = 398600.4415;  // Km^3/s^2; 0.3986004418D15 value on website in m^3/s^2
        
        refMJD = 46431.0;
        dotC20 = 1.1627553400E-11;
        dotC21 = -0.3200000000E-11;
        dotS21 = 1.6200000000E-11;
        
        
        const int MAX = 100;
        char store[MAX]={0};
        std::vector<std::vector<double> > C_norm( m_N + 1, std::vector<double>(m_M + 1, 0.0));
        std::vector<std::vector<double> > S_norm( m_N + 1, std::vector<double>(m_M + 1, 0.0));
        C_norm[0][0] = 1.0;
        int n =0; // degree
        int m =0; //order
        double cnm_normalised, snm_normalised;
        double dummy;
        std::ifstream infile(m_gravityCoefFilePath.c_str());
        if (infile.fail())
        {
            std::cerr << "\nCould not open gravity file. Terminating...\n";
            exit(0);
        }
        bool recorded = false;
        infile.getline(store, MAX); //read past the file header
        while ( !infile.eof() )
        {
            infile >> n; //degree of coefficient
            infile >> m; //order of coefficient
            infile >> cnm_normalised;
            infile >> snm_normalised;
            infile >> dummy; //not used
            infile >> dummy; // not used
            
            if ( (n <= m_N) && (m <= m_M) )
            {
                C_norm[n][m] = cnm_normalised;
                S_norm[n][m] = snm_normalised;
                recorded = true;
            }
            
            if( recorded == false ) // exit reading file
            {
                break;
            }
            
            recorded = false;
            
        } //end of loop reading over gravity field coefficients in file
        
        infile.close();
        
        denorm_coef(C_norm, S_norm);
        
        generate_h_coefs();
        
    }
    
    /// Calculate V' and W' for current position.
    // ecef coordinate should be in meters, Not Kilo meters
    // V and W are different in different epoch time and positions
    void GFMGravity::generate_VW_prime(double* ecef)
    {
        double pos[3]={0.0};  //unit: km
        memcpy(pos, ecef, sizeof(double)*3);
        
        double r2 = (pos[0]*pos[0] + pos[1]*pos[1] + pos[2]*pos[2] ); // unit:km
        
        double r = sqrt(r2);  //unit: km
        
        int n = 0, m = 0;
        int nrow1 = 0, nrow2 = 0, nrow3 = 0, mrow1 = 0, mrow2 = 0;
        //double Rer2 = m_a / state_var.r2;
        double Rer2 = m_R / r2;
        
        double x =0.0, y=0.0, z=0.0;
        
        x = static_cast<double>(pos[0]);
        y = static_cast<double>(pos[1]);
        z = static_cast<double>(pos[2]);
        
        //V[0] = Re / state_var.r;
        V[0] = m_R / r;
        W[0] = 0.0;
        
        V[1] = Rer2 * z * V[0];
        W[1] = 0.0;
        
        nrow1 = 1;
        nrow2 = 1;
        nrow3 = 0;
        
        // Compute the Vn0 terms (all Wn0 terms are 0)
        for ( n = 2; n <= m_N + 2 ; ++n )
        {
            nrow1 += n;
            V[nrow1] = Rer2 * (z * V[nrow2] - h[nrow1] * V[nrow3]);
            nrow3 = nrow2;
            nrow2 = nrow1;
            //std::cout << n << 0 << ", " << V[n][0] << ", " << W[n][0] << std::endl;
        }
        
        mrow1 = 0;
        mrow2 = 0;
        for ( m = 1; m <= m_M + 2; ++m )
        {
            // Compute the Vmm and Wmm terms
            
            mrow1 += m + 1;
            V[mrow1] = Rer2 * (x * V[mrow2] - y * W[mrow2]);
            W[mrow1] = Rer2 * (x * W[mrow2] + y * V[mrow2]);
            //V[m*(m+1)/2+m] = Rer2 * (x*V[(m-1)*m/2+m-1] - y*W[(m-1)*m/2+m-1]);
            //W[m*(m+1)/2+m] = Rer2 * (x*W[(m-1)*m/2+m-1] + y*V[(m-1)*m/2+m-1]);
            
            mrow2 = mrow1;
            
            n = m + 1;
            nrow1 = mrow1 + n;
            nrow2 = nrow1;
            nrow3 = mrow1;
            
            // Compute Vm+1,m and Wm+1,m terms
            if (n <= m_N + 2)
            {
                V[nrow1] = Rer2 * z * V[nrow3];
                W[nrow1] = Rer2 * z * W[nrow3];
                //V[n*(n+1)/2+m] = Rer2 * z*V[(n-1)*n/2+m];
                //W[n*(n+1)/2+m] = Rer2 * z*W[(n-1)*n/2+m]; field l
                for ( n = m + 2; n <= m_N + 1; ++n )
                {
                    // Compute the Vnm and Wnm terms
                    
                    nrow1 += n;
                    
                    V[nrow1] = Rer2 * (z * V[nrow2] - h[nrow1] * V[nrow3]);
                    W[nrow1] = Rer2 * (z * W[nrow2] - h[nrow1] * W[nrow3]);
                    //V[n*(n+1)/2+m] = Rer2 * (z*V[(n-1)*n/2+m] - h[n*(n+1)/2+m]*V[(n-2)*(n-1)/2+m]);
                    //W[n*(n+1)/2+m] = Rer2 * (z*W[(n-1)*n/2+m] - h[n*(n+1)/2+m]*W[(n-2)*(n-1)/2+m]);
                    
                    nrow3 = nrow2;
                    nrow2 = nrow1;
                }
            }
        }
    } // end of function populate_VW_prime
    
    
    void GFMGravity::generate_h_coefs()
    {
        int n = 0, m = 0;
        double n_d=0.0, m_d=0.0;
        for (m = 0; m <= m_M + 1; m++)
        {
            
            n_d = m_d + 2.0;
            for ( n = m + 2; n <= m_N + 1; n++ )
            {
                
                h[n * (n + 1) / 2 + m] =
                ((n_d + m_d - 1.0) * (n_d - m_d - 1.0) * m_R) /
                ((2.0 * n_d - 1.0) * (2.0 * n_d - 3.0));
                
                n_d += 1.0;
            }
            
            m_d += 1.0;
        }
        
    } // End of function generate_VW_coefs
    
    
    std::vector<std::vector<double>> GFMGravity::generate_denorm_factors( int max_n, int max_m) const
    {
        int n = 0, m = 0;
        double n_d = 0, m_d = 0;
        
        std::vector<std::vector<double>> g(max_n + 1, std::vector<double>(max_m + 1));
        
        n_d = 0.0;
        m_d = 0.0;
        
        g[0][0] = 0.25;
        
        if (max_n > 0 && max_m > 0)
        {
            g[1][1] = 0.6L;
        }
        
        for ( m = 0; m <= max_m; m++)
        {
            if ( m > 1 )
            {
                g[m][m] = ((8.0 * m_d + 4.0) * g[m - 1][m - 1]) / (10.0 * m_d);
            }
            
            n_d = m_d + 1.0;
            for (n = m + 1; n <= max_n; n++)
            {
                g[n][m] = (2.0 * n_d - 1.0) * (2.0 * n_d + 1.0) * g[n - 1][m] /
                (5.0 * (n_d - m_d) * (n_d + m_d));
                n_d += 1.0;
            }
            m_d += 1.0;
        }
        
        double n_corr = 1.0;
        double m_corr = 1.0;
        int mod4 = 0;
        
        double n_mult[] = {5.0, 1.0, 1.0, 1.0};
        double g_mult[] = {1.0, 5.0, 25.0, 125.0};
        
        for (m = 0; m <= max_m; m++)
        {
            
            g[m][m] = sqrt(g[m][m] * m_corr);
            
            n_corr = 1.0;
            
            for (n = m + 1; n <= max_n; n++)
            {
                
                mod4 = (n - m) & 3;
                
                n_corr *= n_mult[mod4];
                
                g[n][m] = (sqrt(g[n][m] * g_mult[mod4] * m_corr) * n_corr) * n_corr;
            }
            
            m_corr *= 1.25;
        }
        
        return g;
        
    } // End of function generate_denorm_factors
    
    
    
    /// Populates C and S arrays with enhanced denormalised values (C' and S').
    void GFMGravity::denorm_coef( const std::vector<std::vector<double>> &C_norm, const std::vector<std::vector<double>> &S_norm)
    {
        int n = 0, m = 0;
        
        std::vector<std::vector<double>> g = generate_denorm_factors(m_N, m_M);
        
        for ( m = 0; m <= m_M; m++ )
        {
            for ( n = m; n <= m_N; n++ )
            {
                C[n * (n + 1) / 2 + m] = C_norm[n][m] * g[n][m] * (2 * n + 1);
                S[n * (n + 1) / 2 + m] = S_norm[n][m] * g[n][m] * (2 * n + 1);
            }
        }
        
        //make a copy of C and S (unnormalized coefficients)
        Cp = C;
        Sp = S;
        
    } // End of function denorm_coef
    
    
    // the general relativity correction for the gravity force
    // a :  in ECEF, expected to be in m/s^2
    // reference: satellite orbit, page 111
    void GFMGravity::generalRelativityCorrection( double* pv_eci, double* a )
    {
        double c = GCONST("CLIGHT")/1000.0; //unit: km
        double C2 = c*c;
        double minusGMc2 = -m_GM/C2;  // the unit of m_GM is km^3/s^2 ; the unit of minusGMc2 is meter
        double r2 = pow(pv_eci[0],2.0) + pow(pv_eci[1],2.0) + pow(pv_eci[2],2.0);
        double v2 = pow(pv_eci[3],2.0) + pow(pv_eci[4],2.0) + pow(pv_eci[5],2.0);
        double r = sqrt(r2);
        double minusGMc2r3 = minusGMc2 / (r * r2);
        double pvd =0.0; // unit: km^2
        for( int i = 0 ; i<3; i++ )
        {
            pvd += pv_eci[i]*pv_eci[i+3];
        }
        
        double r_coef = minusGMc2r3 * (4.0 * m_GM / r - v2); //unit: 1.0/s^4
        double v_coef = minusGMc2r3 * 4.0 * pvd;
        
        //double a[3] = {0.0}; // in ECEF, expected to be in km/s^2
        
        a[0] = (r_coef * pv_eci[0] + v_coef * pv_eci[3])*1000.0;
        a[1] = (r_coef * pv_eci[1] + v_coef * pv_eci[4])*1000.0;
        a[2] = (r_coef * pv_eci[2] + v_coef * pv_eci[5])*1000.0;
        
    }
    
    
    /*
     * get the partial derivatives of gravity force
     *
     *
     */
    void GFMGravity::getPartialDerivatives(GVector& pos_ecef)
    {
        m_dadr.resize(3, 3);
        
        
        int n = 0 , m = 0;
        int nrow1 = 0, nrow2 = 0, nrow3 = 0;
        
        n = m_N; //degree
        
        double GMR3 = m_GM/(m_R*m_R*m_R);
        
        double paxpx = 0.0, paxpy =0.0, paxpz =0.0, paypz =0.0,pazpz =0.0;
        double paypy = 0.0;
        int index_nm = 0;
        
        //double t[11] = {0.0};
        
        //sum up all the degree and order
        int degree = m_N;  // should be m_N
        
        for ( n = degree; n >=0 ; --n )  // n==0 is a special situation
        {
            for ( m = n; m >= 0; --m )
            {
                //calculate paz/pz
                pazpz += 2.0*(2*n+3)*( C[n*(n+1)/2+m]*V[(n+2)*(n+3)/2+m] + S[n*(n+1)/2+m]*W[(n+2)*(n+3)/2+m] );
                
                if( m == 0 )
                {
                    paxpx += (2*n+3)*( C[n*(n+1)/2]*V[(n+2)*(n+3)/2+2] - C[n*(n+1)/2]*V[(n+2)*(n+3)/2] );
                    //calculate pax/px
                    //paxpx += t[1];
                    
                    //calculate pax/py
                    paxpy += (2*n+3)*(C[n*(n+1)/2]*W[(n+2)*(n+3)/2+2]);
                    //paxpy += t[2];
                    
                    // calculate pax/pz
                    paxpz += 2.0*(2*n+3)*C[n*(n+1)/2]*V[(n+2)*(n+3)/2+1];
                    //paxpz += t[3];
                    
                    //calculate pay/pz
                    paypz += 2.0*(2*n+3)*C[n*(n+1)/2]*W[(n+2)*(n+3)/2+1];
                    //paypz += t[4];
                    
                }
                else if( m > 0 )
                {
                    paxpz += (2*n+3)*( C[n*(n+1)/2+m]*V[(n+2)*(n+3)/2+m+1] + S[n*(n+1)/2+m]*W[(n+2)*(n+3)/2+m+1]
                                    -C[n*(n+1)/2+m]*V[(n+2)*(n+3)/2+m-1] - S[n*(n+1)/2+m]*W[(n+2)*(n+3)/2+m-1]);
                    //paxpz += t[5];
                    
                    paypz += (2*n+3)*(C[n*(n+1)/2+m]*W[(n+2)*(n+3)/2+m+1] - S[n*(n+1)/2+m]*V[(n+2)*(n+3)/2+m+1]
                                    +C[n*(n+1)/2+m]*W[(n+2)*(n+3)/2+m-1] - S[n*(n+1)/2+m]*V[(n+2)*(n+3)/2+m-1]);
                    //paypz += t[6];
                    
                    if( m==1 )
                    {
                        paxpx += 0.5*(2*n+3)*( C[n*(n+1)/2+1]*V[(n+2)*(n+3)/2+3] + S[n*(n+1)/2+1]*W[(n+2)*(n+3)/2+3]
                                            -C[n*(n+1)/2+1]*V[(n+2)*(n+3)/2+1] - S[n*(n+1)/2+1]*W[(n+2)*(n+3)/2+1]);
                        //paxpx += t[7];
                        
                        paxpy += 0.5*(2*n+3)*( C[n*(n+1)/2+1]*W[(n+2)*(n+3)/2+3] - S[n*(n+1)/2+1]*V[(n+2)*(n+3)/2+3]
                                            -C[n*(n+1)/2+1]*W[(n+2)*(n+3)/2+1] - S[n*(n+1)/2+1]*V[(n+2)*(n+3)/2+1]);
                        //paxpy += t[8];
                        
                    }
                    else if( m > 1 )
                    {
                        paxpx +=  0.5*(2*n+3)*( C[n*(n+1)/2+m]*V[(n+2)*(n+3)/2+m+2] + S[n*(n+1)/2+m]*W[(n+2)*(n+3)/2+m+2]
                                            -2.0*( C[n*(n+1)/2+m]*V[(n+2)*(n+3)/2+m] + S[n*(n+1)/2+m]*W[(n+2)*(n+3)/2+m] )
                                            + C[n*(n+1)/2+m]*V[(n+2)*(n+3)/2+m-2] + S[n*(n+1)/2+m]*W[(n+2)*(n+3)/2+m-2]
                                            );
                        
                        paxpy += 0.5*(2*n+3)*( C[n*(n+1)/2+m]*W[(n+2)*(n+3)/2+m+2] - S[n*(n+1)/2+m]*V[(n+2)*(n+3)/2+m+2]
                                             -C[n*(n+1)/2+m]*W[(n+2)*(n+3)/2+m-2]+ S[n*(n+1)/2+m]*V[(n+2)*(n+3)/2+m-2]
                                             );
                        
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
        
//        //for two body problem
//        double r2 = pos_ecef.norm2();
//        m_dadr(0,0) = 3.0*pos_ecef.x*pos_ecef.x - r2;
//        m_dadr(0,1) = 3.0*pos_ecef.x*pos_ecef.y;
//        m_dadr(0,2) = 3.0*pos_ecef.x*pos_ecef.z;
//        m_dadr(1,0) = 3.0*pos_ecef.x*pos_ecef.y;
//        m_dadr(1,1) = 3.0*pos_ecef.y*pos_ecef.y - r2;
//        m_dadr(1,2) = 3.0*pos_ecef.y*pos_ecef.z;
//        m_dadr(2,0) = 3.0*pos_ecef.x*pos_ecef.z;
//        m_dadr(2,1) = 3.0*pos_ecef.y*pos_ecef.z;
//        m_dadr(2,2) = 3.0*pos_ecef.z*pos_ecef.z - r2;
//        m_dadr = m_dadr*m_GM/(r2*r2*sqrt(r2));
        
     //   printf("origin dadr:\n");
     //   m_dadr.print();
        
        // then transfer the ECEF to ECI
        // according to Montenbruck Satellite Orbit, Page 246, the law of matrix differential
        // dady(ECI) = inv(U)*dady(ECEF)*U , U is the matrix from ECI to ECEF
        double tm[9] = {0.0};
        GSpaceEnv::eop.getECI2ECEFMatrix(tm);
        GMatrix U(tm,3,3);
        
        m_dadr = (~U)*m_dadr*(U);
        
    }  // end of function getPartialDerivatives
    
    
    /*
     corrections to the denormalized coefficients
     */
    void GFMGravity::correctTideCS(GTime ctUTC,double ut1mutc,bool solidFlag, bool poleFlag, bool oceanFlag)
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
        
        //add the corrections to the denormalized C and S coefficients
        double leapYears = ( ctUTC.toDays() - refMJD )/365.25;
        double detC20 = N20*leapYears*dotC20;
        double detC21 = N21*leapYears*dotC21;
        double detS21 = N21*leapYears*dotS21;
        
        //时变重力场???
        C[2*(2+1)/2+0] = Cp[2*(2+1)/2+0] +  detC20;  //C20, Cnm , C[n*(n+1)/2 + m]
        C[2*(2+1)/2+1] = Cp[2*(2+1)/2+1] +  detC21 ; //C21, Cnm , C[n*(n+1)/2 + m]
        S[2*(2+1)/2+1] = Sp[2*(2+1)/2+1] +  detS21;  //S21, Snm , S[n*(n+1)/2 + m]
        
//        C[2*(2+1)/2+0] = C[2*(2+1)/2+0] + detC20;  //C20, Cnm , C[n*(n+1)/2 + m]
//        C[2*(2+1)/2+1] = C[2*(2+1)/2+1] + detC21 ; //C21, Cnm , C[n*(n+1)/2 + m]
//        S[2*(2+1)/2+1] = S[2*(2+1)/2+1] + detS21;  //S21, Snm , S[n*(n+1)/2 + m]
        
        
        GTidalCorrection tidecorrection;
        
        
        if( solidFlag ) // solid earth tide correction
        {
            
            
            // C20 C21 C22 C30 C31 C32 C33 C40 C41 C42
            double dS[10]={0.0},dC[10] = {0.0}; //corrections for C and S
            tidecorrection.getSolidEarthTideCorrection( JD_UT1_I, JD_UT1_F, jdTT.jdt(),
                                                       GSpaceEnv::planetPos_ecef[GJPLEPH::SUN],
                                                       GSpaceEnv::planetPos_ecef[GJPLEPH::MOON], dC, dS);
            
            //C20, Cnm , C[n*(n+1)/2 + m]
            C[2*(2+1)/2+0] = Cp[2*(2+1)/2+0] + N20*dC[0];
            C[2*(2+1)/2+1] = Cp[2*(2+1)/2+1] + N21*dC[1];
            C[2*(2+1)/2+2] = Cp[2*(2+1)/2+2] + N22*dC[2];
            C[3*(3+1)/2+0] = Cp[3*(3+1)/2+0] + N30*dC[3];
            C[3*(3+1)/2+1] = Cp[3*(3+1)/2+1] + N31*dC[4];
            C[3*(3+1)/2+2] = Cp[3*(3+1)/2+2] + N32*dC[5];
            C[3*(3+1)/2+3] = Cp[3*(3+1)/2+3] + N33*dC[6];
            C[4*(4+1)/2+0] = Cp[4*(4+1)/2+0] + N40*dC[7];
            C[4*(4+1)/2+1] = Cp[4*(4+1)/2+1] + N41*dC[8];
            C[4*(4+1)/2+2] = Cp[4*(4+1)/2+2] + N42*dC[9];
            
            S[2*(2+1)/2+0] = Sp[2*(2+1)/2+0] + N20*dS[0]; // S20
            S[2*(2+1)/2+1] = Sp[2*(2+1)/2+1] + N21*dS[1];
            S[2*(2+1)/2+2] = Sp[2*(2+1)/2+2] + N22*dS[2];
            S[3*(3+1)/2+0] = Sp[3*(3+1)/2+0] + N30*dS[3]; // S30
            S[3*(3+1)/2+1] = Sp[3*(3+1)/2+1] + N31*dS[4];
            S[3*(3+1)/2+2] = Sp[3*(3+1)/2+2] + N32*dS[5];
            S[3*(3+1)/2+3] = Sp[3*(3+1)/2+3] + N33*dS[6];
            S[4*(4+1)/2+0] = Sp[4*(4+1)/2+0] + N40*dS[7]; // S40
            S[4*(4+1)/2+1] = Sp[4*(4+1)/2+1] + N41*dS[8];
            S[4*(4+1)/2+2] = Sp[4*(4+1)/2+2] + N42*dS[9];
            
        }
        
        if( poleFlag )  // pole motion tide correction
        {
            double dC21 =0.0, dS21 = 0.0;
            
            GSpaceEnv::eop.getPoleTide(dC21, dS21);
            
            C[2*(2+1)/2+1] = Cp[2*(2+1)/2+1] + N21*dC21;
            S[2*(2+1)/2+1] = Sp[2*(2+1)/2+1] + N21*dS21;
            
        }
        
        if( oceanFlag )  // ocean tide correction
        {
            
            
        }
        
    }  // end of function correctTideCS
    
    
   
    /*
     
     gravity force model with solid earth and pole tide correction
     
     */
    void  GFMGravity::doCompute( GTime ctUTC,double ut1mutc, GVector& pos_ecef  )
    {
        
        if(m_N > 2)
        {
            correctTideCS(ctUTC, ut1mutc,true, true, true);
        }
        
        int n = 0 , m = 0;
        double twice_cn0 = 0.0;
        int nrow1 = 0, nrow2 = 0, nrow3 = 0;
        
        double ecefPosition[3] ={ pos_ecef.x , pos_ecef.y, pos_ecef.z };
        
        n = m_N;
        
        generate_VW_prime(ecefPosition);
        GVector a, a_eci;
        
        // Contributions are summed from smallest to largest to preserve precision
        
        nrow1 = n * (n + 1) / 2 + n;
        nrow2 = (n + 1) * (n + 2) / 2 + n + 1;
        nrow3 = nrow2 - 2;
        
        for (n = m_N; n > 1; --n)
        {
            
            //nrow1 = n*(n+1)/2+n;
            //nrow2 = (n+1)*(n+2)/2+n+1;
            //nrow3 = nrow2-2;
            
            for (m = n; m > 0; --m)
            {
                a.x += C[nrow1] * (V[nrow2] - V[nrow3]) +
                S[nrow1] * (W[nrow2] - W[nrow3]);
                
                a.y += C[nrow1] * (W[nrow2] + W[nrow3]) -
                S[nrow1] * (V[nrow2] + V[nrow3]);
                
                --nrow2;
                --nrow3;
                
                a.z += C[nrow1] * V[nrow2] + S[nrow1] * W[nrow2];
                
                --nrow1;
            }
            
            twice_cn0 = C[nrow1] + C[nrow1]; // Cheaper than multiplying by 2
            
            a.x += (twice_cn0 * V[nrow2]);
            
            a.y += (twice_cn0 * W[nrow2]);
            
            --nrow2;
            
            a.z += (C[nrow1] * V[nrow2] + S[nrow1] * W[nrow2]);
            
            --nrow1;
            nrow2 = nrow3;
            nrow3 -= 2;
            
        }
        
        // Finally add the acceleration for the monopole case:
        twice_cn0 = C[0] + C[0]; // Cheaper than multiplying by 2
        a.x += (twice_cn0 * V[2]);
        a.y += (twice_cn0 * W[2]);
        a.z += (C[0] * V[1] + S[0] * W[1]);
        
        a.z += a.z; // Cheaper than multiplying by 2
        
        double minusGMR2 = -m_GM/ ( m_R * m_R );
        
        a *= minusGMR2*1000.0;  // m/s^2 in ECEF
        
        GSpaceEnv::eop.ECEF2ECI_pos(a, a_eci);
        
        setForce(a_eci);
        
        if( m_hasPartialDerivatives == true)
        {
            getPartialDerivatives(pos_ecef);
        }
        
    }
    
    
}

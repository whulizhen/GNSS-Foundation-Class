
/*! @file Force_earth_gravity.cpp
	@author David Harrison
	@date 16 May 2016
	@brief SGNL OPS implementation of Earth gravity.
 */

#ifndef GFMGravity_hpp
#define GFMGravity_hpp

#include <stdio.h>

#include <fstream>

#include "GForceModel.hpp"

#include "GTidalCorrection.hpp"

namespace gfc
{
    
    
    
    /*gravity force model
     *
     * http://cddis.gsfc.nasa.gov/926/egm96/egm96.html
     * http://earth-info.nga.mil/GandG/wgs84/gravitymod/egm96/egm96.html
     * the gravity force is expressed in ECEF
     */
    class GFMGravity : public GForceModel
    {
        
    public:
        
        GFMGravity(int degree,int order)
        {
            setForceName("GFMGravity");
            
            //setMaxDegree(70);
            //setMaxOrder(70);
            m_N = degree;
            
            m_M = order;
            // For V' and W' function coefficients
            h = std::vector<double>((m_N + 2) * ( m_N + 3) / 2, 0.0);
            
            // For V' and W' functions , should be two higher in degree than C ans S
            V = std::vector<double>((m_N + 3) * (m_N + 4) / 2, 0.0);
            W = std::vector<double>((m_N + 3) * (m_N + 4) / 2, 0.0);
            
            // For denormalised C' and S' coefficients
            C = std::vector<double>((m_N + 1) * (m_N + 2) / 2, 0.0);
            S = std::vector<double>((m_N + 1) * (m_N + 2) / 2, 0.0);
            
            
        }
        
        ~GFMGravity() { }
        
        void loadEGM96( GString filename );
        void loadEGM2008(gfc::GString filename);
        
        void loadGRACE_ggm03c( GString filename );
        
        void correctTideCS(GTime ctUTC,double ut1mutc,bool solidFlag, bool poleFlag, bool oceanFlag);
        
        // must implement this function
        void doCompute(GTime ctUTC,double ut1mutc, GVector& pos_ecef  );
        
        // the general relativity correction
        void generalRelativityCorrection( double* pv_ecef, double* a );
        
        // also, we need the pole, solid earth and ocean tide corrections for the gravity
        
        // get the dadr,
        void getPartialDerivatives(GVector& pos_ecef);
        
    private:
        
        //generate the h coefficient
        void generate_h_coefs();
        
        void generate_VW_prime(double* ecef);
        
        std::vector<std::vector<double>> generate_denorm_factors( int max_n, int max_m) const;
        
        void denorm_coef( const std::vector<std::vector<double>> &C_norm, const std::vector<std::vector<double>> &S_norm);
        
        double m_GM;  // the value of GM constant for earth
        
        double m_R;  // the semi-major axis
        
        int m_M;  // order number
        
        int m_N;  // degree number
        
        //these variables are for solid earth tide corrections
        // ref IERS2010 techn Note36 Page80, equation 6.4
        double refMJD = 46431.0; // reference MJD of the force model
        double dotC20;
        double dotC21;
        double dotC30;
        double dotC40;
        
        double dotS21;
        
        // for the denormalised C' and S' coefficients
        std::vector<double> C;
        std::vector<double> S;
        
        //unnormalized coefficients
        // Cp and Sp are the copy of C and S, without any corrections
        std::vector<double> Cp;
        std::vector<double> Sp;
        
        std::vector<double> h;
        
        // for the V' and W' functions
        std::vector<double> V;
        std::vector<double> W;
        
        
        GString m_gravityCoefFilePath;  // the path for the harmonic spherical coefficients
        
    };
    
}


#endif /* GFMGravity_hpp */

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

/**
* @file IDATA.h
*/

#ifndef GFC_GOBSDATA_H
#define GFC_GOBSDATA_H

#include "GString.h"
#include "GFCCONST.h"
#include "GSensor.h"

namespace gfc
{
    
    NEW_GEXCEPTION_CLASS( frequencyUnexist, gfc::GException );
    
    NEW_GEXCEPTION_CLASS( GOBSTYPEUnexist,  gfc::GException );
    
    // class to manage all the radio frequencies for all the measurements, such as L band , C band , K band , X band  etc.
    // this class cannot be inherited
    class GCarrierFreq
    {
        
    public:
        
        /*
          静态函数，用于对frequencyMap进行初始化
        */
        
        static std::map<const GString, double > Initializer();
        static void RegByName(GString variableName, double variableValue);
        
        static void UnregByName( GString variableName );
        static void dump( std::ostream& s );
        static bool IsValid( GString freqname )
        {
            if( frequencyMap.count(freqname)>0 )
            {
                return true;
            }
            else
            {
                return false;
            }
        }
        
        static double GetByName( GString variableName );
        
        static double Getwavelen(GString variableName);
        
        static const GString* GetPointer(GString variableName);
        
        //GString   m_freqName;  // cf_L1 cf_L2 cf_L3 cf_L4 cf_L5 cf_L6 cf_L7 cf_L8 cf_for cf_GNSS or X1 C1 C2 etc.
        
        
    private:
        
        static std::map<const GString, double >  frequencyMap;  // storing all the frequency information
        
        
        //double m_freqValue;  // the value of frequency, unit(Hz)
        
    };
    
    
    
    // this class should be inherited by different datatype classes for different applications
    // this class should be used as  the father class
    class  GOBSType
    {
        
        public:
        
        GOBSType() {}
        
        // IOBSTYPE(GString datatypename) {m_datatypeName = datatypename;}
        virtual ~GOBSType() {}
        
//        GOBSType( GString obstypeName )
//        {
//            if( IsValid(obstypeName) == true)
//             {
//                 m_dataTypeName = obstypeName;
//             }
//             else
//             {
//                 GOBSTYPEUnexist ir("GOBSTYPE unexist!");
//                 ir.addLocation(FILE_LOCATION);
//                 GFC_THROW( ir );
//             }
//        }
        
        //
//        GOBSType& operator= (const GOBSType& right)  //赋值重载
//        {
//            this->m_dataTypeName = right.m_dataTypeName;
//            return *this;
//        }
//        // copy construction function
//        GOBSType( const GOBSType& right )   //拷贝构造函数
//        {
//            this->m_dataTypeName = right.m_dataTypeName;
//        }
        
        //
        static std::list< GString > Initializer();
        static bool IsValid(GString& variableName);
        
        static void RegByName( GString variableName );
        static void UnregByName( GString variableName );
        static void dump( std::ostream& s );
        
        static GOBSType GetByName(GString sensorSysName, GString obsDataTypeName );
        static const GString* GetPointer(GString variableName);
        
//        void  setOBSDataTypeName(GString dataTypeName)
//        {
//            if( IsValid(dataTypeName) == true )
//            {
//                m_dataTypeName = dataTypeName;
//            }
//            else
//            {
//                GOBSTYPEUnexist ir("GOBSTYPE unexist!");
//                ir.addLocation(FILE_LOCATION);
//                GFC_THROW( ir );
//            }
//        }
        
        //GString  getOBSDataTypeName() {return m_dataTypeName;}
        //GString getSensorSysName() {return m_sensorSystem.getSensorSystemName();}
        
        
        //GString m_dataTypeName;       // name of data type
        
        private:
        // static variable for storing the names of registed datatype.
        // for the inherit purpose, this data should be public
        static std::list< GString >  datatypelist;
        
        
        //GSensorSystem  m_sensorSystem;
        //const GString* m_freq;
        
    };  //end of class IDATATYPE
    
    
    // we need to solve this design problem that diffierent IDATA classes should have different IDATATYPE. It is difficult!
    // IDATA is a template class, different measurements data has different template DT.
    // this template class can also be inherited if it has to be
    //template< class DATA_TYPE >
    class GOBSData
    {
        public:
        
        GOBSData()
        {
            m_data_value = 0.0;   m_std_dev =999.9;
        }
        
        GOBSData( double dataValue, double stdDev )
        {
            //m_typestr = mydt;
            m_data_value = dataValue;
            m_std_dev  = stdDev;
        }
        
        virtual~GOBSData() {};
        
        void setValue( double value, double std_dev= 999.9 )
        {
            //m_typestr = dt;
            m_data_value = value;
            m_std_dev = std_dev;
        }
        
        void       setDataValue(double dataValue) {m_data_value = dataValue;}
        
        void       setDataDev(double dataDev)  {m_std_dev = dataDev;}
        
        //获取观测值类型
        //GString*   getDataType()  { return m_typestr;}
        double     getDataValue() {return m_data_value;}
        double     getDataDev()   {return m_std_dev;}
        
        
        //GString* m_typestr;
        
        double   m_data_value;  // the observable value
        double   m_std_dev;   //the stand deviation of this data
        
    }; // end of class IDATA
    
    
    
}  // end of namespace gfc




#endif








//#ifndef GFC_GDTYPE_HPP
//#define GFC_GDTYPE_HPP
//
//#include <iostream>
//#include <iomanip>
//#include <sstream>
//#include "GString.h"
//#include <map>
//#include "DataTypeBase.h"
//
////#include "RinexObsHeader.hpp"
////#include "Rinex3ObsHeader.hpp"
////#include "RinexObsID.hpp"
//
//namespace gfc
//{
//	
//	/** This class creates an index able to represent any type of observation,
//	*  correction, model parameter or other data value of interest for GNSS
//	*  data processing.
//	*
//	* This class is extensible in run-time, so the programmer may add
//	* indexes on-demand. For instance, in order to create a new GDType
//	* object referring INS-related data, and with "Inertial" as description
//	* string, you may write the following:
//	*
//	* @code
//	*    GDType INS = GDType::newValueType("Inertial");
//	* @endcode
//	*
//	* Or using the constructor:
//	*
//	* @code
//	*    GDType INS(GDType::newValueType("Inertial"));
//	* @endcode
//	*
//	* From now on, you'll be able to use INS as GDType when you need to
//	* refer to inertial system data.
//	*
//	*/
//	
//	class GDType : public DataTypeBase
//	{
//			
//	public:		
//		//”√ªß◊‘∂®“Âµƒπ€≤‚÷µ¿‡–Õµƒ∆µ„
//		const static int USER_OBSTYPE = 100;
//		
//		/*π€≤‚÷µ¿‡–Õ√∂æŸ*/
//		enum OBSTYPE
//		{
//			otUNKNOWN = 0,                //Œ¥÷™π€≤‚÷µ¿‡–Õ
//			otTIME,                   // ±º‰π€≤‚¡ø(µ•Œª:ƒ…√Î)
//			otRANGE,                //Œ±æ‡π€≤‚÷µ(µ•Œª:√◊)
//			otCARRIER,                //‘ÿ≤®π€≤‚÷µ(µ•Œª:÷‹)
//			otDOPPLER,                //∂‡∆’¿’π€≤‚÷µ(µ•Œª:∫’◊»)
//			otSTRENGTH,		        //–≈∫≈«ø∂»π€≤‚÷µ(µ•Œª:∑÷±¥db)
//			otANGLE,                //Ω«∂»π€≤‚÷µ(µ•Œª:ª°∂»)
//			otACCELERATION,           //º”ÀŸ∂»°¢¡¶π€≤‚÷µ°¢—π¡¶°¢∆¯—π°¢÷ÿ¡¶(µ•Œª:m/s2)
//			otVELOCITY,               //ÀŸ∂»π€≤‚÷µ°¢æ‡¿Î“ªΩ◊µº ˝(µ•Œª:m/s)
//			otTEMPERATURE,            //Œ¬∂»π€≤‚÷µ(µ•Œª:ø™∂˚Œƒ)
//			
//			otUserdefined1,      //”√ªß◊‘∂®“Â¿‡–Õ
//			otUserdefined2,      
//			otUserdefined3,      
//			otUserdefined4
//		};
//		
//		/*ŒﬁœﬂµÁ–≈∫≈≤®∂Œ∆µ¬ √∂æŸ£¨≤ªÕ¨Œ¿–«œµÕ≥æﬂ”–œ‡Õ¨µƒ≤®∂Œµ´ «∆µ¬ ≤ªÕ¨*/
//		enum OBSBAND{ FUNKNOWN,FANY, F1, F2, F3, F4, F5, F6, F7, F8 };
//				
//		/*π€≤‚÷µ Ù–‘√∂æŸ*/
//		enum OBSATTRIB{ oaUNKNOWN, oaC,oaS,oaL,oaX,oaP,oaW,oaY,oaM,oaD,oaN,oaI,oaQ,oaA,oaZ,oaUserdefined,oaLast };
//		
//		/// The type of the data value.
//		enum ValueType
//		{
//			Unknown,
//			// Observation-related types
//			C1,        ///< GPS civil code observation in L1 frequency
//			C2,        ///< GPS civil code observation in L2 frequency
//			P1,        ///< GPS precise code observation in L1 frequency
//			P2,        ///< GPS precise code observation in L2 frequency
//			L1,        ///< GPS phase observation in L1 frequency
//			L2,        ///< GPS phase observation in L2 frequency
//			D1,        ///< GPS doppler observation in L1 frequency
//			D2,        ///< GPS doppler observation in L2 frequency
//			S1,        ///< GPS signal strength observation in L1 frequency
//			S2,        ///< GPS signal strength observation in L2 frequency
//			T1,        ///< Transit integrated doppler observation in L1 frequency
//			T2,        ///< Transit integrated doppler observation in L2 frequency
//			SSI1,      ///< Signal strength indicator/index, L1 frequency
//			LLI1,      ///< Loss of Lock Indicator/ lock count, L1 frequency
//			SSI2,      ///< Signal strength indicator/index, L2 frequency
//			LLI2,      ///< Loss of Lock Indicator/ lock count, L2 frequency
//			// v 2.11
//			C5,        ///< GPS L5C-code pseudorange
//			L5,        ///< GPS phase observation in L5 frequency
//			D5,        ///< GPS doppler observation in L5 frequency
//			S5,        ///< GPS signal strength observation in L5 frequency
//			SSI5,      ///< Signal strength indicator/index, L5 frequency
//			LLI5,      ///< Loss of Lock Indicator/ lock count, L5 frequency
//			// Galileo-related
//			C6,        ///< Galileo E6-code pseudorange
//			L6,        ///< Galileo phase observation in L6 frequency
//			D6,        ///< Galileo doppler observation in L6 frequency
//			S6,        ///< Galileo signal strength observation in L6 frequency
//			SSI6,      ///< Signal strength indicator/index, L6 frequency
//			LLI6,      ///< Loss of Lock Indicator/ lock count, L6 frequency
//			C7,        ///< Galileo E5b-code pseudorange
//			L7,        ///< Galileo phase observation in L5b frequency
//			D7,        ///< Galileo doppler observation in L5b frequency
//			S7,        ///< Galileo signal strength observation in L5b frequency
//			SSI7,      ///< Signal strength indicator/index, L5b frequency
//			LLI7,      ///< Loss of Lock Indicator/ lock count, L5b frequency
//			C8,        ///< Galileo E5a+b-code pseudorange
//			L8,        ///< Galileo phase observation in L5a+b frequency
//			D8,        ///< Galileo doppler observation in L5a+b frequency
//			S8,        ///< Galileo signal strength observation in L5a+b frequency
//			SSI8,      ///< Signal strength indicator/index, L5a+b frequency
//			LLI8,      ///< Loss of Lock Indicator/ lock count, L5a+b frequency
//			// Combination-related types
//			PC,        ///< Code-based ionosphere-free combination
//			LC,        ///< Phase-based ionosphere-free combination
//			PI,        ///< Code-based ionospheric combination
//			LI,        ///< Phase-based ionospheric combination
//			Pdelta,    ///< Narrow-lane combination
//			Ldelta,    ///< Wide-lane combination
//			MWubbena,  ///< Melbourne-Wubbena combination
//			GRAPHIC1,  ///< GRoup And PHase Ionospheric Combination in L1
//			GRAPHIC2,  ///< GRoup And PHase Ionospheric Combination in L2
//			GRAPHIC5,  ///< GRoup And PHase Ionospheric Combination in L5
//			GRAPHIC6,  ///< GRoup And PHase Ionospheric Combination in L6
//			GRAPHIC7,  ///< GRoup And PHase Ionospheric Combination in L7
//			GRAPHIC8,  ///< GRoup And PHase Ionospheric Combination in L8
//			WL,        ///< Wide-lane combination(+1*L1-1*L2)
//			WL1,       ///< Wide-lane combination(-1*L1+2*L2)
//			WL2,       ///< Wide-lane combination(-2*L1+3*L2)
//			WL3,       ///< Wide-lane combination(-3*L1+4*L2)
//			WL4,       ///< Wide-lane combination(+4*L1-5*L2)
//			EWL,       ///< Wide-lane combination(-7*L1+9*L2)
//			// Derivatives of observations and combinations
//			L1dot,     ///< GPS L1 phase observation first derivative
//			L1dot2,    ///< GPS L1 phase observation second derivative
//			L2dot,     ///< GPS L2 phase observation first derivative
//			L2dot2,    ///< GPS L2 phase observation second derivative
//			L5dot,     ///< GPS L5 phase observation first derivative
//			L5dot2,    ///< GPS L5 phase observation second derivative
//			P1dot,     ///< GPS P1 precise code observation first derivative
//			P1dot2,    ///< GPS P1 precise code observation second derivative
//			P2dot,     ///< GPS P2 precise code observation first derivative
//			P2dot2,    ///< GPS P2 precise code observation second derivative
//			P5dot,     ///< GPS P5 precise code observation first derivative
//			P5dot2,    ///< GPS P5 precise code observation second derivative
//			L6dot,     ///< Galileo L6 phase observation first derivative
//			L6dot2,    ///< Galileo L6 phase observation second derivative
//			L7dot,     ///< Galileo L7 phase observation first derivative
//			L7dot2,    ///< Galileo L7 phase observation second derivative
//			L8dot,     ///< Galileo L8 phase observation first derivative
//			L8dot2,    ///< Galileo L8 phase observation second derivative
//			LCdot,     ///< Phase-based ionosphere-free combination 1st derivative
//			LCdot2,    ///< Phase-based ionosphere-free combination 2nd derivative
//			LIdot,     ///< Phase-based ionospheric combination 1st derivative
//			LIdot2,    ///< Phase-based ionospheric combination 2nd derivative
//			Ldeltadot, ///< Wide-lane combination 1st derivative
//			Ldeltadot2,///< Wide-lane combination 2nd derivative
//			// Model-related types
//			transmit,  ///< Transmit time of the signal
//			rho,       ///< Geometric distance satellite-receiver
//			rhodot,    ///< First derivative of geometric distance SV-RX
//			rhodot2,   ///< Second derivative of geometric distance SV-RX
//			dtSat,     ///< Satellite clock offset
//			dtSatdot,  ///< Satellite clock offset drift
//			dtSatdot2, ///< Satellite clock offset drift rate
//			rel,       ///< Relativistic delay
//			gravDelay, ///< Gravitational delay
//			tropo,     ///< Vertical tropospheric delay, total
//			dryTropo,  ///< Vertical tropospheric delay, dry component
//			dryMap,    ///< Tropospheric mapping function, dry component
//			wetTropo,  ///< Vertical tropospheric delay, wet component
//			wetMap,    ///< Tropospheric mapping function, wet component
//			tropoSlant, ///< Slant tropospheric delay, total
//			iono,      ///< Vertical ionospheric delay
//			ionoTEC,   ///< Total Electron Content (in TECU), 1TECU = 1e+16 electrons per m**2
//			ionoMap,   ///< Ionospheric mapping function
//			ionoMap2,  ///< Ionospheric mapping function for second order ionospheric delay
//			ionoL1,    ///< Slant ionospheric delay, frequency L1
//			ionoL2,    ///< Slant ionospheric delay, frequency L2
//			ionoL5,    ///< Slant ionospheric delay, frequency L5
//			ionoL6,    ///< Slant ionospheric delay, frequency L6
//			ionoL7,    ///< Slant ionospheric delay, frequency L7
//			ionoL8,    ///< Slant ionospheric delay, frequency L8
//			windUp,    ///< Wind-up effect (in radians)
//			satPCenter,///< Satellite antenna phase center correction
//			satX,      ///< Satellite position, X component
//			satY,      ///< Satellite position, Y component
//			satZ,      ///< Satellite position, Z component
//			satVX,     ///< Satellite velocity, X component
//			satVY,     ///< Satellite velocity, Y component
//			satVZ,     ///< Satellite velocity, Z component
//			satAX,     ///< Satellite acceleration, X component
//			satAY,     ///< Satellite acceleration, Y component
//			satAZ,     ///< Satellite acceleration, Z component
//			satJ2kX,   ///< Satellite position in J2000, X component
//			satJ2kY,   ///< Satellite position in J2000, Y component
//			satJ2kZ,   ///< Satellite position in J2000, Z component
//			satJ2kVX,  ///< Satellite velocity in J2000, X component
//			satJ2kVY,  ///< Satellite velocity in J2000, Y component
//			satJ2kVZ,  ///< Satellite velocity in J2000, Z component
//			satJ2kAX,  ///< Satellite acceleration in J2000, X component
//			satJ2kAY,  ///< Satellite acceleration in J2000, Y component
//			satJ2kAZ,  ///< Satellite acceleration in J2000, Z component
//			elevation, ///< Satellite elevation
//			azimuth,   ///< Satellite azimuth
//			// Cycle slip flags
//			CSL1,      ///< Cycle slip in L1
//			CSL2,      ///< Cycle slip in L2
//			CSL5,      ///< Cycle slip in L5
//			CSL6,      ///< Cycle slip in L6
//			CSL7,      ///< Cycle slip in L7
//			CSL8,      ///< Cycle slip in L8
//			// Satellite 'arcs'
//			satArc,    ///< Satellite arc number
//			// Phase-ambiguity types
//			BL1,       ///< Phase ambiguity in L1
//			BL2,       ///< Phase ambiguity in L2
//			BL5,       ///< Phase ambiguity in L5
//			BL6,       ///< Phase ambiguity in L6
//			BL7,       ///< Phase ambiguity in L7
//			BL8,       ///< Phase ambiguity in L8
//			BLC,       ///< Phase ambiguity in LC
//			BWL,       ///< Phase ambiguity in WL
//			BWL2,      ///< Phase ambiguity in WL2
//			BWL3,      ///< Phase ambiguity in WL3
//			BWL4,      ///< Phase ambiguity in WL4
//			// Multipath-related types
//			mpC1,      ///< Multipath bias, C1
//			mpC2,      ///< Multipath bias, C2
//			mpC5,      ///< Multipath bias, C5
//			mpC6,      ///< Multipath bias, C6
//			mpC7,      ///< Multipath bias, C7
//			mpC8,      ///< Multipath bias, C8
//			mpL1,      ///< Multipath bias, L1
//			mpL2,      ///< Multipath bias, L2
//			mpL5,      ///< Multipath bias, L5
//			mpL6,      ///< Multipath bias, L6
//			mpL7,      ///< Multipath bias, L7
//			mpL8,      ///< Multipath bias, L8
//			// Instrumental delays types
//			instC1,    ///< Instrumental delay, C1
//			instC2,    ///< Instrumental delay, C2
//			instC5,    ///< Instrumental delay, C5
//			instC6,    ///< Instrumental delay, C6
//			instC7,    ///< Instrumental delay, C7
//			instC8,    ///< Instrumental delay, C8
//			instL1,    ///< Instrumental delay, L1
//			instL2,    ///< Instrumental delay, L2
//			instL5,    ///< Instrumental delay, L5
//			instL6,    ///< Instrumental delay, L6
//			instL7,    ///< Instrumental delay, L7
//			instL8,    ///< Instrumental delay, L8
//			// Equation system-related types
//			prefitP1,  ///< Prefit residual, code P1
//			prefitP2,  ///< Prefit residual, code P2
//			prefitL1,  ///< Prefit residual, phase L1
//			prefitL2,  ///< Prefit residual, phase L2
//			postfitP1, ///< Postfit residual, code P1
//			postfitP2, ///< Postfit residual, code P2
//			postfitL1, ///< Postfit residual, phase L1
//			postfitL2, ///< Postfit residual, phase L2
//			prefitC5,  ///< Prefit residual, code C5
//			prefitL5,  ///< Prefit residual, phase L5
//			postfitC5, ///< Postfit residual, code C5
//			postfitL5, ///< Postfit residual, phase L5
//			prefitGRAPHIC1,   ///< Prefit residual, GRAPHIC1
//			prefitGRAPHIC2,   ///< Prefit residual, GRAPHIC2
//			postfitGRAPHIC1,  ///< Postfit residual, GRAPHIC1
//			postfitGRAPHIC2,  ///< Postfit residual, GRAPHIC2
//			prefitMWubbena,   /// Prefit residual, MWubbena
//			prefitWL,  ///< Prefit residual, WL
//			prefitWL2, ///< Prefit residual, WL2
//			prefitWL3, ///< Prefit residual, WL3
//			prefitWL4, ///< Prefit residual, WL4
//			postfitWL, ///< Postfit residual, WL
//			postfitWL2,///< Postfit residual, WL2
//			postfitWL3,///< Postfit residual, WL3
//			postfitWL4,///< Postfit residual, WL4
//			prefitC,   ///< Prefit residual, code
//			prefitL,   ///< Prefit residual, phase
//			postfitC,  ///< Postfit residual, code
//			postfitL,  ///< Postfit residual, phase
//			dx,        ///< In the position domain: Position bias, X component; in the range domain: dx coefficient
//			dy,        ///< In the position domain: Position bias, Y component; in the range domain: dy coefficient
//			dz,        ///< In the position domain: Position bias, Z component; in the range domain: dz coefficient
//			dLat,      ///< Position bias, Latitude component
//			dLon,      ///< Position bias, Longitude component
//			dH,        ///< Position bias, Height component
//			cdt,       ///< In the position domain: Receiver clock offset, meters; in the range domain: cdt coefficient
//			cdtSat,    ///< In the position domain: Satellite clock offset, meters; in the range domain: cdt coefficient
//			dSatX,     ///< dSatX coefficient for satellite position in XYZ
//			dSatY,     ///< dSatY coefficient for satellite position in XYZ
//			dSatZ,     ///< dSatZ coefficient for satellite position in XYZ
//			dSatR,     ///< dSatR coefficient for satellite position in RTN
//			dSatT,     ///< dSatT coefficient for satellite position in RTN
//			dSatN,     ///< dSatN coefficient for satellite position in RTN
//			weight,    ///< Weight assigned to a given observation
//			codeBias,  ///< Code bias by both receiver and satellite
//			cdtC1,     ///< Receiver clock offset of C1
//			cdtP1,     ///< Receiver clock offset of P1
//			cdtC2,     ///< Receiver clock offset of C2
//			cdtP2,     ///< Receiver clock offset of P2
//			cdtC5,     ///< Receiver clock offset of C5
//			cdtP5,     ///< Receiver clock offset of P5
//			cdtL1,     ///< Receiver clock offset of L1
//			cdtL2,     ///< Receiver clock offset of L2
//			cdtL5,     ///< Receiver clock offset of L5
//			cdtPC,     ///< Receiver clock offset of PC
//			cdtLC,     ///< Receiver clock offset of LC
//			cdtWL,     ///< Receiver clock offset of WL
//			cdtWL2,    ///< Receiver clock offset of WL2
//			cdtWL3,    ///< Receiver clock offset of WL3
//			cdtWL4,    ///< Receiver clock offset of WL4
//			cdtMW,     ///< Receiver clock offset of MW
//			cdtSatC1,  ///< Satellite clock offset of C1
//			cdtSatP1,  ///< Satellite clock offset of P1
//			cdtSatC2,  ///< Satellite clock offset of C2
//			cdtSatP2,  ///< Satellite clock offset of P2
//			cdtSatC5,  ///< Satellite clock offset of C5
//			cdtSatP5,  ///< Satellite clock offset of P5
//			cdtSatL1,  ///< Satellite clock offset of L1
//			cdtSatL2,  ///< Satellite clock offset of L2
//			cdtSatL5,  ///< Satellite clock offset of L5
//			cdtSatPC,  ///< Satellite clock offset of PC
//			cdtSatLC,  ///< Satellite clock offset of LC
//			cdtSatWL,  ///< Satellite clock offset of WL
//			cdtSatMW,  ///< Satellite clock offset of MW
//			// Other types
//			recX,      ///< Receiver position, X component
//			recY,      ///< Receiver position, Y component
//			recZ,      ///< Receiver position, Z component
//			recVX,     ///< Receiver velocity, X component
//			recVY,     ///< Receiver velocity, Y component
//			recVZ,     ///< Receiver velocity, Z component
//			recAX,     ///< Receiver acceleration, X component
//			recAY,     ///< Receiver acceleration, Y component
//			recAZ,     ///< Receiver acceleration, Z component
//			recLat,    ///< Receiver position, Latitude component
//			recLon,    ///< Receiver position, Longitude component
//			recH,      ///< Receiver position, Height component
//			recVLat,   ///< Receiver velocity, Latitude component
//			recVLon,   ///< Receiver velocity, Longitude component
//			recVH,     ///< Receiver velocity, Height component
//			recALat,   ///< Receiver acceleration, Latitude component
//			recALon,   ///< Receiver acceleration, Longitude component
//			recAH,     ///< Receiver acceleration, Height component
//			recJ2kX,   ///< Receiver position in J2000, X component
//			recJ2kY,   ///< Receiver position in J2000, Y component
//			recJ2kZ,   ///< Receiver position in J2000, Z component
//			recJ2kVX,  ///< Receiver velocity in J2000, X component
//			recJ2kVY,  ///< Receiver velocity in J2000, Y component
//			recJ2kVZ,  ///< Receiver velocity in J2000, Z component
//			recJ2kAX,  ///< Receiver acceleration in J2000, X component
//			recJ2kAY,  ///< Receiver acceleration in J2000, Y component
//			recJ2kAZ,  ///< Receiver acceleration in J2000, Z component
//			sigma,     ///< Standard deviation
//			iura,      ///< Index User Range Accuracy
//			Action,    ///< Flag for quality control
//			// Handy dummy types for non-standard processing
//			dummy0,    ///< Generic, undefined type #0
//			dummy1,    ///< Generic, undefined type #1
//			dummy2,    ///< Generic, undefined type #2
//			dummy3,    ///< Generic, undefined type #3
//			dummy4,    ///< Generic, undefined type #4
//			dummy5,    ///< Generic, undefined type #5
//			dummy6,    ///< Generic, undefined type #6
//			dummy7,    ///< Generic, undefined type #7
//			dummy8,    ///< Generic, undefined type #8
//			dummy9,    ///< Generic, undefined type #9
//
//			Last,      ///< used to extend this...
//			Placeholder = Last+1000
//		};
//
//
//		/// empty constructor, creates an invalid object
//		GDType()
//			: type(Unknown) {};
//
//
//		/** Explicit constructor
//		*
//		* @param vt   ValueType for the new GDType. If you want to use the
//		*             next available ValueType, generate it using the
//		*             'newValueType()' method, as indicated in the example in
//		*             the documentation.
//		*/
//
//		GDType(ValueType vt)
//			: type(vt) {};
//
//
//		/// Equality requires all fields to be the same
//		virtual bool operator==(const GDType& right) const
//		{ return type==right.type; };
//
//
//		/// This ordering is somewhat arbitrary but is required to be able
//		/// to use an GDType as an index to a std::map. If an application
//		/// needs some other ordering, inherit and override this function.
//		virtual bool operator<(const GDType& right) const
//		{ return type < right.type; };
//
//
//		/// Inequality operator
//		bool operator!=(const GDType& right) const
//		{ return !(operator==(right)); };
//
//
//		/// Greater than operator
//		bool operator>(const GDType& right) const
//		{  return (!operator<(right) && !operator==(right)); };
//
//
//		/// Less than or equal operator
//		bool operator<=(const GDType& right) const
//		{ return (operator<(right) || operator==(right)); };
//
//
//		/// Greater than or equal operator
//		bool operator>=(const GDType& right) const
//		{ return !(operator<(right)); };
//		
//		
//		/// Assignment operator
//		virtual GDType operator=(const GDType& right);
//		
//		
//		/// Convenience output method
//		virtual void dump(std::ostream& s) const;
//		
//		
//		/// Returns true if this is a valid GDType. Basically just
//		/// checks that the enum is defined
//		virtual bool isValid() const;
//
//
//		/// Destructor
//		virtual ~GDType() {};
//
//
//		/** Static method to add new GDType's
//		* @param s      Identifying string for the new GDType
//		*/
//		static ValueType newValueType(const GString& s);
//		
//			
//		/// Type of the value
//		ValueType type;
//		
//		
//		/// Map holding type descriptions
//		static std::map< ValueType, GString > tStrings;
//		
//		
//	public:
//		class Initializer
//		{
//		public:
//			Initializer();
//		};
//
//		static Initializer TypeIDsingleton;
//
//	public:
//
//		virtual GString getClassName(void) const;
//		/** Static method to get the user registered GDType by name string
//		* @param name      Identifying string for the new GDType
//		* @return          The desired GDType
//		*/
//		static GDType byName(GString name)
//			throw(InvalidRequest);
//        
//		/** Static method to add new GDType's by name string
//		* @param name      Identifying string for the new GDType
//		* @param desc      Descriptions of the new GDType
//		* @return          The new GDType
//		*/
//		static GDType regByName( GString name,GString desc );
//		
//		/// unregister a GDType by it's name string
//		static void unregByName( GString name );
//		
//		/// unregister all TypeIDs registered by name string
//		static void unregAll();
//		
//	private:
//		
//		/// Have user deined TypeIDs been registered ?
//		static bool bUserTypeIDRegistered;
//		/// Map holding user defined TypeIDs by a string
//		static std::map<GString,GDType> mapUserTypeID;
//														
//														
//		//lizhen defined
//		unsigned char m_valueType;   // ˝æ›À˘¥˙±ÌµƒŒÔ¿Ì∫¨“Â
//		unsigned char m_frequency;   //ƒø«∞Ωˆ÷ß≥÷ŒﬁœﬂµÁ–≈∫≈£¨…˘≤®–≈∫≈º∞∆‰À˚–≈∫≈£®¥˝¿©’π£©	
//		unsigned char m_attribution;  // ˝æ›µƒ Ù–‘		
//		
//		
//	}; // End of class 'GDType'
//
//
//
//	namespace StringUtils
//	{
//		/// convert this object to a string representation
//		std::string asString(const GDType& p);
//	}
//
//
//
//	/// stream output for GDType
//	std::ostream& operator<<(std::ostream& s, const GDType& p);
//
//
//	// bool IsCarrierPhase(const RinexObsType& rot);
//
//	// int GetCarrierBand(const RinexObsType& rot);
//
//	// int GetCarrierBand(const RinexObsID& roi);
//
//	/* GDType::ValueType ConvertToTypeID(const RinexObsType& rot,
//	const RinexSatID& sat);
//	*/
//	/* GDType::ValueType ConvertToTypeID(const RinexObsID& roi,
//	const RinexSatID& sat);
//	*/
//
//}  // End of namespace gpstk
//
//#endif   // GPSTK_TYPEID_HPP


/*
@file Sensor.h
 class for all the Sensor Systems
*/

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

#ifndef GFC_SENSOR_H
#define GFC_SENSOR_H



#include "GString.h"
#include <algorithm>

namespace gfc
{
	/** @defgroup Sensor System base class file*/
	//@{
	/** This level isn't serving any purpose other than to make
	* the class diagram look nice...
	*/
    
    NEW_GEXCEPTION_CLASS(SensorSystemUnexist, gfc::GException );
    
	class GSensorSystem
	{
		//sensor system  还需要添加其他的一些信息，不只是systemName
        
	public:
        
        GSensorSystem(void) { m_sensorName="ssUKN";}
        GSensorSystem(GString sensorName)
        {
            if( IsValid( sensorName) == true)
            {
                m_sensorName = sensorName;
            }
            else
            {
                SensorSystemUnexist ir("SensorSystem unexist!");
                ir.addLocation(FILE_LOCATION);
                GFC_THROW( ir );
            }
        }
        virtual ~GSensorSystem(void) {};
        
        //
        GSensorSystem& operator= (const GSensorSystem& right)  //赋值重载
        {
            this->m_sensorName = right.m_sensorName;
            return *this;
        }
        // copy construction function
        GSensorSystem( const GSensorSystem& right )   //拷贝构造函数
        {
            this->m_sensorName = right.m_sensorName;
        }
        
        static std::list<GString > Initializer();
        static void RegByName( GString variableName );
        
        static void UnregByName( GString variableName );
        static void dump( std::ostream& s ) ;
        static bool IsValid(GString sensorSysName)
        {
            std::list<GString>::iterator myIterator;
            myIterator = find(sensorSystemList.begin(),sensorSystemList.end(),sensorSysName);
            if( myIterator == sensorSystemList.end() )
            {
                return false;
            }
            else
            {
                return true;
            }
        }
        
        static GSensorSystem GetByName( GString variableName );
        
        // return the name(GString) of certain timesystem, this function need to be rewritten
        GString getSensorSystemName() const
        {
            return m_sensorName;
        }
        
        void     setSensorSystemName(GString sensorName)
        {
           if( IsValid( sensorName) == true)
           {
               m_sensorName = sensorName;
           }
           else
           {
                SensorSystemUnexist ir("SensorSystem unexist!");
                ir.addLocation(FILE_LOCATION);
                GFC_THROW( ir );
           }
        }
        
    private:
        
        // list for all the sensor systems
        static std::list<GString> sensorSystemList;
        
        GString  m_sensorName;
        
    };
	//@}
    
    // a satellite is a kind of sensor, So does a INS and even a totalstation ......
    // class SensorID is a node of the observable TREE ( a tree with multiple types of node and multiple layers).
    class GSensorID
    {
        
    public:
        GSensorID() { m_id = -1;}
        GSensorID(GString sensorName, int sensorid) { m_ss = GSensorSystem::GetByName(sensorName); m_id = sensorid; }
        //deconstruction function must be virtual for the purpose of be derived
        virtual ~GSensorID() {}
        
        void setData(GString sensorName, int sensorid) {m_ss = GSensorSystem::GetByName(sensorName); m_id = sensorid; }
        
        GSensorID& operator= (const GSensorID& right)  //赋值重载
        {
            this->m_ss = right.m_ss;
            this->m_id = right.m_id;
            
            return *this;
        }
        
        
        GSensorID( const GSensorID& right )   //copy construction function
        {
            this->m_ss = right.m_ss;
            this->m_id = right.m_id;
        }
        
        GString getIDString() const
        {
            GString idString = m_ss.getSensorSystemName() + GString(m_id);
            
            return idString;
        }
        
        int getIDnum() {return m_id;}
        
        GString getSystem() { return m_ss.getSensorSystemName(); }
        
       static GSensorID ilrsName2GSensorID(GString spacecraftName);
        
        // need to be rewritten by differernt classes
        // because different classes have different format for the sensor ID
        virtual void setID(GString idString) {}
        
        bool Available();
        
        
    private:
        GSensorSystem  m_ss;
        int         m_id; // if m_ss is a kind of GNSS , then m_id is the prn number of the certain satellite
        
    };
    
    
    bool   operator== (const GSensorID& left,const GSensorID& right) ;
    bool   operator!= (const GSensorID& left,const GSensorID& right) ;
    bool   operator>  (const GSensorID& left,const GSensorID& right) ;
    bool   operator<  (const GSensorID& left,const GSensorID& right) ;
    bool   operator>= (const GSensorID& left,const GSensorID& right) ;
    bool   operator<= (const GSensorID& left,const GSensorID& right) ;
    
}  // end of namespace

#endif



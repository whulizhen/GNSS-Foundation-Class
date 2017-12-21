
/*
@file Sensor.cpp
class for all the  Sensor System,ie. all the measurement sensor Systems
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

#include "GSensor.h"

namespace gfc
{
	/** @defgroup GNSS Sensor System class file*/
	//@{
    
    //类TimeSystem相关的初始化
    std::list<GString> gfc::GSensorSystem::sensorSystemList  = gfc::GSensorSystem::Initializer();
    std::list<GString> GSensorSystem::Initializer()
    {
        std::list<GString>  tmpTab;
        tmpTab.push_back("ssUKN"); // 未知传感器系统
        tmpTab.push_back("ssTEST"); // the test space object
        tmpTab.push_back("ssGPS");
        tmpTab.push_back("ssBDS");
        tmpTab.push_back("ssGLO");
        tmpTab.push_back("ssGAL");
        tmpTab.push_back("ssIRN"); // IRNSS
        tmpTab.push_back("ssQZS"); // QZSS
        tmpTab.push_back("ssSBS"); // SBAS payload
        tmpTab.push_back("ssINS");  //惯导系统
        
        tmpTab.push_back("ssTPS"); // 地面测量系统，如全站仪
        tmpTab.push_back("ssLMS");// 地面水准测量或者重力测量系统
        
        tmpTab.push_back("ssLAGEOS"); // LAGEOS satellte
        tmpTab.push_back("ssJASON");  //jason mission
        
        tmpTab.push_back("ssGRACE") ; //grace mission
        
        
        return tmpTab;
    }
    
    // register user defined sensor system
    void GSensorSystem::RegByName( GString variableName )
    {
        // if variableName does not exist, then register , otherwise throw out exception
        int NumberOfSS(0) ;
        NumberOfSS = static_cast<int>(count(sensorSystemList.begin(), sensorSystemList.end(), variableName ));
        
        if( NumberOfSS == 0 )
        {
            sensorSystemList.push_back(variableName);
        }
        else
        {
            printf("WARNING: SensorSystem : sensor system %s has already existed!\n",variableName.c_str());
        }
    }
    
    void GSensorSystem::UnregByName( GString variableName )
    {
        std::list<GString>::iterator myIterator;
        myIterator = find(sensorSystemList.begin(),sensorSystemList.end(),variableName);
        if( myIterator == sensorSystemList.end() )
        {
            printf("WARNING: SensorSystem : sensor system %s hasn't  existed!\n",variableName.c_str());
        }
        else  // if exist , then erase it
        {
            sensorSystemList.erase(myIterator);
        }
    }
    
    
    
    GSensorSystem GSensorSystem::GetByName( GString variableName )
    {
        std::list<GString>::iterator myIterator;
        myIterator = find(sensorSystemList.begin(),sensorSystemList.end(),variableName);
        if( myIterator == sensorSystemList.end() )
        {
            // throw out exception, timesystem does not exist
            //非法请求
            SensorSystemUnexist ir("SensorSystem unexist!");
            ir.addLocation(FILE_LOCATION);
            GFC_THROW( ir );
            
        }
        else
        {
            GSensorSystem  myts(*myIterator);
            return myts;
        }
    }
    
    
    void GSensorSystem::dump( std::ostream& s )
    {
        s<<"ClassName:    "<<"GSensorSystem"<<std::endl;
        
        //遍历所有的变量
        std::list<GString>::iterator myit = sensorSystemList.begin();
        for(;myit!= sensorSystemList.end(); ++myit)  // here, it must be ++myit, because this is much faster than myit++
        {
            s<<"SensorSystemName:  "<<*myit<<std::endl;
        }
    }
    
	//@}
    
    
    bool GSensorID::Available()
    {
        if( m_ss.getSensorSystemName() == "ssUKN" || m_id == -1)
        {
            return false;
        }
        else
        {
            return true;
        }
    }
    
    //static function for converting ILRS spacecarft name to GSensorID
    GSensorID GSensorID::ilrsName2GSensorID(GString spacecraftName)
    {
        spacecraftName.stripLeading_v();
        spacecraftName.stripTrailing_v();
        
        GSensorID myid;
        
        if( spacecraftName == "compassi3")
        {
            myid.setData("ssBDS", 8);
        }
        if( spacecraftName == "compassg1")
        {
            myid.setData("ssBDS", 1);
        }
        else if(spacecraftName == "compassi5")
        {
            myid.setData("ssBDS", 10);
        }
        else if(spacecraftName == "compassis1")
        {
            
        }
        else if(spacecraftName == "compassis2")
        {
            
        }
        else if(spacecraftName == "compassm3")
        {
             myid.setData("ssBDS", 11);
        }
        else if(spacecraftName == "compassms1")
        {
            
        }
        else if(spacecraftName =="compassms2")
        {
            
        }
        else if(spacecraftName == "galileo101")  //E11
        {
            myid.setData("ssGAL", 11);
        }
        else if(spacecraftName=="galileo102")
        {
            myid.setData("ssGAL", 12);
        }
        else if(spacecraftName=="galileo103")
        {
            myid.setData("ssGAL", 19);
        }
        else if(spacecraftName=="galileo104")
        {
            myid.setData("ssGAL", 20);
        }
        else if(spacecraftName=="galileo203")
        {
            myid.setData("ssGAL", 26);
        }else if(spacecraftName=="galileo204")
        {
            myid.setData("ssGAL", 22);
        }else if(spacecraftName=="galileo205")
        {
            myid.setData("ssGAL", 24);
        }else if(spacecraftName=="galileo206")
        {
            myid.setData("ssGAL", 30);
        }else if(spacecraftName=="galileo208")
        {
            myid.setData("ssGAL", 8);
        }
        else if(spacecraftName=="galileo209")
        {
            myid.setData("ssGAL", 9);
        }
        else if(spacecraftName=="galileo210")
        {
            myid.setData("ssGAL", 1);
        }
        else if(spacecraftName=="galileo211")
        {
            myid.setData("ssGAL", 2);
        }
        else if(spacecraftName == "lageos2" || spacecraftName == "Lageos2")
        {
            myid.setData("ssLAGEOS", 52);
        }
        else if(spacecraftName == "lageos1" || spacecraftName == "Lageos1")
        {
            myid.setData("ssLAGEOS", 51);
        }
        
        
        return myid;
    }
    
    
    
    bool  operator>=  (const GSensorID& left, const GSensorID& right)
    {
        return ( left.getIDString() >= right.getIDString() );
    }
    
    bool  operator<=  (const GSensorID& left,const GSensorID& right)
    {
        return ( left.getIDString() <= right.getIDString() );
    }
    
    bool  operator >  (const GSensorID& left,const GSensorID& right)
    {
        return ( left.getIDString() > right.getIDString() );
    }
    
    bool  operator <  (const GSensorID& left,const GSensorID& right)
    {
        return ( left.getIDString() < right.getIDString() );
    }
    
    bool   operator== (const GSensorID& left,const GSensorID& right)
    {
       return ( left.getIDString() == right.getIDString() );
    }
    
    /*不等号重载*/
    bool operator!=(const GSensorID& left,const  gfc::GSensorID &right)
    {
        return ( left.getIDString() != right.getIDString() );
    }

    
    
    
    
}  //end namespace

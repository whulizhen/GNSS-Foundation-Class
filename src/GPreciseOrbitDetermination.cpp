//
//  GPOD.cpp
//  GFC
//
//  Created by lizhen on 01/04/2017.
//  Copyright © 2017 lizhen. All rights reserved.
//

#include "GPreciseOrbitDetermination.hpp"
namespace gfc
{
    void GPreciseOrbitDetermination::loadObsdata(gfc::GString obsfilename)
    {
        allobsfile = new GRinex[200];
        
        int nfile = 0;
        GRinex myrinex;
        std::fstream  cfgfile;
        char ss[1024]={};
        cfgfile.open(obsfilename);
        
        if( cfgfile.is_open() == false )
        {
            printf("rinex file is unavailable\n");
            return;
        }
        while (!cfgfile.eof())
        {
            cfgfile.getline(ss, 1024);
            GString s(ss);
            if(s=="")
            {
                break;
            }
            
            GString stationName = s.substr(s.size()-12,4);
            GSourceID myid(stationName);
            allobsfile[nfile++].setStation(s, myid);
        }
        
        cfgfile.close();
        
        GgnssStorage mydatastorage;
        for(int i = 0 ; i< nfile; i++ )
        {
            allobsfile[i].nextOBS(mydatastorage);
        }
        
        
        //遍历所有的测站
        GgnssStorage::iterator  iter;
        for(  iter=mydatastorage.begin(); iter!=mydatastorage.end();   iter++)
        {
            //CString a= iter - > first;
            GString stationName = iter->first.getSourceName();
            int nGalileo=0;
            for(int i = 0 ; i< iter->second.m_indicatorlist.size(); i++)
            {
                if(iter->second.m_indicatorlist[i].getSystem()=="ssGAL")
                {
                    nGalileo++;
                }
            }
            
            printf("%s  %d \n",stationName.c_str(),nGalileo);
            
            
            
        }
        
        
        int testc = 0;
        
        
        
    }
    
    
    
    
    
} // end of namespace


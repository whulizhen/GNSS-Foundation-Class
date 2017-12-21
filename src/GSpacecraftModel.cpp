//
//  GSpacecraftModel.cpp
//  GFC
//
//  Created by lizhen on 16/06/2016.
//  Copyright © 2016 lizhen. All rights reserved.
//

#include "GSpacecraftModel.hpp"
namespace gfc
{
    
    std::map<GString, GSpaceCraftModel> GSpacecraftModelMgr::spacecraftModelStore = GSpacecraftModelMgr::initialiseModel();
    
    
    //quit important, may need to change it to configure file
    GString GSpacecraftModelMgr::sensorID2svType(GSensorID svID)
    {
        GString svType;
        
        if( svID.getSystem() == "ssGPS" )
        {
            svType = "GPSIIR";
        }
        else if( svID.getSystem() == "ssBDS")
        {
            if( svID.getIDnum()<=5)
            {
                svType = "BDSGEO";
            }
            else if(svID.getIDnum() > 5 && svID.getIDnum() <=10 )
            {
                svType = "BDSIGSO";
            }
            else
            {
                svType = "BDSMEO";
            }
        }
        else if(svID.getSystem() == "ssGAL")
        {
            //ref: https://www.gsc-europa.eu/system-status/Constellation-Information
            int prn = svID.getIDnum() ;
            if( prn == 11 || prn == 12 || prn == 19 || prn == 20 )
            {
                svType = "GALIOV";
            }
            else
            {
                svType = "GALFOC";
            }
            
        }
        else if(svID.getSystem() == "ssGLS")
        {
            svType = "GLSM"; //GLONASS M satellites
        }
        
        
        return svType;
    }
    
    
    std::map<GString, GSpaceCraftModel> GSpacecraftModelMgr::initialiseModel()
    {
        std::map<GString, GSpaceCraftModel> modelstorage;
        
        //GPS IIR
        GSpaceCraftModel mymodel;
        GSpaceCraftModel::opticalElement element;
        
        mymodel.m_mass = 1086.4512;  // kg !< SVN 46 for March 2004
        mymodel.m_antennaPower = 85.0; // 85.0 watta
        mymodel.m_solarPanelPowerDraw = 0.0;
        
        mymodel.m_bias.set(0.0, 0.0, 0.0);
        
        //front side of solar array
        element.name = "solarArray+";
        element.area = 13.564096;
        element.reflectivity = 0.28;
        element.specularity = 0.85;
        element.thickness = 0.0134876;
        element.absorbtivity = 0.78;
        element.conductivity = 1.86148;
        element.emmisivity = 0.82;
        mymodel.solarArray.push_back(element);
        
        // backside of solar array
        element.name = "solarArray-";
        element.area = 13.564096;
        element.reflectivity = 0.28;
        element.specularity = 0.85;
        element.thickness = 0.0134876;
        element.absorbtivity = 0.78;
        element.conductivity = 1.86148;
        element.emmisivity = 0.82;
        mymodel.solarArray.push_back(element);
        
        //busX+
        element.name = "busX+";
        element.normal.x = 1;element.normal.y = 0; element.normal.z = 0.0;
        element.area = 2.286568;
        element.reflectivity = 0.06;
        element.specularity = 0.0;
        element.thickness = 0.0;
        element.absorbtivity = 0.0;
        element.conductivity = 0.0;
        element.emmisivity = 0.0;
        mymodel.busX.push_back(element);
        
        //busX-
        element.name = "busX-";
        element.normal.x = -1;element.normal.y = 0; element.normal.z = 0.0;
        element.area = 2.286568;
        element.reflectivity = 0.06;
        element.specularity = 0.0;
        element.thickness = 0.0;
        element.absorbtivity = 0.0;
        element.conductivity = 0.0;
        element.emmisivity = 0.0;
        mymodel.busX.push_back(element);
        
        //busY+
        element.name = "busY+";
        element.normal.x = 0;element.normal.y = 1; element.normal.z = 0.0;
        element.area = 2.859672;
        element.reflectivity = 0.06;
        element.specularity = 0.0;
        element.thickness = 0.0;
        element.absorbtivity = 0.0;
        element.conductivity = 0.0;
        element.emmisivity = 0.0;
        mymodel.busY.push_back(element);
        
        //busY-
        element.name = "busY-";
        element.normal.x = 0;element.normal.y = -1; element.normal.z = 0.0;
        element.area = 2.859672;
        element.reflectivity = 0.06;
        element.specularity = 0.0;
        element.thickness = 0.0;
        element.absorbtivity = 0.0;
        element.conductivity = 0.0;
        element.emmisivity = 0.0;
        mymodel.busY.push_back(element);
        
        //busZ+
        element.name = "busZ+";
        element.normal.x = 0;element.normal.y = 0.0; element.normal.z = 1.0;
        element.area = 3.059184;
        element.reflectivity = 0.06;
        element.specularity = 0.0;
        element.thickness = 0.0;
        element.absorbtivity = 0.0;
        element.conductivity = 0.0;
        element.emmisivity = 0.0;
        mymodel.busZ.push_back(element);
        
        //busZ-
        element.name = "busZ-";
        element.normal.x = 0;element.normal.y = 0.0; element.normal.z = -1.0;
        element.area = 3.059184;
        element.reflectivity = 0.06;
        element.specularity = 0.0;
        element.thickness = 0.0;
        element.absorbtivity = 0.0;
        element.conductivity = 0.0;
        element.emmisivity = 0.0;
        mymodel.busZ.push_back(element);
        
        //BDS IGSO/MEO
        
        modelstorage["GPSIIR"] = mymodel;
        
        //BDS GEO
        modelstorage["BDSGEO"] = mymodel;
        modelstorage["BDSIGSO"] = mymodel;
        modelstorage["BDSMEO"] = mymodel;
        
        //GLS
        //GAL
        modelstorage["GALIOV"] = mymodel;
        
        return modelstorage;
    }
    
    /* reading the spacecraft model file */
    void GSpacecraftModelMgr::initialiseModel( GString spacecraftModelFile )
    {
        char tmp[2048]="";
        std::ifstream infile(spacecraftModelFile);
        if( !infile )
        {
            std::cout<<"Cannot open file:"<<spacecraftModelFile<<std::endl;
            return;
        }
        
        GSpaceCraftModel::opticalElement element;
        
        std::vector<GString> split;
        while( !infile.eof() )
        {
            memset(tmp,0,sizeof(char)*2048);
            infile.getline(tmp, 2048);
            
            if( tmp[0] == '{')  // start the spacecraft
            {
                GSpaceCraftModel mymodel;
                GString spacecraftName;
                GString srpgridfile[3];
                while(1)
                {
                    memset(tmp,0,sizeof(char)*2048);
                    infile.getline(tmp, 2048);
                    
                    GString str(tmp);
                    
                    if( tmp[0] == '}')
                    {
                        //printf("spacecraftName:%s\n",spacecraftName.c_str());
                        
                        for( int i = 0 ; i< 3; i++ )
                        {
                            if( srpgridfile[i] != "")
                            {
                                //resize the grid data matrix and read in , this is only for solar radiation
                                mymodel.m_srpgriddata[i].resize(GRadiationGridMgr::grid_rows, GRadiationGridMgr::grid_cols);
                                
                                bool test =GRadiationGridMgr::readGridFile(srpgridfile[i], mymodel.m_nomialMass, mymodel.m_srpgriddata[i]);
                                
                                // the same operation needed for infrared radiation here
                            }
                        }
                        
                        spacecraftModelStore[spacecraftName] = mymodel;
                        
                        break;
                    }
                    
                    
                    
                    if(  strstr(tmp,"spacecraftName") )
                    {
                        split = str.split();
                        spacecraftName = split[1];
                        //printf("2\n");
                    }
                    else if(strstr(tmp,"mass")  )
                    {
                        split = str.split();
                        mymodel.m_mass = split[1].asDOUBLE();
                        //printf("3\n");
                    }
                    else if(strstr(tmp,"COM")  )
                    {
                        split = str.split();
                        mymodel.m_com[0] = split[1].asDOUBLE();
                        mymodel.m_com[1] = split[2].asDOUBLE();
                        mymodel.m_com[2] = split[3].asDOUBLE();
                        //printf("4\n");
                    }
                    else if(strstr(tmp,"LRAoffset")  )
                    {
                        split = str.split();
                        mymodel.m_offsetLRA[0] = split[1].asDOUBLE();
                        mymodel.m_offsetLRA[1] = split[2].asDOUBLE();
                        mymodel.m_offsetLRA[2] = split[3].asDOUBLE();
                        //printf("5\n");
                    }
                    
                    else if(strstr(tmp,"GNSSoffset")  )
                    {
                        split = str.split();
                        mymodel.m_offsetGNSS[0] = split[1].asDOUBLE();
                        mymodel.m_offsetGNSS[1] = split[2].asDOUBLE();
                        mymodel.m_offsetGNSS[2] = split[3].asDOUBLE();
                        //printf("6\n");
                    }
                    else if(strstr(tmp,"nominalMass:"))
                    {
                        split = str.split();
                        mymodel.m_nomialMass = split[1].asDOUBLE();
                        //printf("7\n");
                    }
                    
                    else if(strstr(tmp,"antennaPower") )
                    {
                        split = str.split();
                        mymodel.m_antennaPower = split[1].asDOUBLE();
                        //printf("8\n");
                    }
                    else if(strstr(tmp,"panelPowerDraw") )
                    {
                        split = str.split();
                        mymodel.m_solarPanelPowerDraw = split[1].asDOUBLE();
                        //printf("9\n");
                    }
                    else if(strstr(tmp,"busX+"))
                    {
                        split = str.split();
                        element.name = "busX+";
                        element.normal.x = split[1].asDOUBLE();
                        element.normal.y = split[2].asDOUBLE();
                        element.normal.z = split[3].asDOUBLE();
                        element.area = split[4].asDOUBLE();
                        element.thickness = split[5].asDOUBLE();
                        element.absorbtivity = split[6].asDOUBLE();
                        element.conductivity = split[7].asDOUBLE();
                        element.emmisivity = split[8].asDOUBLE();
                        element.reflectivity = split[9].asDOUBLE();
                        element.specularity = split[10].asDOUBLE();
                        
                        mymodel.busX.push_back(element);
                        
                        //printf("10\n");
                        
                    }
                    else if(strstr(tmp,"busX-"))
                    {
                        split = str.split();
                        element.name = "busX-";
                        element.normal.x = split[1].asDOUBLE();
                        element.normal.y = split[2].asDOUBLE();
                        element.normal.z = split[3].asDOUBLE();
                        element.area = split[4].asDOUBLE();
                        element.thickness = split[5].asDOUBLE();
                        element.absorbtivity = split[6].asDOUBLE();
                        element.conductivity = split[7].asDOUBLE();
                        element.emmisivity = split[8].asDOUBLE();
                        element.reflectivity = split[9].asDOUBLE();
                        element.specularity = split[10].asDOUBLE();
                        
                        mymodel.busX.push_back(element);
                        
                        //printf("11\n");
                        
                    }
                    else if(strstr(tmp,"busY+"))
                    {
                        split = str.split();
                        element.name = "busY+";
                        element.normal.x = split[1].asDOUBLE();
                        element.normal.y = split[2].asDOUBLE();
                        element.normal.z = split[3].asDOUBLE();
                        element.area = split[4].asDOUBLE();
                        element.thickness = split[5].asDOUBLE();
                        element.absorbtivity = split[6].asDOUBLE();
                        element.conductivity = split[7].asDOUBLE();
                        element.emmisivity = split[8].asDOUBLE();
                        element.reflectivity = split[9].asDOUBLE();
                        element.specularity = split[10].asDOUBLE();
                        
                        mymodel.busY.push_back(element);
                        
                        //printf("12\n");
                        
                    }
                    else if(strstr(tmp,"busY-"))
                    {
                        split = str.split();
                        element.name = "busY-";
                        element.normal.x = split[1].asDOUBLE();
                        element.normal.y = split[2].asDOUBLE();
                        element.normal.z = split[3].asDOUBLE();
                        element.area = split[4].asDOUBLE();
                        element.thickness = split[5].asDOUBLE();
                        element.absorbtivity = split[6].asDOUBLE();
                        element.conductivity = split[7].asDOUBLE();
                        element.emmisivity = split[8].asDOUBLE();
                        element.reflectivity = split[9].asDOUBLE();
                        element.specularity = split[10].asDOUBLE();
                        
                        mymodel.busY.push_back(element);
                        
                        //printf("13\n");
                        
                    }
                    else if(strstr(tmp,"busZ+"))
                    {
                        split = str.split();
                        element.name = "busZ+";
                        element.normal.x = split[1].asDOUBLE();
                        element.normal.y = split[2].asDOUBLE();
                        element.normal.z = split[3].asDOUBLE();
                        element.area = split[4].asDOUBLE();
                        element.thickness = split[5].asDOUBLE();
                        element.absorbtivity = split[6].asDOUBLE();
                        element.conductivity = split[7].asDOUBLE();
                        element.emmisivity = split[8].asDOUBLE();
                        element.reflectivity = split[9].asDOUBLE();
                        element.specularity = split[10].asDOUBLE();
                        
                        mymodel.busZ.push_back(element);
                        
                        //printf("14\n");
                        
                    }
                    else if(strstr(tmp,"busZ-"))
                    {
                        split = str.split();
                        element.name = "busZ-";
                        element.normal.x = split[1].asDOUBLE();
                        element.normal.y = split[2].asDOUBLE();
                        element.normal.z = split[3].asDOUBLE();
                        element.area = split[4].asDOUBLE();
                        element.thickness = split[5].asDOUBLE();
                        element.absorbtivity = split[6].asDOUBLE();
                        element.conductivity = split[7].asDOUBLE();
                        element.emmisivity = split[8].asDOUBLE();
                        element.reflectivity = split[9].asDOUBLE();
                        element.specularity = split[10].asDOUBLE();
                        
                        mymodel.busZ.push_back(element);
                        
                        //printf("15\n");
                        
                    }
                    else if(strstr(tmp,"solarArray+"))
                    {
                        split = str.split();
                        element.name = "solarArray+";
                        element.normal.x = split[1].asDOUBLE();
                        element.normal.y = split[2].asDOUBLE();
                        element.normal.z = split[3].asDOUBLE();
                        element.area = split[4].asDOUBLE();
                        element.thickness = split[5].asDOUBLE();
                        element.absorbtivity = split[6].asDOUBLE();
                        element.conductivity = split[7].asDOUBLE();
                        element.emmisivity = split[8].asDOUBLE();
                        element.reflectivity = split[9].asDOUBLE();
                        element.specularity = split[10].asDOUBLE();
                        
                        mymodel.solarArray.push_back(element);
                        
                        //printf("16\n");
                        
                    }
                    else if(strstr(tmp,"solarArray-"))
                    {
                        split = str.split();
                        element.name = "solarArray-";
                        element.normal.x = split[1].asDOUBLE();
                        element.normal.y = split[2].asDOUBLE();
                        element.normal.z = split[3].asDOUBLE();
                        element.area = split[4].asDOUBLE();
                        element.thickness = split[5].asDOUBLE();
                        element.absorbtivity = split[6].asDOUBLE();
                        element.conductivity = split[7].asDOUBLE();
                        element.emmisivity = split[8].asDOUBLE();
                        element.reflectivity = split[9].asDOUBLE();
                        element.specularity = split[10].asDOUBLE();
                        
                        mymodel.solarArray.push_back(element);
                        
                        //printf("17\n");
                        
                    }
                    else if(strstr(tmp,"BFSbias"))
                    {
                        split = str.split();
                        mymodel.m_bias.x = split[1].asDOUBLE()*mymodel.m_mass;
                        mymodel.m_bias.y = split[2].asDOUBLE()*mymodel.m_mass;
                        mymodel.m_bias.z = split[3].asDOUBLE()*mymodel.m_mass;
                        
                        //printf("18\n");
                    }
                    else if(strstr(tmp, "gridX"))
                    {
                        split = str.split();
                        if(split.size() == 2 )
                        {
                           srpgridfile[0] = split[1];
                        }
                        
                        //printf("19\n");
                    }
                    else if(strstr(tmp,"gridY"))
                    {
                        split = str.split();
                        if(split.size() == 2 )
                        {
                           srpgridfile[1] = split[1];
                        }
                        //printf("20\n");
                    }
                    else if(strstr(tmp, "gridZ"))
                    {
                        split = str.split();
                        if(split.size() == 2)
                        {
                          srpgridfile[2] = split[1];
                        }
                        //printf("21\n");
                    }
                    
                }
               
                //printf("test_hello2\n");
                
            }
            
            //printf("test_hello1\n");
           
            
            
        }
        
        //printf("finish reading\n");
        
        infile.close();
        
    }
    
    int testc =0;
    
}

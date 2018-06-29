//
//  GConfigure.cpp
//  GFC
//
//  Created by lizhen on 16/04/2017.
//  Copyright © 2017 lizhen. All rights reserved.
//

#include "GConfigure.hpp"
namespace gfc
{
    void GConfigure::parseCfg(gfc::GString cfgfile)
    {
        GString line, param, value;
        int pos = 0;
        const char delim = '=';
        const char harsh = '#';
        
        // Open configurable parameters file
        std::ifstream infile(cfgfile);
        // Check that the file has opened successfully
        if (infile.fail() || !infile.is_open())
        {
            cout<<"configure file open failed!"<<endl;
            return;
        }
        
        
        // Load parameters from file
        while (getline(infile, line))
        {
            // the Basics block
            if(line == "+Basics")
            {
                while( getline(infile,line))
                {
                    if( line == "-Basics" )
                    {
                        break;
                    }
                    //find the "="
                    pos = line.find(delim, 1);
                    if (pos < line.length())
                    {
                        param = line.substr(0, pos);
                        value = line.substr(pos + 1);
                        //strip all the while space in the end and begining
                       
                        pos = value.find(harsh,1);
                        value = value.substr(0,pos);
                        param.strip_v(" ");
                        value.strip_v(" ");
                        
                        if(param == "Spacecraft models") { config.spacecraftmodel = value; }
                        if(param == "Planet Ephemeris") { config.planetEphemeris = value; }
                        if(param == "Earth Gravity file") { config.earthGravityFile = value; }
                        if(param == "EOP") { config.eopfile = value; }
                        if(param == "Integrator Stepsize(sec)") { config.stepsize = value.asDOUBLE(); }
                        
                    }
                    
                    
                }
            } // end of "Basics" block
            
            
            // the Force models block
            if(line == "+Force models")
            {
                while( getline(infile,line))
                {
                    if( line == "-Force models" )
                    {
                        break;
                    }
                    //find the "="
                    pos = line.find(delim, 1);
                    if (pos < line.length())
                    {
                        param = line.substr(0, pos);
                        value = line.substr(pos + 1);
                        //strip all the while space in the end and begining
                       
                        pos = value.find(harsh,1);
                        value = value.substr(0,pos);
                        param.strip_v(" ");
                        value.strip_v(" ");
                        
                        if(param == "Earth Gravity degree")
                        {
                            std::vector<GString> tmp = value.split();
                            if(tmp[0] == "YES")
                            {
                                config.forcelist.push_back("GFMGravity");
                                config.earthGravityDegree = tmp[1].asINT();
                            }
                        }
                        
                        if(param == "Time Variable Earth Gravity")
                        {
                            if(value == "NO")
                            {
                                config.time_variable_gravity = false;
                            }
                            else if(value == "YES")
                            {
                                config.time_variable_gravity = true;
                            }
                        }
                        
                        if(param == "Third body gravity")
                        {
                            std::vector<GString> tmp = value.split();
                            if(tmp[0] == "YES")
                            {
                                config.forcelist.push_back("GFMNbody");
                                
                                //process the planets
                                for(int i = 1 ; i< tmp.size(); i++)
                                {
                                    if(tmp[i] == "EARTH") { config.planetlist.push_back(GJPLEPH::EARTH); }
                                    if(tmp[i] == "SUN") { config.planetlist.push_back(GJPLEPH::SUN); }
                                    if(tmp[i] == "MOON") { config.planetlist.push_back(GJPLEPH::MOON); }
                                    if(tmp[i] == "MERC") { config.planetlist.push_back(GJPLEPH::MERCURY); }
                                    if(tmp[i] == "VENU") { config.planetlist.push_back(GJPLEPH::VENUS); }
                                    if(tmp[i] == "MARS") { config.planetlist.push_back(GJPLEPH::MARS); }
                                    if(tmp[i] == "JUPI") { config.planetlist.push_back(GJPLEPH::JUPITER); }
                                    if(tmp[i] == "SATU") { config.planetlist.push_back(GJPLEPH::SATURN); }
                                    if(tmp[i] == "URAN") { config.planetlist.push_back(GJPLEPH::URANUS); }
                                    if(tmp[i] == "NEPT") { config.planetlist.push_back(GJPLEPH::NEPTUNE); }
                                    if(tmp[i] == "PLUT") { config.planetlist.push_back(GJPLEPH::PLUTO); }
                                }
                            }
                            
                        }
                        if(param == "Solid Earth Tide")
                        {
                            std::vector<GString> tmp = value.split();
                            if(tmp[0] == "YES")
                            {
                                config.solid_earth_tide = true;
                            }
                        }
                        
                        if(param == "Ocean Tide")
                        {
                            std::vector<GString> tmp = value.split();
                            if(tmp[0] == "YES")
                            {
                                config.ocean_earth_tide = true;
                            }
                        }
                        
                        if(param == "Polar Tide")
                        {
                            std::vector<GString> tmp = value.split();
                            if(tmp[0] == "YES")
                            {
                                config.polar_tide = true;
                            }
                        }

                        if(param == "General Relativity")
                        {
                            std::vector<GString> tmp = value.split();
                            if(tmp[0] == "YES")
                            {
                                config.forcelist.push_back("GFMGR");
                            }
                        }
                        if(param == "Antenna Thrust")
                        {
                            std::vector<GString> tmp = value.split();
                            if(tmp[0] == "YES")
                            {
                                config.forcelist.push_back("GFMANT");
                            }
                        }

                        if(param == "Solar Radiation Pressure")
                        {
                            if(value == "NO")
                            {
                                config.srp_option = 0;
                            }
                            else if(value == "BOXW")
                            {
                                config.srp_option = 1;
                                config.forcelist.push_back("GFMSRP");
                            }
                            else if(value == "GRID")
                            {
                                config.srp_option = 2;
                                config.forcelist.push_back("GFMSRP");
                            }
                        }
                        
                        if(param == "Earth Radiation Pressure")
                        {
                            if(value == "NO")
                            {
                                config.erp_option = 0;
                            }
                            else if(value == "BOXW")
                            {
                                config.erp_option = 1;
                                config.forcelist.push_back("GFMERP");
                            }
                            else if(value == "GRID")
                            {
                                config.erp_option = 2;
                                config.forcelist.push_back("GFMERP");
                                
                            }
                        }
                        
                        if(param == "Thermal Radiation")
                        {
                            if(value != "NO")
                            {
                                config.forcelist.push_back("GFMTRR");
                            }
                        }
                        
                        if(param == "BFSbias")
                        {
                            if(value != "NO")
                            {
                                // Ybias force
                                config.forcelist.push_back("GFMBFSbias");
                            }
                        }

                        if(param == "Empirical Force")
                        {
                            if(value != "NO")
                            {
                                config.forcelist.push_back("GFMEMP");
                            }
                            if(value == "ECOM1")
                            {
                                config.emp_option = 0;
                            }
                            else if(value == "ECOM2")
                            {
                                config.emp_option =1;
                            }
                            else if(value == "DREMT")
                            {
                                config.emp_option = 2;
                            }
                            else if(value == "ECOM3")
                            {
                                config.emp_option = 3;
                            }
                        }
                        
                        if(param == "Atmosphere Drag")
                        {
                            if(value != "NO")
                            {
                                config.forcelist.push_back("GFMDRG");
                            }
                        }
                        
                        if(param == "Earth Radiation Model")
                        {
                            if(value == "NO")
                            {
                                config.earth_radiation_option = 0;
                            }
                            else if(value == "SIMPLE")
                            {
                                config.earth_radiation_option = 1;
                            }
                            else if(value == "CERES_T")
                            {
                                config.earth_radiation_option = 2;
                            }
                            else if(value == "CERES_O")
                            {
                                config.earth_radiation_option = 3;
                            }
                        }
                    }
                
                }
            } // end of "Force models" block
            
            
            // the Basics block
            if(line == "+Directories")
            {
                while( getline(infile,line))
                {
                    if( line == "-Directories" )
                    {
                        break;
                    }
                    
                    //find the "="
                    pos = line.find(delim, 1);
                    if (pos < line.length())
                    {
                        param = line.substr(0, pos);
                        value = line.substr(pos + 1);
                        //strip all the while space in the end and begining
                        
                        pos = value.find(harsh,1);
                        value = value.substr(0,pos);
                        param.strip_v(" ");
                        value.strip_v(" ");
                        
                        if(param == "Earth Radiation file")
                        {
                            config.earthRadiationDirectory = value;
                        }
                        if(param == "Max Earth Triangle Level")
                        {
                            config.max_earth_division = value.asINT();
                        }
                    }
                    
                }
            }
            
            
        } // end of while(1)
        
        infile.close();
        
        
        //rest the force partials
        config.forcePartial.resize(config.forcelist.size(),false);
        GEarthRadiationFlux::triFluxPath = config.earthRadiationDirectory;
        GEarthRadiationFlux::maxLevel = config.max_earth_division;
        
        
//        if( config.earth_radiation_option== 0) // no earth radiation model, no erp force,
//        {
//            config.erp_option = 0;
//        }
        
        if(config.earth_radiation_option == 0)
        {
            GEarthRadiationFlux::simple = false;
            GEarthRadiationFlux::ceres_original = false;
            GEarthRadiationFlux::ceres_tri = false;
        }
        else if(config.earth_radiation_option == 1)
        {
            GEarthRadiationFlux::simple = true;
        }
        else if(config.earth_radiation_option == 2)
        {
            GEarthRadiationFlux::ceres_tri = true;
            //load erp and srp grid files
            GEarthRadiationFlux::populateEMData();
        }
        else if(config.earth_radiation_option == 3)
        {
            GEarthRadiationFlux::ceres_original = true;
            // here need to read in the data for 12 months
            GEarthRadiationFlux::populateCERESGrid();
        }
        
    }
    
    
    void GConfigure::configForceModelMgr(GForceModelMgr& modelMgr)
    {
       
        modelMgr.gravity_file = config.earthGravityFile;
        modelMgr.gravity_degree = config.earthGravityDegree;
        modelMgr.gravity_order = config.earthGravityDegree;
        
        modelMgr.setForceList(config.forcelist, config.forcePartial);
        
        // set the bodies for third body force
        //GFMThirdBody::setBodies(config.planetlist);
        
        // update all the planets list
        GSpaceEnv::planetsUsed = config.planetlist;
        
        // GSpaceEnv has to have sun and moon
        if ( std::find(config.planetlist.begin(), config.planetlist.end(), GJPLEPH::SUN) == config.planetlist.end() )
        {
            GSpaceEnv::planetsUsed.push_back(GJPLEPH::SUN);
        }
        
        if ( std::find(config.planetlist.begin(), config.planetlist.end(), GJPLEPH::MOON) == config.planetlist.end() )
        {
            GSpaceEnv::planetsUsed.push_back(GJPLEPH::MOON);
        }
        
        
        for( int i = 0 ; i< config.forcelist.size(); i++ )
        {
            if(config.forcelist[i] == "GFMGravity")
            {
                ((GFMEarthGravity*)modelMgr.m_forceModels["GFMGravity"])->time_variable = config.time_variable_gravity;
                ((GFMEarthGravity*)modelMgr.m_forceModels["GFMGravity"])->solid_tide = config.solid_earth_tide;
                ((GFMEarthGravity*)modelMgr.m_forceModels["GFMGravity"])->ocean_tide = config.ocean_earth_tide;
                ((GFMEarthGravity*)modelMgr.m_forceModels["GFMGravity"])->polar_tide  = config.polar_tide;
            }
            
            if(config.forcelist[i] == "GFMEMP")
            {
                int num_param = 5;
                if(config.emp_option == 0 )  // ECOM1
                {
                    num_param =5;
                }
                else if(config.emp_option == 1) // ECOM2
                {
                    num_param =7;
                }
                else if(config.emp_option == 2) //DREMT
                {
                    num_param = 6;
                }
                
                else if(config.emp_option == 3) //reduced ECOM2
                {
                    num_param = 7;
                }
                
                ((GFMSolarRadiationPressureEM*)modelMgr.m_forceModels["GFMEMP"])->num_param = num_param;
                ((GFMSolarRadiationPressureEM*)modelMgr.m_forceModels["GFMEMP"])->type_opt = config.emp_option;
                ((GFMSolarRadiationPressureEM*)modelMgr.m_forceModels["GFMEMP"])->parameters.resize(num_param);
                ((GFMSolarRadiationPressureEM*)modelMgr.m_forceModels["GFMEMP"])->m_dadp.resize(3, num_param);
                
                
               
                
                
            }
            
            if( config.forcelist[i] == "GFMSRP" )
            {
                if(config.srp_option == 1)  // box-wing
                {
                    ((GFMSolarRadiationPressure*)modelMgr.m_forceModels["GFMSRP"])->with_grid_on = false;
                }
                else if(config.srp_option == 2) // grid file
                {
                    ((GFMSolarRadiationPressure*)modelMgr.m_forceModels["GFMSRP"])->with_grid_on = true;
                }
                
            }
            
            if( config.forcelist[i] == "GFMERP" )
            {
                
                if(config.erp_option == 1)  // box-wing
                {
                    ((GFMEarthRadiationPressure*)modelMgr.m_forceModels["GFMERP"])->with_grid_on = false;
                }
                else if(config.erp_option == 2) // grid file
                {
                    ((GFMEarthRadiationPressure*)modelMgr.m_forceModels["GFMERP"])->with_grid_on = true;
                }
                
            }
            
            if( config.forcelist[i] == "GFMNbody" )
            {
                ((GFMThirdBody*)modelMgr.m_forceModels["GFMNbody"])->setBodies(config.planetlist);
            }
            
        }
        
        
    }
    
    
    
    
}

//============================================================================
//
//  This file is part of GFC, the GNSS FOUNDATION CLASS.
//
//  The GFC is free software; you can redistribute it and/or modify
//  it under the terms of the GNU Lesser General Public License as published
//  by the Free Software Foundation; either version 3.0 of the License, or
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

//
//  GRadiationForce.cpp
//  GFC
//
//  Created by lizhen on 13/06/2016.
//  Copyright © 2016 lizhen. All rights reserved.
//

#include "GRadiationGrid.hpp"

namespace gfc
{
    //std::map< GString, GMatrix* >  GRadiationGridMgr::griddata;
    //double  GRadiationGridMgr::nominal_mass = 1100.0; //kg
    double  GRadiationGridMgr::nominal_flux = 1368.0; // 1368 w/m^2
    
    bool GRadiationGridMgr::readGridFile( GString file_name,double nominal_mass, GMatrix& data )
    {
        std::ifstream grid_file(file_name);
        
        if (grid_file.fail() || !grid_file.is_open())
        {
            std::cerr << "Error opening SRP grid file" << file_name << std::endl;
            return false;
        }
        
        uint_fast16_t n_cols, n_rows; // no. of horizontal and vertical grid nodes.
        double min_lon, max_lon, min_lat, max_lat; // grid file ranges.
        double min_acc, max_acc;
        
        // First header line:
        std::string scrap;
        getline(grid_file, scrap); // DSAA designation from Surfer6 format.
        
        // Check grid resolution:
        grid_file >> n_cols;
        if (n_cols != grid_cols) {
            std::cerr << "Expected " << grid_cols << " columns and detected "
            << n_cols << " in " << file_name << std::endl;
            return false;
        }
        
        grid_file >> n_rows;
        if (n_rows != grid_rows) {
            std::cerr << "Expected " << grid_rows << " columns and detected "
            << n_rows << " in " << file_name << std::endl;
            return false;
        }
        
        // Check grid value ranges:
        grid_file >> min_lon;
        if (min_lon != min_longitude) {
            std::cerr << "Expected minimum longitude of " << min_longitude
            << " and detected " << min_lon << std::endl;
            return false;
        }
        
        grid_file >> max_lon;
        if (max_lon != max_longitude)
        {
            std::cerr << "Expected maximum longitude of " << max_longitude
            << " and detected " << max_lon << std::endl;
            return false;
        }
        
        grid_file >> min_lat;
        if (min_lat != min_latitude)
        {
            std::cerr << "Expected minimum latitude of " << min_latitude
            << " and detected " << min_lat << std::endl;
            return false;
        }
        
        grid_file >> max_lat;
        if (max_lat != max_latitude)
        {
            std::cerr << "Expected maximum latitude of " << max_latitude
            << " and detected " << max_lat << std::endl;
            return false;
        }
        
        grid_file >> min_acc;
        grid_file >> max_acc;
        
        // std::cout << "Minimum acceleration in " << file_name << " is " << min_acc << std::endl;
        // std::cout << "Maximum acceleration in " << file_name << " is " << max_acc << std::endl;
        // std::cout << "\n" << std::endl;
        
        double data_denorm = nominal_mass / (nominal_flux);
        //cout<<file_name.c_str()<<std::endl;
        // Read in grid node values.
        // Populate grid with decreasing row index - maintaining the Surfer6 format (upside Mercator).
        //data.resize(grid_rows);
        for ( uint_fast16_t i = 0; i < grid_rows; ++i ) //181
        {
            //data[i].resize(grid_cols);
            for ( uint_fast16_t j = 0; j < grid_cols; ++j ) //361
            {
                grid_file >> data(i,j);
                
                data(i,j) *= data_denorm;
               
                //printf("%.16E ", data(i,j) );
            }
            //printf("\n");
        }
        
        //cout<<data;
        
        grid_file.close();
        
        return true;
        
    }
    
//    void GRadiationGridMgr::loadGridFile( GString filename[3] , GString svType )
//    {
//        
//        if( svType == "GPSIIR")
//        {
//            //nominal_mass = 1100.0; // for gpsIIR
//            
//            nominal_flux = 1368.0; // watta
//        }
//        
//        
//        //note that, assign memory here
//        griddata[svType] = new GMatrix[3];
//        
//        for( int i = 0 ; i< 3; i++ )
//        {
//            
//            griddata[svType][i].resize(grid_rows, grid_cols);
//            
//            readGridFile(filename[i], griddata[svType][i]);
//            
//        }
//    }
    
    
    /*
     the biliear interpolation function
     latitude: -90 - 90
     longitude: -180 - 180
     //  https://en.wikipedia.org/wiki/Bilinear_interpolation
     */
    GVector GRadiationGridMgr::mybilinear_interp(GMatrix (&griddata)[3], double longitude, double latitude)
    {
        //first find out the corresponding griddata
        GMatrix& x_grid = griddata[0];
        GMatrix& y_grid = griddata[1];
        GMatrix& z_grid = griddata[2];
        
        GVector acc;
        
        //double min_longitude = -180, max_longitude = 180;
        //double min_latitude = -90, max_latitude = 90;
        
        // Check input for validity.
        if (longitude < min_longitude)
        {
            longitude = min_longitude;
        }
        else if (longitude > max_longitude)
        {
            longitude = max_longitude;
        }
        
        if (latitude < min_latitude)
        {
            latitude = min_latitude;
        }
        else if (latitude > max_latitude)
        {
            latitude = max_latitude;
        }
        
        
        //find the minimum grid first;
        // the grid is 1 by 1 degree
        double u_interval = 1.0;
        double v_interval = 1.0;
        
        int u = floor((latitude - min_latitude)/u_interval);  // latitude
        int v = floor((longitude - min_longitude)/v_interval); // longitude
        
        int u_max = floor((max_latitude - min_latitude)/u_interval);
        int v_max = floor((max_longitude - min_longitude)/v_interval);
        
        //make del_u and del_v both less than 1.0
        double del_u = (latitude - min_latitude -  u*u_interval)/u_interval;
        double del_v = (longitude - min_longitude - v*v_interval)/v_interval;
        
        
        double a_u1 = 0.0, a_u2 =0.0;
        
        double grid_x_u_v =0.0, grid_x_u_1_v =0, grid_x_u_v_1 =0,  grid_x_u_1_v_1 =0;
        double grid_y_u_v =0.0, grid_y_u_1_v =0, grid_y_u_v_1 =0,  grid_y_u_1_v_1 =0;
        double grid_z_u_v =0.0, grid_z_u_1_v =0, grid_z_u_v_1 =0,  grid_z_u_1_v_1 =0;
        
        grid_x_u_v = x_grid(u,v);
        grid_y_u_v = y_grid(u,v);
        grid_z_u_v = z_grid(u,v);
        
        if(u+1<= u_max)
        {
            grid_x_u_1_v = x_grid(u+1,v);
            grid_y_u_1_v = y_grid(u+1,v);
            grid_z_u_1_v = z_grid(u+1,v);
        }
        
        if(v+1 <= v_max)
        {
            grid_x_u_v_1 =x_grid(u,v+1);
            grid_y_u_v_1 =y_grid(u,v+1);
            grid_z_u_v_1 =z_grid(u,v+1);
        }
        if(v+1 <= v_max && u+1<= u_max )
        {
            grid_x_u_1_v_1 = x_grid(u+1,v+1);
            grid_y_u_1_v_1 = y_grid(u+1,v+1);
            grid_z_u_1_v_1 = z_grid(u+1,v+1);
        }
        
        
        a_u1 = (1.0 - del_u)*grid_x_u_v + del_u * grid_x_u_1_v;
        a_u2 = (1.0 - del_u)*grid_x_u_v_1 + del_u * grid_x_u_1_v_1;
        acc.x = (1.0 - del_v)*a_u1  + del_v * a_u2;
        
        a_u1 = (1.0 - del_u)*grid_y_u_v + del_u * grid_y_u_1_v;
        a_u2 = (1.0 - del_u)*grid_y_u_v_1 + del_u *grid_y_u_1_v_1;
        acc.y = (1.0 - del_v)*a_u1  + del_v * a_u2;
        
        a_u1 = (1.0 - del_u)*grid_z_u_v + del_u * grid_z_u_1_v;
        a_u2 = (1.0 - del_u)*grid_z_u_v_1 + del_u *grid_z_u_1_v_1;
        acc.z = (1.0 - del_v)*a_u1  + del_v * a_u2;
        
        return acc;
    }
    
    GVector GRadiationGridMgr::bilinear_interp(GMatrix (&griddata)[3], double longitude, double latitude)
    {
        //first find out the corresponding griddata
        GMatrix& x_grid = griddata[0];
        GMatrix& y_grid = griddata[1];
        GMatrix& z_grid = griddata[2];
        
        GVector acc;
        // Check input for validity.
        if (longitude < min_longitude) {
            longitude = min_longitude;
        } else if (longitude > max_longitude) {
            longitude = max_longitude;
        }
        
        if (latitude < min_latitude) {
            latitude = min_latitude;
        } else if (latitude > max_latitude) {
            latitude = max_latitude;
        }
        
        // Raw horizontal and vertical indexed in array (0 indexing).
        double u_raw = longitude - min_longitude;
        double v_raw = latitude - min_latitude;
        
        double u_floor = std::floor(u_raw);
        double v_floor = std::floor(v_raw);
        
        double del_u = u_raw - u_floor;
        double del_v = v_raw - v_floor;
        
        double tol = 1.0E-9; // ~10^-9
        
        uint_fast16_t i1, i2, j1, j2;
        
        i1 = static_cast<uint_fast16_t>(v_floor);
        i2 = i1 + 1;
        j1 = static_cast<uint_fast16_t>(u_floor);
        j2 = j1 + 1;
        
        // Perform interpolation based on appropriate neighbors:
        
        if (del_u > tol) {
            if (del_v > tol) {
                // Interpolation point falls between grid points in both directions:
                double row1, row2;
                
                // Horizontal interpolation along upper border.
                row1 = (x_grid(i1,j2) - x_grid(i1,j1) ) * del_u + x_grid(i1,j1);
                // Horizontal interpolation along lower border.
                row2 = (x_grid(i2,j2) - x_grid(i2,j1) ) * del_u + x_grid(i2,j1);
                // Vertical interpolation (downward) between row1 and row2.
                acc.x = (row2 - row1) * del_v + row1;
                
                row1 = (y_grid(i1,j2) - y_grid(i1,j1) ) * del_u + y_grid(i1,j1);
                row2 = (y_grid(i2,j2) - y_grid(i2,j1) ) * del_u + y_grid(i2,j1);
                acc.y = (row2 - row1) * del_v + row1;
                
                row1 = (z_grid(i1,j2) - z_grid(i1,j1) ) * del_u + z_grid(i1,j1);
                row2 = (z_grid(i2,j2) - z_grid(i2,j1)) * del_u + z_grid(i2,j1);
                acc.z = (row2 - row1) * del_v + row1;
            } else {
                // Interpolation point (pretty much) falls on horizontal grid line:
                // 1D rightward interpolation (increasing longitude).
                acc.x = (x_grid(i1,j2) - x_grid(i1,j1)) * del_u + x_grid(i1,j1);
                acc.y = (y_grid(i1,j2) - y_grid(i1,j1)) * del_u + y_grid(i1,j1);
                acc.z = (z_grid(i1,j2) - z_grid(i1,j1)) * del_u + z_grid(i1,j1);
            }
        } else {
            if (del_v > tol) {
                // Interpolation point (pretty much) falls on vertical grid line:
                // 1D downward interpolation (increasing latitude).
                acc.x = (x_grid(i2,j1) - x_grid(i1,j1)) * del_v + x_grid(i1,j1);
                acc.y = (y_grid(i2,j1) - y_grid(i1,j1)) * del_v + y_grid(i1,j1);
                acc.z = (z_grid(i2,j1) - z_grid(i1,j1)) * del_v + z_grid(i1,j1);
            } else {
                // Interpolation point falls on grid node exactly (ish):
                acc.x = x_grid(i1,j1);
                acc.y = y_grid(i1,j1);
                acc.z = z_grid(i1,j1);
            }
        }
        
        return acc;
        
    } // end of bilinear_interp
    
    // please note that force is in BFS
    //
    GVector GRadiationGridMgr::getAcc_grid(GMatrix (&griddata)[3], double longitude, double latitude)
    {
        GVector force;
        
        // the results is newton per flux
        //force = bilinear_interp(griddata, longitude, latitude);
        
        force =  mybilinear_interp(griddata, longitude, latitude);
        
        
        return force;
    }
    
    
} // end of namespace gfc

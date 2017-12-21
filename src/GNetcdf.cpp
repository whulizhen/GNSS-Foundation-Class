//
//  GNetcdf.cpp
//  GFC
//
//  Created by lizhen on 16/1/29.
//  Copyright © 2016年 lizhen. All rights reserved.
//



#include "GNetcdf.h"


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


namespace gfc
{
    
    
    GString GNetcdf::DataType[13] = {"NC_NAT","BYTE","CHAR","SHORT","INT","FLOAT","DOUBLE",
                                     "UBYTE","USHORT","UINT","INT64","UINT64","STRING"};
    
    GNetcdf::GNetcdf(GString filename, int mode)
    {
        loadNCFile(filename, mode);
    }
    
    void GNetcdf::handle_error(int status)
    {
         if (status != NC_NOERR)
         {
            fprintf(stderr, "%s\n", nc_strerror(status));
            exit(-1);
         }
    }
    
    void GNetcdf::loadNCFile(GString filename, int mode)
    {
        int status =   nc_open(filename.c_str() , mode, &m_fileID);
        handle_error(status);
        char names[MAXCHAR_NAME]={0};
        m_ncver = nc_inq_libvers();
        // call the nc functions in netcdf library
        status = nc_inq( m_fileID, &m_ndim, &m_nvar, &m_nattGLOBAL, NULL);
        handle_error(status);
        int * unlimitedAttri = new int[m_ndim];
        memset(unlimitedAttri,-1, sizeof(int)*m_ndim);
        status = nc_inq_unlimdim( m_fileID, unlimitedAttri); handle_error(status);
        // get the names of all the dimensions
        
        for( int i = 0 ; i< m_ndim ; i++ )
        {
            memset( names,0,sizeof(char)*MAXCHAR_NAME);
            size_t len =0;
            status = nc_inq_dim(m_fileID, i, names, &len);
            handle_error(status);
            ncDimension mydim;
            mydim.m_dimName = names;
            mydim.m_dimVal  = len;
            mydim.m_ID = i;
            if( unlimitedAttri[i] != -1 )
            {
                mydim.m_isUnlimited = true;
            }
            
            mydim.m_data.resize(len);
            
            status = nc_get_var_double(m_fileID, i, &mydim.m_data[0]);
            handle_error(status);
            if( mydim.m_data.size() > 1 )
            {
                mydim.m_step = mydim.m_data[1] - mydim.m_data[0];
            }
            else if(mydim.m_data.size() == 0 )
            {
                mydim.m_step = 0;
            }
            // the method can be very slow.
            //get the maxmum and minimum in this dimension data
            mydim.m_maxmum =  *std::max_element(mydim.m_data.begin(), mydim.m_data.end()--);
            mydim.m_minimum =  *std::min_element(mydim.m_data.begin(), mydim.m_data.end()--);
            if( mydim.m_isUnlimited == true)
            {
                mydim.m_maxmum = mydim.m_minimum;
            }
            
            m_dimensions.push_back(mydim);
        }
        
        if( unlimitedAttri != NULL) {delete[] unlimitedAttri; unlimitedAttri = NULL;}
        
        //get information on the global attributes
       
        for( int i = 0 ; i< m_nattGLOBAL; i++ )
        {
            memset(names,0,sizeof(char)*MAXCHAR_NAME);
            nc_type mytype= -1 ;
            size_t len = -1 ;
            status = nc_inq_attname( m_fileID, NC_GLOBAL, i, names);
            handle_error(status);
            status = nc_inq_att ( m_fileID, NC_GLOBAL, names, &mytype, &len);
            handle_error(status);
            char attrVal[MAXCHAR_VAL]={0};
            //char * attrVal = new char[len*2];
            status = nc_get_att(m_fileID, NC_GLOBAL, names, attrVal);
            handle_error(status);
            ncAttribute myattri;
            myattri.m_attName = names;
            myattri.m_ID = i;
            myattri.m_dataType = mytype;
            myattri.m_value = attrVal;
            
            m_attriGLOBAL.push_back(myattri);
            
        }
        
        int* ndimids = new int[m_ndim];
        memset(ndimids,-1,sizeof(int)*m_ndim);
        // get information about all the Variables
        for( int i = 0 ; i< m_nvar ; i++ )
        {
            memset(names,0,sizeof(char)*MAXCHAR_NAME);
            nc_type xtype = 0;
            int ndims = 0;
            int natts = 0;
            status =  nc_inq_var( m_fileID,i,names,&xtype,&ndims,ndimids,&natts);
            handle_error(status);
            
            // judge whether current variable is in the dimensions set
            // yes, then add some attributes,
            // no,  then doing nothing
            GString tmpstr(names);
            bool test = false;
            int  index = -1;
            for( int k = 0 ; k< m_ndim ; k++ )
            {
                if( m_dimensions[k].m_dimName == tmpstr )
                {
                    test = true;
                    index = k ;
                    break;
                }
            }
            
            ncVariable myvar;
            for(int j = 0 ; j< m_ndim ; j++)
            {
                if(ndimids[j] != -1)
                {
                    myvar.m_dimen.push_back(ndimids[j]);
                }
            }
            
            myvar.m_varName = names;
            myvar.m_varID = i;
            myvar.m_dataType = xtype;
            
            // get the attributes just belond to this special variable
            for( int j = 0 ; j< natts ; j++ )
            {
                nc_type xtype = -1;
                size_t len = -1 ;
                memset(names,0,sizeof(char)*MAXCHAR_NAME);
                status = nc_inq_attname(m_fileID, i, j, names);
                handle_error(status);
                status = nc_inq_att(m_fileID, i, names, &xtype, &len);
                char attrVal[MAXCHAR_VAL] ={0};
                status = nc_get_att(m_fileID, i, names, attrVal);
                
                ncAttribute myatt;
                myatt.m_attName = names;
                myatt.m_ID = j ;
                myatt.m_value = attrVal;
                myatt.m_value = myatt.m_value.substr(0,len);
                myatt.m_dataType = xtype;
                if( test == true )
                {
                    m_dimensions[index].m_attri.push_back(myatt) ;
                }
                else if( test == false )
                {
                   myvar.m_attribute.push_back(myatt);
                }
            }
            
            if( test == false )
            {
               m_variables.push_back(myvar);
            }
            
        }
        
        if(ndimids != NULL )  {delete[] ndimids; ndimids = NULL;}
        m_nvar = m_variables.size();
        
        
        printf("Data Information:\n");
        printf("Dim_num: %d Att_num: %d Var_num: %d\n",m_ndim, m_nattGLOBAL,m_nvar);
        printf("Dimension Information:\n");
        for(int i =0 ; i< m_ndim ; i++ )
        {
            printf("Dim_Name:%s len: %d Unlimited: %d\n",
                    m_dimensions[i].m_dimName.c_str(),m_dimensions[i].m_dimVal,
                    m_dimensions[i].m_isUnlimited
                  );
            for(int j = 0 ; j<m_dimensions[i].m_attri.size() ; j++ )
            {
                printf("%s : %s\n",m_dimensions[i].m_attri[j].m_attName.c_str(),
                       m_dimensions[i].m_attri[j].m_value.c_str());
            }
        }
        
        printf("****************************************************************************\n");
        printf(" Global Attribute Information :\n");
        for(int i = 0 ; i< m_nattGLOBAL ; i++ )
        {
            printf("%s : %s \n",m_attriGLOBAL[i].m_attName.c_str(),m_attriGLOBAL[i].m_value.c_str());
            //std::cout<<m_attriGLOBAL[i].m_attName<<" : "<<m_attriGLOBAL[i].m_value<<std::endl;
        }
        
        printf("****************************************************************************\n");
        printf(" Variables Information :\n");
        for(int i = 0 ; i< m_variables.size() ; i++ )
        {
            printf("%s dimension: ",m_variables[i].m_varName.c_str());
            for(int j = 0 ; j< m_variables[i].m_dimen.size(); j++ )
            {
                printf("%d ",m_variables[i].m_dimen[j]);
            }
            printf("\n");
            for(int j =0; j< m_variables[i].m_attribute.size() ; j++ )
            {
                printf("%s : %s \n",m_variables[i].m_attribute[j].m_attName.c_str(),m_variables[i].m_attribute[j].m_value.c_str());
            }
        }
        
       // int testc = 0;
    }
    
    
    /*get Attribute Value in the given name*/
    GString GNetcdf::getAttriValue( gfc::GString attriName, void* value)
    {
        GString data_type_string;
        // first for Global attributes
        
        // second for dimensions
        
        
        return data_type_string;
    }
    
    
    std::vector<double> GNetcdf::getDimensionData(gfc::GString varName )
    {
        int index =-1;
        int dimID =-1;
        for( int i = 0 ; i<m_dimensions.size(); i++ )
        {
            if( m_dimensions[i].m_dimName == varName)
            {
                index = i;
                dimID = m_dimensions[i].m_ID;
            }
        }
        
        return m_dimensions[index].m_data;
    }
    
  
    /* 
     * just get one data according to the index
     * index is a array which describling the position in the data
     * index[0], the first dimension , index[1], the second dimension
     * be careful that the lenght of index must equal to the number of dimensions of this data
     */
    double GNetcdf::getData(size_t* indexArray, GString varName)
    {
        double retval =0.0;
        int varid = -1;
        int status =-1;
        int mypos =-1;
        for( int i = 0 ; i< m_variables.size() ; i++ )
        {
            if( m_variables[i].m_varName == varName)
            {
                mypos = i;
                varid = m_variables[i].m_varID;
                break;
            }
        }
        
        if( mypos == -1 )
        {
            printf("can not find the variable: %s\n",varName.c_str());
            return 0.0;
        }
        
        status = nc_get_var1_double( m_fileID, varid,indexArray, &retval);
        handle_error(status);
        return retval;
    }
    
    int GNetcdf::varName2ID(GString varName, int* index)
    {
        if( index != NULL ) { *index = -1;}
        int varid =  -1;
        int status = -1;
        for( int i = 0 ; i< m_variables.size() ; i++ )
        {
            if( m_variables[i].m_varName == varName)
            {
                if( index != NULL)
                {
                    *index = i;
                }
                
                varid = m_variables[i].m_varID;
                
                break;
            }
        }
        
        if( *index == -1 )  { printf("can not find the variable: %s\n",varName.c_str());}
        return varid;
    }
    
    void GNetcdf::getData(GString varName, double* data)
    {
        int varid = varName2ID(varName, NULL);
        int status = nc_get_var_double(m_fileID, varid, data);
        handle_error(status);
    }
    
   /* !!get data that you want!! 
    * start,count are arrays with the same size of number of dimension
    * examples: start={3, 5, 0}, count={1,1,360}; get the 3th dimension data with 1st dimension == 3th and 2nd dimension == 5th
    * the minimum of count is 1
    */
  void GNetcdf::getData( GString varName, size_t* start, size_t* count, double* data)
    {
        int index = -1;
        int varid = -1;
        int status = -1;
        varid = varName2ID( varName, &index );
        status = nc_get_vara_double( m_fileID, varid,start,count, data);
        handle_error(status);
    }
    
    
    //just for test
  void GNetcdf::getData( double *data, gfc::GString varName)
    {
        int index = -1;
        int varid = -1;
        int status = -1;
        varid = varName2ID(varName, &index);
        
        size_t start[3]={0}, count[3]={0};
        
//        count[0] = 1; count[1] = 1;     count[2] = 360;
//        start[0] = 3;   start[1] = 5;     start[2] = 0;
//        
//        double mytestData[360];
//        double toadatatest[182][360];
//        status = nc_get_vara_double( m_fileID, varid,start,count, mytestData);
//        
//        for( int i = 0 ; i< 360; i++ )
//        {
//            start[2] = i;
//            int  status = nc_get_vara_double( m_fileID, varid,start,count, &toadatatest[0][0]);
//            int testc =0;
//        }
        
        count[0] = 1; count[1] = 180; count[2] = 360;
        start[0] = 0; start[1] = 0;   start[2] = 0;
        double toadata[180][360];
        
        CivilTime ct( 2000,3,1,0,0,0 ,"ts"); //starter time
        GTime  mytime;
        mytime.SetFromCivilTime(ct);
        
        std::vector<double> timeData;
        
        double secpday = GCONST("SECPDAY");
        
        struct testData
        {
            double data[180][360];
            testData()
            {
                memset(data,0,sizeof(double)*180*360);
            }
        };
        
        std::vector<testData> totalData;
        totalData.reserve(12);
        
        int    testCount[12]={0};//纪录每个月的次数
        timeData = getDimensionData("time");
        
        FILE* mydatafile = fopen("toa.data", "w+");
        for( int rec = 0 ; rec < 182; rec++ )
        {
            
            start[0] = rec;
            int  status = nc_get_vara_double( m_fileID, varid,start,count, &toadata[0][0]);
            handle_error(status);
            GTime curtime ;
            curtime.SetData(TimeSystem::GetByName("tsUKN"), timeData[rec], 0, 0);
            curtime = mytime + curtime;
            JDTime jt = GTime::GTime2JDTime(curtime);
            CivilTime  myct = GTime::JDTime2CivilTime(jt);
            
            testCount[myct.m_month-1]++;
            for(int i = 0 ; i< 180; i++)
            {
                for( int j =0; j< 360; j++ )
                {
                    totalData[myct.m_month-1].data[i][j]+=toadata[i][j];
                }
            }
            
            fprintf(mydatafile, "time: %04d %02d %02d %02d %02d %4.2f\n",myct.m_year,myct.m_month,myct.m_day,myct.m_hour,myct.m_minute,myct.m_second);
            for( int i = 0 ; i< 360; i++ )
            {
                for( int j =0; j< 180 ; j++ )
                {
                    fprintf(mydatafile,"%8.3f ", toadata[j][i]);
                }
                fprintf(mydatafile, "\n");
            }
            
        }
        fclose(mydatafile);
       
        
        FILE* meanFile = fopen("mean.data","w+");
        for(int k =0 ; k< totalData.size() ; k++ )
        {
            fprintf(meanFile, "%Month: %02d\n",k+1);
            for( int i = 0 ; i<360; i++ )
            {
                for( int j = 0 ; j< 180; j++ )
                {
                    totalData[k].data[j][i] = totalData[k].data[j][i] / testCount[k];
                    fprintf(meanFile, "%8.3f ",totalData[k].data[j][i]);
                }
                fprintf(meanFile, "\n");
            }
        }
        fclose(meanFile);
        
        
        
    }
    
    
    
} // end of namespace gfc



//
//  GOrbitFitting.cpp
//  GFC
//
//  Created by lizhen on 10/07/2016.
//  Copyright © 2016 lizhen. All rights reserved.
//

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


#include "GOrbitFitting.hpp"
namespace gfc
{
    GOrbitFitting::GOrbitFitting()
    {
        m_fittingRMS = 0.0;
    }
    
    GOrbitFitting::GOrbitFitting(GSpaceCraft* spacecraft)
    {
        m_fittingRMS = 0.0;
        m_integrator = NULL;
        m_integrator = new GRungeKuttaFehlberg();
        setSpaceCraft(spacecraft);
    }
    
    void GOrbitFitting::setStepSize(double stepsize)
    {
        m_integrator->setStepsize(stepsize);
    }
    
    void GOrbitFitting::setSpaceCraft(gfc::GSpaceCraft *spacecraft)
    {
        m_spaceCraft = spacecraft;
    }
    
    
    GForceModelMgr GOrbitFitting::getForceManager()
    {
        return m_forceManager;
    }
    
    void GOrbitFitting::setForceManager(GForceModelMgr manager)
    {
        m_forceManager = manager;
        
        std::map< GString, GForceModel* >::iterator   it = m_forceManager.m_forceModels.begin();
        m_forceManager.totalNum_parameters = 0; // recount the num of parameters
        
        // as a orbit fitter, need to calculate all the derivatives and estimating parameters
        for( ;it != m_forceManager.m_forceModels.end(); ++it )
        {
            it->second->hasPartialDerivatives(true);
            
            m_forceManager.totalNum_parameters += it->second->num_param;
        }
        
        m_spaceCraft->getStatePointer()->senMatrix.resize(6, m_forceManager.totalNum_parameters);
        
    }
    
    void GOrbitFitting::setInitialValue(gfc::GTime t, gfc::GVector p, gfc::GVector v)
    {
        m_initialEpoch = t;
        m_initialP = p;
        m_initialV = v;
    }
    
    void GOrbitFitting::setEndEpoch(gfc::GTime &t)
    {
        m_endEpoch = t;
    }
    
    
    // this is the function fiting the orbit with
    // parameters estimated, mainly SRP parameters and satellite mass changes
    // thus, need to use kalman filtering
    
    /* for this kalman filter:
     measurement update:   L = B * Xt
     time        update:   Xt  = PHI(t,t-1)*Xt-1 + Qk
    */
    void GOrbitFitting::fitting()
    {
        
        GVector po_ecf,vo_ecf,po_eci,vo_eci,pc;
        double phiI[36]={
            1,0,0,0,0,0,
            0,1,0,0,0,0,
            0,0,1,0,0,0,
            0,0,0,1,0,0,
            0,0,0,0,1,0,
            0,0,0,0,0,1
        };
        
        int nobs = 50;
        
        GMatrix X(6,1); // the solution, correction to the initial position and velocity
        
        GMatrix P(6,6); // the covariance matrix of state vector X
        P(0,0) = 2*2; P(1,1) = 2*2; P(2,2) = 2*2; // position accurancy: 100 meter
        P(3,3) = 0.01*0.01; P(4,4) = 0.01*0.01; P(5,5) = 0.01*0.01; // velocity accurancy: 1dm
        
        
        GMatrix I(phiI,6,6); // the design matrix for measurement update
        GMatrix B = I;
        
        GMatrix PHI(6,6); // transition matrix from t0 to tn
        PHI = I ;
        
        GMatrix PHIK(6,6);  // the transition matrix for kalman filter
        PHIK = I;
        
        GMatrix R(6,6); // the covariance matrix of observation, 5 cm for position and 1mm for velocity
        GMatrix Q(6,6); // the state noise matrix of time update equation, it should include the mismodelling acceleration error
        double delta0[6] = {0.05, 0.05, 0.05, 0.005,0.005,0.005 };
        R(0,0) = delta0[0]*delta0[0]; R(1,1) = delta0[1]*delta0[1]; R(2,2) = delta0[2]*delta0[2];
        R(3,3) = delta0[3]*delta0[3]; R(4,4) = delta0[4]*delta0[4]; R(5,5) = delta0[5]*delta0[5];
        
        
        m_initialP = m_spaceCraft->getStatePointer()->satpos_eci;
        m_initialV = m_spaceCraft->getStatePointer()->satvel_eci;
        
        m_spaceCraft->getStatePointer()->updateState_eci(m_initialEpoch, m_initialP, m_initialV);
        
        GEarthOrientationParameter eop;
        //*  2004  2 29 21  0  0.00000000
        for( int i = 0 ; i< nobs ; i++ )
        {
            GTime ct_utc = m_initialEpoch + 900*(i+1);
            GTime ct_gps = GTime::UTC2GPST(ct_utc);
            //this is the observation
            //m_spaceCraft->getEphValue(index).getPV(po_ecf,vo_ecf); // these are in ECEF, should be transformed to ECI
            m_spaceCraft->getEphValue(ct_gps).getPV(po_ecf,vo_ecf);
            
            eop.setEpochTime(ct_utc);
            eop.ECEF2ECI_pos(po_ecf, po_eci);
            eop.ECEF2ECI_vel(po_ecf, vo_ecf, vo_eci);
            
            PropagateTo( ct_utc  ); //here ct must be in UTC
            
            //get the calculation value
            pc = m_spaceCraft->getStatePointer()->satpos_eci;
            
            double obs[6] = { po_eci.x - pc.x,
                              po_eci.y - pc.y,
                              po_eci.z - pc.z,
                              vo_eci.x - m_spaceCraft->getStatePointer()->satvel_eci.x,
                              vo_eci.y - m_spaceCraft->getStatePointer()->satvel_eci.y,
                              vo_eci.z - m_spaceCraft->getStatePointer()->satvel_eci.z
                            };
            
            
            //robust obs check !!!
            
            
            GMatrix Z(obs,6,1);
            Z = Z *1000.0;  // unit: meter
            PHI = m_spaceCraft->getStatePointer()->phiMatrix * PHI ;
            
            B = PHI;
            GMatrix X1 = PHIK*X;
            GMatrix P1 = PHIK*P*(~PHIK) + Q;
            GMatrix TM = ( B*P1*(~B) + R);
            GMatrix Kk =   P1*(~B)*(!TM); //增益矩阵
            
            // estimated resulte
            X = X1 + Kk*(Z - B*X1 );
            P = (I - Kk*B)*P1;
            
            JDTime jd_ct = GTime::GTime2JDTime(ct_utc);
            
//            cout<< GTime::JDTime2CivilTime(jd_ct).TimeString()<<","
//                <<(~X)[0]*1000.0<<","<<(~X)[1]*1000.0<<","<<(~X)[2]*1000.0<<","
//                <<(~X)[3]*1000.0<<","<<(~X)[4]*1000.0<<","<<(~X)[5]*1000.0<<","<<endl;
//            
            printf("%s, %.12f, %.12f, %.12f, %.12f, %.12f, %.12f,\n",
                  GTime::JDTime2CivilTime(jd_ct).TimeString().c_str(),
                 (~X)[0]*1000.0,(~X)[1]*1000.0,(~X)[2]*1000.0,
                 (~X)[3]*1000.0,(~X)[4]*1000.0,(~X)[5]*1000.0);
            
           // cout<<P;
            
            
            int testc = 0;
        }
        
        printf("%.12f ,%.12f ,%.12f, %.12f, %.12f, %.12f\n",
              (~X)[0]/1000.0 + m_initialP.x ,
               (~X)[1]/1000.0+ m_initialP.y,
               (~X)[2]/1000.0+ m_initialP.z,
               (~X)[3]/1000.0+ m_initialV.x,
               (~X)[4]/1000.0+ m_initialV.y,
               (~X)[5]/1000.0+ m_initialV.z);
        
        int testc = 0;
        
    }
    
    
    void GOrbitFitting::calculate_SRPModel()
    {
        std::vector<double> phi; // latitude of sun in bfs
        std::vector<double> eta; // the eta angle for ECOM
        
        std::vector<GSpaceCraftAttitude> attitude;
        std::vector<GVector> acc_eci;
        
        //ouput to file
        FILE* fout = fopen(m_outputFile.c_str(), "w+");
        
        if(!fout)
        {
            printf("file %s open failed!\n", m_outputFile.c_str());
        }
        
        GSensorID sensorID = m_spaceCraft->getSpaceCraftID();
        fprintf(fout, "Reference Time: %s \t %s\n",GTime::GTime2CivilTime(m_initialEpoch).TimeString().c_str(),m_initialEpoch.getTimeSystem().getTimeSystemName().c_str());
        
        fprintf(fout, "%s\t%3d\n",sensorID.getSystem().c_str(),sensorID.getIDnum());
        
        GVector po_ecf,vo_ecf,po_eci,vo_eci,pc,vc;
        
        double average_ele= 0.0;
        
        double interval = 900;  // interval :900 seconds, this may need to be changed !
        
        int iteration_num = 0;
        
        //m_spaceCraft->getNumofPreEph();
        
        int np = m_forceManager.totalNum_parameters;
        
        //only estimate the SRP model parameters
        GMatrix X(np,1); // the solution, correction to the initial position and velocity
        
        GMatrix BTPB(np,np);
        GMatrix BTPL(np,1);
        GMatrix PHI(6,6); // transition matrix from t0 to tn
        
        GMatrix PHIpos(3,6);
        
        //sensitivity matrix
        GMatrix S(6,np);
        GMatrix Spos(3,np);
        
        
        GMatrix W(6,6);  // the weighting matrix
        W(0,0) = 1;W(1,1) = 1;W(2,2) = 1;
        //W(3,3) = 1.0/(100000*100000);W(4,4) = 1.0/(100000*100000);W(5,5) = 1.0/(100000*100000);
        W(3,3) = 1;W(4,4) = 1;W(5,5) = 1;
        
        // temp matrice
        GMatrix PHITPHI(6,6),PHITS(6,np),STPHI(np,6),STS(np,np),PHITL(6,1),STL(np,1);
        
        m_initialP = m_spaceCraft->getStatePointer()->satpos_eci;
        m_initialV = m_spaceCraft->getStatePointer()->satvel_eci;
        int iteration = 0;
        double phiI[36]={ 1,0,0,0,0,0,
            0,1,0,0,0,0,
            0,0,1,0,0,0,
            0,0,0,1,0,0,
            0,0,0,0,1,0,
            0,0,0,0,0,1
        };
        
        //double ecomParam[5]={ -0.00000010042745184630757, -0.00000000020611988712060707, 0.0000000016151040633494485,-0.0000000004195205808724743,0.0000000031361947345823236};
        // they are in units of m/s^2
        double ecomParam[10]={0.0};
        fprintf(fout, "++ a priori value\n");
        fprintf(fout, "PX:%20.8f\n",m_initialP.x);
        fprintf(fout, "PY:%20.8f\n",m_initialP.y);
        fprintf(fout, "PZ:%20.8f\n",m_initialP.z);
        fprintf(fout, "VX:%20.8f\n",m_initialV.x);
        fprintf(fout, "VY:%20.8f\n",m_initialV.y);
        fprintf(fout, "VZ:%20.8f\n",m_initialV.z);
        
        fprintf(fout, "D0:%20.8E\n",ecomParam[0]);
        fprintf(fout, "Y0:%20.8E\n",ecomParam[1]);
        fprintf(fout, "B0:%20.8E\n",ecomParam[2]);
        fprintf(fout, "Bc:%20.8E\n",ecomParam[3]);
        fprintf(fout, "Bs:%20.8E\n",ecomParam[4]);
        fprintf(fout, "-- a priori value\n");
        //these two variables used to record the original initial value
        GVector pp = m_initialP;
        GVector vv = m_initialV;
        
        //observation
        GMatrix l(3,1);
        
        //starting calculation
        while(1)
        {
            
            // reset the initial value
            //m_initialP.x += X[0]/1000.0;m_initialP.y += X[1]/1000.0;m_initialP.z += X[2]/1000.0;
            //m_initialV.x += X[3]/1000.0;m_initialV.y += X[4]/1000.0;m_initialV.z += X[5]/1000.0;
            
            for( int i = 0 ; i< np; i++ )
            {
                ecomParam[i] += X[i];
            }
            
            
            if( np > 0)
            {
                m_forceManager.m_forceModels["GFMEMP"]->setModelParameters(np, ecomParam);
            }
            
            
            BTPB.resize(np, np);
            BTPL.resize(np, 1);
            
            PHITPHI.clear();
            PHITL.clear();
            PHITS.clear();
            STPHI.clear();
            STS.clear();
            STL.clear();
            
            
            double LTPL = 0.0;
            
            iteration++;
            
            
            GSpaceEnv::updateSpaceEnvironment(m_initialEpoch);
            m_spaceCraft->getStatePointer()->updateState_eci(m_initialEpoch, m_initialP, m_initialV);
            PHI.setData(phiI, 6, 6); //set phi to be identical matrix
            m_spaceCraft->getStatePointer()->updateTransitionMatrix(PHI);
            m_spaceCraft->getStatePointer()->senMatrix.resize(6, np);
            
            /*
             because there are two types of observations
             the covariance compnents estimation should be applied to get the weight matrix
             */
            
            int iobs = 0;
            average_ele= 0.0;
            phi.clear();
            eta.clear();
            
            acc_eci.clear();
            attitude.clear();
            
            //从初始历元的下一个历元开始
            GTime ct_utc = m_initialEpoch + 900;
            
            GPreciseEphemeris peph;
            //*  2004  2 29 21  0  0.00000000
            for( iobs = 0 ; ct_utc < m_endEpoch ; iobs++,ct_utc = ct_utc + interval )
            {
                
                // if(iobs>=20) break;
                
                //printf("epoch: %d--%s ", iobs,GTime::GTime2CivilTime(ct_utc).TimeString().c_str() );
                
                //GTime ct_utc = m_initialEpoch + 900*(iobs+1);
                
                GTime ct_gps = GTime::UTC2GPST(ct_utc);
                
                //this is the observation
                //m_spaceCraft->getEphValue(index).getPV(po_ecf,vo_ecf); // these are in ECEF, should be transformed to ECI
                //peph = m_spaceCraft->getEphValue(ct_gps);
                
                peph = m_spaceCraft->getEphValue_test(ct_gps);
                
                if( peph.isOK() == false)
                {
                    continue;
                }
                
                peph.getPV(po_ecf,vo_ecf);
                
                // call the propagator, to get the calculation value and state transition matrix
                PropagateTo( ct_utc ); //here ct must be in UTC
                
                //get the calculation value, in km
                pc = m_spaceCraft->getStatePointer()->satpos_eci;
                //vc = m_spaceCraft->getStatePointer()->satvel_eci;
                
                GEarthOrientationParameter myeop(ct_utc);
                myeop.ECEF2ECI_pos(po_ecf, po_eci);
                //myeop.ECEF2ECI_vel(po_ecf, vo_ecf, vo_eci);
                
                
                average_ele += m_spaceCraft->getStatePointer()->attitude_eci.beta;
                
                phi.push_back(m_spaceCraft->getStatePointer()->attitude_eci.phi);
                
                eta.push_back(m_spaceCraft->getStatePointer()->attitude_eci.eta);
                
                attitude.push_back(m_spaceCraft->getStatePointer()->attitude_eci);
                
                // printf("%s, %12.6f %12.6f %12.6f %12.8f %12.8f %12.8f\n",GTime::GTime2CivilTime(ct_utc).TimeString().c_str(),
                //        pc.x, pc.y,pc.z,vc.x,vc.y,vc.z);
                
                //GVector ecom_eci = m_forceManager.m_forceModels["GFMSRP"]->getForce();
                GVector ecom_eci;
                //ecom_eci= m_forceManager.m_forceModels["GFMEMP"]->getForce();
                //acc_eci.push_back(ecom_eci);
                
                
                double obs[3] = { (po_eci.x - pc.x)*1000.0,
                    (po_eci.y - pc.y)*1000.0,
                    (po_eci.z - pc.z)*1000.0
                };
                
                
                //GMatrix l(obs,6,1);
                l.setData(obs, 3, 1);
                
                // be careful the order of multiplication PHi(t2,t0)= PHI(t2,t1)*PHI(t1,t0);
                PHI = m_spaceCraft->getStatePointer()->phiMatrix;  //*PHI;
                S   = m_spaceCraft->getStatePointer()->senMatrix;
                
                for(int i = 0 ; i< 3; i++ )
                {
                    for(int j = 0 ; j< 6; j++)
                    {
                        PHIpos(i,j) = PHI(i,j);
                    }
                    for(int k = 0 ; k< np; k++)
                    {
                        Spos(i,k) = S(i,k);
                    }
                }
                
                //printf("PHI: ");
                //cout<<m_spaceCraft->getStatePointer()->phiMatrix;
                
                
                //printf("l: ");
                //cout<<~l;
                /*
                 //pos+vel version
                 PHITPHI += (~PHI)*PHI;
                 PHITS   += (~PHI)*S;
                 STPHI   += (~S)*PHI;
                 STS     += (~S)*S;
                 PHITL   += (~PHI)*l;
                 STL     += (~S)*l;
                 */
                
                //only pos version
                PHITPHI += (~PHIpos)*PHIpos;
                PHITS   += (~PHIpos)*Spos;
                STPHI   += (~Spos)*PHIpos;
                STS     += (~Spos)*Spos;
                PHITL   += (~PHIpos)*l;
                STL     += (~Spos)*l;
                
                //cout<<"Spos:"<<endl;
                //cout<< Spos;
                
                //calculate LTPL
                for( int j = 0; j< 3; j++ )
                {
                    LTPL += l[j]*l[j];
                }
                
            }
            
            average_ele = average_ele / iobs;
            
            
            
            //cout<<"PHITPHI"<<endl;
            //cout<<PHITPHI<<endl;
            
            //cout<<"PHITS"<<endl;
            //cout<<PHITS<<endl;
            
            //cout<<"STPHI"<<endl;
            //cout<<STPHI<<endl;
            
            cout<<"STS"<<endl;
            cout<<STS<<endl;
            
            
            
            // construct the BTPB and BTPL
            for( int i = 0; i< np; i++ )
            {
                BTPL(i,0) = STL(i,0);
                
                //BTPB
                for( int j = 0 ; j<np; j++)
                {
                    {
                        BTPB(i,j) = STS(i,j);
                    }
                }
            }
            
            GMatrix XX = !BTPB*BTPL;
            
            /*
             cout<<"BTPB"<<endl;
             cout<< BTPB<<endl;
             cout<<"BTPL"<<endl;
             cout<<BTPL<<endl;
             cout<<"X"<<endl;
             cout<<XX<<endl;
             */
            
            //printf("iter:%d\n",iteration);
            //cout << ~XX;
            //cout << ~X;
            
            if( fabs(XX[0]-X[0])<1.0E-9
               &&fabs(XX[1]-X[1])<1.0E-9
               &&fabs(XX[2]-X[2])<1.0E-9
               &&fabs(XX[3]-X[3])<1.0E-9
               &&fabs(XX[4]-X[4])<1.0E-9
               &&fabs(XX[5]-X[5])<1.0E-9
               
               )
            {
                
                //GVector correctionP = (m_initialP - pp)*1000.0; // in meters
                
                //GVector correctionV = (m_initialV - vv)*1000.0; // in meters
                
                //start precision evaluating, VTPV
                double vtv = LTPL;
                vtv = sqrt(vtv/(iobs*3-np));
                
                m_fittingRMS = vtv;
                
                //printf("l: ");
                //cout<<~l<<std::endl;
                
                GMatrix D = (!BTPB);
                
                D = D*vtv*vtv;
                
                //printf("deviation: %.16f,%.16f,%.16f,%.16f,%.16f,%.16f\n",sqrt(D(0,0)),sqrt(D(1,1)),sqrt(D(2,2)),sqrt(D(3,3)),sqrt(D(4,4)),sqrt(D(5,5)));
                GString str_t0 = GTime::GTime2CivilTime(m_initialEpoch).TimeString();
                
                //cout<< D;
                
                //printf("%s,%.6f,%.8f,%.8f,%.8f,%.8f,%.10f,%.10f,%.10f,",str_t0.c_str(), vtv, average_ele,m_initialP.x, m_initialP.y,m_initialP.z,m_initialV.x,m_initialV.y,m_initialV.z);
                
                //printf("%.6E,%.6E,%.6E,%.8E,%.8E,%.8E,%.8E,%.8E,%.8E,%.8E\n", ecomParam[0],ecomParam[1],ecomParam[2],ecomParam[3],ecomParam[4],ecomParam[5],ecomParam[6],ecomParam[7],ecomParam[8],ecomParam[9]);
                
                // printf("%.8E,%.8E,%.8E,%.8E,%.8E\n", ecomParam[0],ecomParam[1],ecomParam[2],ecomParam[3],ecomParam[4]);
                
                /*
                 printf("acc_x, acc_y, acc_z\n");
                 for( int i = 0 ; i< phi.size(); i++ )
                 {
                 double myphi = phi[i];
                 double delta = myphi > 0.0 ?1.0 : -1.0;
                 //double delta = 1.0;
                 
                 GVector acc_xyz;
                 GMatrix T(3,3), temp1(3,1),temp2;
                 T(0,0) = attitude[i].xhat.x; T(0,1) = attitude[i].xhat.y; T(0,2) = attitude[i].xhat.z;
                 T(1,0) = attitude[i].yhat.x; T(1,1) = attitude[i].yhat.y; T(1,2) = attitude[i].yhat.z;
                 T(2,0) = attitude[i].zhat.x; T(2,1) = attitude[i].zhat.y; T(2,2) = attitude[i].zhat.z;
                 temp1(0,0) = acc_eci[i].x/m_spaceCraft->getSpaceCraftGemotry()->m_mass;
                 temp1(1,0) = acc_eci[i].y/m_spaceCraft->getSpaceCraftGemotry()->m_mass;
                 temp1(2,0) = acc_eci[i].z/m_spaceCraft->getSpaceCraftGemotry()->m_mass;
                 
                 temp2 = T*temp1;
                 
                 acc_xyz.x = temp2[0];acc_xyz.y = temp2[1];acc_xyz.z = temp2[2];
                 
                 
                 printf("%.3f,%.8E,%.8E,%.8E\n",myphi,acc_xyz.x,acc_xyz.y,acc_xyz.z);
                 
                 }
                 */
                
                fprintf(fout, "++ estimated\n");
                fprintf(fout, "PX:%20.8f\n",m_initialP.x);
                fprintf(fout, "PY:%20.8f\n",m_initialP.y);
                fprintf(fout, "PZ:%20.8f\n",m_initialP.z);
                fprintf(fout, "VX:%20.8f\n",m_initialV.x);
                fprintf(fout, "VY:%20.8f\n",m_initialV.y);
                fprintf(fout, "VZ:%20.8f\n",m_initialV.z);
                
                fprintf(fout, "D0:%20.8E\n",ecomParam[0]);
                fprintf(fout, "Y0:%20.8E\n",ecomParam[1]);
                fprintf(fout, "B0:%20.8E\n",ecomParam[2]);
                fprintf(fout, "Bc:%20.8E\n",ecomParam[3]);
                fprintf(fout, "Bs:%20.8E\n",ecomParam[4]);
                
                fprintf(fout, "RMS:%20.6f\n",vtv);
                fprintf(fout, "CONVARIANCE:\n");
                for(int i = 0 ; i< np; i++ )
                {
                    for(int j = 0 ; j< np; j++)
                    {
                        fprintf(fout, "%20.12E ", D(i,j));
                    }
                    fprintf(fout, "\n");
                }
                
                fprintf(fout, "-- estimated\n");
                
                
                fprintf(fout, "END of Orbit Fitting Parameters\n");
                
                fclose(fout);
                
                break;
            }
            else
            {
                X = XX;
            }
        }
        
    }
    
    
    
    
    // this is the function with Least Square
    // the main function for calculation
    void GOrbitFitting::calculateInitialError()
    {
        
        std::vector<double> phi; // latitude of sun in bfs
        std::vector<double> eta; // the eta angle for ECOM
        
        std::vector<GSpaceCraftAttitude> attitude;
        std::vector<GVector> acc_eci;
        
        
        
        //ouput to file
        FILE* fout = fopen(m_outputFile.c_str(), "w+");
        
        if(!fout)
        {
           printf("file %s open failed!\n", m_outputFile.c_str());
        }
        
        GSensorID sensorID = m_spaceCraft->getSpaceCraftID();
        fprintf(fout, "Reference Time: %s \t %s\n",GTime::GTime2CivilTime(m_initialEpoch).TimeString().c_str(),m_initialEpoch.getTimeSystem().getTimeSystemName().c_str());
        
        fprintf(fout, "%s\t%3d\n",sensorID.getSystem().c_str(),sensorID.getIDnum());
        
        GVector po_ecf,vo_ecf,po_eci,vo_eci,pc,vc;
        
        double average_ele= 0.0;
        
        double interval = 900;  // interval :900 seconds, this may need to be changed !
        
        int iteration_num = 0;
        
        //m_spaceCraft->getNumofPreEph();
        
        int np = m_forceManager.totalNum_parameters;
        
        GMatrix X(6+np,1); // the solution, correction to the initial position and velocity
        
        GMatrix BTPB(6+np,6+np);
        GMatrix BTPL(6+np,1);
        GMatrix PHI(6,6); // transition matrix from t0 to tn
        
        GMatrix PHIpos(3,6);
        
        //sensitivity matrix
        GMatrix S(6,np);
        GMatrix Spos(3,np);
        
        
        GMatrix W(6,6);  // the weighting matrix
        W(0,0) = 1;W(1,1) = 1;W(2,2) = 1;
        //W(3,3) = 1.0/(100000*100000);W(4,4) = 1.0/(100000*100000);W(5,5) = 1.0/(100000*100000);
        W(3,3) = 1;W(4,4) = 1;W(5,5) = 1;
        
        // temp matrice
        GMatrix PHITPHI(6,6),PHITS(6,np),STPHI(np,6),STS(np,np),PHITL(6,1),STL(np,1);
        
        m_initialP = m_spaceCraft->getStatePointer()->satpos_eci;
        m_initialV = m_spaceCraft->getStatePointer()->satvel_eci;
        int iteration = 0;
        double phiI[36]={ 1,0,0,0,0,0,
                          0,1,0,0,0,0,
                          0,0,1,0,0,0,
                          0,0,0,1,0,0,
                          0,0,0,0,1,0,
                          0,0,0,0,0,1
                        };
        
        //double ecomParam[5]={ -0.00000010042745184630757, -0.00000000020611988712060707, 0.0000000016151040633494485,-0.0000000004195205808724743,0.0000000031361947345823236};
        // they are in units of m/s^2
        double ecomParam[10]={0.0};
        fprintf(fout, "++ a priori value\n");
        fprintf(fout, "PX:%20.8f\n",m_initialP.x);
        fprintf(fout, "PY:%20.8f\n",m_initialP.y);
        fprintf(fout, "PZ:%20.8f\n",m_initialP.z);
        fprintf(fout, "VX:%20.8f\n",m_initialV.x);
        fprintf(fout, "VY:%20.8f\n",m_initialV.y);
        fprintf(fout, "VZ:%20.8f\n",m_initialV.z);
        
        fprintf(fout, "D0:%20.8E\n",ecomParam[0]);
        fprintf(fout, "Y0:%20.8E\n",ecomParam[1]);
        fprintf(fout, "B0:%20.8E\n",ecomParam[2]);
        fprintf(fout, "Bc:%20.8E\n",ecomParam[3]);
        fprintf(fout, "Bs:%20.8E\n",ecomParam[4]);
        fprintf(fout, "-- a priori value\n");
        //these two variables used to record the original initial value
        GVector pp = m_initialP;
        GVector vv = m_initialV;
        
        //observation
        GMatrix l(3,1);
        
        //starting calculation
        while(1)
        {
            // reset the initial value
            m_initialP.x += X[0]/1000.0;m_initialP.y += X[1]/1000.0;m_initialP.z += X[2]/1000.0;
            
            m_initialV.x += X[3]/1000.0;m_initialV.y += X[4]/1000.0;m_initialV.z += X[5]/1000.0;
            
            //mysat.getStatePointer()->updateState_ecef(start_utc, ps, vs);
            
            
            for( int i = 0 ; i< np; i++ )
            {
                ecomParam[i] += X[6+i];
            }
            
            
            if( np > 0)
            {
               m_forceManager.m_forceModels["GFMEMP"]->setModelParameters(np, ecomParam);
            }
            
            
            BTPB.resize(6+np, 6+np);
            BTPL.resize(6+np, 1);
            
            PHITPHI.clear();
            PHITL.clear();
            PHITS.clear();
            STPHI.clear();
            STS.clear();
            STL.clear();
            
            
            double LTPL = 0.0;
            
            iteration++;
            
            
            GSpaceEnv::updateSpaceEnvironment(m_initialEpoch);
            m_spaceCraft->getStatePointer()->updateState_eci(m_initialEpoch, m_initialP, m_initialV);
            PHI.setData(phiI, 6, 6); //set phi to be identical matrix
            m_spaceCraft->getStatePointer()->updateTransitionMatrix(PHI);
            m_spaceCraft->getStatePointer()->senMatrix.resize(6, np);
            
            /*
             because there are two types of observations
             the covariance compnents estimation should be applied to get the weight matrix
             */
            
            int iobs = 0;
            average_ele= 0.0;
            phi.clear();
            eta.clear();
            
            acc_eci.clear();
            attitude.clear();
            
            //从初始历元的下一个历元开始
            GTime ct_utc = m_initialEpoch + 900;
            
            GPreciseEphemeris peph;
            //*  2004  2 29 21  0  0.00000000
            for( iobs = 0 ; ct_utc < m_endEpoch ; iobs++,ct_utc = ct_utc + interval )
            {
                
               // if(iobs>=20) break;
                
                //printf("epoch: %d--%s ", iobs,GTime::GTime2CivilTime(ct_utc).TimeString().c_str() );
                
                //GTime ct_utc = m_initialEpoch + 900*(iobs+1);
                
                GTime ct_gps = GTime::UTC2GPST(ct_utc);
                
                //this is the observation
                //m_spaceCraft->getEphValue(index).getPV(po_ecf,vo_ecf); // these are in ECEF, should be transformed to ECI
                //peph = m_spaceCraft->getEphValue(ct_gps);
                
                peph = m_spaceCraft->getEphValue_test(ct_gps);
                
                if( peph.isOK() == false)
                {
                    continue;
                }
                
                peph.getPV(po_ecf,vo_ecf);
                
                // call the propagator, to get the calculation value and state transition matrix
                PropagateTo( ct_utc ); //here ct must be in UTC
                
                //get the calculation value, in km
                pc = m_spaceCraft->getStatePointer()->satpos_eci;
                //vc = m_spaceCraft->getStatePointer()->satvel_eci;
                
                GEarthOrientationParameter myeop(ct_utc);
                myeop.ECEF2ECI_pos(po_ecf, po_eci);
                //myeop.ECEF2ECI_vel(po_ecf, vo_ecf, vo_eci);
                
                
                average_ele += m_spaceCraft->getStatePointer()->attitude_eci.beta;
                
                phi.push_back(m_spaceCraft->getStatePointer()->attitude_eci.phi);
                
                eta.push_back(m_spaceCraft->getStatePointer()->attitude_eci.eta);
                
                attitude.push_back(m_spaceCraft->getStatePointer()->attitude_eci);
                
               // printf("%s, %12.6f %12.6f %12.6f %12.8f %12.8f %12.8f\n",GTime::GTime2CivilTime(ct_utc).TimeString().c_str(),
               //        pc.x, pc.y,pc.z,vc.x,vc.y,vc.z);
                
                //GVector ecom_eci = m_forceManager.m_forceModels["GFMSRP"]->getForce();
                GVector ecom_eci;
                //ecom_eci= m_forceManager.m_forceModels["GFMEMP"]->getForce();
                //acc_eci.push_back(ecom_eci);
                
                
                double obs[3] = { (po_eci.x - pc.x)*1000.0,
                                  (po_eci.y - pc.y)*1000.0,
                                  (po_eci.z - pc.z)*1000.0
                                };
                
                
                //GMatrix l(obs,6,1);
                l.setData(obs, 3, 1);
                
                // be careful the order of multiplication PHi(t2,t0)= PHI(t2,t1)*PHI(t1,t0);
                PHI = m_spaceCraft->getStatePointer()->phiMatrix;  //*PHI;
                S   = m_spaceCraft->getStatePointer()->senMatrix;
                
                for(int i = 0 ; i< 3; i++ )
                {
                    for(int j = 0 ; j< 6; j++)
                    {
                        PHIpos(i,j) = PHI(i,j);
                    }
                    for(int k = 0 ; k< np; k++)
                    {
                        Spos(i,k) = S(i,k);
                    }
                }
                
                //printf("PHI: ");
                //cout<<m_spaceCraft->getStatePointer()->phiMatrix;
                
                
                //printf("l: ");
                //cout<<~l;
                /*
                //pos+vel version
                PHITPHI += (~PHI)*PHI;
                PHITS   += (~PHI)*S;
                STPHI   += (~S)*PHI;
                STS     += (~S)*S;
                PHITL   += (~PHI)*l;
                STL     += (~S)*l;
                */
                
                //only pos version
                PHITPHI += (~PHIpos)*PHIpos;
                PHITS   += (~PHIpos)*Spos;
                STPHI   += (~Spos)*PHIpos;
                STS     += (~Spos)*Spos;
                PHITL   += (~PHIpos)*l;
                STL     += (~Spos)*l;
                
                
                
                //calculate LTPL
                for( int j = 0; j< 3; j++ )
                {
                    LTPL += l[j]*l[j];
                }
                
            }
            
            average_ele = average_ele / iobs;
            
            
            
             //cout<<"PHITPHI"<<endl;
             //cout<<PHITPHI<<endl;
             
             //cout<<"PHITS"<<endl;
             //cout<<PHITS<<endl;
             
             //cout<<"STPHI"<<endl;
             //cout<<STPHI<<endl;
             
             //cout<<"STS"<<endl;
             //cout<<STS<<endl;
            

            
            // construct the BTPB and BTPL
            for( int i =0; i< 6+np; i++ )
            {
                //BTPL
                if( i<6 )
                {
                    BTPL(i,0) = PHITL(i,0);
                }
                else if( i>=6 )
                {
                    BTPL(i,0) = STL(i-6,0);
                }
                
                //BTPB
                for(int j = 0 ; j<6+np; j++)
                {
                    if( i <6 && j< 6 )
                    {
                        BTPB(i,j) = PHITPHI(i,j);
                    }
                    else if( i<6 && j>=6 )
                    {
                        BTPB(i,j) = PHITS(i,j-6);
                    }
                    else if( i>= 6 && j< 6 )
                    {
                        BTPB(i,j) = STPHI(i-6,j);
                    }
                    else if( i>=6 && j>=6 )
                    {
                        BTPB(i,j) = STS(i-6,j-6);
                    }
                }
            }
            
            
            
            GMatrix XX = !BTPB*BTPL;
            
            /*
            cout<<"BTPB"<<endl;
            cout<< BTPB<<endl;
            cout<<"BTPL"<<endl;
            cout<<BTPL<<endl;
            cout<<"X"<<endl;
            cout<<XX<<endl;
            */
            
            //printf("iter:%d\n",iteration);
            //cout << ~XX;
            //cout << ~X;
            
            if( fabs(XX[0]-X[0])<0.0001
               &&fabs(XX[1]-X[1])<0.0001
               &&fabs(XX[2]-X[2])<0.0001
               &&fabs(XX[3]-X[3])<0.00001
               &&fabs(XX[4]-X[4])<0.00001
               &&fabs(XX[5]-X[5])<0.00001
               
               )
            {
                
                GVector correctionP = (m_initialP - pp)*1000.0; // in meters
                
                GVector correctionV = (m_initialV - vv)*1000.0; // in meters
                
                //start precision evaluating, VTPV
                double vtv = LTPL;
                vtv = sqrt(vtv/(iobs*3-6-np));
                
                m_fittingRMS = vtv;
                
                //printf("l: ");
                //cout<<~l<<std::endl;
                
                GMatrix D = (!BTPB);
                
                D = D*vtv*vtv;
                
                //printf("deviation: %.16f,%.16f,%.16f,%.16f,%.16f,%.16f\n",sqrt(D(0,0)),sqrt(D(1,1)),sqrt(D(2,2)),sqrt(D(3,3)),sqrt(D(4,4)),sqrt(D(5,5)));
                GString str_t0 = GTime::GTime2CivilTime(m_initialEpoch).TimeString();
                
                //cout<< D;
                
                //printf("%s,%.6f,%.8f,%.8f,%.8f,%.8f,%.10f,%.10f,%.10f,",str_t0.c_str(), vtv, average_ele,m_initialP.x, m_initialP.y,m_initialP.z,m_initialV.x,m_initialV.y,m_initialV.z);
                
                //printf("%.6E,%.6E,%.6E,%.8E,%.8E,%.8E,%.8E,%.8E,%.8E,%.8E\n", ecomParam[0],ecomParam[1],ecomParam[2],ecomParam[3],ecomParam[4],ecomParam[5],ecomParam[6],ecomParam[7],ecomParam[8],ecomParam[9]);
                
               // printf("%.8E,%.8E,%.8E,%.8E,%.8E\n", ecomParam[0],ecomParam[1],ecomParam[2],ecomParam[3],ecomParam[4]);
                
                /*
                printf("acc_x, acc_y, acc_z\n");
                for( int i = 0 ; i< phi.size(); i++ )
                {
                    double myphi = phi[i];
                    double delta = myphi > 0.0 ?1.0 : -1.0;
                    //double delta = 1.0;
                    
                    GVector acc_xyz;
                    GMatrix T(3,3), temp1(3,1),temp2;
                    T(0,0) = attitude[i].xhat.x; T(0,1) = attitude[i].xhat.y; T(0,2) = attitude[i].xhat.z;
                    T(1,0) = attitude[i].yhat.x; T(1,1) = attitude[i].yhat.y; T(1,2) = attitude[i].yhat.z;
                    T(2,0) = attitude[i].zhat.x; T(2,1) = attitude[i].zhat.y; T(2,2) = attitude[i].zhat.z;
                    temp1(0,0) = acc_eci[i].x/m_spaceCraft->getSpaceCraftGemotry()->m_mass;
                    temp1(1,0) = acc_eci[i].y/m_spaceCraft->getSpaceCraftGemotry()->m_mass;
                    temp1(2,0) = acc_eci[i].z/m_spaceCraft->getSpaceCraftGemotry()->m_mass;
                    
                    temp2 = T*temp1;
                    
                    acc_xyz.x = temp2[0];acc_xyz.y = temp2[1];acc_xyz.z = temp2[2];
                    
                    
                    printf("%.3f,%.8E,%.8E,%.8E\n",myphi,acc_xyz.x,acc_xyz.y,acc_xyz.z);
                    
                }
                */
                
                fprintf(fout, "++ estimated\n");
                fprintf(fout, "PX:%20.8f\n",m_initialP.x);
                fprintf(fout, "PY:%20.8f\n",m_initialP.y);
                fprintf(fout, "PZ:%20.8f\n",m_initialP.z);
                fprintf(fout, "VX:%20.8f\n",m_initialV.x);
                fprintf(fout, "VY:%20.8f\n",m_initialV.y);
                fprintf(fout, "VZ:%20.8f\n",m_initialV.z);
                
                fprintf(fout, "D0:%20.8E\n",ecomParam[0]);
                fprintf(fout, "Y0:%20.8E\n",ecomParam[1]);
                fprintf(fout, "B0:%20.8E\n",ecomParam[2]);
                fprintf(fout, "Bc:%20.8E\n",ecomParam[3]);
                fprintf(fout, "Bs:%20.8E\n",ecomParam[4]);
                
                fprintf(fout, "RMS:%20.6f\n",vtv);
                fprintf(fout, "CONVARIANCE:\n");
                for(int i = 0 ; i< 6+np; i++ )
                {
                    for(int j = 0 ; j< 6+np; j++)
                    {
                        fprintf(fout, "%20.12E ", D(i,j));
                    }
                    fprintf(fout, "\n");
                }
                
                fprintf(fout, "-- estimated\n");
                
                
                fprintf(fout, "END of Orbit Fitting Parameters\n");
                
                fclose(fout);
                
                
                
                
                
                
                break;
            }
            else
            {
                X = XX;
            }
            
            
            
        }
        
        
        
    }
    
    // get the derivatives including the first derivatives of position , velocity and state transition
    // get dy/dt , dPHI/dt and dS/dt
    void GOrbitFitting::getDerivatives( int n, double x, double *y, double *dydx )
    {
        
        //update p and v from y
        GVector p, v;
        p.x = y[0]; p.y = y[1]; p.z = y[2];
        v.x = y[3]; v.y = y[4]; v.z = y[5];
        
        //need to reset the time each step for the calculation of gravity force
        GTime ct = m_t0 + x;  // for each step, should be in UTC
        //        //need to update the space environment and the motion state at the same time
        GSpaceEnv::updateSpaceEnvironment(ct);  //update the
        //need to update the spacecraft state in every step
        m_spaceCraft->getStatePointer()->updateState_eci(ct,p,v);
        
        //printf("%s\n", GTime::GTime2CivilTime(ct).TimeString().c_str());
        
        //double phi[6*6] = {0.0};
        
        //GMatrix sM = m_spaceCraft->getStatePointer()->senMatrix;
        //GMatrix phiM = m_spaceCraft->getStatePointer()->phiMatrix ;
        
        int np = m_spaceCraft->getStatePointer()->senMatrix.getColNO();
        
        for( int i = 0 ; i< 6; i++ )
        {
            for( int j = 0; j< 6; j++ )
            {
                m_spaceCraft->getStatePointer()->phiMatrix(i,j) = y[6+i*6+j];
            }
            for(int k = 0 ; k<np; k++ )
            {
                m_spaceCraft->getStatePointer()->senMatrix(i,k) = y[42+i*np+k];
            }
        }
        
        //m_spaceCraft->getStatePointer()->updateTransitionMatrix(phiM);
        //m_spaceCraft->getStatePointer()->updateSensitivityMatrix(sM);
        
        
        int m = 18 + 6 + 3*np; // 18 da/df, 6 dr/dt + 3*np da/dp
        
        double *dydx_tmp = new double[m];
        memset(dydx_tmp,0,sizeof(double)*m);
        
        /* the structure of dydx_tmp matrix
         y = (r,v), f=(v,a)
         
       dydt  | u         v       w        ax       ay       az     |
      dax/dr | dax/dx,  dax/dy,  dax/dz,  dax/du,  dax/dv,  dax/dw |
      day/dr | day/dx,  day/dy,  day/dz,  day/du,  day/dv,  day/dw |
      daz/dr | daz/dx,  daz/dy,  daz/dz,  daz/du,  daz/dv,  daz/dw |
      dax/dp | dax/dp1, dax/dp2, dax/dp3, dax/dp4, dax/dp...       |
      day/dp | day/dp1, day/dp2, day/dp3, day/dp4, day/dp...       |
      daz/dp | daz/dp1, daz/dp2, daz/dp3, daz/dp4, daz/dp...       |
         
        */
        
        m_forceManager.getDerivatives(ct, m, x, y, dydx_tmp, m_spaceCraft);
        
        /* Note that, dPHI/dt = dF/dy * PHI
         PHI = dy/dy0 ; PHI(t0) = I
         df/dy is like this :
            |0        I    |
            |              |
            |da/dr   da/dv |
         
         assume da/dv = 0
         
         if considering the sensitivity matrix S:
          S = dy/dp ; S(t0) = 0
         
                                |  0    |
         dS/dt = df/dy * S +    |       |
                                | da/dp |, 6 by np
         
         dPHI/dt = df/dy*PHI
         
                                         |0   0     |
         d(PHI,S)/dt = dF/dy * (PHI,S) + |          |
                                         |0   da/dp |, 6 by (6+np)
         
         */
        
        
        GMatrix dfdy(6,6);
        GMatrix Stmp(6,np);  //
        //set the identity part
        for( int i = 0 ; i< 3; i++ )
        {
            dfdy(i,3+i) = 1.0;
        }
        
        for( int i = 0 ; i< 3; i++ )
        {
            //set da/dr and da/dv part
            for( int j =0; j< 6; j++ )
            {
                dfdy( 3+ i, j) = dydx_tmp[6 + i*6 + j];
            }
            
            //set Stmp, da/dp part
            for( int k = 0 ; k< np; k++ )
            {
                Stmp(3+i,k) = dydx_tmp[24+i*np+k];
            }
        }
        
        // get the dPHIdt matrix
        GMatrix dPHIdt = dfdy * m_spaceCraft->getStatePointer()->phiMatrix;
        
        //get the dS/dt matrix
        GMatrix dSdt = dfdy*m_spaceCraft->getStatePointer()->senMatrix + Stmp;
        
        
       // printf("dfdy:\n");
       // cout<<dfdy<<endl;
//        
//        
       // printf("sM:\n");
       // cout<<sM<<endl;
//        
       // printf("Stmp:\n");
       // cout<<Stmp<<endl;
//        
//        printf("dSdt:\n");
//        cout<<dSdt<<endl;
//        
        // set up dydx, include dy/dt, dPHI/dt and dS/dt
        for( int i = 0 ; i< 6; i++ )
        {
            dydx[i] = dydx_tmp[i];
        }
        
        for( int i = 0 ; i< 6; i++ )
        {
            for( int j = 0 ; j< 6; j++ )
            {
                dydx[6+i*6+j] = dPHIdt(i,j);
            }
            
            for( int k = 0 ; k< np; k++ )
            {
                dydx[42+ i*np +k] = dSdt(i,k);
            }
        }
        
        if(dydx_tmp != NULL)
        {
            delete[] dydx_tmp;
            dydx_tmp = NULL;
        }
        
    }
    
    void GOrbitFitting::PropagateTo(gfc::GTime t)
    {
        m_t0 = m_spaceCraft->getStatePointer()->m_epoch;
        
        GVector p = m_spaceCraft->getStatePointer()->satpos_eci;
        GVector v = m_spaceCraft->getStatePointer()->satvel_eci;
        
        gfc::TimeSystem ts;
        
        long mjd = 0,sod = 0;
        
        double fsod = 0.0;
        
        (t - m_t0).GetData( ts, mjd, sod, fsod );
        
        double sec_end = mjd*GCONST("SECPDAY") + sod + fsod;
        
        //set up the initial value for integration
        //GMatrix myphi = m_spaceCraft->getStatePointer()->phiMatrix;
        
        //std::cout<< myphi;
        
        //this is just for testing
        int np = m_spaceCraft->getStatePointer()->senMatrix.getColNO();
        
        int ndim = 6 + 6*6 + 6*np;  // state vector and the transition matrix
        
        // the structure of y0: 0-5 state , 6-41 state transition matrix, 42-- the sensitiviy matrix( 6 by np )
        double *y0 = new double[ndim];  // the inital value for integration
        double *yend = new double[ndim];  // the state after integration
        memset(y0,0,sizeof(double)*ndim);
        memset(yend,0,sizeof(double)*ndim);
        
        y0[0] = p.x; y0[1] = p.y; y0[2] = p.z; y0[3] = v.x; y0[4] = v.y; y0[5] = v.z;
        
        // PHI matrix is an identity matrix
        // y0[6] = 1.0; y0[13] = 1.0; y0[20] = 1.0; y0[27] = 1.0; y0[34]=1.0; y0[41] = 1.0;
        
        for( int i = 0 ; i< 6; i++ )
        {
            for(int j = 0 ; j< 6; j++)
            {
                y0[6+i*6+j] = m_spaceCraft->getStatePointer()->phiMatrix(i,j);
            }
        }
        
        
        for(int i = 0 ; i< 6; i++ )
        {
           for( int j =0; j< np; j++ )
            {
                y0[42+i*np+j] = m_spaceCraft->getStatePointer()->senMatrix(i,j);
            }
        }
        
        //printf("%f %f %f %f %f %f\n",p.x, p.y, p.z, v.x, v.y, v.z);
        
        GTime testT(57094,55784,0.0,"tsUTC");
        
//        if(t == testT)
//        {
//            int tesc = 0;
//        }
        
        m_integrator->IntegrateTo( this , ndim, 0.0, y0, sec_end, yend);
        
        p.x = yend[0];p.y = yend[1];p.z = yend[2];
        v.x = yend[3];v.y = yend[4];v.z = yend[5];
        
        //printf("%f %f %f %f %f %f\n",p.x, p.y, p.z, v.x, v.y, v.z);
        
        double phi[36]= {0.0};
        double *ss = new double[6*np];
        for( int i = 0 ; i< 36; i++ )
        {
            phi[i] = yend[6+i];
        }
        
        for(int i = 0 ;i< 6*np; i++)
        {
            ss[i] = yend[42+i];
        }
        
        
        GSpaceEnv::updateSpaceEnvironment(t);
        
        m_spaceCraft->getStatePointer()->updateState_eci(t, p, v);
        m_spaceCraft->getStatePointer()->updateTransitionMatrix(phi);
        m_spaceCraft->getStatePointer()->updateSensitivityMatrix(ss);
        
        
        
        //printf("Sen matrix:\n");
        //std::cout<< m_spaceCraft->getStatePointer()->senMatrix<<endl;
        
        
        if( ss != NULL )
        {
            delete[] ss;
            ss = NULL;
        }

        
        if( y0 != NULL )
        {
            delete[] y0;
            y0 = NULL;
        }
        
        if( yend != NULL )
        {
            delete[] yend;
            yend = NULL;
        }
        
        int testc = 0;
        
    }
    
    
} // end of namespace



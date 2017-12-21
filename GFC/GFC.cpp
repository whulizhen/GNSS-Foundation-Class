// GFC.cpp : ∂®“Âøÿ÷∆Ã®”¶”√≥Ã–Úµƒ»Îø⁄µ„°£
//

#include <fstream>
#include <thread>

#include "../../src/GFCCONST.h"
#include "../../src/LogStream.h"
#include "../../src/GOBSData.h"
#include "../../src/GTime.h"
#include "../../src/EllipsoidMgr.h"
#include "../../src/GMatrix.h"
#include "../../src/GSensor.h"
#include "../../src/GNSSData.h"

#include "../../src/GMath.hpp"
#include "../../src/JPLEPH.h"

#include "../../src/GNetcdf.h"

#include <pthread.h>

using namespace std;
using namespace gfc;

void* threadFun(void* param)
{
    while(1)
    {
        cout<<"A"<<endl;
        sleep(1);
    }
    return (void*)(NULL);
}


class A
{
public:
    A() {m_a ="m_a"; p_a="p_a";}
    ~A(){}
    GString m_a;
    static void INIT() {s_a ="s_a";}
    static void setS() {s_a ="setS";}
    static GString getS(){return s_a;}
    GString getM() {return m_a;}
    GString getP() {return p_a;}
private:
    static GString s_a;
    GString p_a;
};


GString A::s_a = "s_a1";

class B : public A
{
public:
    B(){m_b ="B_m_a";}
    void changeM(){m_a = "changed";}
    //static GString getS() {return B::s_a;}
    GString m_b;  //覆盖了父类中的成员变量
    
};


//void HouseholderUpperTri( double* M,int num_of_row, int num_of_column  );
//void CholeskyUUT( double* P,int n, double* U);
//void InverseUpperTriangular(double* R, int n, double* U);

int main( int argc, char* argv[] )
{
    
    GNetcdf  mync("CERES_SYN1deg-Month_Terra-Aqua-MODIS_Ed3A_Subset_200003-201504.nc", GNetcdf::NCMODE::MODE_NOWRITE);
    
    double mydatatest[10];
    mync.getData(mydatatest, "toa_sw_all_mon");
    
    
    
    
    JPLEPH jpleph;
    jpleph.ReadHeader_a("header.405");//("jpleph.405");
    jpleph.ReadRecord_a("ascp2020.405");
    jpleph.EphTest("testpo.405");
    
    
    
    
    //A::INIT();
    A testa;
    B testb;
    GString testS1 = testb.m_a;
    testb.changeM();
    A::setS();
    GString testS2 = testb.m_a;
    GString testS3 = testa.m_a;
    GString  testS4 = B::getS();
    
//    pthread_t  mythread;
//    pthread_mutex_t mymutex;
//    int abcdef[3]={0};
//    
//    pthread_create(&mythread,NULL,threadFun,(void*)abcdef);
    //pthread_join(mythread, NULL);
    //GOBSType gobsdt = GOBSType::GetByName("ot_PHASE");
    
    gnssOBSDataType gdt("ssGPS","ot_RANGE","CA","cf_L1");
    gnssOBSData mydata(gdt,22232324.345,3.5,char(0),char(0));
    //mydata.setValue(gdt,434343.343,23);
    double testvalue = mydata.getDataValue();
    
    //gnssOBSDataType::RegByName("ot_ABC");  //静态函数也一块继承了
    
    //gdt.dump(std::cout);
    //gnssDATA::datatypelist.size();
    //gnssDATA::RegByName();
    
    double t2 = GMath::MAX(5, 3);
    
    //gnssDataType  gdt;
    //GOBSType   odt;
    
    int sizetest = sizeof(GTime);
    
    GCarrierFreq  cf1;
    cf1.setFreq("cf_L1");
    
    GCarrierFreq* p = &cf1;
    long pdata = (long)(&cf1);
    
    GCarrierFreq myp = *(GCarrierFreq*)(pdata);
    
    TimeSystem  ts;
    TimeSystem::dump(std::cout);
    //TimeSystem::RegTimeSystem("tsABC");
    TimeSystem::UnregByName("tsGPS");
    TimeSystem::dump(std::cout);
    
    
	double matrixa[] =
	{
		3,4,-3,-1,-2,2,-1,-1,1
	};
	
	double matrixb[] =
	{
		4.02, -1.56, 9.81,
		6.19,  4.00, -4.09,
		-8.22, -8.67, -4.57,
		-7.57,  1.75, -8.61,
		-3.03,  2.86, 8.99
	};
	
	GMatrix matrixA( matrixa,3,3 );
	GMatrix matrixB( matrixb,5,3 );
	//double detA = matrixA.det();
	cout<<matrixB;
	
	(!matrixA).dump();
	
	//matrixA.testLapack();
	matrixB = matrixA;
    GString   doublestring(12987);
    GString  hexstring = doublestring.d2x();
    
    GString  mystring("     hellohellohello_my_cpp_world");
    GString  string1 = mystring.upperCase();
    GString string2  = string1.lowerCase();
    GString astring("  ");
    GString inputstring ="hello";
    GString outputstring ="HELLO";
    GString teststring = mystring.rightJustify(5);
    mystring.stripLeading_v(' ');
    
	GTime  mytime1;
	GTime  mytime2;
	GTime  mytime3,mytime4, mytime5,mytime6;
	CivilTime ct, ct1;
	
    ct.m_year = 2015; ct.m_month = 4; ct.m_day = 4; ct.m_hour = 23; ct.m_minute = 59; ct.m_second = 59.879876546353;ct.m_ts = GTimeSystem("tsGPS");
	JDTime jdt = GTime::CivilTime2JDTime(ct);
	mytime2 = GTime::JDTime2GTime(jdt);
	JDTime jdt1 = GTime::GTime2JDTime(mytime2);
	ct1 = GTime::JDTime2CivilTime(jdt);		
	mytime1.SetFromCivilTime( ct );	
	
	mytime3 = mytime1 - mytime2;
	mytime5.SetData( GTimeSystem("tsGPS"), 0,0,34434534 );
	mytime4 = mytime1 + mytime5;
	
	mytime6 = mytime4;
	 
	mytime6 = mytime6 + 234.45;
		
	double leapsec = GTime::getLeapSecond( 1980, 5, 1 );
//
	GFCCONST::dump( std::cout );
	double clight = GCONST("CLIGHT");
	long double nanosec = GCONST("CLIGHT")/GCONST("NANO");
	GFCCONST::RegByName("ABC",3);
	double abc = GCONST("ABC");
	GFCCONST::UnregByName("ABC");
	
	
	EllipsoidMgr ellipmgr;
	ellipmgr.dump(std::cout);
	double a = ellipmgr.A_km("cgcs2000");
	
	std::ofstream oflog("test.log",std::ios_base::out);
	ConfigureLOG::Stream() = &oflog; // else use &(std::cout) - the default
	
	ConfigureLOG::ReportLevels() = true;
	ConfigureLOG::ReportTimeTags() = true;
	
	ConfigureLOG::ReportingLevel() = ConfigureLOG::Level("DEBUG");
	
	LOG(INFO)<<"Log Tester"<<endl;
	
	ConfigureLOG::Stream() = &(std::cout); // else use &(std::cout) - the default
	
	LOG(INFO) << "In main(), level is " << ConfigureLOG::ReportingLevel();
	LOG(INFO) << "This is an INFO message in main";
	LOG(WARNING) << "This is a WARNING message in main";
	LOG(ERROR) << "This is an ERROR message in main";
//
		

	return 0;
}


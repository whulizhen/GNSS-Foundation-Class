#include <iostream>
using namespace std;

extern "C"
{
    void erpfboxw_(int* ERM, char* myIPATH, int* LEN_DATPATH,int* ANT, int* GRD, int* REFF,double* YSAT,double* SUN,double* KAPPA,int* MONTH, int* BLKNUM,double* ACCEL);
    void ftest_(int* ii, float* ff, char* DATPATH);
    void myerptest_(int* ERM,int* ANT, int* GRD, int* REFF,double* YSAT,double* SUN,double* KAPPA,int* MONTH, int* BLKNUM,double* ACCEL);
}

int main()
{

	//fortran test
    int ERM = 2;
    int ANT = 0;
    int GRD = 1;
    int REFF = 0;
    int MONTH = 3;
    int BLKNUM = 4;
    double YSAT[6] = { -17662416.430000001,4031322.1279999998,-19465284.702000001,
                        412.19989000000001,-3679.8486600000002,-1122.6964100000001
                     };
    
    double SUN[3] = {139782958138.39048,-45244942785.067253,-19615287917.492874};
    
    double KAPPA[3][3] = {
                            {-0.4052188296198202,0.91421971355857101,0.00012434800857546492},
                            {-0.91421964219783969,-0.40521884732965979,0.00036275125925288264},
                            {0.00038202250902986989,0.000033312248809693642,0.99999992647454571},
                         };


    double accel[3] ={0.0};
    double acc[3] = {0.0};
    //double* acc = new double[3];

    char* DATPATH = "//Users//lizhen//myProject//GFC//GFC_Proj//Earth_Radiation_Pressure//DATA//";
    int len = strlen(DATPATH);
    //erpfboxw_(&ERM, DATPATH, &len, &ANT, &GRD, &REFF, YSAT,SUN, KAPPA[0],&MONTH,&BLKNUM,ACCEL);
    printf("1\n");
    myerptest_(&ERM, &ANT, &GRD, &REFF, YSAT,SUN, KAPPA[0],&MONTH,&BLKNUM,acc);
    

    printf("2\n");

    myerptest_(&ERM, &ANT, &GRD, &REFF, YSAT,SUN, KAPPA[0],&MONTH,&BLKNUM,accel);

   // printf("%e %e %e\n", acc[0], acc[1], acc[2]);

    char testStr[100] = "abcdefg";
    int ii = 5;
    float ff = 5.5;
    ftest_(&ii,&ff,testStr);
    ftest_(&ii,&ff,testStr);
    ftest_(&ii,&ff,testStr);

}


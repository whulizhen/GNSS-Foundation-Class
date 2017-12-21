/*
 some functions for Square Information Root Filter(SRIF)
 reference: factorization method for discrete sequential estimation(Bieman)
 */


#include <math.h>
#include <memory>
#include <iostream>
#include <fstream>
using namespace std;
/*对矩阵A进行三角化（上三角）Householder 正交变换
 reference：Appendix IV.A
 matrix A: m rows, n columns
 |s|  |
 T*A =| |AA|
 |0|  |
 
 T is the Househlder Transformation matrix
 There is no need to know the matrix T
 
 Author: lizhen
 date;2015 11 25
 */
void HouseholderUpperTri( double* M,int num_of_row, int num_of_column )
{
    int colcount = 0;
    int m = num_of_row;
    int n = num_of_column;
    double* A = M ; // use a new pointer, but the memory keep original
    //int* test = num_of_row<=num_of_column?(&m):(&n);  // let test point to the minimum of m and n
    
    int test = num_of_row<=num_of_column?(m-1):(n-1);
    double *u = new double[num_of_row];
    
    while( A != NULL )
    {
        if( colcount >= test+1 )
        {
            break;
        }
        
        double s =0.0 ,s1 =0.0,b =0.0,r = 0.0 ;
        int t =0 , i =0, j = 0 ;
        if( *A >= 0 ) t = 1;
        else if( *A < 0 ) t =-1;
        for( i = 0; i<m;i++ )
        {
            s1+=  (*(A+i*num_of_column ))*(*(A+i*num_of_column ));
        }
        s = -t*sqrt(s1); // or s = -1*sqrt(s1) 保证 s的符号与*A相反，这样可以保证u[0]尽量大
        memset(u,0,sizeof(double)*num_of_row);
        for( i =0 ; i<m; i++ )
        {
            if( i == 0 )
            {
                u[i] = *A - s;
            }
            else
            {
                u[i] = *(A+i*num_of_column);
            }
        }
        
        b = s*u[0];  // variable b can be very small, if b = 1.0/(s*u[0])
        
        for( i = 0; i<m ; i++ )
        {
            *(A+i*num_of_column) = 0.0;
        }
        *A = s;
        
        for( j =1 ; j< n ; j++ )  //the jth column
        {
            s1 = 0.0;
            for( i = 0 ; i< m ; i++ )  // all the rows
            {
                s1+= u[i]*A[i*num_of_column+j];
            }
            r = s1/b;
            for( i = 0 ; i< m ; i++ )  // all the rows
            {
                *(A+i*num_of_column+j) += (r*u[i]);
            }
        }
        
        
        //        printf("iterate count: %d\n",colcount);
        //        for( i  = 0 ; i< num_of_row; i++ )
        //        {
        //            for( j =0; j<num_of_column  ;j++ )
        //            {
        //                printf("%8.4f ",M[i*num_of_column+j]);
        //            }
        //            printf("\n");
        //        }
        
        colcount++;
        m = num_of_row - colcount;
        n = num_of_column - colcount;
        A = M + num_of_column*colcount + colcount;
    }
    
    if( u != NULL )  { delete[] u; u = NULL;}
}


/*
 *
 * back substitution algorithm to solve the Rx=y , where R is upper triangular
 * Author: lizhen
 * date: 2015 12 3
 * reference: Appendix IV.B
 * input: R , R= |R|y|, R should be the matrix including y, for the convinience of SRIF
 * output: x
 */

void backSubstitution(double* R, double* y, double* x, int n )
{
    double total =0.0;
    int j =0, k =0;
    for( j = n-1 ; j>= 0; j-- )
    {
        total = 0.0;
        for(  k =j+1 ; k<n; k++ )
        {
            total+= ( R[j*n+k]*x[k] );
        }
        
        x[j] = (y[j] -total)/R[j*n+j];
    }
}


/* Cholesky Factorization, with P=UUT
 P是实对称矩阵
 U是上三角矩阵
 */
void CholeskyUUT( double* P,int n, double* U)
{
    int i = 0, j = 0, k =0 ;
    for(j =n;j>=2; j-- )
    {
        U[(j-1)*n+j-1] = sqrt(P[(j-1)*n+j-1]);
        for(k = 1; k<=j-1;k++)
        {
            U[(k-1)*n+j-1] = P[(k-1)*n+j-1]/U[(j-1)*n+j-1];
        }
        
        for(k = 1; k<=j-1;k++ )
        {
            for(i=1;i<=k;i++)
            {
                P[(i-1)*n+k-1] = P[(i-1)*n+k-1]-U[(i-1)*n+j-1]*U[(k-1)*n+j-1];
            }
        }
    }
    
    U[0] = sqrt(P[0]);
}


/*Inversion of the upper triangular matrix
 input:
 R-->Upper triangular matrix to be inversed
 n : the column number of the matrix R
 output: U--> the inversion of R matrix
 
 reference: P65
 */

void InverseUpperTriangular(double* R, int n, double* U)
{
    int i =0, j =0, k =0;
    double sum =0.0;
    U[0] = 1.0/R[0];
    for(j=2;j<=n;j++)
    {
        U[(j-1)*n+j-1] = 1.0/R[(j-1)*n+j-1];
        for(k=1;k<=j-1;k++)
        {
            sum =0.0;
            for(i=k; i<=j-1;i++)
            {
                sum = sum + U[(k-1)*n+i-1]*R[(i-1)*n+j-1];
            }
            U[(k-1)*n+j-1] = -sum*U[(j-1)*n+j-1];
        }
    }
}


/* SRIF without process noise, the paraments should NOT be time varying
 *input:
 * A: the current design matrix
 * IM: the current information matrix
 * obsnum: the current obs num
 * xnum:   the current unknown parameter number
 *output: updatd IM
 */
void SRIF0( double* IM,double* obsvalue,int obsnum, int xnum,double* A)
{
    int i = 0, j =0;
    //measurement update, set the Information Matrix
    for( i =0; i< obsnum; i++ )
    {
        for( j=0;j<xnum+1;j++ )
        {
            if( j>=xnum ) // obs value part
            {
                IM[(xnum+i)*(xnum+1)+j] = obsvalue[i];
            }
            else  // design matrix part
            {
                IM[(xnum+i)*(xnum+1)+j] = A[i*xnum + j];
            }
        }
    }
    
    //household transformation
    HouseholderUpperTri(IM, xnum+obsnum, xnum+obsnum);
    
    
}

/*
 * Author: zhen.LI
 * date: 2015 12 25, Merry Christmas!!
 * 该函数是包含有状态转移过程以及状态方程过程噪声的滤波算法
 * reference: P121, equation 2.28 and  2.29
 * equation 2.28 is measuremant update, measruement model: z=Ax+v;
 * equation 2.29 is the time update, state model: x1 = PHI*X0 + G*w0;
 * comments: G: nrow = xnum, ncol =wnum; w0: nrow = wnum, ncol=1;
 * inputs:
 *       wnum           xnum          1
 *      | Rw           Rwx           zw | wnum
 * IM = |                               |
 *      | 0          R(j+1)     z^(j+1) | xnum
 * z^ = R*x0; z^,R are a priori measurement value and covariance(P0 = R^-1*R^-T)
 * IM should be initialized before the estimation process; nrow = xnum+wnum; ncol = xnum+wnum+1
 *
 * A: the design matrix for measurement equation; nrow= obsnum, ncol = xnum
 * z: the current obs value, nrow = obsnum, ncol = 1;
 * PHI: the inversion of the state transformation matrxi(IT MUST BE UNSINGULAR); nrow = xnum, ncol = xnum;
 * G :  the coefficient of process noise ; nrow = xnum, ncol= wnum;
 */
void SRIF( int xnum, int obsnum, int wnum,double* IM,double* A,double* z,
          double* PHI,double* G)
{
    int i = 0, j = 0 , k = 0 ;
    int nrow_DM = xnum+obsnum;
    int ncol_DM = xnum+1;
    double* DM = new double[nrow_DM*ncol_DM];
    memset(DM,0,sizeof(double)*nrow_DM*ncol_DM);
    double sum =0.0, tmp =0.0;
    /*      xnum  1
     *      | R   z^| xnum
     * DM = |       |
     *      | A   z | obsnum
     *
     */
    for( i = 0 ; i< nrow_DM; i++ )
    {
        for( j = 0 ; j< ncol_DM; j++ )
        {
            if( i< xnum && j< xnum )  // R part, data equation
            {
                DM[i*ncol_DM+j] = IM[(i+wnum)*(xnum+wnum+1)+j+wnum];
            }
            else if( i<xnum && j>=xnum ) // z^ part, data equation
            {
                DM[i*ncol_DM+j] = IM[(i+wnum)*(xnum+wnum+1)+j+wnum];
            }
            else if( i>= xnum && j<xnum ) // A part , measurement equation
            {
                DM[i*ncol_DM+j] = A[(i-xnum)*xnum+j];
            }
            else if( i>= xnum && j>=xnum ) // z part, measurement equation
            {
                DM[i*ncol_DM+j] = z[i-xnum];
            }
        }
    }// end of the formation of matrix DM
    
    printf("measurement updata(before) DM matrix:\n");
    for( i = 0 ; i< nrow_DM; i++ )
    {
        for( j = 0 ; j< ncol_DM; j++ )
        {
            
            printf("%6.4f ",DM[i*ncol_DM+j]);
        }
        printf("\n");
    }
    
    /*
     *     |R z^|  |R^ z^^|
     *  Tj*|    |= |      |
     *     |A z |  |0  e  |
     *
     *  the process of Householder transformation
     *  this is also called the measurement updation
     */
    
    HouseholderUpperTri(DM,nrow_DM,ncol_DM);
    
    printf("measurement updata(after) DM matrix:\n");
    for( i = 0 ; i< nrow_DM; i++ )
    {
        for( j = 0 ; j< ncol_DM; j++ )
        {
            
            printf("%6.4f ",DM[i*ncol_DM+j]);
        }
        printf("\n");
    }
    
    /*
     * obtaining the R and z,which are updated by the measurement equation,called measuement updation.
     * then, combining with the process noise w and state transformation matrxi PHI
     * we can doing the state transformation, which is called time updation as well.
     *
     *
     *
     *        | Rw              0            zw  |   | Rw(j+1)    Rwx(j+1)    zw(j+1) |
     * (Tj+1)*|                                  | = |                                |
     *        |-R^*PHI^(-1)*G   R^*PHI^(-1)  z^^ |   | 0          R(j+1)      z^(j+1) |
     *
     * the first row , Rw(j+1) , Rwx(j+1) and zw(j+1) are used in smoothing problem
     * the second row, R(j+1) and z^(j+1) are used in filter problme
     *
     */
    
    /*computing R^*PHI^-1*/
    
    double* RPHI = new double[xnum*xnum];
    memset(RPHI,0,sizeof(double)*xnum*xnum);
    
    for( i = 0 ; i< xnum ; ++i )
    {
        for( k=0;k<xnum;++k )
        {
            sum = DM[i*(xnum+1)+k];
            for( j =0 ;j< xnum; ++j )
            {
                RPHI[i*xnum+j]+= sum*PHI[k*xnum+j];
            }
        }
    }
    
    //computing the R^*PHI^-1*G and apply it into the IM matrix
    for( i = 0 ; i<xnum; i++ )
    {
        for( k =0; k<xnum; k++ )
        {
            tmp = RPHI[i*xnum+k];
            sum = 0.0;
            for( j =0; j<wnum; j++ )
            {
                // 这里在赋值之前需要将IM相应部分清0
                sum+= ( -tmp*G[k*wnum+j] );
            }
            
            IM[(wnum+k)*(wnum+xnum+1)+j]  = sum;
        }
    }
    
    //apply the RPHI^-1 to the IM matrix
    for( i = 0 ; i< xnum ; i++ )
    {
        for( j=0;j<xnum; j++ )
        {
            IM[(wnum+i)*(wnum+xnum+1) + wnum+j] = RPHI[i*xnum+j];
        }
        
        //udpate z^^ from Matrix DM
        IM[(wnum+i)*(wnum+xnum+1) + wnum+xnum] = DM[i*(xnum+1) + xnum];
    }
    
    
    printf("measurement updata(before) IM matrix:\n");
    for( i = 0 ; i< wnum+xnum; i++ )
    {
        for( j = 0 ; j< wnum+xnum+1; j++ )
        {
            
            printf("%6.4f ",IM[i*(wnum+xnum+1)+j]);
        }
        printf("\n");
    }
    
    //Time updating,
    HouseholderUpperTri(IM,wnum+xnum,wnum+xnum+1);
    
    
    printf("measurement updata(after) IM matrix:\n");
    for( i = 0 ; i< wnum+xnum; i++ )
    {
        for( j = 0 ; j< wnum+xnum+1; j++ )
        {
            
            printf("%6.4f ",IM[i*(wnum+xnum+1)+j]);
        }
        printf("\n");
    }
    
    
    
    if( DM != NULL )
    {
        delete[] DM; DM = NULL;
    }
    
    if( RPHI != NULL )
    {
        delete[] RPHI, RPHI = NULL;
    }
    
}


/* generate random number with normal distribution ---------------------------*/
double randn(double myu, double sig)
{
    double PI = 3.14159265357;
    double a,b;
    a=((double)rand()+1.0)/((double)RAND_MAX+1.0);  /* 0<a<=1 */
    b=((double)rand()+1.0)/((double)RAND_MAX+1.0);  /* 0<b<=1 */
    return myu+sqrt(-2.0*log(a))*sin(2.0*PI*b)*sig;
}



/*
 simulate the observables which is time varying
 this kind of paraments can be estimated  withe SRIF() function
 this is a test  main function of the SRIF
 */
void testSRIF1()
{
    double PI = 3.14159265357;
    double error = 0.0;
    double a = 3.2 , b = 1.2, c = 3.0;
    int N = 1000,N0, i =0 , j =0, k = 0;
    int xnum =3;
    int wnum =1;
    int obsnum = 1;
    double w0 = 0.5; // process noise
    double pw0 = 10.0; // the initial variance of w0
    double Rw0 = sqrt(pw0);
    
    double G[3] ={0.0,0.0,0.0}; // the coefficient in front of the processing noise w0
    double IM[4*5] ={0.0}; // the information matrix
    double PHI_[9] ={0.0}; // the state transformation matrix
    //double x0[3] =   {sin(1.0/PI/10.0), cos(1.0/PI/10.0),0.00001};  // initial value , 3 rows and 1 col
    double x0[3] =   {3.2, 1.2, 3.0};  // initial value , 3 rows and 1 col
    
    double D0[9]   = {100,0,0,0,100,0,0,0,100}; // initial cov-variation matrix
    double invR[9] = {0.0};
    double R[9]    = {0.0};
    double y_[3]    = {0.0};
    double yw =  0.0;
    double A[3]={0.0};
    
    double *x = new double[N];
    double *y = new double[N];
    FILE*   myfile = fopen("srif.dat","w+");
    N0 = 0;
    double alfa = 0.005;
    double beita = 10.0;
    for( i = 0 ; i< N ; i++ )
    {
        x[i] = (i-N0)*alfa;
        //a = sin(i/PI/beita);
        a = 3.2;
        //b = cos(i/PI/beita);
        b = 1.2;
        //c = sin(x[i]);
        c = 3.0;
        y[i] = a*x[i]*x[i] + b*x[i] + c + randn(0,0.5);
        
        //fprintf(myfile, "%d %8.4f %8.4f %8.4f %8.4f %8.4f\n",i,a,b,c,x[i],y[i]);
        
    }
    
    //start the processing code
    
    CholeskyUUT(D0, xnum, invR);
    Rw0 = sqrt(pw0);
    InverseUpperTriangular(invR, xnum, R);
    Rw0 = 1.0/Rw0;
    y_[0] = R[0]*x0[0] + R[1]*x0[1] + R[2]*x0[2];
    y_[1] = R[3]*x0[0] + R[4]*x0[1] + R[5]*x0[2];
    y_[2] = R[6]*x0[0] + R[7]*x0[1] + R[8]*x0[2];
    yw = Rw0*w0;
    
    //set up the information matrix, 6 parts
    for( i = 0;  i < wnum+xnum ; i++ ) // row
    {
        for( j=0;j<wnum+xnum+1;j++ )  //col
        {
            if( i<wnum && j<wnum )  // the Rw part
            {
                IM[i*(wnum+xnum+1)+j] = Rw0;
            }
            else if( i<wnum && j>=wnum+xnum) // the zw part
            {
                IM[i*(wnum+xnum+1)+j] = yw;
            }
            else if( i>=wnum && j>=wnum && j<wnum+xnum ) // the Rx part
            {
                IM[i*(wnum+xnum+1)+j] = R[(i-wnum)*(xnum)+j-wnum];
            }
            else if( i>=wnum && j>=wnum+xnum ) // the zx part
            {
                IM[i*(wnum+xnum+1)+j] = y_[i-wnum];
            }
        }
    }
    
    printf("Initial IM:\n");
    for( i = 0 ; i<wnum+xnum ; i++)
    {
        for( j =0; j<wnum+xnum+1;j++ )
        {
            printf("%6.3f ",IM[i*(wnum+xnum+1)+j]);
        }
        printf("\n");
    }
    
    
    // start the filtering
    for( i = 1 ; i< N ; i++ )
    {
        A[0] = x[i]*x[i];
        A[1] = x[i];
        A[2] = 1.0;
        
        double deltaX = 0.0;
        if( i!= 0 )
        {
            deltaX = x[i] - x[i-1];
        }
        
        // reset the PHI matrxi , the state transformation matrix
        //PHI_[0] = 1.0/( 1.0+1.0/alfa/beita/PI/tan((1.0/alfa*x[i]+N0)/beita/PI)*deltaX  );
        //PHI_[4] = 1.0/( 1.0-1.0/alfa/beita/PI*tan((1.0/alfa*x[i]+N0)/beita/PI)*deltaX  );
        //PHI_[8] = 1.0/( 1.0+1.0/tan(x[i])*deltaX );
        
        PHI_[0] = 1.0;  PHI_[4]=1.0; PHI_[8] = 1.0;
        
        printf("state transformation Matrix PHI:\n");
        for( j = 0 ; j< xnum ; j++ )
        {
            for( k = 0 ; k< xnum ; k++ )
            {
                printf("%8.6f ",PHI_[j*xnum+k]);
            }
            printf("\n");
        }
        
        SRIF(xnum, obsnum, wnum, IM, A, y+i, PHI_, G);
        
        // get the solution of abc
        
        double RR[9] ={0.0};
        
        for( j = 0 ; j< xnum ; j++ )
        {
            for( k =0 ; k< xnum ; k++ )
            {
                RR[j*xnum+k] = IM[(j+wnum)*(wnum+xnum+1)+k+wnum];
            }
            y_[j] = IM[(wnum+j)*(wnum+xnum+1)+wnum+xnum];
        }
        
        double x_est[3]={0.0};
        backSubstitution(RR, y_, x_est, xnum);
        
        //double aa = sin(i/PI/beita);
        double aa = 3.2;
        //double bb = cos(i/PI/beita);
        double bb = 1.2;
        double cc = 3.0;//sin(x[i]);
        double testobs = A[0]*x_est[0] + A[1]*x_est[1] + A[2]*x_est[2];
        double trueValue = A[0]*aa + A[1]*bb + A[2]*cc;
        printf("X_est:%8.4f %8.4f %8.4f\n", x_est[0],x_est[1],x_est[2]);
        fprintf(myfile, "%8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f\n", x_est[0],x_est[1],x_est[2],aa,bb,cc,y[i],testobs,trueValue);  //tan((200*x[i]+N0))
        
    }
    
    
    if(myfile != NULL ) { fclose(myfile); myfile = NULL;}
    
    
}




/*
 * test the SRIF Algorithm with simulated data
 * the simulation scenery is the y = ax^2 + bx + c
 */
void testSRIF()
{
    
    double PI = 3.14159265357;
    double error = 0.0;
    double a = 3.2 , b = 1.2, c = 3.0;
    int N = 1000, i =0 , j =0, k = 0;
    double *x = new double[N];
    double *y = new double[N];
    double randvalue = 0.0;  // the random value for simulation
    double testobs =0.0;
    double trueValue = 0.0;
    FILE*   myfile = fopen("srif.dat","w+");
    
    for( i = 0 ; i< N ; i++ )
    {
        x[i] = (i-N/2)*0.005;
        //a = sin(i/PI/10);
        //b = cos(i/PI/10);
        
        //        if( i == N/2)
        //        {
        //            error = 5;
        //        }
        //        else if( fabs(i-N/3)<1 )
        //        {
        //            error = 5;
        //        }
        //        else if( fabs(i-2*N/3)<1 )
        //        {
        //            error = 5 ;
        //        }
        
        y[i] = a*x[i]*x[i] + b*x[i] + c + randn(0,0.2) + error;
        error = 0;
    }
    
    // start the data processing with srif
    // the estimated parameters are a, b and c
    // set the initial information matrix for initial value of abc
    double x0[3] =   {3.2, 1.2,3.0};  // initial value , 3 rows and 1 col
    double D0[9]   = {100,0,0,0,100,0,0,0,100}; // initial cov-variation matrix
    double R[9]    = {0.0};
    double invR[9] = {0.0};
    double A[3]    = {0.0};  // coefficients for the measurement equation
    double IM[16]   = {0.0};  //, 4 rows and 4 colums. the information matrix for the data equation, including the measurement equation and state equation
    double y_[3]    = {0.0};
    double Ay[4]    = {0.0};
    // derive the R0 and Z0 for Z0 = R0*X0 + v; this is called the data equation
    // while R0 = inv(sqrt(D0));
    //first cholesky factorization of the matrix D0, we can get R0^-1
    CholeskyUUT(D0, 3, invR);
    InverseUpperTriangular(invR, 3, R);
    y_[0] = R[0]*x0[0] + R[1]*x0[1] + R[2]*x0[2];
    y_[1] = R[3]*x0[0] + R[4]*x0[1] + R[5]*x0[2];
    y_[2] = R[6]*x0[0] + R[7]*x0[1] + R[8]*x0[2];
    
    // the initial information matrix
    IM[0] = R[0]; IM[1] = R[1]; IM[2] = R[2]; IM[3] = y_[0];
    IM[4] = R[3]; IM[5] = R[4]; IM[6] = R[5]; IM[7] = y_[1];
    IM[8] = R[6]; IM[9] = R[7]; IM[10] = R[8]; IM[11] = y_[2];
    
    //start the estimating process
    for( i = 0 ; i< N ; i++ )
    {
        A[0] = x[i]*x[i];
        A[1] = x[i];
        A[2] = 1.0;
        //Ay[0] = A[0]; Ay[1] = A[1];Ay[2]=A[2];Ay[3]=y[i];
        
        // this is called the measurement update
        //IM[12] = A[0]; IM[13] = A[1]; IM[14]=A[2];IM[15]=y[i];  //complete the Information matrix
        //HouseholderUpperTri(IM, 4, 4);
        SRIF0(IM,y+i,1,3,A);
        //After the HouseholderUpperTri, We can get the TimeUpdate
        // that is to say, we have to update the state equation part of the Information Matrix. the R and y_
        
        //get the solution
        double RR[9] ={0.0};
        RR[0] = IM[0];RR[1] = IM[1];RR[2] = IM[2];
        RR[3] = IM[4];RR[4] = IM[5];RR[5] = IM[6];
        RR[6] = IM[8];RR[7] = IM[9];RR[8] = IM[10];
        y_[0] = IM[3];y_[1] = IM[7];y_[2] = IM[11];
        
        double x_est[3]={0.0};
        backSubstitution(RR, y_, x_est, 3);
        
        testobs = A[0]*x_est[0] + A[1]*x_est[1] + A[2]*x_est[2];
        trueValue = A[0]*a + A[1]*b + A[2]*c;
        
        for(j =0 ; j< 4; j++ )
        {
            for(k = 0; k< 4; k++ )
            {
                printf("%8.6f ",IM[j*4+k]);
            }
            printf("\n");
        }
        
        printf("x_estimation: %.3f %.3f %.3f %.3f\n",x_est[0],x_est[1],x_est[2],IM[15]);
        
        fprintf(myfile, "%8.4f  %8.4f  %8.4f %8.4f %8.4f %8.4f %8.6f\n",trueValue,y[i],testobs,x_est[0],x_est[1],x_est[2],IM[15]);
        
    }
    
    if(myfile != NULL )
    {
        fclose(myfile);
        myfile = NULL;
    }
    
    if(x != NULL)  {delete[] x; x = NULL;}
    if(y != NULL)  {delete[] y; y = NULL;}
    
}








/*
 *
 * Measurement Process Update
 * Author: lizhen
 * date: 2015 12 3
 * reference: Appendix VII.B, page 156-157
 * input: Np, Nx,Ny,m, S(L,Ntot), where L = Np+Nx+m, Ntot = Np + Nx + Ny +1;
 *         S is the a priori information array
 * output: S is the posteriori informaiton array, Lower part of S contains measurement
 */
void  MeasureUpdate(double* S, int Np, int Nx, int Ny, int m)
{
    int j =0, i =0,k =0, L=0, Ntot= 0;
    double delta =0.0,alfa =0.0;
    
    double* v = new double[Np +Nx];
    
    L = Np + Nx + m;
    Ntot = Np + Nx + Ny + 1;
    for(j =0; j< Np+Nx; j++ )
    {
        delta =0.0;
        for( i = j;i<L; i++ )
        {
            v[i] = S[i*Ntot + j];
            S[i*Ntot + j] = 0.0;
            delta = delta + v[i]*v[i];
        }
        if(delta <= 0.0) continue;
        
        delta = sqrt(delta);
        
        if( v[j] > 0.0 ) delta = -delta;
        S[j*Ntot+j] = delta;
        v[j] = v[j] - delta;
        delta  = 1.0/(delta*v[j]);
        for( k = j+1; k< Ntot; k++ )
        {
            alfa = 0.0;
            for( i =j; i< L ; i++ )
            {
                alfa = alfa + S[i*Ntot+k]*v[i];
            }
            alfa = alfa * delta;
            for( i =j; i< L ; i++ )
            {
                S[i*Ntot+k] = S[i*Ntot+k] + alfa*v[i];
            }
        }
    }
    if(v != NULL ) {delete[] v; v = NULL;}
}

/*
 *
 * Measurement Process Update , compute the posteriori of y information array
 * Author: lizhen
 * date: 2015 12 3
 * reference: Appendix VII.B, page 156-157
 * input:  Ny, m, Ry[Ny(Ny+3)/2], Sbar ,where Ny is the number of bias paramenters, m is the number of measurements being processed,
 Ry is the triangular, vector stored SRIF a priori y information array, NRY = Ny(Ny+3)/2
 Sbar is the comprised of the bottom m rows and the last Ny+1 columns of S, this is the y portion of the observation
 *         S is the a priori information array
 * output: S is the posteriori informaiton array, Lower part of S contains measurement
 */
void PostBiasPar(int Ny, int m,double* Ry,double* Sbar )
{
    int kk =0 ,k=0, j =0, i =0, jj =0, L =0;
    double delta =0.0 , alfa =0.0, beta =0.0 ;
    for( j =0; j< Ny; j++ )
    {
        kk = kk + j;
        delta = 0.0;
        for( i =0; i< m; i++ )
        {
            delta = delta + Sbar[i*(Ny+1)+j ]*Sbar[i*(Ny+1)+j ];
        }
        if(delta  == 0.0 ) continue;
        
        delta = sqrt(delta + Ry[kk])*sqrt(delta + Ry[kk]);
        if( Ry[kk] > 0.0 ) delta = -delta;
        alfa  = Ry[kk] - delta;
        Ry[kk] = delta;
        beta = 1.0/(alfa*delta);
        jj = kk;
        L = j;
        for( k=j+1;k<Ny+1;k++ )
        {
            
        }
    }
    
}

/*
 *
 * Time update, SRIF Propagation from Tto T+DT
 * Author: lizhen
 * date: 2015 12 3
 * reference: Appendix VII.B, page 156-157
 * input:  dt : Propagation interval
 Np : Number of colored noise variables
 Nx : Number of dynamic variables
 Ny : Number of bias parameters
 Nd : Number of dynamically related colored noise variables
 tau[Np] : Time constants corresponding to colored noise variables
 Vp[Nx][Nd] : The first Nd columns of the Vp matrix correspond to the dynamic paraments; the last Np-Nd columns are omitted because the are in theory zero
 Rw[Np]: Process noise standard deviation reciprocals
 S[Npx+ Np][Ntot]:  The top Npx rows of S contain the SRIF array corresponding to the p and x variables, the bottom p rows are used to store smoothing-related terms
 where Ntot = Np + Nx + Ny + 1; Npx = Np + Nx; DM[j] = exp(-dt/tau[j]), if tau[j] = 0.0 then DM[j] =0.0 ;
 *
 * output: S[Npx+Np][Ntot]: Time-updated array with smoothing-related terms stored in the bottom portion of S
 SIG[Np] : Smoothing-related coefficients
 */
void TimeUpdate( double dt, int Np, int Nx, int Ny, int Nd, double* tau, double* Vp,double* Rw, double* S ,double* SIG )
{
    int j1 =0,j2 =0, i =0, k =0, L = 0 ;
    double alfa =0.0 , delta =0.0, kappa = 0.0;
    int Ntot = Np + Nx + Ny + 1;
    double* DM = new double[Np];  // remember to delete this pointer ;
    double* v  = new double[Np+Nx];
    double* Zw = new double[Np];  // Colored noise is generally 0 mean; zw , in most cases, is not explicitly included;
    memset(Zw,0,sizeof(double)*Np);
    
    memset(v,0,sizeof(double)*(Np+Nx));
    
    for( i = 0 ; i< Np ; i++ )
    {
        if( fabs(tau[i]) <= 1.0E-8 )
        {
            DM[i] = 0.0;
        }
        else
        {
            DM[i] = exp(-dt/tau[i]);
        }
    }
    
    for(j1 =0; j1< Np ; j1++ )
    {
        if( j1<= Nd )
        {
            for( i =0; i< Np + Nx ; i++ )
            {
                for( k =0; k< Nx; k++ )
                {
                    S[i*Ntot + 0 ] = S[i*Ntot + 0 ] - S[i*Ntot + Np+k]*Vp[k*Nd+j1];
                }
            }
        }  // end of if( j1<= Nd )
        
        alfa = -Rw[j1]*DM[j1];
        alfa = alfa*alfa;
        
        delta = 0.0;
        for( i = 0; i< Np + Nx ; i++ )
        {
            v[i] = S[i*Ntot+0];
            delta = delta + v[i]*v[i];
        }
        delta = sqrt(delta);
        alfa = alfa - delta;
        
        SIG[j1] = delta;
        delta = 1.0/(delta*alfa);
        
        for(j2 =1 ;j2<Ntot; j2++)
        {
            kappa = 0.0;
            // if(j2 == Ntot-1) kappa = alfa*Zw[j1];
            for( i=0; i< Np + Nx ; i++ )
            {
                kappa = kappa + S[i*Ntot+j2]*v[i];
            }
            
            kappa = kappa*delta;
            L = j2 -1;
            if(j2 + 1 > Np ) L = j2;
            S[(Np+Nx+j1)*Ntot+L] = kappa*alfa;
            
            for( i= 0 ; i< Np + Nx ; i++ )
            {
                S[i*Ntot+L] = S[i*Ntot+j2] + kappa*v[i];
            }
            
        } // end of for(j2 =1 ;j2<Ntot; j2++)
        
        // S[(Np+Nx+j1)*Ntot+Ntot-1] = S[(Np+Nx+j1)*Ntot+Ntot-1] + kappa*Zw[j1];
        
        kappa = alfa*Rw[0]*delta;
        S[(Np+Nx+j1)*Ntot+Np-1] = Rw[j1] + kappa*alfa;
        for( i=0 ; i< Np+Nx; i++ )
        {
            S[i*Ntot+Np-1] = kappa*v[i];
        }
    } // end of j1
    
    if(v != NULL )  {delete[] v; v = NULL;}
    if(Zw != NULL ) {delete[] Zw; Zw = NULL;}
    if(DM != NULL ) {delete[] DM; DM = NULL;}
    
    
}




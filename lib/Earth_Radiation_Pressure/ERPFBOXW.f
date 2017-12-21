C*
      SUBROUTINE ERPFBOXW(ERM,DATPATH,ANT,GRD,REFF,YSAT,SUN,KAPPA,MONTH,
     1                    BLKNUM,FORCE)
CC
CC NAME       :  ERPFBOXW
CC
CC PURPOSE    :  COMPUTATION OF EARTH RADIATION PRESSURE ACTING ON A
CC               BOW-WING SATELLITE
CC
CC PARAMETERS :
CC         IN :  ERM       : EARTH RADIATION MODEL
CC                           0 = NONE
CC                           1 = EARTH RADIATION PRESSURE (ANALYTICAL)
CC                           2 = EARTH RADIATION PRESSURE (CERES DATA)
CC               CERES DAT : CERES DATA (EXTERNAL ASCII FILES) -> SET PATH IN DATPATH
CC               ANT       : 0 = NO ANTENNA THRUST
CC                         : 1 = WITH ANTENNA THRUST
CC               GRD       : 1 = 2.5 degree grid
CC                         : 2 = 5.0 degree grid
CC                         : 3 = 10  degree grid
CC               REFF      : ACCELERATION IN REFERENCE FRAME:
CC                           0 = INERTIAL
CC                           1 = BODY FIXED (Z,Y,X)
CC                           2 = SUN FIXED (D,Y,B)
CC                           3 = ORBITAL (RADIAL, ALONG- AND CROSS-TRACK)
CC               YSAT      : SATELLITE POSITION [m] AND VELOCITY [m/s] (INERTIAL), 
CC                           (POSITION,VELOCITY) = (RX,RY,RZ,VX,VY,VZ)
CC               SUN       : SUN POSITION VECTOR [m] (INERTIAL)
CC               KAPPA     : ROTATION MATRIX FROM INERTIAL TO EARTH-FIXED
CC                           REFERENCE FRAME (FOR ERM=2)
CC               MONTH     : MONTH NUMBER (1 ... 12), FOR ERM=2
CC               BLKNUM    : BLOCK NUMBER
CC                             1 = GPS-I
CC                             2 = GPS-II
CC                             3 = GPS-IIA
CC                             4 = GPS-IIR
CC                             5 = GPS-IIR-A
CC                             6 = GPS-IIR-B
CC                             7 = GPS-IIR-M
CC                             8 = GPS-IIF
CC                           101 = GLONASS
CC                           102 = GLONASS-M
CC                           103 = GLONASS-K (not yet available)
CC
CC        OUT : ACCEL      : ACCELERATION VECTOR [m/s^2]
CC
CC AUTHOR     : C.J. RODRIGUEZ-SOLANO
CC              rodriguez@bv.tum.de
CC
CC VERSION    : 1.0 (OCT 2010)
CC
CC CREATED    : 2010/10/18
CC
C*
C
      IMPLICIT NONE
C
      INTEGER*4 BLKNUM,ERM,ANT,GRD,REFF,SBLK,INDB,MONTH
      INTEGER*4 IFIRST,II,JJ,K,LFNLOC,IOSTAT,LENPATH
      INTEGER*4 LATK,LONK,GRDCER
      INTEGER*4 ,save::LATKMX,LONKMX
C
      REAL*8 PI,C,AU,S0,TOA,ALB,MASS
      REAL*8 YSAT(6),SUN(3),FORCE(3),ACCEL(3)
C
C     EARTH AND SATELLITE PROPERTIES,Don't forget the SAVE attribution
      REAL*8 ,save::CERES_R(72,144),CERES_E(72,144)
      REAL*8 ,save::CERGRE(72,144),CERGEM(72,144)
      REAL*8 ,save::D_AREA_ALL(72,144),V_NS_ALL(72,144,3)
      REAL*8 ,save::AREA(4,2,15),REFL(4,2,15),DIFU(4,2,15),ABSP(4,2,15)
      REAL*8 ,save::AREA2(4,2,15),REFL2(4,2,15)
      REAL*8 ,save::DIFU2(4,2,15),ABSP2(4,2,15)
      REAL*8 ,save::REFLIR(4,2,15),DIFUIR(4,2,15),ABSPIR(4,2,15)
      REAL*8 ,save::AREAS(4,2),REFLS(4,2),DIFUS(4,2),ABSPS(4,2)
      REAL*8 ,save::AREA2S(4,2),REFL2S(4,2),DIFU2S(4,2),ABSP2S(4,2)
      REAL*8 ,save::REFLIRS(4,2),DIFUIRS(4,2),ABSPIRS(4,2)
C
C     ATTITUDE
      REAL*8 RADVEC(3),ALGVEC(3),CRSVEC(3),ESUN(3)
      REAL*8 RSUN,ABSPOS,ABSCRS,ABSY0V,Y0SAT(3)
      REAL*8 Z_SAT(3),D_SUN(3),Y_SAT(3),B_SUN(3),X_SAT(3)
      REAL*8 ATTSURF(3,4)
C
      REAL*8 ABSNCFVI,ABSNCFIR,ALBFAC,PHASEVI,PHASEIR
      REAL*8 NCFVEC(3)
      REAL*8 D_ANG,GRDANG,GRDNUM,LATIND,LONIND
      REAL*8 D_AREA,LAT_IN,LON_IN,DIST2,ABSDIST
      REAL*8 COSLAT,PSI
      REAL*8 V_NS(3),V_INS(3),V_DIST(3),V_SAT(3)
      REAL*8 COS_IN,COS_RE,PSIDOT,ABSSUN
      REAL*8 REFL_CF,EMIT_CF,E_REFL,E_EMIT
      REAL*8 FORCE_LL(3),FREF(3)
      REAL*8 ANTFORCE
      REAL*8 KAPPA(3,3)
C
      CHARACTER*2 F_MONTH(12)
      CHARACTER*100 DATPATH,FILREFL,FILEMIT
      CHARACTER*5000 LINE
      CHARACTER*25 ITEM
C
      DATA IFIRST/1/
      DATA (F_MONTH(K),K=1,12)/
     . '01', '02', '03', '04', '05', '06', 
     . '07', '08', '09', '10', '11', '12'/

C DATA FILES PATH AND NAME (CHANGE TO LOCAL PATH)
C     DATPATH = 'your/local/CERES/data/path'

C COMPLETE FILE NAMES
      II = LEN(DATPATH)
      DO WHILE (DATPATH(II:II).EQ.' ')
         II = II-1
      ENDDO
      LENPATH = II

      FILREFL = DATPATH(1:LENPATH) // 'REFLMO' // F_MONTH(MONTH)
      FILEMIT = DATPATH(1:LENPATH) // 'EMITMO' // F_MONTH(MONTH)

C Constants needed
      PI = 4D0*DATAN(1D0)
      C  = 299792458D0
      AU = 149597870691D0
C Solar Constant(W/m2)
      S0 = 1367D0
C Top of Atmosphere for CERES
      TOA = 6371000D0 + 30000D0
C Albedo of the Earth
      ALB = 0.3D0
C Initialization of force vector
      FORCE(1) = 0D0
      FORCE(2) = 0D0
      FORCE(3) = 0D0

      ACCEL(1) = 0D0
      ACCEL(2) = 0D0
      ACCEL(3) = 0D0

      LFNLOC = 999

C      write(*,*) "finished the definitions"
C      write(*,*) "IFIRST,ERM", IFIRST,ERM
C ----------------------------------------
C LOAD SATELLITE PROPERTIES AND CERES DATA
C ----------------------------------------
      IF ((IFIRST==1).AND.(ERM.GT.0)) THEN

C     PROPERTIES FOR ALL SATELLITES BLOCKS
         DO SBLK = 1,105
            IF((SBLK.LE.10).OR.(SBLK.GT.100))THEN 

            CALL PROPBOXW(SBLK,AREAS,REFLS,DIFUS,ABSPS,AREA2S,REFL2S,
     1                    DIFU2S,ABSP2S,REFLIRS,DIFUIRS,ABSPIRS)
            
            IF(SBLK.LE.10)THEN
               INDB = SBLK
            ELSEIF(SBLK.GT.100)THEN
               INDB = SBLK - 90
            ENDIF

            DO II = 1,4
               DO JJ = 1,2
                  AREA(II,JJ,INDB) = AREAS(II,JJ)
                  REFL(II,JJ,INDB) = REFLS(II,JJ)
                  DIFU(II,JJ,INDB) = DIFUS(II,JJ)
                  ABSP(II,JJ,INDB) = ABSPS(II,JJ)

                  AREA2(II,JJ,INDB) = AREA2S(II,JJ)
                  REFL2(II,JJ,INDB) = REFL2S(II,JJ)
                  DIFU2(II,JJ,INDB) = DIFU2S(II,JJ)
                  ABSP2(II,JJ,INDB) = ABSP2S(II,JJ)

                  REFLIR(II,JJ,INDB) = REFLIRS(II,JJ)
                  DIFUIR(II,JJ,INDB) = DIFUIRS(II,JJ)
                  ABSPIR(II,JJ,INDB) = ABSPIRS(II,JJ)
               ENDDO
            ENDDO
            ENDIF
         ENDDO

         IF(ERM.EQ.2)THEN

C     REFLECTIVITY 
            CERES_R = 0
            OPEN(UNIT=LFNLOC,FILE=FILREFL,STATUS='UNKNOWN',
     1           FORM='FORMATTED',IOSTAT=IOSTAT)
            DO II=1,72
               READ(LFNLOC,"(A)")LINE
               DO JJ=1,144
                  ITEM = LINE((JJ-1)*25+1:JJ*25)
                  IF (INDEX(ITEM,'NaN') /= 0) THEN
                     CERES_R(II,JJ) = 0D0
                  ELSE
                     READ(ITEM,*) CERES_R(II,JJ)
                  ENDIF
               ENDDO
            ENDDO
            CLOSE(LFNLOC)
C     EMISSIVITY 
            CERES_E = 0
            OPEN(UNIT=LFNLOC,FILE=FILEMIT,STATUS='UNKNOWN',
     1           FORM='FORMATTED',IOSTAT=IOSTAT)
            DO II=1,72
               READ(LFNLOC,"(A)")LINE
               DO JJ=1,144
                  ITEM = LINE((JJ-1)*25+1:JJ*25)
                  IF (INDEX(ITEM,'NaN') /= 0) THEN
                     CERES_E(II,JJ) = 0D0
                  ELSE
                     READ(ITEM,*) CERES_E(II,JJ)
                  ENDIF
               ENDDO
            ENDDO
            CLOSE(LFNLOC)

C     PRE-INTEGRATION, INDEPENDENT OF SATELLITE POSITION
            GRDANG = 2.5D0
            LATIND = 36.5D0
            LONIND = 72.5D0
            IF(GRD.EQ.1)THEN
               GRDNUM = 1D0
               LATKMX = 72
               LONKMX = 144
            ELSEIF(GRD.EQ.2)THEN
               GRDNUM = 2D0
               GRDCER = 2
               LATKMX = 36
               LONKMX = 72
            ELSEIF(GRD.EQ.3)THEN
               GRDNUM = 4D0
               GRDCER = 4
               LATKMX = 18
               LONKMX = 36
            ENDIF
 
            GRDANG = GRDANG*GRDNUM
            LATIND = (LATIND-0.5D0)/GRDNUM + 0.5D0
            LONIND = (LONIND-0.5D0)/GRDNUM + 0.5D0
              
            D_ANG = (PI*GRDANG/180D0)**2

            DO LATK = 1,LATKMX
               DO LONK = 1,LONKMX

                  LAT_IN = (LATK-LATIND)*GRDANG*(PI/180D0)
                  LON_IN = (LONK-LONIND)*GRDANG*(PI/180D0)

C                 Sphere normal vector and differential of area
                  COSLAT = DCOS(LAT_IN)
                  D_AREA_ALL(LATK,LONK) = (TOA**2)*COSLAT*D_ANG
                  V_NS_ALL(LATK,LONK,1) = COSLAT*DCOS(LON_IN)
                  V_NS_ALL(LATK,LONK,2) = COSLAT*DSIN(LON_IN)
                  V_NS_ALL(LATK,LONK,3) = DSIN(LAT_IN)

C                 New matrix of Reflectivity and Emissivity
                  CERGRE(LATK,LONK) = 0D0
                  CERGEM(LATK,LONK) = 0D0
                  IF(GRD.EQ.1)THEN
                     CERGRE(LATK,LONK) = CERES_R(LATK,LONK)
                     CERGEM(LATK,LONK) = CERES_E(LATK,LONK)
                  ELSEIF((GRD.EQ.2).OR.(GRD.EQ.3))THEN
                     DO II = 0,(GRDCER-1)
                        DO JJ = 0,(GRDCER-1)
                           CERGRE(LATK,LONK) = CERGRE(LATK,LONK) 
     1                     + CERES_R(GRDCER*LATK-II,GRDCER*LONK-JJ)
                           CERGEM(LATK,LONK) = CERGEM(LATK,LONK)
     1                     + CERES_E(GRDCER*LATK-II,GRDCER*LONK-JJ)
                        ENDDO
                     ENDDO
                     CERGRE(LATK,LONK) = CERGRE(LATK,LONK)/(GRDNUM**2)
                     CERGEM(LATK,LONK) = CERGEM(LATK,LONK)/(GRDNUM**2)
                  ENDIF

               ENDDO
            ENDDO
C
          ENDIF
         IFIRST = 0
C         write(*,*) "load properties, this should be loaded one time"
      ENDIF
                  
C      write(*,*) "TEST: AREA" , AREA
C       write(*,*) "TEST: REFL" , REFL
C        write(*,*) "TEST: DIFU" , DIFU
C      write(*,*) "TEST: ABSPS" , ABSPS

C      write(*,*) "TEST: AREA2S" , AREA2S
C      write(*,*) "TEST: REFL2S" , REFL2S
C      write(*,*) "TEST: DIFU2S" , DIFU2S
C      write(*,*) "TEST: ABSP2S" , ABSP2S

C      write(*,*) "TEST: REFLIRS" , REFLIRS
C      write(*,*) "TEST: DIFUIRS" , DIFUIRS
C      write(*,*) "TEST: ABSPIRS" , ABSPIRS
      

C --------------------------
C NOMINAL SATELLITE ATTITUDE
C --------------------------

      ABSPOS = DSQRT(YSAT(1)**2+YSAT(2)**2+YSAT(3)**2)
      DO K=1,3
         RADVEC(K) = YSAT(K)/ABSPOS
      ENDDO
      
      CRSVEC(1) = YSAT(2)*YSAT(6)-YSAT(3)*YSAT(5)
      CRSVEC(2) = YSAT(3)*YSAT(4)-YSAT(1)*YSAT(6)
      CRSVEC(3) = YSAT(1)*YSAT(5)-YSAT(2)*YSAT(4)
      ABSCRS = DSQRT(CRSVEC(1)**2+CRSVEC(2)**2+CRSVEC(3)**2)
      DO K=1,3
         CRSVEC(K) = CRSVEC(K)/ABSCRS
      ENDDO

      ALGVEC(1) = CRSVEC(2)*RADVEC(3)-CRSVEC(3)*RADVEC(2)
      ALGVEC(2) = CRSVEC(3)*RADVEC(1)-CRSVEC(1)*RADVEC(3)
      ALGVEC(3) = CRSVEC(1)*RADVEC(2)-CRSVEC(2)*RADVEC(1) 

      IF(ERM.GT.0)THEN

C        DISTANCE FROM SATELLITE TO SUN
         RSUN = DSQRT((YSAT(1)-SUN(1))**2+(YSAT(2)-SUN(2))**2+
     1           (YSAT(3)-SUN(3))**2)

C        D VECTOR AND Z VECTOR
         DO K=1,3
            D_SUN(K) = (SUN(K)-YSAT(K))/RSUN
            Z_SAT(K) = -YSAT(K)/ABSPOS
         ENDDO

C        Y VECTOR
         Y0SAT(1) = Z_SAT(2)*D_SUN(3)-Z_SAT(3)*D_SUN(2)
         Y0SAT(2) = Z_SAT(3)*D_SUN(1)-Z_SAT(1)*D_SUN(3)
         Y0SAT(3) = Z_SAT(1)*D_SUN(2)-Z_SAT(2)*D_SUN(1)
         ABSY0V = DSQRT(Y0SAT(1)**2 + Y0SAT(2)**2 + Y0SAT(3)**2)
         DO K=1,3
            Y_SAT(K) = Y0SAT(K)/ABSY0V
         ENDDO

C        B VECTOR
         B_SUN(1) = Y_SAT(2)*D_SUN(3) - Y_SAT(3)*D_SUN(2)
         B_SUN(2) = Y_SAT(3)*D_SUN(1) - Y_SAT(1)*D_SUN(3)
         B_SUN(3) = Y_SAT(1)*D_SUN(2) - Y_SAT(2)*D_SUN(1)

C        X VECTOR
         X_SAT(1) = Y_SAT(2)*Z_SAT(3) - Y_SAT(3)*Z_SAT(2)
         X_SAT(2) = Y_SAT(3)*Z_SAT(1) - Y_SAT(1)*Z_SAT(3)
         X_SAT(3) = Y_SAT(1)*Z_SAT(2) - Y_SAT(2)*Z_SAT(1)

         DO K=1,3
            ATTSURF(K,1) = Z_SAT(K)
            ATTSURF(K,2) = Y_SAT(K)
            ATTSURF(K,3) = X_SAT(K)
            ATTSURF(K,4) = D_SUN(K)
         ENDDO
      ENDIF
C      write(*,*) "nomial attitude"
C ---------------------------- 
C OPTICAL PROPERTIES PER BLOCK
C ----------------------------

      IF(ERM.GT.0)THEN
         IF(BLKNUM.LE.10)THEN
            INDB = BLKNUM
         ELSEIF(BLKNUM.GT.100)THEN
            INDB = BLKNUM - 90
         ENDIF

         DO II = 1,4
            DO JJ = 1,2
               AREAS(II,JJ) = AREA(II,JJ,INDB)
               REFLS(II,JJ) = REFL(II,JJ,INDB)
               DIFUS(II,JJ) = DIFU(II,JJ,INDB)
               ABSPS(II,JJ) = ABSP(II,JJ,INDB)

               AREA2S(II,JJ) = AREA2(II,JJ,INDB)
               REFL2S(II,JJ) = REFL2(II,JJ,INDB)
               DIFU2S(II,JJ) = DIFU2(II,JJ,INDB)
               ABSP2S(II,JJ) = ABSP2(II,JJ,INDB)

               REFLIRS(II,JJ) = REFLIR(II,JJ,INDB)
               DIFUIRS(II,JJ) = DIFUIR(II,JJ,INDB)
               ABSPIRS(II,JJ) = ABSPIR(II,JJ,INDB)
            ENDDO
         ENDDO
      ENDIF
C      write(*,*) "optical properties per block"
C ----------------------
C EARTH RADIATION MODELS
C ----------------------

      IF(ERM.GT.0)THEN
         ABSSUN = DSQRT(SUN(1)**2 + SUN(2)**2 + SUN(3)**2)
         DO K=1,3
            ESUN(K) = SUN(K)/ABSSUN
         ENDDO

         PSIDOT = ESUN(1)*RADVEC(1)+ESUN(2)*RADVEC(2)+ESUN(3)*RADVEC(3)
         IF(DABS(PSIDOT).GT.(1D0-1D-6))THEN
            PSI = 0D0
         ELSE
            PSI = DACOS(PSIDOT)
         ENDIF
         S0 = S0*(AU/ABSSUN)**2
      ENDIF

C     ANALYTICAL MODEL
      IF(ERM.EQ.1)THEN

         NCFVEC(1) = RADVEC(1)
         NCFVEC(2) = RADVEC(2)
         NCFVEC(3) = RADVEC(3)
         ALBFAC = (PI*TOA**2)*(S0/C)/(ABSPOS**2)
         PHASEVI = (2*ALB/(3*PI**2))*((PI-PSI)*DCOS(PSI)+DSIN(PSI))
         PHASEIR = (1-ALB)/(4*PI)
         ABSNCFVI = ALBFAC*PHASEVI
         ABSNCFIR = ALBFAC*PHASEIR

         CALL SURFBOXW(AREAS,REFLS,DIFUS,ABSPS,
     1                    AREA2S,REFL2S,DIFU2S,ABSP2S,
     2                    REFLIRS,DIFUIRS,ABSPIRS,
     3                    ABSNCFVI,ABSNCFIR,NCFVEC,ATTSURF,FORCE)



C     NUMERICAL MODEL (CERES DATA)
      ELSEIF(ERM.EQ.2)THEN

         ABSSUN = DSQRT(SUN(1)**2 + SUN(2)**2 + SUN(3)**2)
         DO K=1,3
            ESUN(K) = SUN(K)/ABSSUN
         ENDDO
C         write(*,*) "ERP==2", LATKMX, LONKMX
         DO LATK = 1,LATKMX
            DO LONK = 1,LONKMX

               D_AREA = D_AREA_ALL(LATK,LONK)
               V_NS(1) = V_NS_ALL(LATK,LONK,1)
               V_NS(2) = V_NS_ALL(LATK,LONK,2)
               V_NS(3) = V_NS_ALL(LATK,LONK,3)

               DO II=1,3
                  V_INS(II)=0D0
                  DO JJ=1,3
                     V_INS(II) = V_INS(II) + KAPPA(JJ,II)*V_NS(JJ)
                  ENDDO
               ENDDO

C              Distance and direction from point in the Earth to satellite
               V_DIST(1) = YSAT(1)-TOA*V_INS(1) 
               V_DIST(2) = YSAT(2)-TOA*V_INS(2)
               V_DIST(3) = YSAT(3)-TOA*V_INS(3)
               DIST2 = V_DIST(1)**2 +V_DIST(2)**2 +V_DIST(3)**2
               ABSDIST = DSQRT(DIST2)
               V_SAT(1) = V_DIST(1)/ABSDIST
               V_SAT(2) = V_DIST(2)/ABSDIST
               V_SAT(3) = V_DIST(3)/ABSDIST

C              Cosine of angles of incident and reflected radiation
               COS_IN = ESUN(1)*V_INS(1) + ESUN(2)*V_INS(2) 
     1                + ESUN(3)*V_INS(3)
               COS_RE =V_SAT(1)*V_INS(1) +V_SAT(2)*V_INS(2)
     1                +V_SAT(3)*V_INS(3)

               IF(COS_RE.GE.0)THEN

C              Reflectivity and emissivity coefficients
                  REFL_CF = CERGRE(LATK,LONK)
                  EMIT_CF = CERGEM(LATK,LONK)

C                 Reflected Irradiance
                  IF(COS_IN.GE.0)THEN
                     E_REFL=(REFL_CF/(PI*DIST2))*COS_RE*COS_IN*S0*D_AREA
                  ELSE
                     E_REFL=0D0
                  ENDIF

C                 Emitted Irradiance
                  E_EMIT = (EMIT_CF/(4*PI*DIST2))*COS_RE*S0*D_AREA

C                 Non-conservative force
                  ABSNCFVI = E_REFL/C
                  ABSNCFIR = E_EMIT/C
                  NCFVEC(1) = V_SAT(1)
                  NCFVEC(2) = V_SAT(2)
                  NCFVEC(3) = V_SAT(3)
C                  write(*,*) "LATK, LONK", LATK, LONK   
                  CALL SURFBOXW(AREAS,REFLS,DIFUS,ABSPS,
     1                    AREA2S,REFL2S,DIFU2S,ABSP2S,
     2                    REFLIRS,DIFUIRS,ABSPIRS,
     3                    ABSNCFVI,ABSNCFIR,NCFVEC,ATTSURF,FORCE_LL)
C                  write(*,*) "getPARAM",DIFU2S
C                  , DIFU2S, ABSP2S, REFLIRS,DIFUIRS,ABSPIRS,ABSNCFVI
C                  ,AREAS, REFLS, DIFUS, ABSPS,
C                   ,ABSNCFIR,NCFVEC,ATTSURF
                  FORCE(1) = FORCE(1) + FORCE_LL(1)
                  FORCE(2) = FORCE(2) + FORCE_LL(2)
                  FORCE(3) = FORCE(3) + FORCE_LL(3)
               ENDIF

            ENDDO
         ENDDO
C         write(*,*) "ERP model quit"
      ENDIF

C      write(*,*) "earth radiation model"


C     NAVIGATION ANTENNA THRUST (SIMPLE MODEL)
      IF((ANT.EQ.1).AND.(BLKNUM.LT.100))THEN
         ANTFORCE = 80D0/C
         FORCE(1) = FORCE(1) + ANTFORCE*RADVEC(1)
         FORCE(2) = FORCE(2) + ANTFORCE*RADVEC(2)
         FORCE(3) = FORCE(3) + ANTFORCE*RADVEC(3)
      ENDIF

      IF((REFF.GT.0).AND.((ERM.GT.0).OR.(ANT.EQ.1)))THEN
         DO K=1,3
            FREF(K) = 0D0
         ENDDO

C        FORCE IN BODY-FIXED REFERENCE FRAME
         IF(REFF.EQ.1)THEN
            DO K=1,3
               FREF(1) = FREF(1) + FORCE(K)*Z_SAT(K)
               FREF(2) = FREF(2) + FORCE(K)*Y_SAT(K)
               FREF(3) = FREF(3) + FORCE(K)*X_SAT(K)
            ENDDO

C        FORCE IN SUN-FIXED REFERENCE FRAME
         ELSEIF(REFF.EQ.2)THEN
            DO K=1,3
               FREF(1) = FREF(1) + FORCE(K)*D_SUN(K)
               FREF(2) = FREF(2) + FORCE(K)*Y_SAT(K)
               FREF(3) = FREF(3) + FORCE(K)*B_SUN(K)
            ENDDO

C        FORCE IN ORBITAL REFERENCE FRAME
         ELSEIF(REFF.EQ.3)THEN
            DO K=1,3
               FREF(1) = FREF(1) + FORCE(K)*RADVEC(K)
               FREF(2) = FREF(2) + FORCE(K)*ALGVEC(K)
               FREF(3) = FREF(3) + FORCE(K)*CRSVEC(K)
            ENDDO
         ENDIF

         DO K=1,3
            FORCE(K) = FREF(K)
         ENDDO
      ENDIF

C     MASS OF SATELLITES
      IF(BLKNUM.EQ.1)THEN
         MASS = 500D0
      ELSEIF(BLKNUM.EQ.2)THEN 
         MASS = 880D0
      ELSEIF(BLKNUM.EQ.3)THEN
         MASS = 975D0
      ELSEIF((BLKNUM.GE.4).AND.(BLKNUM.LE.7))THEN
         MASS = 1100D0
      ELSEIF((BLKNUM.EQ.101).OR.(BLKNUM.EQ.102))THEN
         MASS = 1415D0
      ELSEIF(BLKNUM.EQ.103)THEN
         MASS = 935D0
      ENDIF
C      write(*,*) "get mass"
      
C     CONVERSION TO ACCELERATION
      DO K=1,3
         ACCEL(K) = FORCE(K)/MASS
      ENDDO

      END SUBROUTINE

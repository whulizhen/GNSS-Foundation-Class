C*
      SUBROUTINE SURFBOXW(AREA,REFL,DIFU,ABSP,AREA2,REFL2,DIFU2,
     1                    ABSP2,REFLIR,DIFUIR,ABSPIR,ABSNCFVI,ABSNCFIR,
     2                    NCFVEC,ATTSURF,FORCE)
CC
CC NAME       :  SURFBOXW
CC
CC PURPOSE    :  COMPUTES THE FORCE IN THE SATELLITE SURFACES, FOR A GIVEN RADIATION VECTOR
CC
CC PARAMETERS :
CC         IN :  AREA(I,J)   : AREAS OF FLAT SURFACES [m^2]
CC               REFL(I,J)   : REFLEXION COEFFICIENT
CC               DIFU(I,J)   : DIFFUSION COEFFICIENT
CC               ABSP(I,J)   : ABSORPTION COEFFICIENT
CC                             I = 1 +Z (TOWARDS THE EARTH)
CC                             I = 2 +Y (ALONG SOLAR PANELS BEAMS)
CC                             I = 3 +X (ALONG BUS DIRECTION ALWAYS ILLUMINATED BY THE SUN)
CC                             I = 4 SOLAR PANELS
CC                             J = 1 POSITIVE DIRECTION
CC                             J = 2 NEGATIVE DIRECTION
CC               ????2(I,J)  : INDEX 2 INDICATES AREAS AND OPTICAL PROPERTIES OF CYLINDRICAL
CC                             SURFACES, SAME MEANING OF (I,J) AS BEFORE
CC               ????IR(I,J ): OPTICAL PROPERTIES IN THE INFRARED, NO SEPARATION
CC                             BETWEEN FLAT AND CYLINDRICAL SURFACES
CC               ABSNCFVI    : MAGNITUDE OF INCIDENT VISIBLE RADIATION [N/m^2]
CC               ABSNCFIR    : MAGNITUDE OF INCIDENT INFRARED RADIATION [N/m^2]
CC               NCFVEC      : UNITARY INCIDENT RADIATION VECTOR IN INERTIAL REFERENCE FRAME
CC               ATTSURF(K,I): ATTITUDE VECTORS OF SATELLITE SURFACES
CC                             K = 1,2,3 VECTOR COMPONENTS IN INERTIAL REFERENCE FRAME
CC                             I = 1 +Z (TOWARDS THE EARTH)
CC                             I = 2 +Y (ALONG SOLAR PANELS BEAMS)
CC                             I = 3 +X (ALONG BUS DIRECTION ALWAYS ILLUMINATED BY THE SUN)
CC                             I = 4 SOLAR PANELS
CC
CC        OUT : FORCE        : RESULTING FORCE VECTOR [Newtons] IN INERTIAL REFERENCE FRAME
CC
CC AUTHOR     :  C.J. RODRIGUEZ-SOLANO
CC               rodriguez@bv.tum.de
CC
CC VERSION    :  1.0 (OCT 2010)
CC
CC CREATED    :  2010/10/18
CC
C*
      IMPLICIT NONE
C
      INTEGER*4 K,I
      INTEGER*4 SIDE(4)
C
      REAL*8 ABSNCFVI,ABSNCFIR,PI
      REAL*8 FACABDI,FACREFL
      REAL*8 NCFVEC(3),FORCE(3)
      REAL*8 ATTSURF(3,4),NORVEC(3,4)
      REAL*8 AREA(4,2),REFL(4,2),DIFU(4,2),ABSP(4,2)
      REAL*8 AREA2(4,2),REFL2(4,2),DIFU2(4,2),ABSP2(4,2)
      REAL*8 REFLIR(4,2),DIFUIR(4,2),ABSPIR(4,2)
      REAL*8 COSANG(4)
      REAL*8 TABSP(3),TREFL(3),TDIFU(3)
      REAL*8 TABSP2(3),TREFL2(3),TDIFU2(3)
      REAL*8 PRT_FLT(3,3),PRT_CYL(3,3),PRT_FLTIR(3,3),PRT_CYLIR(3,3)
      REAL*8 FVI_FLT(3),FVI_CYL(3),FIR_FLT(3),FIR_CYL(3)
      REAL*8 FVI(3),FIR(3)

C Initialization of variables
      PI = 4D0*DATAN(1D0)

      FORCE(1) = 0D0
      FORCE(2) = 0D0
      FORCE(3) = 0D0

C --------------------------------
C ANGLE BETWEEN FORCE AND SURFACES
C --------------------------------
      SIDE(1) = 1
      SIDE(2) = 1
      SIDE(3) = 1
      SIDE(4) = 1

      DO K=1,4
         COSANG(K) = ATTSURF(1,K)*NCFVEC(1)+ATTSURF(2,K)*NCFVEC(2)
     1            +ATTSURF(3,K)*NCFVEC(3)

         IF(COSANG(K).GE.0D0)THEN
            SIDE(K) = 2
            NORVEC(1,K) = ATTSURF(1,K)
            NORVEC(2,K) = ATTSURF(2,K)
            NORVEC(3,K) = ATTSURF(3,K)
         ELSEIF(COSANG(K).LT.0D0)THEN
            SIDE(K) = 1
            NORVEC(1,K) = -ATTSURF(1,K)
            NORVEC(2,K) = -ATTSURF(2,K)
            NORVEC(3,K) = -ATTSURF(3,K)
         ENDIF
      ENDDO


C --------------------------------------
C FORCE ACTING ON THE SATELLITE SURFACES
C --------------------------------------
      DO K=1,4
         FACABDI = DABS(COSANG(K))
         FACREFL = COSANG(K)**2

         DO I=1,3
            TABSP(I) = FACABDI*NCFVEC(I)
            TABSP2(I)= TABSP(I)

            TREFL(I)   =   (2D0)*FACREFL*NORVEC(I,K)
            TREFL2(I)= (4D0/3D0)*FACREFL*NORVEC(I,K)

            TDIFU(I) = FACABDI*(NCFVEC(I)+(2D0/3D0)*NORVEC(I,K))    
            TDIFU2(I)= FACABDI*(NCFVEC(I) +(PI/6D0)*NORVEC(I,K))

C           SATELLITE BUS, INCLUDES DIFFUSE THERMAL RE-RADIATION
            IF(K.LE.3)THEN
               PRT_FLT(1,I) = TDIFU(I)*(ABSP(K,SIDE(K))+DIFU(K,SIDE(K)))
               PRT_FLT(2,I) = TREFL(I)*REFL(K,SIDE(K))
               PRT_FLT(3,I) = 0D0
               PRT_CYL(1,I) = TDIFU2(I)*(ABSP2(K,SIDE(K))
     1                                  +DIFU2(K,SIDE(K)))
               PRT_CYL(2,I) = TREFL2(I)*REFL2(K,SIDE(K))
               PRT_CYL(3,I) = 0D0

               PRT_FLTIR(1,I) = TDIFU(I)*(ABSPIR(K,SIDE(K))
     1                                   +DIFUIR(K,SIDE(K)))
               PRT_FLTIR(2,I) = TREFL(I)*REFLIR(K,SIDE(K))
               PRT_FLTIR(3,I) = 0D0
               PRT_CYLIR(1,I) = TDIFU2(I)*(ABSPIR(K,SIDE(K))
     1                                    +DIFUIR(K,SIDE(K)))
               PRT_CYLIR(2,I) = TREFL2(I)*REFLIR(K,SIDE(K))
               PRT_CYLIR(3,I) = 0D0

C        SOLAR PANELS, THERMAL RE-RADIATION ALMOST CANCELS (FRONT AND BACK)
            ELSEIF(K.EQ.4)THEN
               PRT_FLT(1,I) = TABSP(I)*ABSP(K,SIDE(K))
               PRT_FLT(2,I) = TREFL(I)*REFL(K,SIDE(K))
               PRT_FLT(3,I) = TDIFU(I)*DIFU(K,SIDE(K))
               PRT_CYL(1,I) = TABSP2(I)*ABSP2(K,SIDE(K))
               PRT_CYL(2,I) = TREFL2(I)*REFL2(K,SIDE(K))
               PRT_CYL(3,I) = TDIFU2(I)*DIFU2(K,SIDE(K))
               PRT_FLTIR(1,I) = TABSP(I)*ABSPIR(K,SIDE(K))
               PRT_FLTIR(2,I) = TREFL(I)*REFLIR(K,SIDE(K))
               PRT_FLTIR(3,I) = TDIFU(I)*DIFUIR(K,SIDE(K))
               PRT_CYLIR(1,I) = TABSP2(I)*ABSPIR(K,SIDE(K))
               PRT_CYLIR(2,I) = TREFL2(I)*REFLIR(K,SIDE(K))
               PRT_CYLIR(3,I) = TDIFU2(I)*DIFUIR(K,SIDE(K))
            ENDIF
            
C           TOTAL FORCE COMPUTATION
            FVI_FLT(I) = PRT_FLT(1,I)  + PRT_FLT(2,I)  + PRT_FLT(3,I)
            FVI_CYL(I) = PRT_CYL(1,I)  + PRT_CYL(2,I)  + PRT_CYL(3,I)
            FIR_FLT(I) = PRT_FLTIR(1,I)+ PRT_FLTIR(2,I)+ PRT_FLTIR(3,I)
            FIR_CYL(I) = PRT_CYLIR(1,I)+ PRT_CYLIR(2,I)+ PRT_CYLIR(3,I)

            FVI_FLT(I) = FVI_FLT(I)*AREA(K,SIDE(K)) 
            FVI_CYL(I) = FVI_CYL(I)*AREA2(K,SIDE(K))
            FIR_FLT(I) = FIR_FLT(I)*AREA(K,SIDE(K))
            FIR_CYL(I) = FIR_CYL(I)*AREA2(K,SIDE(K))
         

            FVI(I) = (FVI_FLT(I)+FVI_CYL(I))*ABSNCFVI
            FIR(I) = (FIR_FLT(I)+FIR_CYL(I))*ABSNCFIR
            

            FORCE(I) = FORCE(I) + FVI(I) + FIR(I)

         ENDDO
      ENDDO

      END SUBROUTINE

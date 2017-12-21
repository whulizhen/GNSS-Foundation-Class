C*
      SUBROUTINE PROPBOXW(BLKNUM,AREA,REFL,DIFU,ABSP,AREA2,REFL2,DIFU2,
     1                    ABSP2,REFLIR,DIFUIR,ABSPIR)
CC
CC NAME       :  PROPBOXW
CC
CC PURPOSE    :  SATELLITE DIMENSIONS AND OPTICAL PROPERTIES FROM ROCK MODELS
CC               BOX-WING MODELS FOR GPS AND GLONASS SATELLITES
CC
CC REFERENCES :  Fliegel H, Gallini T, Swift E (1992) Global Positioning System Radiation
CC                  Force Model for Geodetic Applications. Journal of Geophysical Research
CC                  97(B1): 559-568
CC               Fliegel H, Gallini T (1996) Solar Force Modelling of Block IIR Global
CC                  Positioning System satellites. Journal of Spacecraft and Rockets
CC                  33(6): 863-866
CC               Ziebart M (2001) High Precision Analytical Solar Radiation Pressure
CC                  Modelling for GNSS Spacecraft. PhD Thesis, University of East London
CC               [IGEXMAIL-0086] GLONASS S/C mass and dimension
CC               [IGSMAIL-5104] GLONASS-M dimensions and center-of-mass correction
CC               http://acc.igs.org/orbits/IIF_SV_DimensionsConfiguration.ppt
CC
CC PARAMETERS : 
CC         IN :  BLKNUM     : BLOCK NUMBER
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
CC        OUT :  AREA(I,J)  : AREAS OF FLAT SURFACES [m^2]
CC               REFL(I,J)  : REFLEXION COEFFICIENT
CC               DIFU(I,J)  : DIFFUSION COEFFICIENT
CC               ABSP(I,J)  : ABSORPTION COEFFICIENT
CC                            I = 1 +Z (TOWARDS THE EARTH)
CC                            I = 2 +Y (ALONG SOLAR PANELS BEAMS)
CC                            I = 3 +X (ALONG BUS DIRECTION ALWAYS ILLUMINATED BY THE SUN)
CC                            I = 4 SOLAR PANELS
CC                            J = 1 POSITIVE DIRECTION
CC                            J = 2 NEGATIVE DIRECTION
CC               ????2(I,J) : INDEX 2 INDICATES AREAS AND OPTICAL PROPERTIES OF CYLINDRICAL
CC                            SURFACES, SAME MEANING OF (I,J) AS BEFORE
CC               ????IR(I,J): OPTICAL PROPERTIES IN THE INFRARED (ASSUMED), NO SEPARATION
CC                            BETWEEN FLAT AND CYLINDRICAL SURFACES
CC
CC AUTHOR     :  C.J. RODRIGUEZ-SOLANO
CC               rodriguez@bv.tum.de
CC
CC VERSION    :  1.0 (OCT 2010)
CC
CC CREATED    :  2010/10/18             LAST MODIFIED :  17-MAR-11
CC
CC CHANGES    :  17-MAR-11 : CR: ADD BOX-WING MODELS FOR GLONASS, GLONASS-M, GPS II-F
C*
      IMPLICIT NONE
C
      INTEGER*4 BLKNUM,II,JJ,KK,SS
C
      REAL*8 AREA(4,2),REFL(4,2),DIFU(4,2),ABSP(4,2)
      REAL*8 AREA2(4,2),REFL2(4,2),DIFU2(4,2),ABSP2(4,2)
      REAL*8 REFLIR(4,2),DIFUIR(4,2),ABSPIR(4,2)
      REAL*8 X_SIDE(5,4),Z_SIDE(4,4),S_SIDE(2,4),Y_SIDE(2,4)
      REAL*8 SURFALL(13,4)
      REAL*8 ALLREFL(13),ALLDIFU(13),REFL_AREA(13),DIFU_AREA(13)
      REAL*8 G_AREA1(5),G_REFL1(5),G_DIFU1(5),G_ABSP1(5)
      REAL*8 G_AREA2(5),G_REFL2(5),G_DIFU2(5),G_ABSP2(5)
      REAL*8 BUSFAC

C SATELLITE PROPERTIES FROM ROCK MODEL
C VALUES GIVEN FOR +X, -Z, +Z AND SOLAR PANELS
C X_SIDE(I,J) :   I     DIFFERENT SURFACE COMPONENTS
C                 J = 1 AREA
C                 J = 2 SPECULARITY
C                 J = 3 REFLECTIVITY
C                 J = 4 SHAPE: 1 = flat, 2 = cylindrical

      DO JJ = 1,4
         DO II = 1,5
            X_SIDE(II,JJ) = 0D0
         ENDDO
         DO II = 1,4
            Z_SIDE(II,JJ) = 0D0
         ENDDO
         DO II = 1,2
            S_SIDE(II,JJ) = 0D0
         ENDDO
         DO II = 1,2
            Y_SIDE(II,JJ) = 0D0
         ENDDO
      ENDDO


C     -----------
C     GPS BLOCK I
C     -----------
C     SEE FLIEGEL ET AL (1992)
      IF(BLKNUM.EQ.1)THEN
          
C     +X SIDE
          X_SIDE(1,1) = 1.055D0
          X_SIDE(1,2) = 0.80D0
          X_SIDE(1,3) = 0.50D0
          X_SIDE(1,4) = 1D0

C     ENGINE SIDE
          X_SIDE(2,1) = 0.570D0
          X_SIDE(2,2) = 0.75D0
          X_SIDE(2,3) = 0.86D0
          X_SIDE(2,4) = 2D0

C     TT&C ANTENNA SIDE
          X_SIDE(3,1) = 0.055D0
          X_SIDE(3,2) = 0.05D0
          X_SIDE(3,3) = 0.28D0
          X_SIDE(3,4) = 2D0

C     TT&C ANTENNA TIP
          X_SIDE(4,1) = 0.019D0
          X_SIDE(4,2) = 0.85D0
          X_SIDE(4,3) = 0.28D0
          X_SIDE(4,4) = 2D0

C     EACH NAVIGATIONAL ANTENNA ADAPTER
          X_SIDE(5,1) = 0.029D0
          X_SIDE(5,2) = 0.75D0
          X_SIDE(5,3) = 0.36D0
          X_SIDE(5,4) = 2D0

C     BODY AFT END (-Z)
          Z_SIDE(1,1) = 0.816D0
          Z_SIDE(1,2) = 0.80D0
          Z_SIDE(1,3) = 0.86D0
          Z_SIDE(1,4) = 1D0

C     ENGINE AFT END (-Z)                                                         
          Z_SIDE(2,1) = 0.694D0
          Z_SIDE(2,2) = 0D0
          Z_SIDE(2,3) = 0D0
          Z_SIDE(2,4) = 1D0

C     FORWARD END (+Z)                                                                     
          Z_SIDE(3,1) = 1.510D0
          Z_SIDE(3,2) = 0.75D0
          Z_SIDE(3,3) = 0.86D0
          Z_SIDE(3,4) = 1D0

C     ALL SOLAR PANELS
          S_SIDE(1,1) = 5.583D0
          S_SIDE(1,2) = 0.85D0
          S_SIDE(1,3) = 0.23D0
          S_SIDE(1,4) = 1D0

C     SOLAR PANELS MASTS                                                                                     
          S_SIDE(2,1) = 0.470D0
          S_SIDE(2,2) = 0.85D0
          S_SIDE(2,3) = 0.85D0
          S_SIDE(2,4) = 2D0

C     -----------------
C     GPS BLOCK II, IIA
C     -----------------
C     SEE FLIEGEL ET AL (1992)
      ELSEIF((BLKNUM.EQ.2).OR.(BLKNUM.EQ.3))THEN

C     +X SIDE                                                                                           
          X_SIDE(1,1) = 1.553D0
          X_SIDE(1,2) = 0.20D0
          X_SIDE(1,3) = 0.56D0
          X_SIDE(1,4) = 1D0

C     ENGINE SIDE (INCLUDES PLUME SHILED 0.22*1.84 m^2)                                            
          X_SIDE(2,1) = 1.054D0
          X_SIDE(2,2) = 0.20D0
          X_SIDE(2,3) = 0.56D0
          X_SIDE(2,4) = 2D0

C     TT&C ANTENNA SIDE                                                                       
          X_SIDE(3,1) = 0.105D0
          X_SIDE(3,2) = 0.20D0
          X_SIDE(3,3) = 0.28D0
          X_SIDE(3,4) = 2D0

C     EACH NAVIGATIONAL ANTENNA ADAPTER                                                               
          X_SIDE(5,1) = 0.181D0
          X_SIDE(5,2) = 0.20D0
          X_SIDE(5,3) = 0.36D0
          X_SIDE(5,4) = 2D0

C     BODY AFT END (-Z)                                                                               
          Z_SIDE(1,1) = 2.152D0
          Z_SIDE(1,2) = 0.20D0
          Z_SIDE(1,3) = 0.56D0
          Z_SIDE(1,4) = 1D0

C     ENGINE AFT END (-Z)                                                                          
          Z_SIDE(2,1) = 0.729D0
          Z_SIDE(2,2) = 0D0
          Z_SIDE(2,3) = 0D0
          Z_SIDE(2,4) = 1D0

C     FORWARD END (+Z)                                                                                    
          Z_SIDE(3,1) = 2.881D0
          Z_SIDE(3,2) = 0.20D0
          Z_SIDE(3,3) = 0.56D0
          Z_SIDE(3,4) = 1D0

C     ALL SOLAR PANELS                                                                                     
          S_SIDE(1,1) = 10.886D0
          S_SIDE(1,2) = 0.85D0
          S_SIDE(1,3) = 0.23D0
          S_SIDE(1,4) = 1D0

C     SOLAR PANELS MASTS                                                                    
          S_SIDE(2,1) = 0.985D0
          S_SIDE(2,2) = 0.41D0
          S_SIDE(2,3) = 0.52D0
          S_SIDE(2,4) = 2D0

C     ----------------------------------
C     GPS BLOCK IIR, IIR-A, IIR-B, IIR-M
C     ----------------------------------
C     SEE FLIEGEL AND GALINI (1996)
      ELSEIF((BLKNUM.GE.4).AND.(BLKNUM.LE.7))THEN

C     + AND -X FACES                                                                                
         X_SIDE(1,1) = 3.05D0
         X_SIDE(1,2) = 0D0
         X_SIDE(1,3) = 0.06D0
         X_SIDE(1,4) = 1D0

C     PLUME SHILED                                                                                  
         X_SIDE(2,1) = 0.17D0
         X_SIDE(2,2) = 0D0
         X_SIDE(2,3) = 0.06D0
         X_SIDE(2,4) = 2D0

C     ANTENNA SHROUD                                                                              
         X_SIDE(3,1) = 0.89D0
         X_SIDE(3,2) = 0D0
         X_SIDE(3,3) = 0.06D0
         X_SIDE(3,4) = 2D0

C     -Z FACE                                                                                        
         Z_SIDE(1,1) = 3.75D0
         Z_SIDE(1,2) = 0D0
         Z_SIDE(1,3) = 0.06D0
         Z_SIDE(1,4) = 1D0

C     -Z W-SENSOR
         Z_SIDE(2,1) = 0.50D0
         Z_SIDE(2,2) = 0D0
         Z_SIDE(2,3) = 0.06D0
         Z_SIDE(2,4) = 1D0

C     +Z FACE                                                                                              
         Z_SIDE(3,1) = 3.75D0
         Z_SIDE(3,2) = 0D0
         Z_SIDE(3,3) = 0.06D0
         Z_SIDE(3,4) = 1D0

C     +Z W-SENSOR                                                                                             
         Z_SIDE(4,1) = 0.50D0
         Z_SIDE(4,2) = 0D0
         Z_SIDE(4,3) = 0.06D0
         Z_SIDE(4,4) = 1D0

C     ALL SOLAR PANELS                                                                                
         S_SIDE(1,1) = 13.60D0
         S_SIDE(1,2) = 0.85D0
         S_SIDE(1,3) = 0.28D0
         S_SIDE(1,4) = 1D0

C     SOLAR PANELS BEAMS                                                                     
         S_SIDE(2,1) = 0.32D0
         S_SIDE(2,2) = 0.85D0
         S_SIDE(2,3) = 0.85D0
         S_SIDE(2,4) = 1D0

C     -------------
C     GPS BLOCK IIF
C     -------------
C     SEE IIF PRESENTATION
C     OPTICAL PROPERTIES NOT KNOWN, SAME ASSUMPTIONS AS ZIEBART (2001)
      ELSEIF(BLKNUM.EQ.8)THEN
C     +/- X SIDE
         X_SIDE(1,1) = 5.72D0
         X_SIDE(1,2) = 0.20D0
         X_SIDE(1,3) = 0.56D0
         X_SIDE(1,4) = 1D0

C     +/- Y SIDE
         Y_SIDE(1,1) = 7.01D0
         Y_SIDE(1,2) = 0.20D0
         Y_SIDE(1,3) = 0.56D0
         Y_SIDE(1,4) = 1D0

C     -Z SIDE
         Z_SIDE(1,1) = 5.40D0
         Z_SIDE(1,2) = 0D0
         Z_SIDE(1,3) = 0D0
         Z_SIDE(1,4) = 1D0

C     +Z SIDE
         Z_SIDE(3,1) = 5.40D0
         Z_SIDE(3,2) = 0D0
         Z_SIDE(3,3) = 0D0
         Z_SIDE(3,4) = 1D0

C     SOLAR PANELS
         S_SIDE(1,1) = 22.25D0
         S_SIDE(1,2) = 0.85D0
         S_SIDE(1,3) = 0.23D0
         S_SIDE(1,4) = 1D0

C     -------
C     GLONASS
C     -------
C     SEE ZIEBART (2001)
C     OPTICAL PROPERTIES NOT KNOWN, SAME ASSUMPTIONS AS ZIEBART (2001)
C     ASSUMED SAME SHAPE FOR GLONASS-M, SEE END OF SUBROUTINE
      ELSEIF((BLKNUM.EQ.101).OR.(BLKNUM.EQ.102))THEN
C     +/- X SIDE (FLAT)
         X_SIDE(1,1) = 1.258D0
         X_SIDE(1,2) = 0.20D0
         X_SIDE(1,3) = 0.56D0
         X_SIDE(1,4) = 1D0

C     +/- X SIDE (CYLINDRICAL)
         X_SIDE(2,1) = 2.052D0
         X_SIDE(2,2) = 0.20D0
         X_SIDE(2,3) = 0.56D0
         X_SIDE(2,4) = 2D0

C     +/- Y SIDE (FLAT)
         Y_SIDE(1,1) = 2.591D0
         Y_SIDE(1,2) = 0.20D0
         Y_SIDE(1,3) = 0.56D0
         Y_SIDE(1,4) = 1D0

C     +/- Y SIDE (CYLINDRICAL)
         Y_SIDE(2,1) = 2.532D0
         Y_SIDE(2,2) = 0.20D0
         Y_SIDE(2,3) = 0.56D0
         Y_SIDE(2,4) = 2D0

C     -Z SIDE (BUS)
         Z_SIDE(1,1) = 0.877D0
         Z_SIDE(1,2) = 0.20D0
         Z_SIDE(1,3) = 0.56D0
         Z_SIDE(1,4) = 1D0

C     -Z SIDE (APOGEE ENGINE)
         Z_SIDE(2,1) = 0.785D0
         Z_SIDE(2,2) = 0D0
         Z_SIDE(2,3) = 0D0
         Z_SIDE(2,4) = 1D0

C     +Z SIDE (BUS)
         Z_SIDE(3,1) = 1.412D0
         Z_SIDE(3,2) = 0.20D0
         Z_SIDE(3,3) = 0.56D0
         Z_SIDE(3,4) = 1D0

C     +Z SIDE (RETRO REFLECTOR ARRAY)
         Z_SIDE(4,1) = 0.250D0
         Z_SIDE(4,2) = 1D0
         Z_SIDE(4,3) = 1D0
         Z_SIDE(4,4) = 1D0

C     SOLAR PANELS
         S_SIDE(1,1) = 23.616D0
         S_SIDE(1,2) = 0.85D0
         S_SIDE(1,3) = 0.23D0
         S_SIDE(1,4) = 1D0

      ENDIF


C ---------------------------
C BOX-WING MODEL APROXIMATION
C ---------------------------

C NAVIGATION ANTENNAS
      X_SIDE(5,1) = X_SIDE(5,1)*12

C MATRIX WITH ALL SURFACES
      DO II = 1,13
         DO JJ = 1,4
            IF(II.LE.5)THEN
               SURFALL(II,JJ) = X_SIDE(II,JJ)
            ELSEIF((II.GE.6).AND.(II.LE.9))THEN
               SURFALL(II,JJ) = Z_SIDE(II-5,JJ)
            ELSEIF((II.GE.10).AND.(II.LE.11))THEN
               SURFALL(II,JJ) = S_SIDE(II-9,JJ)
            ELSEIF(II.GE.12)THEN
               SURFALL(II,JJ) = Y_SIDE(II-11,JJ)
            ENDIF
         ENDDO
      ENDDO
      

C COMPUTATION OF FRACTION OF REFLECTED, DIFFUSSED AND ABSORVED PHOTONS
C FROM SPECULARITY AND REFLECTIVITY
      DO II = 1,13
         ALLREFL(II) = SURFALL(II,2)*SURFALL(II,3)
         ALLDIFU(II) = SURFALL(II,3)*(1-SURFALL(II,2))
      ENDDO

C MULTIPLICATION OF OPTICAL PROPERTIES WITH SURFACE AREA
      DO II = 1,13
         REFL_AREA(II) = ALLREFL(II)*SURFALL(II,1)
         DIFU_AREA(II) = ALLDIFU(II)*SURFALL(II,1)
      ENDDO
      
C AVERAGE PER SURFACE (SEPARATE FOR FLAT AND CYLINDRICAL SURFACES)
      DO SS = 1,5
         G_AREA1(SS) = 0D0
         G_REFL1(SS) = 0D0
         G_DIFU1(SS) = 0D0
         G_ABSP1(SS) = 0D0
         G_AREA2(SS) = 0D0
         G_REFL2(SS) = 0D0
         G_DIFU2(SS) = 0D0
         G_ABSP2(SS) = 0D0
      ENDDO


      DO II = 1,13
C        +X SIDE 
         IF(II.LE.5)THEN
            SS = 1
C        -Z SIDE
         ELSEIF((II.GE.6).AND.(II.LE.7))THEN
            SS = 2
C        +Z SIDE
         ELSEIF((II.GE.8).AND.(II.LE.9))THEN
            SS = 3
C        SOLAR PANELS
         ELSEIF((II.GE.10).AND.(II.LE.11))THEN
            SS = 4
C        +/- Y SIDE
         ELSEIF(II.GE.12)THEN
            SS = 5
         ENDIF

         IF(SURFALL(II,4).EQ.1D0)THEN
            G_AREA1(SS) = G_AREA1(SS) + SURFALL(II,1)
            G_REFL1(SS) = G_REFL1(SS) + REFL_AREA(II)
            G_DIFU1(SS) = G_DIFU1(SS) + DIFU_AREA(II)
         ELSEIF(SURFALL(II,4).EQ.2D0)THEN
            G_AREA2(SS) = G_AREA2(SS) + SURFALL(II,1)
            G_REFL2(SS) = G_REFL2(SS) + REFL_AREA(II)
            G_DIFU2(SS) = G_DIFU2(SS) + DIFU_AREA(II)
         ENDIF
      ENDDO

C AVERAGE OF OPTICAL PROPERTIES ACCORDING TO SURFACES AREA
      DO SS = 1,5
         IF(G_AREA1(SS).GT.0D0)THEN
            G_REFL1(SS) = G_REFL1(SS)/G_AREA1(SS)
            G_DIFU1(SS) = G_DIFU1(SS)/G_AREA1(SS)
            G_ABSP1(SS) = 1D0 - (G_REFL1(SS) + G_DIFU1(SS))
         ENDIF
         IF(G_AREA2(SS).GT.0D0)THEN
            G_REFL2(SS) = G_REFL2(SS)/G_AREA2(SS)
            G_DIFU2(SS) = G_DIFU2(SS)/G_AREA2(SS)
            G_ABSP2(SS) = 1D0 -(G_REFL2(SS) + G_DIFU2(SS))
         ENDIF
      ENDDO


C --------------------------------------------
C ARRAYS OF OPTICAL PROPERTIES AND DIMENSIONS
C --------------------------------------------

      DO SS = 1,4
         DO KK = 1,2
            AREA(SS,KK) = 0D0
            REFL(SS,KK) = 0D0
            DIFU(SS,KK) = 0D0
            ABSP(SS,KK) = 0D0
            AREA2(SS,KK) = 0D0
            REFL2(SS,KK) = 0D0
            DIFU2(SS,KK) = 0D0
            ABSP2(SS,KK) = 0D0
         ENDDO
      ENDDO

C     +Z FACE      
      AREA(1,1) = G_AREA1(3)
      REFL(1,1) = G_REFL1(3)
      DIFU(1,1) = G_DIFU1(3)
      ABSP(1,1) = G_ABSP1(3)

C     -Z FACE
      AREA(1,2) = G_AREA1(2)
      REFL(1,2) = G_REFL1(2)
      DIFU(1,2) = G_DIFU1(2)
      ABSP(1,2) = G_ABSP1(2)

C     +X FACE
      AREA(3,1) = G_AREA1(1)
      REFL(3,1) = G_REFL1(1)
      DIFU(3,1) = G_DIFU1(1)
      ABSP(3,1) = G_ABSP1(1)

C     +X CYLINDRICAL
      AREA2(3,1) = G_AREA2(1)
      REFL2(3,1) = G_REFL2(1)
      DIFU2(3,1) = G_DIFU2(1)
      ABSP2(3,1) = G_ABSP2(1)

C     -X FACE (ASSUMED FOR ALL BLOCKS)
      AREA(3,2) = G_AREA1(1)
      REFL(3,2) = G_REFL1(1)
      DIFU(3,2) = G_DIFU1(1)
      ABSP(3,2) = G_ABSP1(1)

C     -X CYLINDRICAL (ASSUMED FOR ALL BLOCKS)                                               
      AREA2(3,2) = G_AREA2(1)
      REFL2(3,2) = G_REFL2(1)
      DIFU2(3,2) = G_DIFU2(1)
      ABSP2(3,2) = G_ABSP2(1)

C     +Y FACE (ASSUMED FOR BLKNUM = 1...7)
      IF((BLKNUM.GE.1).AND.(BLKNUM.LE.7))THEN
C        DUE TO BLOCK IIR W-SENSOR
         IF((BLKNUM.GE.4).AND.(BLKNUM.LE.7))THEN      
            AREA(2,1) = ((AREA(1,1)-0.5D0) + AREA(3,1))/(2D0)
         ELSE
            AREA(2,1) = (AREA(1,1) + AREA(3,1))/(2D0)
         ENDIF
         REFL(2,1) = (REFL(1,1) + REFL(1,2) + REFL(3,1) + REFL(3,2))/4
         DIFU(2,1) = (DIFU(1,1) + DIFU(1,2) + DIFU(3,1) + DIFU(3,2))/4
         ABSP(2,1) = (ABSP(1,1) + ABSP(1,2) + ABSP(3,1) + ABSP(3,2))/4

C        CYLINDRICAL PART EQUAL TO +/- X CYLINDRICAL
         AREA2(2,1) = AREA2(3,1)
         REFL2(2,1) = REFL2(3,1)
         DIFU2(2,1) = DIFU2(3,1)
         ABSP2(2,1) = ABSP2(3,1)

C     FOR OTHER SATELLITE BLOCKS
      ELSE
         AREA(2,1) = G_AREA1(5)
         REFL(2,1) = G_REFL1(5)
         DIFU(2,1) = G_DIFU1(5)
         ABSP(2,1) = G_ABSP1(5)

         AREA2(2,1) = G_AREA2(5)
         REFL2(2,1) = G_REFL2(5)
         DIFU2(2,1) = G_DIFU2(5)
         ABSP2(2,1) = G_ABSP2(5)
      ENDIF

C     -Y FACE (ASSUMED FOR ALL BLOCKS)
      AREA(2,2) = AREA(2,1)
      REFL(2,2) = REFL(2,1)
      DIFU(2,2) = DIFU(2,1)
      ABSP(2,2) = ABSP(2,1)

C     -Y CYLINDRICAL (ASSUMED FOR ALL BLOCKS)
      AREA2(2,2) = AREA2(2,1)
      REFL2(2,2) = REFL2(2,1)
      DIFU2(2,2) = DIFU2(2,1)
      ABSP2(2,2) = ABSP2(2,1)

C     FRONT OF SOLAR PANELS
      AREA(4,1) = G_AREA1(4)
      REFL(4,1) = G_REFL1(4)
      DIFU(4,1) = G_DIFU1(4)
      ABSP(4,1) = G_ABSP1(4)

C     FRONT OF SOLAR PANELS (CYLINDER)
      AREA2(4,1) = G_AREA2(4)
      REFL2(4,1) = G_REFL2(4)
      DIFU2(4,1) = G_DIFU2(4)
      ABSP2(4,1) = G_ABSP2(4)

C     BACK OF SOLAR PANELS (ASSUMED FOR ALL BLOCKS)
      AREA(4,2) = AREA(4,1)
      REFL(4,2) = 0.055D0
      DIFU(4,2) = 0.055D0
      ABSP(4,2) = 0.890D0


C     BACK OF SOLAR PANELS (CYLINDER)
      AREA2(4,2) = G_AREA2(4)
      REFL2(4,2) = G_REFL2(4)
      DIFU2(4,2) = G_DIFU2(4)
      ABSP2(4,2) = G_ABSP2(4)


C     OPTICAL PROPERTIES IN THE INFRARED
C     ASSUMED FOR ALL BLOCKS AND ALL SURFACES
      DO SS=1,4
         DO II=1,2
            REFLIR(SS,II) = 0.10D0
            DIFUIR(SS,II) = 0.10D0
            ABSPIR(SS,II) = 0.80D0
         ENDDO
      ENDDO

C     NOT ENOUGH INFORMATION AVAILABLE FOR GLONASS-M
C     SEE IGSMAIL-5104
      IF(BLKNUM.EQ.102)THEN
         BUSFAC = 4.2D0/3.31D0
         DO SS=1,3
            DO II=1,2
               AREA(SS,II) = BUSFAC*AREA(SS,II)
               AREA2(SS,II) = BUSFAC*AREA2(SS,II)
            ENDDO
         ENDDO
         AREA(4,1) = 30.85D0
         AREA(4,2) = 30.85D0
      ENDIF

C     NO (YET) INFORMATION AVAILABLE FOR GLONASS-K
      IF(((BLKNUM.GE.9).AND.(BLKNUM.LE.100))
     1  .OR.(BLKNUM.GE.103))THEN
         DO SS=1,4
            DO II=1,2
               AREA(SS,II) = 0D0
               REFL(SS,II) = 0D0
               DIFU(SS,II) = 0D0
               ABSP(SS,II) = 0D0

               AREA2(SS,II) = 0D0
               REFL2(SS,II) = 0D0
               DIFU2(SS,II) = 0D0
               ABSP2(SS,II) = 0D0

               REFLIR(SS,II) = 0D0
               DIFUIR(SS,II) = 0D0
               ABSPIR(SS,II) = 0D0
            ENDDO
         ENDDO
      ENDIF

      END SUBROUTINE

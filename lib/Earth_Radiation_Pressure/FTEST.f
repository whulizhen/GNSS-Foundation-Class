      SUBROUTINE FTEST(ii,ff,DATPATH)
      integer ii
      real*4 ff
      CHARACTER*100 DATPATH
      write (6,100 ) ii,ff
 100  format('ii=',i2,'ff=',f6.3)
      write (*,*) DATPATH
      return 
      END SUBROUTINE

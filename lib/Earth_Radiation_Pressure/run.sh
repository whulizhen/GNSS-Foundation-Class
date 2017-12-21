rm *.o
rm *.a
gfortran -c erpTest.f90 ERPFBOXW.f PROPBOXW.f FTEST.f SURFBOXW.f
ar rcs liberp.a erpTest.o ERPFBOXW.o PROPBOXW.o FTEST.o SURFBOXW.o
g++ -L/usr/local/gfortran/lib -o main testC.cpp liberp.a -lgfortran
./main
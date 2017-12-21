# CRD Sample code
#   common_c

#---------------------------------------------------------------------------

WRTSRC= wr_test.c read_crd.c write_crd.c
WRTOBJ= wr_test.o read_crd.o write_crd.o

wr_test: $(WRTOBJ)
	gcc -g $(WRTOBJ) -lm -o wr_test

#---------------------------------------------------------------------------

all:	doytogr.o  getfield.o  grtodoy.o  grtojd.o  jdtogr.o read_crd.o  write_crd.o wr_test

doytogr.o: doytogr.c
	gcc -g -c -o doytogr.o doytogr.c

getfield.o: getfield.c
	gcc -g -c -o getfield.o getfield.c

grtodoy.o: grtodoy.c
	gcc -g -c -o grtodoy.o grtodoy.c

grtojd.o: grtojd.c
	gcc -g -c -o grtojd.o grtojd.c

jdtogr.o: jdtogr.c
	gcc -g -c -o jdtogr.o jdtogr.c

read_crd.o: read_crd.c
	gcc -g -c -o read_crd.o read_crd.c

write_crd.o: write_crd.c
	gcc -g -c -o write_crd.o write_crd.c

clean:
	rm -f *.o


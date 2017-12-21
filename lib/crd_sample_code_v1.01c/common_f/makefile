# CRD Sample code
#   common_f

#---------------------------------------------------------------------------

WRTSRC= wr_testf.c
WRTOBJ= wr_testf.o

wr_testf: $(WRTOBJ)
	f77 -g $(WRTOBJ) \
        read_crdf.o write_crdf.o trimlen.o -lm \
        -o wr_testf

#---------------------------------------------------------------------------
all: read_crdf.o  write_crdf.o trimlen.o wr_testf

read_crdf.o: read_crdf.f
	f77 -g -c -o read_crdf.o read_crdf.f

write_crdf.o: write_crdf.f
	f77 -g -c -o write_crdf.o write_crdf.f

trimlen.o: trimlen.f
	f77 -g -c -o trimlen.o trimlen.f

clean:
	rm -f *.o

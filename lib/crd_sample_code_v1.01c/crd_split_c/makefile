# CRD Sample code
#   crd_split_c

INCLUDES = -I../include -I.

# This variable lists the required libraries.
LIBS=   -lm
CFLAGS = -g -I/usr/include -I. $(INCLUDES)
CC= gcc

.o:     .c

#---------------------------------------------------------------------------

FRSSRC= frd_strip.c

FRSOBJ= frd_strip.o

frs:	 $(FRSOBJ)
	$(CC) $(CFLAGS) $(FRSOBJ) ../common_c/read_crd.o ../common_c/write_crd.o \
	-lm -o frd_strip
#---------------------------------------------------------------------------

CRSSRC= crd_split.c

CRSOBJ= crd_split.o

crs:	 $(CRSOBJ)
	$(CC) $(CFLAGS) $(CRSOBJ) ../common_c/read_crd.o ../common_c/write_crd.o \
	-lm -o crd_split
#---------------------------------------------------------------------------
MCDSRC= merge_crd_daily.c

MCDOBJ= merge_crd_daily.o

mcd:	$(MCDOBJ)
	$(CC) $(CFLAGS) $(MCDOBJ) \
	../common_c/read_crd.o ../common_c/write_crd.o ../common_c/grtojd.o \
	../common_c/doytogr.o ../common_c/jdtogr.o ../common_c/getfield.o \
	-lm -o merge_crd_daily
#---------------------------------------------------------------------------

clean:
	rm -f *.o frd_strip crd_split merge_crd_daily

all:
	make frs
	make crs
	make mcd

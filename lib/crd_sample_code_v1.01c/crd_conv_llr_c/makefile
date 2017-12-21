# CRD Sample code
#   crd_conv_llr_c

INCLUDES = -I../include -I.

# This variable lists the required libraries.
LIBS=   -lm
CFLAGS = -g -I/usr/include -I. $(INCLUDES)
CC= gcc

.o:     .c

#---------------------------------------------------------------------------
LTCSRC= cllr_to_crd.c read_llr.c

LTCOBJ= cllr_to_crd.o read_llr.o

ltc:	 $(LTCOBJ)
	$(CC) $(CFLAGS) $(LTCOBJ) \
	../common_c/write_crd.o ../common_c/doytogr.o ../common_c/grtojd.o \
	../common_c/getfield.o -lm \
	-o cllr_to_crd

#---------------------------------------------------------------------------

clean:
	rm -f *.o cllr_to_crd

all:
	make ltc

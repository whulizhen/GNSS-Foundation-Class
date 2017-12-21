# CRD Sample code
#   crd_chk_c

INCLUDES = -I../include -I.
 
# This variable lists the required libraries.
LIBS=   -lm
CFLAGS = -g -I/usr/include -I. $(INCLUDES) 
CC= gcc

.o:	.c

#---------------------------------------------------------------------------

CCHSRC= crd_chk.c
CCHOBJ= crd_chk.o

crd_chk:	 $(CCHOBJ)
	$(CC) $(CFLAGS) $(CCHOBJ) \
	../common_c/read_crd.o ../common_c/grtodoy.o ../common_c/grtojd.o \
	../common_c/getfield.o ${LIBS} -o crd_chk
#---------------------------------------------------------------------------

CCNSRC= crd_cstg_np_cmp.c
CCNOBJ= crd_cstg_np_cmp.o

crd_cstg_np_cmp:     $(CCNOBJ)
	$(CC) $(CFLAGS) $(CCNOBJ) \
                ../common_c/read_crd.o ../common_c/grtodoy.o \
		../common_c/getfield.o \
		../crd_conv_slr_c/read_cstg.o ${LIBS} \
                -o crd_cstg_np_cmp
#---------------------------------------------------------------------------

CRFSRC= crd_merit_fr_cmp.c
CRFOBJ= crd_merit_fr_cmp.o

crd_merit_fr_cmp:     $(CRFOBJ)
	$(CC) $(CFLAGS) $(CRFOBJ) \
        ../common_c/read_crd.o ../common_c/grtodoy.o \
        ../common_c/getfield.o \
        ../crd_conv_slr_c/read_merit.o ${LIBS} \
        -o crd_merit_fr_cmp

#---------------------------------------------------------------------------

clean:
	rm -f *.o crd_chk crd_cstg_np_cmp crd_merit_fr_cmp

all:
	make crd_chk
	make crd_cstg_np_cmp
	make crd_merit_fr_cmp

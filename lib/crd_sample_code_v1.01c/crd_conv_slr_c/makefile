# CRD Sample code
#   crd_conv_slr_c

INCLUDES = -I../include -I.

# This variable lists the required libraries.
LIBS=   -lm
CFLAGS = -g -I/usr/include -I. $(INCLUDES)
CC= gcc

.o:     .c

#---------------------------------------------------------------------------
CTCSRC= cstg_to_crd.c read_cstg.c

CTCOBJ= cstg_to_crd.o read_cstg.o

ctc:	 $(CTCOBJ)
	$(CC) $(CFLAGS) $(CTCOBJ) \
	../common_c/write_crd.o ../common_c/doytogr.o ../common_c/grtojd.o \
	../common_c/getfield.o -lm \
	-o cstg_to_crd

#---------------------------------------------------------------------------
MTCSRC= merit_to_crd.c read_merit.c

MTCOBJ= merit_to_crd.o read_merit.o

mtc:	 $(MTCOBJ)
	$(CC) $(CFLAGS) $(MTCOBJ) \
	../common_c/write_crd.o ../common_c/doytogr.o ../common_c/grtojd.o \
	../common_c/getfield.o -lm \
	-o merit_to_crd

#---------------------------------------------------------------------------

CFCNSRC= crd_to_cstg_np.c write_cstg.c

CFCNOBJ= crd_to_cstg_np.o write_cstg.o

cfcn:	 $(CFCNOBJ)
	$(CC) $(CFLAGS) $(CFCNOBJ) \
	../common_c/read_crd.o ../common_c/grtodoy.o ../common_c/grtojd.o \
	../common_c/getfield.o -lm \
	-o crd_to_cstg_np
#---------------------------------------------------------------------------

CFCQSRC= crd_to_cstg_ql.c write_cstg.c

CFCQOBJ= crd_to_cstg_ql.o write_cstg.o

cfcq:	 $(CFCQOBJ)
	$(CC) $(CFLAGS) $(CFCQOBJ) \
	../common_c/read_crd.o ../common_c/grtodoy.o ../common_c/grtojd.o \
	../common_c/getfield.o -lm \
	-o crd_to_cstg_ql
#---------------------------------------------------------------------------

CFMSRC= crd_to_merit.c write_merit.c

CFMOBJ= crd_to_merit.o write_merit.o

cfm:	 $(CFMOBJ)
	$(CC) $(CFLAGS) $(CFMOBJ) \
	../common_c/read_crd.o ../common_c/grtodoy.o ../common_c/grtojd.o\
        ../common_c/getfield.o -lm \
	-o crd_to_merit
#---------------------------------------------------------------------------

clean:
	rm -f *.o cstg_to_crd merit_to_crd crd_to_cstg_np \
	crd_to_cstg_ql crd_to_merit

all:
	make ctc
	make mtc
	make cfcn
	make cfcq
	make cfm


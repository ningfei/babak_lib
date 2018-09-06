include ../Makefile

PROJNAME = niftiio

INCFLAGS = $(ZLIB_INC) $(ZNZ_INC)
LIBS = $(ZLIB_LIBS) $(ZNZ_LIBS) $(NIFTI_LIBS)

SRCS=nifti1_io.c
OBJS=nifti1_io.o

SCRIPTS =

depend:
	rm -f depend.mk
	make depend.mk

depend.mk:
	$(CC) $(DEPENDFLAGS) $(INCFLAGS) $(SRCS) >> depend.mk

lib: libniftiio.a

nifti1_io.o: nifti1_io.c nifti1_io.h
	$(CC) -c $(CFLAGS) $(USEZLIB) $(INCFLAGS) $<

libniftiio.a: $(OBJS)
	$(AR) -r libniftiio.a $(OBJS)
	$(RANLIB) $@

include depend.mk

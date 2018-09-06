## Projects
ZNZ		=	znzlib
NIFTI   =	niftilib

## Compiler  defines
CC		=	gcc
AR		=	ar
RANLIB  = 	ranlib
DEPENDFLAGS	=	-MM
GNU_ANSI_FLAGS	= 	-Wall -ansi -pedantic
ANSI_FLAGS	= 	${GNU_ANSI_FLAGS}
CFLAGS		=	$(ANSI_FLAGS)

MACH = $(shell $(CC) -dumpmachine | awk -F - '{ printf "%s", $$2; }')

## Installation
INSTALL_LIB_DIR	=	../$(MACH)

## Zlib defines
ZLIB_INC	=
ZLIB_PATH	=
ZLIB_LIBS 	= 	$(ZLIB_PATH) -lm -lz

## ZNZ defines
ZNZ_INC		=	-I../$(ZNZ)
ZNZ_PATH	=	-L../$(ZNZ)
ZNZ_LIBS	=	$(ZNZ_PATH) -lznz
USEZLIB	    =   -DHAVE_ZLIB

## NIFTI defines
NIFTI_INC	=	-I../$(NIFTI)
NIFTI_PATH	=	-L../$(NIFTI)
NIFTI_LIBS	=	$(NIFTI_PATH) -lniftiio

## Rules

.SUFFIXES: .c .o

.c.o:
	$(CC) -c $(CFLAGS) $(INCFLAGS) $<

## Targets

all:	   znz nifti install

install:   znz_install nifti_install

clean:	   znz_clean nifti_clean

clean_all: clean install_clean

znz:
	(cd $(ZNZ); make lib;)
	@echo ""

nifti: znz
	(cd $(NIFTI); make depend; make lib;)
	@echo ""

znz_install: $(INSTALL_LIB_DIR)
	cp $(ZNZ)/libznz.a $(INSTALL_LIB_DIR)
	$(RANLIB) $(INSTALL_LIB_DIR)/libznz.a
	@echo ""

nifti_install: $(INSTALL_LIB_DIR)
	cp $(NIFTI)/libniftiio.a $(INSTALL_LIB_DIR)
	$(RANLIB) $(INSTALL_LIB_DIR)/libniftiio.a
	@echo ""

install_clean:
	rm -f $(INSTALL_LIB_DIR)/libznz.a $(INSTALL_LIB_DIR)/libniftiio.a

znz_clean:
	(cd $(ZNZ); rm -f *.mk *.o *.a core;)

nifti_clean:
	(cd $(NIFTI); rm -f *.mk *.o *.a core;)

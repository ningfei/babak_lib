## Projects
Z 	 	= 	zlib
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

## Installation
INSTALL_LIB_DIR	=	..

## Zlib defines
ZLIB_INC	=	-I../$(Z)
ZLIB_PATH	=	-L../$(Z)
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

all:	   z znz nifti install

install:   z_install znz_install nifti_install

clean:	   z_clean znz_clean nifti_clean

clean_all: clean install_clean

z:
	(cd $(Z); ./configure; make libz.a;)
	@echo ""

znz: z
	(cd $(ZNZ); make lib;)
	@echo ""

nifti: znz
	(cd $(NIFTI); make depend; make lib;)
	@echo ""

z_install: $(INSTALL_LIB_DIR)
	cp $(Z)/libz.a $(INSTALL_LIB_DIR)
	$(RANLIB) $(INSTALL_LIB_DIR)/libz.a
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

z_clean:
	(cd $(Z); rm -f *.log *.o *.a core;)

znz_clean:
	(cd $(ZNZ); rm -f *.mk *.o *.a core;)

nifti_clean:
	(cd $(NIFTI); rm -f *.mk *.o *.a core;)

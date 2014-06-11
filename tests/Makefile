#LIBS = -L$(HOME)/lib -ltmglib_LINUX -llapack_LINUX -lblas_LINUX -lI77 -lF77
LIBS = -L$(HOME)/lib  -llapack_LINUX -lblas_LINUX -lI77 -lF77
CLIBS = -L$(HOME)/lib -lbabak_lib_linux -lm -lc
INC = -I$(HOME)/include

CC=g++
CFLAGS= -funroll-all-loops -O3

all: statistical_test removeVectorMean_test normalizeVector_test

statistical_test: statistical_test.c
	$(CC) $(CFLAGS) -o statistical_test statistical_test.c $(INC) -ldcdf $(LIBS) $(CLIBS)

removeVectorMean_test: removeVectorMean_test.cxx
	$(CC) $(CFLAGS) -o removeVectorMean_test removeVectorMean_test.cxx $(INC) -ldcdf $(LIBS) $(CLIBS)

normalizeVector_test: normalizeVector_test.cxx
	$(CC) $(CFLAGS) -o normalizeVector_test normalizeVector_test.cxx $(INC) -ldcdf $(LIBS) $(CLIBS)

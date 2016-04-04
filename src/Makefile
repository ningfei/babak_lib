CFLAGS = -funroll-all-loops -O3
CC = g++
#INC = -I$(HOME)/include -I/usr/local/dmp/nifti/include
INC= -I$(HOME)/include -I/usr/local/dmp/include -I/usr/local/dmp/nifti/include
#LIBS = -L$(HOME)/lib -L/usr/local/dmp/nifti/lib -lbabak_lib_linux -lniftiio -lznz -lm -lz -lc
LIBS = -L$(HOME)/lib -L/usr/local/dmp/nifti/lib -L ../lib -lbabak_lib_linux -lniftiio -lznz -lm -lz -lc -L/usr/local/dmp/lib -ldcdf -llevmar -llapack -lblas -lf2c

all: libbabak_lib_linux.a 3dwarper acpcdetect applywarp3d scaleImage avgImage ivf unwarp2d cphdr reslice kaiba oi

# because of functions in matrixops, linking the libbabak_lib will possibly require lapack libs.

set_dim.o: set_dim.c
	$(CC) -O2 -c $(INC) set_dim.c

findMSP.o: findMSP.c
	$(CC) -O2 -c $(INC) findMSP.c

singular_value_decomposition.o: singular_value_decomposition.c
	$(CC) -O2 -c $(INC) singular_value_decomposition.c

DKI.o: DKI.c
	$(CC) -O2 -c $(INC) DKI.c

artlib.o: artlib.c
	$(CC) -O2 -c $(INC) artlib.c

volume.o: volume.c
	$(CC) -O2 -c $(INC) volume.c

binomial.o: binomial.c
	$(CC) -O2 -c $(INC) binomial.c

smooth.o: smooth.c
	$(CC) -O2 -c $(INC) smooth.c

niftiimage.o: niftiimage.c
	$(CC) -O2 -c $(INC) niftiimage.c

utils.o: utils.c
	$(CC) -O2 -c $(INC) utils.c

nifti.o: nifti.c
	$(CC) -O2 -c $(INC) nifti.c

errorMessage.o: errorMessage.c
	$(CC) -O2 -c $(INC) errorMessage.c

maskOps.o: maskOps.c
	$(CC) -O2 -c $(INC) maskOps.c

EMFIT.o: EMFIT.c
	$(CC) -O2 -c $(INC) EMFIT.c

max_cc.o: max_cc.c
	$(CC) -O2 -c $(INC) max_cc.c

statistics.o: statistics.c
	$(CC) -O2 -c $(INC) statistics.c

ginverse.o: ginverse.c
	$(CC) -O2 -c $(INC) ginverse.c

permutation.o: permutation.c
	$(CC) -O2 -c $(INC) permutation.c

hpsort.o: hpsort.c
	$(CC) -O2 -c $(INC) hpsort.c

random.o: random.c
	$(CC) -O2 -c $(INC) random.c

dicomIO.o: dicomIO.c
	$(CC) -O2 -c $(INC) dicomIO.c

nkiIO.o: nkiIO.c
	$(CC) -O2 -c $(INC) nkiIO.c

matrixops.o: matrixops.c
	$(CC) -O2 -c $(INC) matrixops.c

subsets.o: subsets.c 
	$(CC) -O2 -c $(INC) subsets.c 

analyzeio.o: analyzeio.c 
	$(CC) -O2 -c $(INC) analyzeio.c 

getoption.o: getoption.c 
	$(CC) -O2 -c $(INC) getoption.c 

swap.o: swap.c 
	$(CC) -O2 -c $(INC) swap.c 

fileinfo.o: fileinfo.c
	$(CC) -O2 -c $(INC) fileinfo.c

histogram.o: histogram.c
	$(CC) -O2 -c $(INC) histogram.c

resize.o: resize.c
	$(CC) -O2 -c $(INC) resize.c

gaussian_kernel.o: gaussian_kernel.c
	$(CC) -O2 -c $(INC) gaussian_kernel.c

convolution.o: convolution.c
	$(CC) -O2 -c $(INC) convolution.c

resliceImage.o: resliceImage.c
	$(CC) -O2 -c $(INC) resliceImage.c

matrixCom.o: matrixCom.c
	$(CC) -O2 -c $(INC) matrixCom.c

registration.o: registration.c
	$(CC) -O2 -c $(INC) registration.c

legendre.o: legendre.c
	$(CC) -O2 -c $(INC) legendre.c

cubicspline.o: cubicspline.c
	$(CC) -O2 -c $(INC) cubicspline.c

medianfilter.o: medianfilter.c
	$(CC) -O2 -c $(INC) medianfilter.c

directionCode.o: directionCode.c
	$(CC) -O2 -c $(INC) directionCode.c

checkNiftiFileExtension.o: checkNiftiFileExtension.c
	$(CC) -O2 -c $(INC) checkNiftiFileExtension.c

getNiftiImageOrientation.o: getNiftiImageOrientation.c
	$(CC) -O2 -c $(INC) getNiftiImageOrientation.c

isOrientationCodeValid.o: isOrientationCodeValid.c
	$(CC) -O2 -c $(INC) isOrientationCodeValid.c

PILtransform.o: PILtransform.c
	$(CC) -O2 -c $(INC) PILtransform.c

reorientVolume.o: reorientVolume.c
	$(CC) -O2 -c $(INC) reorientVolume.c

rotate.o: rotate.c
	$(CC) -O2 -c $(INC) rotate.c

setLowHigh.o: setLowHigh.c
	$(CC) -O2 -c $(INC) setLowHigh.c

compute_cm.o: compute_cm.c
	$(CC) -O2 -c $(INC) compute_cm.c

standardize.o: standardize.c
	$(CC) -O2 -c $(INC) standardize.c

sph.o: sph.c
	$(CC) -O2 -c $(INC) sph.c

landmarks.o: landmarks.c
	$(CC) -O2 -c $(INC) landmarks.c

getARTHOME.o: getARTHOME.c
	$(CC) -O2 -c $(INC) getARTHOME.c

libbabak_lib_linux.a: set_dim.o findMSP.o singular_value_decomposition.o DKI.o artlib.o volume.o binomial.o smooth.o niftiimage.o utils.o nifti.o errorMessage.o maskOps.o EMFIT.o max_cc.o statistics.o ginverse.o permutation.o hpsort.o random.o dicomIO.o nkiIO.o subsets.o analyzeio.o getoption.o swap.o fileinfo.o histogram.o resize.o gaussian_kernel.o convolution.o resliceImage.o matrixCom.o registration.o legendre.o cubicspline.o medianfilter.o directionCode.o checkNiftiFileExtension.o getNiftiImageOrientation.o isOrientationCodeValid.o PILtransform.o reorientVolume.o rotate.o setLowHigh.o compute_cm.o standardize.o sph.o matrixops.o landmarks.o getARTHOME.o
	ar -ru ../lib/libbabak_lib_linux.a set_dim.o findMSP.o singular_value_decomposition.o DKI.o artlib.o volume.o binomial.o smooth.o niftiimage.o utils.o nifti.o errorMessage.o maskOps.o EMFIT.o max_cc.o statistics.o ginverse.o permutation.o hpsort.o random.o dicomIO.o nkiIO.o subsets.o analyzeio.o getoption.o swap.o fileinfo.o histogram.o resize.o gaussian_kernel.o convolution.o resliceImage.o matrixCom.o registration.o legendre.o cubicspline.o medianfilter.o directionCode.o checkNiftiFileExtension.o getNiftiImageOrientation.o isOrientationCodeValid.o PILtransform.o reorientVolume.o rotate.o setLowHigh.o compute_cm.o standardize.o sph.o matrixops.o landmarks.o getARTHOME.o

clean: 
	rm -f *.o

# compilation of executables
3dwarper: libbabak_lib_linux.a 3dwarper.cxx
	$(CC) $(CFLAGS) -o 3dwarper 3dwarper.cxx $(LIBS) $(INC)
	mv 3dwarper ../bin

acpcdetect: libbabak_lib_linux.a acpcdetect.cxx
	$(CC) $(CFLAGS) -o acpcdetect acpcdetect.cxx $(LIBS) $(INC)
	mv acpcdetect ../bin

applywarp3d: libbabak_lib_linux.a applywarp3d.cxx
	$(CC) $(CFLAGS) -o applywarp3d applywarp3d.cxx $(INC) $(LIBS)
	mv applywarp3d ../bin

scaleImage: libbabak_lib_linux.a scaleImage.cxx
	$(CC) $(CFLAGS) -o scaleImage scaleImage.cxx $(INC) $(LIBS)
	mv scaleImage ../bin

avgImage: libbabak_lib_linux.a avgImage.cxx
	$(CC) $(CFLAGS) -o avgImage avgImage.cxx $(INC) $(LIBS)
	mv avgImage ../bin

ivf: libbabak_lib_linux.a ivf.cxx
	$(CC) $(CFLAGS) -o ivf ivf.cxx $(INC) $(LIBS)
	mv ivf ../bin

unwarp2d: libbabak_lib_linux.a unwarp2d.cxx
	$(CC) $(CFLAGS) -o unwarp2d unwarp2d.cxx $(INC) $(LIBS)
	mv unwarp2d ../bin

cphdr: libbabak_lib_linux.a cphdr.cxx
	$(CC) $(CFLAGS) -o cphdr cphdr.cxx $(INC) $(LIBS)
	mv cphdr ../bin

reslice: libbabak_lib_linux.a reslice.cxx
	$(CC) $(CFLAGS) -o reslice reslice.cxx $(INC) $(LIBS)
	mv reslice ../bin

kaiba: libbabak_lib_linux.a kaiba.cxx
	$(CC) $(CFLAGS) -o kaiba kaiba.cxx $(INC) $(LIBS)
	mv kaiba ../bin

oi: libbabak_lib_linux.a oi.cxx
	$(CC) $(CFLAGS) -o oi oi.cxx $(INC) $(LIBS)
	mv oi ../bin

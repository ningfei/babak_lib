# for makhmal
MPICC = /usr/local/bin/mpic++ -Wall

CXX = g++
CXXFLAGS = -O3

MACH = $(shell $(CXX) -dumpmachine | awk -F - '{ printf "%s", $$2; }')

ifeq ($(MACH), linux)
	SYSTEM = linux
	CXXFLAGS +=  -Wno-unused-result -Wno-format-zero-length
	CXXSTATICFLAGS =  -static
else ifeq ($(MACH), apple)
	SYSTEM = mac
	CXXSTATICFLAGS =
else ifeq ($(MACH), w64)
	SYSTEM = win
	CXXSTATICFLAGS =  -static
endif

INC= -I../include
LIBS =  -L../lib/$(SYSTEM) -lbabak -lniftiio -lznz -lz -lm -llapack -lblas -lgfortran -lquadmath
LESSLIBS =  -L../lib/$(SYSTEM) -lbabak -lniftiio -lznz -lz -lm

ifeq ($(MACH), apple)
	LIBS += -lgcc
endif

all: libbabak.a brainwash 3dwarper acpcdetect applywarp3d scaleImage avgImage one_sample_t ivf unwarp2d cphdr reslice kaiba oi PILtransform atra art2fsl fsl2art yuki analyze2nii multiplyTransformation niftiinfo break4D gambling_vancova_setup medianImage threshold ccstats

# because of functions in matrixops, linking the libbabak will possibly require lapack libs.

leastSquaresAffineTrans.o: leastSquaresAffineTrans.cpp
	$(CXX) $(CXXFLAGS) -c $(INC) leastSquaresAffineTrans.cpp

hist1D_plot.o: hist1D_plot.cpp
	$(CXX) $(CXXFLAGS) -c $(INC) hist1D_plot.cpp

set_dim.o: set_dim.cpp
	$(CXX) $(CXXFLAGS) -c $(INC) set_dim.cpp

findMSP.o: findMSP.cpp
	$(CXX) $(CXXFLAGS) -c $(INC) findMSP.cpp

singular_value_decomposition.o: singular_value_decomposition.cpp
	$(CXX) $(CXXFLAGS) -c $(INC) singular_value_decomposition.cpp

DKI.o: DKI.cpp
	$(CXX) $(CXXFLAGS) -c $(INC) DKI.cpp

artlib.o: artlib.cpp
	$(CXX) $(CXXFLAGS) -c $(INC) artlib.cpp

volume.o: volume.cpp
	$(CXX) $(CXXFLAGS) -c $(INC) volume.cpp

binomial.o: binomial.cpp
	$(CXX) $(CXXFLAGS) -c $(INC) binomial.cpp

smooth.o: smooth.cpp
	$(CXX) $(CXXFLAGS) -c $(INC) smooth.cpp

niftiimage.o: niftiimage.cpp
	$(CXX) $(CXXFLAGS) -c $(INC) niftiimage.cpp

utils.o: utils.cpp
	$(CXX) $(CXXFLAGS) -c $(INC) utils.cpp

nifti.o: nifti.cpp
	$(CXX) $(CXXFLAGS) -c $(INC) nifti.cpp

errorMessage.o: errorMessage.cpp
	$(CXX) $(CXXFLAGS) -c $(INC) errorMessage.cpp

maskOps.o: maskOps.cpp
	$(CXX) $(CXXFLAGS) -c $(INC) maskOps.cpp

EMFIT.o: EMFIT.cpp
	$(CXX) $(CXXFLAGS) -c $(INC) EMFIT.cpp

max_cc.o: max_cc.cpp
	$(CXX) $(CXXFLAGS) -c $(INC) max_cc.cpp

statistics.o: statistics.cpp
	$(CXX) $(CXXFLAGS) -c $(INC) statistics.cpp

ginverse.o: ginverse.cpp
	$(CXX) $(CXXFLAGS) -c $(INC) ginverse.cpp

permutation.o: permutation.cpp
	$(CXX) $(CXXFLAGS) -c $(INC) permutation.cpp

hpsort.o: hpsort.cpp
	$(CXX) $(CXXFLAGS) -c $(INC) hpsort.cpp

random.o: random.cpp
	$(CXX) $(CXXFLAGS) -c $(INC) random.cpp

dicomIO.o: dicomIO.cpp
	$(CXX) $(CXXFLAGS) -c $(INC) dicomIO.cpp

nkiIO.o: nkiIO.cpp
	$(CXX) $(CXXFLAGS) -c $(INC) nkiIO.cpp

matrixops.o: matrixops.cpp
	$(CXX) $(CXXFLAGS) -c $(INC) matrixops.cpp

subsets.o: subsets.cpp
	$(CXX) $(CXXFLAGS) -c $(INC) subsets.cpp

analyzeio.o: analyzeio.cpp
	$(CXX) $(CXXFLAGS) -c $(INC) analyzeio.cpp

getoption.o: getoption.cpp
	$(CXX) $(CXXFLAGS) -c $(INC) getoption.cpp

swap.o: swap.cpp
	$(CXX) $(CXXFLAGS) -c $(INC) swap.cpp

fileinfo.o: fileinfo.cpp
	$(CXX) $(CXXFLAGS) -c $(INC) fileinfo.cpp

histogram.o: histogram.cpp
	$(CXX) $(CXXFLAGS) -c $(INC) histogram.cpp

resize.o: resize.cpp
	$(CXX) $(CXXFLAGS) -c $(INC) resize.cpp

gaussian_kernel.o: gaussian_kernel.cpp
	$(CXX) $(CXXFLAGS) -c $(INC) gaussian_kernel.cpp

convolution.o: convolution.cpp
	$(CXX) $(CXXFLAGS) -c $(INC) convolution.cpp

resliceImage.o: resliceImage.cpp
	$(CXX) $(CXXFLAGS) -c $(INC) resliceImage.cpp

matrixCom.o: matrixCom.cpp
	$(CXX) $(CXXFLAGS) -c $(INC) matrixCom.cpp

registration.o: registration.cpp
	$(CXX) $(CXXFLAGS) -c $(INC) registration.cpp

legendre.o: legendre.cpp
	$(CXX) $(CXXFLAGS) -c $(INC) legendre.cpp

cubicspline.o: cubicspline.cpp
	$(CXX) $(CXXFLAGS) -c $(INC) cubicspline.cpp

medianfilter.o: medianfilter.cpp
	$(CXX) $(CXXFLAGS) -c $(INC) medianfilter.cpp

directionCode.o: directionCode.cpp
	$(CXX) $(CXXFLAGS) -c $(INC) directionCode.cpp

checkNiftiFileExtension.o: checkNiftiFileExtension.cpp
	$(CXX) $(CXXFLAGS) -c $(INC) checkNiftiFileExtension.cpp

getNiftiImageOrientation.o: getNiftiImageOrientation.cpp
	$(CXX) $(CXXFLAGS) -c $(INC) getNiftiImageOrientation.cpp

isOrientationCodeValid.o: isOrientationCodeValid.cpp
	$(CXX) $(CXXFLAGS) -c $(INC) isOrientationCodeValid.cpp

PILtransform.o: PILtransform.cpp
	$(CXX) $(CXXFLAGS) -c $(INC) PILtransform.cpp

reorientVolume.o: reorientVolume.cpp
	$(CXX) $(CXXFLAGS) -c $(INC) reorientVolume.cpp

rotate.o: rotate.cpp
	$(CXX) $(CXXFLAGS) -c $(INC) rotate.cpp

setLowHigh.o: setLowHigh.cpp
	$(CXX) $(CXXFLAGS) -c $(INC) setLowHigh.cpp

compute_cm.o: compute_cm.cpp
	$(CXX) $(CXXFLAGS) -c $(INC) compute_cm.cpp

standardize.o: standardize.cpp
	$(CXX) $(CXXFLAGS) -c $(INC) standardize.cpp

sph.o: sph.cpp
	$(CXX) $(CXXFLAGS) -c $(INC) sph.cpp

landmarks.o: landmarks.cpp
	$(CXX) $(CXXFLAGS) -c $(INC) landmarks.cpp

getARTHOME.o: getARTHOME.cpp
	$(CXX) $(CXXFLAGS) -c $(INC) getARTHOME.cpp

update_qsform.o: update_qsform.cpp
	$(CXX) $(CXXFLAGS) -c $(INC) update_qsform.cpp

libbabak.a: leastSquaresAffineTrans.o hist1D_plot.o set_dim.o findMSP.o singular_value_decomposition.o DKI.o artlib.o volume.o binomial.o smooth.o niftiimage.o utils.o nifti.o errorMessage.o maskOps.o EMFIT.o max_cc.o statistics.o ginverse.o permutation.o hpsort.o random.o dicomIO.o nkiIO.o subsets.o analyzeio.o getoption.o swap.o fileinfo.o histogram.o resize.o gaussian_kernel.o convolution.o resliceImage.o matrixCom.o registration.o legendre.o cubicspline.o medianfilter.o directionCode.o checkNiftiFileExtension.o getNiftiImageOrientation.o isOrientationCodeValid.o PILtransform.o reorientVolume.o rotate.o setLowHigh.o compute_cm.o standardize.o sph.o matrixops.o landmarks.o getARTHOME.o update_qsform.o
	ar -cru ../lib/$(SYSTEM)/libbabak.a leastSquaresAffineTrans.o hist1D_plot.o set_dim.o findMSP.o singular_value_decomposition.o DKI.o artlib.o volume.o binomial.o smooth.o niftiimage.o utils.o nifti.o errorMessage.o maskOps.o EMFIT.o max_cc.o statistics.o ginverse.o permutation.o hpsort.o random.o dicomIO.o nkiIO.o subsets.o analyzeio.o getoption.o swap.o fileinfo.o histogram.o resize.o gaussian_kernel.o convolution.o resliceImage.o matrixCom.o registration.o legendre.o cubicspline.o medianfilter.o directionCode.o checkNiftiFileExtension.o getNiftiImageOrientation.o isOrientationCodeValid.o PILtransform.o reorientVolume.o rotate.o setLowHigh.o compute_cm.o standardize.o sph.o matrixops.o landmarks.o getARTHOME.o update_qsform.o
	ranlib ../lib/$(SYSTEM)/libbabak.a

clean:
	rm -f *.o

# compilation of executables
PILtransform: PILtransform.cxx
	$(CXX) $(CXXSTATICFLAGS) $(CXXFLAGS) -o PILtransform PILtransform.cxx $(INC) $(LIBS)
	mv PILtransform ${ARTHOME}/bin/$(SYSTEM)

brainwash: brainwash.cxx
	$(CXX) $(CXXSTATICFLAGS) $(CXXFLAGS) -o brainwash brainwash.cxx $(LIBS) $(INC)
	mv brainwash ${ARTHOME}/bin/$(SYSTEM)

3dwarper: 3dwarper.cxx
	$(CXX) $(CXXSTATICFLAGS) $(CXXFLAGS) -o 3dwarper 3dwarper.cxx $(LIBS) $(INC)
	mv 3dwarper ${ARTHOME}/bin/$(SYSTEM)

acpcdetect: acpcdetect.cxx
	$(CXX) $(CXXSTATICFLAGS) $(CXXFLAGS) -o acpcdetect acpcdetect.cxx $(LIBS) $(INC)
	mv acpcdetect ${ARTHOME}/bin/$(SYSTEM)

applywarp3d: applywarp3d.cxx
	$(CXX) $(CXXSTATICFLAGS) $(CXXFLAGS) -o applywarp3d applywarp3d.cxx $(INC) $(LESSLIBS)
	mv applywarp3d ${ARTHOME}/bin/$(SYSTEM)

scaleImage: scaleImage.cxx
	$(CXX) $(CXXSTATICFLAGS) $(CXXFLAGS) -o scaleImage scaleImage.cxx $(INC) $(LIBS)
	mv scaleImage ${ARTHOME}/bin/$(SYSTEM)

avgImage: avgImage.cxx
	$(CXX) $(CXXSTATICFLAGS) $(CXXFLAGS) -o avgImage avgImage.cxx $(INC) $(LESSLIBS)
	mv avgImage ${ARTHOME}/bin/$(SYSTEM)

one_sample_t: one_sample_t.cxx
	$(CXX) $(CXXSTATICFLAGS) $(CXXFLAGS) -o one_sample_t one_sample_t.cxx $(INC) $(LIBS)
	mv one_sample_t ${ARTHOME}/bin/$(SYSTEM)

ivf: ivf.cxx
	$(CXX) $(CXXSTATICFLAGS) $(CXXFLAGS) -o ivf ivf.cxx $(INC) $(LIBS)
	mv ivf ${ARTHOME}/bin/$(SYSTEM)

unwarp2d: unwarp2d.cxx
	$(CXX) $(CXXSTATICFLAGS) $(CXXFLAGS) -o unwarp2d unwarp2d.cxx $(INC) $(LIBS)
	mv unwarp2d ${ARTHOME}/bin/$(SYSTEM)

cphdr: cphdr.cxx
	$(CXX) $(CXXSTATICFLAGS) $(CXXFLAGS) -o cphdr cphdr.cxx $(INC) $(LESSLIBS)
	mv cphdr ${ARTHOME}/bin/$(SYSTEM)

reslice: reslice.cxx
	$(CXX) $(CXXSTATICFLAGS) $(CXXFLAGS) -o reslice reslice.cxx $(INC) $(LIBS)
	mv reslice ${ARTHOME}/bin/$(SYSTEM)

kaiba: kaiba.cxx
	$(CXX) $(CXXSTATICFLAGS) $(CXXFLAGS) -o kaiba kaiba.cxx $(INC) $(LIBS)
	mv kaiba ${ARTHOME}/bin/$(SYSTEM)

#Example: $(CXX) $(CXXSTATICFLAGS) $(CXXFLAGS) -o kaiba kaiba.cxx -DBETA_PARAM=0.3 $(INC) $(LIBS)
# If we used the option -DBETA_PARAM=0.3 as in the above example, it will be the same as
# having the following statement in kaiba.cxx
#define BETA_PARAM 0.3

oi: oi.cxx
	$(CXX) $(CXXSTATICFLAGS) $(CXXFLAGS) -o oi oi.cxx $(INC) $(LIBS)
	mv oi ${ARTHOME}/bin/$(SYSTEM)

atra: atra.cxx
	$(CXX) $(CXXSTATICFLAGS) $(CXXFLAGS) -o atra atra.cxx $(INC) $(LIBS)
	mv atra ${ARTHOME}/bin/$(SYSTEM)

art2fsl: art2fsl.cxx
	$(CXX) $(CXXSTATICFLAGS) $(CXXFLAGS) -o art2fsl art2fsl.cxx $(INC) $(LIBS)
	mv art2fsl ${ARTHOME}/bin/$(SYSTEM)

fsl2art: fsl2art.cxx
	$(CXX) $(CXXSTATICFLAGS) $(CXXFLAGS) -o fsl2art fsl2art.cxx $(INC) $(LIBS)
	mv fsl2art ${ARTHOME}/bin/$(SYSTEM)

yuki: yuki.cxx
	$(MPICC) $(CXXFLAGS) -o yuki yuki.cxx $(INC) $(LIBS)
	mv yuki ${ARTHOME}/bin/$(SYSTEM)

analyze2nii: analyze2nii.cxx
	$(CXX) $(CXXSTATICFLAGS) $(CXXFLAGS) -o analyze2nii analyze2nii.cxx $(INC) $(LIBS)
	mv analyze2nii ${ARTHOME}/bin/$(SYSTEM)

multiplyTransformation: multiplyTransformation.cxx
	$(CXX) $(CXXSTATICFLAGS) $(CXXFLAGS) -o multiplyTransformation multiplyTransformation.cxx $(INC) $(LESSLIBS)
	mv multiplyTransformation ${ARTHOME}/bin/$(SYSTEM)

niftiinfo: niftiinfo.cxx
	$(CXX) $(CXXSTATICFLAGS) $(CXXFLAGS) -o niftiinfo niftiinfo.cxx $(INC) $(LESSLIBS)
	mv niftiinfo ${ARTHOME}/bin/$(SYSTEM)

break4D: break4D.cxx
	$(CXX) $(CXXSTATICFLAGS) $(CXXFLAGS) -o break4D break4D.cxx $(INC) $(LIBS)
	mv break4D ${ARTHOME}/bin/$(SYSTEM)

vancova: vancova.cxx
	$(CXX) $(CXXSTATICFLAGS) $(CXXFLAGS) -o vancova vancova.cxx $(INC) $(LIBS) -ldcdf  -lf2c
	mv vancova ${ARTHOME}/bin/$(SYSTEM)

gambling_vancova_setup: gambling_vancova_setup.cxx
	$(CXX) $(CXXSTATICFLAGS) $(CXXFLAGS) -o gambling_vancova_setup gambling_vancova_setup.cxx $(INC) $(LIBS)
	mv gambling_vancova_setup ${ARTHOME}/bin/$(SYSTEM)

medianImage: medianImage.cxx
	$(CXX) $(CXXSTATICFLAGS) $(CXXFLAGS) -o medianImage medianImage.cxx $(INC) $(LIBS)
	mv medianImage ${ARTHOME}/bin/$(SYSTEM)

threshold: threshold.cxx
	$(CXX) $(CXXSTATICFLAGS) $(CXXFLAGS) -o threshold threshold.cxx $(INC) $(LIBS)
	mv threshold ${ARTHOME}/bin/$(SYSTEM)

ccstats: ccstats.cxx
	$(CXX) $(CXXSTATICFLAGS) $(CXXFLAGS) -o ccstats ccstats.cxx $(INC) $(LIBS)
	mv ccstats ${ARTHOME}/bin/$(SYSTEM)

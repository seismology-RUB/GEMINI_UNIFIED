#-------------------------------------------------------------
#   Makefile for geminiUnified
#----------------------------------------------------------------------------
#   Copyright 2016 Wolfgang Friederich (Ruhr-Universitaet Bochum, Germany)
#
#   This file is part of GEMINI_UNIFIED version 1.0.
#
#   GEMINI_UNIFIED version 1.0 is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 2 of the License, or
#   (at your option) any later version.
#
#   GEMINI_UNIFIED version 1.0 is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with GEMINI_UNIFIED version 1.0.  If not, see <http://www.gnu.org/licenses/>.
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#  This is the generic part of the Makefile. Use this
#  template within each project.
#-----------------------------------------------------------------------
#  set the SHELL
#
SHELL = /usr/bin/bash
#
#  General definitions
#
hostname := $(shell hostname)
bindir = ./bin
obsdir = ./objects
moduledir = ./mod
fcversion = 10.2
OS := $(shell uname)
#

CC = gcc
F95 = gfortran

CFLAGS = -O3
FFLAGS = -O3 -J$(moduledir) -I/usr/include -Wunused-variable -Wuninitialized -fcheck=bounds
#-------------------------------------------------------
#  Directory search
#
vpath %.o $(obsdir)
vpath %.f90 f90 episode plotlib stuff fson miniseed specfemCoupling
vpath %.c  stuff miniseed
vpath %.f stuff f77
vpath %.h miniseed
#-------------------------------------------------------
#  Implicit rule to compile .o files from .f files.
#  Output of f2c stays in directory where .f files are.
#  Substitute .f by .c in dependency and compile target.
#  Then remove .c files.
#  Because of vpath, targets and dependencies need not be
#  in the current directory.
#
%.o: %.c
	$(CC) -c $(CFLAGS) $< -o $(obsdir)/$@
%.o: %.f90 
	$(F95) -c -I$(hdfinc) $(FFLAGS) $< -o $(obsdir)/$@
%.o: %.f 
	$(F95) -c $(FFLAGS) -fimplicit-none -ffixed-line-length-132 $< -o $(obsdir)/$@
#--------------------------------------------------------------
#  Object string for linking:
#  Adds object dir as prefix and removes directory part
#  of $^ (all dependencies)
#
obstring = $(addprefix $(obsdir)/,$(notdir $^))
#
#   End of generic part of Makefile
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#  Library paths
#
ifeq ($(hostname),yoga)
  x11 = -L/usr/lib/x86_64-linux-gnu -lX11
  pgplot = -L/usr/lib -lpgplot $(x11)
  lapack = /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3
  hdflib = -L/usr/local/hdf5/lib -lhdf5_fortran -lhdf5
  hdfinc = /usr/local/hdf5/include
else ifeq ($(hostname),espresso)
  x11 = /usr/lib/x86_64-linux-gnu/libX11.so.6
  pgplot = -L/usr/lib -lpgplot $(x11)
  lapack = /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3
  hdflib = -L/usr/lib/x86_64-linux-gnu/hdf5/openmpi -lhdf5_fortran -lhdf5
  hdfinc = /usr/lib/x86_64-linux-gnu/hdf5/openmpi/include
else ifeq ($(OS), Darwin)
	x11 = -L/usr/X11/lib/ -lX11
	lapack = -llapack
  	hdflib = -L/usr/local/opt/hdf5-parallel/lib -lhdf5_fortran -lhdf5
  	hdfinc = /usr/local/opt/hdf5-parallel/include/shared
  	mpilib = -L/usr/local/Cellar/open-mpi/4.0.5/lib/ -lmpi
  	pgplot = -L/usr/local/Cellar/pgplot/5.2.2_9/lib/ -lpgplot $(x11)
else
  x11 = -L/usr/lib/x86_64-linux-gnu -lX11
  pgplot = -L/usr/lib -lpgplot $(x11)
  lapack = -L/usr/lib/lapack -llapack
  hdfinc = /usr/lib/x86_64-linux-gnu/hdf5/openmpi/include 
  hdflib = -L/usr/lib/x86_64-linux-gnu/hdf5/openmpi/lib -lhdf5_fortran -lhdf5
  hdflink = -Wl,-rpath,/usr/lib/x86_64-linux-gnu/hdf5/openmpi/lib
  #hdfinc = /rscratch/minos01/wolle/for-gfortran-$(fcversion)/hdf5-1.10.5/hdf5/include
  #hdflib = -L/rscratch/minos01/wolle/for-gfortran-$(fcversion)/hdf5-1.10.5/hdf5/lib -lhdf5_fortran -lhdf5
  #hdflink = -Wl,-rpath,/rscratch/minos01/wolle/for-gfortran-$(fcversion)/hdf5-1.10.5/hdf5/lib
  mpilib = -L/usr/lib/x86_64-linux-gnu/openmpi/lib -lmpi
endif
#--------------------------------------------------------------
SA = fileStreamAccess.o datasetStreamAccess.o groupStreamAccess.o \
        flexibleType.o simpleString.o kindDefinitions.o primitiveTypeEncoding.o
STUFF = stuff.o gse20.o sunfortran.o sffFree.o sffSource.o sffInfo.o sffHeader.o \
        sffDatablock.o sffDataSection.o gseDatablock.o
FSON = fson.o fson_path_m.o fson_value_m.o fson_string_m.o
LIBMSEED = fileutils.o genutils.o gswap.o lmplatform.o lookup.o msrutils.o \
           pack.o packdata.o traceutils.o unpack.o unpackdata.o logging.o

#-------------------------------------------------------------------------
#  target for file .initdirs
#  to create module, objects and bin directories if file is not present
#
folders:
	if (! -e $(moduledir)) mkdir -p $(moduledir)
	if (! -e $(obsdir)) mkdir -p $(obsdir)
	if (! -e $(bindir)) mkdir -p $(bindir)
#------------------------------------------------------------------------
dependencies:
	./python/makeDepFromUseInclude.py f90 episode plotlib stuff fson specfemCoupling > $@	
#---------------------------------------------------------------
clean:
	-rm -f $(obsdir)/*.o
	-rm -f $(moduledir)/*.mod
	-rm -f dependencies
#----------------------------------------------------------------
#  include dependency file if available
#
-include dependencies
#----------------------------------------------------------------

gemini: folders computeGreenFKSpectraForSynthetics computeGreenFKSpectraForASKI computeSyntheticSeismograms computeTravelTimeTable \
	computeSyntheticsForSpecfemCoupling computeSyntheticsForASKI computeSeismogramsFromSpectraForSpecfemCoupling \
	computeSourceWavelet \
	prepareSpecfemCoupling plotGreenFKSpectra computeStationFilter plotHdfSeismogramGather plotAsciiSeismogramGather \

all: computeGreenFKSpectraForSynthetics computeGreenFKSpectraForASKI computeSyntheticSeismograms \
	computeEventFilter computeStationFilter buildSeismogramGather buildSffGather \
	plotGreenFKSpectra plotSeismogramGather testCompareSeismograms testReadNodeEarthmodel \
	computeKernelWavefield computeKernelGreenTensor
#----------------------------------------------------------------
mpiSupport.o: mpiSupport.f90
	mpif90 -c -J$(moduledir) $< -o $(obsdir)/$@
#-------------------------------------------------------------------------
computeGreenFKSpectraForSynthetics.o: computeGreenFKSpectraForSynthetics.f90
	mpif90 -c -I$(hdfinc) $(FFLAGS) $< -o $(obsdir)/$@
#-------------------------------------------------------------------------
computeGreenFKSpectraForASKI.o: computeGreenFKSpectraForASKI.f90
	mpif90 -c -I$(hdfinc) $(FFLAGS) $< -o $(obsdir)/$@
#-------------------------------------------------------------------------
computeSyntheticsForSpecfemCoupling.o: computeSyntheticsForSpecfemCoupling.f90
	mpif90 -c -I$(hdfinc) $(FFLAGS) $< -o $(obsdir)/$@
#-------------------------------------------------------------------------
hdfWrapper.o: hdfWrapper.f90
	mpif90 -c -I$(hdfinc) $(FFLAGS) $< -o $(obsdir)/$@
#----------------------------------------------------------------
computeGreenFKSpectraForSynthetics: %: %.o hdfWrapper.o reciprocity.o toroidalMinors.o externalRadialNodes.o \
	inputParameter.o unitJumpGreenFunctions.o mpiSupport.o nodeEarthmodel.o locatePoint.o cubicSpline.o mathConstants.o \
	argumentParser.o spheroidalMinors.o sourceTerms.o splineEarthmodelCoefficients.o errorMessage.o \
	geminiIntegrationEnvironment.o propagateOde.o initialValues.o toroidalOdeSystem.o radialOdeSystem.o \
	complexBulirschStep.o spheroidalAdjoint.o string.o constants.o complexElasticConstants.o fson.o realloc.o \
	fluidOdeSystem.o minorsOdeSystem.o baseIntegrationEnvironment.o baseStepEngine.o baseOdeSystem.o \
	integrationStep.o parametersBulirschStep.o solidtOdeSystem.o anyRankArray.o anyRankRealArray.o anyRankIntegerArray.o \
	$(FSON)
	mpif90 -o $(bindir)/$@ $(obstring) $(lapack) $(hdflib) $(hdflink)
#--------------------------------------------------------------------------------------------------
computeGreenFKSpectraForASKI: computeGreenFKSpectraForASKI.o unitJumpGreenFunctions.o solidtOdeSystem.o \
	toroidalOdeSystem.o reciprocity.o nodeEarthmodel.o errorMessage.o argumentParser.o solidOdeSystem.o \
	sourceTerms.o mpiSupport.o geminiIntegrationEnvironment.o splineEarthmodelCoefficients.o inputParameter.o \
	externalRadialNodes.o spheroidalMinors.o toroidalMinors.o spheroidalAdjoint.o radialOdeSystem.o \
	baseOdeSystem.o string.o constants.o realloc.o complexElasticConstants.o locatePoint.o mathConstants.o cubicSpline.o \
	baseIntegrationEnvironment.o complexBulirschStep.o minorsOdeSystem.o fluidOdeSystem.o \
	propagateOde.o initialValues.o parametersBulirschStep.o hdfWrapper.o baseStepEngine.o \
	integrationStep.o dsvDerivatives.o anyRankArray.o anyRankRealArray.o anyRankIntegerArray.o $(FSON)
	mpif90 -o $(bindir)/$@ $(obstring) $(lapack) $(hdflib) $(hdflink)
#------------------------------------------------------------------------------------
computeSyntheticSeismograms: computeSyntheticSeismograms.o singleFrequencyDisplacement.o seismicEvent.o \
	seismicStation.o geo2epi.o fourierTransform.o eventStationGeometry.o harmonicDegreeSums.o \
	pythag.o dateTime.o timeUtils.o locatePoint.o inputParameter.o asciiSynseisIO.o frequencyTime.o \
	errorMessage.o string.o constants.o argumentParser.o greenFKSpectra.o wavenumberIntegrals.o axesRotation.o \
	realloc.o readEventStationFile.o seismicEventList.o seismicNetwork.o asciiDataIO.o eventFilter.o \
	besselFunctions.o legendreFunctions.o hdfWrapper.o synseisHDF.o anyRankArray.o anyRankRealArray.o \
	anyRankIntegerArray.o externalRadialNodes.o butterworthFilter.o displacementSpectra.o travelTimes.o \
	nodeEarthmodel.o complexElasticConstants.o splineEarthmodelCoefficients.o cubicSpline.o $(FSON)
	$(F95) -o $(bindir)/$@ $(obstring) $(hdflib) $(hdflink) $(lapack)
#------------------------------------------------------------------------------------
computeSyntheticsForSpecfemCoupling: computeSyntheticsForSpecfemCoupling.o singleFrequencyDisplacement.o seismicEvent.o \
	fourierTransform.o harmonicDegreeSums.o axesRotation.o seismicStation.o complexElasticConstants.o asciiSynseisIO.o \
	pythag.o dateTime.o timeUtils.o locatePoint.o inputParameter.o frequencyTime.o seismicNetwork.o \
	errorMessage.o string.o constants.o argumentParser.o greenFKSpectra.o externalRadialNodes.o cubicSpline.o \
	realloc.o readEventStationFile.o seismicEventList.o asciiDataIO.o nodeEarthmodel.o heapSort.o \
	legendreFunctions.o hdfWrapper.o anyRankArray.o anyRankRealArray.o splineEarthmodelCoefficients.o \
	anyRankIntegerArray.o synseisHDF.o boundaryPoints.o gll_library.o eventFilter.o mpiSupport.o \
	displacementSpectra.o butterworthFilter.o wavenumberIntegrals.o besselFunctions.o travelTimes.o $(FSON) 
	mpif90 -o $(bindir)/$@ $(obstring) $(hdflib) $(hdflink) $(lapack)
#-------------------------------------------------------------------------------------------------------
computeSeismogramsFromSpectraForSpecfemCoupling: computeSeismogramsFromSpectraForSpecfemCoupling.o seismicEvent.o \
	fourierTransform.o dateTime.o timeUtils.o locatePoint.o inputParameter.o frequencyTime.o \
	errorMessage.o string.o constants.o argumentParser.o realloc.o hdfWrapper.o anyRankArray.o anyRankRealArray.o \
	anyRankIntegerArray.o synseisHDF.o eventFilter.o mpiSupport.o seismicStation.o \
	butterworthFilter.o asciiSynseisIO.o $(FSON) 
	mpif90 -o $(bindir)/$@ $(obstring) $(hdflib) $(hdflink) $(lapack)
#-------------------------------------------------------------------------------------------------------
computeSyntheticsForASKI: computeSyntheticsForASKI.o singleFrequencyDisplacement.o seismicEvent.o \
	fourierTransform.o harmonicDegreeSums.o axesRotation.o seismicStation.o complexElasticConstants.o \
	pythag.o dateTime.o timeUtils.o locatePoint.o inputParameter.o frequencyTime.o seismicNetwork.o \
	errorMessage.o string.o constants.o argumentParser.o greenFKSpectra.o externalRadialNodes.o cubicSpline.o \
	realloc.o readEventStationFile.o seismicEventList.o nodeEarthmodel.o geo2epi.o \
	legendreFunctions.o hdfWrapper.o anyRankArray.o anyRankRealArray.o splineEarthmodelCoefficients.o \
	anyRankIntegerArray.o synseisHDF.o eventFilter.o eventStationGeometry.o  asciiSynseisIO.o\
	displacementSpectra.o butterworthFilter.o wavenumberIntegrals.o besselFunctions.o \
	travelTimes.o $(FSON) 
	$(F95) -o $(bindir)/$@ $(obstring) $(hdflib) $(hdflink) $(lapack)
#-------------------------------------------------------------------------------------------------------
computeSourceWavelet: computeSourceWavelet.o dateTime.o miniSeed.o errorMessage.o string.o constants.o argumentParser.o \
	mathConstants.o constants.o inputParameter.o realloc.o timeSeries.o  timeUtils.o readMSeed.o fBindcUtils.o cBindcUtils.o \
	sourceTimeFunction.o pgPlotWindow.o pgPlotTimeSeries.o hdfWrapper.o anyRankArray.o \
	anyRankRealArray.o anyRankIntegerArray.o synseisHDF.o seismicStation.o seismicEvent.o fourierTransform.o \
	filterCoefficientsDecimate.o recursiveFilterCoefficients.o pythag.o nonNegativeLeastSquares.o \
	readEventStationFile.o seismicEventList.o seismicNetwork.o basePlotGather.o fourierSpectrum.o \
	interactiveBindingsPlotGather.o changeAxesLimitsPgplot.o asciiSynseisIO.o asciiDataIO.o $(LIBMSEED)
	$(F95) -o $(bindir)/$@ $(obstring) $(hdflib) $(hdflink) $(pgplot)
#---------------------------------------------------------------------------------------------------------
computeDyrtSourceWavelet: %: %.o anyRankArray.o anyRankIntegerArray.o anyRankRealArray.o argumentParser.o \
	asciiDataIO.o asciiSynseisIO.o constants.o dateTime.o errorMessage.o filterCoefficientsDecimate.o \
	fourierSpectrum.o fourierTransform.o frequencyTime.o hdfWrapper.o inputParameter.o mathConstants.o \
	miniSeed.o pythag.o readEventStationFile.o realloc.o recursiveFilterCoefficients.o seismicEvent.o \
	seismicEventList.o seismicNetwork.o seismicStation.o string.o timeSeries.o timeUtils.o \
        readMSeed.o fBindcUtils.o cBindcUtils.o $(LIBMSEED)
	$(F95) -o $(bindir)/$@ $(obstring) $(hdflib) $(hdflink)
#-------------------------------------------------------------------------------------------------------
prepareSpecfemCoupling: prepareSpecfemCoupling.o axesRotation.o asciiDataIO.o errorMessage.o realloc.o string.o \
	hdfWrapper.o gll_library.o boundaryPoints.o locatePoint.o anyRankArray.o anyRankRealArray.o anyRankIntegerArray.o \
	argumentParser.o inputParameter.o heapSort.o constants.o
	$(F95) -o $(bindir)/$@ $(obstring) $(hdflib) $(hdflink)
#----------------------------------------------------------------------------
plotGreenFKSpectra: plotGreenFKSpectra.o greenFKSpectra.o hdfWrapper.o externalRadialNodes.o \
	inputParameter.o nodeEarthmodel.o splineEarthmodelCoefficients.o complexElasticConstants.o \
	argumentParser.o string.o constants.o errorMessage.o pgPlotImage.o realloc.o pgColorScale.o pgPlotWindow.o \
	plotUtilities.o zoomBindingsPgPlotSelect.o locatePoint.o anyRankArray.o anyRankRealArray.o \
	anyRankIntegerArray.o cubicSpline.o $(FSON)
	$(F95) -o $(bindir)/$@ $(obstring) $(pgplot) $(hdflib) $(hdflink) $(lapack)
#------------------------------------------------------------------------------
extractSynseisHdfGatherFromSpecfemCoupling: extractSynseisHdfGatherFromSpecfemCoupling.o mathConstants.o \
	synseisHDF.o string.o constants.o argumentParser.o errorMessage.o anyRankArray.o anyRankRealArray.o \
	anyRankIntegerArray.o hdfWrapper.o seismicEvent.o seismicStation.o realloc.o dateTime.o \
	timeUtils.o inputParameter.o heapSort.o
	$(F95) -o $(bindir)/$@ $(obstring) $(hdflib) $(hdflink)
#-----------------------------------------------------------------------------------------------------
computeStationFilter: computeStationFilter.o argumentParser.o inputParameter.o errorMessage.o \
	greenFKSpectra.o string.o constants.o realloc.o readEventStationFile.o seismicEventList.o \
	seismicEvent.o locatePoint.o seismicNetwork.o seismicStation.o dateTime.o \
	timeUtils.o hdfWrapper.o anyRankArray.o anyRankRealArray.o anyRankIntegerArray.o \
	externalRadialNodes.o nodeEarthmodel.o complexElasticConstants.o splineEarthmodelCoefficients.o \
	cubicSpline.o $(FSON)
	$(F95) -o $(bindir)/$@ $(obstring) $(hdflib) $(hdflink) $(lapack)
#------------------------------------------------------------------------------------------------------
plotHdfSeismogramGather: plotHdfSeismogramGather.o hdfSynseisPlotGather.o errorMessage.o pgPlotWindow.o \
	argumentParser.o string.o constants.o basePlotGather.o synseisHDF.o realloc.o fourierSpectrum.o hdfWrapper.o \
	timeSeries.o fourierTransform.o recursiveFilterCoefficients.o interactiveBindingsPlotGather.o \
	filterCoefficientsDecimate.o dateTime.o mathConstants.o timeUtils.o pythag.o changeAxesLimitsPgplot.o \
	seismicStation.o seismicEvent.o anyRankArray.o anyRankRealArray.o anyRankIntegerArray.o asciiDataIO.o
	$(F95) -o $(bindir)/$@ $(obstring) $(pgplot) $(hdflib) $(hdflink)
#--------------------------------------------------------------------------------------------------------
computeTravelTimeTable: computeTravelTimeTable.o argumentParser.o nodeEarthmodel.o externalRadialNodes.o \
	rayIntegrationEnvironment.o splineEarthmodelCoefficients.o inputParameter.o errorMessage.o \
	mpiSupport.o hdfWrapper.o anyRankArray.o anyRankRealArray.o anyRankIntegerArray.o cubicSpline.o \
	hdfWrapper.o string.o constants.o realloc.o complexElasticConstants.o locatePoint.o rayOdeSystem.o \
	realBulirschStep.o propagateOde.o baseIntegrationEnvironment.o parametersBulirschStep.o \
	baseOdeSystem.o baseStepEngine.o integrationStep.o $(FSON)
	mpif90 -o $(bindir)/$@ $(obstring) $(pgplot) $(hdflib) $(hdflink) $(lapack)
#------------------------------------------------------------------------------------
computeKernelWavefield: computeKernelWavefield.o singleFrequencyDisplacement.o seismicEvent.o \
	locatePoint.o inputParameter.o geminiEarthModel.o chunkCubedSphere.o mpiSupport.o \
	fileUnitHandler.o errorMessage.o string.o constants.o argumentParser.o greenFKSpectra.o wavenumberIntegrals.o \
	realloc.o readEventStationFile.o seismicEventList.o coordinateTransform.o asciiDataIO.o \
	geo2epi.o dateTime.o timeUtils.o seismicStation.o seismicNetwork.o \
	besselFunctions.o $(SA)
	mpif90 -o $(bindir)/$@ $(obstring)
#---------------------------------------------------------------------------------------
computeKernelGreenTensor: computeKernelGreenTensor.o singleFrequencyDisplacement.o seismicEvent.o \
	locatePoint.o inputParameter.o geminiEarthModel.o chunkCubedSphere.o mpiSupport.o \
	fileUnitHandler.o errorMessage.o string.o constants.o argumentParser.o greenFKSpectra.o wavenumberIntegrals.o \
	realloc.o readEventStationFile.o seismicEventList.o coordinateTransform.o asciiDataIO.o \
	geo2epi.o dateTime.o timeUtils.o seismicStation.o seismicNetwork.o \
	besselFunctions.o $(SA)
	mpif90 -o $(bindir)/$@ $(obstring)
#--------------------------------------------------------------------------------------------------
computeEventFilter: computeEventFilter.o butterworthFilter.o argumentParser.o inputParameter.o errorMessage.o \
	greenFKSpectra.o string.o constants.o realloc.o readEventStationFile.o seismicEventList.o \
	seismicEvent.o locatePoint.o seismicNetwork.o seismicStation.o dateTime.o \
	timeUtils.o hdfWrapper.o anyRankArray.o anyRankRealArray.o anyRankIntegerArray.o
	$(F95) -o $(bindir)/$@ $(obstring) $(hdflib) $(hdflink)
#------------------------------------------------------------------------------------------------------
buildSeismogramGather: buildSeismogramGather.o argumentParser.o readEventStationFile.o asciiSynseisIO.o \
	asciiDataIO.o errorMessage.o inputParameter.o string.o constants.o realloc.o readEventStationFile.o seismicEventList.o \
	seismicNetwork.o seismicStation.o seismicEvent.o flexibleType.o dateTime.o timeUtils.o \
	simpleString.o
	$(F95) -o $(bindir)/$@ $(obstring)
#------------------------------------------------------------------------------------------------------
buildSffGather: buildSffGather.o argumentParser.o readEventStationFile.o \
	errorMessage.o inputParameter.o string.o constants.o realloc.o readEventStationFile.o seismicEventList.o \
	seismicNetwork.o seismicStation.o realloc.o seismicEvent.o flexibleType.o dateTime.o timeUtils.o \
	simpleString.o asciiSynseisIO.o asciiDataIO.o $(STUFF)
	$(F95) -o $(bindir)/$@ $(obstring)
#------------------------------------------------------------------------------------------------------
plotAsciiSeismogramGather: plotAsciiSeismogramGather.o asciiSynseisPlotGather.o errorMessage.o pgPlotWindow.o \
	argumentParser.o string.o constants.o basePlotGather.o asciiSynseisIO.o realloc.o fourierSpectrum.o asciiDataIO.o \
	timeSeries.o fourierTransform.o recursiveFilterCoefficients.o interactiveBindingsPlotGather.o \
	filterCoefficientsDecimate.o dateTime.o mathConstants.o timeUtils.o pythag.o changeAxesLimitsPgplot.o
	$(F95) -o $(bindir)/$@ $(obstring) $(pgplot)
#------------------------------------------------------------------------------------------------------
plotSpecfemSeismograms: plotSpecfemSeismograms.o specfemSynseisPlotGather.o errorMessage.o pgPlotWindow.o \
	argumentParser.o string.o constants.o basePlotGather.o realloc.o fourierSpectrum.o asciiDataIO.o \
	timeSeries.o fourierTransform.o recursiveFilterCoefficients.o interactiveBindingsPlotGather.o \
	filterCoefficientsDecimate.o dateTime.o mathConstants.o timeUtils.o pythag.o changeAxesLimitsPgplot.o \
	seismicStation.o seismicNetwork.o readEventStationFile.o seismicEvent.o seismicEventList.o
	$(F95) -o $(bindir)/$@ $(obstring) $(pgplot)
#------------------------------------------------------------------------------------------------------
plotNodeEarthmodel: plotNodeEarthmodel.o errorMessage.o pgPlotWindow.o locatePoint.o realloc.o \
	argumentParser.o string.o constants.o realloc.o pgPlotXY.o nodeEarthmodel.o splineEarthmodelCoefficients.o \
	complexElasticConstants.o cubicSpline.o changeAxesLimitsPgplot.o $(FSON)
	$(F95) -o $(bindir)/$@ $(obstring) $(pgplot) $(lapack)
#------------------------------------------------------------------------------------------------------
testCompareSeismograms: testCompareSeismograms.o argumentParser.o readEventStationFile.o asciiSynseisIO.o \
	errorMessage.o inputParameter.o string.o constants.o realloc.o readEventStationFile.o seismicEventList.o \
	seismicNetwork.o seismicStation.o seismicEvent.o flexibleType.o dateTime.o timeUtils.o \
	simpleString.o
	$(F95) -o $(bindir)/$@ $(obstring)
#-------------------------------------------------------------------------------------------------------
testReadNodeEarthmodel: testReadNodeEarthmodel.o nodeEarthmodel.o errorMessage.o argumentParser.o \
	string.o constants.o realloc.o locatePoint.o complexElasticConstants.o checkAssertion.o $(FSON)
	$(F95) -o $(bindir)/$@ $(obstring) $(pgplot)
#-------------------------------------------------------------------------------------------------------
testExternalRadialNodes: testExternalRadialNodes.o externalRadialNodes.o nodeEarthmodel.o errorMessage.o \
	argumentParser.o string.o constants.o realloc.o locatePoint.o complexElasticConstants.o inputParameter.o $(FSON)
	$(F95) -o $(bindir)/$@ $(obstring)
#-------------------------------------------------------------------------------------------------------
testBpIterator: testBpIterator.o axesRotation.o errorMessage.o realloc.o hdfWrapper.o string.o constants.o asciiDataIO.o \
	gll_library.o boundaryPoints.o locatePoint.o anyRankArray.o anyRankRealArray.o anyRankIntegerArray.o \
	heapSort.o
	$(F95) -o $(bindir)/$@ $(obstring) $(hdflib) $(hdflink)
#-------------------------------------------------------------------------------------------------------
testBpUniqueDistances: testBpUniqueDistances.o axesRotation.o errorMessage.o realloc.o hdfWrapper.o string.o asciiDataIO.o \
	gll_library.o boundaryPoints.o locatePoint.o anyRankArray.o anyRankRealArray.o anyRankIntegerArray.o \
	heapSort.o argumentParser.o inputParameter.o constants.o
	$(F95) -o $(bindir)/$@ $(obstring) $(hdflib) $(hdflink)
#-------------------------------------------------------------------------------------------------------
testHeapSort: testHeapSort.o heapSort.o
	$(F95) -o $(bindir)/$@ $(obstring)
#-------------------------------------------------------------------------------------------------------
testArrayPassing: testArrayPassing.o testArrayPassingSubs.o
	$(F95) -o $(bindir)/$@ $(obstring)
#-------------------------------------------------------------------------------------------
testButterworthImpulseResponse: testButterworthImpulseResponse.o errorMessage.o frequencyTime.o \
	readEventStationFile.o seismicEvent.o seismicStation.o realloc.o fourierTransform.o string.o constants.o \
	dateTime.o greenFKSpectra.o hdfWrapper.o anyRankArray.o anyRankRealArray.o anyRankIntegerArray.o \
	timeUtils.o seismicNetwork.o inputParameter.o argumentParser.o seismicEventList.o \
	externalRadialNodes.o synseisHDF.o complexElasticConstants.o splineEarthmodelCoefficients.o \
	locatePoint.o nodeEarthmodel.o cubicSpline.o eventFilter.o butterworthFilter.o $(FSON)
	$(F95) -o $(bindir)/$@ $(obstring) $(hdflib) $(hdflink) $(lapack)
#---------------------------------------------------------------------------------------------------------
testTravelTimes: testTravelTimes.o argumentParser.o errorMessage.o hdfWrapper.o realloc.o anyRankArray.o \
	anyRankRealArray.o anyRankIntegerArray.o locatePoint.o mathConstants.o string.o constants.o \
	travelTimes.o
	$(F95) -o $(bindir)/$@ $(obstring) $(hdflib) $(hdflink)

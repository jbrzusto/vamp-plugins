## NB: on linux, do a make clean whenever there are significant changes to the class structure
## of any plugin.
##
##  Skeleton Makefile for Vamp plugin builds using command-line tools.
##
##  Rename this to Makefile, and edit as appropriate.
##  This Makefile WILL NOT WORK until you have edited it as described
##  below -- the Makefile as supplied does nothing useful at all!
##
##  Various sets of options are provided, commented out -- just uncomment
##  (remove the '#' characters for) the set that most closely resembles
##  your own situation, and adjust to taste.  Then run "make".
##
##  (For Windows builds using MS Visual Studio, start instead with the
##  VampExamplePlugins project found in the build directory of the SDK.)


# Edit this to the base name of your plugin library
#
PLUGIN_LIBRARY_NAME = lotek-plugins

# Where the plugin library should be installed:

PLUGIN_INSTALL_DIR = ~/.vamp

# Edit this to list one .o file for each .cpp file in your plugin project
#
PLUGIN_CODE_OBJECTS = FindPulseFD.o FindPulseFDBatch.o FindPulseTD.o FindPulseTDBatch.o NewFindPulse.o

# Edit this to the location of the Vamp plugin SDK, relative to your
# project directory
#
VAMP_SDK_INCLUDE_DIR = /usr/include/
VAMP_SDK_LIB_DIR = /usr/lib

CXXFLAGS_PRODUCTION = -std=c++0x -I$(VAMP_SDK_INCLUDE_DIR) -Wall -fPIC -O3 -ffast-math  -ftree-vectorize
CXXFLAGS_DEBUG = -std=c++0x -I$(VAMP_SDK_INCLUDE_DIR) -Wall -fPIC -g3 

CXXFLAGS = $(CXXFLAGS_DEBUG)
##CXXFLAGS = $(CXXFLAGS_PRODUCTION)
PLUGIN_EXT = .so
PLUGIN = $(PLUGIN_LIBRARY_NAME)$(PLUGIN_EXT)
LDFLAGS = -shared -Wl,-soname=$(PLUGIN) $(VAMP_SDK_LIB_DIR)/libvamp-sdk.a -L/usr/local/lib -lfftw3f -Wl,--version-script=vamp-plugin.map


$(PLUGIN): $(PLUGIN_CODE_OBJECTS) lotek-plugins.o
	   $(CXX) -o $@ $^ $(LDFLAGS)

lotek-plugins.o: $(PLUGIN_CODE_OBJECTS) lotek-plugins.cpp

FindPulseFD.o: FindPulseFD.h PulseFinder.h

FindPulseFDBatch.o: FindPulseFDBatch.h PulseFinder.h

FindPulseTD.o: FindPulseFD.h ProbPulseFinder.h 

FindPulseTDBatch.o: FindPulseTDBatch.h PulseFinder.h

FindPulse.o: FindPulse.h MovingAverager.h

FindPulseBatch.o: FindPulseBatch.h MovingAverager.h

AveragePower.o: AveragePower.h MovingAverager.h

NewFindPulse.o: NewFindPulse.h SortedFifo.h DynamicQuantile.h PeakFinder.h MovingAverager.h FreqEstimator.h EdgeDetector.h PulseDetector.h KahanSum.h

FreqEstimator.o: FreqEstimator.h

testFreqEstimator.o: FreqEstimator.h testFreqEstimator.cpp

testFreqEstimator: testFreqEstimator.o FreqEstimator.o
	$(CXX) -o $@ $^ -L/usr/local/lib -lfftw3f

testKahanSum: KahanSum.h testKahanSum.cpp
	g++ -$(CXXFLAGS) -o testKahanSum testKahanSum.cpp

testPeakFinder: PeakFinder.h testPeakFinder.cpp SortedFifo.h
	g++ $(CXXFLAGS) -o testPeakFinder testPeakFinder.cpp

testDynamicQuantile: DynamicQuantile.h testDynamicQuantile.cpp
	g++ $(CXXFLAGS) -o testDynamicQuantile testDynamicQuantile.cpp

testEdgeDetector: EdgeDetector.h testEdgeDetector.cpp
	g++ $(CXXFLAGS) -o testEdgeDetector testEdgeDetector.cpp

testPulseDetector: PulseDetector.h testPulseDetector.cpp 
	g++ $(CXXFLAGS) -o testPulseDetector testPulseDetector.cpp

clean:
	rm -f *.o

install: $(PLUGIN)
	cp -f $(PLUGIN) $(PLUGIN_INSTALL_DIR)


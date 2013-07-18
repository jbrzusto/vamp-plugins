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
PLUGIN_CODE_OBJECTS = FindPulseFD.o FindPulseFDBatch.o FindPulseTD.o FindPulseTDBatch.o

# Edit this to the location of the Vamp plugin SDK, relative to your
# project directory
#
VAMP_SDK_INCLUDE_DIR = /usr/include/
VAMP_SDK_LIB_DIR = /usr/lib

CXXFLAGS_PRODUCTION = -I$(VAMP_SDK_INCLUDE_DIR) -Wall -fPIC -O3 -ffast-math -fsingle-precision-constant -ftree-vectorize
CXXFLAGS_DEBUG = -I$(VAMP_SDK_INCLUDE_DIR) -Wall -fPIC -g3 -ffast-math -fsingle-precision-constant -ftree-vectorize
## CXXFLAGS = $(CXXFLAGS_DEBUG)
CXXFLAGS = $(CXXFLAGS_PRODUCTION)
PLUGIN_EXT = .so
PLUGIN = $(PLUGIN_LIBRARY_NAME)$(PLUGIN_EXT)
LDFLAGS = -shared -Wl,-soname=$(PLUGIN) $(VAMP_SDK_LIB_DIR)/libvamp-sdk.a -lfftw3f -Wl,--version-script=vamp-plugin.map


$(PLUGIN): $(PLUGIN_CODE_OBJECTS) lotek-plugins.o
	   $(CXX) -o $@ $^ $(LDFLAGS)

lotek-plugins.o: $(PLUGIN_CODE_OBJECTS) lotek-plugins.cpp

FindPulseFD.o: FindPulseFD.h PulseFinder.h

FindPulseFDBatch.o: FindPulseFDBatch.h PulseFinder.h

FindPulseTD.o: FindPulseFD.h PulseFinder.h

FindPulseTDBatch.o: FindPulseTDBatch.h PulseFinder.h

FindPulse.o: FindPulse.h MovingAverager.h

FindPulseBatch.o: FindPulseBatch.h MovingAverager.h

AveragePower.o: AveragePower.h MovingAverager.h

clean:
	rm -f *.o

install: $(PLUGIN)
	cp -f $(PLUGIN) $(PLUGIN_INSTALL_DIR)


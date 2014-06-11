/* -*- c-basic-offset: 4 indent-tabs-mode: nil -*-  vi:set ts=8 sts=4 sw=4: */

/*
    Lotek-plugins

    VAMP plugins for detecting Lotek telemetry tags.
    Copyright 2012 John Brzustowski

    Based on the examples from the VAMP-SDK
       An API for audio analysis and feature extraction plugins.
       Centre for Digital Music, Queen Mary, University of London.
       Copyright 2006 Chris Cannam.
  
    Permission is hereby granted, free of charge, to any person
    obtaining a copy of this software and associated documentation
    files (the "Software"), to deal in the Software without
    restriction, including without limitation the rights to use, copy,
    modify, merge, publish, distribute, sublicense, and/or sell copies
    of the Software, and to permit persons to whom the Software is
    furnished to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be
    included in all copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
    EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
    MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
    NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS BE LIABLE FOR
    ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF
    CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
    WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

    Except as contained in this notice, the names of the Centre for
    Digital Music; Queen Mary, University of London; and Chris Cannam
    shall not be used in advertising or otherwise to promote the sale,
    use or other dealings in this Software without prior written
    authorization.
*/

#include "vamp/vamp.h"
#include "vamp-sdk/PluginAdapter.h"


#include "FindNarrowBandPulsePlugin.h"
#include "SNRZFindPulse.h"
#include "FindPulseFDBatch.h"
#include "FindPulseFD.h"
#include "FindPulseTDSNRZ.h"

static Vamp::PluginAdapter<FindNarrowBandPulsePlugin> FindNarrowBandPulsePluginAdapater;
static Vamp::PluginAdapter<SNRZFindPulse> SNRZFindPulseAdapater;
static Vamp::PluginAdapter<FindPulseFDBatch> findPulseFDBatchAdapter;
static Vamp::PluginAdapter<FindPulseFD> findPulseFDAdapter;
static Vamp::PluginAdapter<FindPulseTDSNRZ> findPulseTDSNRZAdapter;

const VampPluginDescriptor *vampGetPluginDescriptor(unsigned int version,
                                                    unsigned int index)
{
    if (version < 1) return 0;

    switch (index) {
    case  0: return FindNarrowBandPulsePluginAdapater.getDescriptor();
    case  1: return SNRZFindPulseAdapater.getDescriptor();
    case  2: return findPulseFDBatchAdapter.getDescriptor();
    case  3: return findPulseFDAdapter.getDescriptor();
    case  4: return findPulseTDSNRZAdapter.getDescriptor();
    default: return 0;
    }
}

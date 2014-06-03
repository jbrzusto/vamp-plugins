 /* -*- mode:c++; c-basic-offset: 4 indent-tabs-mode: nil -*-  vi:set ts=8 sts=4 sw=4: 

    vamp-plugins-common.h - common stuff for vamp-plugins 

    Copyright 2014 John Brzustowski

    License: GPL v 2.0 or later.  This is required in order to use fftw.

*/

#ifndef _VAMP_PLUGINS_COMMON_H_
#define _VAMP_PLUGINS_COMMON_H_

#ifdef MINGW
#define exp10(X) powf(10, X)
#define fftw_free(X) fftwf_free(X)
#endif

#define dB(X) (10.0 * log10(X))
#define undB(X) (exp10((X) / 10.0))

#endif /*  _VAMP_PLUGINS_COMMON_H_ */

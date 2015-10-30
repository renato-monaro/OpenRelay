# Copyright (C) 2008-2015 Renato Machado Monaro
# This file is part of OpenRelay.
#
# OpenRelay is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# OpenRelay is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
#ifndef TYPES_H
#define TYPES_H

#include <iostream>
//#include <boost/date_time/posix_time/posix_time_types.hpp>
#include <vector>
#include <complex>
#include <climits>
#include "parameters.h"
#ifdef RTAI
	#include<rtai_lxrt.h>
#endif
#ifndef RTAI
	#include<time.h>
#endif


using namespace std;


enum{ANALOG=0,DIGITAL,TIME,COMPLEX,STRING,NONE};
enum {OUTPUT=0,INPUT};
enum {FIS=0,FCL,FLL};
typedef  bool digital;
typedef double analogic;
//typedef std::complex<double> phasor;
//typedef boost::posix_time::ptime timer;
typedef long long timer;
typedef std::complex<analogic> Complex;

typedef vector<analogic> Vanalogic;

long long orelay_gettime(void);

template <class T> inline T magsq(complex<T> x)
{
  T y;
  y = (x.real()*x.real() + x.imag()*x.imag());
  return y;
};

inline complex<double> expj(double x)
{
  complex<double> y;
  y.real(cos(x));
  y.imag(sin(x));
  return y;
}


#endif /* !TYPES_H */

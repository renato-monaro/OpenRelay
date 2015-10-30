/*
# Copyright (C) 2010-2015 Rapahel Phillipe
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
*/
/*
 * WangSunFrequency.h
 *
 *  Created on: 29/07/2010
 *      Author: rapphil
 */

#ifndef WANGSUNFREQUENCY_H_
#define WANGSUNFREQUENCY_H_

#include <boost/circular_buffer.hpp>

#include "channel.h"
#include "measures.h"
#include "types.h"
namespace orelay{
class WangSunFrequency : public Measure{
public:
	double ERRO;
	WangSunFrequency(Channel<Complex> * ,Channel<analogic > *,Channel<digital> *TR,unsigned ,unsigned, float );
	void Run();
	virtual ~WangSunFrequency();
private:
	Channel<Complex> * CH_IN;
	Channel<analogic> * CH_OUT;
	Channel<digital> * TR;
	boost::circular_buffer<double> phases;
	unsigned M,N,sampleCount,jumps,countRuns;
	double K,K1,K2,K3,phaseAnt;
	float RelaySamplingRate;
	};
}
#endif /* WANGSUNFREQUENCY_H_ */

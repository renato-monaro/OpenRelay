/*
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
*/
#ifndef ENCODER_H
#define ENCODER_H

#include <string>
#include <fstream>
#include <pthread.h>
#include <vector>
#ifdef RTAI
	#ifdef COMEDI
		#include <rtai_comedi.h>
	#endif
#endif
#ifndef RTAI
	#ifdef COMEDI
		#include <comedilib.h>
	#endif
#endif
#include "channel.h"
#include "types.h"
#include "parameters.h"
#include "acquisition.h"

namespace orelay{
#ifdef COMEDI
class Encoder: public Acquisition{
public:
		Encoder(const char *device, Channel<analogic>  *out, unsigned holes, int a, int b, int z);
		Encoder(const char *device, Channel<analogic> *out, unsigned holes, int a, int b);
		~Encoder();
		bool Run(); 
		bool Prepare(float);
		void RefreshTime();
	protected:
		Channel<analogic> *Out;
		int A,B,Z;
		int subDevice;
		double Ratio;
	private:
		comedi_t *comediDev;
		lsampl_t data,data_ant;
};

#endif
}
#endif

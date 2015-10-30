/*
# Copyright (C) 2008-2015 Renato Machado Monaro, Rafael Marsolla
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
#ifndef thd_h
#define thd_h

#include<relay.h>
#include<string>
namespace orelay{
// Compute Power Factor (FP) with or without THD
class FP: public Measure{
	public: 
	FP(Channel<analogic> *CH_V1, Channel<analogic> *CH_I1, Channel<analogic> *CH_FP1, float Freq);
	FP(Channel<analogic> *CH_V1, Channel<analogic> *CH_I1, Channel<analogic> *CH_THD, Channel<analogic> *CH_FP1, float Freq);
	bool Prepare(float Samp);
	void Run();
		
	protected:
	Complex DFT(Channel<analogic> *CH_IN);
	Channel<analogic> *InputV; 
	Channel<analogic> *InputI; 
	Channel<analogic> *Thd; 
	Channel<analogic> *Output; 
	bool WithTHD;
	float SystemFrequency;
	unsigned N;
};



//Compute the Total Harmonic Distortion - THD
class THD: public Measure {
	public:
		THD(Channel<analogic>* CH_IN,Channel<analogic>* CH_OUT, float freq, unsigned Nharm); 
		bool Prepare(float Samp); 
		void Run();
	protected:
		Channel<analogic> *Input;
		Channel<analogic> *Output; 
		float SystemFrequency; 
		unsigned Harmonic; 
		unsigned N;
};
}
#endif

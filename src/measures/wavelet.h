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
#ifndef WAVELET_H
#define WAVELET_H

#include <vector>
#include "measures.h"
#include "channel.h"  
namespace orelay{
enum {haar=0,daub4,daub6,daub8,daub10,daub12,daub14,daub16,daub18,daub20,daub22,daub24,daub26,daub28,\
	daub30,daub32,daub34,daub36,daub38,daub40,daub42,daub44,daub46,daub48,\
	daub50,daub52,daub54,daub56,daub58,daub60,daub62,daub64,daub66,daub68,daub70,\
	daub72,daub74,daub76,sym8,sym16,coif6,coif12,coif18,coif24,coif30,beylkin18,vaidyanathan24\
	};

enum {WAVELET_RMS=0,WAVELET_ENERGY};//,WAVELET_ENTROPY,WAVELET_FILTER};

class Wavelet: public Measure{
	public:
		Wavelet(Channel<analogic> *CH_IN,Channel<analogic> *CH_OUT, unsigned N, unsigned Family,unsigned Met);
		void Leafs(unsigned char ad,unsigned Nad,unsigned Level);
		virtual void Run(); 
		bool Prepare(float);
	protected:
		bool getFilterBank();
		Channel <analogic> *Input;
		Channel <analogic> *Output;
		vector<double> FilterBank;
		unsigned NSamples;
		unsigned FamilyName;
		unsigned Method;
		vector<unsigned char> AD;
		vector<unsigned> Nivel,Pos;
		vector<unsigned> Ai,Di,Af,Df;
		vector<unsigned> PosIni,PosFim;
		vector<bool> Ordem;
		void Transformada_Wavelet(double* f, long n,bool ordem);
	};
}
#endif/* !WAVELET_H */

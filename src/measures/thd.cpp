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
#include "thd.h"
using namespace orelay;
// FP without THD constructor
FP::FP(Channel<analogic> *CH_V1, Channel<analogic> *CH_I1, Channel<analogic> *CH_FP1, float Freq){
	FunctionalDescription_set("Power Factor"); 
	InputV=CH_V1;
	InputI=CH_I1;
	Output=CH_FP1;
	SystemFrequency=Freq;
	WithTHD=false;
}

// FP with THD constructor
FP::FP(Channel<analogic> *CH_V1, Channel<analogic> *CH_I1, Channel<analogic> *CH_THD, Channel<analogic> *CH_FP1, float Freq){
	FunctionalDescription_set("Power Factor"); 
	InputV=CH_V1;
	InputI=CH_I1;
	Thd=CH_THD;
	Output=CH_FP1;
	SystemFrequency=Freq;
	WithTHD=true;
}

// Get SamplingRate
bool FP::Prepare(float Samp){
        RelaySamplingRate=Samp;
        N=(int)(RelaySamplingRate/SystemFrequency);
        return true;
}

//Compute DFT to FP
Complex FP::DFT(Channel<analogic> *CH_IN){
        Complex temp, i(0.0,1.0);
        temp = Complex(0.0,0.0);
        // fundamental - N1
        for(int window=0; window < N ; window++){
                temp += CH_IN->get_Value(window)*exp(((analogic)(2*M_PI/N)*i)*(analogic)(window*1));
        }
        temp/=(analogic)N;
        //temp*=sqrt(2);
	return  temp;
}

//Compute FP
void FP::Run(){

	// Pure Sinoidal 
	analogic FP = cos(arg(FP::DFT(InputV)) - arg(FP::DFT(InputI)));

	// With THD
	if(WithTHD){
		FP /=sqrt(1+ pow(Thd->get_Value(),2));
	}
	Output->insert_Value(FP);	
}

// THD Constructor
THD::THD(Channel<analogic> *CH_IN,Channel<analogic> *CH_OUT,float Freq, unsigned Nharm){	
	FunctionalDescription_set("Discret Fourier Transform"); 
	Input=CH_IN;
	Output=CH_OUT;
	Harmonic=Nharm;
	SystemFrequency=Freq;
	}

// Get SamplingRate
bool THD::Prepare(float Samp){
	RelaySamplingRate=Samp;
	N=(unsigned)(RelaySamplingRate/SystemFrequency);
	return true;
	}

// Compute THD
void THD::Run(){
	analogic tempout = 0, tempfund =0;
	Complex temp, i(0.0,1.0);

	 // Fundamental N1	
	temp = Complex(0.0,0.0);
	for(int window=0; window < N ; window++){
		temp += Input->get_Value(window)*exp(((analogic)(2*M_PI/N)*i)*(analogic)(window*1));
	}
	temp/=(analogic)N;
	//temp*=sqrt(2);
	tempfund = pow(abs(temp),2);
        
	// 2 to N harmonics
	for(int harm =2; harm < Harmonic; harm++){
		temp = Complex(0.0,0.0);
		for(int window=0; window < N ; window++){
			temp += Input->get_Value(window)*exp(((analogic)(2*M_PI/N)*i)*(analogic)(window*harm));
		}
		temp/=(analogic)N;
		//temp*=sqrt(2);

		tempout += pow(abs(temp),2);
	} 
	Output->insert_Value(sqrt(tempout/tempfund));
}

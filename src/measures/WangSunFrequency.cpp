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
/*
 * WangSunFrequency.cpp
 *
 *  Created on: 29/07/2010
 *      Author: rapphil
 */

#include "WangSunFrequency.h"
#include <math.h>
using namespace orelay;
WangSunFrequency::WangSunFrequency(Channel<Complex> * IN, Channel<analogic> * OUT, Channel<digital> * TR, unsigned N, unsigned M, float Samp) 
{
	this->CH_IN=IN;
	this->CH_OUT=OUT;
	this->TR=TR;
	this->M=M;
	this->N=N;
	this->sampleCount=0;
	this->jumps=0;
	this->phases = boost::circular_buffer<double>(M+1);
	this->countRuns=0;
	this->K=M_PI/(N*sin(2*M_PI/double(N)));
	this->K3=2*M_PI*M/(double)N;
	this->RelaySamplingRate=Samp;
	ERRO=1E-6;
}

void WangSunFrequency::Run()
{
	double phase=atan(CH_IN->get_Value().imag()/CH_IN->get_Value().real());
	
	if(fabs(phase-phaseAnt)>0.8*M_PI)
		++jumps;
		
	phaseAnt=phase;
	phases.push_front(phase+M_PI*jumps);
	
	if(this->M <countRuns )
	{
		double deltaEps, deltaEpsAnt=2,frequency,ang1,ang2;
		unsigned count=0;
		
		ang1 = phases.at(0);
		ang2 = phases.at(this->M);
		this->K1=sin(2*M_PI/(double)N -2*ang1);
		this->K2=sin(2*M_PI/(double)N -2*ang2);
		//this->K=M_PI/(N*sin(2*M_PI/double(N)));
		deltaEps=((ang1-ang2)-this->K3)/(this->K3+this->K*(this->K1-this->K2));
		
		while(fabs(deltaEps-deltaEpsAnt)>ERRO && count<1000)
		{
			deltaEpsAnt=deltaEps;
			this->K=M_PI/(this->N*sin((2+deltaEps)/double(this->N)*M_PI));
			deltaEps=((ang1-ang2)-this->K3)/(this->K3+this->K*(this->K1-this->K2));
			count++;
		}
		
		frequency=this->RelaySamplingRate*(1.0+deltaEps)/double(this->N);
		
		if(fabs(frequency)>100)
			this->CH_OUT->insert_Value(0);	
		else
			this->CH_OUT->insert_Value(fabs(frequency));
			
		this->TR->insert_Value(true);
	}
	else
	{
		CH_OUT->insert_Value(0.0);
		this->TR->insert_Value(false);
		countRuns++;
	}
//	}


}

WangSunFrequency::~WangSunFrequency() 
{
	// TODO Auto-generated destructor stub
}

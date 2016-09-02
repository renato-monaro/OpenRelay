/*
# Copyright (C) 2008-2015 Rodolfo Varraschim Rocha
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
#
# DDF is based on: D. Hou, A. Guzman, J. Roberts, Innovative solutions improve transmis-
#   sion line protection, in: 24th Annual Western Protective Relay Confer-
#   ence, Spokane, WA, Vol. 21, Citeseer, 1997.
#
*/

#include "ddf.h"

using namespace orelay;

DDF::DDF(Channel<analogic> *CH_IN,Channel<analogic> *CH_OUT,float Freq){
			FunctionalDescription_set("Double Differentiator Filter");
			Input=CH_IN;
			Output=CH_OUT;
			SystemFrequency=Freq;
        }
bool DDF::Prepare(float Samp){
            float R,I,K;

			RelaySamplingRate=Samp;
			N=(unsigned)(RelaySamplingRate/SystemFrequency);

			R=1-2*cos(2*M_PI/N)+cos(4*M_PI/N);
            I=2*sin(2*M_PI/N)-sin(4*M_PI/N);
            K=sqrt(1/(pow(R,2)+pow(I,2)));
            b[0]=K; b[1]=-2*K; b[2]=K;
            a[0]=1;
            z[0]=0; z[1]=0; z[2]=0;

            cout<<"DDF Compensation Enabled"<<endl;

			return true;
		}
void DDF::Run(){
            float Y;

            Y=b[0]*Input->get_Value()+z[0];
            for (unsigned j=1;j<3;j++){
                    z[j-1]=b[j]*Input->get_Value()+z[j];
            }
            Output->insert_Value(-Y);

            //cout<<Input->get_Value()<<" "<<Output->get_Value()<<endl;
		}


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
#include "switch.h"
using namespace orelay;
Switch::Switch(Channel<analogic> *CH_1,Channel<analogic> *CH_2,Channel<analogic> *CH_OUT, Channel<digital> *CH_S){
	FunctionalDescription_set("Switch"); 
	Input_1=CH_1; 
	Input_2=CH_2;
	Input_D=CH_S;
	Output=CH_OUT;
	withfix=false;
	}
Switch::Switch(double FIX,Channel<analogic> *CH_1,Channel<analogic> *CH_OUT, Channel<digital> *CH_S){
	FunctionalDescription_set("Switch"); 
	Input_1=CH_1; 
	Input_D=CH_S;
	Output=CH_OUT;
	FixValue=FIX;
	withfix=true;
	}
void Switch::Run(){
	if(!Input_D->get_Value()){
		if(withfix)
			Output->insert_Value(FixValue);
		else
			Output->insert_Value(Input_2->get_Value());
		}
	else{
		Output->insert_Value(Input_1->get_Value());
		}
	}


SwitchPos::SwitchPos(Channel<analogic> *CH_OUT, Channel<digital> *CH_S){
	Input_D=CH_S;
	Output=CH_OUT;
	pos=0;
	};

bool SwitchPos::Join(Channel<analogic> *CH_1){
	type.push_back(T_CHANNEL);
	vecChannel.push_back(CH_1);
	index.push_back(vecChannel.size()-1);
	return true;
	}

bool SwitchPos::Join(double fix){
	type.push_back(T_FIXVALUE);
	vecfixvalue.push_back(fix);
	index.push_back(vecfixvalue.size()-1);
	return true;
	}

void SwitchPos::Run(){
	if(Input_D->get_Value()&!Input_D->get_Value(1)){
		pos++;
		if(pos>=type.size())
			pos=0;
		}
	Output->insert_Value(getValue());
	}

double SwitchPos::getValue(){
	if(type[pos]==T_CHANNEL)
		return vecChannel[index[pos]]->get_Value();
	if(type[pos]==T_FIXVALUE)
		return vecfixvalue[index[pos]];
	}


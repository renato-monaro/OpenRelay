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
#include "logic.h"
using namespace orelay;
void Logic::SetClock(Channel<timer> *t){
	MasterClock=t;
	}

void Logic::FunctionalDescription_set(string Desc){
	FunctionalDescription=Desc;
	}

string Logic::FunctionalDescription_get(){
	return FunctionalDescription;
	}

bool Logic::Prepare(float Samp){
	RelaySamplingRate=Samp;
	return true;
	}

void Multiple::Join_Channel(Channel<digital> *In){
	Input.push_back(In);
	}

Not::Not(Channel<digital> *In,Channel<digital> *Out){
	FunctionalDescription_set("Not Logic");
	Input=In;
	Output=Out;
	}
void Not::Run(){	
	Output->insert_Value(!Input->get_Value());
	}

Or::Or(Channel<digital> *Out){
	FunctionalDescription_set("Or Logic");
	Output=Out;
	}

void Or::Run(){
	bool tmp=false;
	for(unsigned i=0;i<Input.size();i++)
		tmp|=Input[i]->get_Value();
	Output->insert_Value(tmp);
	}

And::And(Channel<digital> *Out){
	Output=Out;
	FunctionalDescription_set("And Logic");
	}

void And::Run(){
	bool tmp=true;
	for(unsigned i=0;i<Input.size();i++)
		tmp&=Input[i]->get_Value();
	Output->insert_Value(tmp);
	}

FlipFlop::FlipFlop(Channel<digital> *InR,Channel<digital> *InS,Channel<digital> *Out){
	FunctionalDescription_set("FlipFlop JK Logic");
	InputS=InS;
	InputR=InR;
	Output=Out;
	}
void FlipFlop::Run(){
	if((!InputS->get_Value())&&(!InputR->get_Value())) //0,0=Q
		Output->insert_Value(Output->get_Value());
	if((!InputS->get_Value())&&(InputR->get_Value())) //0,1=1
		Output->insert_Value(false);
	if((InputS->get_Value())&&(!InputR->get_Value())) //1,0=0
		Output->insert_Value(true);
	if((InputS->get_Value())&&(InputR->get_Value())) //1,1=Q
		Output->insert_Value(!Output->get_Value());
	}


bool Timer::Prepare(float Samp){
	RelaySamplingRate=Samp;
	TimeStep=1E9/Samp;
	cout<<"TimerSTep="<<TimeStep<<endl;
	return true;
	}

Timer::Timer(Channel<digital> *In,Channel<digital> *Out,timer Time,timer ResetTime){
	FunctionalDescription_set("Timer Logic");
	Dynamic=false;
	Trigger=true;
	Input=In;
	Output=Out;
	Delay=Time;
	Resetable=true;
	Reset=ResetTime;
	Accumulator=0;
//	TimeStep=1E9/In->get_Sampling();
}

Timer::Timer(Channel<digital> *In,Channel<digital> *Out,timer Time){
	FunctionalDescription_set("Timer Logic");
	Dynamic=false;
	Trigger=true;
	Input=In;
	Output=Out;
	Delay=Time;
	Resetable=false;
	Reset=0;
	Accumulator=0;
//	TimeStep=1E9/In->get_Sampling();
}

Timer::Timer(Channel<digital> *In,Channel<digital> *Out,Channel<timer> *Time,timer ResetTime){
	FunctionalDescription_set("Timer Logic");
	Dynamic=true;
	Trigger=true;
	Input=In;
	Output=Out;
	DynamicDelay=Time;
	Resetable=true;
	Accumulator=0;
	Reset=ResetTime;
//	TimeStep=1E9/In->get_Sampling();
}

Timer::Timer(Channel<digital> *Out,timer Time){
	FunctionalDescription_set("Timer Logic");
	Dynamic=false;
	Trigger=false;
	Output=Out;
	Delay=Time;
	Resetable=false;
	Accumulator=0;
//	TimeStep=1E9/In->get_Sampling();
}

void Timer::Run(){
	if(Trigger){
		if(Input->get_Value())//active timer
			Accumulator+=TimeStep;
		else
			Accumulator=0;
		}
	else{
		Accumulator+=TimeStep;
		}
	if(Dynamic)
		Delay=DynamicDelay->get_Value();
	if(Resetable){
		if((Accumulator>Delay)&&(Accumulator<(Reset+Delay)))
			Output->insert_Value(true);
		else
			Output->insert_Value(false);
		}
	else{
		if(Accumulator>=Delay){
			Output->insert_Value(true);
			}
		else
			Output->insert_Value(false);
		}
}


Selector::Selector(Channel<digital> *In,Channel<string> *Out, vector<string> V){
	Input=In;
	SOutput=Out;
	resetable=false;
	Type=STRING;
	for(unsigned i=0;i<V.size();i++)
		SValues.push_back(V[i]);
	ActualPos=0;
	}

Selector::Selector(Channel<digital> *In,Channel<analogic> *Out, vector<analogic> V){
	Input=In;
	AOutput=Out;
	resetable=false;
	Type=ANALOG;
	for(unsigned i=0;i<V.size();i++)
		AValues.push_back(V[i]);
	ActualPos=0;
	}
Selector::Selector(Channel<digital> *In,Channel<timer> *Out, vector<timer> V){
	Input=In;
	TOutput=Out;
	resetable=false;
	Type=TIME;
	for(unsigned i=0;i<V.size();i++){
		TValues.push_back(V[i]);
		}
	ActualPos=0;
	}
Selector::Selector(Channel<digital> *In,Channel<digital> *R,Channel<analogic> *Out, vector<analogic> V){
	Input=In;
	AOutput=Out;
	Reset=R;
	resetable=true;
	Type=ANALOG;
	for(unsigned i=0;i<V.size();i++)
		AValues.push_back(V[i]);
	ActualPos=0;
	}
Selector::Selector(Channel<digital> *In,Channel<digital> *R,Channel<timer> *Out, vector<timer> V){
	Input=In;
	TOutput=Out;
	Reset=R;
	resetable=true;
	Type=TIME;
	for(unsigned i=0;i<V.size();i++){
		TValues.push_back(V[i]);
		cout<<TValues[i]<<endl;
		}
	ActualPos=0;
	}
void Selector::Run(){
	if(Input->get_Value()&&(!Input->get_Value(1)))
		ActualPos++;
	if(Type==ANALOG)
		if(ActualPos>=AValues.size())
			ActualPos=0;
	if(Type==TIME)
		if(ActualPos>=TValues.size())
			ActualPos=0;
	if(Type==STRING)
		if(ActualPos>=SValues.size())
			ActualPos=0;
	if(resetable)
		if(Reset->get_Value()&&(!Reset->get_Value(1)))
			ActualPos=0;
	if(Type==ANALOG)
		AOutput->insert_Value(AValues[ActualPos]);
	if(Type==STRING)
		SOutput->insert_Value(SValues[ActualPos]);
	if(Type==TIME){
		TOutput->insert_Value(TValues[ActualPos]);
		}
	}

Counter::Counter(Channel<digital> *In, Channel<digital> *Out,unsigned C){
	Input=In;
	Output=Out;
	Count=C;
	ActualPos=0;
	Resetable=false;
	}
Counter::Counter(Channel<digital> *In,Channel<digital> *Rs, Channel<digital> *Out,unsigned C){
	Input=In;
	Output=Out;
	Reset=Rs;
	Count=C;
	ActualPos=0;
	Resetable=true;
	
	}
void Counter::Run(){
	if(Input->get_Value()&&(!Input->get_Value(1)))
		ActualPos++;
	if(ActualPos>=Count){
		ActualPos=0;
		Output->insert_Value(true);
		}
	else
		Output->insert_Value(false);
	if(Resetable)
		if(Reset->get_Value()&&Reset->get_Value(1))
			ActualPos=0;
	}

RiseDown::RiseDown(Channel<digital> *In,Channel<digital> *Out){
	FunctionalDescription_set("RiseDown Logic");
	Input=In;
	Output=Out;
	}
void RiseDown::Run(){
	if((Input->get_Value(1))&&(!Input->get_Value(0)))
		Output->insert_Value(true);
	else
		Output->insert_Value(false);
	}

RiseUp::RiseUp(Channel<digital> *In,Channel<digital> *Out){
	FunctionalDescription_set("RiseUp Logic");
	Input=In;
	Output=Out;
	}
void RiseUp::Run(){
	if((Input->get_Value(0))&&(!Input->get_Value(1)))
		Output->insert_Value(true);
	else
		Output->insert_Value(false);
	}
LockoutCounter::LockoutCounter(Channel<digital> *In, Channel<digital> *Out,unsigned C){
	Input=In;
	Output=Out;
	Count=C;
	ActualPos=0;
	Resetable=false;
	}

LockoutCounter::LockoutCounter(Channel<digital> *In,Channel<digital> *Rs, Channel<digital> *Out,unsigned C){
	Input=In;
	Output=Out;
	Reset=Rs;
	Count=C;
	ActualPos=0;
	Resetable=true;
	}
void LockoutCounter::Run(){
	if(Input->get_Value())
		ActualPos++;
	else
		if(ActualPos>0)
			ActualPos--;
	if(Resetable)
		if(Reset->get_Value())
			ActualPos=0;
	if(ActualPos>=Count){
		ActualPos=Count+1;
		Output->insert_Value(true);
		}
	else{
		Output->insert_Value(false);
		}
	}


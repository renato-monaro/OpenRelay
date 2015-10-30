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
#include "control.h"


using namespace orelay;
void Control::FunctionalDescription_set(string Desc){
	FunctionalDescription=Desc;
	}

string Control::FunctionalDescription_get(){
	return FunctionalDescription;
	}


void Control::SetClock(Channel<timer> *t){
	MasterClock=t;
	}

Fault::Fault(Channel<analogic> *In, Channel<analogic> *Freq,Channel<digital>*Out, Channel<digital>*Ctrl,vector<double>Ang, vector<double>P, vector<double>Q, int Ncicles){
		FunctionalDescription_set("Fault Control"); 
		Angles=Ang;
		ActivePower=P;
		ReactivePower=Q;
		Input=In;
		Frequency=Freq;
		Output=Out;
		Trigger=Ctrl;
		FaultCicles=Ncicles;
		Trigged=false;
		Waiting=false;
		Armed=false;
		Resting=true;
		k=0;
		i=0;
		j=0;
		FaultTime=orelay_gettime();
	}

void Fault::Run(){
	if((Trigger->get_Value(0))&&(!Trigger->get_Value(1))&&(!Trigged)&&(!Waiting)&&(!Armed)&&(!Resting)){
		Armed=true;
		}
	if(Armed){
		if((Input->get_Value()>=0)&&(Input->get_Value(1)<=0)){  //find zero cross
			FaultTime=MasterClock->get_Value(0)+(Input->get_Value(0)*(MasterClock->get_Value(1)-MasterClock->get_Value(0))/(Input->get_Value(1)-Input->get_Value(0)));
			FaultTime=FaultTime+((1E9/Frequency->get_Value())*(Angles[k]/360));//-(1E-3*ACQUISITION_TICK)-(0));
//			FaultTime=FaultTime+boost::posix_time::microsec((1E6/SYSTEM_FREQUENCY)*(Angles[k]/360));//-(1E-3*ACQUISITION_TICK)-(0));
			Waiting=true;
			Armed=false;
			k++;
			}
		}
	if(Waiting)
		if(FaultTime<MasterClock->get_Value()){ //Wait for the right momment to insert fault
			Trigged=true; //Fault!
			Waiting=false;
			FaultTime+=(FaultCicles*(1E9/Frequency->get_Value())+3424*4);
//			FaultTime+=boost::posix_time::microsec(FaultCicles*(1E6/SYSTEM_FREQUENCY)+3424*4);
			}
	if(Trigged){
		if(FaultTime<MasterClock->get_Value()){ //Wait for the right momment to remove fault
			Trigged=false;
			Resting=true;
			FaultTime+=10E9;
			}
		}
	if(Trigged)
		Output->insert_Value(true);
	else
		Output->insert_Value(false);
	if(Resting)
		{
		if(FaultTime<MasterClock->get_Value())
			{
			Resting=false;
			if(k>=Angles.size()){
				k=0;
				i++;
				}
			if(i>=ActivePower.size()){
				i=0;
				}
			cout<<"Ajust P:"<<ActivePower[i]<<" Q:"<<ReactivePower[i]<<" AND Press the Red Button ";
			cout<<"Fault Angle:"<<Angles[k]<<endl;
			}
		}
		
}

/*
PID::PID(Channel<analogic> *In, Channel<analogic>*Out, analogic Set, analogic Kp, analogic Ti, analogic Td){
	FunctionalDescription_set("PID Control"); 
	Dynamic=false;
	Saturable=false;
	T0=In->get_Sampling();
	q0=Kp*(1+Td/T0);
	q1w=-Kp*(1+2*Td/T0);
	if(Ti!=0)
		q1=-Kp*(1+2*Td/T0-T0/Ti);
	else //se for zero considerar tempo de integração infinito
		q1=q1w;
	q2=Kp*Td/T0;
	Input=In;
	Output=Out;
	SetPoint=Set;
	Action=0;
	kp=Kp;
	ti=Ti;
	kp=Kp;
//Error.set_capacity(3);
//Error.resize(3);
	#ifdef DEBUG
		cout<<"q0:"<<q0<<" q1:"<<q1<<" q2:"<<q2<<endl;
	#endif
	}

PID::PID(Channel<analogic> *In, Channel<analogic>*Out, double Set, double Kp, double Ti, double Td,double LowLimit, double UpLimit){
	Dynamic=false;
	FunctionalDescription_set("PID Control"); 
	Saturable=true;
	UpperLimit=UpLimit;
	LowerLimit=LowLimit;
	T0=In->get_Sampling();
	q0=Kp*(1+Td/T0);
	q1w=-Kp*(1+2*Td/T0);
	if(Ti>0.0)
		q1=-Kp*(1+2*Td/T0-T0/Ti);
	else //se for zero considerar tempo de integração infinito
		q1=q1w;
	q2=Kp*Td/T0;
	Input=In;
	Output=Out;
	SetPoint=Set;
	Action=0;
	kp=Kp;
	ti=Ti;
	kp=Kp;
//Error.set_capacity(3);
//Error.resize(3);
	#ifdef DEBUG
		cout<<"q0:"<<q0<<" q1:"<<q1<<" q2:"<<q2<<endl;
	#endif
	}

PID::PID(Channel<analogic> *In, Channel<analogic>*Out, Channel<analogic> *Set, analogic Kp, analogic Ti, analogic Td){
	FunctionalDescription_set("PID Control"); 
	Dynamic=true;
	DynamicSetPoint=Set;
	Saturable=false;
	T0=In->get_Sampling();
	q0=Kp*(1+Td/T0);
	q1w=-Kp*(1+2*Td/T0);
	if(Ti>0.0)
		q1=-Kp*(1+2*Td/T0-T0/Ti);
	else //se for zero considerar tempo de integração infinito
		q1=q1w;
	q2=Kp*Td/T0;
	Input=In;
	Output=Out;
	Action=0;
	kp=Kp;
	ti=Ti;
	kp=Kp;
//Error.set_capacity(3);
//Error.resize(3);
	#ifdef DEBUG
		cout<<"q0:"<<q0<<" q1:"<<q1<<" q2:"<<q2<<endl;
	#endif
	}

PID::PID(Channel<analogic> *In, Channel<analogic>*Out, Channel<analogic> *Set, double Kp, double Ti, double Td,double LowLimit, double UpLimit){
	FunctionalDescription_set("PID Control"); 
	Dynamic=true;
	DynamicSetPoint=Set;
	Saturable=true;
	UpperLimit=UpLimit;
	LowerLimit=LowLimit;
	T0=In->get_Sampling();
	q0=Kp*(1+Td/T0);
	q1w=-Kp*(1+2*Td/T0);
	kp=Kp;
	ti=Ti;
	kp=Kp;
//Error.set_capacity(3);
//Error.resize(3);
	if(Ti>0.0)
		q1=-Kp*(1+2*Td/T0-T0/Ti);
	else //se for zero considerar tempo de integração infinito
		q1=q1w;
	q2=Kp*Td/T0;
	Input=In;
	Output=Out;
	Action=0;
	#ifdef DEBUG
		cout<<"q0:"<<q0<<" q1:"<<q1<<" q2:"<<q2<<endl;
	#endif
	}


void PID::Run(){
	if(Dynamic)
		SetPoint=DynamicSetPoint->get_Value();
	//Error=0;

	Error_2=Error_1;
	Error_1=Error;
	Error=SetPoint-Input->get_Value();



		Action+=q0*Error+q1*Error_1+q2*Error_2;
	Action=kp*Error;
	if(Saturable&&(Action>UpperLimit))
		Action=UpperLimit;	
	if(Saturable&&(Action<LowerLimit))
		Action=LowerLimit;
	Output->insert_Value(Action);
	}*/

PI::PI(Channel<analogic> *In, Channel<analogic>*Out, double Set, double Kp, double Ti, double LowLimit, double UpLimit,bool Inv,double t0){
	FunctionalDescription_set("PI Control"); 
	Dynamic=false;
	Saturable=true;
	UpperLimit=UpLimit;
	LowerLimit=LowLimit;
	T0=t0;
	q0=Kp+(Kp/Ti)*(T0/2);
	q1=(Kp/Ti)*(T0/2)-Kp;
	Input=In;
	Output=Out;
	SetPoint=Set;
	Action=0;
	Invert=Inv;
//Error.set_capacity(3);
//Error.resize(3);
	#ifdef DEBUG
		cout<<"q0:"<<q0<<" q1:"<<q1<<endl;
	#endif
	}

PI::PI(Channel<analogic> *In, Channel<analogic>*Out, Channel<analogic>*Set, double Kp, double Ti, double LowLimit, double UpLimit,bool Inv,double t0){
	FunctionalDescription_set("PI Control"); 
	Dynamic=true;
	DynamicSetPoint=Set;
	Saturable=true;
	UpperLimit=UpLimit;
	LowerLimit=LowLimit;
	Invert=Inv;
	T0=t0;
	q0=Kp+(Kp/Ti)*(T0/2);
	q1=(Kp/Ti)*(T0/2)-Kp;
	Input=In;
	Output=Out;
	Action=0;
//Error.set_capacity(3);
//Error.resize(3);
	#ifdef DEBUG
		cout<<"q0:"<<q0<<" q1:"<<q1<<endl;
	#endif
	}

PI::PI(Channel<analogic> *In, Channel<analogic>*Out, double Set, double Kp, double Ti, double LowLimit, double UpLimit,double t0){
	FunctionalDescription_set("PI Control"); 
	Dynamic=false;
	Saturable=true;
	UpperLimit=UpLimit;
	LowerLimit=LowLimit;
	T0=t0;
	q0=Kp+(Kp/Ti)*(T0/2);
	q1=(Kp/Ti)*(T0/2)-Kp;
	Input=In;
	Output=Out;
	SetPoint=Set;
	Action=0;
	Invert=false;
//Error.set_capacity(3);
//Error.resize(3);
	#ifdef DEBUG
		cout<<"q0:"<<q0<<" q1:"<<q1<<endl;
	#endif
	}

PI::PI(Channel<analogic> *In, Channel<analogic>*Out, Channel<analogic>*Set, double Kp, double Ti, double LowLimit, double UpLimit,double t0){
	FunctionalDescription_set("PI Control"); 
	Dynamic=true;
	DynamicSetPoint=Set;
	Saturable=true;
	UpperLimit=UpLimit;
	LowerLimit=LowLimit;
	Invert=false;
	T0=t0;
	q0=Kp+(Kp/Ti)*(T0/2);
	q1=(Kp/Ti)*(T0/2)-Kp;
	Input=In;
	Output=Out;
	Action=0;
//Error.set_capacity(3);
//Error.resize(3);
	#ifdef DEBUG
		cout<<"q0:"<<q0<<" q1:"<<q1<<endl;
	#endif
	}

void PI::Run(){
	if(Dynamic)
		SetPoint=DynamicSetPoint->get_Value();
	Error_1=Error;
	Error=SetPoint-Input->get_Value();
	if(Invert)
		Error=-Error;
	Action+=q0*Error+q1*Error_1;
	//Action=kp*Error;
	if(Saturable&&(Action>UpperLimit))
		Action=UpperLimit;	
	if(Saturable&&(Action<LowerLimit))
		Action=LowerLimit;
	Output->insert_Value(Action);
	}


TF::TF(Channel<analogic> *In,Channel<analogic> *Out,vector<double> Ak,vector<double>Bk){
	Input=In;
	Output=Out;
	for(unsigned i=0;i<Ak.size();i++)
		A.push_back(Ak[i]);
	for(unsigned i=0;i<Bk.size();i++)
		B.push_back(Bk[i]);
	}
void TF::Run(){
	Out=0;
	for(unsigned k=0;k<B.size();k++)
		Out+=B[k]*Input->get_Value(k);
	for(unsigned k=1;k<A.size();k++)
		Out-=A[k]*Output->get_Value(k-1);
	if(A.size()>0)
		Out/=A[0];
	Output->insert_Value(Out);
	}

PWM::PWM(Channel<analogic> *In,Channel<digital> *Out1,Channel<digital> *Out2, double Amp, double Freq, double dT){
	Input=In;
	Output_1=Out1;
	Output_2=Out2;
	SwFreq=4*Amp*Freq*dT;
	Amplitude=Amp;
	Triang=0.0;
	Cres=true;

	Output_1_int=new Channel<digital>("Out1D",5);
	Output_2_int=new Channel<digital>("Out2D",5);
	T1=new Timer(Output_1_int,Output_1,4*1E9*dT);
	T2=new Timer(Output_2_int,Output_2,4*1E9*dT);

	T1->SetClock(MasterClock);
	T1->Prepare(1/dT);
	T2->SetClock(MasterClock);
	T2->Prepare(1/dT);

	}

void PWM::Run(){
	if(Triang>=Amplitude)
		Cres=false;
	if(Triang<=-Amplitude)
		Cres=true;
	if(Cres)
		Triang+=SwFreq;
	else
		Triang-=SwFreq;
	if(Triang<=Input->get_Value()){
		//Output_1_int->insert_Value(true);
		//Output_2_int->insert_Value(false);
		Output_1->insert_Value(true);
		Output_2->insert_Value(false);
		}
	else{
		//Output_1_int->insert_Value(false);
		//Output_2_int->insert_Value(true);
		Output_1->insert_Value(false);
		Output_2->insert_Value(true);

		}
	//T1->Run();
	//T2->Run();
	}



PWM2::PWM2(Channel<analogic> *Tri){
	Triang=Tri;
	}
	
PWM2::PWM2(Channel<analogic> *Tri,Channel<digital> *B){
	Triang=Tri;
	Blk=B;
	}
	
bool PWM2::Join(Channel<analogic> *Mod,Channel<digital> *Out1,Channel<digital> *Out2){
	vMod.push_back(Mod);
	vOutput_1.push_back(Out1);
	vOutput_2.push_back(Out2);
	return true;
	}
bool PWM2::Join(Channel<analogic> *Mod,Channel<digital> *Out1){
	return Join(Mod,Out1,NULL);
	}
void PWM2::Run(){

	bool Block=false;
	if(Blk!=NULL)
		Block=Blk->get_Value();
	for(unsigned k=0;k<vMod.size();k++){
		if(!Block){
			if(vMod[k]->get_Value()>=Triang->get_Value()){
				vOutput_1[k]->insert_Value(true);
				if(vOutput_2[k]!=NULL)
					vOutput_2[k]->insert_Value(false);
				}
			else{
				vOutput_1[k]->insert_Value(false);
				if(vOutput_2[k]!=NULL)
					vOutput_2[k]->insert_Value(true);
				}
			}
		else{
			vOutput_1[k]->insert_Value(false);
			vOutput_2[k]->insert_Value(false);
			}
		}
	}

TriangularWave::TriangularWave(Channel<analogic> *Out, double Amp,  double Freq, double dT){
	Output=Out;
	SwFreq=4*Amp*Freq*dT;
	Amplitude=Amp;
	Triang=0.0;
	Cres=true;
	}
void TriangularWave::Run(){
	if(Triang>=Amplitude)
		Cres=false;
	if(Triang<=-Amplitude)
		Cres=true;
	if(Cres)
		Triang+=SwFreq;
	else
		Triang-=SwFreq;
	Output->insert_Value(Triang);

	}


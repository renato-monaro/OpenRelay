/*
# Copyright (C) 2008-2015 Renato Machado Monaro, Vinicius de Cillo Moro
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
#include "protection.h"
using namespace orelay;
void Protection::SetClock(Channel<timer> *t){
	MasterClock=t;
	}

void Protection::FunctionalDescription_set(string Desc){
	FunctionalDescription=Desc;
	}

string Protection::FunctionalDescription_get(){
	return FunctionalDescription;
	}

void Protection::ANSICode_set(string Code){
	ANSICode=Code;
	}

string Protection::ANSICode_get(){
	return ANSICode;
	}


InstantaneousOverVoltage::InstantaneousOverVoltage(Channel<analogic>* In, Channel<digital>* Out, analogic Pick) 
	{
	FunctionalDescription_set("OverVoltage Protection"); 
	ANSICode_set("59");
	PickUp=Pick;
	Input_ANAL=In;
	////Sampling=In->get_Sampling();
	Input_Type=false;
	Output=Out;
	Comparation=OVER_MAGNITUDE;
	}

InstantaneousUnderVoltage::InstantaneousUnderVoltage(Channel<analogic>* In, Channel<digital>* Out, analogic Pick)
	{
	FunctionalDescription_set("UnderVoltage Protection");
	ANSICode_set("27");
	PickUp=Pick;
	Input_ANAL=In;
	////Sampling=In->get_Sampling();
	Input_Type=false;
	Output=Out;
	Comparation=UNDER_MAGNITUDE;
	}

InstantaneousOverCurrent::InstantaneousOverCurrent(Channel<analogic>* In, Channel<digital>* Out, analogic Pick)
	{
	FunctionalDescription_set("Instantaneous OverCurrent Protection");
	ANSICode_set("50");
	PickUp=Pick;
	Input_ANAL=In;
	////Sampling=In->get_Sampling();
	Input_Type=false;
	Output=Out;
	Comparation=OVER_MAGNITUDE;
	}

InstantaneousUnderCurrent::InstantaneousUnderCurrent(Channel<analogic>* In, Channel<digital>* Out, analogic Pick)
	{
	FunctionalDescription_set("Instantaneous UnderCurrent Protection");
	ANSICode_set("37");
	PickUp=Pick;
	Input_ANAL=In;
	////Sampling=In->get_Sampling();
	Input_Type=false;
	Output=Out;
	Comparation=UNDER_MAGNITUDE;
	}

InverseTimeOverCurrent::InverseTimeOverCurrent(Channel<analogic>* In, Channel<digital>* Out, analogic Pick, double Dial, unsigned Type){
	FunctionalDescription_set("Inverse Time OverCurrent Protection"); 	
	ANSICode_set("51");
	PickUp=Pick;
	Input_ANAL=In;
	////Sampling=In->get_Sampling();
	Input_Type=false;
	Output=Out;
	TDM=Dial;
	CurveType=Type;
	InverseIntegralTime=0;
	MakeCurve();
	Comparation=OVER_MAGNITUDE;
	Resetable=false;
	}

InverseTimeOverCurrent::InverseTimeOverCurrent(Channel<analogic>* In, Channel<digital>* Out, analogic Pick, double Dial, unsigned Type, bool Reset){
	FunctionalDescription_set("Inverse Time OverCurrent Protection"); 	
	ANSICode_set("51");
	PickUp=Pick;
	Input_ANAL=In;
	////Sampling=In->get_Sampling();
	Input_Type=false;
	Output=Out;
	TDM=Dial;
	CurveType=Type;
	InverseIntegralTime=0;
	MakeCurve();
	Comparation=OVER_MAGNITUDE;
	Resetable=Reset;
	}

InverseTimeUnderCurrent::InverseTimeUnderCurrent(Channel<analogic>* In, Channel<digital>* Out, analogic Pick, double Dial, unsigned Type){
	FunctionalDescription_set("Inverse Time UnderCurrent Protection"); 	
	ANSICode_set("37");
	PickUp=Pick;
	Input_ANAL=In;
	////Sampling=In->get_Sampling();
	Input_Type=false;
	Output=Out;
	TDM=Dial;
	CurveType=Type;
	InverseIntegralTime=0;
	MakeCurve();
	Comparation=UNDER_MAGNITUDE;
	Resetable=false;
	}

InverseTimeOverVoltage::InverseTimeOverVoltage(Channel<analogic>* In, Channel<digital>* Out, analogic Pick, double Dial, unsigned Type){
	FunctionalDescription_set("Inverse Time OverVoltage Protection"); 	
	ANSICode_set("59");
	PickUp=Pick;
	Input_ANAL=In;
	////Sampling=In->get_Sampling();
	Input_Type=false;
	Output=Out;
	TDM=Dial;
	CurveType=Type;
	InverseIntegralTime=0;
	MakeCurve();
	Comparation=OVER_MAGNITUDE;
	Resetable=false;
	}

InverseTimeUnderVoltage::InverseTimeUnderVoltage(Channel<analogic>* In, Channel<digital>* Out, analogic Pick, double Dial, unsigned Type){
	FunctionalDescription_set("Inverse Time UnderVoltage Protection"); 	
	ANSICode_set("27");
	PickUp=Pick;
	Input_ANAL=In;
	////Sampling=In->get_Sampling();
	Input_Type=false;
	Output=Out;
	TDM=Dial;
	CurveType=Type;
	InverseIntegralTime=0;
	MakeCurve();
	Comparation=UNDER_MAGNITUDE;
	Resetable=false;
	}

DefiniteTimeOverCurrent::DefiniteTimeOverCurrent(Channel<analogic>* In, Channel<digital>* Out, analogic Pick, double Dial){
	FunctionalDescription_set("Definite Time OverCurrent Protection"); 	
	ANSICode_set("51");
	PickUp=Pick;
	Input_ANAL=In;
	////Sampling=In->get_Sampling();
	Input_Type=false;
	Output=Out;
	DefiniteTime=Dial;
	DefiniteIntegralTime=0;
	Comparation=OVER_MAGNITUDE;
	}

DefiniteTimeUnderCurrent::DefiniteTimeUnderCurrent(Channel<analogic>* In, Channel<digital>* Out, analogic Pick, double Dial){
	FunctionalDescription_set("Definite Time UnderCurrent Protection"); 	
	ANSICode_set("37");
	PickUp=Pick;
	Input_ANAL=In;
	////Sampling=In->get_Sampling();
	Input_Type=false;
	Output=Out;
	DefiniteTime=Dial;
	DefiniteIntegralTime=0;
	Comparation=UNDER_MAGNITUDE;
	}

DefiniteTimeOverVoltage::DefiniteTimeOverVoltage(Channel<analogic>* In, Channel<digital>* Out, analogic Pick, double Dial){
	FunctionalDescription_set("Definite Time OverVoltage Protection"); 	
	ANSICode_set("59");
	PickUp=Pick;
	Input_ANAL=In;
	////Sampling=In->get_Sampling();
	Input_Type=false;
	Output=Out;
	DefiniteTime=Dial;
	DefiniteIntegralTime=0;
	Comparation=OVER_MAGNITUDE;
	}

DefiniteTimeUnderVoltage::DefiniteTimeUnderVoltage(Channel<analogic>* In, Channel<digital>* Out, analogic Pick, double Dial){
	FunctionalDescription_set("Definite Time UnderVoltage Protection"); 	
	ANSICode_set("27");
	PickUp=Pick;
	Input_ANAL=In;
	////Sampling=In->get_Sampling();
	Input_Type=false;
	Output=Out;
	DefiniteTime=Dial;
	DefiniteIntegralTime=0;
	Comparation=UNDER_MAGNITUDE;
	}



//***

InstantaneousOverVoltage::InstantaneousOverVoltage(Channel<Complex>* In, Channel<digital>* Out, analogic Pick) 
	{
	FunctionalDescription_set("OverVoltage Protection"); 
	ANSICode_set("59");
	PickUp=Pick;
	Input_CPLX=In;
	////Sampling=In->get_Sampling();
	Input_Type=true;
	Output=Out;
	Comparation=OVER_MAGNITUDE;
	}

InstantaneousUnderVoltage::InstantaneousUnderVoltage(Channel<Complex>* In, Channel<digital>* Out, analogic Pick)
	{
	FunctionalDescription_set("UnderVoltage Protection");
	ANSICode_set("27");
	PickUp=Pick;
	Input_CPLX=In;
	////Sampling=In->get_Sampling();
	Input_Type=true;
	Output=Out;
	Comparation=UNDER_MAGNITUDE;
	}

InstantaneousOverCurrent::InstantaneousOverCurrent(Channel<Complex>* In, Channel<digital>* Out, analogic Pick)
	{
	FunctionalDescription_set("Instantaneous OverCurrent Protection");
	ANSICode_set("50");
	PickUp=Pick;
	Input_CPLX=In;
	////Sampling=In->get_Sampling();
	Input_Type=true;
	Output=Out;
	Comparation=OVER_MAGNITUDE;
	}

InstantaneousUnderCurrent::InstantaneousUnderCurrent(Channel<Complex>* In, Channel<digital>* Out, analogic Pick)
	{
	FunctionalDescription_set("Instantaneous UnderCurrent Protection");
	ANSICode_set("37");
	PickUp=Pick;
	Input_CPLX=In;
	////Sampling=In->get_Sampling();
	Input_Type=true;
	Output=Out;
	Comparation=UNDER_MAGNITUDE;
	}

InverseTimeOverCurrent::InverseTimeOverCurrent(Channel<Complex>* In, Channel<digital>* Out, analogic Pick, double Dial, unsigned Type){
	FunctionalDescription_set("Inverse Time OverCurrent Protection"); 	
	ANSICode_set("51");
	PickUp=Pick;
	Input_CPLX=In;
	//Sampling=In->get_Sampling();
	Input_Type=true;
	Output=Out;
	TDM=Dial;
	CurveType=Type;
	InverseIntegralTime=0;
	MakeCurve();
	Comparation=OVER_MAGNITUDE;
	Resetable=false;
	}

InverseTimeOverCurrent::InverseTimeOverCurrent(Channel<Complex>* In, Channel<digital>* Out, analogic Pick, double Dial, unsigned Type, bool Reset){
	FunctionalDescription_set("Inverse Time OverCurrent Protection"); 	
	ANSICode_set("51");
	PickUp=Pick;
	Input_CPLX=In;
	//Sampling=In->get_Sampling();
	Input_Type=true;
	Output=Out;
	TDM=Dial;
	CurveType=Type;
	InverseIntegralTime=0;
	MakeCurve();
	Comparation=OVER_MAGNITUDE;
	Resetable=Reset;
	}

InverseTimeUnderCurrent::InverseTimeUnderCurrent(Channel<Complex>* In, Channel<digital>* Out, analogic Pick, double Dial, unsigned Type){
	FunctionalDescription_set("Inverse Time UnderCurrent Protection"); 	
	ANSICode_set("37");
	PickUp=Pick;
	Input_CPLX=In;
	//Sampling=In->get_Sampling();
	Input_Type=true;
	Output=Out;
	TDM=Dial;
	CurveType=Type;
	InverseIntegralTime=0;
	MakeCurve();
	Comparation=UNDER_MAGNITUDE;
	Resetable=false;
	}

InverseTimeOverVoltage::InverseTimeOverVoltage(Channel<Complex>* In, Channel<digital>* Out, analogic Pick, double Dial, unsigned Type){
	FunctionalDescription_set("Inverse Time OverVoltage Protection"); 	
	ANSICode_set("59");
	PickUp=Pick;
	Input_CPLX=In;
	//Sampling=In->get_Sampling();
	Input_Type=true;
	Output=Out;
	TDM=Dial;
	CurveType=Type;
	InverseIntegralTime=0;
	MakeCurve();
	Comparation=OVER_MAGNITUDE;
	Resetable=false;
	}

InverseTimeUnderVoltage::InverseTimeUnderVoltage(Channel<Complex>* In, Channel<digital>* Out, analogic Pick, double Dial, unsigned Type){
	FunctionalDescription_set("Inverse Time UnderVoltage Protection"); 	
	ANSICode_set("27");
	PickUp=Pick;
	Input_CPLX=In;
	//Sampling=In->get_Sampling();
	Input_Type=true;
	Output=Out;
	TDM=Dial;
	CurveType=Type;
	InverseIntegralTime=0;
	MakeCurve();
	Comparation=UNDER_MAGNITUDE;
	Resetable=false;
	}

DefiniteTimeOverCurrent::DefiniteTimeOverCurrent(Channel<Complex>* In, Channel<digital>* Out, analogic Pick, double Dial){
	FunctionalDescription_set("Definite Time OverCurrent Protection"); 	
	ANSICode_set("51");
	PickUp=Pick;
	Input_CPLX=In;
	//Sampling=In->get_Sampling();
	Input_Type=true;
	Output=Out;
	DefiniteTime=Dial;
	DefiniteIntegralTime=0;
	Comparation=OVER_MAGNITUDE;
	}

DefiniteTimeUnderCurrent::DefiniteTimeUnderCurrent(Channel<Complex>* In, Channel<digital>* Out, analogic Pick, double Dial){
	FunctionalDescription_set("Definite Time UnderCurrent Protection"); 	
	ANSICode_set("37");
	PickUp=Pick;
	Input_CPLX=In;
	//Sampling=In->get_Sampling();
	Input_Type=true;
	Output=Out;
	DefiniteTime=Dial;
	DefiniteIntegralTime=0;
	Comparation=UNDER_MAGNITUDE;
	}

DefiniteTimeOverVoltage::DefiniteTimeOverVoltage(Channel<Complex>* In, Channel<digital>* Out, analogic Pick, double Dial){
	FunctionalDescription_set("Definite Time OverVoltage Protection"); 	
	ANSICode_set("59");
	PickUp=Pick;
	Input_CPLX=In;
	//Sampling=In->get_Sampling();
	Input_Type=true;
	Output=Out;
	DefiniteTime=Dial;
	DefiniteIntegralTime=0;
	Comparation=OVER_MAGNITUDE;
	}

DefiniteTimeUnderVoltage::DefiniteTimeUnderVoltage(Channel<Complex>* In, Channel<digital>* Out, analogic Pick, double Dial){
	FunctionalDescription_set("Definite Time UnderVoltage Protection"); 	
	ANSICode_set("27");
	PickUp=Pick;
	Input_CPLX=In;
	//Sampling=In->get_Sampling();
	Input_Type=true;
	Output=Out;
	DefiniteTime=Dial;
	DefiniteIntegralTime=0;
	Comparation=UNDER_MAGNITUDE;
	}

void InstantaneousMagnitude::Run(){ 
	dt=MasterClock->Sub(0,1);
	if(Input_Type)
		Input=abs(Input_CPLX->get_Value());//Complex
	else
		Input=Input_ANAL->get_Value(); //Analogic  
	if(Comparation==OVER_MAGNITUDE){
		if(Input>=PickUp){   
 			Output->insert_Value(true);
			}
		else
			Output->insert_Value(false);
		}
	if(Comparation==UNDER_MAGNITUDE){
		if(Input<=PickUp){
			Output->insert_Value(true);  
			}
		else
			Output->insert_Value(false);
		} 
	}

bool InstantaneousMagnitude::Prepare(float SamplingFreq){ 
	dt=1/SamplingFreq;
	return true;
	}

void DefiniteTimeMagnitude::Run(){  
	if(Input_Type)
		Input=abs(Input_CPLX->get_Value());//Complex
	else
		Input=Input_ANAL->get_Value(); //Analogic   
	if(Comparation==OVER_MAGNITUDE){
		if(Input>=PickUp){
			DefiniteIntegralTime+=dt;
			if (DefiniteIntegralTime>=DefiniteTime)
				Output->insert_Value(true);
			else 
				Output->insert_Value(false); 
			}
		else{
			DefiniteIntegralTime=0;
			Output->insert_Value(false);
			}
		}
	if(Comparation==UNDER_MAGNITUDE){
		if(Input<=PickUp){
			DefiniteIntegralTime+=dt;
			if (DefiniteIntegralTime>=DefiniteTime)
				Output->insert_Value(true);
			else 
				Output->insert_Value(false); 
			}
		else{
			DefiniteIntegralTime=0;
			Output->insert_Value(false);
			}
		} 
	}   
bool DefiniteTimeMagnitude::Prepare(float SamplingFreq){ 
	dt=1/SamplingFreq;
	return true;
	}
              
void InverseTimeMagnitude::Run(){
	double t,tr;
	if(Input_Type)
		Input=abs(Input_CPLX->get_Value());//Complex
	else
		Input=Input_ANAL->get_Value(); //Analogic  
	if(Comparation==OVER_MAGNITUDE){
		if (Input>=PickUp){
			if((CurveType==IEEE_EI)||(CurveType==IEEE_VI)||(CurveType==IEEE_MI))
				t=TDM*(A/((pow(Input/PickUp,P)-1))+B);
			if((CurveType==IEC_A)||(CurveType==IEC_B)||(CurveType==IEC_C)||(CurveType==IEC_SI))
				t=TDM*(K/(pow(Input/PickUp,E)-1));
			if((CurveType==IAC_EI)||(CurveType==IAC_VI)||(CurveType==IAC_I)||(CurveType==IAC_SI))
				t=TDM*(A+B/((Input/PickUp)-C)+D/pow(Input/PickUp-C,2)+E/pow(Input/PickUp-C,3));
			if(CurveType==I2T)
				t=TDM*(100/pow(Input/PickUp,2));
			InverseIntegralTime+=dt/t;
			if (InverseIntegralTime>1)
				Output->insert_Value(true);
			else 
				Output->insert_Value(false); 
			}	
		else {
			Output->insert_Value(false);
			if(Resetable&&(InverseIntegralTime>0.0)){
				if((CurveType==IEEE_EI)||(CurveType==IEEE_VI)||(CurveType==IEEE_MI)||(CurveType==IEC_A)||(CurveType==IEC_B)||(CurveType==IEC_C)||(CurveType==IEC_SI)||(CurveType==IAC_EI)||(CurveType==IAC_VI)||(CurveType==IAC_I)||(CurveType==IAC_SI))
					tr=TDM*(Tr/(1-pow(Input/PickUp,2)));		
				if(CurveType==I2T)
					tr=TDM*(100/pow(Input/PickUp,-2));
				InverseIntegralTime-=dt/tr;
				}
			else{
				InverseIntegralTime=0; 
				}
			}
		}
	if(Comparation==UNDER_MAGNITUDE){
		if (Input<=PickUp){
			if((CurveType==IEEE_EI)||(CurveType==IEEE_VI)||(CurveType==IEEE_MI))
				t=TDM*(A/((1-pow(Input/PickUp,P)))+B);
				//t=TDM*(A/(pow(Input/PickUp,P)-1)+B);
		
			if((CurveType==IEC_A)||(CurveType==IEC_B)||(CurveType==IEC_C)||(CurveType==IEC_SI))
				t=TDM*(K/(1-pow(Input/PickUp,E)));
				//t=TDM*(K/(pow(Input/PickUp,E)-1));

			if((CurveType==IAC_EI)||(CurveType==IAC_VI)||(CurveType==IAC_I)||(CurveType==IAC_SI))
				t=TDM*(A+B/(C-(Input/PickUp))+D/pow(C-Input/PickUp,2)+E/pow(C-Input/PickUp,3));
				//t=TDM*(A+B/(Input/PickUp-C)+D/pow(Input/PickUp-C,2)+E/pow(Input/PickUp-C,3));
			InverseIntegralTime+=dt/t;
			if (InverseIntegralTime>1)
				Output->insert_Value(true);
			else 
				Output->insert_Value(false); 
			}
		else {
			Output->insert_Value(false);
			if(Resetable&&(InverseIntegralTime>0.0)){
				if((CurveType==IEEE_EI)||(CurveType==IEEE_VI)||(CurveType==IEEE_MI)||(CurveType==IEC_A)||(CurveType==IEC_B)||(CurveType==IEC_C)||(CurveType==IEC_SI)||(CurveType==IAC_EI)||(CurveType==IAC_VI)||(CurveType==IAC_I)||(CurveType==IAC_SI))
					tr=TDM*(Tr/(1-pow(Input/PickUp,2)));		
				if(CurveType==I2T)
					tr=TDM*(100/pow(Input/PickUp,-2));
				InverseIntegralTime-=dt/tr;
				}
			else{
				InverseIntegralTime=0; 
				}
			}
		}
	#ifdef DEBUG
		if(Output->get_Value()&&(!Output->get_Value(1)))
			cout<<"TimedOverMagnitude"<<endl;
	#endif
}

bool InverseTimeMagnitude::Prepare(float SamplingFreq){ 
	dt=1/SamplingFreq;
	return true;
	}

void InverseTimeMagnitude::MakeCurve(){
	switch(CurveType){
		case IEEE_EI:
			A=28.2; B=0.1217; P=2.0; Tr=29.1;
		break;

		case IEEE_VI:
			A=19.61; B=0.491; P=2.0; Tr=21.6;
		break;
		
		case IEEE_MI:
			A=0.0515; B=0.1140; P=0.02; Tr=4.85;
		break;

		case IEC_A:
			K=0.14; E=0.02; Tr=9.7;
		break;

		case IEC_B:
			K=13.5; E=1.0; Tr=43.2;
		break;
		
		case IEC_C:
			K=80.0; E=2.0; Tr=58.2;
		break;

		case IEC_SI:
			K=0.05; E=0.04; Tr=0.5;
		break;

		case IAC_EI:
			A=0.004; B=0.6379; C=0.62; D=1.7872; E=0.2461; Tr=6.008;
		break;

		case IAC_VI:
			A=0.09; B=0.7955; C=0.1; D=-1.2885; E=7.9586; Tr=4.678;
		break;
		
		case IAC_I:
			A=0.2078; B=0.863; C=0.8; D=-0.418; E=0.1947; Tr=0.99;
		break;

		case IAC_SI:
			A=0.0428; B=0.0609; C=0.62; D=-0.001; E=0.0221; Tr=0.222;
		break;
		}	
	}

Differential::Differential(Channel<Complex>* In1, Channel<Complex>* In2, Channel<digital>* Out, double IN,double P1, double B1, double S1, double B2, double S2,double P2){
	FunctionalDescription_set("Differential Protection"); 
 
	ANSICode_set("87");
	PickUp1=P1;
	PickUp2=P2;
	Slope1=S1;
	Slope2=S2;
	BasePoint1=B1;
	BasePoint2=B2;
	Input1=In1;
	Input2=In2;
	Output=Out;
	INorm=IN;
	a1=Slope1/100;
	b1=-a1*BasePoint1;
	Break1=(PickUp1-b1)/a1;
	a2=Slope2/100;
	b2=-a2*BasePoint2;
	Break2=(b2-b1)/(a1-a2);
	Break3=(PickUp2-b2)/a2;
	#ifdef DEBUG
		cout<<"a1:"<<a1<<" b1:"<<b1<<" a2:"<<a2<<" b2:"<<b2;
		cout<<"Break1:"<<Break1<<" Break2:"<<Break2<<" Break3:"<<Break3<<endl;
		PrintCurve();
	#endif
	}

bool Differential::Prepare(float SamplingFreq){ 
	return true;
	}

void Differential::Run(){
	//cout<<Input1->get_Value()<<"\t"<<Input2->get_Value()<<endl;
	bool Cond=false;
	unsigned type=0;
	Ir=(abs(Input1->get_Value())+abs(Input2->get_Value()))/INorm;
	Io=(abs(Input1->get_Value()+Input2->get_Value()))/INorm;
	//	cout<<" Ir:"<<Ir<<"\tIo:"<<Io<<endl;
	if((Ir<Break1)&&(Io>PickUp1)){
		Cond=true;
		type=1;
		//cout<<"Condition PickUp1";
		}
	if((Break2>Ir)&&(Ir>Break1)&&(Io>(a1*Ir+b1))){
		Cond=true;
		type=2;
		//cout<<"Condition Slope1";
		}
	if((Break3>Ir)&&(Ir>Break2)&&(Io>(a2*Ir+b2))){
		Cond=true;
		type=3;
		//cout<<"Condition Slope2";
		}
	if((Ir>Break3)&&(Io>PickUp2)){
		Cond=true;
		type=4;
		//cout<<"Condition PickUp2";
		}
	Output->insert_Value(Cond);
	#ifdef DEBUG
	if(Output->get_Value()&&(!Output->get_Value(1)))
		{
		cout<<"Differential Trip ";
		switch(type){
			case 1: 
				cout<<"Condition PickUp1";
			break;
			case 2: 
				cout<<"Condition Slope1";
			break;
			case 3: 
				cout<<"Condition Slope2";
			break;
			case 4: 
				cout<<"Condition PickUp2";
			break;
			default:
				cout<<"Not Determined";
				
			}
		cout<<" Ir:"<<Ir<<"\tIo:"<<Io<<endl;
		}
	#endif
	}
void Differential::PrintCurve(){
	cout<<"Ir\tIo\n";
	for(unsigned k=0;k<100;k++){
		Ir=k*Break3*2/100;
		if(Ir<=Break1)
			cout<<(Ir)<<"\t"<<(PickUp1)<<endl;
		if((Break2>=Ir)&&(Ir>Break1))
			cout<<(Ir)<<"\t"<<(a1*Ir+b1)<<endl;	
		if((Break3>=Ir)&&(Ir>Break2))
			cout<<(Ir)<<"\t"<<(a2*Ir+b2)<<endl;
		if(Ir>=Break3)
			cout<<(Ir)<<"\t"<<(PickUp2)<<endl;
		}
	}


DifferentialSpline::DifferentialSpline(Channel<Complex>* In1, Channel<Complex>* In2, Channel<digital>* Out, double In, double P1, double B1, double S1, double B2, double S2){
	FunctionalDescription_set("Differential Protection"); 
 
	ANSICode_set("87");
	PickUp1=P1;
	Slope1=S1;
	Slope2=S2;
	BreakPoint1=B1;
	BreakPoint2=B2;
	Input1=In1;
	Input2=In2;
	Output=Out;
	INorm=In;
	a1=Slope1/100;
	b1=0;
	Break1=(PickUp1-b1)/a1;
	a2=Slope2/100;
	b2=0;
	Break2=BreakPoint1;
	Break3=BreakPoint2;
	p0=a1*Break2+b1;
	m0=a1;
	p1=a2*Break3+b2;
	m1=a2;
	#ifdef DEBUG
		cout<<"a1:"<<a1<<" b1:"<<b1<<" a2:"<<a2<<" b2:"<<b2 <<" p0:"<<p0<<" p1:"<<p1<<" m0:"<<m0<<" m1:"<<m1;
		cout<<"Break1:"<<Break1<<" Break2:"<<Break2<<" Break3:"<<Break3<<endl;
		PrintCurve();
	#endif
	}

bool DifferentialSpline::Prepare(float SamplingFreq){ 
	return true;
	}

void DifferentialSpline::Run(){
	//cout<<Input1->get_Value()<<"\t"<<Input2->get_Value()<<endl;
	bool Cond=false;
	unsigned type=0;
	Io=(abs(Input1->get_Value()+Input2->get_Value()))/INorm;
	Ir=(abs(Input1->get_Value())+abs(Input2->get_Value()))/INorm;
	//cout<<"Io:"<<Io<<" Ir:"<<Ir<<endl;  
	if((Ir<Break1)&&(Io>PickUp1)){
		Cond=true;
		type=1;
		//cout<<"Condition PickUp1";
		}
	if((Break2>Ir)&&(Ir>Break1)&&(Io>(a1*Ir+b1))){
		Cond=true;
		type=2;
		//cout<<"Condition Slope1";
		}
	if((Break3>Ir)&&(Ir>Break2)&&(Io>Spline())){
		Cond=true;
		type=3;
		//cout<<"Condition Spline";
		}
	if((Ir>Break3)&&(Io>(a2*Ir+b2))){
		Cond=true;
		type=4;
		//cout<<"Condition Slope2";
		}
	Output->insert_Value(Cond);
	#ifdef DEBUG
	if(Output->get_Value()&&(!Output->get_Value(1)))
		{
		cout<<"Differential Trip ";
		switch(type){
			case 1: 
				cout<<"Condition PickUp";
			break;
			case 2: 
				cout<<"Condition Slope1";
			break;
			case 3: 
				cout<<"Condition Spline";
			break;
			case 4: 
				cout<<"Condition Slope2";
			break;
			default:
				cout<<"Not Determined";
				
			}
		cout<<" Ir:"<<Ir<<"\tIo:"<<Io<<endl;
		}
	#endif
	}

double DifferentialSpline::Spline(){
	t=(Ir-Break2)/(Break3-Break2);
	return ((2*pow(t,3)-3*pow(t,2)+1)*p0+(pow(t,3)-2*pow(t,2)+t)*(Break3-Break2)*m0+(-2*pow(t,3)+3*pow(t,2))*p1+(pow(t,3)-pow(t,2))*(Break3-Break2)*m1);
	}


void DifferentialSpline::PrintCurve(){
	cout<<"Ir\tIo\n";
	for(unsigned k=0;k<100;k++){
		Ir=k*Break3*2/100;
		if(Ir<=Break1)
			cout<<(Ir)<<"\t"<<(PickUp1)<<endl;
		if((Break2>=Ir)&&(Ir>Break1))
			cout<<(Ir)<<"\t"<<(a1*Ir+b1)<<endl;	
		if((Break3>=Ir)&&(Ir>Break2))
			cout<<(Ir)<<"\t"<<(Spline())<<endl;
		if(Ir>=Break3)
			cout<<(Ir)<<"\t"<<(a2*Ir+b2)<<endl;
		}
	}

DefiniteTimeUnderFrequency::DefiniteTimeUnderFrequency(Channel<analogic>* In, Channel<digital>* Out, analogic Pick, double Dial){
	FunctionalDescription_set("Definite Time UnderFrequency Protection"); 	
	ANSICode_set("87");
	PickUp=Pick;
	Input_ANAL=In;
	////Sampling=In->get_Sampling();
	Input_Type=false;
	Output=Out;
	DefiniteTime=Dial;
	DefiniteIntegralTime=0;
	Comparation=UNDER_MAGNITUDE;
	}

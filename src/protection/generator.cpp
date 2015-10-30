/*
# Copyright (C) 2015 Rodrigo Bataglioli
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
#include "generator.h"
using namespace orelay;
IDesbalanceadasG30::IDesbalanceadasG30(Channel<Complex>* In1, Channel<Complex>* In2, Channel<Complex>* In3, Channel<digital>* Out, Channel<analogic>* Out2, double In, double Tmin, double Tmax, double P, double constK){
	FunctionalDescription_set("Differential Protection");

	ANSICode_set("46");
	PickUp=P;
	tmin=Tmin;
	tmax=Tmax;
	K=constK;
	Input1=In1;
	Input2=In2;
	Input3=In3;
	Output=Out;
	Output2=Out2;
	INorm=In;
	IntegralTime=0;
	}

bool IDesbalanceadasG30::Prepare(float SamplingFreq){
    dt=1/SamplingFreq;
	return true;
	}

void IDesbalanceadasG30::Run(){
	bool Cond=false;
	double Iar, Iai, Ibr, Ibi, Icr, Ici, Ibrlinha, Ibilinha, Icrlinha, Icilinha, I2r, I2i;
	Iar=real(Input1->get_Value());
	Iai=imag(Input1->get_Value());
	Ibr=real(Input2->get_Value());
	Ibi=imag(Input2->get_Value());
	Icr=real(Input3->get_Value());
	Ici=imag(Input3->get_Value());
    Ibrlinha=(-0.5)*Ibr-(-sqrt(3)/2)*Ibi;
    Ibilinha=(-0.5)*Ibi+(-sqrt(3)/2)*Ibr;
    Icrlinha=(-0.5)*Icr-(sqrt(3)/2)*Ici;
    Icilinha=(-0.5)*Ici+(sqrt(3)/2)*Icr;
    I2r=(Iar+Ibrlinha+Icrlinha)/3;
	I2i=(Iai+Ibilinha+Icilinha)/3;
	I2calc=sqrt(pow(I2r,2)+pow(I2i,2));
	t=K/pow(I2calc/(INorm),2);
	if(t<tmin){
       t=tmin;
	}
	if(t>tmax){
       t=tmax;
	}
    if(I2calc>=(PickUp*INorm)){
       if(IntegralTime>=t){
                Cond=true;
            }
            IntegralTime+=dt;
		}
    else{
        if(IntegralTime>0){
            IntegralTime-=dt;
        }
    }

	Output->insert_Value(Cond);
	Output2->insert_Value(I2calc);
	}
	
DifferentialG30::DifferentialG30(Channel<Complex>* In1, Channel<Complex>* In2, Channel<digital>* Out, double In, double P1, double B1, double S1, double B2, double S2){
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
	}

bool DifferentialG30::Prepare(float SamplingFreq){
	return true;
	}

void DifferentialG30::Run(){
	bool Cond=false;
	unsigned type=0;
	Io=(abs(Input1->get_Value()+Input2->get_Value()))/INorm;
	Ir=(max(abs(Input1->get_Value()),abs(Input2->get_Value())))/INorm;
	if((Ir<Break1)&&(Io>PickUp1)){
		Cond=true;
		type=1;
		}
	if((Break2>Ir)&&(Ir>Break1)&&(Io>(a1*Ir+b1))){
		Cond=true;
		type=2;
		}
	if((Break3>Ir)&&(Ir>Break2)&&(Io>Spline())){
		Cond=true;
		type=3;
		}
	if((Ir>Break3)&&(Io>(a2*Ir+b2))){
		Cond=true;
		type=4;
		}
	Output->insert_Value(Cond);
	}

double DifferentialG30::Spline(){
	t=(Ir-Break2)/(Break3-Break2);
	return ((2*pow(t,3)-3*pow(t,2)+1)*p0+(pow(t,3)-2*pow(t,2)+t)*(Break3-Break2)*m0+(-2*pow(t,3)+3*pow(t,2))*p1+(pow(t,3)-pow(t,2))*(Break3-Break2)*m1);
	}
	
	
MotorizacaoG30::MotorizacaoG30(Channel<analogic>* In1, Channel<analogic>* In2, Channel<analogic>* In3, Channel<analogic>* In4, Channel<analogic>* In5, Channel<analogic>* In6, Channel<analogic>* Out1, Channel<digital>* Out2, double Snom, double TPprim, double TPsec, double TCsec, double P, double D){
	FunctionalDescription_set("Loss of Field Protection");

	ANSICode_set("32");
	Pickup=P;
    Temp=D;
    SNorm=Snom;
    RTP=TPprim/TPsec;
    RTC=1/TCsec;
	Input1=In1;
	Input2=In2;
	Input3=In3;
	Input4=In4;
	Input5=In5;
	Input6=In6;
	Output1=Out1;
	Output2=Out2;
	IntegralTime=0;
	}

bool MotorizacaoG30::Prepare(float SamplingFreq){
	dt=1/SamplingFreq;
	return true;
	}

void MotorizacaoG30::Run(){
	bool Cond=false;
	PotCalc = (Input1->get_Value())*RTP*(Input4->get_Value())/RTC+(Input2->get_Value())*RTP*(Input5->get_Value())/RTC+(Input3->get_Value())*RTP*(Input6->get_Value())/RTC;
	if(PotCalc<=(Pickup*SNorm*pow(10,6))){
        if(IntegralTime>=Temp){
                Cond=true;
            }
            IntegralTime+=dt;
	}
	else{
        if(IntegralTime>0){
            IntegralTime-=dt;
        }
	}
	Output1->insert_Value(PotCalc);
	Output2->insert_Value(Cond);
	}
	
PerdaCampoG30::PerdaCampoG30(Channel<Complex>* In1, Channel<Complex>* In2, Channel<digital>* Out, double C1, double R1, double C2, double R2, double D){
	FunctionalDescription_set("Loss of Field Protection");

	ANSICode_set("40");
	Centro1=C1;
	Raio1=R1;
	Centro2=C2;
	Raio2=C2;
    Temp=D;
	Input1=In1;
	Input2=In2;
	Output=Out;
	IntegralTime=0;
	}

bool PerdaCampoG30::Prepare(float SamplingFreq){
	dt=1/SamplingFreq;
	return true;
	}

void PerdaCampoG30::Run(){
	bool Cond=false;
	double a, b, c, d, xa, ya, d1, d2;
	a = real(Input1->get_Value());
	b = imag(Input1->get_Value());
	c = real(Input2->get_Value());
	d = imag(Input2->get_Value());
	xa = (a*c+b*d)/(pow(c,2)+pow(d,2));
	ya = (b*c-a*d)/(pow(c,2)+pow(d,2));
	d1 = sqrt(pow(0-xa,2)+pow(-Centro1-ya,2));
	d2 = sqrt(pow(0-xa,2)+pow(-Centro2-ya,2));
	if(d1<Raio1){
        Cond=true;
		}
	if(d2<Raio2){
            if(IntegralTime>=Temp){
                Cond=true;
            }
            IntegralTime+=dt;
		}
    else{
        if(IntegralTime>0){
            IntegralTime-=dt;
        }
    }
	Output->insert_Value(Cond);
	}
	
SobreexcitacaoG30::SobreexcitacaoG30(Channel<analogic>* In1, Channel<Complex>* In2, Channel<analogic>* In3, Channel<Complex>* In4, Channel<analogic>* In5, Channel<Complex>* In6, Channel<digital>* Out, double Vn, double Fn, double P, double D){
	FunctionalDescription_set("Overexcitation Protection");

	ANSICode_set("24");
	PickUp=P;
    Temp=D;
	Input1=In1;
	Input2=In2;
    Input3=In3;
	Input4=In4;
    Input5=In5;
	Input6=In6;
	Output=Out;
	VNorm=Vn;
	FNorm=Fn;
	IntegralTime=0;
	}

bool SobreexcitacaoG30::Prepare(float SamplingFreq){
	dt=1/SamplingFreq;
	return true;
	}

void SobreexcitacaoG30::Run(){
	bool Cond=false;
	double aux1, aux2, aux3;
	aux1 = (abs(Input2->get_Value())/VNorm)/((Input1->get_Value())/FNorm);
	aux2 = (abs(Input4->get_Value())/VNorm)/((Input3->get_Value())/FNorm);
	aux3 = (abs(Input6->get_Value())/VNorm)/((Input5->get_Value())/FNorm);
	VHz = (aux1+aux2+aux3)/3;
	if(VHz>PickUp){
            if(IntegralTime>=Temp){
                Cond=true;
            }
            IntegralTime+=dt;
		}
    else{
        if(IntegralTime>0){
            IntegralTime-=dt;
        }
    }
	Output->insert_Value(Cond);
	}
SobrfrequenciaG30::SobrfrequenciaG30(Channel<analogic>* In1, Channel<Complex>* In2, Channel<analogic>* In3, Channel<Complex>* In4, Channel<analogic>* In5, Channel<Complex>* In6, Channel<digital>* Out, double Vn, double MinV, double P, double D){
	FunctionalDescription_set("Overexcitation Protection");

	ANSICode_set("81");
	PickUp=P;
    Temp=D;
	Input1=In1;
	Input2=In2;
    Input3=In3;
	Input4=In4;
    Input5=In5;
	Input6=In6;
	Output=Out;
	VNorm=Vn;
	VMin=MinV;
	IntegralTime=0;
	}

bool SobrfrequenciaG30::Prepare(float SamplingFreq){
	dt=1/SamplingFreq;
	return true;
	}

void SobrfrequenciaG30::Run(){
	bool Cond=false;
	double aux1, aux2, aux3;
	aux1 = (abs(Input2->get_Value())/VNorm);
	aux2 = (abs(Input4->get_Value())/VNorm);
	aux3 = (abs(Input6->get_Value())/VNorm);
	if((aux1>VMin&&(Input1->get_Value())>PickUp)||(aux2>VMin&&(Input3->get_Value())>PickUp)||(aux3>VMin&&(Input5->get_Value())>PickUp)){
            if(IntegralTime>=Temp){
                Cond=true;
            }
            IntegralTime+=dt;
		}
    else{
        if(IntegralTime>0){
            IntegralTime-=dt;
        }
    }
	Output->insert_Value(Cond);
	}
	
SubfrequenciaG30::SubfrequenciaG30(Channel<analogic>* In1, Channel<Complex>* In2, Channel<analogic>* In3, Channel<Complex>* In4, Channel<analogic>* In5, Channel<Complex>* In6, Channel<digital>* Out, double Vn, double MinV, double P, double D){
	FunctionalDescription_set("Overexcitation Protection");

	ANSICode_set("81");
	PickUp=P;
    Temp=D;
	Input1=In1;
	Input2=In2;
    Input3=In3;
	Input4=In4;
    Input5=In5;
	Input6=In6;
	Output=Out;
	VNorm=Vn;
	VMin=MinV;
	IntegralTime=0;
	}

bool SubfrequenciaG30::Prepare(float SamplingFreq){
	dt=1/SamplingFreq;
	return true;
	}

void SubfrequenciaG30::Run(){
	bool Cond=false;
	double aux1, aux2, aux3;
	aux1 = (abs(Input2->get_Value())/VNorm);
	aux2 = (abs(Input4->get_Value())/VNorm);
	aux3 = (abs(Input6->get_Value())/VNorm);
	if((aux1>VMin&&(Input1->get_Value())<PickUp)||(aux2>VMin&&(Input3->get_Value())<PickUp)||(aux3>VMin&&(Input5->get_Value())<PickUp)){
            if(IntegralTime>=Temp){
                Cond=true;
            }
            IntegralTime+=dt;
		}
    else{
        if(IntegralTime>0){
            IntegralTime-=dt;
        }
    }
	Output->insert_Value(Cond);
	}

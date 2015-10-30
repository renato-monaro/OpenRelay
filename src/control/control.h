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
#ifndef CONTROL_H
#define CONTROL_H

#include <math.h>
#include <vector>
#include "channel.h"
#include "../logic/logic.h"
#include "parameters.h"
#include "types.h"
#include <boost/circular_buffer.hpp>
namespace orelay{
class Control{
	public:
		void SetClock(Channel<timer> *);
		virtual void Run()=0;
		/*! Método que retorna uma string com a descrição funcional do controle.*/
		/** @return Retorna a descrição funcional.*/
		string FunctionalDescription_get(void);
		bool Prepare(float);
	protected:
		void FunctionalDescription_set(string);/*!< Método que atribui a descrição funcional ao Controle.*/   
		string FunctionalDescription;          /*!< Atributo que armazena a descrição funcional.*/
		Channel<timer> *MasterClock;
		float RelaySamplingRate;
	};
class Fault: public Control{
	public:
		Fault(Channel<analogic> *In,Channel<analogic> *Freq, Channel<digital>*Out, Channel<digital>*Ctrl,vector<double>Ang, vector<double>P, vector<double>Q, int Ncicles);
		void Run();
	protected:
		int FaultCicles;
		timer FaultTime;
		vector<double> Angles;
		vector<double> ActivePower;
		vector<double> ReactivePower;
		Channel<analogic> *Input;
		Channel<analogic> *Frequency;
		Channel<digital> *Output;
		Channel<digital> *Trigger;
		bool Trigged,Waiting,Armed,Resting;
		unsigned i,j,k;
};
/*
class PID: public Control{
	public:
		PID(Channel<analogic> *In, Channel<analogic>*Out, double Set, double Kc, double Ti, double Td);
		PID(Channel<analogic> *In, Channel<analogic>*Out, double Set, double Kc, double Ti, double Td,double LowLimit, double UpLimit);
		PID(Channel<analogic> *In, Channel<analogic>*Out, Channel<analogic> *Set, double Kc, double Ti, double Td);
		PID(Channel<analogic> *In, Channel<analogic>*Out, Channel<analogic> *Set, double Kc, double Ti, double Td,double LowLimit, double UpLimit);
		void Run();
	protected:
		Channel<analogic> *Input;
		Channel<analogic> *Output;
		Channel<analogic> *DynamicSetPoint;
		double SetPoint, q0,q1,q1w,q2,T0;
		bool Saturable;
		bool Dynamic;
		double LowerLimit;
		double UpperLimit;
		//boost::circular_buffer<double> Error;
	private:
		//double Action;
		double Error,Error_1,Error_2, Action;
};*/


class PI: public Control{
	public:
		PI(Channel<analogic> *In, Channel<analogic>*Out, double Set, double Kc, double Ti,double LowLimit, double UpLimit,double);
		PI(Channel<analogic> *In, Channel<analogic>*Out, Channel<analogic> *Set, double Kc, double Ti,double LowLimit, double UpLimit,double);
		PI(Channel<analogic> *In, Channel<analogic>*Out, double Set, double Kc, double Ti,double LowLimit, double UpLimit, bool Inv,double);
		PI(Channel<analogic> *In, Channel<analogic>*Out, Channel<analogic> *Set, double Kc, double Ti,double LowLimit, double UpLimit, bool Inv,double);
		void Run();
	protected:
		Channel<analogic> *Input;
		Channel<analogic> *Output;
		Channel<analogic> *DynamicSetPoint;
		double SetPoint, q0,q1,T0;
		bool Saturable;
		bool Dynamic;
		double LowerLimit;
		double UpperLimit;
		bool Invert;
		//boost::circular_buffer<double> Error;
	private:
		//double Action;
		double Error,Error_1, Action;
};

class TF: public Control{
	public:
		TF(Channel<analogic> *In,Channel<analogic> *Out,vector<double> Ak,vector<double> Bk);
		void Run();
	protected:
		Channel<analogic> *Input;
		Channel<analogic> *Output;
		vector<double> A,B;
	private:
		double Out;
};
class PWM: public Control{
	public:
		PWM(Channel<analogic> *In,Channel<digital> *Out1,Channel<digital> *Out2, double Amp, double Freq, double dT);
		void Run();
		bool Prepare(float);
	protected:
		Channel<analogic> *Input;
		Channel<digital> *Output_1;
		Channel<digital> *Output_2;
		double SwFreq;
	private:
		Channel<digital> *Output_1_int;
		Channel<digital> *Output_2_int;
		double Triang;
		bool Cres;
		double Amplitude;
		bool StatOut1,StatOut2;
		bool StatNewOut1,StatNewOut2;
		double CntOut1,CntOut2;
		Timer *T1,*T2;
};


class PWM2: public Control{
	public:
		PWM2(Channel<analogic> *Tri);
		PWM2(Channel<analogic> *Tri,Channel<digital> *Blk);
		bool Join(Channel<analogic> *Mod,Channel<digital> *Out1,Channel<digital> *Out2);
		bool Join(Channel<analogic> *Mod,Channel<digital> *Out1);
		void Run();
	protected:
	vector<Channel<analogic> *> vMod;
		Channel<analogic> *Triang;
		vector<Channel<digital> *> vOutput_1;
		vector<Channel<digital> *> vOutput_2;
		Channel<digital> *Blk;
	private:
		
};


class TriangularWave: public Control{
	public:
		TriangularWave(Channel<analogic> *Out, double Amp, double Freq, double dT);
		void Run();
	protected:
		Channel<analogic> *Output;
		double SwFreq;
	private:
		double Triang;
		double Amplitude;
		bool Cres;
};
}
#endif /*!CONTROL_H*/




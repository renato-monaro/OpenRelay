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
#ifndef LOGIC_H
#define LOGIC_H

#include <vector>
#include "channel.h"
#include "types.h"
namespace orelay{
class Logic{
	public:
		virtual void Run()=0;		
		virtual bool Prepare(float);
		/*! Método que retorna uma string com a descrição funcional da Lógica.*/
		/** @return Retorna a descrição funcional da Lógica.*/
		string FunctionalDescription_get(void);
		void SetClock(Channel<timer> *);
	protected:
		void FunctionalDescription_set(string);/*!< Método que atribui a descrição funcional à Lógica.*/   
		string FunctionalDescription;          /*!< Atributo que armazena a descrição funcional.*/
		Channel<digital> *Output;
		/*! Atributo que irá receber o ponteiro de tempo.*/
		Channel<timer> *MasterClock;
		float RelaySamplingRate;
};

class Single: public Logic{
	protected:
		Channel<digital> *Input;
};

class Multiple: public Logic{
	public:
		void Join_Channel(Channel<digital>*);
	protected:
		vector <Channel <digital>* > Input;
};

class Not: public Single{
	public:
		Not(Channel<digital> *In,Channel<digital> *Out);
		void Run();
};

class Or: public Multiple{
	public:
		Or(Channel<digital>*Out);
		void Run();
};

class And: public Multiple{
	public:
		And(Channel<digital>*Out);
		void Run();
};

class FlipFlop: public Logic{
	public:
		FlipFlop(Channel<digital> *InR,Channel<digital> *InS,Channel<digital> *Out);
		void Run();
	protected:
	Channel<digital> *InputR;
	Channel<digital> *InputS;
	Channel<digital> *Output;
};

class Timer: public Logic{
	public:
		Timer(Channel<digital> *In,Channel<digital> *Out,timer Time);
		//Timer(Channel<digital> *In,Channel<digital> *Out,analogic Time);
		Timer(Channel<digital> *In,Channel<digital> *Out,timer Time,timer Time2);
		//Timer(Channel<digital> *In,Channel<digital> *Out,analogic Time,analogic Time2);
		Timer(Channel<digital> *In,Channel<digital> *Out,Channel<timer> *Time,timer Time2);
		Timer(Channel<digital> *In,Channel<digital> *Out,Channel<analogic> *Time,analogic Time2);
		Timer(Channel<digital> *Out,timer Time);
		void Run();
		bool Prepare(float);
	protected:
		Channel<digital> *Input;
		Channel<digital> *Output;
		Channel<timer> *DynamicDelay;
		timer Delay;
		timer Reset;
		bool Dynamic;
		bool Resetable;
		bool Trigger;
	private:
		timer Accumulator,TimeStep;
};

class Selector: public Logic{
	public:
		Selector(Channel<digital> *In,Channel<analogic> *Out, vector<analogic> Values);
		Selector(Channel<digital> *In,Channel<string> *Out, vector<string> Values);
		Selector(Channel<digital> *In,Channel<digital> *Reset,Channel<analogic> *Out, vector<analogic> Values);
		Selector(Channel<digital> *In,Channel<timer> *Out, vector<timer> Values);
		Selector(Channel<digital> *In,Channel<digital> *Reset,Channel<timer> *Out, vector<timer> Values);
		void Run();
	protected:
		Channel<digital> *Input;
		Channel<digital> *Reset;
		Channel<analogic> *AOutput;
		Channel<string> *SOutput;
		Channel<timer> *TOutput;
		vector<analogic> AValues;
		vector<timer> TValues;
		vector<string> SValues;
		bool resetable;
		unsigned Type;
	private:
		unsigned ActualPos;
};

class Counter: public Logic{
	public:
		Counter(Channel<digital> *In,Channel<digital> *Out,unsigned Count);
		Counter(Channel<digital> *In,Channel<digital> *Rs,Channel<digital> *Out,unsigned Count);
		void Run();
	protected:
		Channel<digital> *Input;
		Channel<digital> *Output;
		Channel<digital> *Reset;
		unsigned Count;
		bool Resetable;
		unsigned ActualPos;
};

class RiseDown: public Single{
	public:
		RiseDown(Channel<digital> *In,Channel<digital> *Out);
		void Run();
};

class RiseUp: public Single{
	public:
		RiseUp(Channel<digital> *In,Channel<digital> *Out);
		void Run();
};

class LockoutCounter: public Logic{
	public:
		LockoutCounter(Channel<digital> *In,Channel<digital> *Out,unsigned Count);
		LockoutCounter(Channel<digital> *In,Channel<digital> *Rs,Channel<digital> *Out,unsigned Count);
		void Run();
	protected:
		Channel<digital> *Input;
		Channel<digital> *Output;
		Channel<digital> *Reset;
		unsigned Count;
		bool Resetable;
		unsigned ActualPos;
};
}
#endif /*!LOGIC_H*/

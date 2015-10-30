/*
# Copyright (C) 2008-2015 Renato Machado Monaro, Hermes Manoel
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
#ifndef CHANNEL_H
#define CHANNEL_H

#include <string>
#include <iostream>
#include <boost/circular_buffer.hpp>
//#include <boost/thread/mutex.hpp>
#include <typeinfo>
#include "parameters.h"
#include "types.h"

using namespace std;
namespace orelay{
//! Gabarito da classe Channel (template). 
template <class T>
//! Classe Channel.
class Channel
	{
	public:

		/*! Construtor da classe Channel.*/
		/** @param name refere-se ao nome do canal interno.
		   @param name refere-se ao nome do canal interno. */
		Channel(string name, unsigned BufferSize);
		Channel(string name);
		/*! Destrutor da classe Channel.*/
		~Channel(){};
		/*! Método para inserir valor no buffer circular.*/
		void insert_Value(T);
		/*! Retorna o valor atual (último valor) que está no buffer circular.*/
		T get_Value();
		/*! Retorna o valor da posição do buffer escolhida.*/
		T get_Value(int);
		/*! Adiciona dois valores contidos em duas posições do buffer.*/
		T Add(int,int);
		/*! Subtrai dois valores contidos em duas posições do buffer.*/
		T Sub(int,int);
		/*! Calcula a média de todos os valores contidos no buffer.*/
		T Average(); 
		/*! Obtém um vetor de tamanho definido do buffer circular começando pela posição mais atual.*/
		/** @param N é o tamanho do vetor desejado.*/
		vector<T> get_Vector(int N); 
		/*! Verifica se o buffer circular está completo.*/
		bool Full();
		/*! Método que retorna uma string com o nome do canal.*/
		string get_Name(void){return Name;};
		/*! Método para impressão do conteúdo do buffer circular.*/
		void Print();
		/*! Método que retorna o tamanho do buffer circular.*/
		unsigned get_Size(){return Size;};
		/*! Método que retorna o tamanho do buffer circular.*/
		void set_Size();
		/*! Método que retorna o tipo do canal.*/
		int get_Type(){return Type;};
		/*! Método para definir o tamanho do buffer.*/
		/** @param BufferSize é taxa o tamannho do Buffer do canal*/
		void set_Size(unsigned BufferSize);
		/*! Método para definir a taxa de amostragem.*/
		/** @param Sampling é taxa de amostragem utilizada no canal*/
		//void set_Sampling(float Sampling);
		/*! Método que retorna a taxa de amostragem.*/		
		//float get_Sampling(){return SamplingRate;};
		/*! Método que bloqueia reserva a escrita no canal a uma determinada classe. Retorna true em caso de sucesso.*/		
		bool UseIt();
		/*! Método que indica se o canal já esta em uso.*/		
		bool InUse(){return inUse;};

//Remoção Futura 

/*
		Channel(string name,bool direction,unsigned char channel,float ratio);
		Channel(string name,bool direction,unsigned char channel,float ratio, float offset);
		Channel(string name,bool direction,unsigned char channel,float ratio,double Low, double Up);
		Channel(string name,bool direction,unsigned char channel,float ratio, float offset,double Low, double Up);
		Channel(string name,bool direction,unsigned char channel);
		Channel(string name);



		float get_Ratio(){return Ratio;};

		float get_Offset(){return Offset;};

		unsigned char get_BoardChannel(){return BoardChannel;};

		int get_Direction(){return Direction;}

		bool is_Saturable(){return Saturable;};

		double get_UpperLimit(){return UpperLimit;};

		double get_LowerLimit(){return LowerLimit;};
*/
	protected:
		/*! Buffer circular.*/
		boost::circular_buffer <T>Value;
		/*! Atributo referente ao tipo do canal.*/
		int Type;
		/*! Atributo referente ao tamanho do buffer. */
		unsigned Size;
		/*! Atributo referente ao nome do canal.*/
		string Name;
		//float SamplingRate;
		bool inUse;
//Remoção Futura 
/*

		void set_Name(string);
		void set_Description(string);
		string Description;
		float Ratio;
		float Offset;
		int Direction;
		unsigned char BoardChannel;		
		double UpperLimit; 
		double LowerLimit; 
		bool Saturable;

*/
	};



template <class T>
Channel<T>::Channel(string name,unsigned size)
{
	Name=name;
	if(typeid(T)==typeid(analogic))
			Type=ANALOG;
	if(typeid(T)==typeid(digital))
			Type=DIGITAL;
	if(typeid(T)==typeid(timer))
			Type=TIME;
	if(typeid(T)==typeid(string))
			Type=STRING;
	if(typeid(T)==typeid(Complex))
			Type=COMPLEX;
	Size=size;
	Value.resize(size);
	//SamplingRate=Sampling;
	inUse=false;
	}
template <class T>
Channel<T>::Channel(string name)
{
	Name=name;
	if(typeid(T)==typeid(analogic))
			Type=ANALOG;
	if(typeid(T)==typeid(digital))
			Type=DIGITAL;
	if(typeid(T)==typeid(timer))
			Type=TIME;
	if(typeid(T)==typeid(string))
			Type=STRING;
	if(typeid(T)==typeid(Complex))
			Type=COMPLEX;
	//SamplingRate=Sampling;
	inUse=false;
	}
template <class T>
bool  Channel<T>::UseIt(){
	if(inUse)
		return false;
	else 
		return true;
	}
template <class T>
void Channel<T>::Print(){
	for (unsigned int i=0;i<Value.size();i++)
	cout<<Value[i]<<endl;
	}

template <class T>
void Channel<T>::insert_Value(T In){
	Value.push_back(In);
	}

template <class T>
T Channel<T>::get_Value(){
	 return Value.back();
	}
template <class T>
T Channel<T>::get_Value(int pos){
	return Value[Size-pos-1];
	}
template <class T>
T Channel<T>::Average(){
	T Average;
	Average=Value[0];
	for(unsigned i=1;i<Size;i++){
		Average+=Value[i];
		}
	return Average/Size;
	}
template <class T>
T Channel<T>::Sub(int start, int end){
	return (Value[start]-Value[end]);
	}

template <class T>
T Channel<T>::Add(int start, int end){
	return (Value[start]+Value[end]);
	}
template <class T>
bool Channel<T>::Full(){
	return Value.full();
	}
template <class T>
vector<T> Channel<T>::get_Vector(int N){
	vector<T> tmp;
	for(int k=0;k<N;k++)
		tmp.push_back(this->get_Value(k));
	tmp.resize(tmp.size());
	return tmp;
	} 

/*template <class T>
void Channel<T>::set_Sampling(float Sampling){
	SamplingRate=Sampling;
	}*/
template <class T>
void Channel<T>::set_Size(unsigned BufferSize){
	Size=BufferSize;
	Value.resize(Size);
	}



//Remoção Futura 
/*
template <class T>
Channel<T>::Channel(string name,bool direction,unsigned char channel,float ratio)
{
	Name=name;
	Direction=direction;
	BoardChannel=channel;
	Ratio=ratio;
	Offset=0.0;
	if(typeid(T)==typeid(analogic))
			Type=ANALOG;
	if(typeid(T)==typeid(digital))
			Type=DIGITAL;
	if(typeid(T)==typeid(timer))
			Type=TIME;
	Saturable=false;
//	Size=BUFFER_SIZE;
//	Value.set_capacity(BUFFER_SIZE);
//	SamplingRate=1E9/(ACQUISITION_TICK);
//	if(Direction==OUTPUT)
//		Value.resize(BUFFER_SIZE);
//	else
//		Value.resize(1);
}

template <class T>
Channel<T>::Channel(string name,bool direction,unsigned char channel,float ratio, float offset)
{
	Name=name;
	Direction=direction;
	BoardChannel=channel;
	Ratio=ratio;
	Offset=offset;
	if(typeid(T)==typeid(analogic))
			Type=ANALOG;
	if(typeid(T)==typeid(digital))
			Type=DIGITAL;
	if(typeid(T)==typeid(timer))
			Type=TIME;
	Saturable=false;
//	Size=BUFFER_SIZE;
//	Value.set_capacity(BUFFER_SIZE);
//	SamplingRate=1E9/(ACQUISITION_TICK);
//	if(Direction==OUTPUT)
//		Value.resize(BUFFER_SIZE);
//	else
//		Value.resize(1);
}

template <class T>
Channel<T>::Channel(string name,bool direction,unsigned char channel,float ratio,double Low, double Up)
{
	Name=name;
	Direction=direction;
	BoardChannel=channel;
	Ratio=ratio;
	Offset=0.0;
	if(typeid(T)==typeid(analogic))
			Type=ANALOG;
	if(typeid(T)==typeid(digital))
			Type=DIGITAL;
	if(typeid(T)==typeid(timer))
			Type=TIME;
	LowerLimit=Low;
	UpperLimit=Up;
	Saturable=true;
//	Size=BUFFER_SIZE;
//	Value.set_capacity(BUFFER_SIZE);
//	SamplingRate=1E9/(ACQUISITION_TICK);
//	if(Direction==OUTPUT)
//		Value.resize(BUFFER_SIZE);
//	else
//		Value.resize(1);
}

template <class T>
Channel<T>::Channel(string name,bool direction,unsigned char channel,float ratio, float offset,double Low, double Up)
{
	Name=name;
	Direction=direction;
	BoardChannel=channel;
	Ratio=ratio;
	Offset=offset;
	if(typeid(T)==typeid(analogic))
			Type=ANALOG;
	if(typeid(T)==typeid(digital))
			Type=DIGITAL;
	if(typeid(T)==typeid(timer))
			Type=TIME;
	LowerLimit=Low;
	UpperLimit=Up;
	Saturable=true;
//	Size=BUFFER_SIZE;
//	Value.set_capacity(BUFFER_SIZE);
//	SamplingRate=1E9/(ACQUISITION_TICK);
//	if(Direction==OUTPUT)
//		Value.resize(BUFFER_SIZE);
//	else
//		Value.resize(1);
}


template <class T>
Channel<T>::Channel(string name,bool direction,unsigned char channel)
{
	Name=name;
	Direction=direction;
	BoardChannel=channel;
	Ratio=1;
	Offset=0.0;
	if(typeid(T)==typeid(analogic))
			Type=ANALOG;
	if(typeid(T)==typeid(digital))
			Type=DIGITAL;
	if(typeid(T)==typeid(timer))
			Type=TIME;
	Saturable=false;
	//Value.resize(Size);
}

template <class T>
Channel<T>::Channel(string name)
{
	Name=name;
	Direction=INTERNAL;
	BoardChannel=-1;
	Ratio=0;
	Offset=0.0;
	if(typeid(T)==typeid(analogic))
			Type=ANALOG;
	if(typeid(T)==typeid(digital))
			Type=DIGITAL;
	if(typeid(T)==typeid(timer))
			Type=TIME;
	Saturable=false;
}

*/
}
#endif /* !CHANNEL_H */

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
#ifndef OSCILLOGRAPHY_H
#define OSCILLOGRAPHY_H

#include <fstream>
#include <vector>
#include <string>
#include <list>
#include <deque>
#include "channel.h"
#include "parameters.h"
#include "types.h"
namespace orelay{
class Oscillogram{
	public:
		//Oscillogram(vector<vector<timer> > time, vector<vector<analogic> > analog, vector<vector<digital> > dig,vector<vector<Complex> > cp, vector<string> n, string F);
		Oscillogram(vector<string> n, string path, string extension, unsigned sizeTimer, unsigned sizeAnalogic, unsigned sizeDigital, unsigned sizeComplex, unsigned numPoints);
		void insert_time(unsigned int, unsigned int, timer);
		void insert_analogic(unsigned int, unsigned int, analogic);
		void insert_digital(unsigned int, unsigned int, digital);
		void insert_complex(unsigned int, unsigned int, Complex);
		void set_ready(bool);
		bool get_ready();
		bool Save();
		void set_file_name(string file_name);
		void set_file_name(int num_file);
		string get_file_name();
	protected:
		vector< vector<timer> > t;
		vector< vector<analogic> > a;
		vector< vector<digital> > d;
		vector< vector<Complex> > ph;
		vector<string> names;
		int num_arq;
		string FileName;
		string extensionArq;
		string caminho;
		bool ready;
};

//!Classe Oscillography.
class SaveData{
	public:	
		bool Prepare(float,unsigned,float);
		/*! Método que adiciona um canal digital à oscilografia.*/
		void Join_Channel(Channel<digital>*);
		/*! Método que adiciona um canal complexo à oscilografia.*/
		void Join_Channel(Channel<Complex>*);
		/*! Método que adiciona um canal analógico à oscilografia.*/
		void Join_Channel(Channel<analogic>*);
		/*! Método que adiciona um canal de tempo à oscilografia.*/
		void Join_Channel(Channel<timer>*);
		/*! Método para monitorar a condição para iniciar a oscilografia.*/
		void Run();
		/*! Definição do ponteiro de tempo para a oscilografia.*/
		void SetClock(Channel<timer> *);
		/*! Definição do ponteiro de tempo para a oscilografia.*/
		void SetControl(bool *);
		/*! A partir da condição satisfeita inicia-se a oscilografia.*/
		//virtual void CreateOscillogram()=0;
		/*! Método para salvar a oscilografia.*/
		/** Quando os dados para a oscilografia estão disponíveis este método é utilizado para salvar em disco no diretório /tmp/ a oscilografia.*/
		/** @param t é o argumento as informações de tempo.
		@param a é o argumento as informações dos canais analógicos.
		@param d é o argumento as informações dos canais digitais.
		@param ph é o argumento as informações dos canais de phasor.
		@param p é o argumento as informações dos canais de potência.
		@param names é o argumento as informações dos nomes dos canais.
		@param N é o nome do arquivo.*/
		void Poll();
		/*! Método que retorna uma string com a descrição funcional da oscilografia.*/
		/** @return Retorna a descrição funcional.*/
		string FunctionalDescription_get(void);
		//unsigned RequiredBufferSize(void);
		unsigned OscillogramSize(){return (PreFaultPoints);};
	protected:
		void FunctionalDescription_set(string);/*!< Método que atribui a descrição funcional à Oscilografia.*/   
		string FunctionalDescription;          /*!< Atributo que armazena a descrição funcional.*/
		/*!Nome da Oscillografia.*/
		string Name;
		/*! Pontos de pré falta.*/
		unsigned int PreFaultPoints;
		/*! Pontos de pós falta.*/
		unsigned int PostFaultPoints;
		/*! Valor que irá inciar a oscilografia.*/
		int Trigged;
		/*! Salva o instante da falta.*/
//		timer FaultInstant;
		/*! Canal para armazenar o valor do disparo da oscilografia.*/
		Channel<digital> *Trigger;
		/*! Vetor que armazena os ponteiros dos canais analógicos que pertencerão à oscilografia. */
		vector <Channel<analogic>* > Analogic;
		/*! Vetor que armazena os ponteiros dos canais digitais que pertencerão à oscilografia*/
		vector <Channel<digital>* > Digital;
		/*! Vetor que armazena os ponteiros dos canais complexos que pertencerão à oscilografia*/
		vector <Channel<Complex>* > Phasors;
		/*! Vetor que armazena o ponteiro do canal de tempo.*/
		vector <Channel<timer>* > Timer;
		/*! Atributo que irá receber o ponteiro de tempo.*/
		Channel<timer> *MasterClock;
		/*! Número de oscilografias solicitadas*/
		long 	NOsc_criadas, 
				NOsc_salvas;
		//boost::circular_buffer<Oscillogram> OscillogramStack;
		vector<Oscillogram> OscillogramStack;
		unsigned max_oscillogram;
		bool *Running;
		float RelaySamplingRate;
		bool Dynamic;
		Channel<string> *FileName;
		vector<string> names;
		string File;
		int BeforeCicles,AfterCicles;
		string extension;
		//unsigned int Pos;
/*------------
		vector< vector<timer> > t;
		vector< vector<analogic> > a;
		vector< vector<digital> > d;
		vector< vector<Complex> > ph;
------------*/
};
class Oscillography : public SaveData{
	public:
		/*! Construtor da classe Oscillography.*/
		/** Quando o relé manda um sinal de trip, faz com que a oscilografia seja iniciada.*/
		/** @param Trig é o parâmetro que define o inicio da oscilografia.
		@param BeforeCicles constituem os ciclos antes da falta.
		@param AfterCicles constituem os ciclos depois da falta.
		@param Frequency é a Frequência nominal de operação do Sistema.
		@param N é o nome da Pasta onde as Oscilografias serão salvas.*/
		Oscillography(Channel<digital> *Trig,int BeforeCicles,int AfterCicles,string N);
		Oscillography(Channel<digital> *Trig,int BeforeCicles,int AfterCicles,string N,unsigned num_max_oscillogram);
		Oscillography(Channel<digital> *Trig,int BeforeCicles,int AfterCicles,string N,Channel<string> *Name);
		Oscillography(Channel<digital> *Trig,int BeforeCicles,int AfterCicles,string N,Channel<string> *Name,unsigned num_max_oscillogram);
		/*! Construtor da classe Oscillography.*/
		/** Quando o relé manda um sinal de trip, faz com que a oscilografia seja iniciada.*/
		/** @param Trig é o parâmetro que define o inicio da oscilografia.
		@param BeforePoints constituem os pontos antes da falta.
		@param AfterPoints constituem os pontos depois da falta.
		@param N é o nome da Pasta onde as Oscilografias serão salvas.*/
	//	Oscillography(Channel<digital> *Trig,int BeforePoints,int AfterPoints,string N);
		//void CreateOscillogram();
		//bool Prepare(float,unsigned,float);
		//void Run();
};
class Evaluation : public SaveData{
	public:
		Evaluation(Channel<digital> *Trig,int BeforeCicles,int AfterCicles,string N);
		Evaluation(Channel<digital> *Trig,int BeforeCicles,int AfterCicles,string N,unsigned num_max_oscillogram);
		//svoid CreateOscillogram();	
	};
}
#endif /*!OSCILLOGRAPHY_H*/


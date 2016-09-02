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
#ifndef RELAY_H
#define RELAY_H

#include <parameters.h>
#include <types.h>
#include <channel.h>
#include <acquisition/acquisition.h>
#include <protection/protection.h>
#include <measures/measures.h>
#include <measures/wavelet.h>
#include <measures/switch.h>
#include <measures/ddf.h>
#include <measures/phasor.h>
#include <control/control.h>
#include <savedata/oscillography.h>
#ifdef NCURSES
	#include <display/nscreen.h>
#endif
#include <filter.h>
#include <logic/logic.h>
#include <acquisition/encoder.h>
#ifdef IEC61850
#include <acquisition/IEC61850.h>
#include <Publisher.h>
#include <Subscriber.h>
#endif
#include <measures/WangSunFrequency.h>
#include <measures/fuzzy.h>
#include <measures/rna.h>
//#include <interface/openrelayxml.h>
//#include <measures/rna.h> 
#include <string>
#include <vector>
#include <boost/thread/thread.hpp>

#ifdef RTAI
	#include <rtai_lxrt.h>
#endif

namespace orelay{
//! Classe Relay.
/** Esta é uma classe que irá agregar e controlar as funções do relé.*/
class Relay{
	public:
		/*!Construtor*/
		/** @param buffer_size é o número de amostras que os canais reterão.*/
		/** @param sampling_rate é a taxa de amostragem do relé em Hz.*/
		Relay(unsigned buffer_size,float sampling_rate, float frequency);
		/*!Destrutor*/
		~Relay();
		/*!Cria um canal analógico no relé.*/
		Channel<analogic>* CreateAnalogicChannel(string name);
		/*!Cria um canal digital no relé.*/
		Channel<digital>* CreateDigitalChannel(string name);
		/*!Cria um canal para String no relé.*/
		Channel<string>* CreateStringChannel(string name);
		/*!Cria um canal de tempo no relé.*/
		Channel<timer>* CreateTimerChannel(string name);
		/*!Cria um canal complexo no relé.*/
		Channel<Complex>* CreateComplexChannel(string name);
		/*!Busca o canal pela sua descrição*/
		/** @param name Designação dado ao Canal no momenento de sua criação.*/
		/** @return O ponteiro do Canal caso este exita, NULL caso contrário.*/
		Channel<analogic>* GetAnalogicChannel(string name);
		/*!Busca o canal pela sua descrição*/
		/** @param name Designação dado ao Canal no momenento de sua criação.*/
		/** @return O ponteiro do Canal caso este exita, NULL caso contrário.*/
		Channel<digital>* GetDigitalChannel(string name);
		/*!Busca o canal pela sua descrição*/
		/** @param name Designação dado ao Canal no momenento de sua criação.*/
		/** @return O ponteiro do Canal caso este exita, NULL caso contrário.*/
		Channel<string>* GetStringChannel(string name);
		/*!Busca o canal pela sua descrição*/
		/** @param name Designação dado ao Canal no momenento de sua criação.*/
		/** @return O ponteiro do Canal caso este exita, NULL caso contrário.*/
		Channel<timer>* GetTimerChannel(string name);
		/*!Busca o canal pela sua descrição*/
		/** @param name Designação dado ao Canal no momenento de sua criação.*/
		/** @return O ponteiro do Canal caso este exita, NULL caso contrário.*/
		Channel<Complex>* GetComplexChannel(string name);
		/*!Atribui o ponteiro para objeto do tipo aquisição ao relé.*/
		void Join(Acquisition*);
		bool Reset();

		/*!Adiciona o ponteiro para objeto do tipo proteção à lista de execução das proteções.*/
		void Join(Protection*);
		/*!Adiciona o ponteiro para objeto do tipo medida à lista de execução das medidas.*/ 
		void Join(Measure*);
		/*!Adiciona o ponteiro para objeto do tipo controle à lista de execução do controle.*/
		void Join(Control*);
		/*!Adiciona o ponteiro para objeto do tipo oscilografia à lista de execução da oscilografia.*/
		void Join(Oscillography*);
		/*!Adiciona o ponteiro para objeto do tipo oscilografia à lista de execução da avaliação de desempenho.*/
		void Join(Evaluation*);
#ifdef NCURSES
		/*!Adiciona o ponteiro para objeto do tipo display à lista de execução do display.*/
		void Join(Display*);
#endif
		/*!Adiciona o ponteiro para objeto do tipo Logic à lista de execução das operações lógicas.*/
		void Join(Logic*);
		/*! Método para listar as funções que serão desenvolvidas no relé.*/
		void List();
		/*! Método para executar o relé.*/
		void Execute();
		bool Step();
		bool Prepare();
		bool Wait();
		void PrepareEvaluation(unsigned );
		void Stop();
		void Start();

		void JoinChannel(Channel<analogic>*);
		void JoinChannel(Channel<digital>*);
		void JoinChannel(Channel<string>*);
		void JoinChannel(Channel<timer>*);
		void JoinChannel(Channel<Complex>*);
		bool Running;
	protected:
		/*! Armazena o endereço da aquisição.*/
		vector<Acquisition*> AcqList;
		/*!Vetor que armazena o ponteiro para os objetos de proteção.*/
		vector<Protection*> ProtectionList; 
		/*!Vetor que armazena o ponteiro para os objetos da medição.*/
		vector<Measure*>  MeasureList;
		/*!Vetor que armazena o ponteiro para os objetos do controle.*/
		vector<Control*>  ControlList;
		/*!Vetor que armazena o ponteiro para os objetos da oscilografia.*/
		vector<Oscillography*> OscillographyList;
		/*!Vetor que armazena o ponteiro para os objetos da oscilografia.*/
		vector<Evaluation*> EvaluationList;
		/*!Vetor que armazena o ponteiro para os objetos de operações lógicas.*/
		vector<Logic*> LogicList;
		/*!Vetor que armazena o ponteiro para os objetos de avaliação de desempenho.*/
		vector<Oscillography*> PerformanceList;
		/*!Vetor que armazena o ponteiro para os objetos do display.*/
		#ifdef NCURSES
		vector<Display*> DisplayList;
		#endif
		/*!< Vetor que armazena os ponteiros dos canais analógicos.*/
		vector <Channel<analogic>* > AnalogicChannels;
		/*!< Vetor que armazena os ponteiros dos canais complexos.*/
		vector <Channel<Complex>* > ComplexChannels;
		/*!< Vetor que armazena os ponteiros dos canais digitais.*/
		vector <Channel<digital>* > DigitalChannels;
		/*!< Vetor que armazena os ponteiros dos canais de tempo.*/
		vector <Channel<timer>* > TimerChannels; 
		/*!< Vetor que armazena os ponteiros dos canais de String.*/
		vector <Channel<string>* > StringChannels; 
		/*!Atributo que armazena o ponteiro de tempo usado no relé.*/
		Channel<timer> *MasterClock;
		/*! Atributo para controlar a execução do relé.*/

		/*! Atributo que armazena o ponteiro para o canal de disparo.*/
		Channel<digital> *Trigger;
		/*! Atributo que armazena o tempo de execução de cada função do relé.*/
		vector <Channel<analogic>* > ExecutionTimeList;
		/*! Atributo que armazena o se ocorreu overrun.*/
		Channel<digital> *ExecutionOverrun;
//		unsigned AfterPoints,BeforPoints;
		bool PerformanceEvaluation;
		boost::thread_group EvalThreadList;
		boost::thread_group OscThreadList;
		/*! Número de oscilografias solicitadas*/
//		long NOsc;
		/*! Valor que irá inciar a oscilografia.*/
//		int Trigged;
		/*!Nome da Oscillografia.*/
//		string Name;	
		unsigned BufferSize;
		float SamplingRate;
		float SystemFrequency;
	private:
		unsigned sizeAcq,sizeMeas,sizeProt,sizeCont,sizeLog,sizeOsc,sizeEval;
		timer T1,T2,Tb;
		boost::thread_group ThreadList;
	#ifdef RTAI
		RT_TASK *tarefa_principal;
		long long unsigned Ov;
	#endif
};
}
#endif /*!RELAY_H*/


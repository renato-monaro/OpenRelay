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
#ifndef ACQUISITION_H
#define ACQUISITION_H

#include <string>
#include <fstream>
#include <boost/circular_buffer.hpp>
#include <boost/thread/thread.hpp>
//#include <boost/thread/mutex.hpp>
//#include <boost/thread/barrier.hpp>
//#include <boost/date_time/posix_time/posix_time.hpp>
#include <pthread.h>
#include <vector>
#ifdef RTAI
	#ifdef COMEDI
		#include <rtai_comedi.h>
	#endif
#endif
#ifndef RTAI
	#ifdef COMEDI
		#include <comedilib.h>
	#endif
#endif
#include "channel.h"
#include "types.h"
#include "parameters.h"
#include "filter.h"
//#include <acquisition/IEC61850.h>

namespace orelay{
//! Classe Acquisition.
/** Esta é uma classe base para as outras classes do algoritmo de aquisição.*/
class  Acquisition{
	public:
		virtual void RefreshTime()=0;
		virtual bool Prepare(float)=0; /*!< Método virtual pura para verificação da coerência dos parâmetros da lista de canais.*/
		virtual bool Run()=0; /*!< Método virtual pura necessária para a execução da aquisição.*/   

		string FunctionalDescription_get(void);

		void SetClock(Channel<timer> *); /*!< Definição do ponteiro de tempo para a aquisição.*/

		void Join_Channel(Channel<digital>*,bool direction,unsigned char channel); /*!< Método que adiciona um canal digital à lista de execução da aquisição.*/

		void Join_Channel(Channel<analogic>*,bool direction,unsigned char channel); /*!< Método que adiciona um canal digital à lista de execução da aquisição.*/
		void Join_Channel(Channel<analogic>*,bool direction,unsigned char channel,float ratio); /*!< Método que adiciona um canal digital à lista de execução da aquisição.*/
		void Join_Channel(Channel<analogic>*,bool direction,unsigned char channel,float ratio, float offset); /*!< Método que adiciona um canal digital à lista de execução da aquisição.*/

//		void SetBufferSize(unsigned);
//		void SetSampling(float);

//		void Join_Channel(Channel<Complex>*); /*!< Método que adiciona um canal complexo à lista de execução da aquisição.*/ 
//		void Join_Channel(Channel<analogic>*); /*!< Método que adiciona um canal analógico à lista de execução da aquisição.*/
//		void Join_Channel(Channel<timer>*); /*!< Método que adiciona um canal de tempo à lista de execução da aquisição.*/
//		void Join_Channel(Channel<string>*); /*!< Método que adiciona um canal de string à lista de execução da aquisição.*/
		
		//float get_Sampling(){return SamplingRate;};  /*!< Método que retorna a taxa de amostragem da aquisição.*/
		//unsigned get_BufferSize(){return BufferSize;}; /*!< Método que retorna o tamnaho do buffer dos canais da aquisição.*/
		/*! Método que retorna uma string com a descrição funcional da aquisição.*/
		/** @return Retorna a descrição funcional da aquisição.*/
	protected:
		void FunctionalDescription_set(string);/*!< Método que atribui a descrição funcional à Aquisição.*/   
		string FunctionalDescription;          /*!< Atributo que armazena a descrição funcional.*/

		vector <Channel<analogic>* > AnalogicOut; /*!< Vetor que armazena os ponteiros dos canais analógicos de saída da aquisição.*/
		vector <Channel<analogic>* > AnalogicIn;  /*!< Vetor que armazena os ponteiros dos canais analógicos de entrada da aquisição.*/
		vector <Channel<digital>* > DigitalOut; /*!< Vetor que armazena os ponteiros dos canais digitais de saída da aquisição.*/
		vector <Channel<digital>* > DigitalIn; /*!< Vetor que armazena os ponteiros dos canais digitais de entrada da aquisição.*/

		vector <unsigned> AnalogicOutBoardChannel; /*!< Vetor que armazena os ponteiros dos canais analógicos de saída da aquisição.*/
		vector <unsigned> AnalogicInBoardChannel;  /*!< Vetor que armazena os ponteiros dos canais analógicos de entrada da aquisição.*/
		vector <unsigned> DigitalOutBoardChannel; /*!< Vetor que armazena os ponteiros dos canais digitais de saída da aquisição.*/
		vector <unsigned> DigitalInBoardChannel; /*!< Vetor que armazena os ponteiros dos canais digitais de entrada da aquisição.*/

		vector <float> AnalogicOutRatio; /*!< Vetor que armazena os ponteiros dos canais analógicos de saída da aquisição.*/
		vector <float> AnalogicInRatio;  /*!< Vetor que armazena os ponteiros dos canais analógicos de entrada da aquisição.*/

		vector <float> AnalogicOutOffset; /*!< Vetor que armazena os ponteiros dos canais analógicos de saída da aquisição.*/
		vector <float> AnalogicInOffset;  /*!< Vetor que armazena os ponteiros dos canais analógicos de entrada da aquisição.*/

		Channel<timer> *MasterClock; /*!< Atributo que irá receber o ponteiro de tempo, o qual servirá de base de tempo para as demais funções do relé.*/
//		vector <Channel<analogic>* > AnalogicInternal; /*!< Vetor que armazena os ponteiros dos canais analógicos internos da aquisição.*/
//		vector <Channel<Complex>* > Phasors; /*!< Vetor que armazena os ponteiros dos canais complexos internos da aquisição.*/
//		vector <Channel<digital>* > DigitalInternal; /*!<Vetor que armazena os ponteiros dos canais digitais internos da aquisição..*/
//		vector <Channel<timer>* > Timer; /*!< Vetor que armazena o ponteiro de tempo da aquisição..*/
//		vector <Channel<string>* > Strings; /*!< Vetor que armazena o ponteiro de tempo da aquisição..*/
//		float SamplingRate; /*! Atributo que armazena a taxa de amostragem empregada na aquisição.*/
//		unsigned BufferSize; /*! Atributo que armazena o tamanho do buffer necessário para os canais*/
	};
#ifdef COMEDI
//! Classe HardwareAcquisition.
/** Esta é uma classe derivada da classe Acquisition, responsável pela aquisição no hardware.*/
class HardwareAcquisition: public Acquisition{
	public:
		/*!Construtor da classe HardwareAcquisition.*/
		/** @param device é o caminho para a placa de aquisição de dados que fornecerá dados para o relé.*/
		HardwareAcquisition(const char *device);
		bool Run();  /*!< Método virtual pura necessária para a execução da aquisição em hardware.*/      
		/** Este método (virtual pura) é necessário para verificar se o tipo de canal solicitado (entrada/saída/interno,analógico/digital) realmente existe na lista de canais, ou seja, verifica se existe coerência entre o canal solicitado e os canais da lista de canais.*/
		bool Prepare(float);
		void RefreshTime();
	private:	
		comedi_t *comediDev;  
		lsampl_t *data_analog_out, *data_analog_in, *data_digital_in, *data_digital_out; 
		timer *t; 
		double sampling_rate; 
		unsigned int range,aref; 
//		comedi_polynomial_t *AnalogicOutCalibration;  
//		comedi_polynomial_t  *AnalogicInCalibration;  
		int subdevice_analog_out, subdevice_analog_in, subdevice_digital_out, subdevice_digital_in;
		lsampl_t *analog_in_maxdata;
		lsampl_t *analog_out_maxdata;
		#ifdef RTAI
			comedi_krange *analog_in_range;
			comedi_krange *analog_out_range;
		#endif
		#ifndef RTAI
			comedi_range **analog_in_range;
			comedi_range **analog_out_range;
		#endif
	};
#endif
//! Classe FileAcquisition.
/** Esta é uma classe derivada da classe Acquisition, responsável pela aquisição em arquivo.
@todo Observar a taxa de amostragem do arquivo de entrada e ver se esta é adequada ao requisitado pelo relé.*/
class FileAcquisition: public Acquisition{
	public:
		/*!Construtor da classe FileAcquisition.*/
		/** @param filename é o caminho para o arquivo que fornecerá dados para o relé.
			@param Bits é o número de bits do Conversor Analógico Digital. Bits=0 implica remoção do CAD.*/
//		FileAcquisition(const char *filename, unsigned Bits);
		/*!Construtor da classe FileAcquisition.*/
		/** @param filename é o caminho para o arquivo que fornecerá dados para o relé.*/
		FileAcquisition(const char *filename);
		FileAcquisition(const char *filename, float Samp);
		virtual bool Run(); /*!< Método virtual pura necessária para a execução da aquisição em arquivo.*/
		/** Este método (virtual pura) é necessário para verificar se o arquivo solicitado (a coluna do arquivo) realmente existe, caso isso não ocorra, haverá um erro. Além disso, a primeira coluna sempre será ocupada pelo MasterClock.*/		
		virtual bool Prepare(float); 
		/*!Método que obtem o valor numérico da coluna especificada.*/
		/** @param Line é uma string que contém uma linha do arquivo.
		@param Index é o número da coluna.*/
		double getValueColumnN(string Line, unsigned int Index); 
		void RefreshTime();
	protected:
		vector<vector<double> > Data; /*!<Atributo que Armazena os dados de entrada previamente reamostrado e filtrado(butterworth segunda ordem corte em SamplingRate/2)*/
		//unsigned DACBits; /*!<Número de Bits do Conversor Analógico-Digital*/
	private:
		ifstream File; 
		unsigned int Columns;
		timer EventDateTime; 
		float FileSamplingRate;
		float SamplingRate;
		unsigned DataCounter;
	};

#ifdef RTAI
#ifdef COMEDI
	double comedi_to_phys(lsampl_t data,comedi_krange *rng,lsampl_t maxdata);
	lsampl_t comedi_from_phys(double data,comedi_krange *rng,lsampl_t maxdata);
#endif
#endif
}
#endif /* !ACQUISITION_H */

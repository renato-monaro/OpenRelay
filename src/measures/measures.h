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
#ifndef MEASURES_H
#define MEASURES_H
#include <math.h>
#include <complex>
#include "channel.h"
#include "../control/control.h"
#include "parameters.h"
#include "types.h"
#include "filter.h"
#include <string>
#include <iostream>
#include <iomanip>
using namespace std;


//using namespace boost::posix_time;
namespace orelay{
enum {DYNAMIC=0,STATIC};
//enum {COMPLEX=0,REAL};

//! Classe Measure.
/** Esta é uma classe que é a base para as outras classes empregadas nos algoritmos de medição.*/
class Measure{
	public:
		virtual void Run()=0;  /*!< Método virtual pura necessária para a execução de medições.*/  
		virtual bool Prepare(float);
		void SetClock(Channel<timer> *); /*!< Definição do ponteiro de tempo para o algoritmo de medição.*/
		/*! Método que retorna uma string com a descrição funcional da medição.*/
		/** @return Retorna a descrição funcional da medição.*/
		string FunctionalDescription_get(void);
		void FunctionalDescription_set(string);/*!< Método que atribui a descrição funcional à medição.*/ 
	protected:
		Channel<timer> *MasterClock; /*!< Atributo que armazena o ponteiro de tempo.*/
		string FunctionalDescription;/*!< Atributo que armazena a descrição funcional.*/
		float RelaySamplingRate;
	};
//! Classe RMS.
/** Esta é uma classe que é derivada da classe Measure.*/
class RMS: public Measure{
	public:
		/*! Construtor da classe RMS.*/
		/** @param CH_IN É a entrada.
		@param CH_F É a frequência do sinal de entrada.
   		@param CH_OUT É a saída contendo o valor eficaz do canal de entrada.*/		
		RMS(Channel<analogic>* CH_IN,Channel<analogic>* CH_F,Channel<analogic>* CH_OUT); 
		/*! Construtor da classe RMS.*/
		/** @param CH_IN É a entrada.
   		@param CH_OUT É a saída contendo o valor eficaz do canal de entrada.*/		
		RMS(Channel<analogic>* CH_IN,Channel<analogic>* CH_OUT);
		//RMS(Channel<Vanalogic>* CH_IN,Channel<analogic>* CH_OUT);
		virtual void Run();  /*!< Método que executa a função de medida.*/
	protected:
		Channel<analogic> *Input; /*!< Canal de entrada o qual será utilizado para realizar a medição.*/ 
		Channel<analogic> *Freq; /*!< Canal de entrada o qual será utilizado para realizar a medição.*/ 
		Channel<analogic> *Output; /*!< Canal de saída o qual irá conter o valor eficaz do canal de entrada.*/ 
		Channel<Vanalogic> *VInput; /*!< Canal de saída o qual irá conter o valor eficaz do canal de entrada.*/ 
		bool vec;
	private:
		float Frequency;
		unsigned Metodology;
};
//! Classe Frequency.
/** Esta é uma classe que é derivada da classe Measure.*/
class Frequency: public Measure{
	public:
		/*! Construtor da classe Frequency.*/	
		/** @param CH_IN é o canal de entrada.
		@param CH_OUT é o canal de saída contendo o valor da frequência.
		@param MET é a metodologia que será escolhida para realizar a medição.*/
		Frequency(Channel<analogic>* CH_IN, Channel<analogic>* CH_OUT,int MET); 
		virtual void Run();  /*!< Método que executa a função de medida.*/
	protected:
		int Methodology; /*!< Atributo que determina a metodologia de cálculo.*/
		double average;
		Channel<analogic> *Input; /*< Canal de entrada, o qual será utilizado para realizar a medição.*/
		Channel<analogic> *Output; /*< Canal de saída, o qual irá conter o valor de frequência.*/
	private:
		timer CrossTime,CrossTime_old;
};
//! Classe Phasor
/** Esta é uma classe que é derivada da classe Measure.*/
class DFT: public Measure{
	public:
		/*! Construtor da classe DFT. */	
		/** @param CH_IN é o canal de entrada.
		@param CH_OUT é o canal de saída contendo os fasores.
		@param Freq é a frequência fundamental do sinal.
		@param Harm é o harmônic desejado.*/ 	
		DFT(Channel<analogic> *CH_IN,Channel<Complex> *CH_OUT,float Freq, unsigned Harm);
		/*! Construtor da classe DFT que calcula Transformada de Fourier da Primeria Harmônica. */	
		/** @param CH_IN é o canal de entrada.
		@param CH_OUT é o canal de saída contendo os fasores.
		@param Freq é a frequência fundamental do sinal.*/ 
		DFT(Channel<analogic> *CH_IN,Channel<Complex> *CH_OUT,float Freq);
		virtual void Run();  /*!< Método que executa a função de medida.*/
		bool Prepare(float);
	protected:	
		Channel<analogic>* Input; /*!< Canal de entrada, o qual será utilizado para realizar a medição.*/
		Channel<Complex>* Output; /*!< Canal de saída, o qual irá conter o conteúdo harmônico solicitado.*/	
		unsigned Harmonic; /*!< Harmônica que será calculada.*/	
		float SystemFrequency;
	private:
		unsigned N;
	};
//! Classe Power
/** Esta é uma classe que é derivada da classe Measure.*/
class Power: public Measure{
	public:	
 		/*! Construtor da classe Power para cálculo de potência monofásica.*/
/** @param CH_I1 é um canal de entrada que irá conter a corrente.
@param CH_V1 é o segundo canal de entrada que irá conter a tensão.
@param CH_P é o canal de saída que irá conter a potência ativa monofásica calculada.	
@param CH_Q é o canal de saída que irá conter a potência reativa monofásica calculada.
@param CH_S é o canal de saída que irá conter a potência aparente monofásica calculada.
@param Freq é a frequência fundamental do sistema. */
Power(Channel<analogic> *CH_I1,Channel<analogic> *CH_V1,Channel<analogic> *CH_S,Channel<analogic> *CH_P,Channel<analogic> *CH_Q, float Freq);
 		/*! Construtor da classe Power para cálculo de potência trifásica.*/	
/**@param CH_V1 é o canal de entrada que irá conter a primeira tensão trifásica.
@param CH_V2 é o canal de entrada que irá conter a segunda tensão trifásica.
@param CH_V3 é o canal de entrada que irá conter a terceira tensao trifásica.
@param CH_I1 é o canal de entrada que irá conter a primeira corrente trifásica.
@param CH_I2 é o canal de entrada que irá conter a segunda corrente trifásica.
@param CH_I3 é o canal de entrada que irá conter a terceira corrente trifásica.
@param CH_P é o canal de saída que irá conter a potência ativa trifásica calculada.	
@param CH_Q é o canal de saída que irá conter a potência reativa trifásica calculada.
@param CH_S é o canal de saída que irá conter a potência aparente trifásica calculada.
@param Freq é a frequência fundamental do sistema. */	
	Power(Channel<analogic> *CH_V1,Channel<analogic> *CH_V2,Channel<analogic> *CH_V3,Channel<analogic> *CH_I1,Channel<analogic> *CH_I2,Channel<analogic> *CH_I3,Channel<analogic> *CH_S,Channel<analogic> *CH_P,Channel<analogic> *CH_Q, float Freq);
		virtual void Run();   /*!< Método que executa a função de medida.*/
		bool Prepare(float);
	protected:
		vector<Channel<analogic>*> Current; /*!< Atribui em um vetor a(s) entrada(s) de tensão (ou corrente).*/
		vector<Channel<analogic>*> Voltage; /*!< Atribui em um vetor a(s) entrada(s) de corrente (ou tensão).*/
		Channel<analogic>* Active; /*!< Atribui o resultado do cálculo da potência ativa (de n fases equivalente ao número de entradas) em um canal de saída. */ 
		Channel<analogic>* Reactive; /*!< Atribui o resultado do cálculo da potência reativa (de n fases equivalente ao número de entradas) em um canal de saída. */ 
		Channel<analogic>* Apparent; /*!< Atribui o resultado do cálculo da potência aparente (de n fases equivalente ao número de entradas) em um canal de saída. */ 
		float SystemFrequency;
	private:
		unsigned N;
	};
	
class PowerDQ: public Measure{
	public:	
	PowerDQ(Channel<analogic> *CH_VD,Channel<analogic> *CH_ID,Channel<analogic> *CH_VQ,Channel<analogic> *CH_IQ,Channel<analogic> *CH_P, Channel<analogic> *CH_Q);
 		virtual void Run();   /*!< Método que executa a função de medida.*/
		bool Prepare(float);
	protected:
		Channel<analogic> *Vd,*Vq,*Iq,*Id; /*!< Atribui em um vetor a(s) entrada(s) de tensão (ou corrente).*/
		Channel<analogic> *P,*Q; /*!< Atribui em um vetor a(s) entrada(s) de corrente (ou tensão).*/
	};
	
class Power2: public Measure{
	public:	
	Power2(Channel<analogic> *CH_VA,Channel<analogic> *CH_VB,Channel<analogic> *CH_VC,Channel<analogic> *CH_IA, Channel<analogic> *CH_IB, Channel<analogic> *CH_IC,Channel<analogic> *CH_P, Channel<analogic> *CH_Q);
 		virtual void Run();   /*!< Método que executa a função de medida.*/
		bool Prepare(float);
	protected:
		Channel<analogic> *Va,*Vb,*Vc,*Ia,*Ib,*Ic; /*!< Atribui em um vetor a(s) entrada(s) de tensão (ou corrente).*/
		Channel<analogic> *P,*Q; /*!< Atribui em um vetor a(s) entrada(s) de corrente (ou tensão).*/
	};


//! Classe Sum.
/** Esta é uma classe que é derivada da classe Measure.*/
class Sum: public Measure{
	public:
		/*! Construtor da classe Sum.*/
		/**@param CH_OUT É a saída contendo a soma dos canais de entrada.*/		
		Sum(Channel<analogic>* CH_OUT); 
		void Join_Channel(Channel<analogic>*,double S); /*!< Método que adiciona um canal a ser somado.*/
		void Join_Channel(Channel<analogic>*); /*!< Método que adiciona um canal a*/
		virtual void Run();  /*!< Método que executa a função de soma.*/
		bool Prepare(float);
	protected:
		vector <Channel <analogic>* > Input; /*!< Vetor quer contém os canais de entrada.*/ 
		vector<analogic> Scalar;
		Channel<analogic> *Output; /*!< Canal de saída o qual irá conter a soma dos canais de entrada.*/ 
};

class Diff: public Sum{
	public:
		Diff(Channel<analogic>* CH_OUT) : Sum (CH_OUT) {};
		void Run();  /*!< Método que executa a função de diferença.*/
};

class Max: public Sum{
	public:
		Max(Channel<analogic>* CH_OUT)  : Sum (CH_OUT){};
		void Run();  /*!< Método que executa a função de diferença.*/
	};



class ZeroCrossDetection: public Measure{
	public:
		ZeroCrossDetection(Channel<analogic>* CH_IN, Channel<digital>* CH_OUT);
		void Run();
	protected:
		Channel<analogic>* Input;
		Channel<digital> *Output; /*!< Canal onde será indicada a passagem por zero.*/ 
};

//! Classe Mean.
/** Esta é uma classe que é derivada da classe Measure.*/
class Mean: public Measure{
	public:
		/*! Construtor da classe Mean.
		Calcula a média móvel para N amostras (N inferior ao tamanho do buffer)*/
		/** @param CH_IN É a entrada.
   		@param CH_OUT É a saída contendo a média móvel do canal de entrada.
	    @param N É o número de amostras na média.*/		
		Mean(Channel<analogic>* CH_IN,Channel<analogic>* CH_OUT,unsigned N); 
		Mean(Channel<analogic>* CH_IN,Channel<analogic>* CH_OUT,unsigned N,double S); 
		/*! Construtor da classe RMS.
		Calcula a média móvel para todo o Buffer do canal de entrada*/
		/** @param CH_IN É a entrada.
   		@param CH_OUT É a saída contendo a média móvel do canal de entrada.*/		
		Mean(Channel<analogic>* CH_IN,Channel<analogic>* CH_OUT);
		virtual void Run();  /*!< Método que executa a função de medida.*/
	protected:
		Channel<analogic> *Input; /*!< Canal de entrada o qual será utilizado para realizar a medição.*/ 
		Channel<analogic> *Output; /*!< Canal de saída o qual irá conter o valor eficaz do canal de entrada.*/ 
		unsigned N; /*!<Tamanho da janela da média móvel*/
		double Scale;
};


class Abs: public Measure{
	public:
	Abs(Channel<Complex> *CH_IN,Channel<analogic> *CH_OUT);	
	Abs(Channel<analogic> *CH_IN,Channel<analogic> *CH_OUT);
	void Run();
	protected:
		Channel<analogic> *Input; /*!< Canal de entrada o qual será utilizado para realizar a medição.*/
		Channel<Complex> *Input_C; /*!< Canal de entrada o qual será utilizado para realizar a medição.*/  
		Channel<analogic> *Output; /*!< Canal de saída o qual irá conter o valor eficaz do canal de entrada.*/ 
		unsigned Type;
};

class Arg: public Measure{
	public:
	Arg(Channel<Complex> *CH_IN,Channel<analogic> *CH_OUT);	
	void Run();
	protected:
		Channel<Complex> *Input_C; /*!< Canal de entrada o qual será utilizado para realizar a medição.*/  
		Channel<analogic> *Output; /*!< Canal de saída o qual irá conter o valor eficaz do canal de entrada.*/ 
};
class Real: public Arg{
	public:
	Real(Channel<Complex> *CH_IN,Channel<analogic> *CH_OUT):Arg(CH_IN,CH_OUT){};
	void Run();
};
class Imag: public Arg{
	public:
	Imag(Channel<Complex> *CH_IN,Channel<analogic> *CH_OUT):Arg(CH_IN,CH_OUT){};
	void Run();
};

class Conj: public Measure{
	public:
	Conj(Channel<Complex> *CH_IN,Channel<Complex> *CH_OUT);
	void Run();
	protected:
		Channel<Complex> *Input_C; /*!< Canal de entrada o qual será utilizado para realizar a medição.*/  
		Channel<Complex> *Output; /*!< Canal de saída o qual irá conter o valor eficaz do canal de entrada.*/ 
};

class DiffArg: public Measure{
	public:
	DiffArg(Channel<Complex> *CH_IN1,Channel<Complex> *CH_IN2,Channel<analogic> *CH_OUT);	
	virtual void Run();
	protected:
		Channel<Complex> *Input_C1; /*!< Canal de entrada o qual será utilizado para realizar a medição.*/
		Channel<Complex> *Input_C2; /*!< Canal de entrada o qual será utilizado para realizar a medição.*/    
		Channel<analogic> *Output; /*!< Canal de saída o qual irá conter o valor eficaz do canal de entrada.*/ 
};


class Scale: public Measure{
	public:
		Scale(Channel<analogic> *CH_IN,Channel<analogic> *CH_OUT, double G, double Off);
		void Run();
	protected:
		Channel<analogic> *Input; /*!< Canal de entrada o qual será utilizado para realizar a medição.*/
		Channel<analogic> *Output; /*!< Canal de saída o qual irá conter o valor eficaz do canal de entrada.*/ 
		double Gain,Offset;
};

class StringComposer: public Measure{
	public:
		StringComposer(Channel<string> *Str, string Sep);
		void Join(Channel<analogic> *ChAna, unsigned Dec,unsigned Width, string Prefix);
		void Join(Channel<string> *ChStr, string Prefix);
		void Run();
	protected:
		Channel<string> *Output;
		string Separator;
		vector<string> AnalogicPrefix, StringPrefix;
		vector<unsigned> AnalogicPosition, StringPosition, AnalogicDecimal, AnalogicWidth;
		vector<Channel <analogic> * > AnalogicList;
		vector<Channel <string> * > StringList;
	private:
		unsigned Position;
		stringstream Out;
};

//class WaveletDetector: public Measure{
//	public:
		/*! Construtor da classe WaveletDetector.
		Monitora um conjunto de canais analógicos e indica mudanças abruptas no sinal na energia do detalhe do último nível da transformada Wavelet. Pode ser usado para identificar faltas. É utilizado a família Daubechies.*/
		/** @param CH_OUT É a saída contendo a média móvel do canal de entrada.
		@param k é a ordem da transformada(k/2 momentos nulos).
		@param N é o nível de decomposição.
		@param SmplSize Tamanho da Janela para a TW.
		@param thr é o lim,iar de energia para a detecção.*/		
//		WaveletDetector(Channel<digital>* CH_OUT,unsigned k, unsigned N, unsigned SmplSize, double thr);
//		~WaveletDetector();
//		void Join_Channel(Channel<analogic>*); /*!< Método que adiciona um canal a ser somado.*/
//		void Run();
//	protected:
//		vector <Channel <analogic>* > Input; /*!< Vetor quer contém os canais de entrada.*/ 
//		Channel<digital> *Output; /*!< Canal onde será indicada a passagem por zero.*/ 
//		unsigned Order;/*!< Ordem da transformada.*/
//		unsigned Level;/*!< Nível de decomposição.*/
//		double	Threshold; /*!< Limiar de Ativação.*/
//		unsigned SampleSize; /*!< Tamanho da Janela para a TW.*/
//	private:
//		gsl_wavelet *w;
//		gsl_wavelet_workspace *work;
//		bool detect;
//		double *data;
//		double Energy;
//};

//! Classe Sum.
/** Esta é uma classe que é derivada da classe Measure.*/
class SumN: public Measure{
	public:
		/*! Construtor da classe Sum.*/
		/**@param CH_OUT É a saída contendo a soma dos canais de entrada.		
		@param CH_IN É a entrada.
		@param begin é a posição para inicio da soma
		@param end é a posiçãi do fim da soma*/		
		SumN(Channel<analogic>* CH_OUT,Channel<analogic>* CH_IN,unsigned begin, unsigned end); 
		/*! Construtor da classe Sum.*/
		/**@param CH_OUT É a saída contendo a soma dos canais de entrada.		
		@param CH_IN É a entrada.
		@param end é a posição do fim da soma que começa na posição 0 (mais atual)*/	
		SumN(Channel<analogic>* CH_OUT,Channel<analogic>* CH_IN, unsigned end); 
		virtual void Run();  /*!< Método que executa a função de soma.*/
	protected:
		unsigned Begin; /*!<Posição Inicial da Soma.*/
		unsigned End;/*!<Posição Final da Soma.*/
		Channel <analogic> *Input; /*!< Canal de entrada.*/ 
		Channel<analogic> *Output; /*!< Canal de saída o qual irá conter a soma dos canais de entrada.*/ 
};

class MeanN: public Measure{
	public:
		/*! Construtor da classe Sum.*/
		/**@param CH_OUT É a saída contendo a soma dos canais de entrada.		
		@param CH_IN É a entrada.
		@param begin é a posição para inicio da soma
		@param end é a posiçãi do fim da soma*/		
		MeanN(Channel<analogic>* CH_OUT,Channel<analogic>* CH_IN,unsigned begin, unsigned end); 
		/*! Construtor da classe Sum.*/
		/**@param CH_OUT É a saída contendo a soma dos canais de entrada.		
		@param CH_IN É a entrada.
		@param end é a posição do fim da soma que começa na posição 0 (mais atual)*/	
		MeanN(Channel<analogic>* CH_OUT,Channel<analogic>* CH_IN, unsigned end); 
		void Run();  /*!< Método que executa a função de soma.*/
	protected:
		unsigned Begin; /*!<Posição Inicial da Soma.*/
		unsigned End;/*!<Posição Final da Soma.*/
		Channel <analogic> *Input; /*!< Canal de entrada.*/ 
		Channel<analogic> *Output; /*!< Canal de saída o qual irá conter a soma dos canais de entrada.*/ 
	};

class Entropy: public Measure{
	public:
		/*! Construtor da classe Sum.*/
		/**@param CH_OUT É a saída contendo a soma dos canais de entrada.		
		@param CH_IN É a entrada.
		@param begin é a posição para inicio da soma
		@param end é a posiçãi do fim da soma*/		
		Entropy(Channel<analogic>* CH_OUT,Channel<analogic>* CH_IN,unsigned begin, unsigned end); 
		/*! Construtor da classe Sum.*/
		/**@param CH_OUT É a saída contendo a soma dos canais de entrada.		
		@param CH_IN É a entrada.
		@param end é a posição do fim da soma que começa na posição 0 (mais atual)*/	
		Entropy(Channel<analogic>* CH_OUT,Channel<analogic>* CH_IN, unsigned end); 
		void Run();  /*!< Método que executa a função de soma.*/
	protected:
		unsigned Begin; /*!<Posição Inicial da Soma.*/
		unsigned End;/*!<Posição Final da Soma.*/
		Channel <analogic> *Input; /*!< Canal de entrada.*/ 
		Channel<analogic> *Output; /*!< Canal de saída o qual irá conter a soma dos canais de entrada.*/ 
	};
//! Classe que calcula a Derivada Númerica.
class Derivative: public Measure{
	public:
		/*! Construtor da classe Derivative.*/
		/**@param CH_IN É a entrada.
		@param CH_OUT É a saída contendo a derivada do canal de entrada.*/		
		Derivative(Channel<analogic> *CH_OUT, Channel<analogic> *CH_IN);
		/*! Construtor da classe Derivative.*/
		/**@param CH_IN É a entrada.
		@param CH_OUT É a saída contendo a derivada do canal de entrada.
		@param N é o passo da Derivada (intervalo entre amostras)  \f$ CH_{OUT}=\frac{CH_{IN}(0)-CH_{IN}(N)}{N\Delta t}\f$.*/		
		Derivative(Channel<analogic> *CH_OUT, Channel<analogic> *CH_IN,unsigned N);
		void Run();
	protected:
		Channel <analogic> *Input; /*!< Canal de entrada.*/ 
		Channel<analogic> *Output; /*!< Canal de saída o qual irá conter a derivada do canal de entrada.*/ 
		unsigned N; /*!< Passo da Derivada*/

	};

//! Classe que calcula a Integral Númerica no tempo.
class Integrate: public Measure{
	public:
		/*! Construtor da classe Integrate.*/
		/**@param CH_IN É a entrada.
		@param CH_OUT É a saída contendo a integral no tempo do canal de entrada.*/		
		Integrate(Channel<analogic> *CH_OUT, Channel<analogic> *CH_IN);
		/*! Construtor da classe Integrate.*/
		/**@param CH_IN É a entrada.
		@param CH_OUT É a saída contendo a integral no tempo do canal de entrada.
		@param N é o passo da Integral (intervalo entre amostras).*/				
		Integrate(Channel<analogic> *CH_OUT, Channel<analogic> *CH_IN, unsigned N);
		void Run();
	protected:
		Channel <analogic> *Input; /*!< Canal de entrada.*/ 
		Channel<analogic> *Output; /*!< Canal de saída o qual irá conter a integral no tempo do canal de entrada.*/ 
		unsigned N; /*!< Passo da Integral*/
		analogic Acc; /*!< Acumulador da Integral*/

	};

class IRestOper: public Measure{
	public:
		IRestOper(Channel<Complex> *CH_IN_A,Channel<Complex> *CH_IN_B, Channel<analogic> *CH_DIFF, Channel<analogic> *CH_OPER);
		IRestOper(Channel<Complex> *CH_IN_A,Channel<Complex> *CH_IN_B, Channel<analogic> *CH_DIFF, Channel<analogic> *CH_OPER, double);
		void Run();
	protected:
		Channel <Complex> *Input_A, *Input_B;
		Channel<analogic> *Output_REST, *Output_OPER;
		double Base;
	};



class Compose: public Measure{
	public:
		Compose(Channel<Vanalogic> *CH_OUT);
		void Join_Channel(Channel<Vanalogic>*); /*!< Método que adiciona um canal a ser somado.*/
		void Run();
	protected:
		vector<Channel<Vanalogic>*> VInputs;
		Channel<Vanalogic> *Output;
	};

class Division: public Measure{
	public:
		Division(Channel<analogic> *CH_1,Channel<analogic> *CH_2, Channel<analogic> *CH_OUT);
		Division(Channel<Complex> *CH_1,Channel<Complex> *CH_2, Channel<Complex> *CH_OUT);
		void Run();
	protected:
		bool comp;
		Channel <analogic> *Input_1; /*!< Canal de entrada.*/
		Channel <analogic> *Input_2; /*!< Canal de entrada.*/  
		Channel<analogic> *Output; /*!< Canal de saída o qual irá conter a derivada do canal de entrada.*/ 

		Channel <Complex> *Input_1_C; /*!< Canal de entrada.*/
		Channel <Complex> *Input_2_C; /*!< Canal de entrada.*/  
		Channel<Complex> *Output_C; /*!< Canal de saída o qual irá conter a derivada do canal de entrada.*/ 
};

class Product: public Measure{
	public:
		Product(Channel<analogic> *CH_1,Channel<analogic> *CH_2, Channel<analogic> *CH_OUT);
		Product(Channel<Complex> *CH_1,Channel<Complex> *CH_2, Channel<Complex> *CH_OUT);
		void Run();
	protected:
		bool comp;
		Channel <analogic> *Input_1; /*!< Canal de entrada.*/
		Channel <analogic> *Input_2; /*!< Canal de entrada.*/  
		Channel<analogic> *Output; /*!< Canal de saída o qual irá conter a derivada do canal de entrada.*/ 

		Channel <Complex> *Input_1_C; /*!< Canal de entrada.*/
		Channel <Complex> *Input_2_C; /*!< Canal de entrada.*/  
		Channel<Complex> *Output_C; /*!< Canal de saída o qual irá conter a derivada do canal de entrada.*/ 
};

class ZeroSequence: public Measure{
public:
		ZeroSequence(Channel<Complex> *CH_1,Channel<Complex> *CH_2,Channel<Complex> *CH_3,Channel<Complex> *CH_OUT);
		void Run();
	protected:
		Channel<Complex> *Input_1; /*!< Canal de entrada.*/
		Channel<Complex> *Input_2; /*!< Canal de entrada.*/
		Channel<Complex> *Input_3; /*!< Canal de entrada.*/    
		Channel<Complex> *Output; /*!< Canal de saída o qual irá conter a derivada do canal de entrada.*/ 
		Complex alpha;
};

class NegativeSequence: public ZeroSequence{
public:
		NegativeSequence(Channel<Complex> *CH_1,Channel<Complex> *CH_2,Channel<Complex> *CH_3,Channel<Complex> *CH_OUT):ZeroSequence(CH_1,CH_2,CH_3,CH_OUT){};
		void Run();
};

class PositiveSequence: public ZeroSequence{
public:
		PositiveSequence(Channel<Complex> *CH_1,Channel<Complex> *CH_2,Channel<Complex> *CH_3,Channel<Complex> *CH_OUT):ZeroSequence(CH_1,CH_2,CH_3,CH_OUT){};
		void Run();
};

class Washout: public Measure{
public:
		Washout(Channel<analogic> *CH_1,Channel<analogic> *CH_OUT,double Tc);
		void Run();
		bool Prepare(float);
	protected:
		Channel<analogic> *Input_1; /*!< Canal de entrada.*/
		Channel<analogic> *Output; /*!< Canal de saída o qual irá conter a derivada do canal de entrada.*/ 
		double a,b,d;
		bool first_run;
};
class AlphaBeta: public Measure{
	public:
		AlphaBeta(Channel<analogic> *CH_1,Channel<analogic> *CH_2,Channel<analogic> *CH_3,Channel<analogic> *CH_4,Channel<analogic> *CH_5,Channel<analogic> *CH_6, bool Direction);
		AlphaBeta(Channel<analogic> *CH_1,Channel<analogic> *CH_2,Channel<analogic> *CH_3,Channel<analogic> *CH_4,Channel<analogic> *CH_5, bool Direction);
		void Run();
		bool Prepare(float);
	protected:
		Channel<analogic> *Input_1; /*!< Canal de entrada.*/
		Channel<analogic> *Input_2; /*!< Canal de entrada.*/
		Channel<analogic> *Input_3; /*!< Canal de entrada.*/    
		Channel<analogic> *Output_1; /*!< Canal de saída o qual irá conter a derivada do canal de entrada.*/ 
		Channel<analogic> *Output_2; /*!< Canal de saída o qual irá conter a derivada do canal de entrada.*/ 
		Channel<analogic> *Output_3; /*!< Canal de saída o qual irá conter a derivada do canal de entrada.*/
		bool WZero;
		bool Direction;
  private:
    double a,b,c,alpha,beta,gamma;
};
class DQ0: public Measure{
	public:
		DQ0(Channel<analogic> *CH_1,Channel<analogic> *CH_2,Channel<analogic> *CH_3,Channel<analogic> *CH_4,Channel<analogic> *CH_5,Channel<analogic> *CH_6,Channel<analogic> *CH_T, bool Direction);
		DQ0(Channel<analogic> *CH_1,Channel<analogic> *CH_2,Channel<analogic> *CH_3,Channel<analogic> *CH_4,Channel<analogic> *CH_5,Channel<analogic> *CH_T, bool Direction);
		void Run();
		bool Prepare(float);
	protected:
		Channel<analogic> *Theta;
		Channel<analogic> *Input_1; /*!< Canal de entrada.*/
		Channel<analogic> *Input_2; /*!< Canal de entrada.*/
		Channel<analogic> *Input_3; /*!< Canal de entrada.*/    
		Channel<analogic> *Output_1; /*!< Canal de saída o qual irá conter a derivada do canal de entrada.*/ 
		Channel<analogic> *Output_2; /*!< Canal de saída o qual irá conter a derivada do canal de entrada.*/ 
		Channel<analogic> *Output_3; /*!< Canal de saída o qual irá conter a derivada do canal de entrada.*/
		bool WZero;
		bool Direction;
  private:
    double a,b,c,d,q,z,ang;
};

class PLLdq: public Measure{
	public:
		PLLdq(Channel<analogic> *CH_1,Channel<analogic> *CH_2,Channel<analogic> *CH_3,Channel<analogic> *CH_T);
		void Run();
		bool Prepare(float);
		~PLLdq();
	protected:
    DQ0 *transformada;
    PI *control;
    Integrate *integral;
		Channel<analogic> *Theta;
		Channel<analogic> *Input_1; /*!< Canal de entrada.*/
		Channel<analogic> *Input_2; /*!< Canal de entrada.*/
		Channel<analogic> *Input_3; /*!< Canal de entrada.*/    
		Channel<analogic> *Vd; /*!< Canal de saída o qual irá conter a derivada do canal de entrada.*/ 
		Channel<analogic> *Vq; /*!< Canal de saída o qual irá conter a derivada do canal de entrada.*/ 
		Channel<analogic> *Ws; /*!< Canal de saída o qual irá conter a derivada do canal de entrada.*/ 
		double Inital_Value;
};

class ButterFilter: public Measure{
	public:
		ButterFilter(Channel<analogic> *CH_1,Channel<analogic> *CH_2, double F, long O, double A);
		ButterFilter(Channel<analogic> *CH_1,Channel<analogic> *CH_2, double F, long O, double A, double S);
		~ButterFilter();
		void Run();
		bool Prepare(float);
	protected:
		Channel<analogic> *Input; /*!< Canal de entrada.*/
		Channel<analogic> *Output; /*!< Canal de entrada.*/
		butterworth<analogic> *Filter;
		long Order;
		double Scale;
		double Cut_Freq, Atennuation;
	};
	
class Limiter: public Measure{
	public:
		Limiter(Channel<analogic> *CH_1,Channel<analogic> *CH_2, double lim_inf,double lim_sup);
		void Run();
		bool Prepare(float);
	protected:
		Channel<analogic> *Input; /*!< Canal de entrada.*/
		Channel<analogic> *Output; /*!< Canal de entrada.*/
		double Limit_Superior,Limit_Inferior;
	};
}
#endif /* !MEASURES_H */


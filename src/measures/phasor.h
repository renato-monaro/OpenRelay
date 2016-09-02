/*
# Copyright (C) 2008-2016 Rodolfo Varraschim Rocha
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

#ifndef PHASOR_H
#define PHASOR_H

#include <relay.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <iostream>
#include <fstream>

using namespace orelay;

//! Classe RDFT.
/** Esta é uma classe que é derivada da classe Measure.*/
class RDFT: public Measure{
	public:
	    /*! Construtor da classe RDFT que calcula o fasor da Harmônica especificada. */
		/** @param CH_IN é o canal de entrada.
		@param CH_OUT é o canal de saída contendo os fasores.
		@param Freq é a frequência fundamental do sinal de entrada.
		@param Harm é o harmônico desejado.*/
		RDFT(Channel<analogic> *CH_IN,Channel<Complex> *CH_OUT,float Freq, unsigned Harm);
		/*! Construtor da classe RDFT que calcula o fasor da Primeria Harmônica. */
		/** @param CH_IN é o canal de entrada.
		@param CH_OUT é o canal de saída contendo os fasores.
		@param Freq é a frequência fundamental do sinal de entrada.*/
		RDFT(Channel<analogic> *CH_IN,Channel<Complex> *CH_OUT,float Freq);
		bool Prepare(float);
		virtual void Run();  /*!< Método que executa a função de medida.*/
	protected:
		Channel<analogic>* Input; /*!< Canal de entrada, o qual será utilizado para realizar a medição.*/
		Channel<Complex>* Output; /*!< Canal de saída, o qual irá conter o conteúdo harmônico solicitado.*/
		unsigned Harmonic; /*!< Harmônica que será calculada.*/
		float SystemFrequency;
	private:
		unsigned N, j;
};

//! Classe LS.
/** Esta é uma classe que é derivada da classe Measure.*/
class LS: public Measure{
	public:
	    /*! Construtor da classe LS que calcula o fasor de uma Harmônica específica dentro de um espectro especificado. */
		/** @param CH_IN é o canal de entrada.
		@param CH_OUT é o canal de saída contendo os fasores.
		@param Freq é a frequência fundamental do sinal de entrada.
		@param Harm é o harmônico desejado.
		@param HarmRange é o Número de harmônicas desejado.*/
		LS(Channel<analogic> *CH_IN,Channel<Complex> *CH_OUT,float Freq, unsigned Harm, unsigned HarmRange);
		/*! Construtor da classe LS que calcula o fasor da Harmônica definida, até a quinta. */
		/** @param CH_IN é o canal de entrada.
		@param CH_OUT é o canal de saída contendo os fasores.
		@param Freq é a frequência fundamental do sinal.
		@param Harm é o harmônico desejado.*/
		LS(Channel<analogic> *CH_IN,Channel<Complex> *CH_OUT,float Freq, unsigned Harm);
		/*! Construtor da classe LS que calcula o fasor da Primeria Harmônica. */
		/** @param CH_IN é o canal de entrada.
		@param CH_OUT é o canal de saída contendo os fasores.
		@param Freq é a frequência fundamental do sinal.*/
		LS(Channel<analogic> *CH_IN,Channel<Complex> *CH_OUT,float Freq);
        bool Prepare(float);
		virtual void Run();  /*!< Método que executa a função LS.*/
	protected:
		Channel<analogic>* Input; /*!< Canal de entrada, o qual será utilizado para realizar a medição.*/
		Channel<Complex>* Output; /*!< Canal de saída, o qual irá conter o conteúdo harmônico solicitado.*/
		unsigned Harmonic; /*!< Harmônica que será calculada.*/
		float SystemFrequency;
		unsigned SpecRange; /*!< Número de harmônicas consideradas.*/
	private:
		unsigned N,M,k,recomp,imcomp,col;
		float nsamp;
	};

//! Classe PLL.
/** Esta é uma classe que é derivada da classe Measure.*/
class PLL: public Measure{
	public:
	    /*! Construtor da classe PLL que calcula o fasor da primeira harmônica e a frequência. */
		/** @param CH_IN é o canal de entrada.
		@param CH_OUT é o canal de saída contendo os fasores.
		@param CH_FREQ é o canal de saída contendo a frequência estimada do sinal.
		@param Freq é a frequência fundamental do sinal. Valor inicial para a estimação*/
		PLL(Channel<analogic> *CH_IN,Channel<Complex> *CH_OUT,Channel<analogic> *CH_FREQ,float Freq);
		bool Prepare(float);
		virtual void Run();  /*!< Método que executa a função PLL.*/
	protected:
		Channel<analogic>* Input; /*!< Canal de entrada, o qual será utilizado para realizar a medição.*/
		Channel<Complex>* Output; /*!< Canal de saída, o qual irá conter o fasor do sinal de entrada.*/
		Channel<analogic>* Frequency; /*!< Canal de saída, o qual irá conter a frequência estimada do sinal de entrada.*/
		float SystemFrequency;
	private:
		unsigned N,j;
		float u1,u2,u3,C,A_old,A_new,W_old,W_new,Fi_old,Fi_new,
		Erro,deltat,Fasor_mod,Fasor_teta,Freq_med,sample,sig_pll;
};

//! Classe RWT.
/** Esta é uma classe que é derivada da classe Measure.*/
class RWT: public Measure{
	public:
	    /*! Construtor da classe RWT que calcula o fasor da frequência fundamental. */
		/** @param CH_IN é o canal de entrada.
		@param CH_OUT é o canal de saída contendo os fasores.
		@param Freq é a frequência fundamental do sinal de entrada.*/
		RWT(Channel<analogic> *CH_IN,Channel<Complex> *CH_OUT,float Freq);
		bool Prepare(float);
		virtual void Run();  /*!< Método que executa a função de medida.*/
	protected:
		Channel<analogic>* Input; /*!< Canal de entrada, o qual será utilizado para realizar a medição.*/
		Channel<Complex>* Output; /*!< Canal de saída, o qual irá conter o fasor do sinal de entrada.*/
		float SystemFrequency;
	private:
		unsigned N,k;
		Complex lam1,lam2,lam3,lam4,lam5,beta1,beta2,beta3,beta4,beta5,beta6,
                alfa,alfa_par,n1,I,I_arg1,I_arg2,I_arg3,W,W_1,W_2,W_3,W_4,W_5,W_6;
        float sig,w0,deltat,f,Am,Fm,temp,temp2;
};

//! Classe RLS.
/** Esta é uma classe que é derivada da classe Measure.*/
class RLS: public Measure{
	public:
	    /*! Construtor da classe RLS que calcula o fasor da primeira Harmônica. */
		/** @param CH_IN é o canal de entrada.
		@param CH_OUT é o canal de saída contendo os fasores.
		@param Freq é a frequência fundamental do sinal de entrada.*/
		RLS(Channel<analogic> *CH_IN,Channel<Complex> *CH_OUT,float Freq);
		bool Prepare(float);
		virtual void Run();  /*!< Método que executa a função de medida.*/
	protected:
		Channel<analogic>* Input; /*!< Canal de entrada, o qual será utilizado para realizar a medição.*/
		Channel<Complex>* Output; /*!< Canal de saída, o qual irá conter o fasor do sinal de entrada.*/
		float SystemFrequency;
	private:
		unsigned N,k;
		float Ampi,Teta,deltat,fii;
		double ForgetingFactor,e;
        double Roi_acum[2][2], aoi_acum[2][1];
};

//! Classe DWT.
/** Esta é uma classe que é derivada da classe Measure.*/
class DWT: public Measure{
    public:
        /*! Construtor da classe DWT que calcula o fasor da frequência fundamental. */
		/** @param CH_IN é o canal de entrada.
		@param N é o número de amostras/ciclo do sinal.
		@param Family é a função wavelet a ser usada.*/
        DWT(Channel<analogic> *CH_IN,Channel<Complex> *CH_OUT, unsigned N, unsigned Family);
        /*! Construtor da subclasse Leafs associada à DWT para calculo do fasor da frequência fundamental. */
		/** @param ad Define se aproximação ou detalhe é desejado.
		@param Nad É a folha desejada para o nivel de aproximação especificado.
		@param Level É o nivel de interesse para decomposição wavelet.*/
        void Leafs(unsigned char ad,unsigned Nad,unsigned Level);
        virtual void Run(); /*!< Método que executa a função de medida.*/
        bool Prepare(float);
    protected:
        bool getFilterBank();
        Channel <analogic> *Input; /*!< Canal de entrada, o qual será utilizado para realizar a medição.*/
        Channel <Complex> *Output; /*!< Canal de saída, o qual conterá os fasores da primeira harmônica.*/
        vector<double> FilterBank;
        vector<double> AR1;
        vector<double> AR2;
        double mAR1, mAR2;
        unsigned NSamples;
        unsigned FamilyName;
        unsigned Method;
        vector<unsigned char> AD;
        vector<unsigned> Nivel,Pos;
        vector<unsigned> Ai,Di,Af,Df;
        vector<unsigned> PosIni,PosFim;
        vector<bool> Ordem;
        bool neworder;
        void Transformada_Wavelet(double* f, long n,bool ordem);
};
#endif

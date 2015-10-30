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

#ifndef GENERATOR_H
#define GENERATOR_H
#include "channel.h"
#include "parameters.h"
#include "types.h"
namespace orelay{
class IDesbalanceadasG30: public Protection{
	public:
	IDesbalanceadasG30(Channel<Complex>* In1, Channel<Complex>* In2, Channel<Complex>* In3, Channel<digital>* Out, Channel<analogic>* Out2, double In, double Tmin, double Tmax, double P, double constK);
	bool Prepare(float SamplingFreq);
	void Run();
		protected:
		double PickUp; /*!< Pick-up*/
		double tmin, tmax; /*!< Tempos de atuação mínimo e máximo (s)*/
		double IntegralTime; /*!< Atributo que armazena o somatório do tempo durante cada execução.*/
		Channel<Complex> *Input1; /*!< Canal de entrada.*/
		Channel<Complex> *Input2; /*!< Canal de entrada.*/
		Channel<Complex> *Input3; /*!< Canal de entrada.*/
		Channel<digital> *Output; /*!< Canal de saída, o qual define o status atual da proteção.*/
		Channel<analogic> *Output2; /*!< Canal de saída, valor de I2.*/
		float dt; /*!< Tempo entre cada amostra.*/
	private:
		double INorm, K, t, I2calc;
	};
	
class DifferentialG30: public Protection{
	public:
DifferentialG30(Channel<Complex>* In1, Channel<Complex>* In2, Channel<digital>* Out, double In, double P1, double B1, double S1, double B2, double S2);
	bool Prepare(float SamplingFreq);
	void Run();
	double Spline();
	protected:
		double PickUp1; /*!< Pick-up baixo (pu)*/
		double PickUp2; /*!< Pick-up alto (pu)*/
		double BreakPoint1; /*!< Ponto base onde a primeira declividade  cruza o eixo ordenado em pu da curva de proteção.*/
		double Slope1; /*!< Primeira declividade da curva de proteção dada em porcentagem (0-100%).*/
		double BreakPoint2; /*!< Ponto base onde a segunda declividade  cruza o eixo ordenado em pu da curva de proteção.*/
		double Slope2; /*!< Segunda declividade da curva de proteção dada em porcentagem (0-100%).*/
		Channel<Complex> *Input1; /*!< Canal de entrada.*/
		Channel<Complex> *Input2; /*!< Canal de entrada.*/
		Channel<digital> *Output; /*!< Canal de saída, o qual define o status atual da proteção.*/
	private:
		double Break1, Break2, Break3, Ir,Io;
		double a1,b1,a2,b2,INorm,t,p0,m0,p1,m1;
};

class MotorizacaoG30: public Protection{
	public:
		MotorizacaoG30(Channel<analogic>* In1, Channel<analogic>* In2, Channel<analogic>* In3, Channel<analogic>* In4, Channel<analogic>* In5, Channel<analogic>* In6, Channel<analogic>* Out1, Channel<digital>* Out2, double Snom, double TPprim, double TPsec, double TCsec, double P, double D);
		bool Prepare(float SamplingFreq);
		void Run();
	protected:
		double Pickup; /*!< Ajustes*/
		double Temp; /*!< Delay (s)*/
		double IntegralTime; /*!< Atributo que armazena o somatório do tempo durante cada execução.*/
		Channel<analogic> *Input1; /*!< Canal de entrada.*/
		Channel<analogic> *Input2; /*!< Canal de entrada.*/
		Channel<analogic> *Input3; /*!< Canal de entrada.*/
		Channel<analogic> *Input4; /*!< Canal de entrada.*/
		Channel<analogic> *Input5; /*!< Canal de entrada.*/
		Channel<analogic> *Input6; /*!< Canal de entrada.*/
		Channel<analogic> *Output1; /*!< Canal de saída, o qual informa a potência ativa calculada.*/
		Channel<digital> *Output2; /*!< Canal de saída, o qual define o status atual da proteção.*/
		float dt; /*!< Tempo entre cada amostra.*/

	private:
	    double SNorm, RTP, RTC, PotCalc;      /*!< Valores nominais.*/
};
class PerdaCampoG30: public Protection{
	public:
PerdaCampoG30(Channel<Complex>* In1, Channel<Complex>* In2, Channel<digital>* Out, double C1, double R1, double C2, double R2, double D);
		bool Prepare(float SamplingFreq);
		void Run();
protected:
		double Centro1, Centro2, Raio1, Raio2; /*!< Ajustes (Ohms)*/
		double Temp; /*!< Delay (s)*/
		double IntegralTime; /*!< Atributo que armazena o somatório do tempo durante cada execução.*/
		Channel<Complex> *Input1; /*!< Canal de entrada.*/
		Channel<Complex> *Input2; /*!< Canal de entrada.*/
		Channel<digital> *Output; /*!< Canal de saída, o qual define o status atual da proteção.*/
		float dt; /*!< Tempo entre cada amostra.*/

	private:
};

class SobreexcitacaoG30: public Protection{
	public:
SobreexcitacaoG30(Channel<analogic>* In1, Channel<Complex>* In2, Channel<analogic>* In3, Channel<Complex>* In4, Channel<analogic>* In5, Channel<Complex>* In6, Channel<digital>* Out, double Vn, double Fn, double P, double D);
		bool Prepare(float SamplingFreq);
		void Run();
protected:
		double PickUp; /*!< Pick-up (pu)*/
		double Temp; /*!< Delay (s)*/
		double IntegralTime; /*!< Atributo que armazena o somatório do tempo durante cada execução.*/
		Channel<analogic> *Input1; /*!< Canal de entrada.*/
		Channel<Complex> *Input2; /*!< Canal de entrada.*/
		Channel<analogic> *Input3; /*!< Canal de entrada.*/
		Channel<Complex> *Input4; /*!< Canal de entrada.*/
		Channel<analogic> *Input5; /*!< Canal de entrada.*/
		Channel<Complex> *Input6; /*!< Canal de entrada.*/
		Channel<digital> *Output; /*!< Canal de saída, o qual define o status atual da proteção.*/
		float dt; /*!< Tempo entre cada amostra.*/

	private:
		double VHz;              /*!< Valor de V/Hz.*/
		double VNorm,FNorm;      /*!< Valores nominais.*/
};

class SobrfrequenciaG30: public Protection{
	public:
SobrfrequenciaG30(Channel<analogic>* In1, Channel<Complex>* In2, Channel<analogic>* In3, Channel<Complex>* In4, Channel<analogic>* In5, Channel<Complex>* In6, Channel<digital>* Out, double Vn, double MinV, double P, double D);
		bool Prepare(float SamplingFreq);
		void Run();
	protected:
		double PickUp; /*!< Pick-up*/
		double Temp; /*!< Delay (s)*/
		double IntegralTime; /*!< Atributo que armazena o somatório do tempo durante cada execução.*/
		Channel<analogic> *Input1; /*!< Canal de entrada.*/
		Channel<Complex> *Input2; /*!< Canal de entrada.*/
		Channel<analogic> *Input3; /*!< Canal de entrada.*/
		Channel<Complex> *Input4; /*!< Canal de entrada.*/
		Channel<analogic> *Input5; /*!< Canal de entrada.*/
		Channel<Complex> *Input6; /*!< Canal de entrada.*/
		Channel<digital> *Output; /*!< Canal de saída, o qual define o status atual da proteção.*/
		float dt; /*!< Tempo entre cada amostra.*/

	private:
		double VHz;              /*!< Valor de V/Hz.*/
		double VNorm,VMin;      /*!< Valores nominais.*/
};

class SubfrequenciaG30: public Protection{
	public:
SubfrequenciaG30(Channel<analogic>* In1, Channel<Complex>* In2, Channel<analogic>* In3, Channel<Complex>* In4, Channel<analogic>* In5, Channel<Complex>* In6, Channel<digital>* Out, double Vn, double MinV, double P, double D);
		bool Prepare(float SamplingFreq);
		void Run();
		protected:
		double PickUp; /*!< Pick-up*/
		double Temp; /*!< Delay (s)*/
		double IntegralTime; /*!< Atributo que armazena o somatório do tempo durante cada execução.*/
		Channel<analogic> *Input1; /*!< Canal de entrada.*/
		Channel<Complex> *Input2; /*!< Canal de entrada.*/
		Channel<analogic> *Input3; /*!< Canal de entrada.*/
		Channel<Complex> *Input4; /*!< Canal de entrada.*/
		Channel<analogic> *Input5; /*!< Canal de entrada.*/
		Channel<Complex> *Input6; /*!< Canal de entrada.*/
		Channel<digital> *Output; /*!< Canal de saída, o qual define o status atual da proteção.*/
		float dt; /*!< Tempo entre cada amostra.*/

	private:
		double VHz;              /*!< Valor de V/Hz.*/
		double VNorm,VMin;      /*!< Valores nominais.*/
};
}
#endif	

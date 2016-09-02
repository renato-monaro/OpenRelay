/*
# Copyright (C) 2008-2015 Rodolfo Varraschim Rocha
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

#ifndef DDF_H
#define DDF_H

#include <relay.h>
#include <iostream>
#include <fstream>

using namespace orelay;

//! Classe DDF.
/** Filtra componentes de decaimento exponencial*/
class DDF: public Measure{
	public:
	    /*! Construtor da classe DDF.*/
		/** @param CH_IN É o canal de entrada entrada.
   		@param CH_OUT É a saída filtrada, sem componente de decaimento, do canal de entrada.
   		@param Freq É a frequência fundamental do sinal de entrada, deve ser fornecida.*/
		DDF(Channel<analogic> *CH_IN,Channel<analogic> *CH_OUT,float Freq);
		bool Prepare(float);
		virtual void Run();  /*!< Método que executa a função de filtragem.*/
	protected:
		Channel<analogic>* Input; /*!< Canal de entrada, o qual será utilizado para realizar a medição.*/
		Channel<analogic>* Output; /*!< Canal de saída, o qual irá conter o canal de entrada filtrado.*/
		float SystemFrequency;
	private:
		unsigned N;
		float b[3],a[1],z[3];
};
#endif

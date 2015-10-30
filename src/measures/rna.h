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
#ifndef RNA_H
#define RNA_H
#ifdef RNA
#include <fann.h>
#include <vector>
#include "measures.h"
#include "channel.h"  
namespace orelay{
class RNA: public Measure{
	public:
		RNA(struct fann *);
		RNA(char*);
		~RNA();
		void Join_InputChannel(Channel<analogic>*); /*!< Método que adiciona um canal de entrada para a Rede neural.*/
		void Join_InputChannel(Channel<analogic>*, analogic scale); /*!< Método que adiciona um canal de entrada para a Rede neural.*/
		void Join_InputChannel(Channel<analogic>*,unsigned); /*!< Método que adiciona um canal de entrada com atraso para a Rede neural.*/
		void Join_InputChannel(Channel<analogic>*,unsigned,analogic scale); /*!< Método que adiciona um canal de entrada com atraso para a Rede neural.*/
		void Join_InputChannel(Channel<analogic>*,unsigned,unsigned); /*!< Método que adiciona um canal de entrada com um conjunto de atrasos para a Rede neural.*/
		void Join_InputChannel(Channel<analogic>*,unsigned,unsigned, analogic scale); /*!< Método que adiciona um canal de entrada com um conjunto de atrasos para a Rede neural.*/
		void Join_OutputChannel(Channel<analogic>*); /*!< Método que adiciona um canal de saida para a Rede neural.*/
		void Join_OutputChannel(Channel<analogic>*, analogic scale); /*!< Método que adiciona um canal de saida para a Rede neural.*/
		void Join_OutputChannel(Channel<digital>*); /*!< Método que adiciona um canal de saida para a Rede neural.*/
		virtual void Run();  /*!< Método que executa a RNA*/
		virtual bool Prepare(float);  /*!< Método que executa a RNA*/
	protected:
		vector <Channel <analogic>* > Input; /*!< Vetor quer contém os canais de entrada.*/ 
		vector <unsigned> Delay;
		vector <analogic> InputScale;
		vector <analogic> OutputScale;
		unsigned OutputType;
		vector <Channel <analogic>* > AnalogicOutput; /*!< Vetor quer contém os canais de entrada.*/ 
		vector <Channel <digital>* > DigitalOutput; /*!< Vetor quer contém os canais de entrada.*/ 
		struct fann *ann;
	private:
		fann_type *output;
    		fann_type *input;
};
}
#endif
#endif/* !RNA_H */

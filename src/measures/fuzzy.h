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
#ifndef FUZZY_H
#define FUZZY_H

#ifdef FUZZYLITE
#include <fl/Headers.h>
#include <vector>
#include "measures.h"
#include "channel.h"  
#include "types.h"  
#include <fstream>
#include <iostream>


namespace orelay{

class Fuzzy: public Measure{
	public:
		Fuzzy(fl::Engine *);
		//Fuzzy(fl::Engine *, unsigned);
		Fuzzy(const char*,unsigned);
		~Fuzzy();
		void Join_InputChannel(Channel<analogic>*); /*!< Método que adiciona um canal de entrada para o sistema Fuzzy.*/
		//void Join_InputChannel(Channel<analogic>*, analogic, analogic); /*!< Método que adiciona um canal de entrada para o sistema Fuzzy.*/
		void Join_OutputChannel(Channel<analogic>*); /*!< Método que adiciona um canal de saida para o sistema Fuzzy.*/
		virtual void Run();  /*!< Método que executa o FIS*/
		bool Prepare(float);  /*!< Método que executa o FIS*/
	protected:
		//bool Discretize();
		vector <Channel <analogic>* > Input; /*!< Vetor quer contém os canais de entrada.*/ 
		vector <Channel <analogic>* > Output; /*!< Vetor quer contém os canais de entrada.*/ 
		fl::Engine *engine;
		vector<fl::InputVariable*> InputL;
		vector<fl::OutputVariable*> OutputL;
		vector<fl::RuleBlock*> RuleBlock;
		unsigned Discretization;
		vector<vector<double> > DiscreteTables;
		vector<analogic> MaxIn;
		vector<analogic> MinIn;
	private:

};
}
#endif
#endif/* !FUZZY_H */

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
#ifndef SWITCH_H
#define SWITCH_H

#include "channel.h"
#include "measures.h"
#include "parameters.h"
#include "types.h"
#include <string>
namespace orelay{
enum {T_CHANNEL,T_FIXVALUE};

class Switch: public Measure{
	public:
		Switch(Channel<analogic> *CH_1,Channel<analogic> *CH_2,Channel<analogic> *CH_OUT, Channel<digital> *CH_S);
		Switch(double,Channel<analogic> *CH_1,Channel<analogic> *CH_OUT, Channel<digital> *CH_S);
		void Run();
	protected:
		Channel<analogic> *Input_1; /*!< Canal de entrada o qual será utilizado para realizar a medição.*/
		Channel<analogic> *Input_2; /*!< Canal de entrada o qual será utilizado para realizar a medição.*/ 
		Channel<digital> *Input_D; /*!< Canal de entrada o qual será utilizado para realizar a medição.*/ 
		Channel<analogic> *Output; /*!< Canal de saída o qual irá conter o valor eficaz do canal de entrada.*/ 
		double FixValue;
		bool withfix;
	};

class SwitchPos: public Measure{
	public:
		SwitchPos(Channel<analogic> *CH_OUT, Channel<digital> *CH_S);
		bool Join(Channel<analogic> *CH_1);
		bool Join(double);
		void Run();
	protected: 
		double getValue();
		Channel<digital> *Input_D; /*!< Canal de entrada o qual será utilizado para realizar a medição.*/ 
		Channel<analogic> *Output; /*!< Canal de saída o qual irá conter o valor eficaz do canal de entrada.*/ 
		vector<unsigned> type;
		vector<unsigned> index;
		vector<Channel<analogic> *> vecChannel;
		vector<double> vecfixvalue;
		unsigned pos;
	};
}
#endif

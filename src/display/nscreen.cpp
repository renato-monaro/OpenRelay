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
#include "nscreen.h"
#ifdef NCURSES
using namespace orelay;
void Display::SetClock(Channel<timer> *t){
	MasterClock=t;
	}
void NScreen::Join_Channel( Channel<analogic> *Ch){
	Analogic.push_back(Ch);
	}
void NScreen::Join_Channel( Channel<digital> *Ch){
	Digital.push_back(Ch);
	}
void NScreen::Join_Channel( Channel<timer> *Ch){
	Timer.push_back(Ch);
	}

void NScreen::Join_Channel( Channel<Complex> *Ch){
	Phasor.push_back(Ch);
	}
void NScreen::Join_Channel( Channel<string> *Ch){
	String.push_back(Ch);
	}
/*
void NScreen::Join_Channel( Channel<power> *Ch){
	Power.push_back(Ch);
	}*/
NScreen::NScreen(){
 	initscr();
 	getmaxyx(stdscr,MaxRow,MaxCol);	
	OldCout = cout.rdbuf(); // save original sbuf
    cout.rdbuf(OutBuffer.rdbuf()); // redirect 'cout' to a 'fout'
	}
NScreen::~NScreen(){
	endwin();	
	cout.rdbuf(OldCout); // restore the original stream buffer
	}

void NScreen::Run(){
	while(1){
		#ifdef DEBUG
			cout<<"NScreen\n" ;
		#endif
//		erase();
		move(0,0);
		printw("Time:%lld\n",MasterClock->get_Value());
		for(unsigned j=0;j<Timer.size();j++)
			printw("%s:%lld\n",Timer[j]->get_Name().c_str(),Timer[j]->get_Value());
		for(unsigned j=0;j<Analogic.size();j++)
			printw("%s:%5.3lf\n",Analogic[j]->get_Name().c_str(),Analogic[j]->get_Value());
		for(unsigned j=0;j<Phasor.size();j++){
			printw("%s:%5.3lf | %5.3lf\n",Phasor[j]->get_Name().c_str(),abs(Phasor[j]->get_Value()),arg(Phasor[j]->get_Value())-arg(Phasor[0]->get_Value()));
			}
		for(unsigned j=0;j<String.size();j++){
			printw("%s:%s\n",String[j]->get_Name().c_str(),String[j]->get_Value().c_str());
			}
/*
		for(unsigned j=0;j<Power.size();j++){
			printw("%s:%5.2lf(VA)\t%5.2lf(W)\t%5.2lf(VAr)\t%5.2lf(VAh)\n",Power[j]->get_Name().c_str(),
				Power[j]->get_Value().get_Apparent(),
				Power[j]->get_Value().get_Active(),
				Power[j]->get_Value().get_Reactive(),
				Power[j]->get_Value().get_Distortion());
			}*/
		for(unsigned j=0;j<Digital.size();j++)
			printw("%s:%d\n",Digital[j]->get_Name().c_str(),Digital[j]->get_Value());

		move(MaxRow-4,0);
		printw("%s",OutBuffer.str().c_str());
		OutBuffer.str("");
		refresh();
		boost::this_thread::sleep(boost::posix_time::milliseconds(REFRESH_TIME));
		}
	}

#endif



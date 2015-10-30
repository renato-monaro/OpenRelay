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
#include "fuzzy.h"


#ifdef FUZZY
using namespace orelay;
int myPow(int x, int p)
{
  if (p == 0) return 1;
  if (p == 1) return x;

  int tmp = myPow(x, p/2);
  if (p%2 == 0) return tmp * tmp;
  else return x * tmp * tmp;
}

//Fuzzy::Fuzzy(char *File){}
/*Fuzzy::Fuzzy(fl::Engine *Eng){
	string tmp;
	tmp="Fuzzy";
	FunctionalDescription_set(tmp); 
	engine=Eng;
	}*/

Fuzzy::Fuzzy(fl::Engine *Eng){
	string tmp;
	tmp="Fuzzy";
	FunctionalDescription_set(tmp); 
	engine=Eng;
	}

Fuzzy::Fuzzy(const char* file, unsigned F){
	engine = new fl::Engine;
	ifstream in;
	in.open(file);
	if (in.fail())
		exit(0);
	std::string fclEngine((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
	in.close();

	fl::Importer* importer;
	switch(F){
		case FIS:
			importer = new fl::FisImporter;
			break;
		case FCL:
			importer = new fl::FclImporter;
			break;
		case FLL:
			importer = new fl::FllImporter;
			break;
		default:
			importer = new fl::FclImporter;	
		}
	engine = importer->fromString(fclEngine);
	}
/*void Fuzzy::Join_InputChannel(Channel<analogic> *In, analogic Max,analogic Min){
	Input.push_back(In);
	MaxIn.push_back(Max);
	MinIn.push_back(Min);
	}*/

void Fuzzy::Join_InputChannel(Channel<analogic> *In){
	Input.push_back(In);
	cout<<"Join In"<<endl;
	}
void Fuzzy::Join_OutputChannel(Channel<analogic> *Out){
	Output.push_back(Out);
	cout<<"Join Out"<<endl;
	}

bool Fuzzy::Prepare(float sampling){
	cout<<"Prepare Fuzzy"<<endl;
	if(!engine){
		cout<<"Fuzzy:FIS Not Allocated."<<endl; 
		return false;		
		}
	cout<<"Fuzzy Name: "<<engine->getName()<<endl;
	if(engine->numberOfOutputVariables()!=Output.size()){
		cout<<"Fuzzy:Number of Outputs ("<<Output.size()<<") Incompatible with FIS("<<engine->numberOfOutputVariables()<<")"<<endl;  
		return false;	
		}
	for(unsigned i=0;i<engine->numberOfOutputVariables();i++)
		OutputL.push_back(engine->getOutputVariable(i));

	if(engine->numberOfInputVariables()!=Input.size()){
		cout<<"Fuzzy:Number of Inputs ("<<Input.size()<<") Incompatible with FIS("<<engine->numberOfInputVariables()<<")"<<endl; 
		return false;
		}
	for(unsigned i=0;i<engine->numberOfInputVariables();i++)
		InputL.push_back(engine->getInputVariable(i));
/*	if(!Discretize())
		return false;*/
        fl::Exporter* exporter =new fl::FclExporter;
	cout<<exporter->toString(engine)<<endl;
	return true;
	}

void Fuzzy::Run(){
	unsigned m=0;
	for(unsigned i=0;i<Input.size();i++)
		InputL[i]->setInputValue(Input[i]->get_Value());
	engine->process();
	for(unsigned i=0;i<Output.size();i++)
		Output[i]->insert_Value(OutputL[i]->getOutputValue());
	}

Fuzzy::~Fuzzy(){
	}

/*bool Fuzzy::Discretize(){
	vector<vector<double> > InputDiscrete;
	InputDiscrete.resize(Input.size());
	DiscreteTables.resize(Output.size());
//Normalização
	for(unsigned k=0;k<Input.size();k++){
		for(unsigned j=0;j<Discretization;j++){
			InputDiscrete[k].push_back(((InputL[k]->maximum()-InputL[k]->minimum())/Discretization)*j);
			}
		}

	for(unsigned i=0;i<myPow(Discretization,Input.size());i++){
		for(unsigned k=0;k<Input.size();k++)
			InputL[k]->setInput(InputDiscrete[k][(i/myPow(Discretization,k))%Discretization]);
		engine->process();
		for(unsigned l=0;l<Output.size();l++)
			DiscreteTables[l].push_back(OutputL[l]->output().defuzzify());
		}
	for(unsigned i=0;i<myPow(Discretization,Input.size());i++)
		cout<<DiscreteTables[0][i]<<" ";
	return true;
	}*/
#endif

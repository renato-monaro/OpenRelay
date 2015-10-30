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
#ifdef RNA
#include "rna.h"
using namespace orelay;
RNA::RNA(char *File){
	FunctionalDescription_set("RNA"); 
	ann=fann_create_from_file(File);
	OutputType=NONE;
	}
RNA::RNA(struct fann *A){
	FunctionalDescription_set("RNA"); 
	ann=A;
	OutputType=NONE;
	}
void RNA::Join_InputChannel(Channel<analogic> *In){
	Input.push_back(In);
	Delay.push_back(0);
	InputScale.push_back(1.0);
	}
void RNA::Join_InputChannel(Channel<analogic> *In,unsigned Del){
	Input.push_back(In);
	Delay.push_back(Del);
	InputScale.push_back(1.0);
	}
void RNA::Join_InputChannel(Channel<analogic> *In,unsigned Beg, unsigned End){
	for(unsigned i=Beg;i<=End;i++){
		Input.push_back(In);
		Delay.push_back(i);
		InputScale.push_back(1.0);
		}
	}

void RNA::Join_InputChannel(Channel<analogic> *In, analogic scale){
	Input.push_back(In);
	Delay.push_back(0);
	InputScale.push_back(scale);
	}
void RNA::Join_InputChannel(Channel<analogic> *In,unsigned Del, analogic scale){
	Input.push_back(In);
	Delay.push_back(Del);
	InputScale.push_back(scale);
	}
void RNA::Join_InputChannel(Channel<analogic> *In,unsigned Beg, unsigned End, analogic scale){
	for(unsigned i=Beg;i<=End;i++){
		Input.push_back(In);
		Delay.push_back(i);
		InputScale.push_back(scale);
		}
	}

void RNA::Join_OutputChannel(Channel<analogic> *Out){
	if((OutputType==NONE)||(OutputType==ANALOG)){
		OutputType=ANALOG;
		}
	else{
		cout<<"RNA: Output Already define as a differet type."<<endl; 
		exit(0);
		}
	AnalogicOutput.push_back(Out);
	OutputScale.push_back(1.0);
	}

void RNA::Join_OutputChannel(Channel<analogic> *Out, analogic scale){
	if((OutputType==NONE)||(OutputType==ANALOG)){
		OutputType=ANALOG;
		}
	else{
		cout<<"RNA: Output Already define as a differet type."<<endl; 
		exit(0);
		}
	AnalogicOutput.push_back(Out);
	OutputScale.push_back(scale);
	}

void RNA::Join_OutputChannel(Channel<digital> *Out){
	if((OutputType==NONE)||(OutputType==DIGITAL)){
		OutputType=DIGITAL;
		}
	else{
		cout<<"RNA: Output Already define as a differet type."<<endl; 
		exit(0);
		}
	fann_set_activation_function_output(ann, FANN_THRESHOLD);
	DigitalOutput.push_back(Out);
	OutputScale.push_back(1.0);
	}

bool RNA::Prepare(float freq){
	cout<<"Prepare RNA"<<endl;
	if(!ann){
		cout<<"RNA:ANN Not Allocated."<<endl; 
		return false;		
		}
	//fann_print_parameters(ann);
	//fann_print_connections(ann);
	input=(fann_type*)malloc(Input.size()*sizeof(fann_type));
	if(fann_get_num_input(ann)!=Input.size()){
		cout<<"RNA:Number of Inputs Incompatible with ANN("<<fann_get_num_input(ann)<<":"<<Input.size()<<")"<<endl; 
		return false;
		}
	output=(fann_type*)malloc((DigitalOutput.size()+AnalogicOutput.size())*sizeof(fann_type));
	if(fann_get_num_output(ann)!=(DigitalOutput.size()+AnalogicOutput.size())){
		cout<<"RNA:Number of Ouputs Incompatible with ANN."<<endl; 
		return false;	
		}
	return true;
	}

void RNA::Run(){
	for(unsigned i=0;i<Input.size();i++){
		input[i]=Input[i]->get_Value(Delay[i])/InputScale[i];
		//cout<<input[i]<<"\t";
		}
	//fann_scale_input(ann, input);
	output=fann_run(ann, input);	
	//fann_descale_input(ann, output);
	for(unsigned i=0;i<(DigitalOutput.size()+AnalogicOutput.size());i++){
		//cout<<"Out:"<<output[i]<<endl;
		//cout<<"Saida:"<<output[i]<<endl;
		if(OutputType==DIGITAL)
			DigitalOutput[i]->insert_Value((bool)output[i]);
		else
			AnalogicOutput[i]->insert_Value(output[i]*OutputScale[i]);
		}
	}

RNA::~RNA(){
	//fann_destroy(ann);
	//free(input);
	//free(output);
	}
#endif

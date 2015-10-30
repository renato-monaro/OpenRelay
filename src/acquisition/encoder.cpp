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
#include "encoder.h"

#ifdef COMEDI

Encoder::Encoder(const char *device, Channel<analogic>  *out, unsigned holes, int a, int b){
	FunctionalDescription_set("Hardware Acquisition"); 
	comediDev=comedi_open(device);
	if (comediDev==NULL){
		cout<<"Error Opening the Comedi Device:"<<(char*)device<<endl;
		exit(1);
		}
	A=a;
	B=b;
	Z=-1;
	subDevice=-1;
	Out=out;
	Ratio=2*M_PI/holes;
	}

Encoder::Encoder(const char *device, Channel<analogic> *out, unsigned holes, int a, int b, int z){
	FunctionalDescription_set("Hardware Acquisition"); 
	comediDev=comedi_open(device);
	if (comediDev==NULL){
		cout<<"Error Opening the Comedi Device:"<<(char*)device<<endl;
		exit(1);
		}
	A=a;
	B=b;
	Z=z;
	subDevice=-1;
	Out=out;
	Ratio=2*M_PI/holes;
	}



bool Encoder::Run()
{
	data_ant=data;
	comedi_data_read(comediDev, subDevice, 0, 0, 0, &data);
//	comedi_data_write(comediDev, subDevice, 0, 0, 0, 0);
//	Out->insert_Value(Ratio*(data-data_ant)/(MasterClock->get_Value(0)-MasterClock->get_Value(1)));
	Out->insert_Value(Ratio*data);
return true;
}

void Encoder::RefreshTime(){
	MasterClock->insert_Value(orelay_gettime());
}

bool Encoder::Prepare(float){
	lsampl_t counter_mode;

	subDevice=comedi_find_subdevice_by_type(comediDev, COMEDI_SUBD_COUNTER,0);
	if(subDevice<0){
		cout<<"Acquisition::Comedi:: Encoder SubDevice Not Found"<<endl;
	    comedi_close(comediDev);
		return false;
		}

	if(comedi_lock(comediDev, subDevice) < 0) {
		cout<< "Comedi lock failed for subdevice "<<subDevice<<endl;
		comedi_close(comediDev);
		return false;
		}

	comedi_insn insn;
	lsampl_t data[3];
	memset(&insn, 0, sizeof(comedi_insn));
	insn.insn = INSN_CONFIG;
	insn.data = data;
	insn.subdev = subDevice;
	insn.chanspec = 0;
	insn.n = 1;
	data[0] = INSN_CONFIG_RESET;

	if (comedi_do_insn(comediDev, &insn) < 0) {
		cout<< "Comedi do_insn failed on instruction "<< data[0]<<endl;
		comedi_unlock(comediDev, subDevice);
		comedi_close(comediDev);
		return false;
		}

	/* set initial counter value by writing to channel 0 */
	comedi_data_write(comediDev, subDevice, 0, 0, 0, 0);

	counter_mode = NI_GPCT_COUNTING_DIRECTION_HW_UP_DOWN_BITS|NI_GPCT_COUNTING_MODE_QUADRATURE_X4_BITS;
	//counter_mode |= (NI_GPCT_INDEX_ENABLE_BIT | NI_GPCT_INDEX_PHASE_HIGH_A_HIGH_B_BITS);

	insn.n=3;
    data[0]=INSN_CONFIG_SET_GATE_SRC;
    data[1]=0;
    data[2]=NI_GPCT_DISABLED_GATE_SELECT;
    if (comedi_do_insn(comediDev, &insn) < 0){
		cout<<"Comedi do_insn failed on instruction "<<data[0]<<endl;
		comedi_unlock(comediDev, subDevice);
		comedi_close(comediDev);
		return false;
		}
	insn.n=3;
    data[0]=INSN_CONFIG_SET_GATE_SRC;
    data[1]=1;
    data[2]=NI_GPCT_DISABLED_GATE_SELECT;
    if (comedi_do_insn(comediDev, &insn) < 0){
		cout<<"Comedi do_insn failed on instruction "<<data[0]<<endl;
		comedi_unlock(comediDev, subDevice);
		comedi_close(comediDev);
		return false;
		}
	insn.n=3;
    data[0]=INSN_CONFIG_SET_OTHER_SRC;
    data[1]=NI_GPCT_SOURCE_ENCODER_A;
    data[2]=NI_GPCT_PFI_OTHER_SELECT(A);
    if (comedi_do_insn(comediDev, &insn) < 0){
		cout<<"Comedi do_insn failed on instruction "<<data[0]<<endl;
		comedi_unlock(comediDev, subDevice);
		comedi_close(comediDev);
		return false;
		}
	insn.n=3;
    data[0]=INSN_CONFIG_SET_OTHER_SRC;
    data[1]=NI_GPCT_SOURCE_ENCODER_B;
    data[2]=NI_GPCT_PFI_OTHER_SELECT(B);
    if (comedi_do_insn(comediDev, &insn) < 0){
		cout<<"Comedi do_insn failed on instruction "<<data[0]<<endl;
		comedi_unlock(comediDev, subDevice);
		comedi_close(comediDev);
		return false;
		}
/*	insn.n=3;
    data[0]=INSN_CONFIG_SET_OTHER_SRC;
    data[1]=NI_GPCT_SOURCE_ENCODER_Z;
    data[2]=NI_GPCT_PFI_OTHER_SELECT(Z);
    if (comedi_do_insn(comediDev, &insn) < 0){
		cout<<"Comedi do_insn failed on instruction "<<data[0]<<endl;
		comedi_unlock(comediDev, subDevice);
		comedi_close(comediDev);
		return false;
		}*/
	insn.n=2;
    data[0]=INSN_CONFIG_SET_COUNTER_MODE;
    data[1]=counter_mode;
    data[2]=0;
    if (comedi_do_insn(comediDev, &insn) < 0){
		cout<<"Comedi do_insn failed on instruction "<<data[0]<<endl;
		comedi_unlock(comediDev, subDevice);
		comedi_close(comediDev);
		return false;
		}
	insn.n=2;
    data[0]=INSN_CONFIG_ARM;
    data[1]=NI_GPCT_ARM_IMMEDIATE;
    data[2]=0;
    if (comedi_do_insn(comediDev, &insn) < 0){
		cout<<"Comedi do_insn failed on instruction "<<data[0]<<endl;
		comedi_unlock(comediDev, subDevice);
		comedi_close(comediDev);
		return false;
		}
	comedi_unlock(comediDev, subDevice);
	return true;
	}
Encoder::~Encoder(){
	comedi_unlock(comediDev, subDevice);
	comedi_close(comediDev);
	}
#endif

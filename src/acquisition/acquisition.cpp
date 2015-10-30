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
#include "acquisition.h"

using namespace orelay;
#ifdef RTAI
#ifdef COMEDI
double comedi_to_phys(lsampl_t data,comedi_krange *rng,lsampl_t maxdata)
{
	double x;

	if(!rng)return NAN;
	if(!maxdata)return NAN;

	x=data;
	x/=maxdata;
	x*=(rng->max-rng->min);
	x+=rng->min;
	x/=1E6;
	return x;
}


lsampl_t comedi_from_phys(double data,comedi_krange *rng,lsampl_t maxdata)
{
	double s;

	if(!rng)return 0;
	if(!maxdata)return 0;
	data*=1E6;
	s=(data-rng->min)/(rng->max-rng->min)*maxdata;
	if(s<0)return 0;
	if(s>maxdata)return maxdata;

	return (lsampl_t)(floor(s+0.5));
}
#endif
#endif

void Acquisition::FunctionalDescription_set(string Desc){
	FunctionalDescription=Desc;
	}

string Acquisition::FunctionalDescription_get(){
	return FunctionalDescription;
	}
void Acquisition::SetClock(Channel<timer> *t){
	MasterClock=t;
	}



void Acquisition::Join_Channel(Channel<digital> *Ch,bool direction,unsigned char channel){
	if(!Ch->InUse()){
		Ch->UseIt();
		if(direction==OUTPUT){
			DigitalOut.push_back(Ch);
			DigitalOutBoardChannel.push_back(channel);
			}
		else{
			DigitalIn.push_back(Ch);
			DigitalInBoardChannel.push_back(channel);
			}
		}
	else{
		cout<<"Acquisition::Error:: Channel" <<Ch->get_Name()<<" Already in Use."<<endl;
		exit(0);
		}

	}
void Acquisition::Join_Channel(Channel<analogic> *Ch,bool direction,unsigned char channel){
	Join_Channel(Ch,direction,channel,1.0,0.0);
	}
void Acquisition::Join_Channel(Channel<analogic> *Ch,bool direction,unsigned char channel,float ratio){
	Join_Channel(Ch,direction,channel,ratio,0.0);
	}
void Acquisition::Join_Channel(Channel<analogic> *Ch,bool direction,unsigned char channel,float ratio, float offset){
	if(!Ch->InUse()){
		Ch->UseIt();
		if(direction==OUTPUT){
			AnalogicOut.push_back(Ch);
			AnalogicOutBoardChannel.push_back(channel);
			AnalogicOutRatio.push_back(ratio);
			AnalogicOutOffset.push_back(offset);
			}
		else{
			AnalogicIn.push_back(Ch);
			AnalogicInBoardChannel.push_back(channel);
			AnalogicInRatio.push_back(ratio);
			AnalogicInOffset.push_back(offset);
			}
		}
	else{
		cout<<"Acquisition::Error:: Channel" <<Ch->get_Name()<<" Already in Use."<<endl;
		exit(0);
		}
	}



#ifdef COMEDI
HardwareAcquisition::HardwareAcquisition(const char *device){
	FunctionalDescription_set("Hardware Acquisition"); 
	range=COMEDI_RANGE;
	aref=COMEDI_REF;
	comediDev=comedi_open(device);
	if (comediDev==NULL){
		cout<<"Error Opening the Comedi Device:"<<(char*)device<<endl;
		exit(1);
		}
	}

bool HardwareAcquisition::Prepare(float Samp){
	int total=0;
//	comedi_calibration_t *calibration;

//TODO: Preparar os Ranges como? Convert to float???
/*
int  comedi_find_range  (comedi_t  *  device,  unsigned  int subdevice,
       unsigned int channel, unsigned int unit, double min, double max);
*/

if(AnalogicIn.size()>0){
	data_analog_in=(lsampl_t*)malloc(AnalogicIn.size()*sizeof(lsampl_t));
	if(data_analog_in==NULL){
		cout<<"Acquisition::Error in Memory Alocation"<<endl;
		return false;
		}
	subdevice_analog_in=comedi_find_subdevice_by_type(comediDev, COMEDI_SUBD_AI,0);
	if(subdevice_analog_in<0){
		cout<<"Acquisition::Comedi:: Analogic Input SubDevice Not Found"<<endl;
		return false;
		}
	if(comedi_lock(comediDev,subdevice_analog_in)<0){
		cout<<"Acquisition::Comedi:: Error Locking the Analogic Input SubDevice"<<endl;
		comedi_unlock(comediDev,subdevice_analog_in);
		return false;
		}
	comedi_unlock (comediDev,subdevice_analog_in);
	if(comedi_get_n_channels (comediDev,subdevice_analog_in)<(int)AnalogicIn.size()){
		cout<<"Acquisition::Comedi:: Insufficient Number of Channel(s) Available in Analogic Input SubDevice"<<endl;
		return false;
		}
	for(unsigned int i=0;i<AnalogicIn.size();i++){
		if(comedi_get_n_channels (comediDev,subdevice_analog_in)<AnalogicInBoardChannel[i]){
			cout<<"Acquisition::Comedi:: Channel "<<(int)AnalogicInBoardChannel[i]<<" Not Available in Analogic Input SubDevice"<<endl;
			return false;
			}
		for(unsigned int j=0;j<i;j++){
			if(AnalogicInBoardChannel[i]==AnalogicInBoardChannel[j]){
				cout<<"Acquisition::Comedi:: Channel "<<(int)AnalogicInBoardChannel[i]<<" Already Taken in Analogic Input SubDevice"<<endl;
				return false;
				}
			}
		}
	}
if(AnalogicOut.size()>0){
	data_analog_out=(lsampl_t*)malloc(AnalogicOut.size()*sizeof(lsampl_t));
	if(data_analog_out==NULL){
		cout<<"Acquisition::Error in Memory Alocation"<<endl;
		return false;
		}
	subdevice_analog_out=comedi_find_subdevice_by_type(comediDev, COMEDI_SUBD_AO,0);
	if(subdevice_analog_out<0){
		cout<<"Acquisition::Comedi:: Analogic Output SubDevice Not Found"<<endl;
		return false;
		}
	if(comedi_lock(comediDev,subdevice_analog_out)<0){
		cout<<"Acquisition::Comedi:: Error Locking the Analogic Output SubDevice"<<endl;
		comedi_unlock(comediDev,subdevice_analog_out);
		return false;
		}
	comedi_unlock (comediDev,subdevice_analog_out);
	if(comedi_get_n_channels (comediDev,subdevice_analog_out)<(int)AnalogicOut.size())
		{
		cout<<"Acquisition::Comedi:: Insufficient Number of Channel(s) Available in Analogic Output SubDevice"<<endl;
		return false;
		}
	for(unsigned int i=0;i<AnalogicOut.size();i++){
		if(comedi_get_n_channels (comediDev,subdevice_analog_out)<AnalogicOutBoardChannel[i]){
			cout<<"Acquisition::Comedi:: Channel "<<(int)AnalogicOutBoardChannel[i]<<" Not Available in Analogic Output SubDevice"<<endl;
			return false;
			}
		for(unsigned int j=0;j<i;j++){
			if(AnalogicOutBoardChannel[i]==AnalogicOutBoardChannel[j]){
				cout<<"Acquisition::Comedi:: Channel "<<(int)AnalogicOutBoardChannel[i]<<" Already Taken in Analogic Output SubDevice"<<endl;
				return false;				
				}
			}
		}
	}
if(DigitalOut.size()>0){
	data_digital_out=(lsampl_t*)malloc(DigitalOut.size()*sizeof(lsampl_t));
	if(data_digital_out==NULL){
		cout<<"Acquisition::Error in Memory Alocation"<<endl;
		return false;
		}
	subdevice_digital_out=comedi_find_subdevice_by_type(comediDev, COMEDI_SUBD_DO,0);
	if(subdevice_digital_out<0){
		cout<<"Acquisition::Comedi:: Digital Output SubDevice Not Found"<<endl;
		cout<<"Trying the Digital Input/Output SubDevice"<<endl;
		subdevice_digital_out=comedi_find_subdevice_by_type(comediDev, COMEDI_SUBD_DIO,0);
		if(subdevice_digital_out<0){
			cout<<"Acquisition::Comedi:: Digital Input/Output SubDevice Not Found"<<endl;
			return false;			
			}
		for(unsigned int i=0;i<DigitalOut.size();i++){
			 if(comedi_dio_config(comediDev,subdevice_digital_out,DigitalOutBoardChannel[i],COMEDI_OUTPUT)<0){
				cout<<"Acquisition::Comedi:: Error Configuring Direction of Digital Input/Output SubDevice"<<endl;
				return false;			
				}			
			}
		}
	if(comedi_lock(comediDev,subdevice_digital_out)<0){
		cout<<"Acquisition::Comedi:: Error Locking the Digital Output SubDevice"<<endl;
		comedi_unlock(comediDev,subdevice_digital_out);
		return false;
		}
	comedi_unlock (comediDev,subdevice_digital_out);
	if(comedi_get_n_channels (comediDev,subdevice_digital_out)<(int)DigitalOut.size()){
		cout<<"Acquisition::Comedi:: Insufficient Number of Channel(s) Available in Digital Output SubDevice"<<endl;
		return false;
		}
	for(unsigned int i=0;i<DigitalOut.size();i++){
		if(comedi_get_n_channels (comediDev,subdevice_digital_out)<DigitalOutBoardChannel[i]){
			cout<<"Acquisition::Comedi:: Channel "<<(int)DigitalOutBoardChannel[i]<<" Not Available in Digital Output SubDevice"<<endl;
			return false;
			}
	for(unsigned int j=0;j<i;j++){
			if(DigitalOutBoardChannel[i]==DigitalOutBoardChannel[j]){
				cout<<"Acquisition::Comedi:: Channel "<<(int)DigitalOutBoardChannel[i]<<" Already Taken in Digital Output SubDevice"<<endl;
				return false;				
				}
			}
		}
	}
if(DigitalIn.size()>0){
	data_digital_in=(lsampl_t*)malloc(DigitalIn.size()*sizeof(lsampl_t));
	total=DigitalIn.size()+DigitalOut.size()+AnalogicIn.size()+AnalogicOut.size();
	if(data_digital_in==NULL){
		cout<<"Acquisition::Error in Memory Alocation"<<endl;
		return false;
		}
	subdevice_digital_in=comedi_find_subdevice_by_type(comediDev, COMEDI_SUBD_DI,0);
	if(subdevice_digital_in<0){
		cout<<"Acquisition::Comedi:: Digital Input SubDevice Not Found"<<endl;
		cout<<"Trying the Digital Input/Output SubDevice"<<endl;
		subdevice_digital_in=comedi_find_subdevice_by_type(comediDev, COMEDI_SUBD_DIO,0);
		if(subdevice_digital_in<0){
			cout<<"Acquisition::Comedi:: Digital Input/Output SubDevice Not Found"<<endl;
			return false;			
			}
		for(unsigned int i=0;i<DigitalIn.size();i++){
			for(unsigned int j=0;j<DigitalOut.size();j++){
				if(DigitalInBoardChannel[i]==DigitalOutBoardChannel[j]){
					cout<<"Acquisition::Comedi:: Channel "<<(int)DigitalInBoardChannel[i]<<" Already Taken in Digital Input/Output SubDevice"<<endl;
					return false;				
					}
			}
			 if(comedi_dio_config(comediDev,subdevice_digital_in,DigitalInBoardChannel[i],COMEDI_INPUT)<0){
				cout<<"Acquisition::Comedi:: Error Configuring Direction of Digital Input/Output SubDevice"<<endl;
				return false;			
				}			
			}
		}
	if(comedi_lock(comediDev,subdevice_digital_in)<0){
		cout<<"Acquisition::Comedi:: Error Locking the Digital Input SubDevice"<<endl;
		comedi_unlock(comediDev,subdevice_digital_in);
		return false;
		}
	comedi_unlock (comediDev,subdevice_digital_in);
	if(comedi_get_n_channels (comediDev,subdevice_digital_in)<(int)DigitalIn.size()){
		cout<<"Acquisition::Comedi:: Insufficient Number of Channel(s) Available in Digital Input SubDevice"<<endl;
		return false;
		}
	for(unsigned int i=0;i<DigitalIn.size();i++){
		if(comedi_get_n_channels (comediDev,subdevice_digital_in)<DigitalInBoardChannel[i]){
			cout<<"Acquisition::Comedi:: Channel "<<(int)DigitalInBoardChannel[i]<<" Not Available in Digital Input SubDevice"<<endl;
			return false;
			}
	for(unsigned int j=0;j<i;j++){
			if(DigitalInBoardChannel[i]==DigitalInBoardChannel[j]){
				cout<<"Acquisition::Comedi:: Channel "<<(int)DigitalInBoardChannel[i]<<" Already Taken in Digital Input SubDevice"<<endl;
				return false;				
				}
			}
		}
	}

/*t=(timer*)malloc(Timer.size()*sizeof(timer));
if(t==NULL){
	cout<<"Acquisition::Error in Memory Alocation"<<endl;
	return false;
	}*/

#ifndef RTAI
analog_in_maxdata=(lsampl_t*)malloc(AnalogicIn.size()*sizeof(lsampl_t));
for(unsigned i=0;i<AnalogicIn.size();i++)
	analog_in_maxdata[i]=comedi_get_maxdata(comediDev,subdevice_analog_in,AnalogicInBoardChannel[i]);
analog_in_range=(comedi_range**)malloc(AnalogicIn.size()*sizeof(comedi_range*));
for(unsigned i=0;i<AnalogicIn.size();i++)
	analog_in_range[i]=comedi_get_range(comediDev,subdevice_analog_in,AnalogicInBoardChannel[i],range);


analog_out_maxdata=(lsampl_t*)malloc(AnalogicOut.size()*sizeof(lsampl_t));
for(unsigned i=0;i<AnalogicOut.size();i++)
	analog_out_maxdata[i]=comedi_get_maxdata(comediDev,subdevice_analog_out,AnalogicOutBoardChannel[i]);
analog_out_range=(comedi_range**)malloc(AnalogicOut.size()*sizeof(comedi_range*));
for(unsigned i=0;i<AnalogicOut.size();i++)
	analog_out_range[i]=comedi_get_range(comediDev,subdevice_analog_out,AnalogicOutBoardChannel[i],range);
#endif


#ifdef RTAI
analog_in_maxdata=(lsampl_t*)malloc(AnalogicIn.size()*sizeof(lsampl_t));
for(unsigned i=0;i<AnalogicIn.size();i++)
	analog_in_maxdata[i]=comedi_get_maxdata(comediDev,subdevice_analog_in,AnalogicInBoardChannel[i]);
analog_in_range=(comedi_krange*)malloc(AnalogicIn.size()*sizeof(comedi_krange));
for(unsigned i=0;i<AnalogicIn.size();i++){
	comedi_get_krange(comediDev,subdevice_analog_in,AnalogicInBoardChannel[i],range,&analog_in_range[i]);
	}


analog_out_maxdata=(lsampl_t*)malloc(AnalogicOut.size()*sizeof(lsampl_t));
for(unsigned i=0;i<AnalogicOut.size();i++)
	analog_out_maxdata[i]=comedi_get_maxdata(comediDev,subdevice_analog_out,AnalogicOutBoardChannel[i]);
analog_out_range=(comedi_krange*)malloc(AnalogicOut.size()*sizeof(comedi_krange));
for(unsigned i=0;i<AnalogicOut.size();i++){
	comedi_get_krange(comediDev,subdevice_analog_out,AnalogicOutBoardChannel[i],range,&analog_out_range[i]);
	}
#endif

/*for(unsigned i=0;i<AnalogicIn.size();i++)
	cout<<"MaxRange:"<<analog_in_maxdata[i]<<" Max:"<<analog_in_range[i].max<<" Min:"<<analog_in_range[i].min<<endl;
for(unsigned i=0;i<AnalogicOut.size();i++)
	cout<<"MaxRange:"<<analog_out_maxdata[i]<<" Max:"<<analog_out_range[i].max<<" Min:"<<analog_out_range[i].min<<endl;*/
/*
calibration_path=comedi_get_default_calibration_path(comediDev);

AnalogicInCalibration=(comedi_polynomial_t*)malloc(AnalogicIn.size()*sizeof(comedi_polynomial_t));
if(comedi_get_subdevice_flags (comediDev,subdevice_analog_in)&SDF_SOFT_CALIBRATED){
	calibration=comedi_parse_calibration_file(calibration_path);
	if(calibration){
		for(unsigned int i=0;i<AnalogicIn.size();i++){
			if(comedi_get_softcal_converter(subdevice_analog_in,AnalogicIn[i]->get_BoardChannel(),range,COMEDI_TO_PHYSICAL, calibration, &AnalogicInCalibration[i])<0){
				cout<<"Acquisition::Comedi:: Soft Calibration Not Avaible for Channel "<<(int) AnalogicIn[i]->get_BoardChannel() <<" In Analogic InPut Subdevice"<<endl;
				comedi_cleanup_calibration(calibration);
				free(calibration_path);
				return true;
				}
			}
		comedi_cleanup_calibration(calibration);
		}
	else{
		cout<<"Acquisition::Comedi:: Unable to get Calibration"<<endl;
		comedi_cleanup_calibration(calibration);
		free(calibration_path);
		return true;
		}		
	}
else{
	for(unsigned int i=0;i<AnalogicIn.size();i++){
		if(comedi_get_hardcal_converter(comediDev,subdevice_analog_in,AnalogicIn[i]->get_BoardChannel(),range,COMEDI_TO_PHYSICAL, &AnalogicInCalibration[i])<0){
			cout<<"Acquisition::Comedi:: Hard Calibration Not Avaible for Channel "<<(int)AnalogicIn[i]->get_BoardChannel() <<" In Analogic InPut Subdevice"<<endl;
			return true;
			}
		}
	}

AnalogicOutCalibration=(comedi_polynomial_t*)malloc(AnalogicOut.size()*sizeof(comedi_polynomial_t));
if(comedi_get_subdevice_flags (comediDev,subdevice_analog_out)&SDF_SOFT_CALIBRATED){
	calibration=comedi_parse_calibration_file(calibration_path);
	if(calibration){
		for(unsigned int i=0;i<AnalogicOut.size();i++){
			if(comedi_get_softcal_converter(subdevice_analog_out,AnalogicOut[i]->get_BoardChannel(),range,COMEDI_FROM_PHYSICAL, calibration, &AnalogicOutCalibration[i])<0){
				cout<<"Acquisition::Comedi:: Soft Calibration Not Avaible for Channel "<<(int)AnalogicOut[i]->get_BoardChannel() <<" In Analogic OutPut Subdevice"<<endl;
				comedi_cleanup_calibration(calibration);
				free(calibration_path);
				return true;
				}
			}
		comedi_cleanup_calibration(calibration);
		}
	else{
		cout<<"Acquisition::Comedi:: Unable to get Calibration"<<endl;
		comedi_cleanup_calibration(calibration);
		free(calibration_path);
		return true;
		}
	}
else{
	for(unsigned int i=0;i<AnalogicOut.size();i++){
		if(comedi_get_hardcal_converter(comediDev,subdevice_analog_out,AnalogicOut[i]->get_BoardChannel(),range,COMEDI_FROM_PHYSICAL, &AnalogicInCalibration[i])<0){
			cout<<"Acquisition::Comedi:: Hard Calibration Not Avaible for Channel "<<(int)AnalogicOut[i]->get_BoardChannel() <<" In Analogic OutPut Subdevice"<<endl;
			return true;
			}
		}
	}
free(calibration_path);
*/
return true; 
}

bool HardwareAcquisition::Run()
{
	#ifdef RTAI
	//	MasterClock->insert_Value(rt_get_time_ns());		
	#endif
//MasterClock->insert_Value(boost::posix_time::microsec_clock::local_time());
/*	if(comedi_do_insnlist(comediDev,&instruction_list)<0){
		cout<<"Error Executing the Acquisition Instruction"<<endl;
		}*/
		#ifdef DEBUG
			cout<<"Hardware Acquisition\n" ;
		#endif
/*	for(unsigned int i=0;i<Timer.size();i++)
		Timer[i]->insert_Value(boost::posix_time::microsec_clock::local_time());*/
	lsampl_t data;
	unsigned int bit;
	for(unsigned int i=0;i<AnalogicIn.size();i++){

		comedi_data_read_delayed(comediDev,subdevice_analog_in,AnalogicInBoardChannel[i],range,aref,&data,1000);

		//AnalogicIn[i]->insert_Value(AnalogicIn[i]->get_Ratio()*comedi_to_physical(data_analog_in[i],&AnalogicInCalibration[i]));
	#ifdef RTAI
		AnalogicIn[i]->insert_Value(AnalogicInRatio[i]*comedi_to_phys(data,&analog_in_range[i],analog_in_maxdata[i])+AnalogicInOffset[i]);
	#endif
	#ifndef RTAI
		AnalogicIn[i]->insert_Value(AnalogicInRatio[i]*comedi_to_phys(data,analog_in_range[i],analog_in_maxdata[i])+AnalogicInOffset[i]);
	#endif
		//AnalogicIn[i]->insert_Value(AnalogicIn[i]->get_Ratio()*comedi_to_physical(data,&AnalogicInCalibration[i]));
		}
	for(unsigned int i=0;i<AnalogicOut.size();i++){
//		data_analog_out[i]=comedi_from_physical (AnalogicOut[i]->get_Value(),&AnalogicOutCalibration[i])/AnalogicOut[i]->get_Ratio();
//		data=comedi_from_physical (AnalogicOut[i]->get_Value()/AnalogicOut[i]->get_Ratio(),&AnalogicOutCalibration[i]);
	#ifdef RTAI
		data=comedi_from_phys (AnalogicOut[i]->get_Value()/AnalogicOutRatio[i]+AnalogicOutOffset[i],&analog_out_range[i],analog_out_maxdata[i]);
	#endif
	#ifndef RTAI
		data=comedi_from_phys (AnalogicOut[i]->get_Value()/AnalogicOutRatio[i]+AnalogicOutOffset[i],analog_out_range[i],analog_out_maxdata[i]);
	#endif
		comedi_data_write(comediDev,subdevice_analog_out,AnalogicOutBoardChannel[i],range,aref,data);
		}
	for(unsigned int i=0;i<DigitalIn.size();i++){
		comedi_dio_read(comediDev,subdevice_digital_in,DigitalInBoardChannel[i],&bit);
		DigitalIn[i]->insert_Value(bit);
		}
	for(unsigned int i=0;i<DigitalOut.size();i++){
		bit=DigitalOut[i]->get_Value();
		comedi_dio_write(comediDev,subdevice_digital_out,DigitalOutBoardChannel[i],bit);
		}
	return true;
}


void HardwareAcquisition::RefreshTime(){
	MasterClock->insert_Value(orelay_gettime());
}	

#endif
/*FileAcquisition::FileAcquisition(const char *filename, unsigned Bits){ 
	FunctionalDescription_set("Software Acquisition"); 
	File.open(filename);
	if (!File.is_open()){
		cout<<"Error Opening the Input File:"<<(char*)filename<<endl;
		exit(1);
		}
	EventDateTime=orelay_gettime();
	DataCounter=0;
	DACBits=Bits;
	}*/

FileAcquisition::FileAcquisition(const char *filename){ 
	FunctionalDescription_set("Software Acquisition"); 
	File.open(filename);
	if (!File.is_open()){
		cout<<"Error Opening the Input File:"<<(char*)filename<<endl;
		exit(1);
		}
	EventDateTime=orelay_gettime();
	DataCounter=0;
	//DACBits=0;
	FileSamplingRate=0;
	}

FileAcquisition::FileAcquisition(const char *filename, float Samp){ 
	FunctionalDescription_set("Software Acquisition"); 
	File.open(filename);
	if (!File.is_open()){
		cout<<"Error Opening the Input File:"<<(char*)filename<<endl;
		exit(1);
		}
	EventDateTime=orelay_gettime();
	DataCounter=0;
	FileSamplingRate=Samp;
	//DACBits=0;
	}
// compute greatest common divisor
// algorithm found on wikipedia.

int gcd (int a, int b) {

	int f,g;
	
	if ( a == 0 || b == 0) 
		return 1;

	g = b;
	f  =  a % b;
	
	while (f != 0) {
		b = f;
		f = g % f;
		g = b;
	}

	return g;
}


bool FileAcquisition::Prepare(float Samp){
	SamplingRate=Samp;
	string Linha;
	unsigned max_column=0;
	Columns=0;
	double tmp;
	if(!File.eof()){
		getline(File,Linha); /*ignore Header*/
		getline(File,Linha);
		getline(File,Linha);
		getline(File,Linha);
		std::stringstream is(Linha);
		while(is>>tmp)
				Columns++;
		#ifdef DEBUG_ALL
			cout<<"File Columns: "<<Columns<<endl;
			for(unsigned k=0;k<Columns;k++)
				cout<<getValueColumnN(Linha,k)<<"\t";
			cout<<endl;
		#endif
		}
	else{
		cout<<"Acquisition::File: Empty InputFile"<<endl;
		return false;
		}
/*	if(Columns<(AnalogicIn.size()+DigitalIn.size())){
		cout<<"Acquisition::File: Insufficient Number of Columns("<<Columns<<") in the InputFile, required ("<<(AnalogicIn.size()+DigitalIn.size())<<")"<<endl;
		return false;
		}*/
	for(unsigned int i=0;i<AnalogicIn.size();i++)
		if(AnalogicInBoardChannel[i]>max_column)
			max_column=AnalogicInBoardChannel[i];
	for(unsigned int i=0;i<DigitalIn.size();i++)
		if(DigitalInBoardChannel[i]>max_column)
			max_column=DigitalInBoardChannel[i];
	max_column+=1;
	if(max_column>=Columns){
		cout<<"Acquisition::File: Column:"<< max_column<<" not Avaible in the InputFile"<<endl;
		return false;
		}
	vector<vector<double> > OriginalData;
	OriginalData.resize(Columns);
	vector<vector<double> > UpSampled;
	UpSampled.resize(Columns);
	Data.resize(Columns);
	while(!File.eof()){
//	while((!File.eof())&&(getValueColumnN(Linha,0)>=0)){
		OriginalData[0].push_back(getValueColumnN(Linha,0));
		for(unsigned int i=1;i<Columns;i++)
			OriginalData[i].push_back(getValueColumnN(Linha,i));
		/*for(unsigned int i=0;i<DigitalIn.size();i++)
			OriginalData[AnalogicIn.size()+i+1].push_back(getValueColumnN(Linha,DigitalInBoardChannel[i]));*/
		getline(File,Linha);
		}
	File.close();
	butterworth<double> BPF(0.5,2,3.0);
	if(FileSamplingRate<=0)
		FileSamplingRate=1/(OriginalData[0][1]-OriginalData[0][0]);
	cout<<"FileSamplingRate:"<<FileSamplingRate<<"("<<(int)ceil(FileSamplingRate)<<")"<<endl;
	cout<<"SamplingRate:"<<SamplingRate<<"("<<(int)ceil(SamplingRate)<<")"<<endl;
	int MDC=gcd((int)ceil(SamplingRate), (int)ceil(FileSamplingRate));
	unsigned M=(unsigned)(ceil(FileSamplingRate)/MDC);
	unsigned L=(unsigned)(SamplingRate/MDC);

	cout<<"L:"<<L<<" M:"<<M<<" MDC:"<<MDC<<endl;

//Upsampling 
cout<<"Upsampling"<<endl;
	for(unsigned i=0;i<OriginalData[0].size()-1;i++){
		for(unsigned j=0;j<OriginalData.size();j++){
			UpSampled[j].push_back(OriginalData[j][i]);
			for(unsigned k=1;k<L;k++)
				UpSampled[j].push_back(BPF.clock(OriginalData[j][i]+k*(OriginalData[j][i+1]-OriginalData[j][i])/L));
			}
		}
cout<<"Downsampling"<<endl;
//DownSampling
	for(unsigned j=0;j<UpSampled.size();j++){
		BPF.reset();		
		for(unsigned i=0;i<UpSampled[0].size();i=i+M)
			Data[j].push_back(UpSampled[j][i]);
	}
	
	Data[0].clear();
	for(unsigned i=0;i<UpSampled[0].size()/M;i++){
		Data[0].push_back(double(i)/Samp);
		//Data[0].push_back(UpSampled[0][i]);
	}


	/*for(unsigned j=0;j<Data[0].size();j++)
		for(unsigned int i=0;i<Columns;i++)
			{
			if(AnalogicIn[i]->is_Saturable()&&(DACBits>0)){
				double Max=AnalogicIn[i]->get_UpperLimit();
				double Min=AnalogicIn[i]->get_LowerLimit();
				unsigned ContMax=(1<<DACBits)-1;
				double DigitalValue;
				double Q=(Max-Min)/(double)ContMax;
				for(unsigned k=0;k<ContMax;k++){
					DigitalValue=Min+k*Q;
					if(Data[i+1][j]<=(DigitalValue+Q/2))
						break;
					}
				Data[i+1][j]=DigitalValue;
				}
			else
			Data[i+1][j]=Data[i+1][j];
			}*/
	cout<<"Saving"<<endl;
	ofstream IMP("/tmp/resposta.dat");
  for (unsigned i=0;i<Data[0].size();i++) {
	for(unsigned j=0;j<Data.size();j++){
		IMP<<Data[j][i]<<"\t";
		}
	IMP<< "\n";
  }
  IMP.close();
	
	return true;
}

double FileAcquisition::getValueColumnN(string Line, unsigned Index){
	double temp;
	stringstream is(Line);
	for(unsigned i=0;i<=Index;i++)
		is>>temp;
	return temp;
}
bool FileAcquisition::Run(){
	if(Data[0].size()>DataCounter){
		#ifdef DEBUG_ALL
			cout<<"File Acquisition\n" ;
		#endif
//		MasterClock->insert_Value(1E9*Data[0][DataCounter]);
		for(unsigned int i=0;i<AnalogicIn.size();i++){
			AnalogicIn[i]->insert_Value(AnalogicInOffset[i]+AnalogicInRatio[i]*Data[AnalogicInBoardChannel[i]+1][DataCounter]);
			}
		for(unsigned int i=0;i<DigitalIn.size();i++){
			DigitalIn[i]->insert_Value(Data[DigitalInBoardChannel[i]+1][DataCounter]>0.5?true:false);
			}
		DataCounter++;
		//cout<<DataCounter;
		return true;
		}
	else{
		#ifdef DEBUG_ALL
		cout<<"End Of File\n";
		#endif
		return false;
		}	
}

void FileAcquisition::RefreshTime(){
	MasterClock->insert_Value(1E9*Data[0][DataCounter]);
}	






/*
void Acquisition::SetBufferSize(unsigned Size){
	BufferSize=Size;
	MasterClock->set_Size(BufferSize);
	for(unsigned i=0;i<AnalogicOut.size();i++)
		AnalogicOut[i]->set_Size(BufferSize);
	for(unsigned i=0;i<AnalogicIn.size();i++)
		AnalogicIn[i]->set_Size(BufferSize);
	for(unsigned i=0;i<AnalogicInternal.size();i++)
		AnalogicInternal[i]->set_Size(BufferSize);
	for(unsigned i=0;i<Phasors.size();i++)
		Phasors[i]->set_Size(BufferSize);
	for(unsigned i=0;i<DigitalOut.size();i++)
		DigitalOut[i]->set_Size(BufferSize);
	for(unsigned i=0;i<DigitalIn.size();i++)
		DigitalIn[i]->set_Size(BufferSize);
	for(unsigned i=0;i<DigitalInternal.size();i++)
		DigitalInternal[i]->set_Size(BufferSize);
	for(unsigned i=0;i<Timer.size();i++)
		Timer[i]->set_Size(BufferSize);
	for(unsigned i=0;i<Strings.size();i++)
		Strings[i]->set_Size(BufferSize);
}
void Acquisition::SetSampling(float Sampling){
	SamplingRate=Sampling;
	MasterClock->set_Sampling(SamplingRate);
	for(unsigned i=0;i<AnalogicOut.size();i++)
		AnalogicOut[i]->set_Sampling(SamplingRate);
	for(unsigned i=0;i<AnalogicIn.size();i++)
		AnalogicIn[i]->set_Sampling(SamplingRate);
	for(unsigned i=0;i<AnalogicInternal.size();i++)
		AnalogicInternal[i]->set_Sampling(SamplingRate);
	for(unsigned i=0;i<Phasors.size();i++)
		Phasors[i]->set_Sampling(SamplingRate);
	for(unsigned i=0;i<DigitalOut.size();i++)
		DigitalOut[i]->set_Sampling(SamplingRate);
	for(unsigned i=0;i<DigitalIn.size();i++)
		DigitalIn[i]->set_Sampling(SamplingRate);
	for(unsigned i=0;i<DigitalInternal.size();i++)
		DigitalInternal[i]->set_Sampling(SamplingRate);
	for(unsigned i=0;i<Timer.size();i++)
		Timer[i]->set_Sampling(SamplingRate);
	for(unsigned i=0;i<Strings.size();i++)
		Strings[i]->set_Sampling(SamplingRate);
}
void Acquisition::Join_Channel( Channel<timer> *Ch)
{
	Timer.push_back(Ch);
}

void Acquisition::Join_Channel( Channel<Complex> *Ch)
{
Phasors.push_back(Ch);
}

void Acquisition::Join_Channel( Channel<string> *Ch)
{
Strings.push_back(Ch);
}

void Acquisition::Join_Channel( Channel<analogic> *Ch)
{
switch(Ch->get_Direction())
	{
	case OUTPUT:
		AnalogicOut.push_back(Ch);
		break;		
	case INPUT:
		AnalogicIn.push_back(Ch);
		break;
	case INTERNAL:
		AnalogicInternal.push_back(Ch);
		break;
	default:
		AnalogicInternal.push_back(Ch);
	}
}

void Acquisition::Join_Channel( Channel<digital> *Ch)
{
switch(Ch->get_Direction())
	{
	case OUTPUT:
		DigitalOut.push_back(Ch);
		break;		
	case INPUT:
		DigitalIn.push_back(Ch);
		break;
	case INTERNAL:
		DigitalInternal.push_back(Ch);
		break;
	default:
		DigitalInternal.push_back(Ch);
	}
}


*/

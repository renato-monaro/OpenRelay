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
#include "relay.h"
using namespace orelay;
Relay::Relay(unsigned Size, float Sampling, float Freq){
	MasterClock= new Channel<timer>("Relay Master Clock",Size);
	BufferSize=Size;
	SamplingRate=Sampling;
	Running=false;
	PerformanceEvaluation=false;
	SystemFrequency=Freq;
	#ifdef RTAI
		Ov=0;
	#endif
}
/*Relay::Relay(Channel<digital> *Trig,unsigned BeforeLoops, unsigned AfterLoops,string N){
	MasterClock= new Channel<timer>("Relay Master Clock");
	Running=false;
	PerformanceEvaluation=true;
	AfterPoints=BeforeLoops;
	BeforPoints=AfterLoops;
	Trigger=Trig;
	Trigged=0;
	Name=N;
	NOsc=0;
	}*/
void Relay::Join(Acquisition *A){
	AcqList.push_back(A);
	A->SetClock(MasterClock);
	//A->SetBufferSize(BufferSize);
	//A->SetSampling(SamplingRate);

	}
void Relay::Join(Protection *A){
	ProtectionList.push_back(A);
	A->SetClock(MasterClock);
	}
void Relay::Join(Measure *A){
	MeasureList.push_back(A);
	A->SetClock(MasterClock);
	}
void Relay::Join(Logic *A){
	LogicList.push_back(A);
	A->SetClock(MasterClock);
	}
void Relay::Join(Control *A){
	ControlList.push_back(A);
	A->SetClock(MasterClock);
	}
void Relay::Join(Oscillography *A){
	OscillographyList.push_back(A);
	A->SetClock(MasterClock);
	A->SetControl(&Running);
	}

void Relay::Join(Evaluation *A){
	PerformanceEvaluation=true;
	EvaluationList.push_back(A);
	A->SetClock(MasterClock);
	A->SetControl(&Running);
	}
#ifdef NCURSES
void Relay::Join(Display *A){
	DisplayList.push_back(A);
	A->SetClock(MasterClock);
	}
#endif
void Relay::List(){
	cout<<"Acquisition Functions:"<<AcqList.size()<<endl;
	cout<<"Protection Functions:"<<ProtectionList.size()<<endl;
	cout<<"Measure Functions:"<<MeasureList.size()<<endl;
	cout<<"Logic Functions:"<<LogicList.size()<<endl;
	cout<<"Control Functions:"<<ControlList.size()<<endl;
	cout<<"Oscillography Functions:"<<OscillographyList.size()<<endl;	
	#ifdef NCURSES
	cout<<"Display Functions:"<<DisplayList.size()<<endl;
	#endif
	}
void Relay::Execute(){

	Start();
	if(!Prepare())
		exit(0);
	while(Running)
		Step();
	Wait();
	#ifdef RTAI
		rt_make_soft_real_time();
		stop_rt_timer();
		rt_thread_delete(tarefa_principal);
	#endif
}
void Relay::PrepareEvaluation(unsigned Size){
	ExecutionOverrun=new Channel<digital>("OverRun",Size);
	ExecutionTimeList.push_back(new Channel<analogic>("Cicle",Size));
	for(unsigned i=0;i<AcqList.size();i++)
		ExecutionTimeList.push_back(new Channel<analogic>(AcqList[i]->FunctionalDescription_get(),Size));
	for(unsigned i=0;i<MeasureList.size();i++)
		ExecutionTimeList.push_back(new Channel<analogic>(MeasureList[i]->FunctionalDescription_get(),Size));
	for(unsigned i=0;i<ProtectionList.size();i++)
		ExecutionTimeList.push_back(new Channel<analogic>(ProtectionList[i]->FunctionalDescription_get(),Size));
	for(unsigned i=0;i<ControlList.size();i++)
		ExecutionTimeList.push_back(new Channel<analogic>(ControlList[i]->FunctionalDescription_get(),Size));
	for(unsigned i=0;i<LogicList.size();i++)
		ExecutionTimeList.push_back(new Channel<analogic>(LogicList[i]->FunctionalDescription_get(),Size));
	for(unsigned i=0;i<OscillographyList.size();i++)
		ExecutionTimeList.push_back(new Channel<analogic>(OscillographyList[i]->FunctionalDescription_get(),Size));
	ExecutionOverrun->set_Size(Size);
	for(unsigned j=0;j<ExecutionTimeList.size();j++)
		ExecutionTimeList[j]->set_Size(Size);
	}

bool Relay::Wait(){
	cout<<"Waiting threads"<<endl; 
	EvalThreadList.join_all();
	OscThreadList.join_all();
	return true;
	}

bool Relay::Prepare(){
	if(PerformanceEvaluation){
		PrepareEvaluation(BufferSize);
		for(unsigned i=0;i<EvaluationList.size();i++){
				EvaluationList[i]->Join_Channel(ExecutionOverrun);
				for(unsigned j=0;j<ExecutionTimeList.size();j++)
					EvaluationList[i]->Join_Channel(ExecutionTimeList[j]);
				if(!EvaluationList[i]->Prepare(SamplingRate,BufferSize,SystemFrequency)){
					cout<<"Error Preparing the Evaluation #"<<i<<"."<<endl;
					return false;
					}
			EvalThreadList.create_thread(boost::bind(&Evaluation::Poll, EvaluationList[i]));
			}
		}	
	if(AcqList.size()==0){
		cout<<"You must create at least one Acquisiton method"<<endl;
		return false;
		}
	for(unsigned i=0;i<AcqList.size();i++){
		if(!AcqList[i]->Prepare(SamplingRate)){
			cout<<"Error Preparing the Acquisition #"<<i<<"."<<endl;
			return false;
			}
		}
	for(unsigned i=0;i<MeasureList.size();i++){
		if(!MeasureList[i]->Prepare(SamplingRate)){
			cout<<"Error Preparing the Measure #"<<i<<"."<<endl;
			return false;
			}
		}
	for(unsigned i=0;i<LogicList.size();i++){
		if(!LogicList[i]->Prepare(SamplingRate)){
			cout<<"Error Preparing the Logic #"<<i<<"."<<endl;
			return false;
			}
		}
	for(unsigned i=0;i<ProtectionList.size();i++){
		if(!ProtectionList[i]->Prepare(SamplingRate)){
			cout<<"Error Preparing the Measure #"<<i<<"."<<endl;
			return false;
			}
		}
/*	for(unsigned i=0;i<ControlList.size();i++){
		if(!ControlList[i]->Prepare(SamplingRate)){
			cout<<"Error Preparing the Control #"<<i<<"."<<endl;
			exit(1);
			}
		}*/
	for(unsigned i=0;i<OscillographyList.size();i++){
		if(!OscillographyList[i]->Prepare(SamplingRate,BufferSize,SystemFrequency)){
			cout<<"Error Preparing the Oscillography #"<<i<<"."<<endl;
			return false;
			}
		}
	for(unsigned i=0;i<OscillographyList.size();i++){
		OscThreadList.create_thread(boost::bind(&Oscillography::Poll, OscillographyList[i]));
		}
#ifdef NCURSES
	for(unsigned i=0;i<DisplayList.size();i++)
		ThreadList.create_thread(boost::bind(&Display::Run, DisplayList[i]));
#endif
	#ifdef RTAI
		long ACQUISITION_TICK=(long)(1E9/SamplingRate);
		rt_allow_nonroot_hrt();
		mlockall(MCL_CURRENT | MCL_FUTURE);
		tarefa_principal = rt_task_init_schmod(nam2num("RELAY"), 2,0, 0, SCHED_FIFO, CPUMAP);
		if (!tarefa_principal){
			printf("Falha ao Iniciar a Tarefa Principal\n");
			return false;
		}
		rt_set_oneshot_mode();
		start_rt_timer(0);
		rt_make_hard_real_time();
		if(rt_task_make_periodic(tarefa_principal,rt_get_time() +100*nano2count(ACQUISITION_TICK), nano2count(ACQUISITION_TICK)))
			cout<<"Error Periodic\n";
		cout<<"Tick:"<<ACQUISITION_TICK<<"Count:"<<nano2count(ACQUISITION_TICK)<<endl;
//		sleep(10);
	#endif

	for(unsigned k=0;k<BufferSize;k++){
		AcqList[0]->RefreshTime();
		for(unsigned i=0;i<AcqList.size();i++)
			Running&=AcqList[i]->Run();
		#ifdef RTAI
			rt_task_wait_period();
		#endif
		}
/*	for(unsigned k=0;k<BufferSize;k++){
		AcqList[0]->RefreshTime();
		for(unsigned i=0;i<AcqList.size();i++)
			Running&=AcqList[i]->Run();
		for(unsigned i=0;i<MeasureList.size();i++)
			MeasureList[i]->Run();
		for(unsigned i=0;i<LogicList.size();i++)
			LogicList[i]->Run();
		#ifdef RTAI
			rt_task_wait_period();
		#endif
		}*/
	Tb=orelay_gettime();
	T1=Tb;
	sizeAcq = AcqList.size();
	sizeMeas = MeasureList.size();
	sizeProt = ProtectionList.size();
	sizeCont = ControlList.size();
	sizeLog = LogicList.size();
	sizeOsc = OscillographyList.size();
	sizeEval = EvaluationList.size();
	cout<<"Relay Ready"<<endl;
	return true; 
	}

void Relay::Stop(){
	Running=false;
	}

void Relay::Start(){
	Running=true;
	}

bool Relay::Reset(){
	for(unsigned k=0;k<BufferSize;k++){
		for(unsigned l=0;l<AnalogicChannels.size();l++)
			AnalogicChannels[l]->insert_Value(0.0);
		for(unsigned l=0;l<ComplexChannels.size();l++)
			ComplexChannels[l]->insert_Value(Complex(0.0,0.0));
		for(unsigned l=0;l<DigitalChannels.size();l++)
			DigitalChannels[l]->insert_Value(false);
		for(unsigned l=0;l<TimerChannels.size();l++)
			TimerChannels[l]->insert_Value(0);
		for(unsigned l=0;l<StringChannels.size();l++)
			StringChannels[l]->insert_Value("");
		MasterClock->insert_Value(0);
		}
	return true;
	}

bool Relay::Step(){
	if(PerformanceEvaluation){
		AcqList[0]->RefreshTime();
	for(unsigned i=0;i<sizeAcq;i++){
		T1=orelay_gettime();			
		Running&=AcqList[i]->Run();
		T2=orelay_gettime();
		ExecutionTimeList[1]->insert_Value((T2-T1)/1E3); //Acquisition
		}
	for(unsigned i=0;i<sizeMeas;i++){
		T1=orelay_gettime();
		MeasureList[i]->Run();
		T2=orelay_gettime();
		ExecutionTimeList[2+i]->insert_Value((T2-T1)/1E3); //Measure
		}
	for(unsigned i=0;i<sizeProt;i++){
		T1=orelay_gettime();
		ProtectionList[i]->Run();
		T2=orelay_gettime();
		ExecutionTimeList[2+i+sizeMeas]->insert_Value((T2-T1)/1E3); //Protection
		}
	for(unsigned i=0;i<sizeCont;i++){
		T1=orelay_gettime();
		ControlList[i]->Run();
		T2=orelay_gettime();
		ExecutionTimeList[2+i+sizeMeas+sizeProt]->insert_Value((T2-T1)/1E3); //Control
		}
	for(unsigned i=0;i<sizeLog;i++){
		T1=orelay_gettime();
		LogicList[i]->Run();
		T2=orelay_gettime();
		ExecutionTimeList[2+i+sizeMeas+sizeProt+sizeCont]->insert_Value((T2-T1)/1E3); //Logic
		}
	for(unsigned i=0;i<sizeOsc;i++){
		T1=orelay_gettime();
		OscillographyList[i]->Run();
		T2=orelay_gettime();
		ExecutionTimeList[2+i+sizeMeas+sizeProt+sizeCont+sizeLog]->insert_Value((T2-T1)/1E3); //OScillography
		}
	for(unsigned i=0;i<sizeEval;i++)
		EvaluationList[i]->Run();
	#ifdef RTAI
		if (rt_task_wait_period()){//Overrun
			Ov++;
			ExecutionOverrun->insert_Value(true);
			cout << "Overrun: " << Ov << endl;
			}
		else{
			ExecutionOverrun->insert_Value(false);
			}
		if(Ov>10000){
			cout<<"OverRun:"<<Ov<<endl;
			Running=false;
			}
	#endif
		T1=orelay_gettime();
		ExecutionTimeList[0]->insert_Value((MasterClock->get_Value(0)-MasterClock->get_Value(1))/1E3); //Cicle Time
		}
	else{
		AcqList[0]->RefreshTime();
		for(unsigned i=0;i<sizeAcq;i++){		
			Running&=AcqList[i]->Run();
			}
		for(unsigned i=0;i<sizeMeas;i++){
			MeasureList[i]->Run();
			}
		for(unsigned i=0;i<sizeProt;i++){
			ProtectionList[i]->Run();
			}
		for(unsigned i=0;i<sizeCont;i++){
			ControlList[i]->Run();
			}
		for(unsigned i=0;i<sizeLog;i++){
			LogicList[i]->Run();
			}
		for(unsigned i=0;i<sizeOsc;i++){
			OscillographyList[i]->Run();
			}
	#ifdef RTAI
		if (rt_task_wait_period()){//Overrun
			Ov++;
			cout << "Overrun: " << Ov << endl;
			}
		if(Ov>10000){
			cout<<"OverRun:"<<Ov<<endl;
			Running=false;
			}
		#endif
		}
	return Running;
	}

Relay::~Relay(){
/*	for(unsigned i=0;i<AcqList.size();i++)
		delete AcqList[i];
	for(unsigned i=0;i<MeasureList.size();i++)
		delete MeasureList[i];
	for(unsigned i=0;i<ProtectionList.size();i++)
		delete ProtectionList[i];
	for(unsigned i=0;i<ControlList.size();i++)
		delete ControlList[i];
	for(unsigned i=0;i<LogicList.size();i++)
		delete LogicList[i];
	for(unsigned i=0;i<OscillographyList.size();i++)
		delete OscillographyList[i];
	for(unsigned i=0;i<EvaluationList.size();i++)
		delete EvaluationList[i];
	for(unsigned i=0;i<AnalogicChannels.size();i++)
		delete AnalogicChannels[i];
	for(unsigned i=0;i<DigitalChannels.size();i++)
		delete DigitalChannels[i];
	for(unsigned i=0;i<StringChannels.size();i++)
		delete StringChannels[i];
	for(unsigned i=0;i<TimerChannels.size();i++)
		delete TimerChannels[i];
	for(unsigned i=0;i<ComplexChannels.size();i++)
		delete ComplexChannels[i];*/
	}

Channel<analogic>* Relay::CreateAnalogicChannel(string name){
	Channel<analogic> *temp=new Channel<analogic>(name,BufferSize);
	AnalogicChannels.push_back(temp);
	return temp;
	}
Channel<digital>* Relay::CreateDigitalChannel(string name){
	Channel<digital> *temp=new Channel<digital>(name,BufferSize);
	DigitalChannels.push_back(temp);
	return temp;
	}
Channel<string>* Relay::CreateStringChannel(string name){
	Channel<string> *temp=new Channel<string>(name,BufferSize);
	StringChannels.push_back(temp);
	return temp;
	}
Channel<timer>* Relay::CreateTimerChannel(string name){
	Channel<timer> *temp=new Channel<timer>(name,BufferSize);
	TimerChannels.push_back(temp);
	return temp;

	TimerChannels.push_back(new Channel<timer>(name,BufferSize));
	}
Channel<Complex>* Relay::CreateComplexChannel(string name){
	Channel<Complex> *temp=new Channel<Complex>(name,BufferSize);
	ComplexChannels.push_back(temp);
	return temp;
	}

void Relay::JoinChannel(Channel<analogic> *temp){
	temp->set_Size(BufferSize);
	AnalogicChannels.push_back(temp);
	}

void Relay::JoinChannel(Channel<digital> *temp){
	temp->set_Size(BufferSize);
	DigitalChannels.push_back(temp);
	}

void Relay::JoinChannel(Channel<string> *temp){
	temp->set_Size(BufferSize);
	StringChannels.push_back(temp);
	}

void Relay::JoinChannel(Channel<timer> *temp){
	temp->set_Size(BufferSize);
	TimerChannels.push_back(temp);
	}

void Relay::JoinChannel(Channel<Complex> *temp){
	temp->set_Size(BufferSize);
	ComplexChannels.push_back(temp);
	}

Channel<analogic>* Relay::GetAnalogicChannel(string name){
	for(unsigned i=0;i<AnalogicChannels.size();i++)
		if(name==AnalogicChannels[i]->get_Name())
			return AnalogicChannels[i];
	return NULL;
	}
Channel<digital>* Relay::GetDigitalChannel(string name){
	for(unsigned i=0;i<DigitalChannels.size();i++)
		if(name==DigitalChannels[i]->get_Name())
			return DigitalChannels[i];
	return NULL;
	}
Channel<string>* Relay::GetStringChannel(string name){
	for(unsigned i=0;i<StringChannels.size();i++)
		if(name==StringChannels[i]->get_Name())
			return StringChannels[i];
	return NULL;
	}
Channel<timer>* Relay::GetTimerChannel(string name){
	for(unsigned i=0;i<TimerChannels.size();i++)
		if(name==TimerChannels[i]->get_Name())
			return TimerChannels[i];
	return NULL;
	}
Channel<Complex>* Relay::GetComplexChannel(string name){
	for(unsigned i=0;i<ComplexChannels.size();i++)
		if(name==ComplexChannels[i]->get_Name())
			return ComplexChannels[i];
	return NULL;
	}


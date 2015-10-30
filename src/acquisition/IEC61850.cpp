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
#ifdef IEC61850
#include "IEC61850.h"



namespace orelay{
IECPublisher::IECPublisher(MacAddress destMac, MacAddress localMac, unsigned short vlanid, unsigned char canonical, unsigned char vlanpri, unsigned short appid, long TMin, long TMax, string gocbref, string datset, string Goid, string Dev){
	goose=new GooseMessage();
	goose->setDestMac(destMac);
	goose->setSrcMac(localMac);
	goose->setVlanId(vlanid);
	goose->setCanonical(canonical);
	goose->setVlanPrio(vlanpri);
	goose->setAppId(appid);
	goose->setTimeAllowedToLive(TMax);
	goose->setGocbRef(gocbref);
	goose->setDatSet(datset);
	goose->setGoId(Goid);

	goose->setTest(false);
	goose->setConfRev(1);
	goose->setNdsCom(false);

	TimetoLiveMin=TMin;
	TimetoLiveMax=TMax;
	MessageType=GOOSE;
	resetTimetoLive();
	sock=0;
	Device=Dev;
	LastSentMsgTime=0;
	channel_cnt=0;

	}

void IECPublisher::Join_Channel(Channel<digital> *Ch){
	if(!Ch->InUse()){
		Ch->UseIt();
		DigitalOut.push_back(Ch);
		channel_cnt++;
		DigitalOutBoardChannel.push_back(channel_cnt);
		goose->addBitStringData();
		channel_cnt++; 
		goose->addBooleanData();//1
		}
	else{
		cout<<"IECPublisher::Error:: Channel" <<Ch->get_Name()<<" Already in Use."<<endl;
		exit(0);
		}
	}


void IECPublisher::Join_Channel(Channel<analogic> *Ch ){
	if(!Ch->InUse()){
		Ch->UseIt();
		AnalogicOut.push_back(Ch);
		channel_cnt++;
		AnalogicOutBoardChannel.push_back(channel_cnt);
		AnalogicOutRatio.push_back(1.0);
		AnalogicOutOffset.push_back(0.0);
		goose->addBitStringData(); 
		channel_cnt++; 
		goose->addFloatData(0.5);//1
		}
	else{
		cout<<"IECPublisher::Error:: Channel" <<Ch->get_Name()<<" Already in Use."<<endl;
		exit(0);
		}
	}


void IECPublisher::RefreshTime(){
	MasterClock->insert_Value(orelay_gettime());
	}

bool IECPublisher::Prepare(float){ 
	return Publisher::Prepare(); 
	};

 bool IECPublisher::Run(){
	for(unsigned int i=0;i<DigitalOut.size();i++){
		goose->setBooleanData(DigitalOutBoardChannel[i],DigitalOut[i]->get_Value());
		}
	for(unsigned int i=0;i<AnalogicOut.size();i++){
		goose->setFloatData(AnalogicOutBoardChannel[i],AnalogicOut[i]->get_Value());
		}
	return Publisher::Execute();
	};


IECSubscriber::IECSubscriber(unsigned type,string Gocb, string Dev){
	GocbRef=Gocb;
	MessageType=type;
	sock=0;
	Device=Dev;
	if(MessageType==GOOSE){
		goose=new GooseMessage();

		}
	else{
		//SV=new SVMessage();
		}
	}

void IECSubscriber::RefreshTime(){
	MasterClock->insert_Value(orelay_gettime());
	}

void IECSubscriber::Join_Channel(Channel<digital> *Ch, unsigned pos){
	if(!Ch->InUse()){
		Ch->UseIt();
		DigitalIn.push_back(Ch);
		DigitalInBoardChannel.push_back(pos);
		}
	else{
		cout<<"IECPublisher::Error:: Channel" <<Ch->get_Name()<<" Already in Use."<<endl;
		exit(0);
		}
	}

void IECSubscriber::Join_Channel(Channel<analogic> *Ch, unsigned pos){
	if(!Ch->InUse()){
		Ch->UseIt();
		AnalogicIn.push_back(Ch);
		AnalogicInBoardChannel.push_back(pos);
		AnalogicInRatio.push_back(1.0);
		AnalogicInOffset.push_back(0.0);
		}
	else{
		cout<<"IECPublisher::Error:: Channel" <<Ch->get_Name()<<" Already in Use."<<endl;
		exit(0);
		}
	}

bool  IECSubscriber::Prepare(float){
	if(!Subscriber::Prepare())
		return false;
	if(MessageType==GOOSE){
		for(unsigned int i=0;i<DigitalIn.size();i++){
			if(!goose->isBooleanData(DigitalInBoardChannel[i]))
				return false;
			}
		for(unsigned int i=0;i<AnalogicIn.size();i++){
			if(!goose->isFloatData(AnalogicInBoardChannel[i]))
				return false;
			}
		}
	return true;
}

bool IECSubscriber::Run(){
	Subscriber::Execute();
	if(MessageType==GOOSE){
		for(unsigned int i=0;i<DigitalIn.size();i++){
			DigitalIn[i]->insert_Value(goose->getBooleanData(DigitalInBoardChannel[i]));
			}
		for(unsigned int i=0;i<AnalogicIn.size();i++){
			AnalogicIn[i]->insert_Value(goose->getFloatData(AnalogicInBoardChannel[i]));
			//cout << "getFloatValue: " << goose->getFloatData(AnalogicInBoardChannel[i]) << endl;
			}
		}
	return true;
}

/*
svMessage sv;
	sv.setDestMac(destMac);
	sv.setSrcMac(localMac);
	sv.setVlanId(0);
	sv.setCanonical(0);
	sv.setVlanPrio(4);
	sv.setAppId(0x4000);
	sv.addASDU(1);
	sv.setsvID(0, "000100MU0444");
	sv.setsmpCnt(0,0);
	sv.setsmpSynch(0,1);
	sv.setConfRef(0,0000000001);


IECPublisher::IECPublisher(MacAddress destMac, MacAddress localMac, unsigned short vlanid, unsigned char canonical, unsigned char vlanpri, unsigned short appid, long TMin, long TMax, string gocbref, string datset, string Goid,unsigned T, string Dev){

	MessageType=SV;
	sock=0;
	Device=Dev;
	LastSentMsgTime=0;
	Skip=T;
//	HouseKeepingRX=false;

	sv.setDestMac(destMac);
	sv.setSrcMac(localMac);
	sv.setVlanId(vlanid);
	sv.setCanonical(canonical);
	sv.setVlanPrio(vlanpri);
	sv.setAppId(appid);
	sv.setsvID(0, "000100MU0444");
	sv.setsmpCnt(0,0);
	sv.setsmpSynch(0,1);
	sv.setConfRef(0,0000000001);
//ressalva
	sv.addASDU(1);
	}
IECPublisher::IECPublisher(string Dev){

	MessageType=sv;
	sock=0;
	Device=Dev;
	LastSentMsgTime=0;
	Skip=0;
	}*/

}
#endif

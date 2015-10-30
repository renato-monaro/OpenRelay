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
/*
 * GooseMessage.cpp
 *
 *  Created on: 07/06/2011
 *      Author: rapphil
 */

#include <string.h>
#include <iostream>
#include <stdio.h>
#include <netinet/in.h>
#include "GooseMessage.h"

using namespace orelay;
GooseMessage::GooseMessage(){
	this->message=NULL;
	this->header.h_proto=VLAN_MESSAGE;
	this->ethrType=GOOSE_MESSAGE;
	this->reserved1=RESERVED_1_VALUE;
	this->reserved2=RESERVED_2_VALUE;
	this->numDatSetsEntries=0;
	this->sqNum=0;
	this->stNum=0;
}
GooseMessage::GooseMessage(unsigned char * outsideMessage,unsigned size) {
	fromNetwork(outsideMessage,size);
}

void GooseMessage::fromNetwork(unsigned char * outsideMessage,unsigned size){
	if(this->message!=NULL){
		delete this->message;
		for (unsigned count=0;count<this->dataPointers.size();count++)
			delete this->dataPointers.at(count);
		this->dataTypes.clear();
		this->dataSizes.clear();
		this->dataPointers.clear();
	}
	this->message = new unsigned char[size];
	memcpy(message,outsideMessage,size);
	unsigned char * auxmessage=message;
	this->size=size;
	this->header=*((struct ethhdr *)auxmessage);
	this->header.h_proto=ntohs(this->header.h_proto);
	auxmessage+=sizeof(struct ethhdr);
	unsigned short vlanAttr=ntohs(*((unsigned short*)auxmessage));	auxmessage+=2; //atributos do vlan Priority network, Canonical e ID
	this->priority=(vlanAttr>>13)&0b111;
	this->canonical=(vlanAttr>>12)&0b1;
	this->idVlan=vlanAttr&0b111111111111;
	this->ethrType=ntohs(*((unsigned short*)auxmessage));  auxmessage+=2;//  tipo da mensagem ethernet IEC 61850 (0x88b8)
	this->appId=ntohs(*((unsigned short*)auxmessage));  auxmessage+=2; // appId
	this->length=ntohs(*((unsigned short*)auxmessage));  auxmessage+=2; // length
	this->reserved1=ntohs(*((unsigned short*)auxmessage));  auxmessage+=2; // reserved 1
	this->reserved2=ntohs(*((unsigned short*)auxmessage));  auxmessage+=3; // reserved 2
	unsigned char sizeAux = *((unsigned char *)auxmessage); auxmessage++;
	if(sizeAux==0x81)
		auxmessage+=2;
	else if (sizeAux==0x82)
		auxmessage+=3;
	else
		auxmessage++;
	sizeAux= *((unsigned char *)auxmessage); auxmessage++;
	memcpy(this->gocbRef,auxmessage,sizeAux); this->gocbRef[sizeAux]='\0'; auxmessage+=sizeAux+1;
    sizeAux = *((unsigned char *)auxmessage); auxmessage++;
    if(sizeAux==2){
    	this->timeAllowedToLive=ntohs(*((unsigned short*)auxmessage));
    	auxmessage+=3;
    }
    else{
    	this->timeAllowedToLive=*auxmessage; auxmessage+=2;
    }// appId
    sizeAux = *((unsigned char *)auxmessage); auxmessage++;
    memcpy(this->datSet,auxmessage,sizeAux); this->datSet[sizeAux]='\0'; auxmessage+=sizeAux+1;
    sizeAux = *((unsigned char *)auxmessage); auxmessage++;
    memcpy(this->goID,auxmessage,sizeAux); this->goID[sizeAux]='\0'; auxmessage+=sizeAux+1;
    sizeAux = *((unsigned char *)auxmessage); auxmessage++;
    memcpy(this->t,auxmessage,sizeAux); this->t[sizeAux]='\0'; auxmessage+=sizeAux+1;
	int *temp=(int *)this->t;
    long seconds=*temp;
	temp=(int*)(this->t+4);
    long nanoseconds=*temp;
    seconds = ntohl(seconds);
    nanoseconds= ntohl(nanoseconds);
    nanoseconds=nanoseconds/4.294967296;

    this->time=seconds*1e9+nanoseconds;

    sizeAux = *((unsigned char *)auxmessage); auxmessage++;
    if(sizeAux==2){
    	this->stNum=ntohs(*((unsigned short*)auxmessage));
    	auxmessage+=3;
    }
    else{
    	this->stNum=*auxmessage; auxmessage+=2;
    }
    sizeAux = *((unsigned char *)auxmessage); auxmessage++;
    if(sizeAux==2){
        	this->sqNum=ntohs(*((unsigned short*)auxmessage));
        	auxmessage+=3;
    }else{
        	this->sqNum=(unsigned char)*auxmessage; auxmessage+=2;
    }
    sizeAux = *((unsigned char *)auxmessage); auxmessage++;
    this->test=*auxmessage; auxmessage+=2;
    sizeAux = *((unsigned char *)auxmessage); auxmessage++;
    if(sizeAux==2){
       	this->confRev=ntohs(*((unsigned short*)auxmessage));
       	auxmessage+=3;
    }else{
       	this->confRev=*auxmessage; auxmessage+=2;
    }
    sizeAux = *((unsigned char *)auxmessage); auxmessage++;
    this->ndsCom=*auxmessage; auxmessage+=2;
    sizeAux = *((unsigned char *)auxmessage); auxmessage++;
    if(sizeAux==2){
       	this->numDatSetsEntries=ntohs(*((unsigned short*)auxmessage));
       	auxmessage+=3;
    }else{
       	this->numDatSetsEntries=*auxmessage; auxmessage+=2;
    }
    sizeAux = *((unsigned char *)auxmessage); auxmessage++;
    if(sizeAux==0x81)
    	auxmessage+=1;
    else if(sizeAux==0x82)
    	auxmessage+=2;
    unsigned char charSize;
    unsigned char * data;
    this->dataSetSize=0;
    for (unsigned count=0;count<this->numDatSetsEntries;count++){
    	this->dataTypes.push_back(*((unsigned char *)auxmessage)); auxmessage++;
    	charSize=*((unsigned char *)auxmessage); auxmessage++;
    	this->dataSizes.push_back(charSize);
    	data = new unsigned char[charSize];
    	this->dataSetSize+=2+charSize;
    	memcpy(data,auxmessage,charSize); auxmessage+=charSize;
    	this->dataPointers.push_back(data);
    }
}

void GooseMessage::printDataSet(){
	char buff[10];

	for (unsigned count =0;count<this->numDatSetsEntries;count++){
		sprintf(buff,"%00x",this->dataTypes.at(count));
		std::cout <<buff<<" ";
		sprintf(buff,"%00x",this->dataSizes.at(count));
		std::cout<<buff<<" ";
		for (int count2=0;count2<this->dataSizes.at(count);count2++){
			sprintf(buff,"%00x",this->dataPointers.at(count)[count2]);
			std::cout << buff << " ";
		}
		std::cout<< std::endl;
	}
}

void GooseMessage::toNetwork(unsigned char * outsideMessage,unsigned * outsideSize){
	this->updateSizes();
	delete message;
	message = new unsigned char[size];
	memset(message,'\0',size);
	unsigned char * auxmessage;
	unsigned short * shptr;
	auxmessage=message;
	struct ethhdr * hdr;
	hdr = (struct ethhdr *)message;
	*hdr = this->header;
	hdr->h_proto=htons(hdr->h_proto);
	auxmessage+=sizeof(struct ethhdr);
	shptr=(unsigned short *) auxmessage;
	*shptr=htons((this->priority<<13)+(this->canonical<<12)+(this->idVlan&0xFFF));	auxmessage+=2; shptr=(unsigned short *) auxmessage; //atributos do vlan Priority network, Canonical e ID
	*shptr=htons(this->ethrType); auxmessage+=2; shptr=(unsigned short *) auxmessage; //eth type goose 0x88b8
	*shptr=htons(this->appId); auxmessage+=2; shptr=(unsigned short *) auxmessage; //app ID
	*shptr=htons(this->length); auxmessage+=2; shptr=(unsigned short *) auxmessage; //length
	*shptr=htons(this->reserved1); auxmessage+=2; shptr=(unsigned short *) auxmessage; // reserved 1
	*shptr=htons(this->reserved2); auxmessage+=2; // reserved 2
	*auxmessage=0x61; auxmessage++;
	if ( this->goosePduSize>0xff){
		*auxmessage=0x82; auxmessage++;
		*((unsigned short *)auxmessage)=htons(this->goosePduSize); auxmessage+=2;
	}else if(this-> goosePduSize>0x7f){
		*auxmessage=0x81; auxmessage++;
		*auxmessage=(unsigned char)this->goosePduSize; auxmessage++;
	}else{
		*auxmessage=(unsigned char)this->goosePduSize; auxmessage++;
	}


	*auxmessage=0x80; auxmessage++;
	*auxmessage=(unsigned char)strlen(this->gocbRef); auxmessage++;
	memcpy(auxmessage,this->gocbRef,strlen(this->gocbRef));auxmessage+=strlen(this->gocbRef);
	*auxmessage=0x81; auxmessage++;//*auxmessage=0x02;auxmessage++;shptr=(unsigned short*)auxmessage;
	*auxmessage=(this->timeAllowedToLive>0x7f)?0x02:0x01;auxmessage++;
	if(this->timeAllowedToLive>0x7f){
		*((unsigned short *)auxmessage)=htons(this->timeAllowedToLive); auxmessage+=2;
	}else{
		*auxmessage=(unsigned char)this->timeAllowedToLive; auxmessage++;
	}
	//*shptr=htons(this->timeAllowedToLive); auxmessage+=2;
	*auxmessage=0x82; auxmessage++;
	*auxmessage=(unsigned char)strlen(this->datSet); auxmessage++;
	memcpy(auxmessage,this->datSet,strlen(this->datSet));auxmessage+=strlen(this->datSet);
	*auxmessage=0x83; auxmessage++;
	*auxmessage=(unsigned char)strlen(this->goID); auxmessage++;
	memcpy(auxmessage,this->goID,strlen(this->goID));auxmessage+=strlen(this->goID);
	*auxmessage=0x84; auxmessage++;
	//*auxmessage=(unsigned char)strlen((char *)this->t); auxmessage++;
	//memcpy(auxmessage,this->t,strlen((char *)this->t));auxmessage+=strlen((char *)this->t);

	*auxmessage=0x08; auxmessage++;
	memcpy(auxmessage,this->t,0x08); auxmessage+=0x08;

	*auxmessage=0x85; auxmessage++;
	*auxmessage=(this->stNum>0x7f)?0x02:0x01;auxmessage++;
	if(this->stNum>0x7f){
		*((unsigned short *)auxmessage)=htons(this->stNum); auxmessage+=2;
	}else{
		*auxmessage=(unsigned char)this->stNum; auxmessage++;
	}
	*auxmessage=0x86; auxmessage++;
	*auxmessage=(this->sqNum>0x7f)?0x02:0x01;auxmessage++;
	if(this->sqNum>0x7f){
		*((unsigned short *)auxmessage)=htons(this->sqNum); auxmessage+=2;
	}else{
		*auxmessage=(unsigned char)this->sqNum; auxmessage++;
	}
	*auxmessage=0x87; auxmessage++;*auxmessage=0x01;auxmessage++;
	*auxmessage=this->test; auxmessage++;
	*auxmessage=0x88; auxmessage++;
	*auxmessage=(this->confRev>0x7f)?0x02:0x01;auxmessage++;
	if(this->confRev>0x7f){
		*((unsigned short *)auxmessage)=htons(this->confRev); auxmessage+=2;
	}else{
		*auxmessage=(unsigned char)this->confRev; auxmessage++;
	}
	*auxmessage=0x89; auxmessage++;*auxmessage=0x01;auxmessage++;
	*auxmessage=this->ndsCom; auxmessage++;
	*auxmessage=0x8a; auxmessage++;
	*auxmessage=(this->numDatSetsEntries>0x7f)?0x02:0x01;auxmessage++;
	if(this->numDatSetsEntries>0x7f){
		*((unsigned short *)auxmessage)=htons(this->numDatSetsEntries); auxmessage+=2;
	}else{
		*auxmessage=(unsigned char)this->numDatSetsEntries; auxmessage++;
	}
	*auxmessage=0xab; auxmessage++;


	if ( this->dataSetSize>0xff){
		*auxmessage=0x82; auxmessage++;
		*((unsigned short *)auxmessage)=htons(this->dataSetSize); auxmessage+=2;
	}else if(this->dataSetSize>0x7f){
		*auxmessage=0x81; auxmessage++;
		*auxmessage=(unsigned char)this->dataSetSize; auxmessage++;
	}else{
		*auxmessage=(unsigned char)this->dataSetSize; auxmessage++;
	}


	for (unsigned count=0;count<this->dataSizes.size();count++){
		*auxmessage=this->dataTypes.at(count); auxmessage++;
		*auxmessage=this->dataSizes.at(count); auxmessage++;
		for (unsigned count2=0;count2<this->dataSizes.at(count);count2++){
			*auxmessage=this->dataPointers.at(count)[count2];
			auxmessage++;
		}

	}
	memcpy(outsideMessage,message,size);
	*outsideSize=size;
}

void GooseMessage::printMessage(){
	unsigned count,contbytes=0;
	char buff[6];
	for (count=0;count<size;count++){
		sprintf(buff," %2.2x ",message[count]);
		std::cout<< buff;
		contbytes++;
		if(contbytes==16){
			std::cout<<std::endl;
			contbytes=0;
		}
	}
	std::cout<<std::endl;

}

void GooseMessage::updateSizes(){
	this->dataSetSize=0;
	for(unsigned count =0;count<this->numDatSetsEntries;count++)
		this->dataSetSize+=2+this->dataSizes.at(count);
	this->goosePduSize=this->dataSetSize;
	if(this->dataSetSize>0xff)
		this->goosePduSize+=4; // tamanho dados + 0xab 0x82 0xtt 0xtt
	else if(this->dataSetSize>0x7f)
		this->goosePduSize+=3;  // tamanho dados + 0xab 0x81 0xtt
	else
		this->goosePduSize+=2; // tamanho dados + 0xab 0xtt

	this->goosePduSize+=(this->numDatSetsEntries>0x7f)? 4: 3;// (numDatSetEntries) 0x8a 0x02 0xtt 0xtt | 0x8a 0x01 0xtt
	this->goosePduSize+=6; // (ndsCom test) 0x89 0x01 0xtt | 0x87 0x01 0xtt
	this->goosePduSize+=(this->confRev>0x7f)?4:3; //(confRev) 0x88 0x02 0xtt 0xtt | 0x88 0x01 0xtt
	this->goosePduSize+=(this->sqNum>0x7f)?4:3; //(confRev) 0x86 0x02 0xtt 0xtt | 0x86 0x01 0xtt
	this->goosePduSize+=(this->stNum>0x7f)?4:3; //(confRev) 0x86 0x02 0xtt 0xtt | 0x86 0x01 0xtt
	this->goosePduSize+=10; // (t) 0x84 0x08 0xtt 0xtt 0xtt 0xtt 0xtt 0xtt 0xtt 0xtt
	this->goosePduSize+=strlen(this->goID)+2; // (goId) length (goid) + 0x83 0xtt
	this->goosePduSize+=strlen(this->datSet)+2; // (datSet) length (datSet) + 0x82 + 0xtt
	this->goosePduSize+=(this->timeAllowedToLive>0x7f)?4:3; //(timeAllowedToLive) 0x81 0x02 0xtt 0xtt | 0x81 0x01 0xtt
	this->goosePduSize+=strlen(this->gocbRef)+2; // (gocbRef) length (gocbRef) + 0x80  0xtt
	this->length=8; // APPID LENGHT REserverd 1 Reserved 2
	if(this->goosePduSize>0xff)
		this->length+=4 + goosePduSize; // 0x61 0x82 0xtt 0xtt
	else if(this->goosePduSize>0x7f)
		this->length+=3+ goosePduSize; // 0x61 0x81 0xtt
	else
		this->length+=2+goosePduSize; // 0x61 0xtt
	this->size = this->length+ 18;  // MACDEST MACSOURCE 0x8100(VLAN) ATTRVLAN 0x88b8(IEC 61850)
}



void GooseMessage::memcpyinv(unsigned char * dest, unsigned char * source, unsigned short t){
	for (unsigned short count=0;count<t;count++)
		dest[t-1-count]=source[count];
}






/*
void GooseMessage::setVlanId(unsigned short id){
	this->idVlan=id;
}*/
/*
void GooseMessage::setVlanPrio(unsigned char prio){
	this->priority=prio;
}*/
/*
char* GooseMessage::getdataSet(){
	return datSet;
}*/




void GooseMessage::setSrcMac(unsigned char * mac){
	memcpy(&(header.h_source),mac,ETH_ALEN);
}

void GooseMessage::setDestMac(unsigned char * mac){
	memcpy(&(header.h_dest),mac,ETH_ALEN);
}


/*
unsigned short GooseMessage::getVlanId(){
	return idVlan;
}
unsigned char GooseMessage::getVlanPrio(){
	return priority;
}*/



void GooseMessage::setIdVlan(unsigned short id){
	this->idVlan=id;
}

unsigned short GooseMessage::getIdVlan(){
	return this->idVlan;
}

void GooseMessage::setCanonical(unsigned char canonical){
	this->canonical=canonical;
}

unsigned char GooseMessage::getCanonical(){
	return this->canonical;
}

void GooseMessage::setPriority(unsigned char priority){
	this->priority = priority;
}

unsigned char GooseMessage::getPriority(){
	return this->priority;
}

void GooseMessage::setAppId(unsigned short appId){
	this->appId=appId;
}

unsigned short GooseMessage::getAppId(){
	return appId;
}

void GooseMessage::setGocbRef(char * gocbRef){
	strcpy(this->gocbRef,gocbRef);
}

char* GooseMessage::getGocbRef(){
	return gocbRef;
}

void GooseMessage::setTimeAllowedToLive(unsigned short time){
	this->timeAllowedToLive= time;
}

unsigned short GooseMessage::getTimeAllowedToLive(){
	return this->timeAllowedToLive;
}

void GooseMessage::setDatSet(char * datSet){
	strcpy(this->datSet,datSet);
}
char* GooseMessage::getDatSet(){
	return this->datSet;
}


void GooseMessage::setGoId(char *Id){
	strcpy(goID,Id);
}

char* GooseMessage::getGoId(){
	return goID;
}

void GooseMessage::setTime(long long tm){
	long seconds,nanoseconds;
	seconds = tm/1e9;
	nanoseconds=tm-1e9*seconds;
	nanoseconds= nanoseconds*(0x100000000/1000000000);
	int *temp;
	temp=(int *)t;
	*temp=htonl(seconds);
	temp=(int *)(t+4);
	*temp=htonl(nanoseconds);
//	*(int *)t=htons(seconds);
//	*(int *)(t+4)=htons(nanoseconds);
	this->time=tm;
}

long long GooseMessage::getTime(){
	return this->time;
}

void GooseMessage::setStNum(unsigned short num){
	this->stNum=num;
}

unsigned short GooseMessage::getStNum(){
	return stNum;
}


void GooseMessage::setSqNum(unsigned short num ){
	this->sqNum=num;
}

unsigned short GooseMessage::getSqNum(){
	return this->sqNum;
}

void GooseMessage::incSqNum(){
	this->sqNum++;
}

void GooseMessage::incStNum(){
	this->stNum++;
}


void GooseMessage::resetSqNum(){
	this->sqNum=0;
}


void GooseMessage::setTest(bool test){
	this->test=test;
}
bool GooseMessage::getTest(){
	return test;
}

void GooseMessage::setConfRev(unsigned char confRev){
	this->confRev=confRev;
}
unsigned char GooseMessage::getConfRev(){
	return this->confRev;
}

void GooseMessage::setNdsCom(bool ndsCom){
	this->ndsCom=ndsCom;
}
bool GooseMessage::getNdsCom(){
	return this->ndsCom;
}

unsigned GooseMessage::getNumDatSetEntries(){
	return this->numDatSetsEntries;
}

unsigned char GooseMessage::getDataType(unsigned position){
	if(position>=this->dataTypes.size())
		return NO_DATA_TYPE;
	return this->dataTypes.at(position);
}

void GooseMessage::addBooleanData(){
	this->dataSizes.push_back(1);
	this->dataTypes.push_back(BOOL_DATA_TYPE);
	unsigned char * data = new unsigned char[1];
//	*data=value;
	this->dataPointers.push_back(data);
	this->numDatSetsEntries++;
}

void GooseMessage::addFloatData(){
	this->dataSizes.push_back(5);
	this->dataTypes.push_back(FLOAT_DATA_TYPE);
	unsigned char * data = new unsigned char[5];
	this->dataPointers.push_back(data);
//	data[0]=0x8;
//	data++;
//	float auxv=11.5;
//	unsigned char *p=(unsigned char *) &auxv;
//	memcpyinv(data,p,sizeof(float));
	this->numDatSetsEntries++;

}

void GooseMessage::setFloatData(unsigned position,float value){
	if(this->dataTypes.at(position)==FLOAT_DATA_TYPE){
		unsigned char * data = this->dataPointers.at(position);
		data[0]=0x8;
		data++;
		float auxv=value;
		unsigned char *p=(unsigned char *) &auxv;
		memcpyinv(data,p,sizeof(float));
	}
}

void GooseMessage::setBooleanData(unsigned position,bool value){
	if(this->dataTypes.at(position)==BOOL_DATA_TYPE){
		unsigned char * data = this->dataPointers.at(position);
		//data++;
		unsigned char auxv=value;
		unsigned char *p=(unsigned char *) &auxv;
		memcpyinv(data,p,sizeof(unsigned char));
		
	}
}

bool GooseMessage::getBooleanData(unsigned position){
	unsigned type = this->dataTypes.at(position);
	if(type==BOOL_DATA_TYPE){
		bool returnValue;
		returnValue = (bool)(*this->dataPointers.at(position));
		return returnValue;
	}
	return false;
}

float GooseMessage::getFloatData(unsigned position){
	unsigned type = this->dataTypes.at(position);
	if(type==FLOAT_DATA_TYPE){
		unsigned char * data = this->dataPointers.at(position);
		data++;
		float returnValue;
		unsigned char * p = (unsigned char *) &returnValue;
		memcpyinv(p,data,sizeof(float));
		return returnValue;
	}

	return 0;
}

GooseMessage::~GooseMessage() {
}


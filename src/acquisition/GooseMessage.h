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
 * GooseMessage.h
 *
 *  Created on: 07/06/2011
 *      Author: rapphil
 */

#include <linux/if_ether.h>
#include <time.h>
#include <vector>

#ifndef GOOSEMESSAGE_H_
#define GOOSEMESSAGE_H_


#define BOOL_DATA_TYPE 0x83
#define BIT_STRING_DATA_TYPE 0x84
#define FLOAT_DATA_TYPE 0x87
#define NO_DATA_TYPE 0x0
#define GOOSE_MESSAGE 0x88b8
#define VLAN_MESSAGE 0x8100

#define RESERVED_1_VALUE 0x0000
#define RESERVED_2_VALUE 0x0000
namespace orelay{
class GooseMessage {
public:
	GooseMessage();
	GooseMessage(unsigned char *,unsigned);
	void setDestMac(unsigned char *);
	void setSrcMac(unsigned char *);

	void setIdVlan(unsigned short);
	unsigned short getIdVlan();

	void setCanonical(unsigned char);
	unsigned char getCanonical();

	void setPriority(unsigned char);
	unsigned char getPriority();

	void setAppId(unsigned short);
	unsigned short getAppId();

	void setGocbRef(char *);
	char* getGocbRef();

	void setTimeAllowedToLive(unsigned short);
	unsigned short getTimeAllowedToLive();

	void setDatSet(char *);
	char* getDatSet();

	char* getGoId();
	void setGoId(char *);

	long long getTime();
	void setTime(long long);

	void setStNum(unsigned short);	
	unsigned short getStNum();
		
	void setSqNum(unsigned short);
	unsigned short getSqNum();
	void incSqNum();
	void incStNum();
	void resetSqNum();

	void setTest(bool);
	bool getTest();

	void setConfRev(unsigned char);
	unsigned char getConfRev();

	void setNdsCom(bool);
	bool getNdsCom();

	void fromNetwork(unsigned char *,unsigned);
	void toNetwork(unsigned char *,unsigned *);

	void printDataSet();
	void printMessage();
	
	unsigned getNumDatSetEntries();
	void addBooleanData();
	void addFloatData();

	unsigned char getDataType(unsigned);
	float getFloatData(unsigned);
	bool getBooleanData(unsigned);

	void setFloatData(unsigned,float);
	void setBooleanData(unsigned,bool);

	virtual ~GooseMessage();

private:
	void memcpyinv(unsigned char * , unsigned char * , unsigned short);
	void updateSizes();
	unsigned size;
	unsigned char *message;
	struct ethhdr header;
	unsigned char priority;
	unsigned char canonical;
	unsigned short idVlan;
	unsigned short ethrType;
	unsigned short appId;
	unsigned short length;
	unsigned short reserved1;
	unsigned short reserved2;
	char gocbRef[60];
	unsigned short timeAllowedToLive;
	char datSet[60];
	char goID[60];
	unsigned char t[60];
	long long time;
	unsigned short stNum;
	unsigned short sqNum;
	bool test;
	unsigned char confRev;
	bool ndsCom;
	unsigned numDatSetsEntries;
	unsigned char *allData;
	unsigned char * types;
	unsigned char ** addrs;
	unsigned dataSetSize;
	unsigned goosePduSize;
	std::vector<unsigned char>dataSizes;
	std::vector<unsigned char *> dataPointers;
	std::vector<unsigned char> dataTypes;
};
}

#endif /* GOOSEMESSAGE_H_ */

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
#include "acquisitiongoose.h"

#ifdef RTAI
#ifdef RTNET
GooseAcquisition::GooseAcquisition(char *OutMac,char *goid, unsigned short VId, unsigned short VPrio, unsigned short TimetoLive, char *Device){
	FunctionalDescription_set("Goose Acquisition"); 
 	struct sockaddr_ll addr;
	unsigned char *OwnMac;
	OwnMac=(unsigned char*)malloc(MAC_ADDR_LEN);
	BroadcastMac=(unsigned char*)malloc(MAC_ADDR_LEN);
	rt_eth_aton(BroadcastMac, OutMac);
	strcpy(IED,goid);
		
	Sock=rt_dev_socket(PF_PACKET, SOCK_RAW, htons(ETH_P_ALL));
	if(Sock<0){
		cout<<"Acquisition::Goose::Error Creating Socket"<<endl;
		exit(0);
		}
	struct ifreq ifr;
	strncpy(ifr.ifr_name, Device, IFNAMSIZ);
	if (rt_dev_ioctl(Sock, SIOCGIFINDEX, &ifr) < 0) {
		cout<<"Acquisition::Goose::Error Getting Device"<<endl;
		close(Sock);
		exit(0);
		}
	addr.sll_family= AF_PACKET;
	addr.sll_protocol=htons(ETH_P_ALL);
	addr.sll_ifindex= ifr.ifr_ifindex;
//	addr.sll_ifindex= 1;
	addr.sll_halen=6;
	if(rt_dev_ioctl(Sock,SIOCGIFHWADDR,&ifr)<0){
	cout<<"Can not get interface address"<<endl;
	close(Sock);
	exit(0);
	}

    if (rt_dev_bind(Sock, (struct sockaddr *)&addr, sizeof(addr)) < 0) {
		cout<<"Acquisition::Goose::Error Binding Socket"<<endl;
		close(Sock);
		exit(0);
		}
	Packet = (unsigned char*)malloc(MAX_FRAME_SIZE);
	if(Packet<0){
		cout<<"Acquisition::Goose::Error Malloc Buffer"<<endl;
		close(Sock);
		exit(0);
		}
	memset(Packet, '\0', MAX_FRAME_SIZE);
	EtherHeader=(struct ether_header*)Packet;

    if (rt_dev_ioctl(Sock, SIOCGIFHWADDR, &ifr) < 0) {
       cout<<"Acquisition::Goose::Error Getting Interface Address"<<endl;
        close(Sock);
        exit(0);
    }

	memcpy(OwnMac, ifr.ifr_hwaddr.sa_data, MAC_ADDR_LEN);
	cout<<"OwnMac:"; 
	for(unsigned i=0;i<6;i++)
		cout<<hex<<ifr.ifr_hwaddr.sa_data[i];
	cout<<endl;
	Delay=MIN_PUBLISHING_TIME;
	OutputMessage.resetSqNum();
	OutputMessage.setDestMac(BroadcastMac);
	OutputMessage.setSrcMac(OwnMac);
	OutputMessage.setIdVlan(VId);
	OutputMessage.setPriority(VPrio);
	OutputMessage.setCanonical(1);
	OutputMessage.setAppId(3);
	OutputMessage.setGocbRef("PC104");
	OutputMessage.setTimeAllowedToLive(TimetoLive);
	OutputMessage.setDatSet("Dataset1");
	OutputMessage.setGoId("PC104");

	OutputMessage.setTest(false);
	OutputMessage.setConfRev(1);
	OutputMessage.setNdsCom(false);

}


void GooseAcquisition::Join_Channel(Channel<digital> *Ch,bool direction, unsigned channel){
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
void GooseAcquisition::Join_Channel(Channel<analogic> *Ch,bool direction, unsigned channel,analogic Min, analogic Max, analogic Dead){
	if(!Ch->InUse()){
		Ch->UseIt();
		if(direction==OUTPUT){
			AnalogicOut.push_back(Ch);
			AnalogicOutBoardChannel.push_back(channel);
			AnalogicOutRatio.push_back(1.0);
			AnalogicOutOffset.push_back(0.0);
			AnalogicOutputDeadBand.push_back(abs((Max-Min)*(Dead/100)));
			}
		else{
			AnalogicIn.push_back(Ch);
			AnalogicInBoardChannel.push_back(channel);
			AnalogicInRatio.push_back(1.0);
			AnalogicInOffset.push_back(0.0);
			}
		}
	else{
		cout<<"Acquisition::Error:: Channel" <<Ch->get_Name()<<" Already in Use."<<endl;
		exit(0);
		}
	}

bool GooseAcquisition::ReciveMessage(){

	unsigned char dmac[6]={0x01,0x0C,0xCD,0x01,0x01,0xFF};
	unsigned len=0;
	goose=htons(GOOSE);
	len = rt_dev_recv(Sock, Packet, ETH_FRAME_LEN, MSG_DONTWAIT);
	if ((len > 0)&&(len<MAX_FRAME_SIZE)){
//		if(!memcmp(EtherHeader->ether_dhost, BroadcastMac, MAC_ADDR_LEN))
		if(1)
		{
		if(!memcmp(EtherHeader->ether_dhost, dmac, MAC_ADDR_LEN)){								
			if(!memcmp(Packet+ETHERTYPE_LEN+(2*MAC_ADDR_LEN)+2,&goose,2)){
				InputMessage.fromNetwork(Packet,len);

				if(!strcmp(InputMessage.getGoId(),IED))
					return true;	
				}
			}
		}
	}
	//if((MasterClock->get_Value()-InputMessage.getTime())>(10*InputMessage.getTimeAllowedToLive()*1E6))
		//cout<<"Acquisition::Goose::No Goose Packets for" << (MasterClock->get_Value()-InputMessage.getTime())/1E9<<" seconds."<<endl;
	
	return false;
}
bool GooseAcquisition::SendMessage(){
	unsigned len=0,sentBytes=0;
	if(OutputMessage.getSqNum()==0)
		Delay=MIN_PUBLISHING_TIME;

	if((MasterClock->get_Value()-LastMessageSentTime)>=Delay){
		LastMessageSentTime=MasterClock->get_Value();
		OutputMessage.setTime(MasterClock->get_Value());
		OutputMessage.toNetwork(Packet,&len);
//		OutputMessage.printMessage();
		sentBytes = rt_dev_send(Sock, Packet, len, 0);
		if (sentBytes < len){
			cout<<"Acquisition::Goose::Error Sending Packet"<<endl;
			return false;
			}
		OutputMessage.incSqNum();
		if((1E-6*Delay)<OutputMessage.getTimeAllowedToLive())
			Delay<<=1;
		else{
			Delay=1E6*OutputMessage.getTimeAllowedToLive();
			}
		}
	return true;
}
bool GooseAcquisition::Prepare(float SFreq){
	unsigned i=0;
	if((DigitalIn.size()>0)||(AnalogicIn.size()>0)){
		while((!ReciveMessage())&&(i<INIT_TRIES)){
			rt_sleep(100000);
			i++;
			}
		if(i>=INIT_TRIES){
			cout<<"Acquisition::Goose::Network not Ready"<<endl;
			return false;
			}
		for(i=0;i<DigitalIn.size();i++)
			if(InputMessage.getDataType(DigitalInBoardChannel[i])!=BOOL_DATA_TYPE){
				cout<<"Acquisition::Goose::Boolean Channel Solicited: "<<DigitalInBoardChannel[i]<<" not Available."<<endl;
				return false;
				}
		for(i=0;i<AnalogicIn.size();i++)
			if(InputMessage.getDataType(AnalogicInBoardChannel[i])!=FLOAT_DATA_TYPE){
				cout<<"Acquisition::Goose::Float Channel Solicited: "<<AnalogicInBoardChannel[i]<<" not Available."<<endl;
				return false;
				}


	}
	for(unsigned int i=0;i<DigitalOut.size();i++)
		OutputMessage.addBooleanData();
	for(unsigned int i=0;i<AnalogicOut.size();i++)
		OutputMessage.addFloatData();

	cout<<"Acquisition::Goose::Prepare Ok."<<endl;
	return true;
}

bool GooseAcquisition::Run(){
		ReciveMessage();
		for(unsigned int i=0;i<DigitalIn.size();i++){
			DigitalIn[i]->insert_Value(InputMessage.getBooleanData(DigitalInBoardChannel[i]));
			}
		for(unsigned int i=0;i<AnalogicIn.size();i++){
			AnalogicIn[i]->insert_Value(AnalogicInRatio[i]*InputMessage.getFloatData(AnalogicInBoardChannel[i])+AnalogicInOffset[i]);
			}
	for(unsigned int i=0;i<DigitalOut.size();i++){
		OutputMessage.setBooleanData(i,DigitalOut[i]->get_Value());
		if(DigitalOut[i]->get_Value()!=DigitalOut[i]->get_Value(1)){
			OutputMessage.resetSqNum();
			OutputMessage.incStNum();
			}
		}
	for(unsigned int i=0;i<AnalogicOut.size();i++){
		OutputMessage.setFloatData(DigitalOut.size()+i,AnalogicOut[i]->get_Value()/AnalogicOutRatio[i]+AnalogicOutOffset[i]);
		if(abs(AnalogicOut[i]->get_Value()-AnalogicOut[i]->get_Value(1))>=AnalogicOutputDeadBand[i]){
			OutputMessage.resetSqNum();
			OutputMessage.incStNum();
			}

		}
	if((DigitalOut.size()>0)||(AnalogicOut.size()>0))
		SendMessage();
	return true;
}

void GooseAcquisition::RefreshTime(){
	MasterClock->insert_Value(orelay_gettime());
}

int rt_eth_aton(unsigned char *addr_buf, const char *mac)
{
    int i = 0;
    int nibble;


    while (1) {
        if (*mac == 0)
            return -EINVAL;

        if ((nibble = hex2int(*mac++)) < 0)
            return nibble;
        *addr_buf = nibble << 4;

        if (*mac == 0)
            return -EINVAL;

        if ((nibble = hex2int(*mac++)) < 0)
            return nibble;
        *addr_buf++ |= nibble;

        if (++i == 6)
            break;

        if ((*mac == 0) || (*mac++ != ':'))
            return -EINVAL;

    }
    return 0;
}


int hex2int(unsigned char hex_char)
{
    if ((hex_char >= '0') && (hex_char <= '9'))
        return hex_char - '0';
    else if ((hex_char >= 'a') && (hex_char <= 'f'))
        return hex_char - 'a' + 10;
    else if ((hex_char >= 'A') && (hex_char <= 'F'))
        return hex_char - 'A' + 10;
    else
        return -EINVAL;
}

#endif
#endif

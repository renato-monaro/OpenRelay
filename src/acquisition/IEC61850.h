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
#ifndef IEC61850_H
#define IEC61850_H

#ifdef IEC61850
#include <Subscriber.h> 
#include <Publisher.h> 
#include "channel.h"
#include "types.h"
#include "parameters.h"
#include "acquisition.h"


namespace orelay{
using namespace Open61850;



class IECSubscriber : public Acquisition, public Subscriber
{
    public:
    IECSubscriber(unsigned type,string Goc, string dev);
	void RefreshTime();
	bool Prepare(float);
	/*IECPublisher(unsigned T, string Dev);
	IECPublisher(string Dev);*/
	void Join_Channel(Channel<digital>*, unsigned channel);
	void Join_Channel(Channel<analogic>*, unsigned channel);
    	bool Run();
};

class IECPublisher : public Acquisition, public Publisher
{
    public:
    IECPublisher(MacAddress destMac, MacAddress localMac, unsigned short vlanid, unsigned char canonical, unsigned char vlanpri, unsigned short appid, long TMin, long TMax, string gocbref, string datset, string Goid, string Dev);
	void RefreshTime();
	bool Prepare(float);
	/*IECPublisher(unsigned T, string Dev);
	IECPublisher(string Dev);*/
	void Join_Channel(Channel<digital>*);
	void Join_Channel(Channel<analogic>*);
    bool Run();
	private:
	unsigned channel_cnt;
};
}
#endif /* !IEC61850_H */
#endif




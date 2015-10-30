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
#ifndef ACQUISITIONGOOSE_H
#define ACQUISITIONGOOSE_H

#ifdef RTAI
#ifdef RTNET
#include <stdio.h>
#include <stdlib.h>
#include <signal.h>
#include <fcntl.h>
#include <sched.h>
#include <sys/types.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <netdb.h>
#include <arpa/inet.h>
#include <netinet/in.h>

#include <net/if.h>

#include <linux/net.h>
#include <linux/if_packet.h>
#include <netinet/if_ether.h>

#include <string>
#include <fstream>

#include <sys/ioctl.h>


#include <rtnet.h>



#include "channel.h"
#include "types.h"
#include "parameters.h"
#include "acquisition.h"
#include "GooseMessage.h"

#define GOOSE 			0x88B8
#define MAC_ADDR_LEN	6
#define ETHERTYPE_LEN	2
#define INIT_TRIES	1E9
#define MAX_FRAME_SIZE	1521
#define MIN_PUBLISHING_TIME	4E6 /*ns*/

namespace orelay{
class GooseAcquisition: public Acquisition{
	public:
		GooseAcquisition(char *OutMac,char *IED, unsigned short VId, unsigned short VPrio, unsigned short TimetoLive, char *Device);
		virtual bool Run();
		bool Prepare(float); 
		void RefreshTime();
		void Join_Channel(Channel<digital>*,bool direction, unsigned channel); /*!< Método que adiciona um canal digital à lista de execução da aquisição.*/
		void Join_Channel(Channel<analogic>*,bool direction, unsigned channel, analogic Min, analogic Max, analogic Dead); /*!< Método que adiciona um canal digital à lista de execução da aquisição.*/
	protected:
		bool ReciveMessage();
		bool SendMessage();
		unsigned char *BroadcastMac;
		char IED[30];
		int Sock;
		GooseMessage InputMessage,OutputMessage;
		unsigned char *Packet;
		struct ether_header *EtherHeader;
		vector<analogic> AnalogicOutputDeadBand;
	private:
		unsigned goose;
		timer Delay;
		timer LastMessageSentTime;
	};

int rt_eth_aton(unsigned char *addr_buf, const char *mac);
int hex2int(unsigned char hex_char);
}
#endif
#endif
#endif /* !ACQUISITIONGOOSE_H */


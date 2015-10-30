/*
# Copyright (C) 2008-2015 Renato Machado Monaro, Hermes Manoel
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
#include "types.h"


#ifdef RTAI
long long orelay_gettime(void){
//	return rt_get_cpu_time_ns();
	return rt_get_real_time_ns();
	}
#endif
#ifndef RTAI
long long orelay_gettime(void){
	struct timespec x;
	long long r;
	clock_gettime(CLOCK_REALTIME, &x);
	r=(long long)(x.tv_sec*1E9+x.tv_nsec);
	return r;
	}
#endif


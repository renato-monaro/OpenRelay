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
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <libxml/xmlmemory.h>
#include <libxml/parser.h>
#include <libxml/xmlschemas.h>
//#include "relay.h"
using namespace orelay;
int is_valid(const xmlDocPtr doc, const char *schema_filename)
{
	xmlDocPtr schema_doc = xmlReadFile(schema_filename, NULL, XML_PARSE_NONET);

	if (schema_doc == NULL) {
		fprintf(stderr,"!!!ERRO!!!  1\n");
		return -1;
	}

	xmlSchemaParserCtxtPtr parser_ctxt = xmlSchemaNewDocParserCtxt(schema_doc);

	if (parser_ctxt == NULL) {
		fprintf(stderr,"!!!ERRO!!!  2\n");
		xmlFreeDoc(schema_doc);
		return -2;
	}

	xmlSchemaPtr schema = xmlSchemaParse(parser_ctxt);

	if (schema == NULL) {
		fprintf(stderr,"!!!ERRO!!!  3\n");
		xmlSchemaFreeParserCtxt(parser_ctxt);
		xmlFreeDoc(schema_doc);
		return -3;
	}

	xmlSchemaValidCtxtPtr valid_ctxt = xmlSchemaNewValidCtxt(schema);

	if (valid_ctxt == NULL) {
		fprintf(stderr,"!!!ERRO!!!  4\n");
		xmlSchemaFree(schema);
		xmlSchemaFreeParserCtxt(parser_ctxt);
		xmlFreeDoc(schema_doc);
		return -4; 
	}

	int is_valid = (xmlSchemaValidateDoc(valid_ctxt, doc) == 0);
	xmlSchemaFreeValidCtxt(valid_ctxt);
	xmlSchemaFree(schema);
	xmlSchemaFreeParserCtxt(parser_ctxt);
	xmlFreeDoc(schema_doc);
	return is_valid ? 1 : 0;
}

void parseC_RELAY (xmlDocPtr doc) {

	xmlNodePtr cur;
	xmlChar *c_name;
	cur = xmlDocGetRootElement(doc);
	cur = cur->xmlChildrenNode;

	while (cur != NULL) {
		if ((!xmlStrcmp(cur->name, (const xmlChar *)"CHANNELS"))) {
			printf("\n	- CHANNELS\n\n");
			cur = cur->xmlChildrenNode;
			do {
				if ((!xmlStrcmp(cur->name, (const xmlChar *)"ANALOGIC"))) {
					printf("\n		- ANALOGIC\n");
					cur=cur->xmlChildrenNode;
					do{
						if ((!xmlStrcmp(cur->name, (const xmlChar *)"CHANNEL"))) {
							xmlChar *c_name, *c_direction, *c_channel, *c_ratio, *c_offset, *c_low, *c_up;
							c_name = xmlGetProp(cur, (const xmlChar*) "name");
							c_direction = xmlGetProp(cur, (const xmlChar*) "direction");
							c_channel = xmlGetProp(cur, (const xmlChar*) "channel");
							c_ratio = xmlGetProp(cur, (const xmlChar*) "ratio");
							c_offset = xmlGetProp(cur, (const xmlChar*) "offset");
							c_low = xmlGetProp(cur, (const xmlChar*) "low");
							c_up = xmlGetProp(cur, (const xmlChar*) "up");
						if((c_name!=NULL)&&(c_direction=NULL)&&(c_channel=NULL)&&(c_ratio=NULL)&&(c_offset=NULL)&&(c_low=NULL)&&(c_up=NULL)){
							printf("		Construtor 1 (%s)\n", c_name);
							}
						if((c_name!=NULL)&&(c_direction!=NULL)&&(c_channel!=NULL)&&(c_ratio=NULL)&&(c_offset=NULL)&&(c_low=NULL)&&(c_up=NULL)){
							printf("		Construtor 2 (%s, %d, %d)\n", c_name, (unsigned)atoi(c_direction), (unsigned)atoi(c_channel));
							}
						if((c_name!=NULL)&&(c_direction!=NULL)&&(c_channel!=NULL)&&(c_ratio!=NULL)&&(c_offset=NULL)&&(c_low=NULL)&&(c_up=NULL)){
							printf("		Construtor 3 (%s, %d, %d, %.2f)\n", c_name, (unsigned)atoi(c_direction), (unsigned)atoi(c_channel), atof(c_ratio));
							}
						if((c_name!=NULL)&&(c_direction!=NULL)&&(c_channel!=NULL)&&(c_ratio!=NULL)&&(c_offset!=NULL)&&(c_low=NULL)&&(c_up=NULL)){
							printf("		Construtor 4 (%s, %d, %d, %.2f, %.2f)\n", c_name, (unsigned)atoi(c_direction), (unsigned)atoi(c_channel), atof(c_ratio), atof(c_offset));
							}
						if((c_name!=NULL)&&(c_direction!=NULL)&&(c_channel!=NULL)&&(c_ratio!=NULL)&&(c_offset=NULL)&&(c_low!=NULL)&&(c_up!=NULL)){
							printf("		Construtor 5 (%s, %d, %d, %f, %.2lf, %.2lf)\n", c_name, (unsigned)atoi(c_direction), (unsigned)atoi(c_channel), atof(c_ratio), strtod(c_low,NULL), strtod(c_up,NULL));
							}
						if((c_name!=NULL)&&(c_direction!=NULL)&&(c_channel!=NULL)&&(c_ratio!=NULL)&&(c_offset!=NULL)&&(c_low!=NULL)&&(c_up!=NULL)){
							printf("		Construtor 6 (%s, %d, %d, %f, %.2f, %.2lf, %.2lf)\n", c_name, (unsigned)atoi(c_direction), (unsigned)atoi(c_channel), atof(c_ratio), atof(c_offset),strtod(c_low,NULL), strtod(c_up,NULL));
							}
						xmlFree(c_name);
						xmlFree(c_direction);
						xmlFree(c_channel);
						xmlFree(c_ratio);
						xmlFree(c_offset);
						xmlFree(c_low);
						xmlFree(c_up);
						}
						cur = cur->next;
						}while(cur->next!=NULL);
					cur = cur->parent;		
					}
				if ((!xmlStrcmp(cur->name, (const xmlChar *)"DIGITAL"))) {
					printf("\n		- DIGITAL\n");
					cur=cur->xmlChildrenNode;
					do{
						if ((!xmlStrcmp(cur->name, (const xmlChar *)"CHANNEL"))) {
							xmlChar *c_name, *c_direction, *c_channel, *c_ratio, *c_offset, *c_low, *c_up;
							c_name = xmlGetProp(cur, (const xmlChar*) "name");
							c_direction = xmlGetProp(cur, (const xmlChar*) "direction");
							c_channel = xmlGetProp(cur, (const xmlChar*) "channel");
							c_ratio = xmlGetProp(cur, (const xmlChar*) "ratio");
							c_offset = xmlGetProp(cur, (const xmlChar*) "offset");
							c_low = xmlGetProp(cur, (const xmlChar*) "low");
							c_up = xmlGetProp(cur, (const xmlChar*) "up");
						if((c_name!=NULL)&&(c_direction=NULL)&&(c_channel=NULL)&&(c_ratio=NULL)&&(c_offset=NULL)&&(c_low=NULL)&&(c_up=NULL)){
							printf("		Construtor 1 (%s)\n", c_name);
							}
						if((c_name!=NULL)&&(c_direction!=NULL)&&(c_channel!=NULL)&&(c_ratio=NULL)&&(c_offset=NULL)&&(c_low=NULL)&&(c_up=NULL)){
							printf("		Construtor 2 (%s, %d, %d)\n", c_name, (unsigned)atoi(c_direction), (unsigned)atoi(c_channel));
							}
						if((c_name!=NULL)&&(c_direction!=NULL)&&(c_channel!=NULL)&&(c_ratio!=NULL)&&(c_offset=NULL)&&(c_low=NULL)&&(c_up=NULL)){
							printf("		Construtor 3 (%s, %d, %d, %.2f)\n", c_name, (unsigned)atoi(c_direction), (unsigned)atoi(c_channel), atof(c_ratio));
							}
						if((c_name!=NULL)&&(c_direction!=NULL)&&(c_channel!=NULL)&&(c_ratio!=NULL)&&(c_offset!=NULL)&&(c_low=NULL)&&(c_up=NULL)){
							printf("		Construtor 4 (%s, %d, %d, %.2f, %.2f)\n", c_name, (unsigned)atoi(c_direction), (unsigned)atoi(c_channel), atof(c_ratio), atof(c_offset));
							}
						if((c_name!=NULL)&&(c_direction!=NULL)&&(c_channel!=NULL)&&(c_ratio!=NULL)&&(c_offset=NULL)&&(c_low!=NULL)&&(c_up!=NULL)){
							printf("		Construtor 5 (%s, %d, %d, %f, %.2lf, %.2lf)\n", c_name, (unsigned)atoi(c_direction), (unsigned)atoi(c_channel), atof(c_ratio), strtod(c_low,NULL), strtod(c_up,NULL));
							}
						if((c_name!=NULL)&&(c_direction!=NULL)&&(c_channel!=NULL)&&(c_ratio!=NULL)&&(c_offset!=NULL)&&(c_low!=NULL)&&(c_up!=NULL)){
							printf("		Construtor 6 (%s, %d, %d, %f, %.2f, %.2lf, %.2lf)\n", c_name, (unsigned)atoi(c_direction), (unsigned)atoi(c_channel), atof(c_ratio), atof(c_offset),strtod(c_low,NULL), strtod(c_up,NULL));
							}
						xmlFree(c_name);
						xmlFree(c_direction);
						xmlFree(c_channel);
						xmlFree(c_ratio);
						xmlFree(c_offset);
						xmlFree(c_low);
						xmlFree(c_up);
							}
						cur = cur->next;
						}while(cur->next!=NULL);
					cur = cur->parent;		
					}
				if ((!xmlStrcmp(cur->name, (const xmlChar *)"COMPLEX"))) {
					printf("\n		- COMPLEX\n");
					cur=cur->xmlChildrenNode;
					do{
						if ((!xmlStrcmp(cur->name, (const xmlChar *)"CHANNEL"))) {
							xmlChar *c_name, *c_direction, *c_channel, *c_ratio, *c_offset, *c_low, *c_up;
							c_name = xmlGetProp(cur, (const xmlChar*) "name");
							c_direction = xmlGetProp(cur, (const xmlChar*) "direction");
							c_channel = xmlGetProp(cur, (const xmlChar*) "channel");
							c_ratio = xmlGetProp(cur, (const xmlChar*) "ratio");
							c_offset = xmlGetProp(cur, (const xmlChar*) "offset");
							c_low = xmlGetProp(cur, (const xmlChar*) "low");
							c_up = xmlGetProp(cur, (const xmlChar*) "up");
						if((c_name!=NULL)&&(c_direction=NULL)&&(c_channel=NULL)&&(c_ratio=NULL)&&(c_offset=NULL)&&(c_low=NULL)&&(c_up=NULL)){
							printf("		Construtor 1 (%s)\n", c_name);
							}
						if((c_name!=NULL)&&(c_direction!=NULL)&&(c_channel!=NULL)&&(c_ratio=NULL)&&(c_offset=NULL)&&(c_low=NULL)&&(c_up=NULL)){
							printf("		Construtor 2 (%s, %d, %d)\n", c_name, (unsigned)atoi(c_direction), (unsigned)atoi(c_channel));
							}
						if((c_name!=NULL)&&(c_direction!=NULL)&&(c_channel!=NULL)&&(c_ratio!=NULL)&&(c_offset=NULL)&&(c_low=NULL)&&(c_up=NULL)){
							printf("		Construtor 3 (%s, %d, %d, %.2f)\n", c_name, (unsigned)atoi(c_direction), (unsigned)atoi(c_channel), atof(c_ratio));
							}
						if((c_name!=NULL)&&(c_direction!=NULL)&&(c_channel!=NULL)&&(c_ratio!=NULL)&&(c_offset!=NULL)&&(c_low=NULL)&&(c_up=NULL)){
							printf("		Construtor 4 (%s, %d, %d, %.2f, %.2f)\n", c_name, (unsigned)atoi(c_direction), (unsigned)atoi(c_channel), atof(c_ratio), atof(c_offset));
							}
						if((c_name!=NULL)&&(c_direction!=NULL)&&(c_channel!=NULL)&&(c_ratio!=NULL)&&(c_offset=NULL)&&(c_low!=NULL)&&(c_up!=NULL)){
							printf("		Construtor 5 (%s, %d, %d, %f, %.2lf, %.2lf)\n", c_name, (unsigned)atoi(c_direction), (unsigned)atoi(c_channel), atof(c_ratio), strtod(c_low,NULL), strtod(c_up,NULL));
							}
						if((c_name!=NULL)&&(c_direction!=NULL)&&(c_channel!=NULL)&&(c_ratio!=NULL)&&(c_offset!=NULL)&&(c_low!=NULL)&&(c_up!=NULL)){
							printf("		Construtor 6 (%s, %d, %d, %f, %.2f, %.2lf, %.2lf)\n", c_name, (unsigned)atoi(c_direction), (unsigned)atoi(c_channel), atof(c_ratio), atof(c_offset),strtod(c_low,NULL), strtod(c_up,NULL));
							}
						xmlFree(c_name);
						xmlFree(c_direction);
						xmlFree(c_channel);
						xmlFree(c_ratio);
						xmlFree(c_offset);
						xmlFree(c_low);
						xmlFree(c_up);
							}
						cur = cur->next;
						}while(cur->next!=NULL);
					cur = cur->parent;		
					}
				cur = cur->next;
				}while(cur->next!=NULL);
			cur = cur->parent;		
			}

		if ((!xmlStrcmp(cur->name, (const xmlChar *)"ACQUISITION"))) {
				printf("\n	- ACQUISITION\n\n");
				cur = cur->xmlChildrenNode;
				do{
					if ((!xmlStrcmp(cur->name, (const xmlChar *)"FILE_ACQUISITION"))) {
						xmlChar *c_filename,*c_bits;
						printf("\n		- FILE_ACQUISITION\n\n");
						c_filename = xmlGetProp(cur, "filename");
						c_bits = xmlGetProp(cur, "bits");
						if((c_bits!=NULL)&&(c_filename!=NULL)){
				    			printf("		Construtor 1 (%s,%d)\n", c_filename,(unsigned)atoi(c_bits));
							}
						if((c_bits==NULL)&&(c_filename!=NULL)){
				    			printf("		Construtor 2 (%s)\n", c_filename);
							}						
						cur = cur->xmlChildrenNode;
							do{
								if ((!xmlStrcmp(cur->name, (const xmlChar *)"CHANNEL"))) {
									c_name = xmlGetProp(cur, "name");
			    						printf("		Channel_Name: %s\n", c_name);
									xmlFree(c_name);
									}
								cur = cur->next;
							}while(cur->next!=NULL);
						cur = cur->parent;
						xmlFree(c_filename);
						xmlFree(c_bits);
						}
					if ((!xmlStrcmp(cur->name, (const xmlChar *)"HARDWARE_ACQUISITION"))) {
						printf("\n		- HARDWARE_ACQUISITION\n\n");
						c_name = xmlGetProp(cur, "device");
							if(c_name!=NULL){
				    				printf("		Construtor 1 (%s)\n", c_name);
								}
			    			printf("		Device: %s\n", c_name);
						xmlFree(c_name);
						cur = cur->xmlChildrenNode;
							do{
								if ((!xmlStrcmp(cur->name, (const xmlChar *)"CHANNEL"))) {
									c_name = xmlGetProp(cur, "name");
			    						printf("		Channel_Name: %s\n", c_name);
									xmlFree(c_name);
									}
								cur = cur->next;
							}while(cur->next!=NULL);
						cur = cur->parent;
						}
				cur = cur->next;
				}while(cur->next!=NULL);
			cur = cur->parent;
			}

		if ((!xmlStrcmp(cur->name, (const xmlChar *)"MEASURE"))) {
			printf("\n	- MEASURE\n\n");
			cur = cur->xmlChildrenNode;
			do {
				if ((!xmlStrcmp(cur->name, (const xmlChar *)"ABS"))) {
					printf("\n		- ABS\n");
					xmlChar *c_ch_in, *c_ch_out;
					c_ch_in = xmlGetProp(cur, "ch_in");
					c_ch_out = xmlGetProp(cur, "ch_out");
						if((c_ch_in!=NULL)&&(c_ch_out!=NULL)){
							printf("		Construtor 1 (%s, %s)\n", c_ch_in, c_ch_out);
							}
					xmlFree(c_ch_in);
					xmlFree(c_ch_out);
					}
				if ((!xmlStrcmp(cur->name, (const xmlChar *)"DFT"))) {
					printf("\n		- DFT\n");
					xmlChar *c_ch_in, *c_ch_out, *c_freq, *c_harm;
					c_ch_in = xmlGetProp(cur, "ch_in");
					c_ch_out = xmlGetProp(cur, "ch_out");
					c_freq = xmlGetProp(cur, "freq");
					c_harm = xmlGetProp(cur, "harm");
						if((c_ch_in!=NULL)&&(c_ch_out!=NULL)&&(c_freq!=NULL)&&(c_harm!=NULL)){
							printf("		Construtor 1 (%s, %s, %.2f, %d)\n", c_ch_in, c_ch_out, atof(c_freq), (unsigned)atoi(c_harm));
							}
						if((c_ch_in!=NULL)&&(c_ch_out!=NULL)&&(c_freq!=NULL)&&(c_harm=NULL)){
							printf("		Construtor 2 (%s, %s, %.2f)\n", c_ch_in, c_ch_out, atof(c_freq));
							}
					xmlFree(c_ch_in);
					xmlFree(c_ch_out);
					xmlFree(c_freq);
					xmlFree(c_harm);	
					}
				if ((!xmlStrcmp(cur->name, (const xmlChar *)"FREQUENCY"))) {
					printf("\n		- FREQUENCY\n");
					xmlChar *c_ch_in, *c_ch_out, *c_met;
					c_ch_in = xmlGetProp(cur, "ch_in");
					c_ch_out = xmlGetProp(cur, "ch_out");
					c_met = xmlGetProp(cur, "met");
						if((c_ch_in!=NULL)&&(c_ch_out!=NULL)&&(c_met!=NULL)){
							printf("		Construtor 1 (%s, %s, %d)\n", c_ch_in, c_ch_out, atoi(c_met));
							}
					xmlFree(c_ch_in);
					xmlFree(c_ch_out);
					xmlFree(c_met);	
					}
				if ((!xmlStrcmp(cur->name, (const xmlChar *)"MEAN"))) {
					printf("\n		- MEAN\n");
					xmlChar *c_ch_in, *c_ch_out, *c_n;
					c_ch_in = xmlGetProp(cur, "ch_in");
					c_ch_out = xmlGetProp(cur, "ch_out");
					c_n = xmlGetProp(cur, "n_samples");
						if((c_ch_in!=NULL)&&(c_ch_out!=NULL)&&(c_n!=NULL)){
							printf("		Construtor 1 (%s, %s, %d)\n", c_ch_in, c_ch_out, (unsigned)atoi(c_n));
							}
						if((c_ch_in!=NULL)&&(c_ch_out!=NULL)&&(c_n=NULL)){
							printf("		Construtor 2 (%s, %s)\n", c_ch_in, c_ch_out);
							}
					xmlFree(c_ch_in);
					xmlFree(c_ch_out);
					xmlFree(c_n);	
					}
				if ((!xmlStrcmp(cur->name, (const xmlChar *)"POWER"))) {
					printf("\n		- POWER\n");
					xmlChar *c_ch_in_i_1, *c_ch_in_i_2, *c_ch_in_i_3, *c_ch_in_v_1, *c_ch_in_v_2, *c_ch_in_v_3, *c_ch_out_s, *c_ch_out_p, *c_ch_out_q, *c_freq;
					c_ch_in_i_1 = xmlGetProp(cur, "ch_in_i_1");
					c_ch_in_i_2 = xmlGetProp(cur, "ch_in_i_2");
					c_ch_in_i_3 = xmlGetProp(cur, "ch_in_i_3");
					c_ch_in_v_1 = xmlGetProp(cur, "ch_in_v_1");
					c_ch_in_v_2 = xmlGetProp(cur, "ch_in_v_2");
					c_ch_in_v_3 = xmlGetProp(cur, "ch_in_v_3");
					c_ch_out_s = xmlGetProp(cur, "ch_out_s");
					c_ch_out_p = xmlGetProp(cur, "ch_out_p");
					c_ch_out_q = xmlGetProp(cur, "ch_out_q");
					c_freq = xmlGetProp(cur, "freq");
						if((c_ch_in_i_1!=NULL)&&(c_ch_in_i_2=NULL)&&(c_ch_in_i_3=NULL)&&(c_ch_in_v_1!=NULL)&&(c_ch_in_v_2=NULL)&&(c_ch_in_v_3=NULL)&&(c_ch_out_s!=NULL)&&(c_ch_out_p!=NULL)&&(c_ch_out_q!=NULL)&&(c_freq!=NULL)){
							printf("		Construtor 1 (%s, %s, %s, %s, %s, %f)\n", c_ch_in_i_1, c_ch_in_v_1, c_ch_out_s, c_ch_out_p, c_ch_out_q,atof(c_freq));
							}
						if((c_ch_in_i_1!=NULL)&&(c_ch_in_i_2!=NULL)&&(c_ch_in_i_3!=NULL)&&(c_ch_in_v_1!=NULL)&&(c_ch_in_v_2!=NULL)&&(c_ch_in_v_3!=NULL)&&(c_ch_out_s!=NULL)&&(c_ch_out_p!=NULL)&&(c_ch_out_q!=NULL)&&(c_freq!=NULL)){
							printf("		Construtor 1 (%s, %s, %s, %s, %s, %s, %s, %s, %s, %f)\n", c_ch_in_i_1, c_ch_in_i_2, c_ch_in_i_3, c_ch_in_v_1, c_ch_in_v_2, c_ch_in_v_3, c_ch_out_s, c_ch_out_p, c_ch_out_q,atof(c_freq));
							}
					xmlFree(c_ch_in_i_1);
					xmlFree(c_ch_in_i_2);
					xmlFree(c_ch_in_i_3);
					xmlFree(c_ch_in_v_1);
					xmlFree(c_ch_in_v_2);
					xmlFree(c_ch_in_v_3);
					xmlFree(c_ch_out_s);
					xmlFree(c_ch_out_p);
					xmlFree(c_ch_out_q);
					xmlFree(c_freq);
					}
				if ((!xmlStrcmp(cur->name, (const xmlChar *)"RMS"))) {
					printf("\n		- RMS\n");
					xmlChar *c_ch_in, *c_ch_out, *c_freq;
					c_ch_in = xmlGetProp(cur, "ch_in");
					c_ch_out = xmlGetProp(cur, "ch_out");
					c_freq = xmlGetProp(cur, "freq");
						if((c_ch_in!=NULL)&&(c_ch_out!=NULL)&&(c_freq!=NULL)){
							printf("		Construtor 1 (%s, %s, %s)\n", c_ch_in, c_ch_out, c_freq);
							}
						if((c_ch_in!=NULL)&&(c_ch_out!=NULL)&&(c_freq=NULL)){
							printf("		Construtor 2 (%s, %s)\n", c_ch_in, c_ch_out);
							}
					xmlFree(c_ch_in);
					xmlFree(c_ch_out);
					xmlFree(c_freq);
					}
				if ((!xmlStrcmp(cur->name, (const xmlChar *)"SCALE"))) {
					printf("\n		- SCALE\n");
					xmlChar *c_ch_in, *c_ch_out, *c_gain, *c_offset;
					c_ch_in = xmlGetProp(cur, "ch_in");
					c_ch_out = xmlGetProp(cur, "ch_out");
					c_gain = xmlGetProp(cur, "gain");
					c_offset = xmlGetProp(cur, "offset");
						if((c_ch_in!=NULL)&&(c_ch_out!=NULL)&&(c_gain!=NULL)&&(c_offset!=NULL)){
							printf("		Construtor 1 (%s, %s, %.2lf, %.2lf)\n", c_ch_in, c_ch_out, strtod(c_gain,NULL), strtod(c_offset,NULL));
							}
					xmlFree(c_ch_in);
					xmlFree(c_ch_out);
					xmlFree(c_gain);
					xmlFree(c_offset);
					}
				if ((!xmlStrcmp(cur->name, (const xmlChar *)"STRING_COMPOSER"))) {
					printf("\n		- STRING_COMPOSER\n");
					xmlChar *c_str, *c_sep;
					c_str = xmlGetProp(cur, "str");
					c_sep = xmlGetProp(cur, "sep");
						if((c_str!=NULL)&&(c_sep!=NULL)){
							printf("		Construtor 1 (%s, %s)\n", c_str, c_sep);
							}
					xmlFree(c_str);
					xmlFree(c_sep);
					}
				if ((!xmlStrcmp(cur->name, (const xmlChar *)"SUM"))) {
					printf("\n		- SUM\n");
					xmlChar *c_ch_out;
					c_ch_out = xmlGetProp(cur, "ch_out");
						if(c_ch_out!=NULL){
							printf("		Construtor 1 (%s)\n", c_ch_out);
							}
					xmlFree(c_ch_out);
					cur = cur->xmlChildrenNode;
						do{
							if ((!xmlStrcmp(cur->name, (const xmlChar *)"CHANNEL"))) {
								c_name = xmlGetProp(cur, "name");
			    					printf("		Channel_Name: %s\n", c_name);
								xmlFree(c_name);
								}
							cur = cur->next;
						}while(cur->next!=NULL);
					cur = cur->parent;
					}
				if ((!xmlStrcmp(cur->name, (const xmlChar *)"WANG_SUN_FREQUENCY"))) {
					printf("\n		- WANG_SUN_FREQUENCY\n");
					xmlChar *c_ch_in, *c_ch_out, *c_ch_tr, *c_total_length, *c_step;
					c_ch_in = xmlGetProp(cur, "ch_in");
					c_ch_out = xmlGetProp(cur, "ch_out");
					c_ch_tr = xmlGetProp(cur, "ch_tr");
					c_total_length = xmlGetProp(cur, "total_length");
					c_step = xmlGetProp(cur, "step");
						if((c_ch_in!=NULL)&&(c_ch_out!=NULL)&&(c_ch_tr!=NULL)&&(c_total_length!=NULL)&&(c_step!=NULL)){
							printf("		Construtor 1 (%s, %s, %s, %d, %d)\n", c_ch_in, c_ch_out, c_ch_tr, (unsigned)atof(c_total_length), (unsigned)atof(c_step));
							}
					xmlFree(c_ch_in);
					xmlFree(c_ch_out);
					xmlFree(c_ch_tr);
					xmlFree(c_total_length);
					xmlFree(c_step);
					}
				if ((!xmlStrcmp(cur->name, (const xmlChar *)"ZERO_CROSS_DETECTION"))) {
					printf("\n		- ZERO_CROSS_DETECTION\n");
					xmlChar *c_ch_in, *c_ch_out;
					c_ch_in = xmlGetProp(cur, "ch_in");
					c_ch_out = xmlGetProp(cur, "ch_out");
						if((c_ch_in!=NULL)&&(c_ch_out!=NULL)){
							printf("		Construtor 1 (%s, %s)\n", c_ch_in, c_ch_out);
							}
					xmlFree(c_ch_in);
					xmlFree(c_ch_out);
					}
			cur = cur->next;

			}while(cur->next!=NULL);
			cur = cur->parent;
		}
		if ((!xmlStrcmp(cur->name, (const xmlChar *)"PROTECTION"))) {
			printf("\n	- PROTECTION\n\n");
			cur = cur->xmlChildrenNode;
			do {
				if ((!xmlStrcmp(cur->name, (const xmlChar *)"DIFFERENTIAL"))) {
					printf("\n		- DIFFERENTIAL\n");
					xmlChar *c_ch_in_1,*c_ch_in_2, *c_ch_out, *c_in, *c_pick_up_1, *c_base_point_1, *c_slope_1, *c_pick_up_2, *c_base_point_2, *c_slope_2;
					c_ch_in_1 = xmlGetProp(cur, "ch_in_1");
					c_ch_in_2 = xmlGetProp(cur, "ch_in_2");
					c_ch_out = xmlGetProp(cur, "ch_out");
					c_in = xmlGetProp(cur, "in");
					c_pick_up_1 = xmlGetProp(cur, "pick_up_1");
					c_base_point_1 = xmlGetProp(cur, "base_point_1");
					c_slope_1 = xmlGetProp(cur, "slope_1");
					c_pick_up_2 = xmlGetProp(cur, "pick_up_2");
					c_base_point_2 = xmlGetProp(cur, "base_point_2");
					c_slope_2 = xmlGetProp(cur, "slope_2");
						if((c_ch_in_1!=NULL)&&(c_ch_in_2!=NULL)&&(c_ch_out!=NULL)&&(c_in!=NULL)&&(c_pick_up_1!=NULL)&&(c_base_point_1!=NULL)&&(c_slope_1!=NULL)&&(c_pick_up_2!=NULL)&&(c_base_point_2!=NULL)&&(c_slope_2!=NULL)){
							printf("		Construtor 1 (%s, %s, %s, %.2lf, %.2lf, %.2lf, %.2lf, %.2lf, %.2lf)\n", c_ch_in_1, c_ch_in_2, c_ch_out, strtod(c_in,NULL), strtod(c_pick_up_1,NULL), strtod(c_base_point_1,NULL), strtod(c_slope_1,NULL), strtod(c_pick_up_2,NULL), strtod(c_base_point_2,NULL), strtod(c_slope_2,NULL));
							}
					xmlFree(c_ch_in_1);
					xmlFree(c_ch_in_2);
					xmlFree(c_ch_out);
					xmlFree(c_in);
					xmlFree(c_pick_up_1);
					xmlFree(c_base_point_1);
					xmlFree(c_slope_1);
					xmlFree(c_pick_up_2);
					xmlFree(c_base_point_2);
					xmlFree(c_slope_2);
					}
				if ((!xmlStrcmp(cur->name, (const xmlChar *)"DIFFERENTIAL_SPLINE"))) {
					printf("\n		- DIFFERENTIAL_SPLINE\n");
					xmlChar *c_ch_in_1,*c_ch_in_2, *c_ch_out, *c_in, *c_pick_up, *c_break_point_1, *c_slope_1, *c_break_point_2, *c_slope_2;
					c_ch_in_1 = xmlGetProp(cur, "ch_in_1");
					c_ch_in_2 = xmlGetProp(cur, "ch_in_2");
					c_ch_out = xmlGetProp(cur, "ch_out");
					c_in = xmlGetProp(cur, "in");
					c_pick_up = xmlGetProp(cur, "pick_up");
					c_break_point_1 = xmlGetProp(cur, "break_point_1");
					c_slope_1 = xmlGetProp(cur, "slope_1");
					c_break_point_2 = xmlGetProp(cur, "break_point_2");
					c_slope_2 = xmlGetProp(cur, "slope_2");
						if((c_ch_in_1!=NULL)&&(c_ch_in_2!=NULL)&&(c_ch_out!=NULL)&&(c_in!=NULL)&&(c_pick_up!=NULL)&&(c_break_point_1!=NULL)&&(c_slope_1!=NULL)&&(c_break_point_2!=NULL)&&(c_slope_2!=NULL)){
							printf("		Construtor 1 (%s, %s, %s, %.2lf, %.2lf, %.2lf, %.2lf, %.2lf)\n", c_ch_in_1, c_ch_in_2, c_ch_out, strtod(c_in,NULL), strtod(c_pick_up,NULL), strtod(c_break_point_1,NULL), strtod(c_slope_1,NULL), strtod(c_break_point_2,NULL), strtod(c_slope_2,NULL));
							}
					xmlFree(c_ch_in_1);
					xmlFree(c_ch_in_2);
					xmlFree(c_ch_out);
					xmlFree(c_in);
					xmlFree(c_pick_up);
					xmlFree(c_break_point_1);
					xmlFree(c_slope_1);
					xmlFree(c_break_point_2);
					xmlFree(c_slope_2);
					}
				if ((!xmlStrcmp(cur->name, (const xmlChar *)"DEFINITE_TIME_OVER_CURRENT"))) {
					printf("\n		- DEFINITE_TIME_OVER_CURRENT\n");
					xmlChar *c_ch_in, *c_ch_out, *c_pick, *c_dial;
					c_ch_in = xmlGetProp(cur, "ch_in");
					c_ch_out = xmlGetProp(cur, "ch_out");
					c_pick = xmlGetProp(cur, "pick");
					c_dial = xmlGetProp(cur, "dial");
						if((c_ch_in!=NULL)&&(c_ch_out!=NULL)&&(c_pick!=NULL)&&(c_dial!=NULL)){
							printf("		Construtor 1 (%s, %s, %.2lf, %.2lf)\n", c_ch_in, c_ch_out, strtod(c_pick,NULL), strtod(c_dial,NULL));
							}
					xmlFree(c_ch_in);
					xmlFree(c_ch_out);
					xmlFree(c_pick);
					xmlFree(c_dial);
					}
				if ((!xmlStrcmp(cur->name, (const xmlChar *)"DEFINITE_TIME_OVER_VOLTAGE"))) {
					printf("\n		- DEFINITE_TIME_OVER_VOLTAGE\n");
					xmlChar *c_ch_in, *c_ch_out, *c_pick, *c_dial;
					c_ch_in = xmlGetProp(cur, "ch_in");
					c_ch_out = xmlGetProp(cur, "ch_out");
					c_pick = xmlGetProp(cur, "pick");
					c_dial = xmlGetProp(cur, "dial");
						if((c_ch_in!=NULL)&&(c_ch_out!=NULL)&&(c_pick!=NULL)&&(c_dial!=NULL)){
							printf("		Construtor 1 (%s, %s, %.2lf, %.2lf)\n", c_ch_in, c_ch_out, strtod(c_pick,NULL), strtod(c_dial,NULL));
							}
					xmlFree(c_ch_in);
					xmlFree(c_ch_out);
					xmlFree(c_pick);
					xmlFree(c_dial);
					}
				if ((!xmlStrcmp(cur->name, (const xmlChar *)"DEFINITE_TIME_UNDER_CURRENT"))) {
					printf("\n		- DEFINITE_TIME_UNDER_CURRENT\n");
					xmlChar *c_ch_in, *c_ch_out, *c_pick, *c_dial;
					c_ch_in = xmlGetProp(cur, "ch_in");
					c_ch_out = xmlGetProp(cur, "ch_out");
					c_pick = xmlGetProp(cur, "pick");
					c_dial = xmlGetProp(cur, "dial");
						if((c_ch_in!=NULL)&&(c_ch_out!=NULL)&&(c_pick!=NULL)&&(c_dial!=NULL)){
							printf("		Construtor 1 (%s, %s, %.2lf, %.2lf)\n", c_ch_in, c_ch_out, strtod(c_pick,NULL), strtod(c_dial,NULL));
							}
					xmlFree(c_ch_in);
					xmlFree(c_ch_out);
					xmlFree(c_pick);
					xmlFree(c_dial);
					}
				if ((!xmlStrcmp(cur->name, (const xmlChar *)"DEFINITE_TIME_UNDER_VOLTAGE"))) {
					printf("\n		- DEFINITE_TIME_UNDER_VOLTAGE\n");
					xmlChar *c_ch_in, *c_ch_out, *c_pick, *c_dial;
					c_ch_in = xmlGetProp(cur, "ch_in");
					c_ch_out = xmlGetProp(cur, "ch_out");
					c_pick = xmlGetProp(cur, "pick");
					c_dial = xmlGetProp(cur, "dial");
						if((c_ch_in!=NULL)&&(c_ch_out!=NULL)&&(c_pick!=NULL)&&(c_dial!=NULL)){
							printf("		Construtor 1 (%s, %s, %.2lf, %.2lf)\n", c_ch_in, c_ch_out, strtod(c_pick,NULL), strtod(c_dial,NULL));
							}
					xmlFree(c_ch_in);
					xmlFree(c_ch_out);
					xmlFree(c_pick);
					xmlFree(c_dial);
					}
				if ((!xmlStrcmp(cur->name, (const xmlChar *)"INSTANTANEOUS_OVER_CURRENT"))) {
					printf("\n		- INSTANTANEOUS_OVER_CURRENT\n");
					xmlChar *c_ch_in, *c_ch_out, *c_pick;
					c_ch_in = xmlGetProp(cur, "ch_in");
					c_ch_out = xmlGetProp(cur, "ch_out");
					c_pick = xmlGetProp(cur, "pick");
						if((c_ch_in!=NULL)&&(c_ch_out!=NULL)&&(c_pick!=NULL)){
							printf("		Construtor 1 (%s, %s, %.2lf)\n", c_ch_in, c_ch_out, strtod(c_pick,NULL));
							}
					xmlFree(c_ch_in);
					xmlFree(c_ch_out);
					xmlFree(c_pick);
					}
				if ((!xmlStrcmp(cur->name, (const xmlChar *)"INSTANTANEOUS_OVER_VOLTAGE"))) {
					printf("\n		- INSTANTANEOUS_OVER_VOLTAGE\n");
					xmlChar *c_ch_in, *c_ch_out, *c_pick;
					c_ch_in = xmlGetProp(cur, "ch_in");
					c_ch_out = xmlGetProp(cur, "ch_out");
					c_pick = xmlGetProp(cur, "pick");
						if((c_ch_in!=NULL)&&(c_ch_out!=NULL)&&(c_pick!=NULL)){
							printf("		Construtor 1 (%s, %s, %.2lf)\n", c_ch_in, c_ch_out, strtod(c_pick,NULL));
							}
					xmlFree(c_ch_in);
					xmlFree(c_ch_out);
					xmlFree(c_pick);
					}
				if ((!xmlStrcmp(cur->name, (const xmlChar *)"INSTANTANEOUS_UNDER_CURRENT"))) {
					printf("\n		- INSTANTANEOUS_UNDER_CURRENT\n");
					xmlChar *c_ch_in, *c_ch_out, *c_pick;
					c_ch_in = xmlGetProp(cur, "ch_in");
					c_ch_out = xmlGetProp(cur, "ch_out");
					c_pick = xmlGetProp(cur, "pick");
						if((c_ch_in!=NULL)&&(c_ch_out!=NULL)&&(c_pick!=NULL)){
							printf("		Construtor 1 (%s, %s, %.2lf)\n", c_ch_in, c_ch_out, strtod(c_pick,NULL));
							}
					xmlFree(c_ch_in);
					xmlFree(c_ch_out);
					xmlFree(c_pick);
					}
				if ((!xmlStrcmp(cur->name, (const xmlChar *)"INSTANTANEOUS_UNDER_VOLTAGE"))) {
					printf("\n		- INSTANTANEOUS_UNDER_VOLTAGE\n");
					xmlChar *c_ch_in, *c_ch_out, *c_pick;
					c_ch_in = xmlGetProp(cur, "ch_in");
					c_ch_out = xmlGetProp(cur, "ch_out");
					c_pick = xmlGetProp(cur, "pick");
						if((c_ch_in!=NULL)&&(c_ch_out!=NULL)&&(c_pick!=NULL)){
							printf("		Construtor 1 (%s, %s, %.2lf)\n", c_ch_in, c_ch_out, strtod(c_pick,NULL));
							}
					xmlFree(c_ch_in);
					xmlFree(c_ch_out);
					xmlFree(c_pick);
					}
				if ((!xmlStrcmp(cur->name, (const xmlChar *)"INVERSE_TIME_OVER_CURRENT"))) {
					printf("\n		- INVERSE_TIME_OVER_CURRENT\n");
					xmlChar *c_ch_in, *c_ch_out, *c_pick, *c_dial, *c_type, *c_reset;
					c_ch_in = xmlGetProp(cur, "ch_in");
					c_ch_out = xmlGetProp(cur, "ch_out");
					c_pick = xmlGetProp(cur, "pick");
					c_dial = xmlGetProp(cur, "dial");
					c_type = xmlGetProp(cur, "type");
					c_reset = xmlGetProp(cur, "reset");
						if((c_ch_in!=NULL)&&(c_ch_out!=NULL)&&(c_pick!=NULL)&&(c_dial!=NULL)&&(c_type!=NULL)&&(c_reset!=NULL)){
							printf("		Construtor 1 (%s, %s, %.2lf, %.2lf, %d)\n", c_ch_in, c_ch_out, strtod(c_pick,NULL), strtod(c_dial,NULL), (unsigned)atoi(c_type));
							}
						if((c_ch_in!=NULL)&&(c_ch_out!=NULL)&&(c_pick!=NULL)&&(c_dial!=NULL)&&(c_type!=NULL)&&(c_reset=NULL)){
							printf("		Construtor 2 (%s, %s, %.2lf, %.2lf, %d)\n", c_ch_in, c_ch_out, strtod(c_pick,NULL), strtod(c_dial,NULL), (unsigned)atoi(c_type));
							}
					xmlFree(c_ch_in);
					xmlFree(c_ch_out);
					xmlFree(c_pick);
					xmlFree(c_dial);
					xmlFree(c_type);
					xmlFree(c_reset);
					}
				if ((!xmlStrcmp(cur->name, (const xmlChar *)"INVERSE_TIME_OVER_VOLTAGE"))) {
					printf("\n		- INVERSE_TIME_OVER_VOLTAGE\n");
					xmlChar *c_ch_in, *c_ch_out, *c_pick, *c_dial, *c_type;
					c_ch_in = xmlGetProp(cur, "ch_in");
					c_ch_out = xmlGetProp(cur, "ch_out");
					c_pick = xmlGetProp(cur, "pick");
					c_dial = xmlGetProp(cur, "dial");
					c_type = xmlGetProp(cur, "type");
						if((c_ch_in!=NULL)&&(c_ch_out!=NULL)&&(c_pick!=NULL)&&(c_dial!=NULL)&&(c_type!=NULL)){
							printf("		Construtor 1 (%s, %s, %.2lf, %.2lf, %d)\n", c_ch_in, c_ch_out, strtod(c_pick,NULL), strtod(c_dial,NULL), (unsigned)atoi(c_type));
							}
					xmlFree(c_ch_in);
					xmlFree(c_ch_out);
					xmlFree(c_pick);
					xmlFree(c_dial);
					xmlFree(c_type);
					}
				if ((!xmlStrcmp(cur->name, (const xmlChar *)"INVERSE_TIME_UNDER_CURRENT"))) {
					printf("\n		- INVERSE_TIME_UNDER_CURRENT\n");
					xmlChar *c_ch_in, *c_ch_out, *c_pick, *c_dial, *c_type;
					c_ch_in = xmlGetProp(cur, "ch_in");
					c_ch_out = xmlGetProp(cur, "ch_out");
					c_pick = xmlGetProp(cur, "pick");
					c_dial = xmlGetProp(cur, "dial");
					c_type = xmlGetProp(cur, "type");
						if((c_ch_in!=NULL)&&(c_ch_out!=NULL)&&(c_pick!=NULL)&&(c_dial!=NULL)&&(c_type!=NULL)){
							printf("		Construtor 1 (%s, %s, %.2lf, %.2lf, %d)\n", c_ch_in, c_ch_out, strtod(c_pick,NULL), strtod(c_dial,NULL), (unsigned)atoi(c_type));
							}
					xmlFree(c_ch_in);
					xmlFree(c_ch_out);
					xmlFree(c_pick);
					xmlFree(c_dial);
					xmlFree(c_type);
					}
				if ((!xmlStrcmp(cur->name, (const xmlChar *)"INVERSE_TIME_UNDER_VOLTAGE"))) {
					printf("\n		- INVERSE_TIME_UNDER_VOLTAGE\n");
					xmlChar *c_ch_in, *c_ch_out, *c_pick, *c_dial, *c_type;
					c_ch_in = xmlGetProp(cur, "ch_in");
					c_ch_out = xmlGetProp(cur, "ch_out");
					c_pick = xmlGetProp(cur, "pick");
					c_dial = xmlGetProp(cur, "dial");
					c_type = xmlGetProp(cur, "type");
						if((c_ch_in!=NULL)&&(c_ch_out!=NULL)&&(c_pick!=NULL)&&(c_dial!=NULL)&&(c_type!=NULL)){
							printf("		Construtor 1 (%s, %s, %.2lf, %.2lf, %d)\n", c_ch_in, c_ch_out, strtod(c_pick,NULL), strtod(c_dial,NULL), (unsigned)atoi(c_type));
							}
					xmlFree(c_ch_in);
					xmlFree(c_ch_out);
					xmlFree(c_pick);
					xmlFree(c_dial);
					xmlFree(c_type);
					}
			cur = cur->next;

			}while(cur->next!=NULL);
			cur = cur->parent;
	
		}
		if ((!xmlStrcmp(cur->name, (const xmlChar *)"DISPLAY"))) {
			printf("\n	- DISPLAY\n\n");
			cur = cur->xmlChildrenNode;
				do{
					if ((!xmlStrcmp(cur->name, (const xmlChar *)"NSCREEN"))) {
						printf("\n		- NSCREEN\n\n");
						printf("		Construtor 1 ()\n");
						cur = cur->xmlChildrenNode;
							do{
								if ((!xmlStrcmp(cur->name, (const xmlChar *)"CHANNEL"))) {
									c_name = xmlGetProp(cur, "name");
			    						printf("		Channel_Name: %s\n", c_name);
									xmlFree(c_name);
									}
								cur = cur->next;
							}while(cur->next!=NULL);
						cur = cur->parent;
						}
				cur = cur->next;
				}while(cur->next!=NULL);
			cur = cur->parent;
		}
		if ((!xmlStrcmp(cur->name, (const xmlChar *)"OSCILLOGRAPHY"))) {
			printf("\n	- OSCILLOGRAPHY\n\n");
			xmlChar *c_trig, *c_before_cicles, *c_after_cicles, *c_frequency, *c_save_folder;
			c_trig = xmlGetProp(cur, "trig");
			c_before_cicles = xmlGetProp(cur, "before_cicles");
			c_after_cicles = xmlGetProp(cur, "after_cicles");
			c_frequency = xmlGetProp(cur, "frequency");
			c_save_folder = xmlGetProp(cur, "save_folder");
				if((c_trig!=NULL)&&(c_before_cicles!=NULL)&&(c_after_cicles!=NULL)&&(c_frequency!=NULL)&&(c_save_folder!=NULL)){
					printf("		Construtor 1 (%s, %d, %d, %.2lf, %s)\n", c_trig, atoi(c_before_cicles), atoi(c_after_cicles), strtod(c_frequency,NULL), c_save_folder);
					}
			xmlFree(c_trig);
			xmlFree(c_before_cicles);
			xmlFree(c_after_cicles);
			xmlFree(c_frequency);
			xmlFree(c_save_folder);
			cur = cur->xmlChildrenNode;
				do{
					if ((!xmlStrcmp(cur->name, (const xmlChar *)"CHANNEL"))) {
						c_name = xmlGetProp(cur, "name");
			    			printf("		Channel_Name: %s\n", c_name);
						xmlFree(c_name);
						}
					cur = cur->next;
				}while(cur->next!=NULL);
			cur = cur->parent;		
			}
	cur = cur->next;
	}

printf("\n");
return;
}




int main(int argc, char **argv) 
{

	char *docname;
	char *n_schema;
	xmlDocPtr doc;
		
	if (argc <= 2) {
		printf("!!!ERRO!!!  INDICAR ARQUIVO XML\n", argv[0]);
		return(0);
		}

	docname = argv[1];
	n_schema = argv[2];
	doc = xmlParseFile(docname);

	if(is_valid (doc, n_schema)){
		fprintf(stderr,"!!!OPERACAO REALIZADA COM EXITO!!!  \n");
		}
	else
		fprintf(stderr,"!!!ERRO NA VALIDACAO!!!  \n");

	printf("\n- RELAY\n\n");
	parseC_RELAY (doc);
	return (1);
}


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
#include "openrelayxml.h"
#define VALIDATOR "open_relay.xsd"
using namespace orelay;
int main(int argc , char** argv){
	OpenRelayXML relay;
	bool verbose=true;
	if(argc<2){
		cout<<"Modo de Uso: openrelay [-v] file.xml [validator.xsd]"<<endl;
		return 1;
		}
	if(relay.read_xml(argv[1])){
		cout<<"Arquivo XML Carregado"<<endl;
		if(argc>2)
			if(relay.read_xsd(argv[2])){
				cout<<"Arquivo XSD Carregado"<<endl;
				}
			else{
				cout<<"Arquivo XSD Não Carregado"<<endl;	
				return 1;
				}
		else
			if(relay.read_xsd((char*)VALIDATOR)){
				cout<<"Arquivo XSD Carregado"<<endl;
				}
			else{
				cout<<"Arquivo XSD Não Carregado"<<endl;
				return 1;
				}
		if(relay.validate_xml()){
			cout<<"Arquivo Validado"<<endl;
			relay.print();
			}
		else{
			cout<<"Arquivo Inválido"<<endl;
			}
/*			if(verbose){
				cout<<"Conjunto do Relé:"<<endl;
				relay.print();
				}
			if(!relay.build()){
				cout<<"Falha na montagem do relé"<<endl;
				return 1;
				}
			relay.run();
			}
		else{
			cout<<"Arquivo: "<<argv[1]<<"Não está conforme padrão."<<endl;
			return 1;
			}*/
		}
	else{
		cout<<"Arquivo: "<<argv[1]<<" Não encontrado."<<endl;
		return 1;
		}
	return 0;
	}

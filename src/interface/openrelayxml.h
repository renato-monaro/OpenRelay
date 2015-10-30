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
#ifndef OPENRELAYXML_H
#define OPENRELAYXML_H

#include<vector>
#include<string>
#include<relay.h>
#include <libxml++/libxml++.h>
using namespace std;

namespace orelay{
//! Esta classe cria um relé dinamicamente a partir de um arquivo XML .
/** Esta classe cria um relé dinamicamente a partir de um arquivo XML .....*/
class OpenRelayXML{
	public:
		OpenRelayXML(); 
		~OpenRelayXML();
		/*! Executa a Leitura de um arquivo XML.*/
		/** @param file É o caminho do arquivo a ser lido.*/
		/** @return Retorna TRUE em caso de sucesso.*/
		bool read_xml(char* file);
		bool read_xsd(char* xsd_file); /*!< Executa a Leitura de um arquivo XSD a ser utilizado como validador do arquivo XML.*/
		bool validate_xml();
		void print();
		bool build();
		bool run();
	protected:
		void print_channels();
		void print_acquisition();
		void print_measure();
		void print_protection();
		void print_save_data();
		void print_display();
		Relay *Target; /*!< Armazena o relé montado.*/
		xmlpp::DomParser parser;
		xmlpp::SchemaValidator validator;
		xmlpp::Node *node;
		xmlpp::Node::NodeList *nodelist;
};
}
#endif /* !OPENRELAY_H */


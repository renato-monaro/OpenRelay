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
using namespace orelay;
struct indent {
	int depth_;
	indent(int depth): depth_(depth) {};
};

ostream & operator<<(std::ostream & o, indent const & in) {
	for(int i = 0; i != in.depth_; ++i) {
		o << "  ";
		}
	return o;
}

OpenRelayXML::OpenRelayXML(){
	}
OpenRelayXML::~OpenRelayXML(){
	}

bool OpenRelayXML::read_xml(char* file){
	try {
		parser.parse_file(string(file));	
		}
	catch(const exception& error){
		cout <<  error.what() << endl;
		return false;
		}
return true;
}

bool OpenRelayXML::read_xsd(char* file){
	try{
		validator.parse_file(string(file));
		}
	catch(const xmlpp::parse_error& error){
		cout <<  error.what() << endl;
		return false;
		}
return true;
}

bool OpenRelayXML::validate_xml(){
		try{
			validator.validate( parser.get_document() );
			}
		catch( const xmlpp::validity_error& error){
			#ifdef LIBXMLCPP_EXCEPTIONS_ENABLED
			cout << error.what();
			#endif	
			return false;	
			}
return true;
}

void OpenRelayXML::print(){
	try {
		node = parser.get_document()->get_root_node();
		print_channels();
	}	
	catch(const std::exception& error){
		cout << error.what() << endl;
	}
}

void OpenRelayXML::print_channels(){

}

void OpenRelayXML::print_acquisition(){
	try {
		cout<<node->get_name()<<endl;
		nodelist = new xmlpp::Node::NodeList(node->get_children("ACQUISITION"));
		nodelist = new xmlpp::Node::NodeList(nodelist->front()->get_children("FILE_ACQUISITION"));
		cout<<"Tamanho:"<<nodelist->size()<<endl;
		cout<<nodelist->front()->get_name()<<endl;
		cout<<nodelist->back()->get_name()<<endl;
		nodelist = new xmlpp::Node::NodeList(nodelist->front()->get_children("HARDWARE_ACQUISITION"));
		cout<<"Tamanho:"<<nodelist->size()<<endl;
		cout<<nodelist->front()->get_name()<<endl;
		cout<<nodelist->back()->get_name()<<endl;
		nodelist = new xmlpp::Node::NodeList(nodelist->front()->get_children("GOOSE_ACQUISITION"));
		cout<<"Tamanho:"<<nodelist->size()<<endl;
		cout<<nodelist->front()->get_name()<<endl;
		cout<<nodelist->back()->get_name()<<endl;

/*
	while (cur != NULL) {
		while(node. && ){
				cout << "*** ELEMENTO ***" << endl;
				cout << "Nome: " << node.get_name() << endl;
				cout << "Posição na árvore: " << node.get_depth() << endl;

				if(reader.has_attributes())
					{
					cout << "Atributos: " << endl;
					node.move_to_first_attribute();
					do
						{
						cout << indent(depth) << "  " << node.get_name() << ": " << node.get_value() << endl;
					} while(node.move_to_next_attribute());
					node.move_to_element();
				}
				else
					{
					cout << "Não há atributos!" << endl;
				}
			}*/
		}
	catch(const std::exception& error){
		cout << error.what() << endl;
	}
}

void OpenRelayXML::print_measure(){
	try {


		}
	catch(const std::exception& error){
		cout << error.what() << endl;
	}
}

void OpenRelayXML::print_protection(){
	try {


		}
	catch(const std::exception& error){
		cout << error.what() << endl;
	}
}

void OpenRelayXML::print_save_data(){
	try {


		}
	catch(const std::exception& error){
		cout << error.what() << endl;
	}
}

void OpenRelayXML::print_display(){
	try {


		}
	catch(const std::exception& error){
		cout << error.what() << endl;
	}
}

bool OpenRelayXML::build(){
	try {


		}
	catch(const std::exception& error){
		cout << error.what() << endl;
	}
return true;
}

bool OpenRelayXML::run(){
	try {


		}
	catch(const std::exception& error){
		cout << error.what() << endl;
	}
return true;
}




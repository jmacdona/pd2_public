//
// (c)  JAMES T. MACDONALD 2010
// This file is part of the PRODART 2 software
// suite.
//
// Copyright (c) 2010, James T. MacDonald
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
//
// Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
// Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
// Neither the name of the James T. MacDonald nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// For more information please contact by email: j.t.macdonald+prodart@gmail.com
// or visit the sclp (Simple Command Line Parser) project home: http://code.google.com/p/sclp/
//
/*
 * options_store.cpp
 *
 *  Created on: 25 Feb 2010
 *      Author: jmacdona
 */
#include "options_store.h"
using boost::any_cast;
using boost::any;
using std::cout;
using std::endl;
using std::cerr;
using std::string;
using std::stringstream;
using std::vector;
using std::map;
using boost::trim;
using boost::lexical_cast;
using boost::split;
using boost::is_any_of;

namespace PRODART {
namespace UTILS {

typedef vector<string> string_vector;
typedef map<string, string> string_string_map;

//std::string str_example;


/*
template <class T>
option_entry::option_entry(const T& val, const std::string descr){
	this->init();
	value = val;
	this->description = descr;
	set_type_info();
	if (valid_flag == false){
		cerr << "ERROR: option_entry: option type is invalid"
			 << endl;
	}
}


option_entry::option_entry(const char* val, const std::string descr, const bool verbose){
	this->option_entry(std::string(val), descr, verbose);
}
*/

option_entry::option_entry() : valid_flag(false),
		set_flag(false), is_hidden_flag(false), type(opt_invalid), description("") {
	this->init();
}

void option_entry::init(){
	valid_flag = false;
	set_flag = false;
	is_hidden_flag = false;
	type = opt_invalid;
}

bool option_entry::is_hidden() const{
	return this->is_hidden_flag;
}

void option_entry::set_hidden(const bool val){
	this->is_hidden_flag = val;
}

bool option_entry::is_valid() const{
	return valid_flag;
}

bool option_entry::is_set() const{
	return this->set_flag;
}

option_entry::option_type option_entry::get_type() const{
	return type;
}

std::string option_entry::get_description() const{
	return this->description;
}

/*
template <class T>
const T option_entry::get_value() const{
	const int ival = any_cast<T>(value);
	return ival;
}
*/



option_entry::option_type option_entry::get_type_info(const boost::any& val) const{
	if (val.type() == typeid(int)){
		return opt_int;
	}
	else if (val.type() == typeid(long)){
		return opt_long;
	}
	else if (val.type() == typeid(unsigned int)){
		return opt_uint;
	}
	else if (val.type() == typeid(unsigned long)){
		return opt_ulong;
	}
	else if (val.type() == typeid(std::string)){
		return opt_string;
	}
	else if (val.type() == typeid(double)){
		return opt_double;
	}
	else if (val.type() == typeid(bool)){
		return opt_bool;
	}
	else if (val.type() == typeid(char)){
		return opt_char;
	}
	else {
		return opt_invalid;
	}
}

std::string option_entry::type_to_string(option_type this_type){
	if (this_type == opt_int){
		return string("int");
	}
	else if (this_type == opt_long){
		return string("long");
	}
	else if (this_type == opt_uint){
		return string("unsigned int");
	}
	else if (this_type == opt_ulong){
		return string("unsigned long");
	}
	else if (this_type == opt_string){
		return string("string");
	}
	else if (this_type == opt_double){
		return string("double");
	}
	else if (this_type == opt_bool){
		return string("bool");
	}
	else if (this_type == opt_char){
		return string("char");
	}
	else {
		return string("invalid");
	}
}

std::string option_entry::value_to_string() const {
	stringstream sstrm (stringstream::in | stringstream::out);
	if (type == opt_int){
		sstrm << any_cast<int>(value);
	}
	else if (type == opt_long){
		sstrm << any_cast<long>(value);
	}
	else if (type == opt_uint){
		sstrm << any_cast<unsigned int>(value);
	}
	else if (type == opt_ulong){
		sstrm << any_cast<unsigned long>(value);
	}
	else if (type == opt_string){
		sstrm << any_cast<string>(value);
	}
	else if (type == opt_double){
		sstrm << any_cast<double>(value);
	}
	else if (type == opt_bool){
		sstrm << any_cast<bool>(value);
	}
	else if (type == opt_char){
		sstrm << any_cast<char>(value);
	}
	else {
		sstrm << "invalid";
	}
	return sstrm.str();
}

bool option_entry::set_value_auto_cast(const std::string value, const bool verbose ){
	try{
	if (type == opt_int){
		set_value(lexical_cast<int>(value), verbose);
	}
	else if (type == opt_long){
		set_value(lexical_cast<long>(value), verbose);
	}
	else if (type == opt_uint){
		set_value(lexical_cast<unsigned int>(value), verbose);
	}
	else if (type == opt_ulong){
		set_value(lexical_cast<unsigned long>(value), verbose);
	}
	else if (type == opt_string){
		set_value(lexical_cast<string>(value), verbose);
	}
	else if (type == opt_double){
		set_value(lexical_cast<double>(value), verbose);
	}
	else if (type == opt_bool){
		if (boost::to_upper_copy(value).compare("TRUE") == 0 || value.compare("T") == 0){
			set_value(true, verbose);
		}
		else if (boost::to_upper_copy(value).compare("FALSE") == 0 || value.compare("F") == 0){
			set_value(false, verbose);
		}
		else{
			set_value(lexical_cast<bool>(value), verbose);
		}
	}
	else if (type == opt_char){
		set_value(lexical_cast<char>(value), verbose);
	}
	else {
		return false;
	}
	}
	catch (std::exception &e){
		throw options_manager_bad_value_exception();
	}
	return true;
}

std::ostream &operator<<( std::ostream& output,  option_entry::option_type type ){
	output << option_entry::type_to_string(type);
	return output;
}

void option_entry::set_type_info(){

	type = this->get_type_info(value);

	if(type == opt_invalid ){
		valid_flag = false;
	}
	else {
		valid_flag = true;
	}

}

/*
int option_entry::get_int() const{
	const int ival = any_cast<int>(value);
	return ival;
}
*/


options_manager::options_manager(){

	option_entry new_opt(string(""), "catch all for arguments without a flag", true);
	string key("");
	new_opt.set_hidden(true);
	entries[key] = new_opt;



}


option_entry::option_type options_manager::get_option_type(const std::string key) const{
	std::map<std::string, option_entry>::const_iterator iter = entries.find(key);
	if (iter != entries.end()){
		option_entry::option_type type = iter->second.get_type();
		return type;
	}
	else {
		std::cerr << "\nERROR: options_manager: this option key does not exist: '" << key << "'"
				  << std::endl;
		throw options_manager_bad_key_exception();
	}

	return option_entry::opt_invalid;
}

bool options_manager::is_set(const std::string key) const{
	std::map<std::string, option_entry>::const_iterator iter = entries.find(key);
	if (iter != entries.end()){
		bool isset = iter->second.is_set();
		return isset;
	}
	else {
		std::cerr << "\nERROR: options_manager: this option key does not exist: '" << key << "'"
				  << std::endl;
		throw options_manager_bad_key_exception();
	}

	return false;
}

bool options_manager::is_registered(const std::string key) const{
	std::map<std::string, option_entry>::const_iterator iter = entries.find(key);
	if (iter != entries.end()){
		return true;
	}
	else {
		return false;
	}

	return false;
}

void options_manager::hide_option(const std::string key){
	std::map<std::string, option_entry>::iterator iter = entries.find(key);
	if (iter != entries.end()){
		iter->second.set_hidden(true);
		return;
	}
	else {
		std::cerr << "\nERROR: options_manager: this option key does not exist: '" << key << "'"
				  << std::endl;
		throw options_manager_bad_key_exception();
	}

	return;
}

void options_manager::unhide_option(const std::string key){
	std::map<std::string, option_entry>::iterator iter = entries.find(key);
	if (iter != entries.end()){
		iter->second.set_hidden(false);
		return;
	}
	else {
		std::cerr << "\nERROR: options_manager: this option key does not exist: '" << key << "'"
				  << std::endl;
		throw options_manager_bad_key_exception();
	}

	return;
}


std::ostream& options_manager::list_options_full_format(std::ostream& output) const{


	output << cmd << endl;

	std::map<std::string, option_entry>::const_iterator iter;

	output << "\nColumn description:" << endl;
	output << "\toption_key" << "\t|"
		   << "value type" << "\t|"
		   << "current value" << "\t|"
		   << "description" << "\t|"
		   << "\n";

	output << "\nEntries:" << endl;
	for (iter = this->entries.begin(); iter != entries.end(); iter++){
		if (!iter->second.is_hidden()){
			output << "\t--" << iter->first << "\t|"
				<< iter->second.get_type() << "\t|"
				<< iter->second.value_to_string() << "\t|"
				<< iter->second.get_description() << "\t|"
				<< endl;
		}
	}

	output << endl;

	return output;
}

bool options_manager::set_option_value_auto_cast(const std::string key, const std::string value){
	std::map<std::string, option_entry>::iterator iter = entries.find(key);
	bool set_result = false;
	if (iter != entries.end()){
		const option_entry::option_type type = iter->second.get_type();
		try {
			set_result = iter->second.set_value_auto_cast(value);
			if (set_result == false){
				throw options_manager_bad_value_exception();
			}
		}
		catch (std::exception& e){
			std::cerr << "\nERROR: options_manager: setting wrong option value type for option key: '" << key << "'"
					  << " with value type: '" << type << "'"
					  << std::endl;
			throw options_manager_bad_value_exception();
		}

	}
	else {
		std::cerr << "\nERROR: options_manager: this option key does not exist: '" << key << "'"
				  << std::endl;
		throw options_manager_bad_key_exception();
	}
	return set_result;
}

bool options_manager::find_args(std::map<std::string, std::string> &parsed_args, std::vector<std::string> &args){

	bool first_flag_found = false;

	const int argc = args.size();
	for (int i = 0; i < argc; i++){
		if ((args[i].substr(0,2).compare("--") == 0 && args[i].size() > 2)
				|| (args[i].substr(0,1).compare("-") == 0 && args[i].size() == 2 && args[i].compare("--") != 0)){

			first_flag_found = true;

			string key = args[i];
			if (args[i].substr(0,2).compare("--") == 0){
				key.replace(0,2,"");
			}
			else if (args[i].substr(0,1).compare("-") == 0){
				key.replace(0,1,"");
			}

			string value;

			int j = i + 1;
			bool get_out = false;
			while (j < argc && get_out != true){
				if ((args[j].substr(0,1).compare("-") != 0)){
					if (value.size() != 0) value.append("\t");
					value.append(args[j]);
				}
				else {
					get_out = true;
				}
				j++;
			}

			trim(value);
			/*
			if (value.substr(value.size() - 1 , 1).compare("\t") == 0){
				value.replace(value.size() - 1,1, "");
			}

			if (value.size() == 0){
				throw options_manager_no_value_exception();
			}
			*/
			//cout << key << "\t" << value << endl;
			parsed_args[key] = value;

		}
		else if (args[i].compare("-") == 0){
			cerr << "\nERROR: options_manager: can not parse command line with a floating '-'\n" << endl;
			return false;
		}
		else if (args[i].compare("--") == 0){
			cerr << "\nERROR: options_manager: can not parse command line with a floating '--'\n" << endl;
			return false;
		}
		else if (args[i].substr(0,1).compare("-") == 0){
			cerr << "\nERROR: options_manager: short command line options with a single '-' and multiple options are not supported\n" << endl;
			return false;
		}
		else if (!first_flag_found){
			string value;
			string key("");
			first_flag_found = true;

			int j = i + 1;
			bool get_out = false;
			while (j < argc && get_out != true){
				if ((args[j].substr(0,1).compare("-") != 0)){
					if (value.size() != 0) value.append("\t");
					value.append(args[j]);
				}
				else {
					get_out = true;
				}
				j++;
			}

			trim(value);
			/*
			if (value.substr(value.size() - 1 , 1).compare("\t") == 0){
				value.replace(value.size() - 1,1, "");
			}

			if (value.size() == 0){
				throw options_manager_no_value_exception();
			}
			*/
			//cout << key << "\t" << value << endl;
			parsed_args[key] = value;
		}
	}

	return true;
}

bool options_manager::parse_cmd_line(int argc, char *argv[]){

	full_cmd_line.clear();

	string_vector args(argc);
	for (int i = 0; i < argc; i ++){
		//cout << i << "\t" << argv[i] << endl;
		args[i] = string(argv[i]);
		trim(args[i]);
		full_cmd_line.append(args[i]);
		full_cmd_line.append(" ");
	}
	cmd.assign(argv[0]);
	string_string_map parsed_args;
	if(!this->find_args(parsed_args, args)){
		cerr << "\nERROR: options_manager: command line args could not be parsed" << endl;
		return false;
	}
	auto_complete_args(parsed_args, true);

	string_string_map::const_iterator iter;
	for (iter = parsed_args.begin(); iter != parsed_args.end(); iter++){
		try {
			const option_entry::option_type type = this->get_option_type(iter->first);

			if (iter->second.size() == 0 && type != option_entry::opt_bool && iter->first.size() != 0){
				cerr << "ERROR: options_manager: parse_cmd_line: no value given" << endl;
				throw options_manager_no_value_exception();
			}
			else if (iter->second.size() == 0 && type == option_entry::opt_bool){
				this->set_option_value(iter->first, true);
			}
			else if (iter->first.size() == 0 && iter->second.size() == 0){
				// do nothing if it previously had empty value
				if (this->get_option_value<std::string>("").size() != 0){
					this->set_option_value_auto_cast(iter->first, iter->second);
				}
			}
			else {
				this->set_option_value_auto_cast(iter->first, iter->second);
			}
		}
		catch (std::exception &e){
			cerr << "ERROR: options_manager: parse_cmd_line: bad arguments: "
				 << iter->first << " = "
				 << iter->second
				 << endl;
			return false;
		}
	}



	return true;
}

bool options_manager::lazy_match_keys(const std::string key, const std::string search_key) const{

	//cout << "matching: '" << key << "' with '" << search_key << "'" << endl;

	if (key.compare(search_key) == 0){
		//cout << "db: a" << endl;
		return true;
	}

	if (!(key.compare("") == 0 || search_key.compare("") == 0)){
		string_vector key_split, search_split;

		split( key_split, key, is_any_of(":") );
		split( search_split, search_key, is_any_of(":") );

		//cout << "db: b" << endl;


		if (search_split.size() >= key_split.size()){
			return false;
		}

		//cout << "db: c" << endl;

		for (int i = 0; i < (int)search_split.size(); i++){
			//cout << "db: d " << i << endl;
			if (key_split[key_split.size() - i - 1].compare(search_split[search_split.size() - i - 1]) != 0){
				//cout << "db: e " << i << endl;
				return false;
			}
		}
		//cout << "db: f" << endl;
		return true;
	}
	return false;
}

std::string options_manager::auto_complete_arg(const std::string arg, const bool verbose) const{
	std::map<std::string, option_entry>::const_iterator iter;
	int match_count = 0;
	std::string matched_key;
	for (iter = this->entries.begin(); iter != entries.end(); iter++){
		if (this->lazy_match_keys(iter->first, arg)){
			//cout << arg << "\t" << iter->first << endl;
			matched_key = iter->first;
			match_count++;
		}
	}

	if (match_count == 1){
		return matched_key;
	}
	else if (verbose && match_count > 1){
		cerr << "ERROR: can not find unique match for argument: '" << arg << "'" << endl;
		return arg;
	}
	else {
		return arg;
	}
	return arg;
}

void options_manager::auto_complete_args(std::map<std::string, std::string> &parsed_args,
		const bool verbose) const{

	std::map<std::string, std::string> new_parsed_args;

	std::map<std::string, std::string>::iterator iter;
	for (iter = parsed_args.begin(); iter != parsed_args.end(); iter++){
		std::string new_arg = this->auto_complete_arg(iter->first, verbose);
		new_parsed_args[new_arg] = iter->second;

	}
	parsed_args = new_parsed_args;
}

//! strangely the constructor doesn't work with the templates unless these are here or it is inlined:
void test(){
	option_entry opt1((double)1.0, std::string("test"));
	option_entry opt2(std::string("blah blah"), std::string("test"));
	option_entry opt3((int)1, "test");
	option_entry opt4((long)1, std::string("test"));
}

}
}



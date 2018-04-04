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
 * options_store.h
 *
 *  Created on: 25 Feb 2010
 *      Author: jmacdona
 */

#ifndef OPTIONS_STORE_H_
#define OPTIONS_STORE_H_

#include <map>
#include <vector>
#include <string>
#include <boost/any.hpp>
#include <iostream>
#include <typeinfo>
#include <exception>
#include <sstream>

#include<boost/algorithm/string.hpp>
#include<boost/lexical_cast.hpp>
#include<boost/algorithm/string/split.hpp>


namespace PRODART {
namespace UTILS {

//! stores an option entry
class option_entry {



public:

	//! allowed option types
	enum option_type{
		opt_invalid = -1,
		opt_string = 0,
		opt_int = 1,
		opt_long = 2,
		opt_double = 3,
		opt_ulong = 4,
		opt_uint = 5,
		opt_bool = 6,
		opt_char = 7
	};

	//template <class T>
	//option_entry(const T& val, const std::string descr);
	template <class T>
	option_entry(const T& val, const std::string descr, const bool verbose = true);
	//option_entry(const char* val, const std::string descr, const bool verbose = true);

	option_entry();

	bool is_hidden() const;
	void set_hidden(const bool val);
	bool is_valid() const;
	bool is_set() const;
	option_type get_type() const;
	std::string get_description() const;

	template <class T>
	bool set_value(const T val, const bool verbose = true);

	bool set_value_auto_cast(const std::string value, const bool verbose = true);

	template <class T>
	const T get_value() const;

	std::string value_to_string() const;

	static std::string type_to_string(option_type this_type);



private:

	void init();

	option_type get_type_info(const boost::any& ) const;
	void set_type_info();

	boost::any value;
	bool valid_flag , set_flag, is_hidden_flag;
	option_type type;
	std::string description;

};

std::ostream &operator<<( std::ostream&,  option_entry::option_type );


//! stores and manages all option entries
class options_manager {

private:
	std::map<std::string, option_entry> entries;
	std::string cmd;
	std::string full_cmd_line;

	bool lazy_match_keys(const std::string key, const std::string search_key) const;
	std::string auto_complete_arg(const std::string arg, const bool verbose = true) const;
	void auto_complete_args(std::map<std::string, std::string> &parsed_args , const bool verbose = true) const;


	bool find_args(std::map<std::string, std::string> &parsed_args, std::vector<std::string> &args);

public:

	options_manager();

	template <class T>
	bool register_option(const std::string key,
			const T default_value,
			const std::string description,
			const bool hidden = false,
			const bool verbose = true);

	template <class T>
	T get_option_value(const std::string key) const;

	template <class T>
	bool set_option_value(const std::string key, const T new_value);

	bool set_option_value_auto_cast(const std::string key, const std::string value);

	option_entry::option_type get_option_type(const std::string key) const;
	bool is_set(const std::string key) const;
	bool is_registered(const std::string key) const;
	void hide_option(const std::string key);
	void unhide_option(const std::string key);

	std::ostream& list_options_full_format(std::ostream& output) const;

	//! command line parser - seems to be working
	//!TODO there is a problem with parsing negative double values due to '-' sign
	bool parse_cmd_line(int argc, char *argv[]);

	std::string get_full_cmd_line() const{
		return full_cmd_line;
	}

	std::string get_cmd() const{
		return cmd;
	}


};




class options_manager_bad_key_exception: public std::exception{
  virtual const char* what() const throw()
  {
    return "options_manager_bad_key_exception";
  }
};

class options_manager_bad_value_exception: public std::exception{
  virtual const char* what() const throw()
  {
    return "options_manager_bad_value_exception";
  }
};

class options_manager_no_value_exception: public std::exception{
  virtual const char* what() const throw()
  {
    return "options_manager_no_value_exception";
  }
};






/********************inlined functions****************************
 * ***************************************************************
 */



template <class T>
inline bool options_manager::register_option(const std::string key,
		const T default_value,
		const std::string description,
		const bool hidden,
		const bool verbose){

	if (entries.find(key) == entries.end() && key.size() > 0){
		option_entry new_opt(default_value, description, verbose);
		if (new_opt.is_valid()){
			new_opt.set_hidden(hidden);
			entries[key] = new_opt;
			return true;
		}
		else{
			std::cerr << "\noptions_manager: option value type for key: '" << key << "' is not a valid type"
					  << std::endl;
			throw options_manager_bad_value_exception();
		}
	}
	else{
		std::cerr << "\noptions_manager: option key: '" << key << "' already registered"
				  << std::endl;
	}
	return false;
}

template <class T>
inline T options_manager::get_option_value(const std::string key) const{
	std::map<std::string, option_entry>::const_iterator iter = entries.find(key);
	if (iter != entries.end()){
		const option_entry::option_type type = iter->second.get_type();
		try {
			const T value = iter->second.get_value<T>();
			return value;
		}
		catch (std::exception& e){
			std::cerr << "\nERROR: options_manager: requested wrong option value type for option key: '" << key << "'"
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
}

template <class T>
inline bool options_manager::set_option_value(const std::string key, const T new_value){
	std::map<std::string, option_entry>::iterator iter = entries.find(key);
	bool set_result = false;
	if (iter != entries.end()){
		const option_entry::option_type type = iter->second.get_type();
		try {
			set_result = iter->second.set_value(new_value);
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


template <class T>
inline option_entry::option_entry(const T& val, const std::string descr, const bool verbose){
	this->init();
	value = val;
	this->description = descr;
	set_type_info();
	if (valid_flag == false && verbose){
		std::cerr << "ERROR: option_entry: option value type: " <<   typeid(T).name()
			 <<" is not permitted"
			 << std::endl;
	}
}


template <class T>
inline bool option_entry::set_value(const T val, const bool verbose){
	const boost::any new_val = val;
	option_type new_type = get_type_info(new_val);
	if (new_type != type || new_type == opt_invalid){
		if (verbose){
			std::cerr << "ERROR: option_entry: option value type: '" << this->get_type_info(boost::any(val)) //typeid(T).name()
				<<"' does not match previous value type: '" <<  this->type //value.type().name()
				<< "' or is of an non-permitted type"
				<< std::endl;
		}
		return false;
	}
	else {
		value = new_val;
		set_flag = true;
		return true;
	}
}

template <class T>
inline const T option_entry::get_value() const{
	const T ival = boost::any_cast<T>(value);
	return ival;
}


}
}
#endif /* OPTIONS_STORE_H_ */

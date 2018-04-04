//
// (c)  JAMES T. MACDONALD 2010
// This file is part of the PRODART 2 software
// suite and is made available under license.
//
// For more information please contact by email: j.t.macdonald+prodart@gmail.com
//
/*
 * prodart_env.h
 *
 *  Created on: 25 Feb 2010
 *      Author: jmacdona
 */

#ifndef PRODART_ENV_H_
#define PRODART_ENV_H_

// macro for debug messages
#define PRINT_EXPR(a) std::cout << __FILE__ << ":" << __LINE__ << ": EXPR: `" << #a << "` = `" << (a) << "`" << std::endl


#include <map>
#include <vector>
#include <string>
#include <iostream>
#include <exception>
#include <sstream>
#include <boost/shared_ptr.hpp>
#include<boost/algorithm/string.hpp>
#include<boost/lexical_cast.hpp>
#include<boost/algorithm/string/split.hpp>
#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/thread.hpp>
#include <boost/timer.hpp>
#include "utils/options_store.h"
#include "utils/MersenneTwister.h"
#include <boost/progress.hpp>
#include <cstdlib>
#include <fstream>

#include <execinfo.h>
#include <signal.h>

#ifdef PD2_MPI_MODE
#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>
#endif

namespace PRODART {
namespace ENV {


namespace {
boost::mutex options_mutex;
}
	
	template <class T>
	bool register_option(const std::string key,
								const T default_value,
								const std::string description,
								const bool hidden = false,
								const bool verbose = true);

//! singleton class
class prodart_env{

	typedef std::map<boost::thread::id, MTRand::MTRand_shared_ptr> thread_id_MTRAND_ptr_map;

private:
	prodart_env(const prodart_env&);

	static void Init();

	PRODART::UTILS::options_manager options;

	MTRand::uint32 rand_num_incr;
	thread_id_MTRAND_ptr_map rand_gen_store;

	bool command_line_parsed;



protected:
	prodart_env();
	virtual ~prodart_env();

	//! before full command line parse retrieve --database option value
	bool get_database_path_option(int argc, char *argv[]);
	bool parse_paths_dat();
	void register_common_options();

	//! The use of get_random_num_gen() is preferred
	//! return random number generator and increments the seed by one each time - each thread should use a different random number generator to be thread safe
	MTRand::MTRand_shared_ptr get_new_random_num_gen();

	PRODART::UTILS::options_manager& get_options_manager(){
		return this->options;
	}
	const PRODART::UTILS::options_manager& get_options_manager() const{
		return this->options;
	}

public:

	//! should be thread safe
	static prodart_env* Instance();

	bool init(int argc, char *argv[]);

	//! should be thread safe now
	MTRand::MTRand_shared_ptr get_random_num_gen(boost::thread::id tid);

	//! should be thread safe now
	MTRand::MTRand_shared_ptr get_random_num_gen(){
		return this->get_random_num_gen(boost::this_thread::get_id());
	}

	bool is_command_line_parsed() const;


private:

	//thread safe friends
	
	friend  MTRand::MTRand_shared_ptr get_random_num_gen();

	template <class T>
	friend T get_option_value(const std::string key);

	template <class T>
	friend bool register_option(const std::string key,
			const T default_value,
			const std::string description,
			const bool hidden,
			const bool verbose);

	
	template <class T>
	friend bool set_option_value(const std::string key, const T new_value);
	friend bool set_option_value_auto_cast(const std::string key, const std::string value);
	friend PRODART::UTILS::option_entry::option_type get_option_type(const std::string key);
	friend bool is_set(const std::string key) ;
	friend bool is_registered(const std::string key) ;
	friend void hide_option(const std::string key);
	friend void unhide_option(const std::string key);
	friend std::string get_full_cmd_line();
	friend std::ostream& list_options_full_format(std::ostream& output);
	friend bool is_command_line_parsed();
	


};


void backtrace_handler(int sig);
void interrupt_handler(int sig);


// Convenience functions:

inline MTRand::MTRand_shared_ptr get_random_num_gen(){
	return prodart_env::Instance()->get_random_num_gen();
}



template <class T>
inline T get_option_value(const std::string key){
	boost::mutex::scoped_lock lock(options_mutex);
	return prodart_env::Instance()->get_options_manager().get_option_value<T>(key);
}

template <class T>
inline bool register_option(const std::string key,
		const T default_value,
		const std::string description,
		const bool hidden,
		const bool verbose){
	boost::mutex::scoped_lock lock(options_mutex);
	return prodart_env::Instance()->get_options_manager().register_option(key, default_value, description, hidden, verbose);
}

template <class T>
inline bool set_option_value(const std::string key, const T new_value){
	boost::mutex::scoped_lock lock(options_mutex);
	return prodart_env::Instance()->get_options_manager().set_option_value(key, new_value);
}

inline bool set_option_value_auto_cast(const std::string key, const std::string value){
	boost::mutex::scoped_lock lock(options_mutex);
	return prodart_env::Instance()->get_options_manager().set_option_value_auto_cast(key, value);
}

inline PRODART::UTILS::option_entry::option_type get_option_type(const std::string key) {
	boost::mutex::scoped_lock lock(options_mutex);
	return prodart_env::Instance()->get_options_manager().get_option_type(key);
}
inline bool is_set(const std::string key) {
	boost::mutex::scoped_lock lock(options_mutex);
	return prodart_env::Instance()->get_options_manager().is_set(key);
}
inline bool is_registered(const std::string key) {
	boost::mutex::scoped_lock lock(options_mutex);
	return prodart_env::Instance()->get_options_manager().is_registered(key);
}
inline void hide_option(const std::string key){
	boost::mutex::scoped_lock lock(options_mutex);
	prodart_env::Instance()->get_options_manager().hide_option(key);
}
inline void unhide_option(const std::string key){
	boost::mutex::scoped_lock lock(options_mutex);
	prodart_env::Instance()->get_options_manager().unhide_option(key);
}
inline std::string get_full_cmd_line(){
	boost::mutex::scoped_lock lock(options_mutex);
	return prodart_env::Instance()->get_options_manager().get_full_cmd_line();
}
inline std::ostream& list_options_full_format(std::ostream& output){
	boost::mutex::scoped_lock lock(options_mutex);
	return prodart_env::Instance()->get_options_manager().list_options_full_format(output);
}

inline bool is_command_line_parsed(){
	boost::mutex::scoped_lock lock(options_mutex);
	return prodart_env::Instance()->is_command_line_parsed();
}


}
}
#endif /* PRODART_ENV_H_ */

/*
 * mpi_test.cpp
 *
 *  Created on: 27 Sep 2013
 *      Author: jmacdona
 */


#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>
#include <iostream>
#include <iostream>
#include <fstream>
#include "prodart_env/prodart_env.h"

namespace mpi = boost::mpi;
using namespace PRODART::ENV;
using namespace PRODART;


using std::cout;
using std::endl;
using std::cerr;
using std::ios;
using std::ifstream;
using std::ofstream;
using std::istream;
using std::strcpy;
using std::strcat;
using std::strcmp;
using std::pow;
using std::sqrt;
using std::asin;
using std::string;

int main(int argc, char* argv[])
{
	mpi::environment env(argc, argv);
	mpi::communicator world;

	PRODART::ENV::register_option("help",bool(false),"print options");


	if (!prodart_env::Instance()->init(argc, argv)){
		cerr << "initialisation errors" << endl;
		return -1;
	}

	if (world.rank() == 0){


		if (PRODART::ENV::get_option_value<bool>("help") == true){
			PRODART::ENV::list_options_full_format(cout);
		}
	}

	world.barrier();

	MTRand::MTRand_shared_ptr rnd =  ENV::get_random_num_gen();

	world.barrier();


	std::cout << "I am process " << world.rank() << " of " << world.size()
            		<< ". Rand: " << rnd->randInt() << std::endl;

	world.barrier();

	return 0;
}



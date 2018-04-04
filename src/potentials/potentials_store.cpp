/*
 * potentials_store.cpp
 *
 *  Created on: 6 Aug 2010
 *      Author: jmacdona
 */
#include "potentials_store.h"

using namespace std;

namespace PRODART {
namespace POSE {
namespace POTENTIALS {



std::ostream& potentials_energies_map::print_headers(std::ostream& output) const{

	potentials_energies_map::const_iterator iter;

	output << "total_energy";

	for (iter = this->begin(); iter != this->end(); iter++ ){
		output << "\t" << iter->first.get_label();
	}

	output << "\n";

	return output;
}


std::ostream& potentials_energies_map::print_unweighted_components(std::ostream& output) const{

	potentials_energies_map::const_iterator iter;

	output << total_energy;

	for (iter = this->begin(); iter != this->end(); iter++ ){
		output << "\t" << iter->second.energy;
	}

	output << "\n";

	return output;
}

std::ostream& potentials_energies_map::print_weighted_components(std::ostream& output) const{

	potentials_energies_map::const_iterator iter;

	output << total_energy;

	for (iter = this->begin(); iter != this->end(); iter++ ){
		output << "\t" << iter->second.energy * iter->second.weight;
	}

	output << "\n";

	return output;
}

std::ostream& potentials_energies_map::print_weights(std::ostream& output) const{

	potentials_energies_map::const_iterator iter;

	output << 1.0;

	for (iter = this->begin(); iter != this->end(); iter++ ){
		output << "\t" << iter->second.weight;
	}

	output << "\n";

	return output;

}

potentials_energies_map potentials_energies_map::get_diff(const potentials_energies_map other) const{
	potentials_energies_map::const_iterator iter;

	potentials_energies_map diff_map;

	diff_map.total_energy = total_energy - other.total_energy;

	for (iter = this->begin(); iter != this->end(); iter++ ){
		//output << "\t" << iter->second.weight;
		if (other.find(iter->first) != other.end()){
			diff_map[iter->first] = (iter->second - other.find(iter->first)->second);
		}
		else {
			// not in other
			diff_map[iter->first] = (iter->second) - potentials_store(0,0);
		}
	}

	for (iter = other.begin(); iter != other.end(); iter++ ){
		if (diff_map.find(iter->first) == diff_map.end()){
			// not in merged map
			diff_map[iter->first] = potentials_store(0,0) - (iter->second);
		}
	}



	return diff_map;
}


std::ostream& potentials_energies_map::print_unweighted_diff(std::ostream& output, const potentials_energies_map other) const{

	potentials_energies_map::const_iterator iter;

	potentials_energies_map diff_map = get_diff(other);




	//output << diff_map.total_energy << endl;


	for (iter = diff_map.begin(); iter != diff_map.end(); iter++ ){
		if (iter->second.energy != 0 ||  iter->second.weight != 0 ){
			output << "COMPONENT:\t" << iter->first.get_label() << "\t"
					<< iter->second.energy << "\t"
					<< iter->second.weight << "\t"
					<< endl;
		}
	}

	//output << "\n";


	return output;
}


}
}
}



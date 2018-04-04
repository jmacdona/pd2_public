/*
 * potentials_factory.cpp
 *
 *  Created on: 5 Aug 2010
 *      Author: jmacdona
 */
#include "potentials_factory.h"

using namespace boost::filesystem;
using namespace std;

using std::string;
using std::ostream;
using std::istream;

using std::ofstream;

using std::ios;

using std::cout;
using std::endl;

using std::cerr;

using std::vector;
using std::cin;
using std::strcpy;
using std::strcat;
using std::strcmp;
using std::strncmp;
using std::exp;
using std::sqrt;
using boost::trim;
using boost::lexical_cast;
using boost::split;
using boost::is_any_of;

namespace PRODART {
namespace POSE {
namespace POTENTIALS {

namespace {
potentials_factory *instance = NULL;
boost::once_flag once_flag = BOOST_ONCE_INIT;
}


void potentials_factory::Init(){
	if (!instance){
		instance = new potentials_factory;
	}
}

potentials_factory::potentials_factory(){


	cout << "potentials_factory: Initialising..." << endl;

	this->init_pot_vec();
	this->init_presets();

}

potentials_factory::~potentials_factory(){
	delete instance;
	instance = NULL;
}

potentials_factory* potentials_factory::Instance(){

	/*
	if (!instance){
		instance = new potentials_factory;
	}
	*/

	boost::call_once(&potentials_factory::Init, once_flag);
	return instance;
}


void potentials_factory::init_pot_vec(){

	all_potentials.clear();

	all_potentials.push_back(new_bond_potential());
	all_potentials.push_back(new_as_geom_pot());
	all_potentials.push_back(new_dihedral_restraint());
	all_potentials.push_back(new_sec_struct_restraint());
	all_potentials.push_back(new_strand_pairing_rst());
	all_potentials.push_back(new_angle_restraint());
	all_potentials.push_back(new_coord_restraint());
	all_potentials.push_back(new_sec_struct_dih_ang_rst());
	all_potentials.push_back(new_sse_axis_restraint());
	all_potentials.push_back(new_residue_contact_restraint());

	all_potentials.push_back(CA::new_ca_bump());
	all_potentials.push_back(CA::new_ca_bond());
	all_potentials.push_back(CA::new_ca_tau_omega());
	all_potentials.push_back(CA::new_ca_frag_freq());
	all_potentials.push_back(CA::new_ca_hbond());
	all_potentials.push_back(CA::new_ca_radgyr());
	all_potentials.push_back(CA::new_ca_pseudo_cb_pos());
	all_potentials.push_back(CA::new_ca_frag_pac());
	all_potentials.push_back(CA::new_ca_GO());
	all_potentials.push_back(CA::new_ca_density());



	all_potentials.push_back(BB::new_bb_bond());
	all_potentials.push_back(BB::new_bb_bump());
	all_potentials.push_back(BB::new_bb_angle());
	all_potentials.push_back(BB::new_bb_dihedral());
	all_potentials.push_back(BB::new_bb_phi_psi_tether());
	all_potentials.push_back(BB::new_bb_frag3_mq());
	all_potentials.push_back(BB::new_bb_forbidden_phi_psi());
	all_potentials.push_back(BB::new_bb_pp_bias_correction());
	all_potentials.push_back(BB::new_bb_old_hbond());
	all_potentials.push_back(BB::new_bb_14_15lj());
	all_potentials.push_back(BB::new_bb_strict_omega());

	potential_shared_ptr_vector::iterator pot_iter;
	for (pot_iter = all_potentials.begin(); pot_iter != all_potentials.end(); pot_iter++){
		if (!(*pot_iter)->init()){
			std::cerr << "ERROR: potentials_factory: problem initialising potentials" << endl;
			throw potentials_factory_init_exception();
		}
	}


}


potential_shared_ptr potentials_factory::get_potential(const potentials_name& name){
	potential_shared_ptr_vector::iterator pot_iter;
	for (pot_iter = all_potentials.begin(); pot_iter != all_potentials.end(); pot_iter++){
		if ((*pot_iter)->provides(name)){
			if ((*pot_iter)->is_disposable()){
				potential_shared_ptr new_pot((*pot_iter)->get_new_instance());
				new_pot->init();
				return new_pot;
			}
			return (*pot_iter);
		}
	}
	return potential_shared_ptr();
}


potentials_container_shared_ptr potentials_factory::make_potentials_container(const potentials_name_double_map& name_weights_list){

	potentials_container_shared_ptr ret_ptr = new_potentials_container();

	potentials_energies_map& nrg_map = ret_ptr->default_energies_map;

	potentials_name_double_map::const_iterator iter;

	for (iter = name_weights_list.begin(); iter != name_weights_list.end(); iter++){

		potential_shared_ptr pot_ptr = this->get_potential(iter->first);

		if (pot_ptr){
			ret_ptr->add_potential(pot_ptr);
			nrg_map[iter->first].weight = iter->second;
			nrg_map[iter->first].energy = 0;
		}
		else {
			// TODO throw some exception

			cerr << "potentials_factory: potential label not found" << endl;
			throw potentials_factory_init_exception();

		}

	}


	return ret_ptr;
}

potentials_container_shared_ptr potentials_factory::make_preset_potentials_container(const std::string preset_label){

	if (presets.find(preset_label) != presets.end()){

		potentials_name_double_map this_preset = presets[preset_label];
		return make_potentials_container(this_preset);
	}
	else {
		cerr << "\nERROR: potentials_factory: preset " << preset_label << " does not exist\n" << endl;
		throw potentials_factory_init_exception();
		return potentials_container_shared_ptr();
	}

}


void potentials_factory::init_presets(){

	cout << "potentials_factory: loading preset weights" << endl;

	const string weightDirPath = PRODART::ENV::get_option_value<string>("database:path:potentials_weights");
	const path dir(weightDirPath);

	if (exists(dir)){
		directory_iterator end;
		for( directory_iterator iter(dir) ; iter != end ; ++iter )
			if ( is_directory( *iter ) )
			{

			}
			else {
				path wts_path(*iter);
				cout << "potentials_factory: Opening weights:\t" << wts_path << endl;
				load_preset(wts_path);

			}
	}

}

double potentials_factory::get_preset_potentials_weight(const std::string preset_label, const potentials_name& name) const{
	if (presets.find(preset_label) != presets.end()){
		const potentials_name_double_map& pres_map = presets.find(preset_label)->second;
		if (pres_map.find(name) != pres_map.end()){
			return pres_map.find(name)->second;
		}
		else {
			return 0;
		}
	}
	else {
		return 0;
	}
}

bool potentials_factory::load_preset(const boost::filesystem::path& wts_path){

	const std::string wts_stem = wts_path.stem().string();

	std::ifstream input(wts_path.string().c_str(), ios::in);

	if (!input.is_open()){
		return false;
	}


    string lineStr;

    long length, lineNum = 0 ;

    vector<string> SplitVec;

    potentials_name_double_map wts_setting;

    while ( !input.eof() ) {
            getline(input, lineStr);
            string resStr;
            lineNum++;

            length = lineStr.length();

            //cout << endl << lineNum << " " << length << " ";

            if (length > 0) {
                    split( SplitVec, lineStr, is_any_of("\t ") );
                    if ( SplitVec[0].substr(0,1).compare("#") != 0
                     && SplitVec[0].substr(0,1).compare("") != 0 && SplitVec.size() >= 2 ){

                        string potName = SplitVec[0];
                        trim(potName);
                        const potentials_name _label(potName);
                        const double weight = lexical_cast<double>(SplitVec[1]);
                        wts_setting[_label] = weight;

                    }

            }

    }

    presets[wts_stem] = wts_setting;

    return true;
}


}
}
}



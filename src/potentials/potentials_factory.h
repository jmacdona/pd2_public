/*
 * potentials_factory.h
 *
 *  Created on: 5 Aug 2010
 *      Author: jmacdona
 */

#ifndef POTENTIALS_FACTORY_H_
#define POTENTIALS_FACTORY_H_

#include <map>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <exception>
#include <sstream>
#include <boost/shared_ptr.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>

#include <boost/thread.hpp>


#include "prodart_env/prodart_env.h"

#include "potential_interface.h"
#include "potentials_container.h"

#include "bond_potential.h"
#include "as_geom_pot.h"
#include "dihedral_restraint.h"
#include "sec_struct_restraint.h"
#include "strand_pairing_rst.h"
#include "angle_restraint.h"
#include "coord_restraint.h"
#include "sec_struct_dih_ang_rst.h"
#include "sse_axis_restraint.h"
#include "residue_contact_restraint.h"

#include "ca_potentials/ca_bump.h"
#include "ca_potentials/ca_bond.h"
#include "ca_potentials/ca_tau_omega.h"
#include "ca_potentials/ca_frag_freq.h"
#include "ca_potentials/ca_hbond.h"
#include "ca_potentials/ca_radgyr.h"
#include "ca_potentials/ca_pseudo_cb_pos.h"
#include "ca_potentials/ca_frag_pac.h"
#include "ca_potentials/ca_GO.h"
#include "ca_potentials/ca_density.h"

#include "bb_potentials/bb_bond.h"
#include "bb_potentials/bb_bump.h"
#include "bb_potentials/bb_angle.h"
#include "bb_potentials/bb_dihedral.h"
#include "bb_potentials/bb_phi_psi_tether.h"
#include "bb_potentials/bb_frag3_mq.h"
#include "bb_potentials/bb_forbidden_phi_psi.h"
#include "bb_potentials/bb_pp_bias_correction.h"
#include "bb_potentials/bb_old_hbond.h"
#include "bb_potentials/bb_14_15lj.h"
#include "bb_potentials/bb_strict_omega.h"

namespace PRODART {
namespace POSE {
namespace POTENTIALS {


typedef std::map<std::string, double> string_double_map;
typedef std::map<std::string, potentials_name_double_map> string_potentials_name_double_map;


class potentials_factory{


private:
	potentials_factory(const potentials_factory&);
	static void Init();

	potential_shared_ptr_vector all_potentials;
	string_potentials_name_double_map presets;

	void init_pot_vec();

	void init_presets();
	bool load_preset(const boost::filesystem::path&);

protected:
	potentials_factory();
	virtual ~potentials_factory();


public:

	static potentials_factory* Instance();

	potentials_container_shared_ptr make_potentials_container(const potentials_name_double_map& name_weights_list);
	potentials_container_shared_ptr make_preset_potentials_container(const std::string preset_label);

	potential_shared_ptr get_potential(const potentials_name& name);

	double get_preset_potentials_weight(const std::string preset_label, const potentials_name& name) const;



};



class potentials_factory_init_exception: public std::exception{
  virtual const char* what() const throw()
  {
    return "potentials_factory_init_exception";
  }
};


}
}
}


#endif /* POTENTIALS_FACTORY_H_ */

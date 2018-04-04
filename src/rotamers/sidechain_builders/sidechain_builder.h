/*
 * sidechain_builder.h
 *
 *  Created on: 12 Jan 2011
 *      Author: jmacdona
 */

#ifndef SIDECHAIN_BUILDER_H_
#define SIDECHAIN_BUILDER_H_

#include "pose/pose.h"
#include "pose_utils/pose_utils.h"
#include "pose_utils/residue_reconstructor.h"
#include "prodart_env/prodart_env.h"
#include "rotamers/sidechain_rotamer.h"

#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>
#include <cstring>
#include <iomanip>
#include <map>

#include<boost/algorithm/string.hpp>
#include<boost/lexical_cast.hpp>
#include<boost/algorithm/string/split.hpp>

#include <boost/thread.hpp>
#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/shared_ptr.hpp>



namespace PRODART {
namespace ROTAMERS {





//! singleton class
class sidechain_builder{

	typedef std::map<POSE::residue_type, POSE_UTILS::residue_reconstructor_shared_ptr> res_type_resconstructor_map;

private:
	sidechain_builder(const sidechain_builder&);
	static void Init();


protected:
	sidechain_builder();
	virtual ~sidechain_builder();

	res_type_resconstructor_map res_rescontr_map;

	void init_resdefs();
	void add_def(boost::filesystem::path def_path);


public:

	static const sidechain_builder* Instance();

	POSE_UTILS::const_residue_reconstructor_shared_ptr get_reconstructor( POSE::residue_type type) const;




};




}
}










#endif /* SIDECHAIN_BUILDER_H_ */

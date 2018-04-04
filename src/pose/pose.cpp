//
// (c)  JAMES T. MACDONALD 2010 
// This file is part of the PRODART 2 software
// suite and is made available under license.
//
// For more information please contact by email: j.t.macdonald+prodart@gmail.com
//
/*
 * pose.cpp
 *
 *  Created on: Jan 31, 2010
 *      Author: jmacdon
 */

#include "pose.h"



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

using boost::lexical_cast;
using boost::split;
using boost::is_any_of;
using boost::trim;

using std::ios;
using std::setiosflags;
using std::resetiosflags;
using std::setprecision;
using std::setw;

using PRODART::UTILS::vector3d;

typedef std::vector< std::string > string_vector;


namespace PRODART {
namespace POSE {

const int AstralStrArrSize = 10;
const char* const AstralStrArr[ AstralStrArrSize ] = {};
    boost::shared_ptr<pose> new_pose()
    {
        boost::shared_ptr<pose> npose(new pose);
        npose->_this = npose;
        return npose;
    }

    //const prot_backbone_map* const pose::bb_map(prot_backbone_map::Instance());
    //const residue_type_map* const pose::rt_map(residue_type_map::Instance());
    unsigned long pose::instance_count = 0;
    unsigned long pose::get_instance_count()
    {
        return instance_count;
    }

    int pose::get_num_protein_backbone_atom_types()
    {
        return prot_backbone_map::Instance()->get_num_protein_backbone_atoms();
    }

    /*
const_pose_shared_ptr pose::get_instance(unsigned long index){
	return instance_register[index].lock();
}
*/
    pose::pose()
    {
        this->init();
        //instance_register[instance_num] = _this;
    }

    pose::~pose()
    {
        //cout << "\tdeleting instance: " << instance_num << endl;
		//cout << "deleting pose: " << this->get_label() << endl;
        this->clear();
        instance_count--;
    }

    unsigned long pose::get_instance_num() const
    {
        return instance_num;
    }

    int pose::get_chain_count() const
    {
        return this->chain_vec.size();
    }

    chain_shared_ptr pose::get_chain(const int chain_num)
    {
        return this->chain_vec[chain_num];
    }

    const_chain_shared_ptr pose::get_chain(const int chain_num) const
    {
        return this->chain_vec[chain_num];
    }

    chain_shared_ptr pose::get_chain(const char chainID)
    {
        chain_shared_ptr_vector::const_iterator iter;
        for(iter = this->chain_vec.begin();iter != chain_vec.end();iter++){
            if((*iter)->getChainID() == chainID){
                return *iter;
            }
        }

        return chain_shared_ptr();
    }

    const_chain_shared_ptr pose::get_chain(const char chainID) const
    {
        chain_shared_ptr_vector::const_iterator iter;
        for(iter = this->chain_vec.begin();iter != chain_vec.end();iter++){
            if((*iter)->getChainID() == chainID){
                return *iter;
            }
        }

        return chain_shared_ptr();
    }

    int pose::get_residue_count() const
    {
        return residue_vec.size();
    }

    const int pose::get_bb_atom_count() const
    {
        return this->prot_backbone_atom_vec.size();
    }

    residue_shared_ptr pose::get_residue(std::string pdb_residue_num, const char chainID, const char iCode)
    {
        //residue_shared_ptr_vector::const_iterator iter;
        const_chain_shared_ptr the_chain = this->get_chain(chainID);
        if(the_chain){
            std::string searchStr = pdb_residue_num;
            trim(searchStr);
            const int start_index = the_chain->get_first_internal_residue_index();
            const int end_index = the_chain->get_last_internal_residue_index();
            if(start_index >= 0 && end_index >= 0){
                for(int resNum = start_index;resNum <= end_index;resNum++){
                    string resNumStr = this->residue_vec[resNum]->get_pdb_residue_index();
                    trim(resNumStr);
                    if(searchStr.compare(resNumStr) == 0 && residue_vec[resNum]->get_icode() == iCode){
                        return residue_vec[resNum];
                    }
                }

            }

        }

        // return NULL
        return residue_shared_ptr();
    }

    const_residue_shared_ptr pose::get_residue(std::string pdb_residue_num, const char chainID, const char iCode) const
    {
        //residue_shared_ptr_vector::const_iterator iter;
        const_chain_shared_ptr the_chain = this->get_chain(chainID);
        if(the_chain){
            std::string searchStr = pdb_residue_num;
            trim(searchStr);
            const int start_index = the_chain->get_first_internal_residue_index();
            const int end_index = the_chain->get_last_internal_residue_index();
            if(start_index >= 0 && end_index >= 0){
                for(int resNum = start_index;resNum <= end_index;resNum++){
                    string resNumStr = this->residue_vec[resNum]->get_pdb_residue_index();
                    trim(resNumStr);
                    if(searchStr.compare(resNumStr) == 0 && residue_vec[resNum]->get_icode() == iCode){
                        return residue_vec[resNum];
                    }
                }

            }

        }

        // return NULL
        return residue_shared_ptr();
    }

    atom_shared_ptr pose::add_new_atom(const PRODART::UTILS::vector3d & vec, const atom_type type, const int residue_num)
    {
        const int rel_pos = bb_map->get_relative_location(type);
        atom_shared_ptr currAtom;
        if(rel_pos < bb_map->get_num_protein_backbone_atoms()){
            currAtom = prot_backbone_atom_vec[(residue_num * this->bb_map->get_num_protein_backbone_atoms()) + rel_pos];
            currAtom->set_type(type);
            currAtom->set_coords(vec);
        }else{
            currAtom = this->residue_vec[residue_num]->get_sidechain()->get_atom(type);
            if(!currAtom){
                currAtom = this->residue_vec[residue_num]->get_sidechain()->add_new_atom(vec, type);
            }else{
                currAtom->set_coords(vec);
            }
        }

        currAtom->setSet(true);
        currAtom->setActive(true);
        return currAtom;
    }

    atom_shared_ptr pose::get_atom(const atom_type type, const int residue_num)
    {
        const int rel_pos = bb_map->get_relative_location(type);
        if(rel_pos < bb_map->get_num_protein_backbone_atoms()){
            atom_shared_ptr currAtom = prot_backbone_atom_vec[(residue_num * this->bb_map->get_num_protein_backbone_atoms()) + rel_pos];
            if(currAtom->isSet()){
                return currAtom;
            }else{
                return atom_shared_ptr();
            }
        }
        else{
            return this->residue_vec[residue_num]->get_sidechain()->get_atom(type);
        }
    }

    const_atom_shared_ptr pose::get_atom(const atom_type type, const int residue_num) const
    {
        const int rel_pos = bb_map->get_relative_location(type);
        if(rel_pos < bb_map->get_num_protein_backbone_atoms()){
            atom_shared_ptr currAtom = prot_backbone_atom_vec[(residue_num * this->bb_map->get_num_protein_backbone_atoms()) + rel_pos];
            if(currAtom->isSet()){
                return currAtom;
            }else{
                return atom_shared_ptr();
            }
        }
        else{
            return this->residue_vec[residue_num]->get_sidechain()->get_atom(type);
        }
    }

    PRODART::UTILS::vector3d pose::get_atom_coords(const atom_type type, const int residue_num) const
    {
        const_atom_shared_ptr currAtom = this->get_atom(type, residue_num);
        if(currAtom){
            return currAtom->get_coords();
        }else{
            return PRODART::UTILS::vector3d();
        }
    }

    atom_shared_ptr pose::get_atom(const atom_type type, std::string pdb_residue_num, const char chainID, const char iCode)
    {
        const_residue_shared_ptr currRes = get_residue(pdb_residue_num, chainID, iCode);
        if(currRes){
            const int resNum = currRes->get_internal_residue_index();
            return this->get_atom(type, resNum);
        }
        //return NULL
        return atom_shared_ptr();
    }

    const_atom_shared_ptr pose::get_atom(const atom_type type, std::string pdb_residue_num, const char chainID, const char iCode) const
    {
        const_residue_shared_ptr currRes = get_residue(pdb_residue_num, chainID, iCode);
        if(currRes){
            const int resNum = currRes->get_internal_residue_index();
            return this->get_atom(type, resNum);
        }
        //return NULL
        return atom_shared_ptr();
    }

    PRODART::UTILS::vector3d pose::get_atom_coords(const atom_type type, std::string pdb_residue_num, const char chainID, const char iCode) const
    {
        const_atom_shared_ptr currAtom = this->get_atom(type, pdb_residue_num, chainID, iCode);
        if(currAtom){
            return currAtom->get_coords();
        }
        return vector3d();
    }

    void pose::set_atom_coords(const PRODART::UTILS::vector3d & vec, const atom_type type, const int residue_num)
    {
        atom_shared_ptr currAtom = this->get_atom(type, residue_num);
        if(currAtom){
            currAtom->set_coords(vec);
        }
    }

    void pose::set_atom_coords(const PRODART::UTILS::vector3d & vec, const atom_type type, std::string pdb_residue_num, const char chainID, const char iCode)
    {
        atom_shared_ptr currAtom = this->get_atom(type, pdb_residue_num, chainID, iCode);
        if(currAtom){
            currAtom->set_coords(vec);
        }
    }

    atom_shared_ptr4_tuple pose::get_phi_atoms(const int residue_num)
    {
        /*

	(|)
	 |

	C'(-1) - N(0) - Ca(0) - C'(0)

	 */
        //
        const_residue_shared_ptr currRes = this->get_residue(residue_num);
        const_residue_shared_ptr prevRes = currRes->get_prev_residue();
        atom_shared_ptr N_0Ptr = this->get_bb_atom(N, residue_num);
        atom_shared_ptr Ca_0Ptr = this->get_bb_atom(CA, residue_num);
        atom_shared_ptr C_0Ptr = this->get_bb_atom(C, residue_num);
        if(!prevRes || residue_num < 1){
            // no previous atom
            atom_shared_ptr4_tuple tup(atom_shared_ptr(), N_0Ptr, Ca_0Ptr, C_0Ptr);
            return tup;
        }
        atom_shared_ptr C_minus1Ptr = this->get_bb_atom(C, residue_num - 1);
        atom_shared_ptr4_tuple tup(C_minus1Ptr, N_0Ptr, Ca_0Ptr, C_0Ptr);
        return tup;
    }

    atom_shared_ptr4_tuple pose::get_psi_atoms(const int residue_num)
    {
        /*

	\|/
	 |

	 N(0) - Ca(0) - C'(0) - N(+1)

	 */
        const_residue_shared_ptr currRes = this->get_residue(residue_num);
        const_residue_shared_ptr nextRes = currRes->get_next_residue();
        atom_shared_ptr N_0Ptr = this->get_bb_atom(N, residue_num);
        atom_shared_ptr Ca_0Ptr = this->get_bb_atom(CA, residue_num);
        atom_shared_ptr C_0Ptr = this->get_bb_atom(C, residue_num);
        if(!nextRes || residue_num >= this->get_residue_count()){
            // no next res
            atom_shared_ptr4_tuple tup(N_0Ptr, Ca_0Ptr, C_0Ptr, atom_shared_ptr());
            return tup;
        }
        atom_shared_ptr N_plus1Ptr = this->get_bb_atom(N, residue_num + 1);
        atom_shared_ptr4_tuple tup(N_0Ptr, Ca_0Ptr, C_0Ptr, N_plus1Ptr);
        return tup;
        //const double PsiAng = PRODART::UTILS::dihedral(N_0, Ca_0, C_0, N_plus1);
    }

    atom_shared_ptr4_tuple pose::get_omega_atoms(const int residue_num)
    {
        /*

	\|/

	 Ca(-1) - C'(-1) - N(0) - Ca(0)

	 */
        const_residue_shared_ptr currRes = this->get_residue(residue_num);
        const_residue_shared_ptr prevRes = currRes->get_prev_residue();
        atom_shared_ptr N_0Ptr = this->get_bb_atom(N, residue_num);
        atom_shared_ptr Ca_0Ptr = this->get_bb_atom(CA, residue_num);
        if(!prevRes || residue_num < 1){
            //const double OmegaToPreviousAng = PRODART::UTILS::dihedral(Ca_minus1, C_minus1, N_0, Ca_0);
            // no previous atom
            atom_shared_ptr4_tuple tup(atom_shared_ptr(), atom_shared_ptr(), N_0Ptr, Ca_0Ptr);
            return tup;
        }
        atom_shared_ptr Ca_minus1Ptr = this->get_bb_atom(CA, residue_num - 1);
        atom_shared_ptr C_minus1Ptr = this->get_bb_atom(C, residue_num - 1);
        atom_shared_ptr4_tuple tup(Ca_minus1Ptr, C_minus1Ptr, N_0Ptr, Ca_0Ptr);
        return tup;
    }

    double pose::get_phi(const int residue_num) const
    {
        /*

	(|)
	 |

	C'(-1) - N(0) - Ca(0) - C'(0)

	 */
        //
        const_residue_shared_ptr currRes = this->get_residue(residue_num);
        const_residue_shared_ptr prevRes = currRes->get_prev_residue();
        if(!prevRes || residue_num < 1){
            // no previous atom
            return 0;
        }
        vector3d C_minus1, N_0, Ca_0, C_0;
        const_atom_shared_ptr C_minus1Ptr = this->get_bb_atom(C, residue_num - 1);
        const_atom_shared_ptr N_0Ptr = this->get_bb_atom(N, residue_num);
        const_atom_shared_ptr Ca_0Ptr = this->get_bb_atom(CA, residue_num);
        const_atom_shared_ptr C_0Ptr = this->get_bb_atom(C, residue_num);
        if(C_minus1Ptr->isSet() && N_0Ptr->isSet() && Ca_0Ptr->isSet() && C_0Ptr->isSet()){
            C_minus1 = C_minus1Ptr->get_coords();
            N_0 = N_0Ptr->get_coords();
            Ca_0 = Ca_0Ptr->get_coords();
            C_0 = C_0Ptr->get_coords();
        }else{
            return 0;
        }
        const double PhiAng = PRODART::UTILS::dihedral(C_minus1, N_0, Ca_0, C_0);
        return PhiAng;
    }

    double pose::get_psi(const int residue_num) const
    {
        /*

	\|/
	 |

	 N(0) - Ca(0) - C'(0) - N(+1)

	 */
        const_residue_shared_ptr currRes = this->get_residue(residue_num);
        const_residue_shared_ptr nextRes = currRes->get_next_residue();
        if(!nextRes || residue_num >= this->get_residue_count()){
            // no previous atom
            return 0;
        }
        vector3d N_0, Ca_0, C_0, N_plus1;
        const_atom_shared_ptr N_0Ptr = this->get_bb_atom(N, residue_num);
        const_atom_shared_ptr Ca_0Ptr = this->get_bb_atom(CA, residue_num);
        const_atom_shared_ptr C_0Ptr = this->get_bb_atom(C, residue_num);
        const_atom_shared_ptr N_plus1Ptr = this->get_bb_atom(N, residue_num + 1);
        if(N_plus1Ptr->isSet() && N_0Ptr->isSet() && Ca_0Ptr->isSet() && C_0Ptr->isSet()){
            N_0 = N_0Ptr->get_coords();
            Ca_0 = Ca_0Ptr->get_coords();
            C_0 = C_0Ptr->get_coords();
            N_plus1 = N_plus1Ptr->get_coords();
        }else{
            return 0;
        }
        const double PsiAng = PRODART::UTILS::dihedral(N_0, Ca_0, C_0, N_plus1);
        return PsiAng;
    }

    double pose::get_omega_to_prev(const int residue_num) const
    {
        /*

	\|/

	 Ca(-1) - C'(-1) - N(0) - Ca(0)

	 */
        const_residue_shared_ptr currRes = this->get_residue(residue_num);
        const_residue_shared_ptr prevRes = currRes->get_prev_residue();
        if(!prevRes || residue_num < 1){
            // no previous atom
            return 0;
        }
        vector3d Ca_minus1, C_minus1, N_0, Ca_0;
        const_atom_shared_ptr Ca_minus1Ptr = this->get_bb_atom(CA, residue_num - 1);
        const_atom_shared_ptr C_minus1Ptr = this->get_bb_atom(C, residue_num - 1);
        const_atom_shared_ptr N_0Ptr = this->get_bb_atom(N, residue_num);
        const_atom_shared_ptr Ca_0Ptr = this->get_bb_atom(CA, residue_num);
        if(Ca_minus1Ptr->isSet() && N_0Ptr->isSet() && Ca_0Ptr->isSet() && C_minus1Ptr->isSet()){
            Ca_minus1 = Ca_minus1Ptr->get_coords();
            C_minus1 = C_minus1Ptr->get_coords();
            N_0 = N_0Ptr->get_coords();
            Ca_0 = Ca_0Ptr->get_coords();
        }else{
            return 0;
        }
        const double OmegaToPreviousAng = PRODART::UTILS::dihedral(Ca_minus1, C_minus1, N_0, Ca_0);
        return OmegaToPreviousAng;
    }

    double pose::get_omega_to_next(const int residue_num) const
    {
        /*

	\|/

	 Ca(0) - C'(0) - N(+1) - Ca(+1)

	 */
        const_residue_shared_ptr currRes = this->get_residue(residue_num);
        const_residue_shared_ptr nextRes = currRes->get_next_residue();
        if(!nextRes || residue_num >= this->get_residue_count()){
            // no previous atom
            return 0;
        }
        vector3d Ca_0, C_0, N_plus1, Ca_plus1;
        const_atom_shared_ptr Ca_0Ptr = this->get_bb_atom(CA, residue_num);
        const_atom_shared_ptr C_0Ptr = this->get_bb_atom(C, residue_num);
        const_atom_shared_ptr N_plus1Ptr = this->get_bb_atom(N, residue_num + 1);
        const_atom_shared_ptr Ca_plus1Ptr = this->get_bb_atom(CA, residue_num + 1);
        if(Ca_0Ptr != 0 && C_0Ptr != 0 && N_plus1Ptr != 0 && Ca_plus1Ptr != 0){
            Ca_0 = Ca_0Ptr->get_coords();
            C_0 = C_0Ptr->get_coords();
            N_plus1 = N_plus1Ptr->get_coords();
            Ca_plus1 = Ca_plus1Ptr->get_coords();
        }else{
            return 0;
        }
        const double OmegaToNextAng = PRODART::UTILS::dihedral(Ca_0, C_0, N_plus1, Ca_plus1);
        return OmegaToNextAng;
    }

    void pose::add_remark(const std::string & str)
    {
        string space(" ");
        space.append(str);
        this->REMARK_vec.push_back(space);
    }

    void pose::add_remark_with_prodart_label(const std::string & str)
    {
        string space(" PRODART2 ");
        space.append(str);
        this->REMARK_vec.push_back(space);
    }

    void pose::clear_remarks()
    {
        this->REMARK_vec.clear();
    }

    void pose::suppress_remark_output()
    {
        this->suppress_REMARK_output_flag = true;
    }

    void pose::unsuppress_remark_output()
    {
        this->suppress_REMARK_output_flag = false;
    }

    //! add line to append to output
    void pose::add_appendline(std::string str)
    {
        trim(str);
        appendlines_vec.push_back(str);
    }

    //! clear REMARK records
    void pose::clear_appendlines()
    {
        this->appendlines_vec.clear();
    }

    //! suppress REMARK record output
    void pose::suppress_appendlines_output()
    {
        this->suppress_appendlines_output_flag = true;
    }

    //! unsupress REMARK record output
    void pose::unsuppress_appendlines_output()
    {
        this->suppress_appendlines_output_flag = false;
    }

    void pose::init()
    {
        bb_map = prot_backbone_map::Instance();
        rt_map = residue_type_map::Instance();
        prot_backbone_atom_vec.reserve(1000 * 11);
        residue_vec.reserve(1000);
        chain_vec.reserve(10);
        all_atoms_vec.reserve(1000 * 40);
        instance_num = instance_count;
        suppress_REMARK_output_flag = false;
        suppress_appendlines_output_flag = false;
        //CAonly=false;				//MIS
        this->set_version_date_strings();
        this->reset_cryst_record();
        model_number = 1;
        instance_count++;
    }

    void pose::clear()
    {
        //cout << "deleting backbone atoms:" << endl;
        this->prot_backbone_atom_vec.clear();
        prot_backbone_atom_vec.reserve(1000 * 11);
        //cout << "deleting residues:" << endl;
        this->residue_vec.clear();
        residue_vec.reserve(1000);
        //cout << "deleting chains:" << endl;
        this->chain_vec.clear();
        chain_vec.reserve(10);
		//cout << "deleting all remaining atoms:" << endl;
        all_atoms_vec.clear();
        all_atoms_vec.reserve(1000 * 40);
        this->clear_remarks();
        this->reset_cryst_record();
        model_number = 1;
        /*
	all_atom_vec.clear();
	all_atom_vec.reserve(1000*30);
*/
    }

    void pose::reset_cryst_record()
    {
        cryst_a = cryst_b = cryst_c = 1.0;
        cryst_alpha = cryst_beta = cryst_gamma = 90.0;
        cryst_space_group = "P 1";
        cryst_z_value = 1;
    }

    void pose::set_version_date_strings()
    {
        version_string.assign("PRODART2 version: ");
        version_string.append(ENV::get_option_value<string>("prodart_version")); //_PRODART_VERSION_);
        if (ENV::get_option_value<string>("hg_version_number").compare("") != 0){
        	version_string.append(" build: ");
        	version_string.append(ENV::get_option_value<string>("hg_version_number")); //SVN_REV);
        }
        version_string.append(" compile date: ");
        version_string.append(ENV::get_option_value<string>("compile_date")); //_COMPILE_DATE_);
        run_time_string.assign("PRODART2 run date: ");
        time_t rawtime;
        struct tm *timeinfo;
        time(&rawtime);
        timeinfo = localtime(&rawtime);
        string tmp_time_str(asctime(timeinfo));
        trim(tmp_time_str);
        run_time_string.append(tmp_time_str);
    }

    void pose::outputREMARKs(std::ostream & output) const
    {
        if(!suppress_REMARK_output_flag){
            string_vector::const_iterator str_vec_iter;
            for(str_vec_iter = REMARK_vec.begin();str_vec_iter != REMARK_vec.end();str_vec_iter++){
                output << "REMARK" << *str_vec_iter << "\n";
            }
            output << "REMARK " << version_string << "\n" << "REMARK " << run_time_string << "\n";
        }

    }

    void pose::outputCRYST1(std::ostream & output) const
    {
        /*
	Contains:	 unit cell parameters, space group, and Z value
	Notes:
	If the structure was not determined by crystallographic means, simply defines a unit cube
	(a = b =c = 1.0, alpha = beta = gamma = 90
	degrees, space group = P 1, and Z = 1)
	The Hermann-Mauguin space group symbol is given without parenthesis,
	e.g., P 21 21 2 and using the full symbol, e.g., C 1 2 1 instead of C 2.
	The screw axis is described as a two digit number.
	For a rhombohedral space group in the hexagonal setting, the lattice type symbol used is H.
	The Z value is the number of polymeric chains in a unit cell. In the case of heteropolymers,
	Z is the number of occurrences of the most populous chain.
	In the case of a polycrystalline fiber diffraction study,
	CRYST1 and SCALE contain the normal unit cell data.
	The unit cell parameters are used to calculate SCALE.
	COLUMNS       DATA TYPE      CONTENTS
	--------------------------------------------------------------------------------
	 	 1 -  6       Record name    "CRYST1"
	 	 7 - 15       Real(9.3)      a (Angstroms)
	16 - 24       Real(9.3)      b (Angstroms)
	25 - 33       Real(9.3)      c (Angstroms)
	34 - 40       Real(7.2)      alpha (degrees)
	41 - 47       Real(7.2)      beta (degrees)
	48 - 54       Real(7.2)      gamma (degrees)
	56 - 66       LString        Space group
	67 - 70       Integer        Z value
	 */
        output << setw(6) << "CRYST1" << resetiosflags(ios::left) << setiosflags(ios::fixed) << setw(9) << setprecision(3) << cryst_a << setw(9) << setprecision(3) << cryst_b << setw(9) << setprecision(3) << cryst_c << setw(7) << setprecision(2) << cryst_alpha << setw(7) << setprecision(2) << cryst_beta << setw(7) << setprecision(2) << cryst_gamma << " " << setiosflags(ios::left) << setw(11) << cryst_space_group << resetiosflags(ios::left) << setw(4) << cryst_z_value << "\n";
    }

    void pose::output_appendlines(std::ostream & output) const
    {
        if(!suppress_appendlines_output_flag){
            string_vector::const_iterator str_vec_iter;
            for(str_vec_iter = appendlines_vec.begin();str_vec_iter != appendlines_vec.end();str_vec_iter++){
                output << *str_vec_iter << "\n";
            }
        }

    }

    void pose::renumber_residues(const int chain_num, const int start_res_num)
    {
        this->index();
        if(chain_num >= this->get_chain_count()){
            //error out of range
            return;
        }
        const chain_shared_ptr currChain = this->get_chain(chain_num);
        const int start_index = currChain->get_first_internal_residue_index();
        const int end_index = currChain->get_last_internal_residue_index();
        int res_num = start_res_num;
        for(int i = start_index;i <= end_index;i++){
            const residue_shared_ptr currRes = this->get_residue(i);
            currRes->set_pdb_residue_index(res_num);
            res_num++;
        }
    }

    void pose::renumber_residues(const char chainID, const int start_res_num)
    {
        this->index();
        const chain_shared_ptr currChain = this->get_chain(chainID);
        if(currChain){
            const int start_index = currChain->get_first_internal_residue_index();
            const int end_index = currChain->get_last_internal_residue_index();
            int res_num = start_res_num;
            for(int i = start_index;i <= end_index;i++){
                const residue_shared_ptr currRes = this->get_residue(i);
                currRes->set_pdb_residue_index(res_num);
                res_num++;
            }
        }
        else{
            cerr << "pose: ERROR: renumber_residues: can not find chain\t" << chainID << endl;
        }
    }

    void pose::renumber_residues(const int start_res_num)
    {
        const int chain_count = this->get_chain_count();
        for(int i = 0;i < chain_count;i++){
            this->renumber_residues(i, start_res_num);
        }
    }

	//! renumber (external) residue numbers in a chain
	void pose::cyclic_renumber_residues(const int chain_num, const int start_res_num, const int lower_bound, const int upper_bound, const int period){
        this->index();
        if(chain_num >= this->get_chain_count()){
            //error out of range
            return;
        }
        const chain_shared_ptr currChain = this->get_chain(chain_num);
        const int start_index = currChain->get_first_internal_residue_index();
        const int end_index = currChain->get_last_internal_residue_index();
        int res_num = start_res_num;
        for(int i = start_index;i <= end_index;i++){
            const residue_shared_ptr currRes = this->get_residue(i);
            currRes->set_pdb_residue_index(res_num);
            res_num++;
            if (res_num > upper_bound){
            	res_num = lower_bound;
            }
        }
	}
	//! renumber (external) residue numbers for each chain starting from 1
	void pose::cyclic_renumber_residues(const int start_res_num, const int lower_bound, const int upper_bound, const int period){
        const int chain_count = this->get_chain_count();
        for(int i = 0;i < chain_count;i++){
            this->cyclic_renumber_residues(i, start_res_num, lower_bound, upper_bound, period);
        }
	}

	void pose::renumber_chainIDs(){
		const int chain_count = this->get_chain_count();
		std::string alphabet = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789";
		for(int i = 0;i < chain_count;i++){
			if (i < alphabet.size()){
				this->get_chain(i)->setChainID(alphabet[i]);
			}
			else {
				this->get_chain(i)->setChainID(' ');
			}
		}
	}

	void pose::auto_assign_missing_chainIDs(){
		const int chain_count = this->get_chain_count();
		std::string alphabet = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789";
		for(int i = 0;i < chain_count;i++){
			if (this->get_chain(i)->getChainID() == ' '){
				char char_to_use = ' ';
				for (unsigned int jj = 0; jj < alphabet.size(); jj++){
					if(!this->get_chain(alphabet[jj])){
						char_to_use = alphabet[jj];
						break;
					}
				}
				this->get_chain(i)->setChainID(char_to_use);
			}
		}
	}

    std::istream & pose::loadPdb(std::istream & input, const bool check_prev_atomid)
    {
        this->clear();
        //pose_weak_ptr this_pose = pose_weak_ptr(this);
        string lineStr;
        string record;
        string atomID, prevAtomID;
        string atomName;
        string altLoc;
        string resName;
        string chainID, prevChainID;
        string resID, prevResID;
        string iCode, prevICode;
        string xcoor, ycoor, zcoor;
        string occupancy;
        string tempFactor;
        string segID;
        string element;
        string charge;
        int lineLength, colonPosition;
        int internal_residue_index = -1;
        string astral, tempStr;
        atom_shared_ptr currAtom;
        residue_shared_ptr currResidue;
        chain_shared_ptr currChain;
        bool TER_rec = false;
        //char tempCStr[20];
        while(!input.eof()){
            getline(input, lineStr);
            lineLength = lineStr.length();
            record = lineStr.substr(0, 6);
            string record_trim = record;
            trim(record_trim);
            if(record_trim.compare("TER") == 0){
            	TER_rec = true;
            }
            if(record.compare("ENDMDL") == 0){
                //newpdb.reindex();
                this->index();
                return input;
            }
            if(record.compare("REMARK") == 0){
                const string this_remark = lineStr.substr(6, lineLength - 6);
                //cout << this_remark << endl;
                REMARK_vec.push_back(this_remark);
                if(lineLength > 18){
                    astral = lineStr.substr(11, 6);
                    if(astral.compare("ASTRAL") == 0){
                        colonPosition = lineStr.find(":");
                        tempStr = lineStr.substr(18, colonPosition - 18);
                        int type = 999;
                        string valueStr;
                        valueStr = lineStr.substr(colonPosition + 2, lineLength - colonPosition + 2);
                        for(int n = 0;n < AstralStrArrSize;n++){
                            if(AstralStrArr[n] != 0){
                                if(std::strcmp(tempStr.c_str(), AstralStrArr[n]) == 0){
                                    type = n;
                                }
                            }

                        }

                        switch (type){
                            case 0:
                                ASTRAL_version = valueStr;
                                //cout << newpdb.ASTRAL_version << endl;
                                break;
                            case 1:
                                SCOP_sid = valueStr;
                                //cout << newpdb.ASTRAL_version << endl;
                                break;
                            case 2:
                                SCOP_sun = valueStr;
                                //cout << newpdb.ASTRAL_version << endl;
                                break;
                            case 3:
                                SCOP_sccs = valueStr;
                                //cout << newpdb.ASTRAL_version << endl;
                                break;
                            case 4:
                                Source_PDB = valueStr;
                                pdbCode = Source_PDB;
                                //cout << newpdb.ASTRAL_version << endl;
                                break;
                            case 5:
                                Source_PDB_REVDAT = valueStr;
                                //cout << newpdb.ASTRAL_version << endl;
                                break;
                            case 6:
                                Region = valueStr;
                                //cout << newpdb.ASTRAL_version << endl;
                                break;
                            case 7:
                                ASTRAL_SPACI = valueStr;
                                //cout << newpdb.ASTRAL_version << endl;
                                break;
                            case 8:
                                ASTRAL_AEROSPACI = valueStr;
                                //cout << newpdb.ASTRAL_version << endl;
                                break;
                            case 9:
                                Data_updated_release = valueStr;
                                //cout << newpdb.ASTRAL_version << endl;
                                break;
                            default:
                                break;
                        }
                        //cout << "!" << tempStr << "! !" << valueStr << "!\n";
                    }

                }

            }
            else
                if(record.compare("CRYST1") == 0){
                    if(lineLength >= 70){
                        /*
				Contains:	 unit cell parameters, space group, and Z value
				Notes:
				If the structure was not determined by crystallographic means, simply defines a unit cube
				(a = b =c = 1.0, alpha = beta = gamma = 90
				degrees, space group = P 1, and Z = 1)
				The Hermann-Mauguin space group symbol is given without parenthesis,
				e.g., P 21 21 2 and using the full symbol, e.g., C 1 2 1 instead of C 2.
				The screw axis is described as a two digit number.
				For a rhombohedral space group in the hexagonal setting, the lattice type symbol used is H.
				The Z value is the number of polymeric chains in a unit cell. In the case of heteropolymers,
				Z is the number of occurrences of the most populous chain.
				In the case of a polycrystalline fiber diffraction study,
				CRYST1 and SCALE contain the normal unit cell data.
				The unit cell parameters are used to calculate SCALE.
				COLUMNS       DATA TYPE      CONTENTS
				--------------------------------------------------------------------------------
 	 	 	 	 1 -  6       Record name    "CRYST1"
 	 	 	 	 7 - 15       Real(9.3)      a (Angstroms)
				16 - 24       Real(9.3)      b (Angstroms)
				25 - 33       Real(9.3)      c (Angstroms)
				34 - 40       Real(7.2)      alpha (degrees)
				41 - 47       Real(7.2)      beta (degrees)
				48 - 54       Real(7.2)      gamma (degrees)
				56 - 66       LString        Space group
				67 - 70       Integer        Z value
				 */
                        string aStr = lineStr.substr(6, 9);
                        trim(aStr);
                        string bStr = lineStr.substr(15, 9);
                        trim(bStr);
                        string cStr = lineStr.substr(24, 9);
                        trim(cStr);
                        string alphaStr = lineStr.substr(33, 7);
                        trim(alphaStr);
                        string betaStr = lineStr.substr(40, 7);
                        trim(betaStr);
                        string gammaStr = lineStr.substr(47, 7);
                        trim(gammaStr);
                        string sgStr = lineStr.substr(55, 11);
                        trim(sgStr);
                        string zStr = lineStr.substr(66, 4);
                        trim(zStr);
                        try {
                            cryst_a = lexical_cast<double>(aStr);
                            cryst_b = lexical_cast<double>(bStr);
                            cryst_c = lexical_cast<double>(cStr);
                            cryst_alpha = lexical_cast<double>(alphaStr);
                            cryst_beta = lexical_cast<double>(betaStr);
                            cryst_gamma = lexical_cast<double>(gammaStr);
                            cryst_space_group = sgStr;
                            cryst_z_value = lexical_cast<int>(zStr);
                        }catch(boost::bad_lexical_cast & e){
                            std::cerr << "\nWARNING: PDB CRYST1 field has formatting problems" << endl;
                        }
                    }
                    else{
                        cerr << "\nWARNING: PDB contains bad CRYST1 field" << endl;
                    }
                }


            if((record.compare("ATOM  ") == 0) || record.compare("HETATM") == 0){
                resName = lineStr.substr(17, 3);
                trim(resName);
            }
            if((record.compare("ATOM  ") == 0) || (record.compare("HETATM") == 0 && (resName.compare("MSE") == 0 || resName.compare("CSS") == 0))){
                //cout << "Its an ATOM!"  << endl;
                atomID = lineStr.substr(6, 5);
                atomName = lineStr.substr(12, 4);
                trim(atomName);
                altLoc = lineStr.substr(16, 1);
                resName = lineStr.substr(17, 3);
                if(resName.compare("MSE") == 0){
                    resName.assign("MET");
                    if (atomName.compare("SE") == 0){
                    	atomName.assign("SD");
                    }
                }
                if(resName.compare("CSS") == 0){
                    resName.assign("CYS");
                }
                chainID = lineStr.substr(21, 1);
                resID = lineStr.substr(22, 5);
                //trim( resID );
                iCode = lineStr.substr(26, 1);
                xcoor = lineStr.substr(30, 8);
                trim(xcoor);
                ycoor = lineStr.substr(38, 8);
                trim(ycoor);
                zcoor = lineStr.substr(46, 8);
                trim(zcoor);
                double occ_dbl = 1.0;
                double b_factor_dbl = 0;
                if(lineLength >= 60){
                    occupancy = lineStr.substr(54, 6);
                    trim(occupancy);
                    try {
                        occ_dbl = lexical_cast<double>(occupancy);
                    }catch(boost::bad_lexical_cast & e){
                        std::cerr << "\nWARNING!!! bad occupancy field in PDB '" << occupancy << "'" << endl;
                        occ_dbl = 1.0;
                    }
                }

                if(lineLength >= 66){
                    tempFactor = lineStr.substr(60, 6);
                    trim(tempFactor);
                    try {
                        b_factor_dbl = lexical_cast<double>(tempFactor);
                    }catch(boost::bad_lexical_cast & e){
                        std::cerr << "\nWARNING!!! bad B-factor field in PDB '" << tempFactor << "'" << endl;
                        b_factor_dbl = 0.0;
                    }
                }

                if(lineLength >= 76){
                    segID = lineStr.substr(72, 4);
                }
                if(lineLength >= 78){
                    element = lineStr.substr(76, 2);
                }
                if(lineLength >= 80){
                    charge = lineStr.substr(78, 2);
                }
                if((chainID.compare(prevChainID) != 0) || TER_rec){
                    currChain = new_chain();
                    currChain->setChainID(chainID[0]);
                    currChain->set_pose(_this.lock());
                    currChain->setIsPeptide(false);
                    chain_vec.push_back(currChain);
                    TER_rec = false;
                }
                if(resID.compare(prevResID) != 0 || iCode.compare(prevICode) != 0){
                    internal_residue_index++;
                    //cout << "RESNAME::" << resName << "::" << endl;
                    //strcpy( tempCStr,resName.c_str() );
                    currResidue = new_residue();
                    currResidue->set_pdb_residue_index(resID);
                    currResidue->set_type(resName);
                    currResidue->set_icode(iCode[0]);
                    currResidue->set_pose(_this.lock());
                    currResidue->set_chain(currChain);
                    for(int i = 0;i < bb_map->get_num_protein_backbone_atoms();i++){
                        atom_shared_ptr natptr = new_atom();
                        natptr->set_pose(_this.lock());
                        natptr->set_residue(currResidue);
                        natptr->set_chain(currChain);
                        prot_backbone_atom_vec.push_back(natptr);
                    }
                    residue_vec.push_back(currResidue);
                }

                // checking for atomID compared to last here
                if(atomID.compare(prevAtomID) != 0 || check_prev_atomid == false){
                    //strcpy( tempCStr,atomName.c_str() );
                    const double x = lexical_cast<double>(xcoor);
                    const double y = lexical_cast<double>(ycoor);
                    const double z = lexical_cast<double>(zcoor);
                    const int rel_pos = bb_map->get_relative_location(atomName);
                    const int atom_index = (internal_residue_index * bb_map->get_num_protein_backbone_atoms()) + rel_pos;
                    //check if is bb atom
                    if(rel_pos < bb_map->get_num_protein_backbone_atoms()){
                        currAtom = prot_backbone_atom_vec[atom_index];
                        if(!currAtom->isSet()){
                            currAtom->set_type(atomName);
                            currAtom->set_coords(vector3d(x, y, z));
                            currAtom->set_b_factor(b_factor_dbl);
                            currAtom->set_occupancy(occ_dbl);
                            currAtom->setSet(true);
                            currAtom->setActive(true);
                            currAtom->set_pose(_this.lock());
                            currAtom->set_chain(currChain);
                            currAtom->set_residue(currResidue);
                            currChain->setIsPeptide(true);
                        }
                    }
                    else{
                        //is sidechain atom
                        currAtom = currResidue->get_sidechain()->get_atom(atomName);
                        if(!currAtom){
                            currAtom = currResidue->get_sidechain()->add_new_atom(vector3d(x, y, z), atomName);
                            currAtom->set_b_factor(b_factor_dbl);
                            currAtom->set_occupancy(occ_dbl);
                            currAtom->setSet(true);
                            currAtom->setActive(true);
                            currAtom->set_pose(_this.lock());
                            currAtom->set_chain(currChain);
                            currAtom->set_residue(currResidue);
                        }
                    }

                }
                else {

                }

                prevChainID = chainID;
                prevResID = resID;
                prevAtomID = atomID;
                prevICode = iCode;
            }

            //cout << lineStr << endl;
        }

        //cout << endl << "end of file reached." << endl;
        //cout << "Number of chains: " << newpdb.getChainCount() << endl;
        //cout << "Sequence of current chain:\n" << currChain->get3LetterSeq(" ") << endl;
        this->index();
        return input;
    }

    void pose::outputPdb_ATOM_line(std::ostream & output, const atom_shared_ptr currAtom, const residue_shared_ptr currResidue, const chain_shared_ptr currChain, const int atom_seq_num) const
    {
    	string atom_seq_num_str = boost::lexical_cast<string>(atom_seq_num);
    	if (atom_seq_num_str.length() > 5){
    		atom_seq_num_str = atom_seq_num_str.substr(atom_seq_num_str.length() - 5);
    	}
        output << setw(6) << "ATOM  " << setw(5) << atom_seq_num_str << " " << setw(4) << setiosflags(ios::left) << currAtom->get_type().get_pdb_formatted_label() << " " << resetiosflags(ios::left) << setw(3) << currResidue->get_type().get_pdb_formatted_label() << " " << setw(1) << currChain->getChainID() << setw(4) << currResidue->get_pdb_residue_index() << currResidue->get_icode() << "   " << setiosflags(ios::fixed) << setw(8) << setprecision(3) << currAtom->get_coords().x << setw(8) << setprecision(3) << currAtom->get_coords().y << setw(8) << setprecision(3) << currAtom->get_coords().z << setw(6) << setprecision(2) << currAtom->get_occupancy() << setw(6) << setprecision(2) << currAtom->get_b_factor() << "\n";
    }

    void pose::outputPdb(std::ostream & output, bool outputCAonly) const
    {
        const bool CAonly = outputCAonly;
        const int a_model_num = model_number;
        this->outputREMARKs(output);
        this->outputCRYST1(output);
        output << "MODEL " << "    " << resetiosflags(ios::left) << setw(4) << a_model_num << "\n";
        atom_shared_ptr currAtom;
        residue_shared_ptr currResidue;
        chain_shared_ptr currChain, prevChain;
        sidechain_shared_ptr currSC;
        int atom_seq_num = 1;
        const int resCount = this->get_residue_count();
        for(int resNum = 0;resNum < resCount;resNum++){
            currResidue = residue_vec[resNum];
            currChain = currResidue->get_chain();
            if(prevChain != currChain && atom_seq_num != 1)
                output << "TER   \n";

            for(int rel_bb_at_num = 0;rel_bb_at_num < this->bb_map->get_num_protein_backbone_atoms();rel_bb_at_num++){
                currAtom = prot_backbone_atom_vec[this->bb_map->get_output_order(rel_bb_at_num) + currResidue->first_bb_atom_index];
                // MIS quick hack
                atom_type T = currAtom->get_type();
                atom_type CA("CA");
                if(!(T == CA) && CAonly)
                    continue;

                if(currAtom->isSet() && currAtom->isActive()){
                    //if (prevChain != currChain && atom_seq_num != 1) output << "TER   \n";
                    outputPdb_ATOM_line(output, currAtom, currResidue, currChain, atom_seq_num);
                    atom_seq_num++;
                }
            }

            if(!CAonly){
                /*
			prevChain=currChain;
			continue;
		}
			 */
                currSC = currResidue->get_sidechain();
                const int sc_atom_count = currSC->get_atom_count();
                for(int sc_atom_num = 0;sc_atom_num < sc_atom_count;sc_atom_num++){
                    currAtom = currSC->get_atom(sc_atom_num);
                    if(currAtom->isSet() && currAtom->isActive()){
                        //if (prevChain != currChain && atom_seq_num != 1) output << "TER   \n";
                        outputPdb_ATOM_line(output, currAtom, currResidue, currChain, atom_seq_num);
                        atom_seq_num++;
                    }
                }

            }

            prevChain = currChain;
        }

        output << "TER   \n";
        output << "ENDMDL" << endl;
        this->output_appendlines(output);
    }

    int pose::get_all_atom_count() const
    {
        return all_atoms_vec.size();
    }

    pose::pose(const pose & old_pose)
    {
        this->init();
        instance_num = instance_count;
        this->set_version_date_strings();
        instance_count++;
        pdbCode = old_pose.pdbCode;
        REMARK_vec = old_pose.REMARK_vec;
        suppress_REMARK_output_flag = old_pose.suppress_REMARK_output_flag;
        appendlines_vec = old_pose.appendlines_vec;
        suppress_appendlines_output_flag = old_pose.suppress_appendlines_output_flag;
        ASTRAL_version = old_pose.ASTRAL_version;
        SCOP_sid = old_pose.SCOP_sid;
        SCOP_sun = old_pose.SCOP_sun;
        SCOP_sccs = old_pose.SCOP_sccs;
        Source_PDB = old_pose.Source_PDB;
        Source_PDB_REVDAT = old_pose.Source_PDB_REVDAT;
        Region = old_pose.Region;
        ASTRAL_SPACI = old_pose.ASTRAL_SPACI;
        ASTRAL_AEROSPACI = old_pose.ASTRAL_AEROSPACI;
        Data_updated_release = old_pose.Data_updated_release;
    }

    pose_shared_ptr pose::clone() const
    {
        atom_shared_ptr currAtom, prevAtom;
        residue_shared_ptr currResidue, prevResidue;
        chain_shared_ptr currChain, prevChain;
        sidechain_shared_ptr currSC;
        this->index();
        boost::shared_ptr<pose> npose(new pose(*this));
        npose->_this = npose;
        //const int resCount = this->get_residue_count();
        const int chainCount = this->chain_vec.size();
        for(int chainNum = 0;chainNum < chainCount;chainNum++){
            currChain = chain_vec[chainNum]->clone();
            npose->chain_vec.push_back(currChain);
            currChain->set_pose(npose);
            const int start_res = chain_vec[chainNum]->get_first_internal_residue_index();
            const int end_res = chain_vec[chainNum]->get_last_internal_residue_index();
            for(int resnum = start_res;resnum <= end_res;resnum++){
                currResidue = residue_vec[resnum]->clone();
                npose->residue_vec.push_back(currResidue);
                currResidue->set_pose(npose);
                currResidue->set_chain(currChain);
                currSC = residue_vec[resnum]->sc->clone();
                currResidue->sc = currSC;
                currSC->set_residue(currResidue);
                const int sc_atom_count = currSC->get_atom_count();
                for(int sc_atom_num = 0;sc_atom_num < sc_atom_count;sc_atom_num++){
                    currAtom = currSC->get_atom(sc_atom_num);
                    currAtom->set_pose(npose);
                    currAtom->set_residue(currResidue);
                    currAtom->set_chain(currChain);
                    currAtom->set_sidechain(currSC);
                }
                const int start_bb_at = residue_vec[resnum]->get_first_bb_atom_index();
                for(int bb_atom_num = start_bb_at;bb_atom_num < start_bb_at + this->bb_map->get_num_protein_backbone_atoms();bb_atom_num++){
                    currAtom = prot_backbone_atom_vec[bb_atom_num]->clone();
                    npose->prot_backbone_atom_vec.push_back(currAtom);
                    currAtom->set_pose(npose);
                    currAtom->set_residue(currResidue);
                    currAtom->set_chain(currChain);
                }
            }

        }

        npose->index();
        return npose;
    }

    void pose::index() const
    {
        atom_shared_ptr currAtom, prevAtom;
        residue_shared_ptr currResidue, prevResidue;
        chain_shared_ptr currChain, prevChain;
        sidechain_shared_ptr currSC;
        const int resCount = this->get_residue_count();
        all_atoms_vec.resize(0);
        all_atoms_vec.reserve(resCount * 50);
        const int chainCount = this->chain_vec.size();
        for(int chainNum = 0;chainNum < chainCount;chainNum++){
            currChain = chain_vec[chainNum];
            currChain->set_first_internal_residue_index(-1);
            currChain->set_first_bb_atom_index(-1);
            currChain->set_last_internal_residue_index(-1);
            currChain->set_last_bb_atom_index(-1);
        }
        int resNum = 0;
        //const int resCount = residue_vec.size();
        for(resNum = 0;resNum < resCount;resNum++){
            currResidue = residue_vec[resNum];
            currChain = currResidue->get_chain();
            currResidue->set_prev_residue(residue_shared_ptr());
            currResidue->set_next_residue(residue_shared_ptr());
            if(currChain && currChain != prevChain){
                currChain->set_first_internal_residue_index(resNum);
                currChain->set_first_bb_atom_index(resNum * this->bb_map->get_num_protein_backbone_atoms());
                if(prevChain){
                    prevChain->set_last_internal_residue_index(resNum - 1);
                    prevChain->set_last_bb_atom_index((resNum * this->bb_map->get_num_protein_backbone_atoms()) - 1);
                }
            }

            if(currChain == prevChain){
                currResidue->set_prev_residue(prevResidue);
                if(prevResidue){
                    prevResidue->set_next_residue(currResidue);
                }
            }

            currResidue->set_internal_residue_index(resNum);
            currResidue->set_first_bb_atom_index(resNum * this->bb_map->get_num_protein_backbone_atoms());
            prevChain = currChain;
            prevResidue = currResidue;
        }

        if(currChain){
            currChain->set_last_internal_residue_index(resNum - 1);
            currChain->set_last_bb_atom_index((resNum * this->bb_map->get_num_protein_backbone_atoms()) - 1);
        }
        int atom_seq_num = 0;
        for(int resNum = 0;resNum < resCount;resNum++){
            currResidue = residue_vec[resNum];
            currChain = currResidue->get_chain();
            for(int rel_bb_at_num = 0;rel_bb_at_num < this->bb_map->get_num_protein_backbone_atoms();rel_bb_at_num++){
                currAtom = prot_backbone_atom_vec[this->bb_map->get_output_order(rel_bb_at_num) + currResidue->first_bb_atom_index];
                if(currAtom->isSet()){
                    currAtom->set_seq_num(atom_seq_num);
                    all_atoms_vec.push_back(currAtom);
                    atom_seq_num++;
                }
            }

            currSC = currResidue->get_sidechain();
            const int sc_atom_count = currSC->get_atom_count();
            for(int sc_atom_num = 0;sc_atom_num < sc_atom_count;sc_atom_num++){
                currAtom = currSC->get_atom(sc_atom_num);
                if(currAtom->isSet()){
                    currAtom->set_seq_num(atom_seq_num);
                    all_atoms_vec.push_back(currAtom);
                    atom_seq_num++;
                }
            }

        }

    }

    /* NOT SAFE
    void pose::backup_coords()
    {
        atom_shared_ptr_vector::iterator iter;
        for(iter = all_atoms_vec.begin();iter != all_atoms_vec.end();iter++){
            (*iter)->backup_coords();
        }
    }
    */

    atom_shared_ptr_vector3d_map pose::get_all_coords(){
    	atom_shared_ptr_vector3d_map rtn_vec;
    	for (unsigned int i = 0; i < (unsigned int)this->get_all_atom_count(); i++ ){
    		rtn_vec[this->get_atom(i)] = this->get_atom(i)->get_coords();
    	}
    	return rtn_vec;
    }


    /* NOT SAFE
    void pose::restore_backed_up_coords()
    {
        atom_shared_ptr_vector::iterator iter;
        for(iter = all_atoms_vec.begin();iter != all_atoms_vec.end();iter++){
            (*iter)->restore_backed_up_coords();
        }
    }
    */

	bool pose::restore_coords(atom_shared_ptr_vector3d_map& backup ){
		atom_shared_ptr_vector3d_map::const_iterator it;

		bool isOK = true;

		pose_shared_ptr this_locked = _this.lock();

		for (it = backup.begin(); it != backup.end(); it++){
			if (it->first->get_pose() == this_locked){
				it->first->set_coords(it->second);
			}
			else {
				cerr << "WARNING: pose::restore_coords: backing up coords to wrong pose" << endl;
				isOK = false;
			}
		}

		return isOK;
	}

    std::string pose::get_label() const
    {
        return this->label;
    }

    void pose::set_label(std::string new_label)
    {
        this->label = new_label;
    }

    //! testing don't use yet
    residue_shared_ptr pose::add_nterm_residue(const residue_type type, const int chain_num)
    {
        chain_shared_ptr currChain = this->get_chain(chain_num);
        if(currChain){
            residue_shared_ptr currResidue = new_residue();
            residue_shared_ptr old_n_term = residue_shared_ptr();
            const int old_n_term_index = currChain->get_first_internal_residue_index();
            if(old_n_term_index != -1){
                old_n_term = residue_vec[old_n_term_index];
                string old_resID = old_n_term->get_pdb_residue_index();
                trim(old_resID);
                string new_resID;
                try {
                    int old_num = lexical_cast<int>(old_resID);
                    new_resID = lexical_cast<string>(old_num - 1);
                }catch(boost::bad_lexical_cast & e){
                    std::cerr << "\nWARNING!!! can't cast old N-terminal residue ID as an integer :'" << old_resID << "' " << "adding residue ID as :'" << 0 << endl;
                    new_resID = lexical_cast<string>(0);
                }
                currResidue->set_pdb_residue_index(new_resID);
            }

            currResidue->set_type(type);
            currResidue->set_icode(' ');
            currResidue->set_pose(_this.lock());
            currResidue->set_chain(currChain);
            const int first_bb_num = currChain->get_first_bb_atom_index();
            for(int i = 0;i < bb_map->get_num_protein_backbone_atoms();i++){
                atom_shared_ptr natptr = new_atom();
                natptr->set_pose(_this.lock());
                natptr->set_residue(currResidue);
                natptr->set_chain(currChain);
                prot_backbone_atom_vec.insert(prot_backbone_atom_vec.begin() + first_bb_num, natptr);
                //prot_backbone_atom_vec.push_back(natptr);
            }
            residue_vec.insert(residue_vec.begin() + old_n_term_index, currResidue);
            //residue_vec.push_back(currResidue);
            this->index();
            return currResidue;
        }
        else{
            std::cerr << "ERROR: pose: can not add residue - chain_num does not exist" << endl;
            return residue_shared_ptr();
        }
        //std::cerr << "ERROR: pose: can not add residue - not implemented yet" << endl;
        //return residue_shared_ptr();
    }

    //! testing don't use yet
    residue_shared_ptr pose::add_cterm_residue(const residue_type type, const int chain_num)
    {
        chain_shared_ptr currChain = this->get_chain(chain_num);
        if(currChain){
            return this->append_residue(type, currChain->getChainID());
        }else{
            std::cerr << "ERROR: pose: can not add residue - chain_num does not exist" << endl;
            return residue_shared_ptr();
        }
        /*
	if (currChain){
		residue_shared_ptr currResidue = new_residue();
		residue_shared_ptr old_c_term = residue_shared_ptr();
		const int old_c_term_index = currChain->get_last_internal_residue_index();
		if (old_c_term_index != -1){
			old_c_term = residue_vec[currChain->get_last_internal_residue_index()];
			string old_resID = old_c_term->get_pdb_residue_index();
			trim(old_resID);
			string new_resID;
			try {
				int old_num = lexical_cast<int>(old_resID);
				new_resID = lexical_cast<string>(old_num+1);
			}
			catch(boost::bad_lexical_cast& e) {
				std::cerr << "\nWARNING!!! can't cast old C-terminal residue ID as an integer :'" << old_resID << "' "
						<< "adding residue ID as :'" << this->get_residue_count() + 1
						<< endl;
				new_resID = lexical_cast<string>(this->get_residue_count() + 1);
			}
			currResidue->set_pdb_residue_index( new_resID );
		}
		currResidue->set_type(type);
		currResidue->set_icode(' ');

		currResidue->set_pose(_this.lock());
		currResidue->set_chain(currChain);

		const int last_bb_num = currChain->get_last_bb_atom_index();

		for (int i = 0; i < bb_map->get_num_protein_backbone_atoms(); i++){
			atom_shared_ptr natptr= new_atom();
			natptr->set_pose(_this.lock());
			natptr->set_residue(currResidue);
			natptr->set_chain(currChain);
			//prot_backbone_atom_vec.push_back(natptr);
			prot_backbone_atom_vec.insert(prot_backbone_atom_vec.begin()+last_bb_num+1, natptr);
		}

		//residue_vec.push_back(currResidue);
		residue_vec.insert(residue_vec.begin()+old_c_term_index+1, currResidue);

		this->index();

		return currResidue;
	}
	else {
		std::cerr << "ERROR: pose: can not add residue - chain_num does not exist" << endl;
		return residue_shared_ptr();
	}
	*/
    }

    residue_shared_ptr pose::insert_residue(const residue_type type, const int internal_residue_num)
    {
        if(internal_residue_num >= this->get_residue_count() || internal_residue_num < 0){
            std::cerr << "ERROR: pose: can not insert residue - internal_residue_num out of bounds:\t" << internal_residue_num << endl;
            return residue_shared_ptr();
        }
        residue_shared_ptr ins_res = this->get_residue(internal_residue_num);
        chain_shared_ptr currChain = ins_res->get_chain();
        if(currChain){
            residue_shared_ptr currResidue = new_residue();
            residue_shared_ptr old_n_term = residue_shared_ptr();
            const int old_n_term_index = internal_residue_num; //currChain->get_first_internal_residue_index();
            if(old_n_term_index != -1){
                old_n_term = ins_res; //residue_vec[old_n_term_index];
                string old_resID = old_n_term->get_pdb_residue_index();
                trim(old_resID);
                string new_resID;
                try {
                    int old_num = lexical_cast<int>(old_resID);
                    new_resID = lexical_cast<string>(old_num - 1);
                }catch(boost::bad_lexical_cast & e){
                    std::cerr << "\nWARNING!!! can't cast insert point residue ID as an integer :'" << old_resID << "' " << "adding residue ID as :'" << 0 << endl;
                    new_resID = lexical_cast<string>(0);
                }
                currResidue->set_pdb_residue_index(new_resID);
            }

            currResidue->set_type(type);
            currResidue->set_icode(' ');
            currResidue->set_pose(_this.lock());
            currResidue->set_chain(currChain);
            const int first_bb_num = old_n_term->get_first_bb_atom_index(); //currChain->get_first_bb_atom_index();
            for(int i = 0;i < bb_map->get_num_protein_backbone_atoms();i++){
                atom_shared_ptr natptr = new_atom();
                natptr->set_pose(_this.lock());
                natptr->set_residue(currResidue);
                natptr->set_chain(currChain);
                prot_backbone_atom_vec.insert(prot_backbone_atom_vec.begin() + first_bb_num, natptr);
                //prot_backbone_atom_vec.push_back(natptr);
            }
            residue_vec.insert(residue_vec.begin() + old_n_term_index, currResidue);
            //residue_vec.push_back(currResidue);
            this->index();
            return currResidue;
        }
        else{
            std::cerr << "ERROR: pose: can not add residue - chain_num does not exist" << endl;
            return residue_shared_ptr();
        }
    }

    residue_shared_ptr pose::append_residue(const residue_type type, const char chainID)
    {
        chain_shared_ptr currChain = this->get_chain(chainID);
        if(!currChain){
            cerr << "pose: ERROR: chainID: " << chainID << " does not exist" << endl;
            return residue_shared_ptr();
        }
        const int old_c_term_res_num = currChain->get_last_internal_residue_index();
        residue_shared_ptr old_c_term = residue_shared_ptr(); //
        if(old_c_term_res_num >= 0 && old_c_term_res_num < this->get_residue_count()){
            old_c_term = this->get_residue(old_c_term_res_num);
        }
        residue_shared_ptr currResidue = new_residue();
        int new_res_num = 1;
        if(old_c_term){
            try {
                string str = old_c_term->get_pdb_residue_index();
                trim(str);
                int old_num = lexical_cast<int>(str);
                new_res_num = (old_num + 1);
            }catch(boost::bad_lexical_cast & e){
                new_res_num = 1;
                std::cerr << "\nWARNING!!! can't cast old C-term residue ID as an integer :'" << old_c_term->get_pdb_residue_index() << "' " << "adding residue ID as :'" << new_res_num << endl;
                //new_resID = lexical_cast<string>(new_res_num);
            }
        }

        string new_resID = lexical_cast<string>(new_res_num);
        currResidue->set_pdb_residue_index(new_resID);
        currResidue->set_type(type);
        currResidue->set_icode(' ');
        currResidue->set_pose(_this.lock());
        currResidue->set_chain(currChain);
        const int first_bb_num = old_c_term ? old_c_term->get_first_bb_atom_index() + bb_map->get_num_protein_backbone_atoms() : prot_backbone_atom_vec.size();
        for(int i = 0;i < bb_map->get_num_protein_backbone_atoms();i++){
            atom_shared_ptr natptr = new_atom();
            natptr->set_pose(_this.lock());
            natptr->set_residue(currResidue);
            natptr->set_chain(currChain);
            prot_backbone_atom_vec.insert(prot_backbone_atom_vec.begin() + first_bb_num, natptr);
            //prot_backbone_atom_vec.push_back(natptr);
        }
        if(old_c_term){
            residue_vec.insert(residue_vec.begin() + old_c_term_res_num + 1, currResidue);
        }else{
            residue_vec.push_back(currResidue);
        }
        this->index();
        return currResidue;
    }

    chain_shared_ptr pose::add_new_chain(const char chainID, bool peptide_chain)
    {
        if(this->get_chain(chainID)){
            cerr << "pose: ERROR: chainID: " << chainID << " is already assigned" << endl;
            return chain_shared_ptr();
        }
        chain_shared_ptr currChain = new_chain();
        currChain->setChainID(chainID);
        currChain->set_pose(_this.lock());
        currChain->setIsPeptide(peptide_chain);
        chain_vec.push_back(currChain);
        this->index();
        return currChain;
    }

    //!
    chain_shared_ptr pose::add_new_chain(bool peptide_chain)
    {
    	std::string alphabet = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789";
        char char_to_use = 'Z';
        for(unsigned int ii = 0; ii < alphabet.size(); ii++){
            if(!this->get_chain(alphabet[ii])){
                char_to_use = alphabet[ii];
                break;
            }
        }

        return add_new_chain(char_to_use, peptide_chain);
    }

    // TODO: debug this as it seems to be adding duplicate H atoms
    chain_shared_ptr pose::add_duplicated_chain(const_chain_shared_ptr chain_to_copy)
    {
        const_pose_shared_ptr p_to_cp = chain_to_copy->get_pose();
        chain_shared_ptr added_ch = get_chain(chain_to_copy->getChainID()) ? add_new_chain(chain_to_copy->isPeptide()) : add_new_chain(chain_to_copy->getChainID(), chain_to_copy->isPeptide());
        for(int i = chain_to_copy->get_first_internal_residue_index();i <= chain_to_copy->get_last_internal_residue_index();i++){
            const_residue_shared_ptr r_to_cp = p_to_cp->get_residue(i);
            residue_shared_ptr new_r = r_to_cp->clone();
            new_r->set_pose(_this.lock());
            new_r->set_chain(added_ch);
            for(int rn = 0;rn < bb_map->get_num_protein_backbone_atoms();rn++){
                atom_shared_ptr natptr = p_to_cp->get_bb_atom(bb_map->get_BBAtomType_from_rel_loc(rn), i)->clone();
                natptr->set_pose(_this.lock());
                natptr->set_residue(new_r);
                natptr->set_chain(added_ch);
                prot_backbone_atom_vec.push_back(natptr);
            }
            new_r->get_sidechain()->copy_sidechain(r_to_cp->get_sidechain());
            residue_vec.push_back(new_r);
        }

        this->index();
        return added_ch;
    }

    void pose::delete_residue(const unsigned int internal_residue_num)
    {
        if(internal_residue_num < residue_vec.size() ){
            residue_shared_ptr_vector::iterator del_res_it = residue_vec.begin() + internal_residue_num;
            const int first_bb_atom = (*del_res_it)->get_first_bb_atom_index();
            const int last_bb_atom = first_bb_atom + bb_map->get_num_protein_backbone_atoms() - 1;
            residue_vec.erase(residue_vec.begin() + internal_residue_num);
            prot_backbone_atom_vec.erase(prot_backbone_atom_vec.begin() + first_bb_atom, prot_backbone_atom_vec.begin() + last_bb_atom + 1);
            this->index();
        }else{
            std::cerr << "ERROR: pose: can not delete non-existent residue " << internal_residue_num << endl;
        }
    }

    bool pose::coords_numerically_ok() const
    {
        const int atom_count = this->get_all_atom_count();
        bool ok = true;
        //cout << "pose: coords_numerically_ok: checking" << endl;
        for(int i = 0;i < atom_count;i++){
            const_atom_shared_ptr atm = this->get_atom(i);
            vector3d v3d = atm->get_coords();
            for(int r = 0;r < 3;r++){
                if(!boost::math::isfinite(v3d[r])){
                    //(isnan(v3d[r]) || isinf(v3d[r])){
                    const_residue_shared_ptr res = atm->get_residue();
                    ok = false;
                    cerr << "pose: problem at atom " << i << " from residue " << res->get_pdb_residue_index() << " its value is " << v3d << endl;
                }
            }

        }

        return ok;
    }

    int pose::get_model_number() const
    {
        return model_number;
    }

    void pose::set_model_number(int modelNumber)
    {
        model_number = modelNumber;
}



}
}

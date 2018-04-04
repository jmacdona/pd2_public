//
// (c)  JAMES T. MACDONALD 2010 
// This file is part of the PRODART 2 software
// suite and is made available under license.
//
// For more information please contact by email: j.t.macdonald+prodart@gmail.com
//
/*
 * backbone_map.cpp
 *
 *  Created on: Jan 31, 2010
 *      Author: jmacdon
 */

#include "backbone_map.h"

using namespace std;

namespace PRODART {
namespace POSE {

namespace {
const prot_backbone_map *instance = NULL;
boost::once_flag once_flag = BOOST_ONCE_INIT;
}

void prot_backbone_map::Init(){
	if (!instance){
		instance = new const prot_backbone_map;
	}
}

const prot_backbone_map* prot_backbone_map::Instance(){
	/*
	if (!instance){
		instance = new const prot_backbone_map;
	}
	*/

	boost::call_once(&prot_backbone_map::Init, once_flag);

	return instance;
}

template <class T>
T **MatInit(const int rows, const int cols)
{
    int             i;
    T        **matrix = NULL;
    T         *matspace = NULL;

    matspace = (T *) calloc((rows * cols), sizeof(T));
    if (matspace == NULL)
    {
        perror("\n ERROR");
        printf("\n ERROR: Failure to allocate matrix space in MatInit(): (%d x %d)\n", rows, cols);
        exit(EXIT_FAILURE);
    }

    /* allocate room for the pointers to the rows */
    matrix = (T **) malloc(rows * sizeof(T *));
    if (matrix == NULL)
    {
        perror("\n ERROR");
        printf("\n ERROR: Failure to allocate room for row pointers in MatInit(): (%d)\n", rows);
        exit(EXIT_FAILURE);
    }

    /*  now 'point' the pointers */
    for (i = 0; i < rows; i++)
        matrix[i] = matspace + (i * cols);

    return(matrix);
}

template <class T>
void MatDestroy(T ***matrix_ptr)
{
    T **matrix = *matrix_ptr;

    if (matrix != NULL)
    {
        if (matrix[0] != NULL)
        {
            free(matrix[0]);
            matrix[0] = NULL;
        }

        free(matrix);
        *matrix_ptr = NULL;
    }
}

prot_backbone_map::prot_backbone_map(){

	//cerr << "prot_backbone_map: initialising..." << endl;

	assert(ENV::is_command_line_parsed());


	at_map.clear();
	this->at_enum_vec.clear();
	at_enum_vec.resize(20, num_protein_backbone_atoms);
	inv_at_enum_vec.resize(num_protein_backbone_atoms+1,num_protein_backbone_atoms);
	/*enum AtomType {UNDEF_ATOM,
	 * N,	0
	 * H, 	1
	 * CA, 	2
	 * 1HA, 3
	 * 2HA, 4
	 * CB, 	5
	 * 1HB, 6
	 * 2HB, 7
	 * 3HB, 8
	 * C, 	9
	 * O, 	10
	 * H2	11
	 * H3	12
	 * OXT	13
	 * END_ATOMTYPE_ENUM};

	 */

	at_map[atom_type("N")] = 0;
	at_map[atom_type("H")] = 1;
	at_map[atom_type("HN")] = 1;
	at_map[atom_type("CA")] = 2;
	at_map[atom_type("HA")] = 3;
	at_map[atom_type("1HA")] = 3;
	at_map[atom_type("2HA")] = 4;
	at_map[atom_type("HA1")] = 3;
	at_map[atom_type("HA2")] = 4;
	at_map[atom_type("CB")] = 5;
	at_map[atom_type("1HB")] = 6;
	at_map[atom_type("2HB")] = 7;
	at_map[atom_type("3HB")] = 8;
	at_map[atom_type("HB1")] = 6;
	at_map[atom_type("HB2")] = 7;
	at_map[atom_type("HB3")] = 8;
	at_map[atom_type("C")] = 9;
	at_map[atom_type("O")] = 10;
	at_map[atom_type("1H")] = 1;
	at_map[atom_type("2H")] = 11;
	at_map[atom_type("3H")] = 12;
	at_map[atom_type("OXT")] = 13;

	at_map[atom_type(" N")] = 0;
	at_map[atom_type(" H")] = 1;
	at_map[atom_type(" HN")] = 1;
	at_map[atom_type(" CA")] = 2;
	at_map[atom_type(" HA")] = 3;
	at_map[atom_type(" 1HA")] = 3;
	at_map[atom_type(" 2HA")] = 4;
	at_map[atom_type(" HA1")] = 3;
	at_map[atom_type(" HA2")] = 4;
	at_map[atom_type(" CB")] = 5;
	at_map[atom_type(" 1HB")] = 6;
	at_map[atom_type(" 2HB")] = 7;
	at_map[atom_type(" 3HB")] = 8;
	at_map[atom_type(" HB1")] = 6;
	at_map[atom_type(" HB2")] = 7;
	at_map[atom_type(" HB3")] = 8;
	at_map[atom_type(" C")] = 9;
	at_map[atom_type(" O")] = 10;
	at_map[atom_type(" 1H")] = 1;
	at_map[atom_type(" 2H")] = 11;
	at_map[atom_type(" 3H")] = 12;
	at_map[atom_type(" OXT")] = 13;



	//enum BBAtomType { N, H, CA, HA1, HA2, CB, HB1, HB2, HB3, C,  O, HA};


	at_enum_vec[static_cast<int>(N)] = 0;
	at_enum_vec[static_cast<int>(H)] = 1;
	at_enum_vec[static_cast<int>(CA)] = 2;
	at_enum_vec[static_cast<int>(HA1)] = 3;
	at_enum_vec[static_cast<int>(HA2)] = 4;
	at_enum_vec[static_cast<int>(CB)] = 5;
	at_enum_vec[static_cast<int>(HB1)] = 6;
	at_enum_vec[static_cast<int>(HB2)] = 7;
	at_enum_vec[static_cast<int>(HB3)] = 8;
	at_enum_vec[static_cast<int>(C)] = 9;
	at_enum_vec[static_cast<int>(O)] = 10;
	at_enum_vec[static_cast<int>(HA)] = 3;
	at_enum_vec[static_cast<int>(H1)] = 1;
	at_enum_vec[static_cast<int>(H2)] = 11;
	at_enum_vec[static_cast<int>(H3)] = 12;
	at_enum_vec[static_cast<int>(OXT)] = 13;
	at_enum_vec[static_cast<int>(END_BBAtomType)] = num_protein_backbone_atoms;

	for (int i = 0; i < num_protein_backbone_atoms; i++){
		int type = num_protein_backbone_atoms;
		for (unsigned int v = 0; v <  at_enum_vec.size(); v++){
			if (at_enum_vec[v] == i){
				type = i;
				break;
			}
		}
		inv_at_enum_vec[i] = type;
	}



	pdb_output_order_vec.resize(num_protein_backbone_atoms, 0);

	pdb_output_order_vec[0] = this->get_relative_location(N);
	pdb_output_order_vec[1] = this->get_relative_location(CA);
	pdb_output_order_vec[2] = this->get_relative_location(C);
	pdb_output_order_vec[3] = this->get_relative_location(O);
	pdb_output_order_vec[4] = this->get_relative_location(OXT);
	pdb_output_order_vec[5] = this->get_relative_location(CB);
	pdb_output_order_vec[6] = this->get_relative_location(H);
	pdb_output_order_vec[7] = this->get_relative_location(H2);
	pdb_output_order_vec[8] = this->get_relative_location(H3);
	pdb_output_order_vec[9] = this->get_relative_location(HA);
	pdb_output_order_vec[10] = this->get_relative_location(HA2);
	pdb_output_order_vec[11] = this->get_relative_location(HB1);
	pdb_output_order_vec[12] = this->get_relative_location(HB2);
	pdb_output_order_vec[13] = this->get_relative_location(HB3);


	small_bb_atom_list.clear();
	small_bb_atom_list.reserve(8);
	small_bb_atom_list.push_back(N);
	small_bb_atom_list.push_back(CA);
	small_bb_atom_list.push_back(C);
	small_bb_atom_list.push_back(O);
	small_bb_atom_list.push_back(CB);
	small_bb_atom_list.push_back(H);

	smaller_bb_atom_list.clear();
	smaller_bb_atom_list.reserve(8);
	smaller_bb_atom_list.push_back(N);
	smaller_bb_atom_list.push_back(CA);
	smaller_bb_atom_list.push_back(C);
	smaller_bb_atom_list.push_back(O);
	smaller_bb_atom_list.push_back(CB);
	//smaller_bb_atom_list.push_back(H);

	bb_hydrogen_atom_list.clear();
	bb_hydrogen_atom_list.push_back(H);
	bb_hydrogen_atom_list.push_back(HA1);
	bb_hydrogen_atom_list.push_back(HA2);
	bb_hydrogen_atom_list.push_back(HB1);
	bb_hydrogen_atom_list.push_back(HB2);
	bb_hydrogen_atom_list.push_back(HB3);
	bb_hydrogen_atom_list.push_back(HA);
	bb_hydrogen_atom_list.push_back(H1);
	bb_hydrogen_atom_list.push_back(H2);
	bb_hydrogen_atom_list.push_back(H3);





	full_bb_atom_list.clear();
	for (int i = 0 ; i < this->get_num_protein_backbone_atoms(); i++){
		const BBAtomType atty  = this->get_BBAtomType_from_rel_loc(i);
		full_bb_atom_list.push_back(atty);
	}

	this->init_bond_sep_data();

	//cerr << "prot_backbone_map: init ok" << endl;

}

int prot_backbone_map::get_output_order(const int index) const{
	if (index < num_protein_backbone_atoms){
		return pdb_output_order_vec[index];
	}
	else {
		return num_protein_backbone_atoms;
	}
}

prot_backbone_map::~prot_backbone_map(){

	MatDestroy(&internal_residue_bond_sep_mat);

	delete instance;
	instance = 0;
}


// TODO add rest of BB atoms
void prot_backbone_map::init_bond_sep_data(){
	internal_residue_bond_sep_mat = MatInit<int>(num_protein_backbone_atoms,num_protein_backbone_atoms);
	bonds_to_N.resize(num_protein_backbone_atoms,0);
	bonds_to_C.resize(num_protein_backbone_atoms,0);

	//const int init_val = static_cast<int>(N);
	for (BBAtomType at1 = N; at1 < END_BBAtomType; at1=static_cast<BBAtomType>(at1+1)){
		const int int_at1 = get_relative_location(at1);//static_cast<int>(at1) - init_val;
		for (BBAtomType at2 = N; at2 < END_BBAtomType; at2=static_cast<BBAtomType>(at2+1)){
			const int int_at2 = get_relative_location(at2);//static_cast<int>(at2) - init_val;


			//std::cout << int_at1 << "\t" << int_at2 << std::endl;

			if (get_relative_location(at1) == get_relative_location(at2)){
				internal_residue_bond_sep_mat[int_at1][int_at2] = 0;
				internal_residue_bond_sep_mat[int_at2][int_at1] = 0;

			}
			else if (get_relative_location(at1)==get_relative_location(H)){
				if (get_relative_location(at2) == get_relative_location(N)){
					internal_residue_bond_sep_mat[int_at1][int_at2] = 1;
					internal_residue_bond_sep_mat[int_at2][int_at1] = 1;
				}
				else if (get_relative_location(at2) == get_relative_location(CA)){
					internal_residue_bond_sep_mat[int_at1][int_at2] = 2;
					internal_residue_bond_sep_mat[int_at2][int_at1] = 2;
				}
				else if (get_relative_location(at2) == get_relative_location(CB)){
					internal_residue_bond_sep_mat[int_at1][int_at2] = 3;
					internal_residue_bond_sep_mat[int_at2][int_at1] = 3;
				}
				else if (get_relative_location(at2) == get_relative_location(C)){
					internal_residue_bond_sep_mat[int_at1][int_at2] = 3;
					internal_residue_bond_sep_mat[int_at2][int_at1] = 3;
				}
				else if (get_relative_location(at2) == get_relative_location(O)){
					internal_residue_bond_sep_mat[int_at1][int_at2] = 4;
					internal_residue_bond_sep_mat[int_at2][int_at1] = 4;
				}
			}
			else if (get_relative_location(at1) == get_relative_location(N)){
				if (get_relative_location(at2) == get_relative_location(H)){
					internal_residue_bond_sep_mat[int_at1][int_at2] = 1;
					internal_residue_bond_sep_mat[int_at2][int_at1] = 1;
				}
				else if (get_relative_location(at2) == get_relative_location(CA)){
					internal_residue_bond_sep_mat[int_at1][int_at2] = 1;
					internal_residue_bond_sep_mat[int_at2][int_at1] = 1;
				}
				else if (get_relative_location(at2) == get_relative_location(CB)){
					internal_residue_bond_sep_mat[int_at1][int_at2] = 2;
					internal_residue_bond_sep_mat[int_at2][int_at1] = 2;
				}
				else if (get_relative_location(at2) == get_relative_location(C)){
					internal_residue_bond_sep_mat[int_at1][int_at2] = 2;
					internal_residue_bond_sep_mat[int_at2][int_at1] = 2;
				}
				else if (get_relative_location(at2) == get_relative_location(O)){
					internal_residue_bond_sep_mat[int_at1][int_at2] = 3;
					internal_residue_bond_sep_mat[int_at2][int_at1] = 3;
				}
			}
			else if (get_relative_location(at1) == get_relative_location(CA)){
				if (get_relative_location(at2) == get_relative_location(CB)){
					internal_residue_bond_sep_mat[int_at1][int_at2] = 1;
					internal_residue_bond_sep_mat[int_at2][int_at1] = 1;
				}
				else if (get_relative_location(at2) == get_relative_location(C)){
					internal_residue_bond_sep_mat[int_at1][int_at2] = 1;
					internal_residue_bond_sep_mat[int_at2][int_at1] = 1;
				}
				else if (get_relative_location(at2) == get_relative_location(O)){
					internal_residue_bond_sep_mat[int_at1][int_at2] = 2;
					internal_residue_bond_sep_mat[int_at2][int_at1] = 2;
				}
			}
			else if (get_relative_location(at1) == get_relative_location(CB)){
				if (get_relative_location(at2) == get_relative_location(C)){
					internal_residue_bond_sep_mat[int_at1][int_at2] = 2;
					internal_residue_bond_sep_mat[int_at2][int_at1] = 2;
				}
				else if (get_relative_location(at2) == get_relative_location(O)){
					internal_residue_bond_sep_mat[int_at1][int_at2] = 3;
					internal_residue_bond_sep_mat[int_at2][int_at1] = 3;
				}
			}
			else if (get_relative_location(at1) == get_relative_location(C)){
				if (get_relative_location(at2) == get_relative_location(O)){
					internal_residue_bond_sep_mat[int_at1][int_at2] = 1;
					internal_residue_bond_sep_mat[int_at2][int_at1] = 1;
				}
			}

		}
	}

	for (BBAtomType at1 = N; at1 < END_BBAtomType; at1=static_cast<BBAtomType>(at1+1)){
		const int int_at1 = get_relative_location(at1);// static_cast<int>(at1) - init_val;
		const int int_n =  get_relative_location(N) ;//static_cast<int>(N) - init_val;
		const int int_c =  get_relative_location(C) ;//static_cast<int>(C) - init_val;

		bonds_to_N[int_at1] = internal_residue_bond_sep_mat[int_n][int_at1];
		bonds_to_C[int_at1] = internal_residue_bond_sep_mat[int_at1][int_c];


	}
}

int prot_backbone_map::get_bond_sep(const int seq_sep, const BBAtomType lower_at, const BBAtomType upper_at) const{
	//const int init_val = get_relative_location(N) ;//static_cast<int>(N);
	const int int_lower = get_relative_location(lower_at) ;//static_cast<int>(lower_at) - init_val;
	const int int_upper = get_relative_location(upper_at) ;//static_cast<int>(upper_at) - init_val;

	if (int_lower != num_protein_backbone_atoms && int_upper != num_protein_backbone_atoms){
		if (seq_sep == 0){
			return internal_residue_bond_sep_mat[int_lower][int_upper];
		}
		else{
			const int middle_bonds = 1 + ((abs(seq_sep) - 1) * 3);
			const int internal_lower = bonds_to_C[int_lower];
			const int internal_upper = bonds_to_N[int_upper];
			return middle_bonds + internal_lower + internal_upper;
		}
	}

	return 999;
}

int prot_backbone_map::get_bond_sep(const int seq_sep, const atom_type lower_at, const atom_type upper_at) const{
	const int int_lower = get_relative_location(lower_at) ;//static_cast<int>(lower_at) - init_val;
	const int int_upper = get_relative_location(upper_at) ;//static_cast<int>(upper_at) - init_val;

	if (int_lower != num_protein_backbone_atoms && int_upper != num_protein_backbone_atoms){
		if (seq_sep == 0){
			return internal_residue_bond_sep_mat[int_lower][int_upper];
		}
		else{
			const int middle_bonds = 1 + ((abs(seq_sep) - 1) * 3);
			const int internal_lower = bonds_to_C[int_lower];
			const int internal_upper = bonds_to_N[int_upper];
			return middle_bonds + internal_lower + internal_upper;
		}
	}
	return 999;
}






}
}
